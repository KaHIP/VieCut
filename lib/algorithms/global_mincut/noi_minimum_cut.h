/******************************************************************************
 * noi_minimum_cut.h
 *
 * Source of VieCut
 *
 ******************************************************************************
 * Copyright (C) 2017-2018 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#pragma once

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "algorithms/global_mincut/minimum_cut.h"
#include "algorithms/global_mincut/minimum_cut_helpers.h"
#include "algorithms/multicut/multicut_problem.h"
#include "common/configuration.h"
#include "common/definitions.h"
#include "data_structure/graph_access.h"
#include "data_structure/mutable_graph.h"
#include "data_structure/priority_queues/bucket_pq.h"
#include "data_structure/priority_queues/fifo_node_bucket_pq.h"
#include "data_structure/priority_queues/maxNodeHeap.h"
#include "data_structure/priority_queues/node_bucket_pq.h"
#include "data_structure/priority_queues/vecMaxNodeHeap.h"
#include "tools/random_functions.h"
#include <type_traits>
#include "tools/timer.h"

#ifdef PARALLEL
#include "parallel/coarsening/contract_graph.h"
#include "parallel/coarsening/contraction_tests.h"
#include "parallel/coarsening/pr_cascade.h"
#include "parallel/data_structure/union_find.h"
#else
#include "coarsening/contract_graph.h"
#include "coarsening/contraction_tests.h"
#include "data_structure/union_find.h"
#endif

template <class GraphPtr>
class noi_minimum_cut : public minimum_cut {
 public:
    typedef GraphPtr GraphPtrType;
    noi_minimum_cut() { }
    ~noi_minimum_cut() { }
    static constexpr bool debug = false;
    static constexpr bool timing = true;

    EdgeWeight perform_minimum_cut(GraphPtr G) {
        return perform_minimum_cut(G, false);
    }

    EdgeWeight perform_minimum_cut(GraphPtr G, bool indirect) {
        if (!G) {
            return -1;
        }

        std::vector<GraphPtr> graphs;
        timer t;
        EdgeWeight mincut = G->getMinDegree();
        graphs.push_back(G);
        minimum_cut_helpers<GraphPtr>::setInitialCutValues(graphs);

        int pr_cascade_misses = 0;
        int pr_cascade_budget = 8;
        (void)pr_cascade_misses; (void)pr_cascade_budget;
        while (graphs.back()->number_of_nodes() > 2 && mincut > 0) {
            // The same cascading Padberg-Rinaldi reductions the parallel loop
            // applies, folded into this round's contraction. This loop is what
            // the pre-run uses to finish every instance, so without them a
            // small bikewheel takes four seconds at one thread where the
            // reductions solve it in milliseconds. Guards as in the parallel
            // loop: sparse graphs, three unproductive calls retire it.
            union_find seq_pre(graphs.back()->number_of_nodes());
            bool have_seq_pre = false;
#ifdef PARALLEL
            if constexpr (std::is_same_v<GraphPtr, mutableGraphPtr>) {
            if (pr_cascade_misses < 3 && pr_cascade_budget > 0 &&
                graphs.back()->number_of_nodes() < (1u << 17) &&
                graphs.back()->number_of_edges() <=
                    size_t{32} * graphs.back()->number_of_nodes()) {
                EdgeWeight pr_bound = mincut;
                const size_t got =
                    pr_cascade::run(graphs.back(), &pr_bound, &seq_pre);
                mincut = std::min(mincut, pr_bound);
                if (got > 0) have_seq_pre = true;
                if (got * 100 < graphs.back()->number_of_nodes()) {
                    ++pr_cascade_misses;
                }
                if (got * 4 < graphs.back()->number_of_nodes()) {
                    --pr_cascade_budget;
                }
            }
            }
#endif
            auto uf = modified_capforest(graphs.back(), mincut);
            if (have_seq_pre) {
                for (NodeID v = 0; v < graphs.back()->number_of_nodes(); ++v) {
                    const NodeID r = seq_pre.Find(v);
                    if (r != v) uf.Union(v, r);
                }
            }
            // Copy the graph once and then contract that copy in place, rather
            // than building a whole new graph in every round and keeping all of
            // them. The first round still copies, because the caller keeps
            // using the graph it handed us, and with save_cut the history is
            // needed to map the cut back onto the input.
            const bool save_cut = configuration::getConfig()->save_cut;
            const bool own_copy =
                save_cut ||
                (graphs.size() == 1 &&
                 configuration::getConfig()->num_iterations > 1);
            GraphPtr next =
                contraction::fromUnionFind(graphs.back(), &uf, own_copy);
            if (own_copy) {
                graphs.emplace_back(next);
            } else {
                graphs.back() = next;
            }
            mincut = minimum_cut_helpers<GraphPtr>::updateCut(graphs, mincut);
        }

        if (!indirect && configuration::getConfig()->save_cut)
            minimum_cut_helpers<GraphPtr>::retrieveMinimumCut(graphs);

        return mincut;
    }

    static priority_queue_interface * selectPq(GraphPtr G,
                                               EdgeWeight mincut) {
        priority_queue_interface* pq;

        size_t max_deg = G->getMaxDegree();

        if (configuration::getConfig()->disable_limiting) {
            mincut = max_deg;
        }

        if (configuration::getConfig()->pq == "default") {
            if (mincut > 10000 && mincut > G->number_of_nodes()) {
                pq = new vecMaxNodeHeap(G->number_of_nodes());
            } else {
                pq = new fifo_node_bucket_pq(G->number_of_nodes(), mincut);
            }
        } else {
            if (configuration::getConfig()->pq == "bqueue") {
                pq = new fifo_node_bucket_pq(G->number_of_nodes(), mincut);
            } else {
                if (configuration::getConfig()->pq == "heap") {
                    pq = new vecMaxNodeHeap(G->number_of_nodes());
                } else {
                    if (configuration::getConfig()->pq == "bstack") {
                        pq = new node_bucket_pq(G->number_of_nodes(), mincut);
                    } else {
                        std::cerr << "unknown pq type "
                                  << configuration::getConfig()->pq
                                  << std::endl;
                        exit(1);
                        pq = new maxNodeHeap();
                    }
                }
            }
        }
        return pq;
    }

    union_find modified_capforest(GraphPtr G,
                                  EdgeWeight mincut) {
        union_find uf(G->number_of_nodes());

        priority_queue_interface* pq = selectPq(G, mincut);

        std::vector<bool> visited(G->number_of_nodes(), false);
        std::vector<bool> seen(G->number_of_nodes(), false);
        std::vector<EdgeWeight> r_v(G->number_of_nodes(), 0);

        NodeID starting_node = random_functions::next() % G->number_of_nodes();

        NodeID current_node = starting_node;

        pq->insert(current_node, 0);

        while (!pq->empty()) {
            current_node = pq->deleteMax();
            visited[current_node] = true;
            if (!configuration::getConfig()->disable_limiting) {
                for (EdgeID e : G->edges_of(current_node)) {
                    NodeID tgt = G->getEdgeTarget(current_node, e);
                    if (!visited[tgt]) {
                        bool increase = false;

                        if (r_v[tgt] < mincut || mincut == 0) {
                            increase = true;
                            if ((r_v[tgt] + G->getEdgeWeight(
                                     current_node, e)) >= mincut) {
                                uf.Union(current_node, tgt);
                            }
                        }

                        r_v[tgt] += G->getEdgeWeight(current_node, e);

                        size_t new_rv = std::min(r_v[tgt], mincut);

                        if (seen[tgt]) {
                            if (increase && !visited[tgt]) {
                                pq->increaseKey(tgt, new_rv);
                            }
                        } else {
                            seen[tgt] = true;
                            pq->insert(tgt, new_rv);
                        }
                    }
                }
            } else {
                for (EdgeID e : G->edges_of(current_node)) {
                    NodeID tgt = G->getEdgeTarget(current_node, e);

                    if (!visited[tgt]) {
                        if (r_v[tgt] < mincut) {
                            if ((r_v[tgt] + G->getEdgeWeight(
                                     current_node, e)) >= mincut) {
                                uf.Union(current_node, tgt);
                            }
                        }

                        r_v[tgt] += G->getEdgeWeight(current_node, e);

                        if (seen[tgt]) {
                            if (!visited[tgt]) {
                                pq->increaseKey(tgt, r_v[tgt]);
                            }
                        } else {
                            seen[tgt] = true;
                            pq->insert(tgt, r_v[tgt]);
                        }
                    }
                }
            }
        }
        delete pq;
        return uf;
    }
};
