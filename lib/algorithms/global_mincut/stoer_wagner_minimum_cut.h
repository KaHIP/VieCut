/******************************************************************************
 * stoer_wagner_minimum_cut.h
 *
 * Source of VieCut.
 *
 ******************************************************************************
 * Copyright (C) 2017 Alexander Noe <alexander.noe@univie.ac.at>
 * Correctness fixes by Bora Ucar <bora.ucar@ens-lyon.fr>.
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#pragma once

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <memory>
#include <vector>

#include "algorithms/global_mincut/minimum_cut.h"
#include "common/configuration.h"
#include "data_structure/adjlist_graph.h"
#include "data_structure/graph_access.h"
#include "data_structure/mutable_graph.h"
#include "data_structure/priority_queues/bucket_pq.h"
#include "data_structure/priority_queues/priority_queue_interface.h"
#include "data_structure/priority_queues/vecMaxNodeHeap.h"
#ifdef PARALLEL
#include "parallel/data_structure/union_find.h"
#else
#include "data_structure/union_find.h"
#endif
#include "tools/random_functions.h"

class stoer_wagner_minimum_cut : public minimum_cut {
 public:
    stoer_wagner_minimum_cut() { }
    virtual ~stoer_wagner_minimum_cut() { }

    using minimum_cut::perform_minimum_cut;

    // Entry points of the minimum_cut interface (mincut <graph> sw). With
    // save_cut set they also mark the side of the cut (setNodeInCut).
    EdgeWeight perform_minimum_cut(graphAccessPtr G) override {
        return solveAndMark(G);
    }

    EdgeWeight perform_minimum_cut(mutableGraphPtr G) override {
        return solveAndMark(G);
    }

    /**
     * Stoer-Wagner: n - 1 maximum adjacency phases, each of which yields a
     * cut-of-the-phase and merges the last two vertices it added. Returns 0
     * for a disconnected graph and UINT64_MAX for a graph with fewer than
     * two vertices. If side is given, it is set to one side of the cut,
     * indexed by vertex.
     */
    uint64_t perform_minimum_cut(adjlist_graph G,
                                 std::vector<bool>* side = nullptr) {
        uint64_t global_mincut = UINT64_MAX;
        union_find uf(G.number_of_nodes());
        if (side) {
            side->assign(G.number_of_nodes(), false);
        }
        for (uint64_t i = G.number_of_nodes(); i >= 2; --i) {
            global_mincut = contract(G, i, uf, global_mincut, side);
            if (global_mincut == 0) {
                break;
            }
        }
        return global_mincut;
    }

 private:
    template <class GraphPtr>
    EdgeWeight solveAndMark(GraphPtr G) {
        if (!G) {
            return -1;
        }
        if (!configuration::getConfig()->save_cut) {
            return perform_minimum_cut(toAdjlist(G));
        }
        std::vector<bool> side;
        EdgeWeight cut = perform_minimum_cut(toAdjlist(G), &side);
        for (NodeID n : G->nodes()) {
            G->setNodeInCut(n, side[n]);
        }
        return cut;
    }

    template <class GraphPtr>
    static adjlist_graph toAdjlist(GraphPtr G) {
        adjlist_graph A;
        A.start_construction(G->number_of_nodes(), G->number_of_edges());
        for (NodeID n : G->nodes()) {
            A.new_node();
            for (EdgeID e : G->edges_of(n)) {
                auto [tgt, wgt] = G->getEdge(n, e);
                A.new_edge(n, tgt, wgt);
            }
        }
        A.finish_construction();
        return A;
    }

    // Sets side to the vertices whose supervertex satisfies in_side.
    template <class InSide>
    static void markSide(union_find& uf, std::vector<bool>* side,
                         InSide in_side) {
        for (NodeID v = 0; v < side->size(); ++v) {
            (*side)[v] = in_side(uf.Find(v));
        }
    }

    // One phase on the graph with size (super)vertices left. G and uf are
    // updated in place, so that the merge carries over to the next phase.
    // Returns the smaller of best and the cut of the phase, and records the
    // side of the cut in side if the phase improved on best.
    uint64_t contract(adjlist_graph& G, uint64_t size, union_find& uf,
                      uint64_t best, std::vector<bool>* side) {
        assert(size >= 2);
        // A bucket queue holds one bucket per attainable key, so it is only
        // used while the keys (bounded by the maximum degree) stay below n.
        std::unique_ptr<priority_queue_interface> pq;
        if (G.getMaxDegree() < G.number_of_nodes()) {
            pq = std::make_unique<bucket_pq>(G.getMaxDegree());
        } else {
            pq = std::make_unique<vecMaxNodeHeap>(G.number_of_nodes());
        }

        std::vector<bool> visited(G.number_of_nodes(), false);
        NodeID starting_node =
            random_functions::nextInt(0, G.number_of_nodes() - 1);
        // The start vertex enters with key 0; it is the first one removed.
        pq->insert(uf.Find(starting_node), 0);

        NodeID current_node = UNDEFINED_NODE;
        for (uint64_t i = 0; i < size - 1; ++i) {
            if (pq->empty()) {
                // Not every vertex is reachable: the graph is disconnected,
                // and the reachable part is a side of a cut of weight 0.
                if (side) {
                    markSide(uf, side, [&](NodeID r) { return visited[r]; });
                }
                return 0;
            }
            current_node = pq->deleteMax();
            assert(uf.Find(current_node) == current_node);
            visited[current_node] = true;

            for (EdgeID e = G.get_first_edge(current_node);
                 e < G.get_first_invalid_edge(current_node); ++e) {
                NodeID tgt = uf.Find(G.getEdgeTarget(e));
                if (!visited[tgt]) {
                    assert(G.getEdgeWeight(e) > 0);
                    if (pq->contains(tgt)) {
                        pq->increaseKey(tgt, pq->getKey(tgt)
                                        + G.getEdgeWeight(e));
                    } else {
                        pq->insert(tgt, G.getEdgeWeight(e));
                    }
                }
            }
        }

        if (pq->empty()) {
            if (side) {
                markSide(uf, side, [&](NodeID r) { return visited[r]; });
            }
            return 0;
        }

        // last_node, the only vertex left in the queue, is the last one added;
        // the cut of the phase separates it from the rest. current_node was
        // added just before it, and the two are merged.
        NodeID last_node = pq->maxElement();
        uint64_t cut_of_the_phase = pq->maxValue();
        if (side && cut_of_the_phase < best) {
            markSide(uf, side, [&](NodeID r) { return r == last_node; });
        }
        uf.Union(current_node, last_node);
        if (uf.Find(current_node) == current_node) {
            G.merge_nodes(current_node, last_node, uf);
        } else {
            G.merge_nodes(last_node, current_node, uf);
        }
        return std::min(best, cut_of_the_phase);
    }
};
