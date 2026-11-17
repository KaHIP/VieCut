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

#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <memory>
#include <vector>

#include "data_structure/adjlist_graph.h"
#include "data_structure/priority_queues/bucket_pq.h"
#include "data_structure/priority_queues/priority_queue_interface.h"
#include "data_structure/priority_queues/vecMaxNodeHeap.h"
#ifdef PARALLEL
#include "parallel/data_structure/union_find.h"
#else
#include "data_structure/union_find.h"
#endif
#include "tools/random_functions.h"

class stoer_wagner_minimum_cut {
 public:
    stoer_wagner_minimum_cut() { }
    virtual ~stoer_wagner_minimum_cut() { }

    /**
     * Stoer-Wagner: n - 1 maximum adjacency phases, each of which yields a
     * cut-of-the-phase and merges the last two vertices it added. Returns 0
     * for a disconnected graph and UINT64_MAX for a graph with fewer than
     * two vertices.
     */
    uint64_t perform_minimum_cut(adjlist_graph G) {
        uint64_t global_mincut = UINT64_MAX;
        union_find uf(G.number_of_nodes());
        for (uint64_t i = G.number_of_nodes(); i >= 2; --i) {
            uint64_t mincut = contract(G, i, uf);
            if (mincut < global_mincut) {
                global_mincut = mincut;
            }
            if (global_mincut == 0) {
                break;
            }
        }
        return global_mincut;
    }

 private:
    // One phase on the graph with size (super)vertices left. G and uf are
    // updated in place, so that the merge carries over to the next phase.
    uint64_t contract(adjlist_graph& G, uint64_t size, union_find& uf) {
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
                // Not every vertex is reachable: the graph is disconnected.
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
            return 0;
        }

        // Merge the last two vertices of the phase: current_node was removed
        // last, last_node is the only one left in the queue.
        NodeID last_node = pq->maxElement();
        uint64_t cut_of_the_phase = pq->maxValue();
        uf.Union(current_node, last_node);
        if (uf.Find(current_node) == current_node) {
            G.merge_nodes(current_node, last_node, uf);
        } else {
            G.merge_nodes(last_node, current_node, uf);
        }
        return cut_of_the_phase;
    }
};
