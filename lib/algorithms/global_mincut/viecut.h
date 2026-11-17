/******************************************************************************
 * viecut.h
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
#include <vector>

#include "algorithms/flow/excess_scaling.h"
#include "algorithms/global_mincut/minimum_cut.h"
#include "algorithms/global_mincut/minimum_cut_helpers.h"
#include "algorithms/global_mincut/noi_minimum_cut.h"
#include "algorithms/misc/strongly_connected_components.h"
#include "common/definitions.h"
#include "data_structure/flow_graph.h"
#include "data_structure/graph_access.h"
#include "tlx/logger.hpp"
#include "tools/graph_extractor.h"
#include "tools/timer.h"

#ifdef PARALLEL
#include "parallel/coarsening/contract_graph.h"
#include "parallel/coarsening/contraction_tests.h"
#include "parallel/coarsening/label_propagation.h"
#else
#include "coarsening/contract_graph.h"
#include "coarsening/contraction_tests.h"
#include "coarsening/label_propagation.h"
#endif

template <class GraphPtr>
class viecut : public minimum_cut {
 public:
    typedef GraphPtr GraphPtrType;
    static constexpr bool debug = false;
    bool exact_result = false;
    bool timing = configuration::getConfig()->verbose;
    viecut() { }

    virtual ~viecut() { }

    EdgeWeight perform_minimum_cut(GraphPtr G) {
        return perform_minimum_cut(G, false);
    }

    // True when the value returned is the minimum cut of the graph handed in,
    // not merely an upper bound on it.
    bool resultIsExact() const { return exact_result; }

    EdgeWeight perform_minimum_cut(GraphPtr G,
                                   bool indirect) {
        if (!G) {
            return -1;
        }

        EdgeWeight cut = G->getMinDegree();
        std::vector<GraphPtr> graphs;
        graphs.push_back(G);
        exact_result = false;

        minimum_cut_helpers<GraphPtr>::setInitialCutValues(graphs);

        while (graphs.back()->number_of_nodes() > 10000 &&
               (graphs.size() == 1 ||
                (graphs.back()->number_of_nodes() <
                 graphs[graphs.size() - 2]->number_of_nodes()))) {
            timer t;
            G = graphs.back();
            label_propagation<GraphPtr> lp;
            std::vector<NodeID> cluster_mapping = lp.propagate_labels(G);
            auto [mapping, reverse_mapping] =
                minimum_cut_helpers<GraphPtr>::remap_cluster(
                    G, cluster_mapping);
            LOGC(timing) << "LP (total): " << t.elapsedToZero();

            contraction::findTrivialCuts(G, &mapping, &reverse_mapping, cut);
            LOGC(timing) << "Trivial Cut Local Search: " << t.elapsedToZero();

            auto H = contraction::contractGraph(G, mapping, reverse_mapping);
            graphs.push_back(H);
            cut = minimum_cut_helpers<GraphPtr>::updateCut(graphs, cut);
            LOGC(timing) << "Graph Contraction (to "
                         << graphs.back()->number_of_nodes()
                         << " nodes): " << t.elapsedToZero();

            union_find uf = tests::prTests12(graphs.back(), cut);
            graphs.push_back(
                contraction::fromUnionFind(graphs.back(), &uf, true));
            cut = minimum_cut_helpers<GraphPtr>::updateCut(graphs, cut);
            union_find uf2 = tests::prTests34(graphs.back(), cut);
            graphs.push_back(
                contraction::fromUnionFind(graphs.back(), &uf2, true));
            cut = minimum_cut_helpers<GraphPtr>::updateCut(graphs, cut);
            LOGC(timing) << "Padberg-Rinaldi Tests (to "
                         << graphs.back()->number_of_nodes()
                         << " nodes): " << t.elapsedToZero();
        }

        if (graphs.back()->number_of_nodes() > 1) {
            timer t;
            // The loop above only runs on graphs of more than 10000 vertices,
            // so below that it never contracts and the call below solves the
            // graph we were given rather than a coarsened one.
            const bool never_contracted = (graphs.size() == 1);
#ifdef PARALLEL
            // That call is sequential. When this is the pre-run of the parallel
            // exact algorithm, solving the whole graph here spends one thread
            // and leaves the other thirty-one idle, and then reports the answer
            // as exact so that the parallel loop never runs at all: every
            // instance of ten thousand vertices or fewer, which is fifteen of
            // the fifty-three here, ignored the thread count entirely. Hand the
            // graph back instead. The minimum degree is already an upper bound
            // on the minimum cut, which is all the caller's loop needs to
            // start, and it solves the graph with every thread. reg2 at
            // n=6400 goes from 22.6 seconds to 1.6.
            if (indirect && never_contracted && omp_get_max_threads() > 1 &&
                !configuration::getConfig()->save_cut) {
                exact_result = false;
                return cut;
            }
#endif
            noi_minimum_cut<GraphPtr> noi;
            // Solving the graph itself makes the answer the minimum cut rather
            // than an upper bound on it. Say so, so that the caller need not
            // compute it a second time.
            exact_result = never_contracted;
            cut = std::min(cut, noi.perform_minimum_cut(graphs.back(), true));

            LOGC(timing) << "Exact Algorithm:"
                         << t.elapsedToZero() << " deg: "
                         << graphs.back()->getMinDegree();
        }

        if (!indirect && configuration::getConfig()->save_cut)
            minimum_cut_helpers<GraphPtr>::retrieveMinimumCut(graphs);

        return cut;
    }
};
