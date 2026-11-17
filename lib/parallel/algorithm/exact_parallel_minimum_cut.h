/******************************************************************************
 * exact_parallel_minimum_cut.h
 *
 * Source of VieCut.
 *
 ******************************************************************************
 * Copyright (C) 2018 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#pragma once

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <deque>
#include <functional>
#include <memory>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include "algorithms/global_mincut/minimum_cut_helpers.h"
#include "algorithms/global_mincut/noi_minimum_cut.h"
#include "algorithms/global_mincut/viecut.h"
#include "common/configuration.h"
#include "common/definitions.h"
#include "data_structure/graph_access.h"
#include "data_structure/priority_queues/fifo_node_bucket_pq.h"
#include "data_structure/priority_queues/list_bucket_pq.h"
#include "data_structure/priority_queues/maxNodeHeap.h"
#include "data_structure/priority_queues/node_bucket_pq.h"
#include "data_structure/priority_queues/vecMaxNodeHeap.h"
#include "tools/random_functions.h"
#include "tools/timer.h"

#ifdef PARALLEL

#include "parallel/coarsening/contract_graph.h"
#include "parallel/coarsening/contraction_tests.h"
#include "parallel/coarsening/pr_cascade.h"
#include "parallel/coarsening/sparsify.h"
#include "parallel/data_structure/union_find.h"

#else
#include "coarsening/contract_graph.h"
#include "coarsening/contraction_tests.h"
#include "coarsening/sparsify.h"
#include "data_structure/union_find.h"
#endif



template <class GraphPtr>
class exact_parallel_minimum_cut : public minimum_cut {
 public:
    typedef GraphPtr GraphPtrType;
    exact_parallel_minimum_cut() { }
    ~exact_parallel_minimum_cut() { }

    static constexpr bool debug = false;
    bool timing = configuration::getConfig()->verbose;

    EdgeWeight perform_minimum_cut(GraphPtr G) {
        return perform_minimum_cut(G, false);
    }

    EdgeWeight perform_minimum_cut(GraphPtr G,
                                   bool indirect) {
        if (!G) {
            return -1;
        }

        std::vector<GraphPtr> graphs;
        timer t;
        EdgeWeight mincut = G->getMinDegree();
#ifdef PARALLEL
        // On a small sparse graph the pre-run buys nothing the minimum degree
        // does not already give: on every sparse family here the trivial
        // bound is exact, and the heuristic's fixed cost exceeds the whole
        // solve. The 1997 comparison made this visible: its ttime on a small
        // bikewheel is under a millisecond, and our floor was the pre-run.
        const bool skip_prerun =
            G->number_of_edges() < (1u << 19) &&
            !configuration::getConfig()->save_cut;
        if (!skip_prerun) {
        viecut<GraphPtr> heuristic_mc;
        mincut = heuristic_mc.perform_minimum_cut(G, true);
        LOGC(timing) << "VieCut found cut " << mincut
                     << " [Time: " << t.elapsed() << "s]";
        // On a graph small enough that the heuristic never contracted, it ran
        // the exact algorithm on the graph itself, so the bound it returned is
        // already the minimum cut and the loop below can only rediscover it.
        // With save_cut the loop is still needed, since the sides of the cut
        // are recovered from the contraction history.
        if (heuristic_mc.resultIsExact() &&
            !configuration::getConfig()->save_cut) {
            return mincut;
        }
        }
#endif

        graphs.push_back(G);

        // if PARALLEL is set, NodeInCut are already set to the result of viecut
        // This is what we want.
#ifndef PARALLEL
        minimum_cut_helpers<GraphPtr>::setInitialCutValues(graphs);
#endif

        int passes = omp_get_max_threads();
        // Sharing one visited array across the passes is work sharing when
        // the passes skip each other's vertices, and cache line ping-pong when
        // they do not. Which of the two it is follows the size of the graph.
        //
        // Every pass walks the whole graph on its own when the arrays are
        // separate, so the round reads the edges once per thread. While the
        // edges fit in the last level cache that costs nothing and the sharing
        // only buys contention, and independent arrays win by a factor of two
        // to twenty on every generated family. Once the edges are several
        // times the cache the same round reads tens of gigabytes from memory
        // per pass, and dividing the walk between the passes is the only thing
        // that matters: com-orkut's first round takes 10.6 seconds apart and
        // about 4 shared, and twitter-2010 four times over.
        //
        // This machine has 256 MiB of L3, and the two groups sit far apart:
        // the largest graph that prefers separate arrays needs 659 MiB and the
        // smallest that prefers a shared one needs 5134 MiB. The threshold is
        // 1 GiB, in the gap.
        //
        // It has to be decided before the first round rather than probed,
        // because on these graphs the first round is the whole run: it leaves
        // com-orkut with 3743 of its 2.4 million vertices, and a probe over
        // the rounds after it measures a graph that no longer costs anything.
        // 45 million directed entries, at 24 bytes each, is the 1 GiB.
        int pr_cascade_misses = 0;
        int pr_cascade_budget = 8;
        int pr_batch_misses = 0;
        NodeID last_round_merges = 0;
        bool share_visited =
            graphs.back()->number_of_edges() > (size_t{45} << 20);
        int round_idx = 0;
        double rate_shared = 0.0, rate_alone = 0.0;
        while (graphs.back()->number_of_nodes() > 2 && mincut > 0) {
            GraphPtr curr_g = graphs.back();

            // Reductions collect merges here and the round applies them in
            // the same contraction as the traversal's. Materialising each
            // reduction separately rebuilt the whole graph per reduction per
            // round, which is what the previous version paid: on irreg, which
            // needs hundreds of rounds and where the reductions find a pair
            // or two each, the extra rebuilds cost a fifth of the run.
            union_find pre_merges(curr_g->number_of_nodes());
            bool have_pre = false;

            // Cascade the Padberg-Rinaldi reductions before the traversal.
            // They find nothing on the raw input of these families and
            // thousands of contractions once the capforest has merged a few
            // vertices, so they run at the top of every round rather than
            // once. The cascade itself is serial, which the size guard pays
            // for; the merges it returns are applied through the same
            // parallel contraction as the capforest's. It stops asking after
            // a few empty rounds, since a graph on which it never fires only
            // pays the extraction.
            // The cascade must earn its call. On the families it exists
            // for it removes a large fraction of the graph every time; on
            // reg2, where it finds a handful of pairs per round, its per-round
            // extraction summed to sixteen times the whole solve. A call that
            // merges less than a hundredth of the graph is a strike, and three
            // strikes retire it for the run. The cascade is implemented for
            // mutable_graph only.
            constexpr bool is_mutable =
                std::is_same_v<GraphPtr, mutableGraphPtr>;
            if (is_mutable && pr_cascade_misses < 3 && pr_cascade_budget > 0 &&
                curr_g->number_of_nodes() < (1u << 17) &&
                curr_g->number_of_edges() <=
                    size_t{32} * curr_g->number_of_nodes()) {
                EdgeWeight pr_bound = mincut;
                union_find pruf(curr_g->number_of_nodes());
                size_t got = 0;
                if constexpr (is_mutable) {
                    got = pr_cascade::run(curr_g, &pr_bound, &pruf);
                }
                mincut = std::min(mincut, pr_bound);
                if (got > 0) {
                    pre_merges = std::move(pruf);
                    have_pre = true;
                }
                if (got * 100 < curr_g->number_of_nodes()) {
                    ++pr_cascade_misses;
                }
                // A call that collapses a quarter of the graph is the mode the
                // cascade exists for and stays free. A call that merely trims
                // draws on a fixed budget: on irreg, which needs hundreds of
                // rounds and where every call trims a little, per-round
                // re-extraction summed to twenty-five times the whole solve,
                // and the yield rule alone never fired.
                if (got * 4 < curr_g->number_of_nodes()) {
                    --pr_cascade_budget;
                }
            }

            // The dense complement of the cascade: the batch tests, parallel,
            // at the same equality boundary the 1997 code uses. On the dense
            // families the contractions the reductions find are mostly
            // disjoint pairs, which is exactly what a batch pass collects,
            // and the pass is a parallel sweep over the edges with no scratch
            // structure to build. Equality is sound here for the same reason
            // as in the cascade: every degree of the graph being tested is
            // already folded into the bound, by updateCut on the previous
            // round's contraction, and for the first round by the fold below.
            // The batch pass stays alive while the graph keeps changing.
            // It fires on fresh degrees, and the first rounds are exactly the
            // ones where there are none: on random the 1997 code's initial
            // pass also finds nothing, and its 872 thousand contractions all
            // come from passes made during contraction. Retiring after two
            // empty rounds, as this first did, threw that away.
            if (last_round_merges * 100 >= curr_g->number_of_nodes() &&
                pr_batch_misses < 8) {
                pr_batch_misses = 0;
            }
            if (pr_batch_misses < 2 &&
                curr_g->number_of_edges() < (size_t{1} << 25) &&
                curr_g->number_of_nodes() > 2) {
                curr_g->computeDegrees();
                EdgeWeight mindeg = mincut;
#pragma omp parallel for reduction(min : mindeg)
                for (NodeID v = 0; v < curr_g->number_of_nodes(); ++v) {
                    const EdgeWeight d = curr_g->getWeightedNodeDegree(v);
                    if (d > 0 && d < mindeg) mindeg = d;
                }
                mincut = std::min(mincut, mindeg);
                union_find uf12 =
                    tests::prTests12(curr_g, mincut, false, true);
                if (uf12.n() == curr_g->number_of_nodes()) {
                    // A pass that found nothing costs a sweep of the edges
                    // every round; two such rounds retire it.
                    ++pr_batch_misses;
                } else if (have_pre) {
                    for (NodeID v = 0; v < curr_g->number_of_nodes(); ++v) {
                        const NodeID r = uf12.Find(v);
                        if (r != v) pre_merges.Union(v, r);
                    }
                } else {
                    pre_merges = std::move(uf12);
                    have_pre = true;
                }
            }
            timer ts;
#ifdef PARALLEL

            noi_minimum_cut<GraphPtr> noi;

            // Each pass of the capforest is a full traversal of the graph
            // from its own starting vertex, and the passes are independent:
            // the work is replicated across threads, not divided among them.
            // That pays only when the extra starting vertices find extra
            // contractible edges. On the families built to resist contraction
            // a round merges one pair however many passes look for it, so
            // running one pass per thread costs a multiple of the work for the
            // same progress, and the algorithm is slower at 32 threads than at
            // one. So let the previous round decide: ask for about as many
            // passes as that round produced merges, never fewer than one and
            // never more than the threads available.
            const NodeID nodes_before = curr_g->number_of_nodes();
            // Round zero runs the default, and the two probe rounds follow
            // it. The first round works on the largest graph the run will ever
            // see, so it is the worst round to spend on a setting that might
            // turn out to be the bad one: on random it alone costs five
            // seconds of the eleven.
            const bool probing = (round_idx == 1 || round_idx == 2);
            const bool alone_now = probing ? (round_idx == 2) : !share_visited;
            timer t_cf;
            auto uf = parallel_modified_capforest(curr_g, mincut, alone_now,
                                                  passes);
            const double cf_seconds = t_cf.elapsed();
            if (have_pre) {
                // The reductions' merges ride the traversal's contraction, so
                // the round rebuilds the graph once whatever found the pairs.
                for (NodeID v = 0; v < curr_g->number_of_nodes(); ++v) {
                    const NodeID r = pre_merges.Find(v);
                    if (r != v) uf.Union(v, r);
                }
            }
            if (uf.n() == curr_g->number_of_nodes()) {
                uf = noi.modified_capforest(curr_g, mincut);
                LOGC(timing) << "seq capforest needed";
            }

#else
            LOG1 << "Error: Running exact_parallel_minimum_cut without PARALLEL"
                 << " Using normal noi_minimum_cut instead!";

            noi_minimum_cut noi;
            auto uf = noi.modified_capforest(curr_g, mincut);
#endif

            if (uf.n() > 1) {
                // Copy the graph once, then keep contracting that same copy in
                // place. Every round used to build a full new graph, which is
                // O(n+m) even when the round merges a single pair of vertices,
                // and on the families designed to resist contraction there is
                // one round per vertex: that is the Theta(n^2) the profile
                // shows, and it is also what makes the run retain one graph per
                // round.
                //
                // The first round must still copy, because mincut.cpp loads the
                // graph once and hands us the same object again for every -i
                // iteration; contracting it in place would corrupt the runs
                // that follow. Every later round works on a graph we made and
                // may modify. With save_cut the mapping back onto the input
                // needs the whole history, so there each round keeps its graph.
                const bool save_cut = configuration::getConfig()->save_cut;
                // The first round used to copy unconditionally, because the
                // driver reuses the input graph across repetitions. When there
                // is one repetition and the sides of the cut are not wanted,
                // nothing reads the input again, and the copy is a full
                // rebuild of the largest graph of the run: over a second on
                // noi1 at one thread.
                const bool own_copy =
                    save_cut ||
                    (graphs.size() == 1 &&
                     configuration::getConfig()->num_iterations > 1);
                GraphPtr next =
                    contraction::fromUnionFind(curr_g, &uf, own_copy);
                if (own_copy) {
                    graphs.push_back(next);
                } else {
                    graphs.back() = next;
                }

                mincut = minimum_cut_helpers<GraphPtr>::updateCut(
                    graphs, mincut);

                // Every pass runs on its own thread and they all run at once,
                // so a pass that finds nothing costs a core rather than any
                // wall clock time, while a pass that finds a merge removes a
                // whole round. The passes therefore stay at the thread count.
                //
                // Deriving the count from the merges the last round achieved,
                // which is what this did before, drives the very quantity it
                // measures: fewer passes start from fewer vertices, so they
                // find fewer merges, which asks for fewer passes again. On
                // bikewheel that loop settles at one pass and one merge per
                // round, and the run needs one round per vertex; at the full
                // thread count the same instance finishes in 0.6 seconds
                // rather than 15. The instances that prefer fewer passes are
                // the small dense ones, where the traversals compete for
                // memory bandwidth, and they lose at most about a seventh.
                const NodeID merges = nodes_before - uf.n();
                last_round_merges = merges;
                if (probing) {
                    // The fraction of the graph a round removes per second,
                    // not the number of vertices per second. The probe
                    // alternates the two settings over consecutive rounds, and
                    // consecutive rounds see a smaller graph than the one
                    // before, so counting vertices hands whichever setting
                    // draws the earlier rounds a lead it did not earn. On
                    // dblcyc that bias picked the shared array and cost a
                    // factor of four.
                    // Timed over the whole round, not over the traversal.
                    // Every round also rebuilds the contracted graph, which
                    // costs the same whether the round merged one vertex or a
                    // thousand, and on the cycle families that fixed cost is
                    // most of the round. Judging by the traversal alone picks
                    // the setting with the cheapest traversal even when it
                    // merges half as much per round and so needs twice as many
                    // rounds: on dblcyc that is a factor of four.
                    const double round_seconds = t_cf.elapsed();
                    const double rate =
                        (round_seconds > 0 && nodes_before > 0)
                            ? (static_cast<double>(merges) / nodes_before
                               / round_seconds)
                            : 0.0;
                    if (alone_now) {
                        rate_alone = rate;
                    } else {
                        rate_shared = rate;
                    }
                    if (round_idx == 2) {
                        share_visited = (rate_shared > rate_alone);
                    }
                }
                ++round_idx;
            } else {
                break;
            }
        }

        if (!indirect && configuration::getConfig()->save_cut)
            minimum_cut_helpers<GraphPtr>::retrieveMinimumCut(graphs);

        return mincut;
    }

    std::vector<NodeID> randomStartNodes(GraphPtr G, int passes) {
        std::vector<NodeID> start_nodes;
        for (int i = 0; i < passes; ++i)
            start_nodes.push_back(
                random_functions::next() % G->number_of_nodes());

        return start_nodes;
    }

    std::vector<NodeID> bfsStartNodes(GraphPtr G) {
        NodeID starting_node = random_functions::next() % G->number_of_nodes();
        std::vector<NodeID> start_nodes;
        start_nodes.push_back(starting_node);

        for (int i = 1; i < omp_get_max_threads(); ++i) {
            std::deque<NodeID> bfs;
            std::vector<bool> nodes(G->number_of_nodes(), false);
            size_t found = i;

            for (auto el : start_nodes) {
                bfs.push_back(el);
                nodes[el] = true;
            }

            while (!bfs.empty() && found < G->number_of_nodes()) {
                NodeID no = bfs.front();
                bfs.pop_front();
                for (EdgeID e : G->edges_of(no)) {
                    NodeID tgt = G->getEdgeTarget(e);
                    if (!nodes[tgt]) {
                        found++;
                        nodes[tgt] = true;
                        bfs.push_back(tgt);
                        if (found == G->number_of_nodes()) {
                            start_nodes.push_back(tgt);
                            break;
                        }
                    }
                }
            }
        }
        return start_nodes;
    }

    // The bucket priority queue is sized by the cut value: it holds
    // 2*(mincut+1)+1 std::deques, and libstdc++ gives every deque a 512 byte
    // node as soon as it is constructed. On an instance whose minimum cut is
    // large that is hundreds of millions of allocations for a graph with a few
    // thousand vertices, and it is paid again in every round of the outer loop.
    // A heap costs O(log n) per operation instead of O(1), but it is O(n) in
    // total and allocates nothing per unit of cut. The sequential
    // Nagamochi-Ibaraki code makes exactly this switch in
    // noi_minimum_cut::selectPq; this is the same rule.
    //
    // The choice cannot change the result. The forest only ever contracts edges
    // that provably do not cross a minimum cut, so a different order of equal
    // priorities contracts a different set of safe edges and still leaves the
    // minimum cut intact.
    union_find parallel_modified_capforest(
        GraphPtr G,
        const EdgeWeight mincut,
        const bool disable_blacklist = false,
        int passes = 0) {
        if (passes <= 0) passes = omp_get_max_threads();
        if (mincut > 10000 && mincut > G->number_of_nodes()) {
            return capforest<vecMaxNodeHeap>(G, mincut, disable_blacklist,
                                             passes);
        }
        // Between the two bucket queues, choose by how crowded a bucket will
        // be. A deque holds its contents contiguously, so when the cut is small
        // and a handful of buckets hold thousands of vertices each it is the
        // cache-friendly choice. When the cut is larger the buckets hold a few
        // vertices apiece, the deque's own machinery costs more than the
        // vertices it holds, and intrusive lists win: on reg2 they are a third
        // faster and on dblcyc nearly twice.
        const size_t buckets = 2 * static_cast<size_t>(mincut) + 3;
        const size_t per_bucket = G->number_of_nodes() / (buckets ? buckets : 1);
        if (per_bucket > 2000) {
            return capforest<fifo_node_bucket_pq>(G, mincut, disable_blacklist,
                                                  passes);
        }
        return capforest<list_bucket_pq>(G, mincut, disable_blacklist,
                                         passes);
    }

    // Templated on the queue so that the inner loop keeps calling it directly:
    // going through priority_queue_interface would put a virtual call on the
    // hottest path in the algorithm.
    template <typename PQ>
    union_find capforest(
        GraphPtr G,
        const EdgeWeight mincut,
        const bool disable_blacklist,
        int passes) {
        union_find uf(G->number_of_nodes());
        timer t;

        timer timer2;
        std::vector<NodeID> start_nodes = randomStartNodes(G, passes);

        // std::vector<bool> would be bad for thread-safety
        std::vector<uint8_t> visited(G->number_of_nodes(), false);

#pragma omp parallel for num_threads(passes)
        for (int i = 0; i < passes; ++i) {
            PQ pq(G->number_of_nodes(), mincut + 1);
            // Bytes rather than std::vector<bool>: these are probed once per
            // edge visit on the hottest path in the algorithm, and a bit vector
            // turns each probe into a load, a shift and a mask. The extra
            // memory is one byte per vertex.
            std::vector<NodeID> r_v(G->number_of_nodes(), 0);
            std::vector<uint8_t> local_visited(G->number_of_nodes(), 0);

            NodeID starting_node = start_nodes[i];
            NodeID current_node = starting_node;

            pq.insert(current_node, 0);

            timer t;
            size_t elements = 0;

            while (!pq.empty()) {
                current_node = pq.deleteMax();
                elements++;
                local_visited[current_node] = 1;

                if (!disable_blacklist) {
                    if (visited[current_node]) {
                        continue;
                    } else {
                        visited[current_node] = true;
                    }
                }

                const auto& adj = G->adjacencyOf(current_node);
                for (const auto& edge : adj) {
                    const NodeID tgt = edge.target;
                    const EdgeWeight wgt = edge.weight;
                    if (!local_visited[tgt]) {
                        if (r_v[tgt] < mincut) {
                            if ((r_v[tgt] + wgt) >= mincut) {
                                // The blacklist this used to test is the set of
                                // vertices already taken from the queue, which
                                // is what local_visited records, and it was
                                // written in the same place. Inside this branch
                                // the vertex is not in it, so the test was
                                // always true: an array of a byte per vertex
                                // allocated and cleared per traversal, and a
                                // random probe and a branch per edge visit, all
                                // of which decided nothing.
                                uf.Union(current_node, tgt);
                            }

                            if (!visited[tgt]) {
                                size_t new_rv =
                                    std::min(r_v[tgt] + wgt,
                                             mincut);
                                r_v[tgt] = new_rv;
                                // Both !visited[tgt] and !local_visited[tgt]
                                // are already known here, from this branch and
                                // from the one above it, so the queue update is
                                // unconditional; re-testing them cost two more
                                // random probes on every edge visit.
                                pq.insertOrIncrease(tgt, new_rv);
                            }
                        }
                    }
                }
            }
        }
        return uf;
    }
};
