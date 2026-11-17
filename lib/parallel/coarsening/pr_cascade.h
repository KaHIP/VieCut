/******************************************************************************
 * pr_cascade.h
 *
 * The Padberg-Rinaldi reductions, applied so that they cascade.
 *
 * The tests only become able to fire once contraction has raised the degrees
 * they compare against: on every family here their initial pass finds nothing
 * and the passes made while the graph contracts find thousands. The 1997
 * implementation of Chekuri et al. rides that chain by contracting on the spot
 * and re-reading the degrees; on random graphs of a million vertices its
 * internal passes remove 87 percent of the graph. This is that mechanism, with
 * a representation chosen so that it stays cheap at that size.
 *
 * The graph is held as one plain vector of (target, weight) entries per
 * supervertex. Entries are allowed to go stale: an entry pointing at a merged
 * vertex is redirected through the union-find when it is read, and parallel
 * entries are combined only when their owner is tested, by a stamped
 * aggregation pass over its list. A merge splices the shorter list onto the
 * longer and updates one degree arithmetically; it touches neither the
 * neighbours nor the longer list. An earlier version kept one hash map per
 * supervertex instead, and the maps were the cost: admitting a graph of a
 * million vertices meant a million allocations per call, and the run was
 * slower than leaving the reductions off entirely.
 *
 * Two soundness points, both learnt the hard way.
 *
 * The tests fire at equality, as the 1997 code has them, and equality is only
 * sound if the degree of every vertex the cascade removes is first recorded as
 * a cut: contracting u into v because twice the edge weight reaches deg(u)
 * argues that u can cross the cut, and that argument fails for the cut that
 * isolates u. The bound this maintains is part of the correctness.
 *
 * The fourth test sums min(w2, w3) over the common neighbourhood, and the lazy
 * representation can hold the weight between w and a common neighbour split
 * across several stale entries. Taking the minimum per entry overstates the
 * sum: min(3,5) + min(3,5) = 6 where the true min(6,5) = 5, and an overstated
 * sum contracts an edge the test does not justify. The inner neighbourhood is
 * therefore aggregated before any minimum is taken.
 *
 * The cascade itself is sequential and near-linear; the extraction is
 * parallel, and the merges it returns are applied through the caller's
 * parallel contraction.
 *****************************************************************************/

#pragma once

#include <algorithm>
#include <utility>
#include <vector>

#include "common/definitions.h"
#include "data_structure/mutable_graph.h"
#include "parallel/data_structure/union_find.h"

class pr_cascade {
 public:
    // Reduces G, recording the merges in *out and the strongest upper bound
    // seen in *bound. Returns the number of merges. G is not modified.
    static size_t run(mutableGraphPtr G, EdgeWeight* bound, union_find* out) {
        const NodeID n = G->number_of_nodes();
        if (n <= 2) return 0;

        std::vector<std::vector<std::pair<NodeID, EdgeWeight> > > adj(n);
        std::vector<EdgeWeight> deg(n, 0);

#pragma omp parallel for schedule(dynamic, 512)
        for (NodeID v = 0; v < n; ++v) {
            const auto& row = G->adjacencyOf(v);
            adj[v].reserve(row.size());
            EdgeWeight d = 0;
            for (const auto& e : row) {
                if (e.target == v) continue;
                adj[v].emplace_back(e.target, e.weight);
                d += e.weight;
            }
            deg[v] = d;
        }

        // Every degree is a cut; the smallest is a bound before any merge.
        for (NodeID v = 0; v < n; ++v) {
            if (deg[v] > 0 && deg[v] < *bound) *bound = deg[v];
        }

        std::vector<NodeID> parent(n);
        for (NodeID v = 0; v < n; ++v) parent[v] = v;
        auto find = [&parent](NodeID v) {
            while (parent[v] != v) {
                parent[v] = parent[parent[v]];
                v = parent[v];
            }
            return v;
        };

        // Stamped aggregation scratch: one pair of arrays for the vertex under
        // test, a second pair for the neighbour the triangle tests walk, since
        // both neighbourhoods are needed at once.
        std::vector<uint32_t> st1(n, 0), st2(n, 0);
        std::vector<EdgeWeight> ac1(n, 0), ac2(n, 0);
        uint32_t ep = 0;
        std::vector<std::pair<NodeID, EdgeWeight> > nb2;

        // Combines a list's entries by current representative, dropping self
        // loops; writes the distinct neighbours back into the list.
        auto compact = [&](NodeID v, std::vector<uint32_t>& st,
                           std::vector<EdgeWeight>& ac, uint32_t epoch,
                           std::vector<std::pair<NodeID, EdgeWeight> >* into) {
            into->clear();
            for (const auto& [t, w] : adj[v]) {
                const NodeID x = find(t);
                if (x == v) continue;
                if (st[x] != epoch) {
                    st[x] = epoch;
                    ac[x] = w;
                    into->emplace_back(x, 0);
                } else {
                    ac[x] += w;
                }
            }
            for (auto& e : *into) e.second = ac[e.first];
        };

        std::vector<uint8_t> queued(n, 1);
        std::vector<NodeID> work(n);
        for (NodeID v = 0; v < n; ++v) work[v] = n - 1 - v;
        std::vector<std::pair<NodeID, EdgeWeight> > nb1;

        size_t merges = 0;
        while (!work.empty()) {
            NodeID v = work.back();
            work.pop_back();
            queued[v] = 0;
            if (find(v) != v || adj[v].empty()) continue;

            // Aggregate v's neighbourhood once; st1/ac1 stay valid for the
            // triangle tests below.
            if (++ep == 0) {
                std::fill(st1.begin(), st1.end(), 0);
                std::fill(st2.begin(), st2.end(), 0);
                ep = 1;
            }
            compact(v, st1, ac1, ep, &nb1);
            adj[v].assign(nb1.begin(), nb1.end());

            NodeID partner = UNDEFINED_NODE;
            EdgeWeight w_link = 0;
            const EdgeWeight dv = deg[v];

            for (const auto& [x, w] : nb1) {
                // Test 1: the edge alone reaches the bound.
                // Test 2: the edge carries at least half of either endpoint.
                if ((w >= *bound) ||
                    ((2 * w >= dv || 2 * w >= deg[x]))) {
                    partner = x;
                    w_link = w;
                    break;
                }
            }

            if (partner == UNDEFINED_NODE && nb1.size() <= 32) {
                for (const auto& [x, w1] : nb1) {
                    if (adj[x].size() > 256) continue;
                    // Aggregate x's side before taking any minimum: split
                    // stale entries would overstate the shared sum.
                    compact(x, st2, ac2, ep, &nb2);
                    EdgeWeight shared = w1;
                    bool hit = false;
                    for (const auto& [y, w2] : nb2) {
                        if (y == v || st1[y] != ep) continue;
                        const EdgeWeight w3 = ac1[y];   // v -- y
                        shared += std::min(w2, w3);
                        // Test 3: the triangle covers half of both endpoints.
                        if (2 * (w1 + w3) >= dv &&
                            2 * (w1 + w2) >= deg[x]) {
                            hit = true;
                            break;
                        }
                    }
                    // Test 4: the edge plus the shared neighbourhood reaches
                    // the bound.
                    if (hit || shared >= *bound) {
                        partner = x;
                        w_link = w1;
                        break;
                    }
                }
            }
            if (partner == UNDEFINED_NODE) continue;

            // Merge: splice the shorter list onto the longer. Nothing else is
            // touched; entries elsewhere go stale and are redirected when
            // their owner is next tested.
            NodeID big = v, small = partner;
            if (adj[big].size() < adj[small].size()) std::swap(big, small);
            adj[big].insert(adj[big].end(), adj[small].begin(),
                            adj[small].end());
            adj[small].clear();
            adj[small].shrink_to_fit();
            deg[big] = deg[v] + deg[partner] - 2 * w_link;
            parent[small] = big;
            ++merges;
            if (deg[big] > 0 && deg[big] < *bound) *bound = deg[big];
            if (!queued[big]) { queued[big] = 1; work.push_back(big); }
        }

        if (merges > 0) {
            for (NodeID v = 0; v < n; ++v) {
                if (find(v) != v) out->Union(v, find(v));
            }
        }
        return merges;
    }
};
