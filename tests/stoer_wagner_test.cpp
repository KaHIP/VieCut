/******************************************************************************
 * stoer_wagner_test.cpp
 *
 * Source of VieCut.
 *
 ******************************************************************************
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#include <cstdint>
#include <random>
#include <string>
#include <tuple>
#include <vector>

#include "algorithms/global_mincut/stoer_wagner_minimum_cut.h"
#include "common/definitions.h"
#include "data_structure/adjlist_graph.h"
#include "data_structure/graph_access.h"
#include "gtest/gtest.h"
#include "io/graph_io.h"
#include "tools/random_functions.h"

typedef std::vector<std::tuple<NodeID, NodeID, EdgeWeight> > EdgeList;

static graph_access make_graph(NodeID n, const EdgeList& edges) {
    std::vector<std::vector<std::pair<NodeID, EdgeWeight> > > adj(n);
    for (auto [u, v, w] : edges) {
        adj[u].emplace_back(v, w);
        adj[v].emplace_back(u, w);
    }
    graph_access G;
    G.start_construction(n, 2 * edges.size());
    for (NodeID u = 0; u < n; ++u) {
        G.new_node();
        for (auto [v, w] : adj[u]) {
            G.new_edge(u, v, w);
        }
    }
    G.finish_construction();
    return G;
}

static uint64_t stoer_wagner(const graph_access& G) {
    adjlist_graph A;
    A.construct_from_graph(G);
    stoer_wagner_minimum_cut sw;
    return sw.perform_minimum_cut(A);
}

// Minimum over all bipartitions; independent of every solver in VieCut.
static uint64_t brute_force(NodeID n, const EdgeList& edges) {
    uint64_t best = UINT64_MAX;
    for (uint64_t side = 1; side < (uint64_t{1} << (n - 1)); ++side) {
        uint64_t cut = 0;
        for (auto [u, v, w] : edges) {
            if (((side >> u) & 1) != ((side >> v) & 1)) {
                cut += w;
            }
        }
        best = std::min(best, cut);
    }
    return best;
}

TEST(StoerWagnerTest, GraphsFromFile) {
    for (auto [file, cut] : { std::pair<std::string, uint64_t>{"small", 2},
                              {"small-wgt", 3} }) {
        auto G = graph_io::readGraphWeighted<graph_access>(
            std::string(VIECUT_PATH) + "/graphs/" + file + ".metis");
        ASSERT_EQ(stoer_wagner(*G), cut) << file;
    }
}

TEST(StoerWagnerTest, TwoVertices) {
    EdgeList edges = { { 0, 1, 7 } };
    ASSERT_EQ(stoer_wagner(make_graph(2, edges)), 7);
}

TEST(StoerWagnerTest, Disconnected) {
    EdgeList two_components = { { 0, 1, 3 }, { 2, 3, 2 } };
    ASSERT_EQ(stoer_wagner(make_graph(4, two_components)), 0);
    // Only the last vertex added in the phase is unreachable.
    EdgeList isolated_vertex = { { 0, 1, 1 }, { 1, 2, 1 }, { 0, 2, 1 } };
    for (int seed = 0; seed < 20; ++seed) {
        random_functions::setSeed(seed);
        ASSERT_EQ(stoer_wagner(make_graph(4, isolated_vertex)), 0);
    }
}

// Random weighted graphs, sparse and dense, light and heavy weights. Heavy
// weights make the maximum degree exceed n, which selects the binary heap
// instead of the bucket queue.
TEST(StoerWagnerTest, RandomGraphsAgainstBruteForce) {
    std::mt19937 rng(42);
    for (int round = 0; round < 500; ++round) {
        NodeID n = 2 + rng() % 11;
        EdgeWeight max_w = (round % 2) ? 1000 : 3;
        EdgeList edges;
        for (NodeID u = 0; u < n; ++u) {
            for (NodeID v = u + 1; v < n; ++v) {
                if (rng() % 3 == 0) {
                    edges.emplace_back(u, v, 1 + rng() % max_w);
                }
            }
        }
        random_functions::setSeed(round);
        ASSERT_EQ(stoer_wagner(make_graph(n, edges)), brute_force(n, edges))
            << "round " << round << ", n = " << n;
    }
}
