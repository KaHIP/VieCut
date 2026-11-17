/******************************************************************************
 * mutable_graph.h
 *
 * Source of VieCut.
 *
 ******************************************************************************
 * Copyright (C) 2018-2019 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#pragma once

#include <algorithm>
#include <iostream>
#include <memory>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "common/definitions.h"
#include "data_structure/graph_access.h"
#include "tlx/logger.hpp"

struct RevEdge {
    NodeID     target;
    EdgeWeight weight;
    EdgeID     reverse_edge;
#ifndef COMPACT_EDGES
    // Only the flow algorithms read these. The minimum cut solver scans this
    // record once per edge visit on its hottest loop, so builds that do not
    // need flow leave them out and the record shrinks from 40 to 24 bytes.
    FlowType   flow;
    size_t     problem_id;
#endif

    RevEdge() { }

    explicit RevEdge(NodeID p_target)
        : target(p_target), weight(1)
#ifndef COMPACT_EDGES
          , flow(0), problem_id(0)
#endif
    { }

    RevEdge(NodeID p_target, EdgeWeight p_wgt)
        : target(p_target), weight(p_wgt)
#ifndef COMPACT_EDGES
          , flow(0), problem_id(0)
#endif
    { }

    RevEdge(NodeID p_target, EdgeWeight p_wgt, EdgeID p_rev)
        : target(p_target),
          weight(p_wgt),
          reverse_edge(p_rev)
#ifndef COMPACT_EDGES
          , flow(0), problem_id(0)
#endif
    { }

    RevEdge(NodeID p_target, EdgeWeight p_wgt, EdgeID p_rev, FlowType p_flow)
        : target(p_target),
          weight(p_wgt),
          reverse_edge(p_rev)
#ifdef COMPACT_EDGES
    { (void)p_flow; }
#else
          , flow(p_flow),
          problem_id(0) { }
#endif
};

class mutable_graph {
 public:
    static constexpr bool debug = false;

    virtual ~mutable_graph() { }

    mutable_graph(const mutable_graph& copy)
        : partition_index(copy.partition_index),
          vertices(copy.vertices),
          weighted_degree(copy.weighted_degree),
          node_in_cut(copy.node_in_cut),
          current_position(copy.current_position),
          contained_in_this(copy.contained_in_this),
          last_node(copy.last_node),
          num_edges(copy.num_edges),
          partition_count(copy.partition_count),
          original_nodes(copy.original_nodes) { }

    mutable_graph() : last_node(0), num_edges(0), partition_count(0) { }

    void start_construction(NodeID n, EdgeID = 0) {
        vertices.resize(n);
        weighted_degree.resize(n, 0);
        partition_index.resize(n, 0);
        node_in_cut.resize(n, 0);
        current_position.resize(n);
        original_nodes = n;

        for (size_t i = 0; i < n; ++i) {
            contained_in_this.emplace_back(1, i);
            current_position[i] = i;
        }
    }

    void resizePositions(NodeID n) {
        current_position.resize(n, 0);
        original_nodes = n;
    }

    NodeID new_node() {
        if (last_node == vertices.size()) {
            contained_in_this.emplace_back(last_node);
            current_position.emplace_back();
            current_position[last_node] = vertices.size();
            vertices.emplace_back();
            partition_index.emplace_back();
            node_in_cut.emplace_back();
            weighted_degree.emplace_back();
        }

        return last_node++;
    }

    NodeID getLastNode() {
        return last_node;
    }

    NodeID new_empty_node() {
        if (last_node == vertices.size()) {
            contained_in_this.emplace_back();
            vertices.emplace_back();
            partition_index.emplace_back();
            node_in_cut.emplace_back();
            weighted_degree.emplace_back();
        }

        return last_node++;
    }

    // Give a vertex room for its edges before they are added. Every edge
    // emplaces into the adjacency of both endpoints, so without this each of
    // those vectors grows geometrically and the records are copied again at
    // every step; on a graph with hundreds of millions of edges that copying
    // is most of the time it takes to build it.
    void reserveNodeDegree(NodeID node, size_t degree) {
        vertices[node].reserve(degree);
    }

    EdgeID new_edge(NodeID source, NodeID target) {
        return new_edge(source, target, 1);
    }

    EdgeID new_edge(NodeID source, NodeID target, EdgeWeight wgt) {
        if (source < target) {
            EdgeID tgt_id = vertices[target].size();
            EdgeID src_id = vertices[source].size();
            vertices[source].emplace_back(target, wgt, tgt_id);
            vertices[target].emplace_back(source, wgt, src_id);
            num_edges += 2;
            weighted_degree[source] += wgt;
            weighted_degree[target] += wgt;
            return src_id;
        } else {
            return -1;
        }
    }

    void new_edge_order(NodeID source, NodeID target, EdgeWeight wgt) {
        if (source < target) {
            new_edge(source, target, wgt);
        } else {
            if (source > target) {
                new_edge(target, source, wgt);
            }
        }
    }

    void computeDegrees() {
        // do nothing, this exists for compatibility with graph_access in
        // templates
    }

    auto nodes() const {
        return iterator<NodeID>(0, number_of_nodes());
    }

    auto edges_of(NodeID n) const {
        return iterator<EdgeID>(get_first_edge(n), get_first_invalid_edge(n));
    }

    // The adjacency of a vertex as a contiguous range. Reading it through
    // getEdge(n, e) makes the compiler reload the row for every edge, since it
    // cannot prove that the queue updates in the loop leave it alone.
    const std::vector<RevEdge>& adjacencyOf(NodeID n) const {
        return vertices[n];
    }

    auto edges_of_starting_at(NodeID n, EdgeID e) const {
        return iterator<EdgeID>(e, get_first_invalid_edge(n));
    }

    bool isEmpty(NodeID n) const {
        return (contained_in_this[n].size() == 0);
    }

    NodeID getEdgeSource(NodeID node, EdgeID edge) const {
        RevEdge e = vertices[node][edge];

        return vertices[e.target][e.reverse_edge].target;
    }

    void setEdgeFlow(NodeID node, EdgeID edge, FlowType flow, size_t f_prob) {
#ifdef COMPACT_EDGES
        (void)node; (void)edge; (void)flow; (void)f_prob;
#else
        vertices[node][edge].flow = flow;
        vertices[node][edge].problem_id = f_prob;
#endif
    }

    void addEdgeFlow(NodeID node, EdgeID edge, FlowType addFlow, size_t f_pro) {
#ifdef COMPACT_EDGES
        (void)node; (void)edge; (void)addFlow; (void)f_pro;
#else
        if (vertices[node][edge].problem_id == f_pro) {
            vertices[node][edge].flow += addFlow;
        } else {
            vertices[node][edge].flow = addFlow;
            vertices[node][edge].problem_id = f_pro;
        }
#endif
    }

    FlowType getEdgeFlow(NodeID node, EdgeID edge, size_t f_problem) {
#ifdef COMPACT_EDGES
        (void)node; (void)edge; (void)f_problem;
        return 0;
#else
        if (vertices[node][edge].problem_id == f_problem) [[unlikely]] {
            return vertices[node][edge].flow;
        } else {
            vertices[node][edge].flow = 0;
            vertices[node][edge].problem_id = f_problem;
            return 0;
        }
#endif
    }

    FlowType getEdgeFlow(NodeID node, EdgeID edge) {
#ifdef COMPACT_EDGES
        (void)node; (void)edge;
        return 0;
#else
        return vertices[node][edge].flow;
#endif
    }

    void finish_construction() {
        last_node = vertices.size();
    }

    // Install one vertex's neighbour list directly, and set the edge count
    // once at the end.
    //
    // The parallel contraction computes every contracted vertex's neighbours
    // independently and can therefore fill them concurrently, which new_edge
    // cannot: it writes both endpoints of an edge, so two threads working on
    // different vertices still collide.
    void setAdjacency(NodeID node, std::vector<RevEdge>&& adj) {
        EdgeWeight wd = 0;
        for (const auto& e : adj) {
            wd += e.weight;
        }
        weighted_degree[node] = wd;
        vertices[node] = std::move(adj);
    }

    void setEdgeCount(EdgeID m) {
        num_edges = m;
    }

    /* ============================================================= */
    /* graph access methods */
    /* ============================================================= */
    NodeID number_of_nodes() const {
        return vertices.size();
    }

    NodeID n() const {
        return number_of_nodes();
    }

    EdgeID number_of_edges() const {
        return num_edges;
    }

    EdgeID m() const {
        return number_of_edges();
    }

    EdgeID get_first_edge(NodeID) const {
        return 0;
    }

    EdgeID get_first_invalid_edge(NodeID node) const {
        return vertices[node].size();
    }

    PartitionID get_partition_count() const {
        return partition_count;
    }

    void set_partition_count(PartitionID count) {
        partition_count = count;
    }

    bool getNodeInCut(NodeID n) {
        return node_in_cut[n];
    }

    void setNodeInCut(NodeID n, bool b) {
        node_in_cut[n] = b;
    }

    PartitionID getPartitionIndex(NodeID node) {
        return partition_index[node];
    }

    void setPartitionIndex(NodeID node, PartitionID id) {
        partition_index[node] = id;
    }

    EdgeWeight getUnweightedNodeDegree(NodeID node) const {
        return get_first_invalid_edge(node);
    }

    EdgeWeight getWeightedNodeDegree(NodeID node) const {
        return weighted_degree[node];
    }

    EdgeWeight getMaxDegree() {
        return *std::max_element(weighted_degree.begin(),
                                 weighted_degree.end());
    }

    EdgeID getMaxUnweightedDegree() {
        EdgeID max_size = 0;
        for (NodeID n : nodes()) {
            EdgeID curr_size = vertices[n].size();
            max_size = std::max(max_size, curr_size);
        }
        return max_size;
    }

    NodeID getOriginalNodes() {
        return original_nodes;
    }

    EdgeWeight sumOfEdgeWeights() {
        EdgeWeight sum_weights = 0;
        for (NodeID n : nodes()) {
            sum_weights += weighted_degree[n];
        }
        return sum_weights;
    }

    void setOriginalNodes(NodeID n) {
        original_nodes = n;
        if (n > current_position.size()) {
            current_position.resize(n);
        }
    }

    EdgeWeight getMinDegree() {
        return *std::min_element(weighted_degree.begin(),
                                 weighted_degree.end());
    }

    EdgeWeight getEdgeWeight(NodeID node, EdgeID edge) const {
        return vertices[node][edge].weight;
    }

    size_t numContainedVertices(NodeID node) const {
        return contained_in_this[node].size();
    }

    std::vector<NodeID> containedVertices(NodeID node) const {
        return contained_in_this[node];
    }

    // containedVertices returns by value, so reading a single element with
    // containedVertices(n)[0] allocates and copies the whole vector. The
    // contraction does that once per vertex per round, which is quadratic over
    // a run on the families that contract one pair at a time.
    NodeID firstContainedVertex(NodeID node) const {
        return contained_in_this[node][0];
    }

    void setContainedVertices(NodeID node, std::vector<NodeID> v) {
        contained_in_this[node] = v;
    }

    void addContainedVertex(NodeID node, NodeID v) {
        contained_in_this[node].emplace_back(v);
    }

    void setEdgeWeight(NodeID node, EdgeID edge, EdgeWeight weight) {
        RevEdge e = vertices[node][edge];

        EdgeWeight curr_weight = e.weight;

        NodeID target = e.target;
        EdgeID rev = e.reverse_edge;

        vertices[node][edge].weight = weight;
        vertices[target][rev].weight = weight;

        weighted_degree[node] = weighted_degree[node] + weight - curr_weight;
        weighted_degree[target] =
            weighted_degree[target] + weight - curr_weight;
    }

    std::pair<NodeID, EdgeWeight> getEdge(NodeID node, EdgeID edge) const {
        return { getEdgeTarget(node, edge), getEdgeWeight(node, edge) };
    }

    NodeID getEdgeTarget(NodeID node, EdgeID edge) const {
        return vertices[node][edge].target;
    }

    EdgeID getReverseEdge(NodeID node, EdgeID edge) {
        return vertices[node][edge].reverse_edge;
    }

    NodeID getCurrentPosition(NodeID node) {
        return current_position[node];
    }

    void setCurrentPosition(NodeID node, NodeID pos) {
        current_position[node] = pos;
    }

    void deleteVertex(NodeID node) {
        last_node--;
        NodeID back = number_of_nodes() - 1;

        if (node < back) {
            for (NodeID n : contained_in_this[back]) {
                current_position[n] = node;
            }

            for (NodeID n : contained_in_this[node]) {
                current_position[n] = back;
            }

            contained_in_this[node] = contained_in_this[back];
        }

        for (EdgeID e : edges_of(node)) {
            NodeID tgt = getEdgeTarget(node, e);
            EdgeID rev = getReverseEdge(node, e);
            EdgeWeight wgt = getEdgeWeight(node, e);
            internalDeleteEdge(tgt, rev);
            num_edges -= 2;
            weighted_degree[tgt] -= wgt;
        }

        if (node < back) {
            for (EdgeID e : edges_of(back)) {
                NodeID tgt = getEdgeTarget(back, e);
                EdgeID rev = getReverseEdge(back, e);
                vertices[tgt][rev].target = node;
            }
        }

        vertices[node] = std::move(vertices.back());
        weighted_degree[node] = std::move(weighted_degree.back());
        partition_index[node] = std::move(partition_index.back());
        node_in_cut[node] = std::move(node_in_cut.back());
        vertices.pop_back();
        weighted_degree.pop_back();
        partition_index.pop_back();
        node_in_cut.pop_back();
    }

    // Contraction methods
    void deleteEdge(NodeID node, EdgeID edge) {
        RevEdge e = vertices[node][edge];
        NodeID target = e.target;
        EdgeWeight weight = e.weight;

        if (!vertices[node].size() || !vertices[target].size()) {
            LOG1 << "ERROR: Deleting edge from empty vertex! Exiting!";
            exit(3);
        }

        if (get_first_invalid_edge(node) > edge + 1) {
            vertices[node][edge] =
                std::move(vertices[node][vertices[node].size() - 1]);
            vertices[vertices[node][edge].target]
            [vertices[node][edge].reverse_edge].reverse_edge = edge;
        }
        vertices[node].pop_back();

        if (get_first_invalid_edge(target) > e.reverse_edge + 1) {
            vertices[target][e.reverse_edge] =
                std::move(vertices[target][vertices[target].size() - 1]);
            vertices[vertices[target][e.reverse_edge].target]
            [vertices[target][e.reverse_edge].reverse_edge].reverse_edge
                = e.reverse_edge;
        }
        vertices[target].pop_back();

        weighted_degree[node] -= weight;
        weighted_degree[target] -= weight;

        num_edges -= 2;
    }

    void mergeEdgeSparse(NodeID node, NodeID target, EdgeID ed) {
        NodeID e_target = getEdgeTarget(target, ed);
        EdgeWeight e_weight = getEdgeWeight(target, ed);
        EdgeID del_rev = getReverseEdge(target, ed);
        bool edge_found = false;
        for (EdgeID src_edge : edges_of(node)) {
            if (getEdgeTarget(node, src_edge) == e_target) {
                EdgeID rev = getReverseEdge(node, src_edge);
                vertices[node][src_edge].weight += e_weight;
                vertices[e_target][rev].weight += e_weight;
                weighted_degree[node] += e_weight;
                num_edges -= 2;
                internalDeleteEdge(e_target, del_rev);
                edge_found = true;
                break;
            }
        }
        if (!edge_found) {
            // create new edge and map reverse edge
            vertices[e_target][del_rev].reverse_edge
                = vertices[node].size();
            vertices[e_target][del_rev].target = node;
            vertices[node].emplace_back(
                e_target, e_weight, del_rev);
            weighted_degree[node] += e_weight;
        }
    }

    NodeID contractSparseTargetNoEdge(NodeID node, NodeID target) {
        if (node == target + 1) {
            std::unordered_set<NodeID> nodes = { node, target };
            contractVertexSet(nodes);
            return target;
        }

        last_node--;
        for (EdgeID ed : edges_of(target)) {
            mergeEdgeSparse(node, target, ed);
        }

        for (NodeID n : contained_in_this[target]) {
            contained_in_this[node].emplace_back(n);
            current_position[n] = node;
        }

        contained_in_this[target].clear();

        for (NodeID n : contained_in_this[vertices.size() - 1]) {
            contained_in_this[target].emplace_back(n);
            current_position[n] = target;
        }

        vertices[target] = std::move(vertices.back());
        weighted_degree[target] = std::move(weighted_degree.back());
        partition_index[target] = std::move(partition_index.back());
        node_in_cut[target] = std::move(node_in_cut.back());
        vertices.pop_back();
        weighted_degree.pop_back();
        partition_index.pop_back();
        node_in_cut.pop_back();
        contained_in_this.pop_back();

        // remap all reverse edges of vertex that was now moved to 'target'
        for (EdgeID ed : edges_of(target)) {
            RevEdge e = vertices[target][ed];
            vertices[e.target][e.reverse_edge].target = target;
        }

        return target;
    }

    NodeID contractEdgeSparseTarget(NodeID node, EdgeID edge) {
        RevEdge e = vertices[node][edge];
        NodeID target = e.target;

        if (node == target + 1) {
            return contractEdge(node, edge);
        }
        last_node--;
        for (EdgeID ed : edges_of(target)) {
            if (ed != e.reverse_edge) {
                mergeEdgeSparse(node, target, ed);
            }
        }

        for (NodeID n : contained_in_this[target]) {
            contained_in_this[node].emplace_back(n);
            current_position[n] = node;
        }

        contained_in_this[target].clear();

        for (NodeID n : contained_in_this[vertices.size() - 1]) {
            contained_in_this[target].emplace_back(n);
            current_position[n] = target;
        }

        // delete edge of 'node' to 'target', remap reverse
        internalDeleteEdge(node, edge);
        num_edges -= 2;
        weighted_degree[node] -= e.weight;

        vertices[target] = std::move(vertices.back());
        weighted_degree[target] = std::move(weighted_degree.back());
        partition_index[target] = std::move(partition_index.back());
        node_in_cut[target] = std::move(node_in_cut.back());
        vertices.pop_back();
        weighted_degree.pop_back();
        partition_index.pop_back();
        node_in_cut.pop_back();
        contained_in_this.pop_back();

        // remap all reverse edges of vertex that was now moved to 'target'
        for (EdgeID ed : edges_of(target)) {
            RevEdge e = vertices[target][ed];
            vertices[e.target][e.reverse_edge].target = target;
        }

        return target;
    }

    NodeID contractEdge(NodeID node, EdgeID edge) {
        last_node--;

        if (getEdgeTarget(node, edge) < node) {
            EdgeID re = getReverseEdge(node, edge);
            node = getEdgeTarget(node, edge);
            edge = re;
        }

        RevEdge e = vertices[node][edge];
        NodeID target = e.target;
        EdgeID del_id = 0;
        EdgeWeight del_wgt = 0;

        // The shared neighbourhood used to be found through an
        // unordered_map built per contraction. This runs once per merge on
        // the hottest contraction path, and hashing was most of its constant:
        // on the dense families every merge paid a map build over a
        // neighbourhood of thousands. A stamp array indexed by vertex does
        // the same lookup at the cost of an array read; the epoch makes
        // clearing free, and the arrays persist with the graph.
        if (ctr_stamp_.size() < vertices.size()) {
            ctr_stamp_.assign(vertices.size(), 0);
            ctr_edge_.resize(vertices.size());
        }
        if (++ctr_epoch_ == 0) {          // wrapped: reset lazily
            std::fill(ctr_stamp_.begin(), ctr_stamp_.end(), 0);
            ctr_epoch_ = 1;
        }
        for (EdgeID ed : edges_of(node)) {
            if (ed != edge) {
                RevEdge del_ed = vertices[node][ed];
                const NodeID t = getEdgeTarget(node, ed);
                ctr_stamp_[t] = ctr_epoch_;
                ctr_edge_[t] = std::make_tuple(del_ed.target,
                                               del_ed.reverse_edge, ed);
            } else {
                del_wgt = vertices[node][ed].weight;
                del_id = ed;
            }
        }

        for (EdgeID ed : edges_of(target)) {
            if (ed != e.reverse_edge) {
                RevEdge del_edge = vertices[target][ed];
                NodeID e_target = getEdgeTarget(target, ed);
                // neighbour of 'target' also in neighbourhood of 'node'
                if (ctr_stamp_[e_target] == ctr_epoch_) {
                    // sum up edge weights
                    auto [tgt, rev, ed2] = ctr_edge_[e_target];
                    vertices[node][ed2].weight += del_edge.weight;
                    vertices[tgt][rev].weight += del_edge.weight;

                    weighted_degree[node] += del_edge.weight;

                    num_edges -= 2;

                    internalDeleteEdge(del_edge.target, del_edge.reverse_edge);
                } else {
                    // create new edge and map reverse edge
                    vertices[e_target][del_edge.reverse_edge].reverse_edge
                        = vertices[node].size();
                    vertices[e_target][del_edge.reverse_edge].target = node;
                    vertices[node].emplace_back(e_target, del_edge.weight,
                                                del_edge.reverse_edge);
                    weighted_degree[node] += del_edge.weight;
                }
            }
        }

        for (NodeID n : contained_in_this[target]) {
            contained_in_this[node].emplace_back(n);
            current_position[n] = node;
        }

        contained_in_this[target].clear();

        for (NodeID n : contained_in_this[vertices.size() - 1]) {
            contained_in_this[target].emplace_back(n);
            current_position[n] = target;
        }

        vertices[target] = std::move(vertices.back());
        weighted_degree[target] = std::move(weighted_degree.back());
        partition_index[target] = std::move(partition_index.back());
        node_in_cut[target] = std::move(node_in_cut.back());
        vertices.pop_back();
        weighted_degree.pop_back();
        partition_index.pop_back();
        node_in_cut.pop_back();
        contained_in_this.pop_back();

        // remap all reverse edges of vertex that was now moved to 'target'
        for (EdgeID ed : edges_of(target)) {
            RevEdge e = vertices[target][ed];
            vertices[e.target][e.reverse_edge].target = target;
        }

        // delete edge of 'node' to 'target', remap reverse
        internalDeleteEdge(node, del_id);
        num_edges -= 2;
        weighted_degree[node] -= del_wgt;

        return target;
    }

    void contractVertexSet(const std::unordered_set<NodeID>& vertex_set) {
        std::vector<NodeID> vertex_set_vec;
        NodeID first = UNDEFINED_NODE;
        size_t idx = 0;

        for (auto i : vertex_set) {
            if (i >= vertices.size())
                continue;

            if (i < first) {
                first = i;
                idx = vertex_set_vec.size();
            }
            vertex_set_vec.emplace_back(i);
        }

        if (vertex_set_vec.size() <= 1)
            return;

        vertex_set_vec.erase(vertex_set_vec.begin() + idx);

        std::vector<RevEdge> edges;
        // Membership and the shared-neighbourhood positions both go through
        // stamp arrays rather than the hash containers this used: those cost
        // a hash per edge of every member, once per merge, on the hottest
        // contraction path, and on the dense families that was most of the
        // running time of a round.
        if (mem_stamp_.size() < vertices.size()) {
            mem_stamp_.assign(vertices.size(), 0);
        }
        if (++mem_epoch_ == 0) {
            std::fill(mem_stamp_.begin(), mem_stamp_.end(), 0);
            mem_epoch_ = 1;
        }
        for (NodeID i : vertex_set) {
            if (i < vertices.size()) mem_stamp_[i] = mem_epoch_;
        }
        if (ctr_stamp_.size() < vertices.size()) {
            ctr_stamp_.assign(vertices.size(), 0);
            ctr_edge_.resize(vertices.size());
        }
        if (++ctr_epoch_ == 0) {
            std::fill(ctr_stamp_.begin(), ctr_stamp_.end(), 0);
            ctr_epoch_ = 1;
        }

        for (EdgeID e : edges_of(first)) {
            NodeID target = getEdgeTarget(first, e);
            if (mem_stamp_[target] != mem_epoch_) {
                RevEdge edge_to_set = vertices[first][e];
                edges.push_back(vertices[first][e]);
                ctr_stamp_[target] = ctr_epoch_;
                ctr_edge_[target] = std::make_tuple(
                    static_cast<EdgeWeight>(edge_to_set.target),
                    edge_to_set.reverse_edge,
                    static_cast<EdgeID>(edges.size() - 1));
                vertices[edge_to_set.target]
                [edge_to_set.reverse_edge].reverse_edge = edges.size() - 1;
            } else {
                weighted_degree[first] -= vertices[first][e].weight;
                num_edges--;
            }
        }

        vertices[first].swap(edges);

        for (NodeID n : vertex_set_vec) {
            if (n == first)
                continue;

            for (EdgeID e : edges_of(n)) {
                NodeID target = getEdgeTarget(n, e);
                if (mem_stamp_[target] != mem_epoch_) {
                    RevEdge del_edge = vertices[n][e];
                    // neighbour of 'target' also in neighbourhood of 'node'
                    if (ctr_stamp_[target] == ctr_epoch_) {
                        // sum up edge weights
                        auto[tgt, rev, ed2] = ctr_edge_[target];

                        VIECUT_ASSERT_EQ(vertices[first][ed2].weight,
                                         vertices[tgt][rev].weight);
                        vertices[first][ed2].weight += del_edge.weight;
                        vertices[tgt][rev].weight += del_edge.weight;

                        weighted_degree[first] += del_edge.weight;
                        num_edges -= 2;

                        if (vertices[tgt].size() - 1 == rev) {
                            // we need to move position of reverse edge in 'map'
                            ctr_edge_[target] = std::make_tuple(
                                static_cast<EdgeWeight>(tgt),
                                del_edge.reverse_edge, ed2);
                        }

                        internalDeleteEdge(del_edge.target,
                                           del_edge.reverse_edge);
                    } else {
                        VIECUT_ASSERT_EQ(
                            del_edge.weight,
                            vertices[target][del_edge.reverse_edge].weight);

                        // create new edge and map reverse edge
                        vertices[target][del_edge.reverse_edge].reverse_edge =
                            vertices[first].size();
                        vertices[target][del_edge.reverse_edge].target = first;

                        vertices[first].emplace_back(target, del_edge.weight,
                                                     del_edge.reverse_edge);
                        weighted_degree[first] += del_edge.weight;
                        ctr_stamp_[target] = ctr_epoch_;
                        ctr_edge_[target] = std::make_tuple(
                            static_cast<EdgeWeight>(target),
                            del_edge.reverse_edge,
                            static_cast<EdgeID>(vertices[first].size() - 1));
                    }
                } else {
                    --num_edges;
                }
            }
        }

        std::sort(vertex_set_vec.begin(), vertex_set_vec.end(),
                  [](auto v1, auto v2) {
                      return v1 > v2;
                  });

        for (NodeID vtx : vertex_set_vec) {
            for (NodeID n : contained_in_this[vtx]) {
                contained_in_this[first].emplace_back(n);
                current_position[n] = first;
            }
            contained_in_this[vtx].clear();

            if (vtx < vertices.size() - 1) {
                for (NodeID n : contained_in_this[vertices.size() - 1]) {
                    contained_in_this[vtx].emplace_back(n);
                    current_position[n] = vtx;
                }
                vertices[vtx] = std::move(vertices.back());
                weighted_degree[vtx] = std::move(weighted_degree.back());
                partition_index[vtx] = std::move(partition_index.back());
                node_in_cut[vtx] = std::move(node_in_cut.back());

                // remap all rev edges of vertex that was now moved to 'target'
                for (EdgeID ed : edges_of(vtx)) {
                    RevEdge e = vertices[vtx][ed];
                    VIECUT_ASSERT_EQ(vertices[e.target][e.reverse_edge].target,
                                     vertices.size() - 1);
                    vertices[e.target][e.reverse_edge].target = vtx;
                }
            }

            vertices.pop_back();
            weighted_degree.pop_back();
            partition_index.pop_back();
            node_in_cut.pop_back();
            contained_in_this.pop_back();
            last_node--;
        }
    }

    // Graph class translation
    static mutableGraphPtr from_graph_access(
        graphAccessPtr G) {
        mutableGraphPtr m_G = std::make_shared<mutable_graph>();
        std::vector<NodeID> last_incident(G->number_of_nodes(), UNDEFINED_NODE);
        m_G->start_construction(G->number_of_nodes());
        for (NodeID n : G->nodes()) {
            m_G->new_node();
            for (EdgeID e : G->edges_of(n)) {
                NodeID t = G->getEdgeTarget(e);
                if (last_incident[t] == n) {
                    // double edge
                    // LOG1 << "double edge";

                    for (EdgeID e2 : m_G->edges_of(n)) {
                        if (m_G->getEdgeTarget(n, e2) == t) {
                            m_G->setEdgeWeight(n, e2,
                                               m_G->getEdgeWeight(n, e2)
                                               + G->getEdgeWeight(e));
                            break;
                        }
                    }
                } else {
                    last_incident[t] = n;
                    m_G->new_edge(n, G->getEdgeTarget(e), G->getEdgeWeight(e));
                }
            }
        }
        m_G->finish_construction();
        return m_G;
    }

    graphAccessPtr to_graph_access() {
        graphAccessPtr G = std::make_shared<graph_access>();
        G->start_construction(number_of_nodes(), number_of_edges());
        for (NodeID n : G->nodes()) {
            G->new_node();
            for (EdgeID e : edges_of(n)) {
                EdgeID edge = G->new_edge(n, getEdgeTarget(n, e));
                G->setEdgeWeight(edge, getEdgeWeight(n, e));
            }
        }
        G->finish_construction();
        return G;
    }

    static mutableGraphPtr deserialize(std::vector<uint64_t> v) {
        mutableGraphPtr G = std::make_shared<mutable_graph>();
        uint64_t num_nodes = v[0];
        size_t num_edges = v[1];
        uint64_t original_nodes = v[3];
        G->start_construction(num_nodes);
        G->set_partition_count(v[2]);
        G->setOriginalNodes(original_nodes);
        size_t n = 0;
        size_t deserial = 4;
        while (n < num_nodes) {
            uint64_t next = v[deserial++];
            if (next >= num_edges) {
                if (next == UNDEFINED_EDGE) {
                    break;
                }
                n = next - num_edges;
                continue;
            }
            uint64_t wgt = v[deserial++];
            G->new_edge(n, next, wgt);
        }

        for (size_t i = 0; i < num_nodes; ++i) {
            G->setPartitionIndex(i, v[deserial++]);
            G->setContainedVertices(i, { });
        }

        for (size_t i = 0; i < original_nodes; ++i) {
            NodeID pos = static_cast<NodeID>(v[deserial++]);
            G->setCurrentPosition(i, pos);
            G->addContainedVertex(pos, i);
        }

        G->finish_construction();

        return G;
    }

    std::vector<uint64_t> serialize() {
        std::vector<uint64_t> serial(5 + 2 * n() + m() + original_nodes);
        serial[0] = static_cast<uint64_t>(vertices.size());
        serial[1] = static_cast<uint64_t>(num_edges);
        serial[2] = static_cast<uint64_t>(partition_count);
        serial[3] = static_cast<uint64_t>(original_nodes);
        size_t next = 4;

        for (NodeID n : nodes()) {
            serial[next++] = static_cast<uint64_t>(num_edges + n);
#ifdef COMPACT_EDGES
            for (const auto& [t, w, r] : vertices[n]) {
#else
            for (const auto& [t, w, r, f, l] : vertices[n]) {
#endif
                // I am deeply sorry for this ugly code, but structured bindings
                // seem to not work in combination with maybe_unused to suppress
                // unintended unused warnings.
                (void)r;
#ifndef COMPACT_EDGES
                (void)f;
                (void)l;
#endif
                if (t > n) {
                    serial[next++] = static_cast<uint64_t>(t);
                    serial[next++] = static_cast<uint64_t>(w);
                }
            }
        }

        serial[next++] = UNDEFINED_EDGE;

        for (const auto& p : partition_index) {
            serial[next++] = static_cast<uint64_t>(p);
        }

        for (const auto& c : current_position) {
            serial[next++] = static_cast<uint64_t>(c);
        }

        return serial;
    }

    mutableGraphPtr simplify() {
        mutableGraphPtr G = std::make_shared<mutable_graph>();
        G->start_construction(number_of_nodes());
        for (NodeID n : G->nodes()) {
            G->new_node();
            for (EdgeID e : edges_of(n)) {
                G->new_edge(n, getEdgeTarget(n, e), getEdgeWeight(n, e));
            }
            G->setPartitionIndex(n, getPartitionIndex(n));
        }
        G->finish_construction();
        return G;
    }

    void resetContainedvertices() {
        original_nodes = vertices.size();
        for (NodeID n : nodes()) {
            current_position[n] = n;
            contained_in_this[n] = { n };
        }
    }

 private:
    void internalDeleteEdge(NodeID n, EdgeID e) {
        if (vertices[n].size() > e + 1) {
            vertices[n][e] = std::move(vertices[n][vertices[n].size() - 1]);
            vertices[n].pop_back();

            NodeID tgt = vertices[n][e].target;
            EdgeID rev = vertices[n][e].reverse_edge;
            VIECUT_ASSERT_EQ(vertices[tgt][rev].reverse_edge,
                             vertices[n].size());
            vertices[tgt][rev].reverse_edge = e;
        } else {
            if (vertices[n].size()) {
                vertices[n].pop_back();
            } else {
                LOG1 << "ERROR: Trying to delete edge from empty vertex!";
                exit(-1);
            }
        }
    }

    std::vector<PartitionID> partition_index;
    std::vector<std::vector<RevEdge> > vertices;
    // scratch for contractEdge's shared-neighbourhood lookup
    std::vector<uint32_t> mem_stamp_;
    uint32_t mem_epoch_ = 0;
    std::vector<uint32_t> ctr_stamp_;
    std::vector<std::tuple<EdgeWeight, EdgeID, EdgeID> > ctr_edge_;
    uint32_t ctr_epoch_ = 0;
    std::vector<EdgeWeight> weighted_degree;
    std::vector<bool> node_in_cut;
    std::vector<NodeID> current_position;
    std::vector<std::vector<NodeID> > contained_in_this;

    NodeID last_node;
    EdgeID num_edges;
    PartitionID partition_count;
    NodeID original_nodes;
};

[[maybe_unused]] static std::string toStringUnweighted(
    mutableGraphPtr G) {
    std::ostringstream oss;

    for (NodeID n : G->nodes()) {
        oss << n << ": [";
        for (EdgeID e : G->edges_of(n)) {
            oss << G->getEdgeTarget(n, e);
            if (e < G->get_first_invalid_edge(n) - 1) {
                oss << ",";
            }
        }
        oss << "] ";
    }
    return oss.str();
}

[[maybe_unused]] static std::string toStringCompact(
    mutableGraphPtr G) {
    std::ostringstream oss;

    for (NodeID n : G->nodes()) {
        oss << " [";
        for (EdgeID e : G->edges_of(n)) {
            oss << G->getEdgeTarget(n, e);
            if (e < G->get_first_invalid_edge(n) - 1) {
                oss << ",";
            }
        }
        oss << "]";
    }
    return oss.str();
}

[[maybe_unused]] static std::string toStringWeighted(
    mutableGraphPtr G) {
    std::ostringstream oss;

    oss << "Graph n " << G->n() << " m " << G->m() << "\n";
    for (NodeID n : G->nodes()) {
        oss << n << ": [";
        for (EdgeID e : G->edges_of(n)) {
            oss << G->getEdgeTarget(n, e) << "("
                << G->getEdgeWeight(n, e) << ")";
            if (e < G->get_first_invalid_edge(n) - 1) {
                oss << ",";
            }
        }
        oss << "]\n";
    }
    return oss.str();
}

[[maybe_unused]] static std::string toStringContains(
    mutableGraphPtr G) {
    std::ostringstream oss;

    for (NodeID n : G->nodes()) {
        oss << n << ": [";
        for (NodeID v : G->containedVertices(n)) {
            oss << v;
            oss << ",";
        }
        oss << "]\n";
    }
    return oss.str();
}

[[maybe_unused]] static std::ostream& operator << (
    std::ostream& os, mutableGraphPtr G) {
    // return os << toStringUnweighted(G);
    return os << toStringWeighted(G);
    // return os << toStringCompact(G);
    // return os << toStringContains(G);
}
