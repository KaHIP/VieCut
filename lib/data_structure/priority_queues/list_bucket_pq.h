/******************************************************************************
 * list_bucket_pq.h
 *
 * A bucket priority queue whose buckets are intrusive linked lists.
 *
 * The bucket queue this replaces holds a std::deque per bucket. The capforest
 * inserts or raises a key once per edge visit, and every one of those is a
 * deque operation that walks the deque's map array and can allocate. With the
 * lists below a bucket is one index, and an insert or a removal is a handful of
 * array writes into arrays that are sized by the graph, not by the cut.
 *
 * Keys are the accumulated weight of a vertex and lie in [0, span].
 *****************************************************************************/

#pragma once

#include <limits>
#include <vector>

#include "common/definitions.h"
#include "data_structure/priority_queues/priority_queue_interface.h"

class list_bucket_pq {
 public:
    static constexpr NodeID kNone = std::numeric_limits<NodeID>::max();

    list_bucket_pq(const NodeID& num_nodes, const EdgeWeight& gain_span)
        : m_elements(0),
          m_max_idx(0),
          m_head(static_cast<size_t>(gain_span) + 2, kNone),
          m_next(num_nodes, kNone),
          m_prev(num_nodes, kNone),
          m_key(num_nodes, kNone) { }

    // The capforest drains the queue every round, so when a round ends every
    // key is already kNone and every bucket head is already kNone: the queue is
    // clean without being cleared. That makes it reusable across rounds, which
    // matters because the arrays are sized by the graph and allocating them per
    // round costs more than the queue operations themselves on the instances
    // whose cut is small.
    void reuse(const NodeID& num_nodes, const EdgeWeight& gain_span) {
        if (m_head.size() < static_cast<size_t>(gain_span) + 2) {
            m_head.assign(static_cast<size_t>(gain_span) + 2, kNone);
        }
        if (m_next.size() < num_nodes) {
            m_next.assign(num_nodes, kNone);
            m_prev.assign(num_nodes, kNone);
            m_key.assign(num_nodes, kNone);
        }
        m_elements = 0;
        m_max_idx = 0;
    }

    NodeID size() const { return m_elements; }
    bool empty() const { return m_elements == 0; }
    bool contains(NodeID node) const { return m_key[node] != kNone; }

    void insert(NodeID node, Gain gain) {
        const NodeID k = static_cast<NodeID>(gain);
        m_key[node] = k;
        link(node, k);
        if (k > m_max_idx) m_max_idx = k;
        ++m_elements;
    }

    void increaseKey(NodeID node, Gain gain) {
        const NodeID k = static_cast<NodeID>(gain);
        unlink(node, m_key[node]);
        m_key[node] = k;
        link(node, k);
        if (k > m_max_idx) m_max_idx = k;
    }

    // One lookup instead of contains() followed by insert or increaseKey.
    void insertOrIncrease(NodeID node, Gain gain) {
        if (m_key[node] == kNone) {
            insert(node, gain);
        } else {
            increaseKey(node, gain);
        }
    }

    NodeID deleteMax() {
        while (m_head[m_max_idx] == kNone && m_max_idx > 0) {
            --m_max_idx;
        }
        const NodeID node = m_head[m_max_idx];
        unlink(node, m_max_idx);
        m_key[node] = kNone;
        --m_elements;
        return node;
    }

 private:
    void link(NodeID node, NodeID bucket) {
        const NodeID h = m_head[bucket];
        m_next[node] = h;
        m_prev[node] = kNone;
        if (h != kNone) m_prev[h] = node;
        m_head[bucket] = node;
    }

    void unlink(NodeID node, NodeID bucket) {
        const NodeID nx = m_next[node], pv = m_prev[node];
        if (pv != kNone) {
            m_next[pv] = nx;
        } else {
            m_head[bucket] = nx;
        }
        if (nx != kNone) m_prev[nx] = pv;
    }

    NodeID m_elements;
    NodeID m_max_idx;
    std::vector<NodeID> m_head;
    std::vector<NodeID> m_next;
    std::vector<NodeID> m_prev;
    std::vector<NodeID> m_key;
};
