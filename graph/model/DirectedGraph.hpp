/**
 * @file Graph.hpp
 * @author bwu
 * @brief basic graph models
 * @version 0.1
 * @date 2024-11-29
 */
#pragma once
#include "generic/utils/Index.hpp"
#include "generic/utils/LinearMap.hpp"
#include <optional>
#include <set>
namespace generic::graph::model {

struct NodeIdTag;
struct EdgeIdTag;
struct LevelIdTag;

using namespace generic::utils;
using NodeId = Index<NodeIdTag>;
using EdgeId = Index<EdgeIdTag>;
using LevelId = Index<LevelIdTag>;

template<typename T>
class Range
{
public:
    Range(T b, T e): m_begin(b), m_end(e) {}

    T Begin() const { return m_begin; }
    T End() const { return m_end; }
    bool Empty() const { return m_begin == m_end; }
    size_t Size() const { return std::distance(m_begin, m_end); }
    std::pair<T, T> toStdPair() const { return {m_begin, m_end}; }

    // std-style interface
    T begin() const { return Begin(); }
    T end() const { return End(); }
private:
    T m_begin, m_end;
};

template<typename T>
inline Range<T> makeRange(T b, T e) { return Range<T>(b, e); }

struct GraphIdMaps
{
    LinearMap<NodeId, NodeId> nodeMap;
    LinearMap<EdgeId, EdgeId> edgeMap;
    GraphIdMaps(LinearMap<NodeId, NodeId> nodeMap, LinearMap<EdgeId, EdgeId> edgeMap)
     : nodeMap(std::move(nodeMap)), edgeMap(std::move(edgeMap))
    {
    }
};

template <typename T>
class Recycler
{
public:
    size_t Size() const { return m_objs.size(); }
    bool Empty() const { return m_objs.empty(); }
    void Add(T t) { m_objs.emplace_back(t); }
    void Clear() { m_objs.clear(); }

    T Take();
    void ResetLast();
    std::optional<T> LastTake() const;
private:
    std::vector<T> m_objs;
    std::optional<T> m_last{std::nullopt};
};

class DirectedGraph
{
public:
    using NodeIterator = typename LinearMap<NodeId, NodeId>::ConstIterator;
    using EdgeIterator = typename LinearMap<EdgeId, EdgeId>::ConstIterator;
    using NodeRange = Range<NodeIterator>;
    using EdgeRange = Range<EdgeIterator>;

    using OutEdgeIterator = typename std::set<EdgeId>::iterator;
    using OutEdgeRange = Range<OutEdgeIterator>;

    using InEdgeIterator = typename std::set<EdgeId>::iterator;
    using InEdgeRange = Range<InEdgeIterator>;

    virtual ~DirectedGraph() = default;

    virtual typename NodeId::SizeType Nodes() const { return DirectedGraph::Nodes() - m_nodeRecycler.Size(); }
    virtual typename EdgeId::SizeType Edges() const { return DirectedGraph::Edges() - m_edgeRecycler.Size(); }

    virtual NodeRange AllNodes() const { return makeRange(m_nodes.Begin(), m_nodes.End()); }
    virtual EdgeRange AllEdges() const { return makeRange(m_edges.Begin(), m_edges.End()); }

    virtual NodeId Source(const EdgeId id) const { return m_source[id]; }
    virtual NodeId Target(const EdgeId id) const { return m_target[id]; }

    virtual bool isValid(const NodeId id) const { return size_t(id) < m_nodes.Size(); }
    virtual bool isValid(const EdgeId id) const { return size_t(id) < m_edges.Size(); }

    virtual InEdgeRange InEdges(const NodeId id) const { return makeRange(m_inEdges[id].begin(), m_inEdges[id].end()); }
    virtual InEdgeRange OutEdges(const NodeId id) const { return makeRange(m_inEdges[id].begin(), m_inEdges[id].end()); }

    virtual EdgeId FindEdge(NodeId source, NodeId target) const;

// mutable
    virtual NodeId AddNode();

    virtual EdgeId AddEdge(NodeId source, NodeId target);

    virtual void RemoveNode(const NodeId nodeId);

    virtual void RemoveEdge(const NodeId source, const NodeId target);

    virtual void RemoveEdge(const EdgeId edgeId);
protected:
    LinearMap<NodeId, NodeId> m_nodes;
    LinearMap<EdgeId, EdgeId> m_edges;
    LinearMap<EdgeId, NodeId> m_source;
    LinearMap<EdgeId, NodeId> m_target;

    LinearMap<NodeId, std::set<EdgeId>> m_inEdges;
    LinearMap<NodeId, std::set<EdgeId>> m_outEdges;

protected:
    Recycler<NodeId> m_nodeRecycler;
    Recycler<EdgeId> m_edgeRecycler;
};

/// Implementation

template <typename T>
T Recycler<T>::Take()
{
    T t = m_objs.back();
    m_objs.pop_back();
    m_last = t;
    return t;
}

template <typename T>
void Recycler<T>::ResetLast()
{
    m_last = std::nullopt;
}

template <typename T>
std::optional<T> Recycler<T>::LastTake() const
{
    return m_last;
}

EdgeId DirectedGraph::FindEdge(NodeId source, NodeId target) const
{
    for (auto edge : OutEdges(source)) {
        if (target == Target(edge))
            return edge;
    }
    return EdgeId::Invalid();
}

NodeId DirectedGraph::AddNode()
{
    NodeId id;
    if (m_nodeRecycler.Empty()) {
        id = NodeId(m_nodes.Size());
        m_nodes.Append(id);
        m_inEdges.Append(std::set<EdgeId>{});
        m_outEdges.Append(std::set<EdgeId>{});
    }
    else {
        id = m_nodeRecycler.Take();
        m_nodes[id] = id;
        m_inEdges[id] = std::set<EdgeId>{};
        m_outEdges[id] = std::set<EdgeId>{};
    }
    GENERIC_ASSERT(m_inEdges.Size() == m_nodes.Size());
    GENERIC_ASSERT(m_outEdges.Size() == m_nodes.Size());
    return id;
}

EdgeId DirectedGraph::AddEdge(NodeId source, NodeId target)
{
    GENERIC_ASSERT(isValid(source));
    GENERIC_ASSERT(isValid(target));

    EdgeId id;
    if (m_edgeRecycler.Empty()) {
        id = EdgeId(m_edges.Size());
        m_edges.Append(id);
        m_source.Append(source);
        m_target.Append(target);
    }
    else {
        id = m_edgeRecycler.Take();
        m_edges[id] = id;
        m_source[id] = source;
        m_target[id] = target;
    }

    GENERIC_ASSERT(m_source.Size() == m_target.Size());

    m_inEdges[target].insert(id);
    m_outEdges[source].insert(id);

    GENERIC_ASSERT(Source(id) == source);
    GENERIC_ASSERT(Target(id) == target);

    return id;
}

void DirectedGraph::RemoveNode(const NodeId nodeId)
{
    GENERIC_ASSERT(isValid(nodeId));
    for (auto edge : InEdges(nodeId)) {
        if (edge) RemoveEdge(edge);
    }

    for (auto edge : OutEdges(nodeId)) {
        if (edge) RemoveEdge(edge);
    }
    
    m_nodes[nodeId].makeInvalid();
    m_nodeRecycler.Add(nodeId);
}

void DirectedGraph::RemoveEdge(const NodeId source, const NodeId target)
{
    while (auto edgeId = FindEdge(source, target)) {
        if (not isValid(edgeId)) break;
        RemoveEdge(edgeId);
    }
}

void DirectedGraph::RemoveEdge(const EdgeId edgeId)
{
    GENERIC_ASSERT(isValid(edgeId));
    
    auto source = Source(edgeId);
    m_outEdges[source].erase(edgeId);

    auto target = Target(edgeId);
    m_inEdges[target].erase(edgeId);

    m_edges[edgeId].makeInvalid();
    m_edgeRecycler.Add(edgeId);
}

} // namespace generic::graph::model
