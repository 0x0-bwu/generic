/**
 * @file BoostGraphTraits.hpp
 * @author bwu
 * @brief boost graph adapter
 * @version 0.1
 * @date 2024-11-29
 */
#pragma once
#include "generic/graph/model/DirectedGraph.hpp"
#include <boost/graph/graph_traits.hpp>
#include <tuple>
namespace boost {

using namespace generic::graph::model;

struct GenericDirectedGraphTraversalCategory
 : public virtual bidirectional_graph_tag
 , public virtual vertex_list_graph_tag
 , public virtual edge_list_graph_tag
{
};

template <>
struct graph_traits<DirectedGraph>
{
    using vertex_descriptor = NodeId;
    using edge_descriptor = EdgeId;

    using vertex_iterator = typename DirectedGraph::NodeIterator;
    using edge_iterator = typename DirectedGraph::EdgeIterator;

    using vertices_size_type = typename NodeId::SizeType;
    using degree_size_type = typename EdgeId::SizeType;
    using edges_size_type = typename EdgeId::SizeType;

    using edge_parallel_category = allow_parallel_edge_tag;
    using traversal_category = GenericDirectedGraphTraversalCategory;

    using directed_category = directed_tag ;
    using in_edge_iterator = typename DirectedGraph::InEdgeIterator;
    using out_edge_iterator = typename DirectedGraph::OutEdgeIterator;

/**
 * @brief Brief description of null_vertex.
 * @return static vertex_descriptor
 */
    static vertex_descriptor null_vertex() { return NodeId(); }
/**
 * @brief Brief description of null_edge.
 * @return static edge_descriptor
 */
    static edge_descriptor null_edge() { return EdgeId(); }
};


// IncidenceGraph APIs
inline std::pair<typename graph_traits<DirectedGraph>::vertex_iterator,
                 typename graph_traits<DirectedGraph>::vertex_iterator>
vertices(const DirectedGraph & g)
{
    return g.AllNodes().toStdPair();
}

inline typename graph_traits<DirectedGraph>::vertices_size_type
num_vertices(const DirectedGraph & g)
{
    return g.Nodes();
}

inline std::pair<typename graph_traits<DirectedGraph>::edge_iterator,
                 typename graph_traits<DirectedGraph>::edge_iterator>
edges(const DirectedGraph & g)
{
    return g.AllEdges().toStdPair();
}

inline typename graph_traits<DirectedGraph>::edges_size_type
num_edges(const DirectedGraph & g)
{
    return g.Edges();
}

inline typename graph_traits<DirectedGraph>::vertex_descriptor
source(typename graph_traits<DirectedGraph>::edge_descriptor e, const DirectedGraph & g)
{
    return g.Source(e);
}

inline typename graph_traits<DirectedGraph>::vertex_descriptor
target(typename graph_traits<DirectedGraph>::edge_descriptor e, const DirectedGraph & g)
{
    return g.Target(e);
}

inline std::pair<typename graph_traits<DirectedGraph>::out_edge_iterator, 
                 typename graph_traits<DirectedGraph>::out_edge_iterator>
out_edges(typename graph_traits<DirectedGraph>::vertex_descriptor v, const DirectedGraph & g)
{
    return g.OutEdges(v).toStdPair();
}

inline typename graph_traits<DirectedGraph>::degree_size_type
out_degree(typename graph_traits<DirectedGraph>::vertex_descriptor v, const DirectedGraph & g)
{
    return g.OutEdges(v).Size();
}

// BidirectionalGraph APIs
inline std::pair<typename graph_traits<DirectedGraph>::in_edge_iterator, 
                 typename graph_traits<DirectedGraph>::in_edge_iterator>
in_edges(typename graph_traits<DirectedGraph>::vertex_descriptor v, const DirectedGraph & g)
{
    return g.InEdges(v).toStdPair();
}

inline typename graph_traits<DirectedGraph>::degree_size_type
in_degree(typename graph_traits<DirectedGraph>::vertex_descriptor v, const DirectedGraph & g)
{
    return g.InEdges(v).Size();
}

inline typename graph_traits<DirectedGraph>::degree_size_type
degree(typename graph_traits<DirectedGraph>::vertex_descriptor v, const DirectedGraph & g)
{
    return g.InEdges(v).Size() + g.OutEdges(v).Size();
}

// MutableGraph APIs
inline std::pair<typename graph_traits<DirectedGraph>::edge_descriptor, bool>
add_edge(typename graph_traits<DirectedGraph>::vertex_descriptor u, typename graph_traits<DirectedGraph>::vertex_descriptor v, DirectedGraph & g)
{
    return {g.AddEdge(u, v), true};
}

inline void remove_edge(typename graph_traits<DirectedGraph>::vertex_descriptor u, typename graph_traits<DirectedGraph>::vertex_descriptor v, DirectedGraph & g)
{
    return g.RemoveEdge(u, v);
}

inline void remove_edge(typename graph_traits<DirectedGraph>::edge_descriptor e, DirectedGraph & g)
{
    return g.RemoveEdge(e);
}

inline void remove_edge(typename graph_traits<DirectedGraph>::edge_iterator e, DirectedGraph & g)
{
    return g.RemoveEdge(*e);
}

inline typename graph_traits<DirectedGraph>::vertex_descriptor
add_vertex(DirectedGraph & g)
{
    return g.AddNode();
}

inline void clear_vertex(typename graph_traits<DirectedGraph>::vertex_descriptor u, DirectedGraph & g)
{
    for (auto edge : g.OutEdges(u)) {
        if (g.isValid(edge)) g.RemoveEdge(edge);
    }
}

inline void remove_vertex(typename graph_traits<DirectedGraph>::vertex_descriptor u, DirectedGraph & g)
{
    return g.RemoveNode(u);
}

} // namespace boost
