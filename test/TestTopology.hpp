/**
 * @file TestTopology.hpp
 * @author bwu
 * @brief Unit test of classes and functions in namespace topology
 * @version 0.1
 * @date 2025-10-29
 */
#pragma once
#include "generic/test/TestCommon.hpp"
#include "generic/topology/IndexGraph.hpp"

using namespace boost::unit_test;
using namespace generic::topology;

void t_topology_undirected_index_edge()
{
    // Test default constructor
    UndirectedIndexEdge e1;
    BOOST_CHECK_EQUAL(e1.v1(), noIndex);
    BOOST_CHECK_EQUAL(e1.v2(), noIndex);
    
    // Test parameterized constructor
    UndirectedIndexEdge e2(5, 3);
    BOOST_CHECK_EQUAL(e2.v1(), 3);  // Should be min
    BOOST_CHECK_EQUAL(e2.v2(), 5);  // Should be max
    
    UndirectedIndexEdge e3(2, 7);
    BOOST_CHECK_EQUAL(e3.v1(), 2);
    BOOST_CHECK_EQUAL(e3.v2(), 7);
    
    // Test SetVertices
    UndirectedIndexEdge e4;
    e4.SetVertices(10, 1);
    BOOST_CHECK_EQUAL(e4.v1(), 1);
    BOOST_CHECK_EQUAL(e4.v2(), 10);
    
    // Test equality
    UndirectedIndexEdge e5(3, 5);
    BOOST_CHECK(e2 == e5);
    BOOST_CHECK(e2 != e3);
    
    // Test hasVertex
    BOOST_CHECK(e2.hasVertex(3));
    BOOST_CHECK(e2.hasVertex(5));
    BOOST_CHECK(!e2.hasVertex(7));
}

void t_topology_undirected_index_edge_hash()
{
    UndirectedIndexEdgeHash hasher;
    
    UndirectedIndexEdge e1(3, 5);
    UndirectedIndexEdge e2(5, 3);  // Same edge, different order
    UndirectedIndexEdge e3(2, 7);  // Different edge
    
    // Same edge should have same hash
    BOOST_CHECK_EQUAL(hasher(e1), hasher(e2));
    
    // Different edges should (likely) have different hashes
    BOOST_CHECK(hasher(e1) != hasher(e3));
}

void t_topology_undirected_index_edge_compare()
{
    UndirectedIndexEdgeCompare compare;
    
    UndirectedIndexEdge e1(3, 5);
    UndirectedIndexEdge e2(5, 3);
    UndirectedIndexEdge e3(2, 7);
    
    BOOST_CHECK(compare(e1, e2));
    BOOST_CHECK(!compare(e1, e3));
}

void t_topology_undirected_index_edge_set()
{
    UndirectedIndexEdgeSet edgeSet;
    
    edgeSet.insert(UndirectedIndexEdge(1, 2));
    edgeSet.insert(UndirectedIndexEdge(2, 3));
    edgeSet.insert(UndirectedIndexEdge(3, 4));
    
    BOOST_CHECK_EQUAL(edgeSet.size(), 3);
    
    // Try inserting duplicate (same edge, different vertex order)
    edgeSet.insert(UndirectedIndexEdge(2, 1));
    BOOST_CHECK_EQUAL(edgeSet.size(), 3);  // Should not add duplicate
    
    // Check if edge exists
    BOOST_CHECK(edgeSet.find(UndirectedIndexEdge(1, 2)) != edgeSet.end());
    BOOST_CHECK(edgeSet.find(UndirectedIndexEdge(5, 6)) == edgeSet.end());
}

void t_topology_sparse_index_graph_basic()
{
    // Test creating empty graph
    SparseIndexGraph graph;
    BOOST_CHECK_EQUAL(boost::num_vertices(graph), 0);
    BOOST_CHECK_EQUAL(boost::num_edges(graph), 0);
}

void t_topology_sparse_index_graph_add_edges()
{
    SparseIndexGraph graph(5);  // Create graph with 5 vertices
    
    AddEdge(0, 1, graph);
    AddEdge(1, 2, graph);
    AddEdge(2, 3, graph);
    AddEdge(3, 4, graph);
    
    BOOST_CHECK_EQUAL(boost::num_vertices(graph), 5);
    BOOST_CHECK_EQUAL(boost::num_edges(graph), 4);
}

void t_topology_sparse_index_graph_edges()
{
    SparseIndexGraph graph(4);
    
    AddEdge(0, 1, graph);
    AddEdge(1, 2, graph);
    AddEdge(2, 3, graph);
    
    auto [begin, end] = Edges(graph);
    
    size_t edgeCount = 0;
    for (auto iter = begin; iter != end; ++iter) {
        auto src = Source(*iter, graph);
        auto tgt = Target(*iter, graph);
        BOOST_CHECK(src < 4);
        BOOST_CHECK(tgt < 4);
        edgeCount++;
    }
    
    BOOST_CHECK_EQUAL(edgeCount, 3);
}

void t_topology_make_sparse_index_graph()
{
    // Create adjacency list
    std::vector<std::set<int>> connection(4);
    connection[0].insert(1);
    connection[1].insert(0);
    connection[1].insert(2);
    connection[2].insert(1);
    connection[2].insert(3);
    connection[3].insert(2);
    
    auto graph = makeSparseIndexGraph(connection);
    
    BOOST_CHECK(graph != nullptr);
    BOOST_CHECK_EQUAL(boost::num_vertices(*graph), 4);
    BOOST_CHECK_EQUAL(boost::num_edges(*graph), 3);
}

void t_topology_connected_component()
{
    SparseIndexGraph graph(6);
    
    // Create two separate components
    // Component 1: 0-1-2
    AddEdge(0, 1, graph);
    AddEdge(1, 2, graph);
    
    // Component 2: 3-4-5
    AddEdge(3, 4, graph);
    AddEdge(4, 5, graph);
    
    // Test ConnectedComponent for vertex 0
    std::vector<index_t> component1;
    ConnectedComponent(graph, 0, component1);
    BOOST_CHECK_EQUAL(component1.size(), 3);
    
    // Test ConnectedComponent for vertex 3
    std::vector<index_t> component2;
    ConnectedComponent(graph, 3, component2);
    BOOST_CHECK_EQUAL(component2.size(), 3);
    
    // Test with std::list
    std::list<index_t> componentList;
    ConnectedComponent(graph, 0, componentList);
    BOOST_CHECK_EQUAL(componentList.size(), 3);
}

void t_topology_connected_components()
{
    SparseIndexGraph graph(8);
    
    // Create three separate components
    // Component 1: 0-1-2
    AddEdge(0, 1, graph);
    AddEdge(1, 2, graph);
    
    // Component 2: 3-4
    AddEdge(3, 4, graph);
    
    // Component 3: 5-6-7
    AddEdge(5, 6, graph);
    AddEdge(6, 7, graph);
    
    std::vector<std::vector<index_t>> components;
    ConnectedComponents(graph, components);
    
    BOOST_CHECK_EQUAL(components.size(), 3);
    
    // Check component sizes
    bool found3 = false, found2 = false, found3b = false;
    for (const auto& comp : components) {
        if (comp.size() == 3) {
            if (!found3) found3 = true;
            else found3b = true;
        }
        else if (comp.size() == 2) found2 = true;
    }
    
    BOOST_CHECK(found3);
    BOOST_CHECK(found2);
    BOOST_CHECK(found3b);
}

void t_topology_edge_set()
{
    SparseIndexGraph graph(5);
    
    AddEdge(0, 1, graph);
    AddEdge(1, 2, graph);
    AddEdge(2, 3, graph);
    AddEdge(3, 4, graph);
    
    auto edgeSet = EdgeSet(graph);
    
    BOOST_CHECK_EQUAL(edgeSet.size(), 4);
    
    // Check if edges are in the set
    BOOST_CHECK(edgeSet.find(UndirectedIndexEdge(0, 1)) != edgeSet.end());
    BOOST_CHECK(edgeSet.find(UndirectedIndexEdge(1, 2)) != edgeSet.end());
    BOOST_CHECK(edgeSet.find(UndirectedIndexEdge(2, 3)) != edgeSet.end());
    BOOST_CHECK(edgeSet.find(UndirectedIndexEdge(3, 4)) != edgeSet.end());
    
    // Check for non-existing edge
    BOOST_CHECK(edgeSet.find(UndirectedIndexEdge(0, 4)) == edgeSet.end());
}

test_suite * create_topology_test_suite()
{
    test_suite * topology_suite = BOOST_TEST_SUITE("s_topology");
    //
    topology_suite->add(BOOST_TEST_CASE(&t_topology_undirected_index_edge));
    topology_suite->add(BOOST_TEST_CASE(&t_topology_undirected_index_edge_hash));
    topology_suite->add(BOOST_TEST_CASE(&t_topology_undirected_index_edge_compare));
    topology_suite->add(BOOST_TEST_CASE(&t_topology_undirected_index_edge_set));
    topology_suite->add(BOOST_TEST_CASE(&t_topology_sparse_index_graph_basic));
    topology_suite->add(BOOST_TEST_CASE(&t_topology_sparse_index_graph_add_edges));
    topology_suite->add(BOOST_TEST_CASE(&t_topology_sparse_index_graph_edges));
    topology_suite->add(BOOST_TEST_CASE(&t_topology_make_sparse_index_graph));
    topology_suite->add(BOOST_TEST_CASE(&t_topology_connected_component));
    topology_suite->add(BOOST_TEST_CASE(&t_topology_connected_components));
    topology_suite->add(BOOST_TEST_CASE(&t_topology_edge_set));
    //
    return topology_suite;
}
