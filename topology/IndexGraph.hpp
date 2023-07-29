/**
 * @file IndexGraph.hpp
 * @author bwu
 * @brief Model of index graph concept and graph related algorithms
 * @version 0.1
 * @date 2022-02-22
 */
#pragma once
#include "generic/common/Exception.hpp"
#include "generic/common/Archive.hpp"
#include "generic/common/Traits.hpp"
#include "Common.hpp"
#include <boost/graph/connected_components.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/functional/hash.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/visitors.hpp>
#include <unordered_set>

namespace generic  {
namespace topology {

using namespace common;

///@brief represents model of undirected index edge concept
struct UndirectedIndexEdge
{
    ///@brief constructs an invalid edge
    UndirectedIndexEdge()
    {
        m_vertices = std::make_pair(noIndex, noIndex);
    }

    ///@brief constructs an undirected edge with vertex `iv1` and `iv2`
    UndirectedIndexEdge(index_t iv1, index_t iv2)
    {
       SetVertices(iv1, iv2);
    }

    ///@brief checks if this edge equals to `e`
    bool operator==(const UndirectedIndexEdge & e) const { return m_vertices == e.m_vertices; }
    ///@brief checks if this edge not equals to `e`
    bool operator!=(const UndirectedIndexEdge & e) const { return !(*this == e); }

    ///@brief set two vertices index of this edges
    void SetVertices(index_t iv1, index_t iv2) { m_vertices.first = std::min(iv1, iv2); m_vertices.second = std::max(iv1, iv2); }

    ///@brief returns the smaller vertex index of this edge
    index_t v1() const { return m_vertices.first ; }
    ///@brief returns the bigger vertex index of this edge
    index_t v2() const { return m_vertices.second; }

    ///@brief chechks if this edge contains vertex `iv`
    bool hasVertex(index_t iv) const { return iv == m_vertices.first || iv == m_vertices.second; }

private:
    std::pair<index_t, index_t> m_vertices;

#if BOOST_SERIALIZATION_SUPPORT
private:
    friend class boost::serialization::access;
    template <typename Archive>
    void serialize(Archive & ar, const unsigned int)
    {
        ar & m_vertices;
    }
#endif
};

struct UndirectedIndexEdgeHash
{
    size_t operator() (const UndirectedIndexEdge & edge) const noexcept
    {
        size_t seed(0);
        boost::hash_combine(seed, edge.v1());
        boost::hash_combine(seed, edge.v2());
        return seed;
    }
};

struct UndirectedIndexEdgeCompare
{
    bool operator() (const UndirectedIndexEdge & e1, const UndirectedIndexEdge & e2) const noexcept { return e1 == e2; }
};

using UndirectedIndexEdgeSet = std::unordered_set<UndirectedIndexEdge, UndirectedIndexEdgeHash, UndirectedIndexEdgeCompare>;

template <typename T>
using UndirectedIndexEdgeMap = std::unordered_map<UndirectedIndexEdge, T, UndirectedIndexEdgeHash, UndirectedIndexEdgeCompare>;

///@brief represents model of sparse undirect index graph concept
using SparseIndexGraph = boost::adjacency_list<
                         boost::setS,
                         boost::vecS,
                         boost::undirectedS,
                         boost::property<boost::vertex_index_t, index_t >
                         >;

using SIGEdge = boost::graph_traits<SparseIndexGraph>::edge_descriptor; 
using SIGEdgeIter = boost::graph_traits<SparseIndexGraph>::edge_iterator;
using SIGVertex = boost::graph_traits<SparseIndexGraph>::vertex_descriptor;

inline void AddEdge(index_t i, index_t j, SparseIndexGraph & g)
{
    boost::add_edge(i, j, g);
}

inline std::pair<SIGEdgeIter, SIGEdgeIter> Edges(const SparseIndexGraph & g)
{
    return boost::edges(g);
}

inline SIGVertex Source(const SIGEdge & e, const SparseIndexGraph & g)
{
    return boost::source(e, g);
}

inline SIGVertex Target(const SIGEdge & e, const SparseIndexGraph & g)
{
    return boost::target(e, g);
}

///@brief construct a sparse index grahp from adjecent list formed by vector of set
inline std::unique_ptr<SparseIndexGraph> makeSparseIndexGraph(const std::vector<std::set<int> > & connection)
{
    size_t vertex = connection.size();
    auto graph = std::make_unique<SparseIndexGraph>();
    for(size_t v = 0; v < vertex; ++v) {
        const auto & adjecent = connection[v];
        for(int w : adjecent) {
            if(w < 0 || w >= vertex) continue;
                AddEdge(v, static_cast<index_t>(w), *graph);
        }
    }

    return graph;
}

/**
 * @brief gets connected component by BFS
 * 
 * @param[in] g the connection graph
 * @param[in] v the source vertex
 * @param[out] c contains the vertices connected to source v
 * @note `c` should be one of the std containers
 */
template <template <typename, typename> class container,
          template <typename> class allocator = std::allocator>
inline void ConnectedComponent(const SparseIndexGraph & g, const index_t v, container<index_t, allocator<index_t> > & c)
{
    using namespace boost;
    using namespace common;
    GENERIC_ASSERT(v < num_vertices(g))

    using Component = container<index_t, allocator<index_t> >;
    class BFSVisitor : public default_bfs_visitor
    {
    public:
        Component & visited;
        BFSVisitor(Component & _visited) : visited(_visited) {}
        void discover_vertex(index_t s, const SparseIndexGraph &) { append_to_std_container(visited, s, 0); }
    };

    c.clear();
    BFSVisitor vis(c);
    breadth_first_search(g, v, visitor(vis));
}

/**
 * @brief gets all connected components from graph
 * 
 * @param[in] g the connection graph
 * @param[out] cc connected components
 * @note the value type of `cc` should be one of the std containers
 */
template <template <typename, typename> class container,
          template <typename> class allocator = std::allocator>
inline void ConnectedComponents(const SparseIndexGraph & g, std::vector<container<index_t, allocator<index_t> > > & cc)
{
    using namespace boost;
    std::vector<index_t> c(num_vertices(g));
    auto numComp = connected_components(g, make_iterator_property_map(c.begin(), get(vertex_index, g)));

    cc.clear();
    cc.resize(numComp);
    for(size_t i = 0; i < c.size(); ++i)
        append_to_std_container(cc[c[i]], i, 0);
}

inline UndirectedIndexEdgeSet EdgeSet(const SparseIndexGraph & g)
{
    UndirectedIndexEdgeSet edgeSet;
    auto [begin, end] = topology::Edges(g);
    for(auto iter = begin; iter != end; ++iter){
        auto u = topology::Source(*iter, g);
        auto v = topology::Target(*iter, g);
        edgeSet.insert(UndirectedIndexEdge{u, v});
    }
    return edgeSet;
}

}//namespace topology
}//namespace generic