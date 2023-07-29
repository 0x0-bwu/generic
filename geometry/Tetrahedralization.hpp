/**
 * @file Tetrahedralization.hpp
 * @author bwu
 * @brief Model of tetrahedralization concept
 * @version 0.1
 * @date 2022-02-22
 */
#pragma once
#include "generic/topology/IndexGraph.hpp"
#include "Point.hpp"
#include <vector>
#include <array>
#include <list>

#if BOOST_SERIALIZATION_SUPPORT
#include "Serialization.hpp"
#endif

namespace generic {
namespace geometry {
///@brief tetrahedralization related classes and functions
namespace tet {

using IndexEle4 = std::array<size_t,  4>;
using IndexEle10 = std::array<size_t, 10>;

template <typename point_t>
class PiecewiseLinearComplex
{
public:
    using IdxPolyLine = std::vector<size_t>;
    struct Surface
    {
        std::vector<IdxPolyLine> faces;
        std::vector<point_t> holes;
    };

    std::vector<point_t> points;
    std::vector<Surface> surfaces;
    
    void Clear();
};

template <typename point_t>
void PiecewiseLinearComplex<point_t>::Clear()
{
    points.clear();
    surfaces.clear();
}

struct IndexVertex;
struct IndexTetrahedron;
using index_t = topology::index_t;
using PosIdx = index_t;
using VerIdx = index_t;
using TetIdx = index_t;

template <typename point_t>
using PointVec = std::vector<point_t>;
using VertexVec = std::vector<IndexVertex>;
using TetIdxVec = std::vector<TetIdx>;
using VerIdxSet = std::unordered_set<VerIdx>;
using TetIdxSet = std::unordered_set<TetIdx>;
using TetrahedronVec = std::vector<IndexTetrahedron>;
using IndexEdge = topology::UndirectedIndexEdge;
using IndexEdgeHash = topology::UndirectedIndexEdgeHash;
using IndexEdgeCompare = topology::UndirectedIndexEdgeCompare;
using IndexFace = std::vector<index_t>;

using topology::noIndex;
inline static constexpr VerIdx noVertex = noIndex;
inline static constexpr TetIdx noNeighbor = noIndex;

struct IndexVertex
{
    PosIdx index;
    TetIdxSet tetrahedrons;

    void RemoveTetrahedron(TetIdx it) { tetrahedrons.erase(it); }

    bool operator != (const IndexVertex & v) const { return index != v.index || tetrahedrons != v.tetrahedrons; }
    bool operator == (const IndexVertex & v) const { return !(*this != v); }

    friend std::ostream & operator<< (std::ostream & os, const IndexVertex & v)
    {
        os << "point index: " << v.index << ", tetrahedrons index:";
        std::for_each(v.tetrahedrons.begin(), v.tetrahedrons.end(), [&](TetIdx i) mutable { os << " " << i; });
        return os;
    }

    void Clear() { tetrahedrons.clear(); }

    static bool isShareEdge(const IndexVertex & a, const IndexVertex & b)
    {
        if(a.tetrahedrons.size() && b.tetrahedrons.size()){
            for(TetIdx it : a.tetrahedrons){
                if(b.tetrahedrons.count(it)) return true;
            }
        }
        return false;
    }

#if BOOST_SERIALIZATION_SUPPORT
private:
    friend class boost::serialization::access;
    template <typename Archive>
    void serialize(Archive & ar, const unsigned int)
    {
        ar & index;
        ar & tetrahedrons;
    }
#endif
};

struct IndexTetrahedron
{
    using IdxEdge = IndexEdge;
    std::array<VerIdx, 4> vertices = {noVertex, noVertex, noVertex, noVertex};
    std::array<TetIdx, 4> neighbors = {noNeighbor, noNeighbor, noNeighbor, noNeighbor};

    void Reset()
    {
        vertices[0] = vertices[1] = vertices[2] = vertices[3] = noVertex;
        neighbors[0] = neighbors[1] = neighbors[2] = neighbors[3] = noNeighbor;
    }

#if BOOST_SERIALIZATION_SUPPORT
private:
    friend class boost::serialization::access;
    template <typename Archive>
    void serialize(Archive & ar, const unsigned int)
    {
        ar & vertices;
        ar & neighbors;
    }
#endif
};

template <typename point_t>
struct Tetrahedralization
{
    using Edge = IndexEdge;
    using Vertex = IndexVertex;
    using Tetrahedron = IndexTetrahedron;
    using EdgeSet = topology::UndirectedIndexEdgeSet;

    EdgeSet fixedEdges;
    VertexVec vertices;
    PointVec<point_t> points;
    TetrahedronVec tetrahedrons;
    
    void Clear()
    {
        points.clear();
        vertices.clear();
        fixedEdges.clear();
        tetrahedrons.clear();
    }

#if BOOST_SERIALIZATION_SUPPORT
private:
    friend class boost::serialization::access;
    template <typename Archive>
    void serialize(Archive & ar, const unsigned int)
    {
        ar & tetrahedrons;
        ar & fixedEdges;
        ar & vertices;
        ar & points;
    }
#endif
};

}//namespace tet
}//namespace geometry
}//namespace generic