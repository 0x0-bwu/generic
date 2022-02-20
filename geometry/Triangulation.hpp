/**
 * @file Triangulation.hpp
 * @author bwu
 * @brief Model of triangulation concept
 * @version 0.1
 * @date 2022-02-22
 */
#ifndef GENERIC_GEOMETRY_TRI_TRIANGULATION_HPP
#define GENERIC_GEOMETRY_TRI_TRIANGULATION_HPP
#include "generic/topology/IndexGraph.hpp"
#include "generic/math/MathUtility.hpp"
#include "generic/common/Traits.hpp"
#include "GeometryIO.hpp"
#include "Utility.hpp"
#include "Common.hpp"
#include "Point.hpp"
#include "Box.hpp"
#include <unordered_set>
#include <utility>
#include <vector>
#include <array>

#ifdef BOOST_SERIALIZATION_SUPPORT
#include "Serialization.hpp"
#endif

namespace generic  {
namespace geometry {
namespace tri {
struct IndexVertex;
struct IndexTriangle;
using index_t = topology::index_t;
using PosIdx = index_t;
using VerIdx = index_t;
using TriIdx = index_t;

template <typename point_t>
using PointVec = std::vector<point_t>;
using VertexVec = std::vector<IndexVertex>;
using TriIdxVec = std::vector<TriIdx>;
using VerIdxSet = std::unordered_set<VerIdx>;
using TriIdxSet = std::unordered_set<TriIdx>;
using TriangleVec = std::vector<IndexTriangle>;
using IndexEdge = topology::UndirectedIndexEdge;
using IndexEdgeHash = topology::UndirectedIndexEdgeHash;
using IndexEdgeCompare = topology::UndirectedIndexEdgeCompare;

using topology::noIndex;
inline static constexpr VerIdx noVertex = noIndex;
inline static constexpr TriIdx noNeighbor = noIndex;

struct IndexVertex
{
    PosIdx index;
    TriIdxSet triangles;

    void RemoveTriangle(TriIdx it) { triangles.erase(it); }

    bool operator != (const IndexVertex & v) const { return index != v.index || triangles != v.triangles; }
    bool operator == (const IndexVertex & v) const { return !(*this != v); }

    friend std::ostream & operator<< (std::ostream & os, const IndexVertex & v)
    {
        os << "point index: " << v.index << ", triangle index:";
        std::for_each(v.triangles.begin(), v.triangles.end(), [&](TriIdx i) mutable { os << " " << i; });
        return os;
    }

    void Clear() { triangles.clear(); }

    static bool isShareEdge(const IndexVertex & a, const IndexVertex & b)
    {
        if(a.triangles.size() && b.triangles.size()){
            for(TriIdx it : a.triangles){
                if(b.triangles.count(it)) return true;
            }
        }
        return false;
    }

#ifdef BOOST_SERIALIZATION_SUPPORT
private:
    friend class boost::serialization::access;
    template <typename Archive>
    void serialize(Archive & ar, const unsigned int)
    {
        ar & index;
        ar & triangles;
    }
#endif
};

struct IndexTriangle
{
    using IdxEdge = IndexEdge;
    std::array<VerIdx, 3> vertices = {noVertex, noVertex, noVertex};
    std::array<TriIdx, 3> neighbors = {noNeighbor, noNeighbor, noNeighbor};

    static index_t cw(const index_t i)  { return (i + 2) % 3; }
    static index_t ccw(const index_t i) { return (i + 1) % 3; }
    ///@brief opposed triangle neighbor index from vertex index iv
    static index_t iVeOpiTn(const index_t iv) { return (iv + 1) % 3; }
    ///@brief opposed vertex index from triangle neighbor index in
    static index_t iTnOpiVe(const index_t in) { return (in + 2) % 3; }
    ///@brief opposed vertex index from edge index ie
    static index_t iEgOpiVe(const index_t ie) { return (ie + 2) % 3; }
    ///@brief opposed edge index from vertex index iv
    static index_t iVeOpiEg(const index_t iv) { return (iv + 1) % 3; }
    index_t iVe(const VerIdx v) const
    {
        for(index_t i = 0; i < 3; ++i){ if(v == vertices[i]) return i; }
        throw std::runtime_error("Could not find vertex index in triangle");
    }

    index_t iEg(const IdxEdge & e) const
    {
        for(index_t i = 0; i < 3; ++i) { if(IdxEdge(vertices[i], vertices[ccw(i)]) == e) return i; }
        throw std::runtime_error("Could not find edge index in triangle");
    }

    index_t iTn(const TriIdx t) const
    {
        for(index_t i = 0; i < 3; ++i){ if(t == neighbors[i]) return i; }
        throw std::runtime_error("Could not find neighbor triangle index");
    }

    index_t iTn(const VerIdx v1, const VerIdx v2) const
    {
        for(index_t i = 0; i < 3; ++i){ if(v1 != vertices[i] && v2 != vertices[i]) return iVeOpiTn(i); }
        throw std::runtime_error("Could not find opposed-to-edge triangle index");
    }

    index_t VeOpiTn(const VerIdx vOp) const { return iVeOpiTn(iVe(vOp)); }
    index_t TnOpiVe(const TriIdx tOp) const { return iTnOpiVe(iTn(tOp)); }
    index_t VeOpiEg(const VerIdx vOp) const { return iVeOpiEg(iVe(vOp)); }
    index_t EgOpiVe(const IdxEdge & e) const { return iEgOpiVe(iEg(e)); }
    
    VerIdx TnOpVe(const TriIdx tOp) const { return vertices[TnOpiVe(tOp)]; }
    VerIdx EgOpVe(const IdxEdge & e) const {return vertices[EgOpiVe(e)]; }
    VerIdx NextVertex(const VerIdx v, WindingDirection dir = WindingDirection::CounterClockwise) const
    { return vertices[dir == WindingDirection::Clockwise ? cw(iVe(v)) : ccw(iVe(v))]; }

    TriIdx Tn(const IdxEdge & e) const { return Tn(e.v1(), e.v2()); }
    TriIdx Tn(const VerIdx v1, const VerIdx v2) const { return neighbors[iTn(v1, v2)]; }
    TriIdx VeOpTn(const VerIdx vOp) const { return neighbors[VeOpiTn(vOp)]; }
    IdxEdge VeOpEg(const VerIdx vOp) const { return Edge(VeOpiEg(vOp)); }
    IdxEdge VeCwEdge(const VerIdx v) const { return Edge(iVeOpiEg(ccw(iVe(v)))); }
    IdxEdge VeCcwEdge(const VerIdx v) const { return Edge(iVeOpiEg(cw(iVe(v)))); }
    IdxEdge Edge(const index_t i) const { return IdxEdge(vertices[i], vertices[ccw(i)]); }
    std::pair<VerIdx, VerIdx> VeOpVes(const VerIdx v) const { auto i = iVe(v); return std::make_pair(vertices[cw(i)], vertices[ccw(i)]); }//clock wise
    std::pair<TriIdx, TriIdx> EgOpTns(const IdxEdge & e) const { auto i = EgOpiVe(e); return std::make_pair(neighbors[iVeOpiTn(cw(i))], neighbors[iVeOpiTn(ccw(i))]); }

    void ChangeVertex(const VerIdx from, const VerIdx to) { vertices[iVe(from)] = to; }
    void ChangeNeighbor(const TriIdx from, const TriIdx to) { neighbors[iTn(from)] = to; }
    void ChangeNeighbor(const VerIdx v1, const VerIdx v2, const TriIdx to) { neighbors[iTn(v1, v2)] = to; }

    bool hasNeighbor(const TriIdx neighbor) const
    {
        for(TriIdx n : neighbors){
            if(n == neighbor) return true;
        }
        return false;
    }

    bool operator != (const IndexTriangle & t) const { return vertices != t.vertices || neighbors != t.neighbors; }
    bool operator == (const IndexTriangle & t) const { return !(*this != t); }

    friend std::ostream & operator<< (std::ostream & os, const IndexTriangle & t)
    {
        os << "vertex:";
        std::for_each(t.vertices.begin(), t.vertices.end(),[&](VerIdx i) mutable { os << " " << i; });
        os << ", neighbor:";
        std::for_each(t.neighbors.begin(), t.neighbors.end(),[&](TriIdx n) mutable { os << " " << n; });
        return os;
    }

    void Reverse()
    {
        std::swap(vertices[1], vertices[2]);
        std::swap(neighbors[0], neighbors[2]);
    }

    void Reset()
    {
        vertices[0] = vertices[1] = vertices[2] = noVertex;
        neighbors[0] = neighbors[1] = neighbors[2] = noNeighbor;
    }

#ifdef BOOST_SERIALIZATION_SUPPORT
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
struct Triangulation
{
    using Edge = IndexEdge;
    using Vertex = IndexVertex;
    using Triangle = IndexTriangle;
    using EdgeSet = topology::UndirectedIndexEdgeSet;

    EdgeSet fixedEdges;
    VertexVec vertices;
    TriangleVec triangles;
    PointVec<point_t> points;

    const point_t & Point(const index_t i) const { return points[i]; }
    const point_t & VerIdxPoint(const VerIdx iv) const { return VertexPoint(vertices[iv]); }
    const point_t & VertexPoint(const Vertex & v) const { return Point(v.index); }
    const point_t & TriIdxPoint(const TriIdx it, const index_t i) const { return TrianglePoint(triangles[it], i); }
    const point_t & TrianglePoint(const Triangle & t, const index_t i) const { return VerIdxPoint(t.vertices[i % 3]); }

    void AddAdjacentTriangle(const VerIdx iv, const TriIdx it)
    {
        vertices[iv].triangles.insert(it);
    }

    void RemoveAdjacentTriangle(const VerIdx iv, const TriIdx it)
    {
        vertices[iv].triangles.erase(it);
    }

    void ChangeNeighbor(const TriIdx it, const TriIdx from, const TriIdx to){ if(it != noNeighbor) triangles[it].ChangeNeighbor(from, to); }
    void ChangeNeighbor(const TriIdx it, const VerIdx iv1, const VerIdx iv2, const TriIdx neighbor) { if(it != noNeighbor) triangles[it].ChangeNeighbor(iv1, iv2, neighbor); }
    bool hasNeighbor(const TriIdx it, const TriIdx neighbor) const { return triangles[it].hasNeighbor(neighbor); }
    bool hasFixedEdge(const Edge & e) const { return fixedEdges.count(e); }
    void FlipEdge(const TriIdx t, const TriIdx tOp)
    {
        Triangle & triangle   = triangles[t  ];
        Triangle & triangleOp = triangles[tOp];

        index_t iv  = triangle.TnOpiVe(tOp);
        index_t ivOp = triangleOp.TnOpiVe(t);

        const auto & triNs   = triangle.neighbors;
        const auto & triOpNs = triangleOp.neighbors;
        const auto & triVs   = triangle.vertices;
        const auto & triOpVs = triangleOp.vertices;

        auto v1 = triVs[iv];
        auto v2 = triVs[Triangle::ccw(iv)];
        auto n1 = triNs[iv];
        auto n3 = triNs[Triangle::cw(iv)];
        auto v3 = triOpVs[ivOp];
        auto v4 = triOpVs[Triangle::ccw(ivOp)];
        auto n4 = triOpNs[ivOp];
        auto n2 = triOpNs[Triangle::cw(ivOp)];

        triangle.vertices  = std::array<VerIdx, 3>{v4, v1, v3};
        triangle.neighbors = std::array<TriIdx, 3>{n3, tOp, n4};

        triangleOp.vertices  = std::array<VerIdx, 3>{v2, v3, v1};
        triangleOp.neighbors = std::array<TriIdx, 3>{n2, t, n1};

        ChangeNeighbor(n1, t, tOp);
        ChangeNeighbor(n4, tOp, t);
        AddAdjacentTriangle(v1, tOp);
        AddAdjacentTriangle(v3, t);
        RemoveAdjacentTriangle(v2, t);
        RemoveAdjacentTriangle(v4, tOp);
    }

    std::list<Edge> ConcentricFixedEdges(const VerIdx & iv) const
    {
        std::list<Edge> edges;
        const auto & v = vertices[iv];
        for(TriIdx it : v.triangles){
            const auto & triangle = triangles[it];
            Edge e = triangle.VeCwEdge(iv);
            if(fixedEdges.count(e))
                edges.emplace_back(std::move(e));
        }
        return edges;
    }

    std::pair<TriIdx, TriIdx> GetTriangles(const Edge & e) const
    {
        auto tris = std::make_pair(noNeighbor, noNeighbor);
        const auto & v1 = vertices[e.v1()];
        const auto & v2 = vertices[e.v2()];
        for(TriIdx it1 : v1.triangles){
            for(TriIdx it2 : v2.triangles){
                if(it1 == it2){
                    if(noNeighbor == tris.first)
                        tris.first = it1;
                    else tris.second = it2;
                }
            }
        }
        return tris;
    }

    void Clear()
    {
        points.clear();
        vertices.clear();
        triangles.clear();
        fixedEdges.clear();
    }

#ifdef BOOST_SERIALIZATION_SUPPORT
private:
    friend class boost::serialization::access;
    template <typename Archive>
    void serialize(Archive & ar, const unsigned int)
    {
        ar & fixedEdges;
        ar & vertices;
        ar & triangles;
        ar & points;
    }
#endif
};

using geometry::segment_type;
using geometry::triangle_type;
template <typename point_t>
class TriangulationUtility
{
    using coor_t = typename point_t::coor_t;
    using float_t = float_type<typename point_t::coor_t>;
public:
    static std::array<point_t, 3> GetVertexPoints(const Triangulation<point_t> & tri, const TriIdx it)
    {
        return { tri.TriIdxPoint(it, 0), tri.TriIdxPoint(it, 1), tri.TriIdxPoint(it, 2) };
    }

    static point_t GetVertexPoint(const Triangulation<point_t> & tri, const VerIdx iv)
    {
        return tri.VerIdxPoint(iv);
    }

    static coor_t GetEdgeLenSq(const Triangulation<point_t> & tri, const IndexEdge & e)
    {
        return geometry::DistanceSq(GetVertexPoint(tri, e.v1()), GetVertexPoint(tri, e.v2()));
    }

    static coor_t GetEdgeLenSq(const Triangulation<point_t> & tri, const TriIdx it, index_t ie)
    {
        return GetEdgeLenSq(tri, tri.triangles[it].Edge(ie));
    }

    static float_t GetEdgeLength(const Triangulation<point_t> & tri, const IndexEdge & e)
    {
        return std::sqrt(GetEdgeLenSq(tri, e));
    }

    static float_t GetEdgeLength(const Triangulation<point_t> & tri, const TriIdx it, index_t ie)
    {
        return std::sqrt(GetEdgeLenSq(tri, it, ie));
    }

    static float_t GetAverageEdgeLength(const Triangulation<point_t> & tri, const TriIdx it)
    {
        float_t len = 0, inv_3 = 1.0 / 3;
        for(index_t ie = 0; ie < 3; ++ie)
            len += inv_3 * GetEdgeLength(tri, it, ie);
        return len;
    }

    static std::array<float_t, 3> GetInteriorAngles(const Triangulation<point_t> & tri, const TriIdx it)
    {
        std::array<point_t, 3> pts = GetVertexPoints(tri, it);
        return geometry::InteriorAngles(pts[0], pts[1], pts[2]);
    }

    static float_t GetMinimumAngle(const Triangulation<point_t> & tri, const TriIdx it)
    {
        auto angles = GetInteriorAngles(tri, it);
        float_t min = std::min(angles[0], angles[1]);
        return std::min(min, angles[2]);
    }

    static float_t GetMaximumAngle(const Triangulation<point_t> & tri, const TriIdx it)
    {
        auto angles = GetInteriorAngles(tri, it);
        float_t max = std::max(angles[0], angles[1]);
        return std::max(max, angles[2]);
    }

    static IndexEdge GetMaxLenEdge(const Triangulation<point_t> & tri, const TriIdx it)
    {
        index_t max = 0;
        coor_t maxDistSq = 0;
        for(size_t i = 0; i < 3; ++i){
            auto distSq = GetEdgeLenSq(tri, it, i);
            if(distSq > maxDistSq){
                maxDistSq = distSq;
                max = i;
            }
        }
        return tri.triangles[it].Edge(max);
    }

    static IndexEdge GetMinLenEdge(const Triangulation<point_t> & tri, const TriIdx it)
    {
        index_t min = 0;
        coor_t minDistSq = std::numeric_limits<coor_t>::max();
        for(size_t i = 0; i < 3; ++i){
            auto distSq = GetEdgeLenSq(tri, it, i);
            if(distSq < minDistSq){
                minDistSq = distSq;
                min = i;
            }
        }
        return tri.triangles[it].Edge(min);
    }

    static coor_t GetMaxEdgeLenSq(const Triangulation<point_t> & tri, const TriIdx it)
    {
        return GetEdgeLenSq(tri, GetMaxLenEdge(tri, it));
    }

    static coor_t GetMinEdgeLenSq(const Triangulation<point_t> & tri, const TriIdx it)
    {
        return GetEdgeLenSq(tri, GetMinLenEdge(tri, it));
    }

    static float_t GetMaxEdgeLen(const Triangulation<point_t> & tri, const TriIdx it)
    {
        return std::sqrt(GetMaxEdgeLenSq(tri, it));
    }

    static float_t GetMinEdgeLen(const Triangulation<point_t> & tri, const TriIdx it)
    {
        return std::sqrt(GetMinEdgeLenSq(tri, it));
    }

    static segment_type<point_t> GetSegment(const Triangulation<point_t> & tri, const IndexEdge & e)
    {
        return segment_type<point_t>(GetVertexPoint(tri, e.v1()), GetVertexPoint(tri, e.v2()));
    }

    template <typename std::enable_if<point_t::dim == 2, bool>::type = true>
    static Circle<float_t> GetDiametralCircle(const Triangulation<point_t> & tri, const IndexEdge & e)
    {
        return geometry::DiametralCircle(GetSegment(tri, e));
    }

    template <typename std::enable_if<point_t::dim == 2, bool>::type = true>
    static Circle<float_t> GetCircumCircle(const Triangulation<point_t> & tri, const TriIdx it)
    {
        auto [p1, p2, p3] = GetVertexPoints(tri, it);
        return CircumCircle(p1, p2, p3);
    }

    template <typename std::enable_if<point_t::dim == 2, bool>::type = true>
    static bool isInCircumCircle(const Triangulation<point_t> & tri, const TriIdx it, const point_t & p)
    {
        auto [p1, p2, p3] = GetVertexPoints(tri, it);
        return geometry::isInCircumCircle(p1, p2, p3, p, false);
    }

    static triangle_type<point_t> GetTriangle(const Triangulation<point_t> & tri, const TriIdx it)
    {
        auto [p1, p2, p3] = GetVertexPoints(tri, it);
        return triangle_type<point_t>(p1, p2, p3);
    }

    static float_t GetTriangleArea(const Triangulation<point_t> & tri, const TriIdx it)
    {
        return GetTriangle(tri, it).Area();
    }
    
    template <typename std::enable_if<point_t::dim == 2, bool>::type = true>
    static bool isCCW(const Triangulation<point_t> & tri, const TriIdx it)
    {
        const auto & t = tri.triangles[it];
        return isCCW(tri, t.vertices[0], t.vertices[1], t.vertices[2]);
    }

    template <typename std::enable_if<point_t::dim == 2, bool>::type = true>
    static bool isCCW(const Triangulation<point_t> & tri, const VerIdx iv1, const VerIdx iv2, const VerIdx iv3)
    {
        const auto & p1 = GetVertexPoint(tri, iv1);
        const auto & p2 = GetVertexPoint(tri, iv2);
        const auto & p3 = GetVertexPoint(tri, iv3);
        return triangle_type<point_t>::isCCW(p1, p2, p3);
    }

    template <typename std::enable_if<point_t::dim == 2, bool>::type = true>
    static bool isLocallyDelaunay(const Triangulation<point_t> & tri, const IndexEdge & e)
    {
        auto [it1, it2] = tri.GetTriangles(e);
        if(noNeighbor = it2) return true;
        return isLocallyDelaunay(it1, it2);
    }

    template <typename std::enable_if<point_t::dim == 2, bool>::type = true>
    static bool isLocallyDelaunay(const Triangulation<point_t> & tri, const TriIdx it, const TriIdx itOp)
    {
        GENERIC_ASSERT(tri.hasNeighbor(it, itOp))

        const auto & t = tri.triangles[it];
        const auto & tOp = tri.triangles[itOp];

        VerIdx iv = t.TnOpVe(itOp);
        VerIdx ivOp = tOp.TnOpVe(it);
        if(isInCircumCircle(tri, it, tri.VerIdxPoint(ivOp))) return false;
        if(isInCircumCircle(tri, itOp, tri.VerIdxPoint(iv))) return false;
        return true;
    }

    static point_f<point_t> GetCenter(const Triangulation<point_t> & tri, const TriIdx it)
    {
        return GetTriangle(tri, it).Center();
    }

    template <typename std::enable_if<point_t::dim == 2, bool>::type = true>
    static point_f<point_t> GetCircumCenter(const Triangulation<point_t> & tri, const TriIdx it)
    {
        return GetCircumCircle(tri, it).o;
    }

    template <typename point_t2, typename std::enable_if<point_t::dim == 2 && point_t2::dim == 2, bool>::type = true>
    static PointTriangleLocation GetPointTriangleLocation(const Triangulation<point_t> & tri, const TriIdx it, const point_t2 & p)
    {
        auto [p1, p2, p3] = GetVertexPoints(tri, it);
        return geometry::GetPointTriangleLocation(p, p1, p2, p3);
    }

    template <typename std::enable_if<point_t::dim == 2, bool>::type = true>
    static float_t CircumRadius2ShortestEdgeRatioSq(const Triangulation<point_t> & tri, const TriIdx it)
    {
        const auto & p1 = tri.TriIdxPoint(it, 0);
        const auto & p2 = tri.TriIdxPoint(it, 1);
        const auto & p3 = tri.TriIdxPoint(it, 2);
        return geometry::CircumRadius2ShortestEdgeRatioSq(p1, p2, p3);
    }

    template <typename std::enable_if<point_t::dim == 2, bool>::type = true>
    static void GetNeighborVerticesCCW(const Triangulation<point_t> & tri, const VerIdx iv, std::vector<VerIdx> & neighbors)
    {
        const auto & vertex = tri.vertices[iv];
        const auto & triangles = vertex.triangles;

        neighbors.clear();
        neighbors.reserve(triangles.size());
        
        //start from boundary edge if it's not inner vertex
        TriIdx start;
        for(auto it : triangles){
            start = it;
            const auto & startT = tri.triangles.at(start);
            auto nextV = startT.NextVertex(iv);
            if(noNeighbor == startT.Tn(iv, nextV)) break;
        }

        auto currT = tri.triangles[start];
        while(neighbors.size() != triangles.size()){
            auto nextV = currT.NextVertex(iv);
            neighbors.push_back(nextV);
            
            nextV = currT.NextVertex(nextV);
            currT = tri.triangles[currT.Tn(iv, nextV)];
        }
    }

#ifdef BOOST_SERIALIZATION_SUPPORT
    static bool Write(const Triangulation<point_t> & tri, const std::string & archive)
    {
        std::ofstream ofs(archive);
        if(!ofs.is_open()) return false;
        boost::archive::text_oarchive oa(ofs);
        oa << tri;
        return true;
    }

    static bool Read(const std::string & archive, Triangulation<point_t> & tri)
    {
        std::ifstream ifs(archive);
        if(!ifs.is_open()) return false;
        boost::archive::text_iarchive ia(ifs);
        tri.Clear();
        ia >> tri;
        return true;
    }
#endif
};

class VertexGraph
{
public:
    using Edge = IndexEdge;
    using Edges = std::list<Edge>;

    VertexGraph() = default;

    template <typename point_t>
    explicit VertexGraph(const Triangulation<point_t> & tri, const TriIdxSet & ignored = {});
    void AddEdge(VerIdx i, VerIdx j);
    void GetEdges(Edges & edges) const;

private:
    topology::SparseIndexGraph m_graph;
};

template <typename point_t>
inline VertexGraph::VertexGraph(const Triangulation<point_t> & tri, const TriIdxSet & ignored)
{
    m_graph.clear();
    for(VerIdx i = 0; i < tri.vertices.size(); ++i){
        const IndexVertex & vertex =  tri.vertices[i];
        for(TriIdx it : vertex.triangles){
            if(ignored.count(it)) continue;
            const IndexTriangle & triangle = tri.triangles[it];
            auto vertices = triangle.VeOpVes(i);
            topology::AddEdge(i, vertices.first, m_graph);
            topology::AddEdge(i, vertices.second, m_graph);
       }
    }
}

inline void VertexGraph::AddEdge(VerIdx i, VerIdx j)
{
    topology::AddEdge(i, j, m_graph);
}

inline void VertexGraph::GetEdges(Edges & edges) const
{
    edges.clear();
    auto [begin, end] = topology::Edges(m_graph);
    for(auto iter = begin; iter != end; ++iter){
        auto u = topology::Source(*iter, m_graph);
        auto v = topology::Target(*iter, m_graph);
        edges.emplace_back(Edge(u, v));
    }
}
}//namespace tri
}//namespace geometry
}//namespace generic
#endif//GENERIC_GEOMETRY_TRI_TRIANGULATION_HPP
