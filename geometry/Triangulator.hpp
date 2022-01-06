#ifndef GENERIC_GEOMETRY_TRI_TRIANGULATOR_HPP
#define GENERIC_GEOMETRY_TRI_TRIANGULATOR_HPP
#include <boost/geometry/index/rtree.hpp>
#include <boost/geometry/geometry.hpp>
#include "TriangulationOperator.hpp"
#include <unordered_map>
#include <unordered_set>
#include <limits>
namespace generic  {
namespace geometry {
namespace tri      {
using generic::common::float_type;

template <typename num_type>
class RuppertRefinement2D;

template <typename num_type>
class ChewSecondRefinement2D;

template <typename num_type>
class Triangulator2D
{
public:
    using Depth = unsigned short;
    using Point = Point2D<num_type>;
    using Box = Box2D<num_type>;
    using Edge = IndexEdge;
    using Vertex = IndexVertex;
    using Triangle = IndexTriangle;

    using EdgeSet = typename Triangulation<Point>::EdgeSet;
    using TriIdxMap = std::unordered_map<TriIdx, TriIdx>;

    Triangulation<Point> & tri;
    TriangulationOperator<Point> op;
    Triangulator2D(Triangulation<Point> & t) : tri(t), op(t) {}

    template <typename VertexIterator,
              typename CoorXGetter,
              typename CoorYGetter>
    void InsertVertices(VertexIterator begin, VertexIterator end, CoorXGetter xGetter, CoorYGetter yGetter);

    template <typename EdgeIterator,
              typename sVerIdxGetter,
              typename eVerIdxGetter>
    void InsertEdges(EdgeIterator begin, EdgeIterator end, sVerIdxGetter sGetter, eVerIdxGetter eGetter);

    void EraseSuperTriangle();
    void EraseOuterTriangles();
    void EraseOuterTrianglesAndHoles();

    void Clear();
private:
    template <typename VertexIterator,
                typename CoorXGetter,
                typename CoorYGetter>
    Box Extent(VertexIterator begin, VertexIterator end, CoorXGetter xGetter, CoorYGetter yGetter);

    void AddSuperTriangle(const Box & box);
    
    void InsertOneEdge(const Edge & edge);
    void InsertOneVertex(const Point & pos);
    
    VerIdx InsertPointInTriangle(const Point & pos, TriIdx i0, std::stack<TriIdx> & out);
    VerIdx InsertPointOnSharedEdge(const Point & pos, TriIdx i1, TriIdx i2, std::stack<TriIdx> & out);
    VerIdx InsertPointOnBoundaryEdge(const Point & pos, TriIdx it, const Edge & e, std::stack<TriIdx> & out);

    void Delaunay(std::stack<TriIdx> & triangls, VerIdx iv);
    bool isNeedFlip(VerIdx iv, TriIdx it, TriIdx itOp) const;
    bool isInCircumCircle(const Point & pos, const Point & v1, const Point & v2, const Point & v3) const;

    std::tuple<TriIdx, VerIdx, VerIdx> IntersectedTriangle(VerIdx iA,
                                                            const TriIdxSet & candidates,
                                                            const Point & a, const Point & b) const;      
    void RemoveSuperTriangleVertices();
    template <typename TriIdxIterator>
    void RemoveTriangles(TriIdxIterator begin, TriIdxIterator end);

    TriIdx TriangulatePseudoPolygon(VerIdx ia, VerIdx ib, const std::vector<VerIdx> & points);
    TriIdx PesudoPolyOuterTriangle(VerIdx ia, VerIdx ib) const;
    VerIdx FindDelaunayPoint(VerIdx ia, VerIdx ib, const std::vector<VerIdx> & points) const;
    std::pair<std::vector<VerIdx>, std::vector<VerIdx> >
    SplitPsedoPolygon(VerIdx vi, const std::vector<VerIdx> & points);

    TriIdxSet Grow2Boundary(std::stack<TriIdx> & seeds) const;
    std::vector<Depth> CalculateTriangleDepths(const std::vector<Vertex> & vertices,
                                                const TriangleVec & triangles,
                                                const EdgeSet & fixedEdges);
    TriIdxSet PeelLayer(std::stack<TriIdx> seeds,
                        const TriangleVec & triangles,
                        const EdgeSet & fixedEdges,
                        Depth depth,
                        std::vector<Depth> & triDepths);
    
    void RebuildRTree();
};

template <typename num_type>
template <typename VertexIterator, typename CoorXGetter, typename CoorYGetter>
inline void Triangulator2D<num_type>::InsertVertices(VertexIterator begin, VertexIterator end, CoorXGetter xGetter, CoorYGetter yGetter)
{
    if(tri.vertices.empty()){
        auto bounds = Extent(begin, end, xGetter, yGetter);
        AddSuperTriangle(bounds);
    }
    
    tri.vertices.reserve(tri.vertices.size() + std::distance(begin, end));
    for(auto iter = begin; iter != end; ++iter){
        InsertOneVertex(Point(xGetter(*iter), yGetter(*iter)));
    }
}

template <typename num_type>
template <typename EdgeIterator, typename sVerIdxGetter, typename eVerIdxGetter>
inline void Triangulator2D<num_type>::InsertEdges(EdgeIterator begin, EdgeIterator end, sVerIdxGetter sGetter, eVerIdxGetter eGetter)
{
    for(auto iter = begin; iter != end; ++iter){
        InsertOneEdge(Edge(sGetter(*iter) + 3, eGetter(*iter) + 3));
    }
    op.ReallocateTriangulation();
}

template <typename num_type>
inline void Triangulator2D<num_type>::EraseSuperTriangle()
{
    TriIdxVec toErase;
    for(TriIdx it = 0; it < tri.triangles.size(); ++it){
        auto & t = tri.triangles[it];
        if(t.vertices[0] < 3 || t.vertices[1] < 3 || t.vertices[2] < 3)
            toErase.push_back(it);
    }
    RemoveTriangles(toErase.begin(), toErase.end());
    RemoveSuperTriangleVertices();
}

template <typename num_type>
inline void Triangulator2D<num_type>::EraseOuterTriangles()
{
    std::stack<TriIdx> seed(std::deque<TriIdx>(1, *(tri.vertices[0].triangles.begin())));
    TriIdxSet toErase = Grow2Boundary(seed);
    RemoveTriangles(toErase.begin(), toErase.end());
    RemoveSuperTriangleVertices();
}

template <typename num_type>
inline void Triangulator2D<num_type>::EraseOuterTrianglesAndHoles()
{
    auto triDepths = CalculateTriangleDepths(tri.vertices, tri.triangles, tri.fixedEdges);

    TriIdxVec toErase;
    toErase.reserve(tri.triangles.size());
    for(size_t it = 0; it != tri.triangles.size(); ++it){
        if(triDepths[it] % 2 == 0)
            toErase.push_back(it);
    }
    RemoveTriangles(toErase.begin(), toErase.end());
    RemoveSuperTriangleVertices();
}

template <typename num_type>
template <typename VertexIterator, typename CoorXGetter, typename CoorYGetter>
inline Box2D<num_type> Triangulator2D<num_type>::Extent(VertexIterator begin, VertexIterator end, CoorXGetter xGetter, CoorYGetter yGetter)
{
    Box box;
    for(auto iter = begin; iter != end; ++iter)
        box |= Point(xGetter(*iter), yGetter(*iter));
    return box;
}

template <typename num_type>
inline void Triangulator2D<num_type>::AddSuperTriangle(const Box & box)
{
    auto center = box.Center();
    auto l = box.Length();
    auto w = box.Width();
    l *= l; w *= w;
    auto r = 0.5 * 1.1 * std::sqrt(float_type<num_type>(l + w));
    auto sx = r * std::sqrt(3);
    VerIdx iv1 = op.AddOneVertex(Point(center[0] - sx, center[1] - r), {0});
    VerIdx iv2 = op.AddOneVertex(Point(center[0] + sx, center[1] - r), {0});
    VerIdx iv3 = op.AddOneVertex(Point(center[0], center[1] + 2 * r), {0});
    Triangle super{{iv1, iv2, iv3}, {noNeighbor, noNeighbor, noNeighbor}};
    op.AddOneTriangle(super);
}

template <typename num_type>
inline void Triangulator2D<num_type>::InsertOneVertex(const Point & pos)
{
    VerIdx iv;
    std::stack<TriIdx> tris;
    auto [it1, it2] = op.TraversalTriangleAt(pos);
    if(noNeighbor == it2) iv = op.InsertPointInTriangle(pos, it1, tris);
    else iv = op.InsertPointOnSharedEdge(pos, it1, it2, tris);
    Delaunay(tris, iv);
}

template <typename num_type>
inline void Triangulator2D<num_type>::InsertOneEdge(const Edge & edge)
{
    auto is = edge.v1();
    auto ie = edge.v2();
    if(is == ie) return;

    const Vertex & vs = tri.vertices[is];
    const Vertex & ve = tri.vertices[ie];

    if(Vertex::isShareEdge(vs, ve)){
        tri.fixedEdges.insert(Edge(is, ie));
        return;
    }

    TriIdx it;
    VerIdx ivLeft, ivRight;
    std::tie(it, ivLeft, ivRight) = IntersectedTriangle(is, vs.triangles, tri.VertexPoint(vs), tri.VertexPoint(ve));
    if(it == noNeighbor){
        tri.fixedEdges.insert(Edge(is, ivLeft));
        return InsertOneEdge(Edge(ivLeft, ie));
    }

    std::vector<TriIdx> intersected(1, it);
    std::vector<VerIdx> ptsLeft(1, ivLeft);
    std::vector<VerIdx> ptsRight(1, ivRight);
    VerIdx iv = is;
    Triangle t = tri.triangles[it];
    const auto & vrts = t.vertices;
    while(std::find(vrts.begin(), vrts.end(), ie) == vrts.end()){
        
        auto itOp = t.VeOpTn(iv);
        const Triangle & tOp = tri.triangles[itOp];
        auto ivOp = tOp.TnOpVe(it);
        auto vOp = tri.vertices[ivOp];

        intersected.push_back(itOp);
        it = itOp;
        t = tri.triangles[it];

        PointLineLocation loc = GetPointLineLocation(tri.VertexPoint(vOp), tri.VertexPoint(vs), tri.VertexPoint(ve));
        if(loc == PointLineLocation::Left){
            ptsLeft.push_back(ivOp);
            iv = ivLeft;
            ivLeft = ivOp;
        }
        else if(loc == PointLineLocation::Right){
            ptsRight.push_back(ivOp);
            iv = ivRight;
            ivRight = ivOp;
        }
        else { ie = ivOp; }
    }

    auto iter = intersected.begin();
    for(; iter != intersected.end(); ++iter)
        op.RemoveOneTriangle(*iter);
    
    auto itLeft = TriangulatePseudoPolygon(is, ie, ptsLeft);
    std::reverse(ptsRight.begin(), ptsRight.end());
    auto itRight = TriangulatePseudoPolygon(ie, is, ptsRight);
    
    tri.ChangeNeighbor(itLeft, noNeighbor, itRight);
    tri.ChangeNeighbor(itRight, noNeighbor, itLeft);

    tri.fixedEdges.insert(Edge(is, ie));
    if(ie != edge.v2())
        return InsertOneEdge(Edge(ie, edge.v2()));
}

template <typename num_type>
inline void Triangulator2D<num_type>::Delaunay(std::stack<TriIdx> & triangles, VerIdx iv)
{
    while(!triangles.empty()){
        TriIdx it = triangles.top();
        triangles.pop();

        const auto & triangle = tri.triangles[it];
        TriIdx itOp = triangle.VeOpTn(iv);
        if(itOp == noNeighbor) continue;
        if(isNeedFlip(iv, it, itOp)){
            tri.FlipEdge(it, itOp);
            triangles.push(it);
            triangles.push(itOp);
        }
    }
}

template <typename num_type>
inline bool Triangulator2D<num_type>::isNeedFlip(VerIdx iv, TriIdx it, TriIdx itOp) const
{
    const Triangle & triOp = tri.triangles[itOp];
    size_t i = triOp.TnOpiVe(it);
    VerIdx ivOp = triOp.vertices[i];
    if(iv < 3 && ivOp < 3) return false;

    VerIdx ivCw  = triOp.vertices[Triangle:: cw(i)];
    VerIdx ivCcw = triOp.vertices[Triangle::ccw(i)];
    const Point & p0 = tri.VerIdxPoint(iv);
    const Point & p1 = tri.VerIdxPoint(ivCw);
    const Point & p2 = tri.VerIdxPoint(ivOp);
    const Point & p3 = tri.VerIdxPoint(ivCcw);
    if(ivCw < 3) return GetPointLineLocation(p1, p2, p3) == GetPointLineLocation(p0, p2, p3);
    else if(ivCcw < 3) return GetPointLineLocation(p3, p1, p2) == GetPointLineLocation(p0, p1, p2);
    else return isInCircumCircle(p0, p1, p2, p3);
}

template <typename num_type>
inline bool Triangulator2D<num_type>::isInCircumCircle(const Point & p, const Point & p1, const Point & p2, const Point & p3) const
{
    return geometry::isInCircumCircle(p1, p2, p3, p, false);
}

template <typename num_type>
inline std::tuple<TriIdx, VerIdx, VerIdx>
Triangulator2D<num_type>::IntersectedTriangle(VerIdx iA,
                                                const TriIdxSet & candidates,
                                                const Point & a, const Point & b) const
{
    for(TriIdx it : candidates){
        const Triangle & triangle = tri.triangles[it];
        auto i = triangle.iVe(iA);
        auto iv1 = triangle.vertices[Triangle::cw(i)];
        auto iv2 = triangle.vertices[Triangle::ccw(i)];
        auto locP1 = GetPointLineLocation(tri.VerIdxPoint(iv1), a, b);
        auto locP2 = GetPointLineLocation(tri.VerIdxPoint(iv2), a, b);
        if(locP2 == PointLineLocation::Right){
            if(locP1 == PointLineLocation::OnLine)
                return std::make_tuple(noNeighbor, iv1, iv2);
            else if(locP1 == PointLineLocation::Left)
                return std::make_tuple(it, iv1, iv2);
        }
    }
    throw std::runtime_error("Could not find vertex triangle intersected by "
                             "edge. Note: can be caused by duplicate points.");
}

template <typename num_type>
inline void Triangulator2D<num_type>::RemoveSuperTriangleVertices()
{
    for(VerIdx iv = 0; iv < 3; ++iv)
        op.RemoveOneVertex(iv);
    op.ReallocateTriangulation();
}

template <typename num_type>
template <typename TriIdxIterator>
void Triangulator2D<num_type>::RemoveTriangles(TriIdxIterator begin, TriIdxIterator end)
{
    for(auto iter = begin; iter != end; ++iter)
        op.RemoveOneTriangle(*iter);
    op.ReallocateTriangulation();
}

template <typename num_type>
inline TriIdx Triangulator2D<num_type>::TriangulatePseudoPolygon(VerIdx ia, VerIdx ib, const std::vector<VerIdx> & points)
{
    if(points.empty())
        return PesudoPolyOuterTriangle(ia, ib);
    
    VerIdx ic = FindDelaunayPoint(ia, ib, points);
    auto split = SplitPsedoPolygon(ic, points);
    auto it2 = TriangulatePseudoPolygon(ic, ib, split.second);
    auto it1 = TriangulatePseudoPolygon(ia, ic, split.first);
    Triangle t = {{ia, ib, ic}, {noNeighbor, it2, it1}};
    TriIdx it = op.AddOneTriangle(t);
    if(it1 != noNeighbor){
        if(split.first.empty())
            tri.ChangeNeighbor(it1, ia, ic, it);
        else
            tri.triangles[it1].neighbors[0] = it;
    }
    if(it2 != noNeighbor){
        if(split.second.empty())
            tri.ChangeNeighbor(it2, ic, ib, it);
        else
            tri.triangles[it2].neighbors[0] = it;
    }
    tri.AddAdjacentTriangle(ia, it);
    tri.AddAdjacentTriangle(ib, it);
    tri.AddAdjacentTriangle(ic, it);
    return it;
}

template <typename num_type>
inline TriIdx Triangulator2D<num_type>::PesudoPolyOuterTriangle(VerIdx ia, VerIdx ib) const
{
    const auto & aTris = tri.vertices[ia].triangles;
    const auto & bTris = tri.vertices[ib].triangles;
    for(auto iter = aTris.begin(); iter != aTris.end(); ++iter){
        if(std::find(bTris.begin(), bTris.end(), *iter) != bTris.end())
            return *iter;
    }
    return noNeighbor;
}

template <typename num_type>
inline VerIdx Triangulator2D<num_type>::FindDelaunayPoint(VerIdx ia, VerIdx ib, const std::vector<VerIdx> & vs) const
{
    assert(!vs.empty());
    const auto & a = tri.VerIdxPoint(ia);
    const auto & b = tri.VerIdxPoint(ib);
    VerIdx ic = vs.front();
    auto c = tri.VerIdxPoint(ic);
    for(auto iter = vs.begin(); iter != vs.end(); ++iter){
        const auto & v = tri.VerIdxPoint(*iter);
        if(!isInCircumCircle(v, a, b, c)) continue;
        ic = *iter;
        c = tri.VerIdxPoint(ic);
    }
    return ic;
}

template <typename num_type>
inline std::pair<std::vector<VerIdx>, std::vector<VerIdx> >
Triangulator2D<num_type>::SplitPsedoPolygon(VerIdx vi, const std::vector<VerIdx> & points)
{
    std::pair<std::vector<VerIdx>, std::vector<VerIdx> > out;
    auto iter = points.begin();
    for(; vi != *iter; ++iter)
        out.first.push_back(*iter);
    for(++iter; iter != points.end(); ++iter)
        out.second.push_back(*iter);
    return out;
}

template <typename num_type>
inline TriIdxSet Triangulator2D<num_type>::Grow2Boundary(std::stack<TriIdx> & seeds) const
{
    TriIdxSet traversed;
    while(!seeds.empty()){
        auto it = seeds.top();
        seeds.pop();
        traversed.insert(it);
        const auto & triange = tri.triangles[it];
        for(auto i = 0; i < 3; ++i){
            Edge opEdge = Edge(
                        triange.vertices[Triangle::ccw(i)],
                        triange.vertices[Triangle::cw(i)]);
            if(tri.fixedEdges.count(opEdge)) continue;
            TriIdx in = triange.neighbors[Triangle::iVeOpiTn(i)];
            if(in != noNeighbor && !traversed.count(in))
                seeds.push(in);
        }
    }
    return traversed;
}

template <typename num_type>
inline std::vector<typename Triangulator2D<num_type>::Depth>
Triangulator2D<num_type>::CalculateTriangleDepths(const std::vector<Vertex> & vertices,
                                                    const TriangleVec & triangles,
                                                    const EdgeSet & fixedEdges)
{
    std::vector<Depth> triDepths(triangles.size(), std::numeric_limits<Depth>::max());

    using TriDeque = std::deque<TriIdx>;
    using TriStack = std::stack<TriIdx>;
    TriStack seeds(TriDeque(1, *(tri.vertices[0].triangles.begin())));
    Depth layerDepth = 0;
    do{
        auto newSeeds = PeelLayer(seeds, triangles, fixedEdges, layerDepth++, triDepths);
        seeds = TriStack(TriDeque(newSeeds.begin(), newSeeds.end()));
    } while(!seeds.empty());

    return triDepths;
}

template <typename num_type>
inline TriIdxSet Triangulator2D<num_type>::PeelLayer(std::stack<TriIdx> seeds,
                                        const TriangleVec & triangles,
                                        const EdgeSet & fixedEdges,
                                        Depth depth,
                                        std::vector<Depth> & triDepths)
{
    TriIdxSet behindBoundary;
    while(!seeds.empty()){
        auto it = seeds.top();
        seeds.pop();
        triDepths[it] = depth;
        behindBoundary.erase(it);
        const auto & triangle = triangles[it];
        for(auto i = 0; i < 3; ++i){
            Edge opEdge(triangle.vertices[Triangle::ccw(i)],
                        triangle.vertices[Triangle:: cw(i)]);
            TriIdx in = triangle.neighbors[Triangle::iVeOpiTn(i)];
            if(in == noNeighbor || triDepths[in] <= depth) continue;
            if(fixedEdges.count(opEdge)){
                behindBoundary.insert(in);
                continue;
            }
            seeds.push(in);
        }
    }
    return behindBoundary;
}

template <typename num_type>
inline void Triangulator2D<num_type>::Clear()
{
    op.Clear();
}

template <typename num_type>
class ConvexTriangulator2D
{
    struct DoublyLinkedIndex { VerIdx prev, curr, next; };
public:
    using float_t = float_type<num_type>;
    using Point = Point2D<num_type>;
    using Edge = IndexEdge;
    using Vertex = IndexVertex;
    using Triangle = IndexTriangle;
    using Utility = TriangulationUtility<Point>;
    using EdgeSet = typename Triangulation<Point>::EdgeSet;
    using TriIdxMap = std::unordered_map<TriIdx, TriIdx>;
    

    Triangulation<Point> & tri;
    TriangulationOperator<Point> op;
    ConvexTriangulator2D(Triangulation<Point> & t) : tri(t), op(t) { op.Clear(); }

    template <typename VertexIterator,
              typename CoorXGetter,
              typename CoorYGetter>
    void Triangulate(VertexIterator begin, VertexIterator end, CoorXGetter xGetter, CoorYGetter yGetter)
    {
        const size_t size = std::distance(begin, end);
        assert(size >= 3);
        
        size_t i = 0;
        std::vector<DoublyLinkedIndex> vertices(size);
        for(auto iter = begin; iter != end; ++iter){
            vertices[i++].curr = op.AddOneVertex(Point(xGetter(*iter), yGetter(*iter)));
        }

        for(size_t i = 0; i < size; ++i){
            vertices[i].prev = vertices[(i + size - 1) % size].curr;
            vertices[i].next = vertices[(i + 1) % size].curr;
        }
        
        for(size_t i = size - 1; i > 3; --i){
            vertices[i - 1].next = vertices[i].next;
        }

        TriIdx it = AddTriangle(vertices[0].curr, vertices[1].curr, vertices[2].curr);
        for(size_t i = 3; i < size; ++i)
            ConvexInsertVertex(vertices[i].curr, Edge(vertices[i].prev, vertices[i].next));
        
        op.ReallocateTriangulation();

    }

private:

    TriIdx AddTriangle(VerIdx iv1, VerIdx iv2, VerIdx iv3)
    {
        return op.AddOneTriangle(iv1, iv2, iv3);
    }

    void ConvexInsertVertex(VerIdx iv, const Edge & e)
    {
        TriIdx it;
        std::tie(it, std::ignore) = tri.GetTriangles(e);
        VerIdx iu = e.v1(), iw = e.v2();
        if(noNeighbor == it){
            AddTriangle(iv, iu, iw);
            return;
        }
        VerIdx ix = tri.triangles[it].EgOpVe(e);
        if(Utility::isInCircumCircle(tri, it, tri.VerIdxPoint(iv))){
            op.RemoveOneTriangle(it);
            ConvexInsertVertex(iv, Edge(iu, ix));
            ConvexInsertVertex(iv, Edge(ix, iw));
        }
        else AddTriangle(iv, iu, iw);
    }
};

struct DuplicatesInfo
{
    std::vector<size_t> mapping;
    std::unordered_set<size_t> duplicates;
};

template <typename num_type>
inline DuplicatesInfo FindDuplicates(const std::vector<Point2D<num_type> > & points, num_type tolerance)
{
    using Grid = std::pair<size_t, size_t>;
    struct GridHash
    {
        size_t operator() (const Grid & grid) const noexcept
        {
            size_t seed(0);
            boost::hash_combine(seed, grid.first);
            boost::hash_combine(seed, grid.second);
            return seed;
        }
    };

    struct GridCmp
    {
        bool operator() (const Grid & g1, const Grid & g2) const noexcept
        {
            return g1 == g2;
        }
    };

    using GridIdxMap = std::unordered_map<Grid, size_t, GridHash, GridCmp>;
    const size_t shift = std::numeric_limits<size_t>::max() / 2;
    auto toGrid = [tolerance, shift](const Point2D<num_type> & p)
    {
        size_t x = static_cast<size_t>(p[0] / tolerance + shift);
        size_t y = static_cast<size_t>(p[1] / tolerance + shift);
        return std::make_pair(x, y);
    };

    GridIdxMap gridIdxMap;
    DuplicatesInfo dup { std::vector<size_t>(points.size()), std::unordered_set<size_t>() };
    for(size_t iIn = 0, iOut = iIn; iIn < points.size(); ++iIn){
        bool isUnique;
        typename GridIdxMap::const_iterator iter;
        std::tie(iter, isUnique) = gridIdxMap.insert(std::make_pair(toGrid(points[iIn]), iOut));
        if(isUnique){
            dup.mapping[iIn] = iOut++;
            continue;
        }
        dup.mapping[iIn] = iter->second;
        dup.duplicates.insert(iIn);
    }
    return dup;
}

inline void RemapEdges(std::list<IndexEdge> & edges, const std::vector<size_t> & mapping)
{
    for(auto iter = edges.begin(); iter != edges.end();){
        iter->SetVertices(mapping[iter->v1()], mapping[iter->v2()]);
        if(iter->v1() == iter->v2())
            iter = edges.erase(iter);
        else ++iter;
    }
}

template <typename num_type>
inline void RemoveDuplicates(std::vector<Point2D<num_type> > & points, const std::unordered_set<size_t> & duplicates)
{
    for(size_t i = 0, iNew = i; i < points.size(); ++i){
        if(duplicates.count(i)) continue;
        points[iNew] = points[i];
        iNew++;
    }
    points.erase(points.end() - duplicates.size(), points.end());
}

template <typename num_type>
inline DuplicatesInfo RemoveDuplicatesAndRemapEdges(std::vector<Point2D<num_type> > & points, std::list<IndexEdge> & edges, num_type tolerance)
{
    DuplicatesInfo dup = FindDuplicates(points, tolerance);
    RemoveDuplicates(points, dup.duplicates);
    RemapEdges(edges, dup.mapping);
    return dup;
}

}//namespace tri
}//namespace geometry
}//namespace generic
#endif//GENERIC_GEOMETRY_TRI_TRIANGULATOR_HPP
