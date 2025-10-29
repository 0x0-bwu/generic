/**
 * @file Mesh2D.hpp
 * @author bwu
 * @brief 2d mesh flow
 * @version 0.1
 * @date 2024-09-20
 */
#pragma once
#include "../TriangulationRefinement.hpp"
#include "generic/tree/QuadTreeUtilityMT.hpp"

namespace generic::geometry::mesh2d {

using coor_t = int64_t;
using float_t = float_type<coor_t>;
using IndexEdge = tri::IndexEdge;
using IndexEdgeList = std::list<IndexEdge>;
using Point2DContainer = std::vector<Point2D<coor_t> >;
using PolygonContainer = std::vector<Polygon2D<coor_t> >;
using Segment2DContainer = std::vector<Segment2D<coor_t> >;

inline void ExtractSegment(const Polygon2D<coor_t> & polygon, Segment2DContainer & segments)
{
    size_t size = polygon.Size();
    for(size_t i = 0; i < size; ++i){
        size_t j = (i + 1) % size;
        segments.emplace_back(Segment2D<coor_t>(polygon[i], polygon[j]));
    }
}

inline void ExtractSegments(const PolygonContainer & polygons, Segment2DContainer & segments)
{
    std::for_each(polygons.begin(), polygons.end(), [&](const auto & p){ ExtractSegment(p, segments); });
}

inline void ExtractIntersections(const Segment2DContainer & segments, Segment2DContainer & intersections)
{
    intersections.clear();
    boost::polygon::intersect_segments(intersections, segments.begin(), segments.end());
}

inline void ExtractTopology(const Segment2DContainer & segments, Point2DContainer & points, IndexEdgeList & edges)
{
    edges.clear();
    points.clear();
    using EdgeSet = topology::UndirectedIndexEdgeSet;
    using PointIdxMap = std::unordered_map<Point2D<coor_t>, size_t>;
    
    EdgeSet edgeSet;
    PointIdxMap pointIdxMap;

    auto getIndex = [&pointIdxMap, &points](const Point2D<coor_t> & p) mutable
    {
        if(!pointIdxMap.count(p)){
            pointIdxMap.insert(std::make_pair(p, points.size()));
            points.push_back(p);
        }
        return pointIdxMap.at(p);
    };

    points.reserve(2 * segments.size());
    for(const auto & segment : segments){
        IndexEdge e(getIndex(segment[0]), getIndex(segment[1]));
        if(edgeSet.count(e)) continue;
        edgeSet.insert(e);
        edges.emplace_back(std::move(e));
    }
}

inline void MergeClosePointsAndRemapEdge(Point2DContainer & points, IndexEdgeList & edges, coor_t tolerance)
{
    if(tolerance != 0) tri::RemoveDuplicatesAndRemapEdges(points, edges, tolerance);
}

inline void SplitOverlengthEdges(Point2DContainer & points, IndexEdgeList & edges, coor_t maxLength)
{
    if (0 >= maxLength) return;
    auto maxLenSq = maxLength * maxLength;

    auto lenSq = [&points](const IndexEdge & e) { return DistanceSq(points[e.v1()], points[e.v2()]); };
    auto split = [&points](const IndexEdge & e) mutable
    {
        size_t index = points.size();
        points.push_back((points[e.v1()] + points[e.v2()]) * 0.5);
        return std::make_pair(IndexEdge(e.v1(), index), IndexEdge(index, e.v2()));
    };

    IndexEdgeList tmp;
    while (edges.size()) {
        IndexEdge e = edges.front();
        edges.pop_front();
        if (maxLenSq < lenSq(e)) {
            auto added = split(e);
            edges.emplace_front(std::move(added.first));
            edges.emplace_front(std::move(added.second));
        }
        else tmp.emplace_back(std::move(e));
    }
    std::swap(edges, tmp);
}

inline void AddPointsFromBalancedQuadTree(const Polygon2D<coor_t> & outline, Point2DContainer & points, size_t threshold, size_t threads)
{
    struct PointExtent
    {
        Box2D<coor_t> operator()(const Point2D<coor_t> & point) const
        {
            return Box2D<coor_t>(point, point);
        }
    };
    
    std::list<Point2D<coor_t> * > objs;
    for(size_t i = 0; i < points.size(); ++i)
        objs.push_back(&points[i]);
    
    threshold = std::max(size_t(1), threshold);
    Box2D<coor_t> bbox = Extent(outline);
    using Tree = tree::QuadTree<coor_t, Point2D<coor_t>, PointExtent>;
    using Node = typename Tree::QuadNode;
    using TreeBuilder = tree::QuadTreeBuilderMT<Point2D<coor_t>, Tree>;
    Tree tree(bbox);
    TreeBuilder builder(tree, threads);
    builder.Build(objs, threshold);

    tree.Balance();

    std::list<Node * > leafNodes;
    Tree::GetAllLeafNodes(&tree, leafNodes);

    for(auto node : leafNodes){
        if(node->GetObjs().size() > 0) continue;
        const auto & box = node->GetBBox();
        Point2D<coor_t> ct = box.Center().Cast<coor_t>();
        if(Contains(outline, ct))
            points.emplace_back(std::move(ct));
    }
}

inline void TriangulatePointsAndEdges(const Point2DContainer & points, const IndexEdgeList & edges, tri::Triangulation<Point2D<coor_t> > & tri)
{
    tri.Clear();
    try {
        tri::Triangulator2D<coor_t> triangulator(tri);
        triangulator.InsertVertices(points.begin(), points.end(), [](const Point2D<coor_t> & p){ return p[0]; }, [](const Point2D<coor_t> & p){ return p[1]; });
        triangulator.InsertEdges(edges.begin(), edges.end(), [](const IndexEdge & e){ return e.v1(); }, [](const IndexEdge & e){ return e.v2(); });
        triangulator.EraseOuterTriangles();
    }
    catch (...) {
        tri.Clear();
        ThrowException("failed to generate mesh");
    }
}

template <class Refinement = tri::JonathanRefinement2D<coor_t>>
inline void TriangulationRefinement(tri::Triangulation<Point2D<coor_t> > & triangulation, float_t minAlpha, coor_t minLen, coor_t maxLen, size_t iteration)
{
    Refinement refinement(triangulation);
    refinement.SetParas(minAlpha, minLen, maxLen);
    refinement.Refine(iteration);
    refinement.ReallocateTriangulation();
}

} // namespace generic::geometry::mesh2d