/**
 * @file Mesh2D.hpp
 * @author bwu
 * @brief 2d mesh flow
 * @version 0.1
 * @date 2024-09-20
 */
#pragma once
#include "TriangulationRefinement.hpp"
#include "generic/tree/QuadTreeUtilityMT.hpp"

namespace generic::geometry::mesh2d {

using Coord = int64_t;
using Float = float_type<Coord>;
using Edge = tri::IndexEdge;
using Edges = std::list<Edge>;
using Point = Point2D<Coord>;
using Segment = Segment2D<Coord>;
using Polygon = Polygon2D<Coord>;
using Points = std::vector<Point>;
using Segments = std::vector<Segment>;
using Polygons = std::vector<Polygon>;

inline void ExtractSegments(const Polygon & polygon, Segments & segments)
{
    size_t size = polygon.Size();
    for(size_t i = 0; i < size; ++i){
        size_t j = (i + 1) % size;
        segments.emplace_back(Segment(polygon[i], polygon[j]));
    }
}
    
inline void ExtractSegments(const Polygons & polygons, Segments & segments)
{
    std::for_each(polygons.begin(), polygons.end(), [&](const auto & p){ ExtractSegments(p, segments); });
}

inline void ExtractIntersections(const Segments & segments, Segments & intersections)
{
    intersections.clear();
    boost::polygon::intersect_segments(intersections, segments.begin(), segments.end());
}

inline void ExtractTopology(const Segments & segments, Points & points, Edges & edges)
{
    edges.clear();
    points.clear();
    using EdgeSet = topology::UndirectedIndexEdgeSet;
    using PointIdxMap = std::unordered_map<Point, size_t>;
    
    EdgeSet edgeSet;
    PointIdxMap pointIdxMap;

    auto getIndex = [&pointIdxMap, &points](const Point & p) mutable
    {
        if (not pointIdxMap.count(p)) {
            pointIdxMap.insert(std::make_pair(p, points.size()));
            points.push_back(p);
        }
        return pointIdxMap.at(p);
    };

    points.reserve(2 * segments.size());
    for(const auto & segment : segments){
        Edge e(getIndex(segment[0]), getIndex(segment[1]));
        if(edgeSet.count(e)) continue;
        edgeSet.insert(e);
        edges.emplace_back(std::move(e));
    }
    points.shrink_to_fit();
}

inline void MergeClosePointsAndRemapEdge(Points & points, Edges & edges, Coord tolerance)
{
    if (tolerance != 0) tri::RemoveDuplicatesAndRemapEdges(points, edges, tolerance);
}

inline void SplitOverlengthEdges(Points & points, Edges & edges, Coord maxLength)
{
    if (0 >= maxLength) return;
    auto maxLenSq = maxLength * maxLength;

    auto lenSq = [&points](const Edge & e) { return DistanceSq(points[e.v1()], points[e.v2()]); };
    auto split = [&points](const Edge & e) mutable
    {
        size_t index = points.size();
        points.push_back((points[e.v1()] + points[e.v2()]) * 0.5);
        return std::make_pair(Edge(e.v1(), index), Edge(index, e.v2()));
    };

    Edges tmp;
    while (edges.size()) {
        auto e = edges.front();
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

inline void AddPointsFromBalancedQuadTree(const Polygon & outline, Points & points, size_t threshold, size_t threads)
{
    struct PointExtent
    {
        Box2D<Coord> operator()(const Point & point) const
        {
            return Box2D<Coord>(point, point);
        }
    };
    
    std::list<Point * > objs;
    for(size_t i = 0; i < points.size(); ++i)
        objs.push_back(&points[i]);
    
    threshold = std::max(size_t(1), threshold);
    Box2D<Coord> bbox = Extent(outline);
    using Tree = tree::QuadTree<Coord, Point, PointExtent>;
    using Node = typename Tree::QuadNode;
    using TreeBuilder = tree::QuadTreeBuilderMT<Point, Tree>;
    Tree tree(bbox);
    TreeBuilder builder(tree, threads);
    builder.Build(objs, threshold);

    tree.Balance();

    std::list<Node * > leafNodes;
    Tree::GetAllLeafNodes(&tree, leafNodes);

    for(auto node : leafNodes){
        if(node->GetObjs().size() > 0) continue;
        const auto & box = node->GetBBox();
        auto ct = box.Center().Cast<Coord>();
        if(Contains(outline, ct))
            points.emplace_back(std::move(ct));
    }
}

inline void TriangulatePointsAndEdges(const Points & points, const Edges & edges, tri::Triangulation<Point> & tri)
{
    tri.Clear();
    try {
        tri::Triangulator2D<Coord> triangulator(tri);
        triangulator.InsertVertices(points.begin(), points.end(), [](const Point & p){ return p[0]; }, [](const Point & p){ return p[1]; });
        triangulator.InsertEdges(edges.begin(), edges.end(), [](const Edge & e){ return e.v1(); }, [](const Edge & e){ return e.v2(); });
        triangulator.EraseOuterTriangles();
    }
    catch (...) {
        tri.Clear();
        ThrowException("failed to generate mesh");
    }
}

template <class Refinement = tri::JonathanRefinement2D<Coord>>
inline void TriangulationRefinement(tri::Triangulation<Point> & triangulation, float_t minAlpha, Coord minLen, Coord maxLen, size_t iteration)
{
    Refinement refinement(triangulation);
    refinement.SetParas(minAlpha, minLen, maxLen);
    refinement.Refine(iteration);
    refinement.ReallocateTriangulation();
}

} // namespace generic::geometry::mesh2d
