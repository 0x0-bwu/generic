/**
 * @file Mesh2D.hpp
 * @author bwu
 * @brief 2d mesh flow
 * @version 0.1
 * @date 2024-09-20
 */
#pragma once
#include "TriangulationRefinement.hpp"

namespace generic::geometry{
namespace mesh2d  {

using coor_t = int64_t;
using float_t = float_type<coor_t>;
using IndexEdge = tri::IndexEdge;
using IndexEdgeList = std::list<IndexEdge>;
using Point2DContainer = std::vector<Point2D<coor_t> >;
using PolygonContainer = std::vector<Polygon2D<coor_t> >;
using Segment2DContainer = std::vector<Segment2D<coor_t> >;

inline void ExtractIntersections(const PolygonContainer & polygons, Segment2DContainer & segments)
{
    segments.clear();
    std::list<Segment2D<coor_t> > temp;
    for(const auto & polygon : polygons){
        size_t size = polygon.Size();
        for(size_t i = 0; i < size; ++i){
            size_t j = (i + 1) % size;
            temp.emplace_back(Segment2D<coor_t>(polygon[i], polygon[j]));
        }
    }
    boost::polygon::intersect_segments(segments, temp.begin(), temp.end());
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

} // namespace mesh2d
} // namespace generic::geometry


