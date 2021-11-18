#ifndef GENERIC_GEOMETRY_UTILITY_HPP
#define GENERIC_GEOMETRY_UTILITY_HPP
#include "BoostGeometryRegister.hpp"
#include "BoostPolygonRegister.hpp"
#include "generic/math/MathUtility.hpp"
#include "generic/common/Version.hpp"
#include "generic/math/Numbers.hpp"
#include "GeometryTraits.hpp"
#include "Trapezoid.hpp"
#include "Sphere.hpp"
#include "Plane.hpp"
#include "Line.hpp"
#include <algorithm>
#include <cassert>
#include <deque>
#include <cmath>
namespace generic {
namespace geometry{  

using generic::common::float_type;

template<typename vector_t, typename std::enable_if<traits::is_vector_t<vector_t>::value, bool>::type = true>
inline auto Inverse(const vector_t & vec) -> vector_f<vector_t>;

template<typename vector_t, typename std::enable_if<traits::is_vector_t<vector_t>::value, bool>::type = true>
inline auto SafeInverse(const vector_t & vec) -> vector_f<vector_t>;

template <typename vector_t, typename std::enable_if<traits::is_vector_t<vector_t>::value, bool>::type = true>
inline typename vector_t::coor_t DotProduct(const vector_t & a, const vector_t & b) { return a.Dot(b); }

template <typename num_type>
inline num_type CrossProduct(const Vector2D<num_type> & a, const Vector2D<num_type> & b) { return a.CrossProduct(b); }

template <typename num_type>
inline Vector3D<num_type> CrossProduct(const Vector3D<num_type> & a, const Vector3D<num_type> & b) { return a.CrossProduct(b); }

template <typename point_t>
inline typename point_t::coor_t DistanceSq(const point_t & a, const point_t & b) { return (a - b).NormSquare(); }

template <typename point_t>
inline coor_f<point_t> Distance(const point_t & a, const point_t & b) { return std::sqrt(DistanceSq(a, b)); }

template <typename num_type>
inline Polygon2D<num_type> toPolygon(const Triangle2D<num_type> & tri);

template <typename num_type>
inline Polygon2D<num_type> toPolygon(const Box2D<num_type> & box);

template <typename num_type>
inline Polygon2D<num_type> toPolygon(const Polyline2D<num_type> & polyline, num_type width);

template <typename num_type>
inline Polygon2D<num_type> InscribedPolygon(const Circle<num_type> & c, size_t div);

template <typename num_type>
inline Polygon2D<num_type> CircumscribedPolygon(const Circle<num_type> & c, size_t div);

template <typename num_type>
inline Circle<float_type<num_type> > DiametralCircle(const Segment2D<num_type> & s) { return DiametralCircle(s[0], s[1]); }

template <typename num_type>
inline Circle<float_type<num_type> > DiametralCircle(const Point2D<num_type> & p1, const Point2D<num_type> & p2);

template <typename num_type>
inline Circle<float_type<num_type> > CircumCircle(const Triangle2D<num_type> & tri) { return CircumCircle(tri[0], tri[1], tri[2]); }

template <typename num_type>
inline Circle<float_type<num_type> > CircumCircle(const Point2D<num_type> & p1, const Point2D<num_type> & p2, const Point2D<num_type> & p3);

template <typename num_type>//return coor of circle
inline Point2D<float_type<num_type> > CircumCircle(const Point2D<num_type> & p1, const Point2D<num_type> & p2, const Point2D<num_type> & p3, float_type<num_type> & radius2);

template <typename point_t, typename std::enable_if<traits::is_2d_point_t<point_t>::value, bool>::type = true>//todo point3d
inline coor_f<point_t> CircumRadius2ShortestEdgeRatioSq(const point_t & p1, const point_t & p2, const point_t & p3);

template <typename point_t, typename std::enable_if<traits::is_2d_point_t<point_t>::value, bool>::type = true>//todo point3d
inline coor_f<point_t> CircumRadius2ShortestEdgeRatio(const point_t & p1, const point_t & p2, const point_t & p3);

template <typename vector_t, typename std::enable_if<traits::is_point_t<vector_t>::value, bool>::type = true>
inline coor_f<vector_t> Angle(const vector_t & a, const vector_t & b);

template <typename point_t, typename std::enable_if<traits::is_point_t<point_t>::value, bool>::type = true>//angle of vec(p, s) and vec(p, e)
inline coor_f<point_t> Angle(const point_t & s, const point_t & p, const point_t & e) { return Angle(s - p, e - p); }

template <size_t vertex, typename triangle_t, typename std::enable_if<traits::is_triangle_t<triangle_t>::value, bool>::type = true>
inline coor_f<triangle_t> InteriorAngle(const triangle_t & t);

template <typename triangle_t, typename std::enable_if<traits::is_triangle_t<triangle_t>::value, bool>::type = true>
inline coor_f<triangle_t> InteriorAngle(const triangle_t & t, const size_t vertex);

template <typename point_t, typename std::enable_if<traits::is_point_t<point_t>::value, bool>::type = true>
inline std::array<coor_f<point_t>, 3> InteriorAngles(const point_t & p1, const point_t & p2, const point_t & p3);

template <typename num_type>
inline Point3D<float_type<num_type> > ClosestPointOnPlane(const Point3D<num_type> & p, const Plane<num_type> & plane);

template <typename num_type>
inline float_type<num_type> PointPlaneDistance(const Point3D<num_type> & p, const Plane<num_type> & plane);

template <typename point_t, typename std::enable_if<traits::is_point_t<point_t>::value, bool>::type = true>//return projection dist^2
inline coor_f<point_t> PointLineDistanceSq(const point_t & p, const point_t & a, const point_t & b);

template <typename point_t, typename segment_t,
          typename std::enable_if<traits::is_same_dim_t<point_t, segment_t>::value &&
                                  traits::is_same_coor_t<point_t, segment_t>::value, bool>::type = true>
inline point_t ClosestPointOnSegment(const point_t & p, const segment_t & s);

template <typename point_t, typename segment_t,
          typename std::enable_if<traits::is_same_dim_t<point_t, segment_t>::value &&
                                  traits::is_same_coor_t<point_t, segment_t>::value, bool>::type = true>//return dist^2
inline coor_f<point_t> PointSegmentDistanceSq(const point_t & p, const segment_t & s);

template <typename point_t, typename box_t,
          typename std::enable_if<traits::is_same_dim_t<point_t, box_t>::value &&
                                  traits::is_same_coor_t<point_t, box_t>::value, bool>::type = true>
inline point_t ClosestPointInBox(const point_t & p, const box_t & b);

template <typename point_t, typename box_t,
          typename std::enable_if<traits::is_same_dim_t<point_t, box_t>::value &&
                                  traits::is_same_coor_t<point_t, box_t>::value, bool>::type = true>//return dist^2
inline typename point_t::coor_t PointBoxDistanceSq(const point_t & p, const box_t & b);

template <typename point_t, typename segment_t, 
          typename std::enable_if<traits::is_same_dim_t<point_t, segment_t>::value, bool>::type = true>
inline PointLineLocation GetPointSegmentLocation(const point_t & p, const segment_t & seg);

template <typename point_t, typename line_t, 
          typename std::enable_if<traits::is_same_dim_t<point_t, line_t>::value, bool>::type = true>
inline PointLineLocation GetPointLineLocation(const point_t & p, const line_t & line);

template <typename point_t1, typename point_t2,
          typename std::enable_if<traits::is_same_dim_t<point_t1, point_t2>::value, bool>::type = true>
inline PointLineLocation GetPointLineLocation(const point_t1 & p, const point_t2 & v1, const point_t2 & v2);

template <typename point_t1, typename point_t2,
          typename std::enable_if<traits::is_same_dim_t<point_t1, point_t2>::value, bool>::type = true>
inline PointTriangleLocation GetPointTriangleLocation(const point_t1 & p, const point_t2 & v1, const point_t2 & v2, const point_t2 & v3);

template <typename point_t, typename triangle_t,
          typename std::enable_if<traits::is_same_dim_t<point_t, triangle_t>::value, bool>::type = true>
inline PointTriangleLocation GetPointTriangleLocation(const point_t & p, const triangle_t & tri);

template <typename segment_t1, typename segment_t2,
          typename std::enable_if<traits::is_2d_geometry_t<segment_t1>::value && 
                                  traits::is_2d_geometry_t<segment_t2>::value, bool>::type = true>
inline bool Intersects(const segment_t1 & s1, const segment_t2 & s2, bool considerTouch = true);

template <typename num_type>
inline bool Intersection(const Segment2D<num_type> & s1, const Segment2D<num_type> & s2, std::vector<Point2D<num_type> > & points);

template <typename num_type>
inline bool Intersection(const Segment2D<num_type> & s, const Line2D<num_type> & line, std::vector<Point2D<num_type> > & points);

template <typename geometry_t1, typename geometry_t2,
          typename std::enable_if<traits::is_geometry_t<geometry_t1>::value &&
                                  traits::is_geometry_t<geometry_t2>::value, bool>::type = true>
inline bool Contains(const geometry_t1 & g1, const geometry_t2 & g2, bool considerTouch = true);

template <typename point_t, typename std::enable_if<traits::is_point_t<point_t>::value, bool>::type = true>
inline point_t Extent(const point_t & p) { return p; }

template <typename segment_t, typename std::enable_if<traits::is_segment_t<segment_t>::value, bool>::type = true>
inline box_type<segment_t> Extent(const segment_t & segment);

template <typename triangle_t, typename std::enable_if<traits::is_triangle_t<triangle_t>::value, bool>::type = true>
inline box_type<triangle_t> Extent(const triangle_t & triangle);

template <typename box_t, typename std::enable_if<traits::is_box_t<box_t>::value, bool>::type = true>
inline box_t Extent(const box_t & box) { return box; }

template <typename polygon_t, typename std::enable_if<traits::is_polygon_t<polygon_t>::value ||
                                                       traits::is_polygon_with_holes_t<polygon_t>::value, bool>::type = true>
inline box_type<polygon_t> Extent(const polygon_t & polygon);

template <typename iterator,
          typename std::enable_if<traits::is_geometry_t<
          typename std::iterator_traits<iterator>::value_type>::value, bool>::type = true>
inline box_type<typename std::iterator_traits<iterator>::value_type> Extent(iterator begin, iterator end);

template <typename num_type>
inline Box2D<num_type> Extent(const Polyline2D<num_type> & polyline);

#if GENERIC_CURRENT_BOOST_LIBRARY_VER >= 165
template <typename polygon_t, template <typename, typename> class container, template <typename> class allocator = std::allocator>
inline polygon_t ConvexHull(const container<polygon_t, allocator<polygon_t> > & polygons);
#endif

}//namespace geometry
}//namespace generic

#include "Utility.ipp"
#endif //GENERIC_GEOMETRY_UTILITY_HPP
