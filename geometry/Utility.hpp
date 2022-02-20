/**
 * @file Utility.hpp
 * @author bwu
 * @brief Geometry related utilities
 * @version 0.1
 * @date 2022-02-22
 */
#ifndef GENERIC_GEOMETRY_UTILITY_HPP
#define GENERIC_GEOMETRY_UTILITY_HPP
#include "generic/math/MathUtility.hpp"
#include "generic/common/Version.hpp"
#include "generic/math/Numbers.hpp"
#include "BoostGeometryRegister.hpp"
#include "BoostPolygonRegister.hpp"
#include "GeometryTraits.hpp"
#include "Predicates.hpp"
#include "Trapezoid.hpp"
#include "Curves.hpp"
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

/**
 * @brief inverse of vector, will get +-INF or NAN if zero input
 * @tparam vector_t vector type, could be Vector2D, Vector3D or VectorN
 * @param[in] vec input vector
 * @return vector_f<vector_t> inversed vector with floating points type
 */
template<typename vector_t, typename std::enable_if<traits::is_vector_t<vector_t>::value, bool>::type = true>
inline auto Inverse(const vector_t & vec) -> vector_f<vector_t>;

/**
 * @brief saves inverse of vector, will get 1/epsilon if zero input
 * @tparam vector_t vector type, could be Vector2D, Vector3D or VectorN
 * @param[in] vec input vector
 * @return vector_f<vector_t> inversed vector with floating points type
 */
template<typename vector_t, typename std::enable_if<traits::is_vector_t<vector_t>::value, bool>::type = true>
inline auto SafeInverse(const vector_t & vec) -> vector_f<vector_t>;

/**
 * @brief dot product of two vectors
 * @tparam vector_t vector type, could be Vector2D, Vector3D or VectorN
 * @param[in] a one of the input vectors
 * @param[in] b one of the input vectors
 * @return dot product result
 */
template <typename vector_t, typename std::enable_if<traits::is_vector_t<vector_t>::value, bool>::type = true>
inline typename vector_t::coor_t DotProduct(const vector_t & a, const vector_t & b) { return a.Dot(b); }

/**
 * @brief cross product of two Vector2D
 * @param[in] a vector2d a
 * @param[in] b vector2d b
 * @return scalar result of a x b
 */
template <typename num_type>
inline num_type CrossProduct(const Vector2D<num_type> & a, const Vector2D<num_type> & b) { return a.CrossProduct(b); }

/**
 * @brief cross product of two Vector3D
 * @param[in] a vector3d a
 * @param[in] b vector3d b
 * @return vector3d result of a x b
 */
template <typename num_type>
inline Vector3D<num_type> CrossProduct(const Vector3D<num_type> & a, const Vector3D<num_type> & b) { return a.CrossProduct(b); }

/**
 * @brief squared distance of two points
 * @tparam point_t point type, could be Point2D or Point3D
 * @param[in] a one of the input points
 * @param[in] b one of the input points
 * @return square of a-b 's euclidean distance
 */
template <typename point_t>
inline typename point_t::coor_t DistanceSq(const point_t & a, const point_t & b) { return (a - b).NormSquare(); }

/**
 * @brief distance of two points
 * @tparam point_t point type, could be Point2D or Point3D
 * @param[in] a one of the input points
 * @param[in] b one of the input points
 * @return floating points type of a-b 's euclidean distance
 */
template <typename point_t>
inline coor_f<point_t> Distance(const point_t & a, const point_t & b) { return std::sqrt(DistanceSq(a, b)); }

///@brief convert an arc3 to arc
template<typename num_type>
inline Arc<num_type> toArc(const Arc3<num_type> & arc3);

/**
 * @brief converts an arc to polyline
 * @param[in] arc the input arc
 * @param[in] div circle divide number
 * @return the converted polyline
 */
template <typename num_type>
inline Polyline2D<num_type> toPolyline(const Arc<num_type> & arc, size_t div);

/**
 * @brief converts an arc3 to polyline
 * @param[in] arc3 the input arc3
 * @param[in] div circle divide number
 * @return the converted polyline
 */
template <typename num_type>
inline Polyline2D<num_type> toPolyline(const Arc3<num_type> & arc3, size_t div);

///@brief converts triangle2d to polygon2d
template <typename num_type>
inline Polygon2D<num_type> toPolygon(const Triangle2D<num_type> & tri);

///@brief converts box2d to polygon2d
template <typename num_type>
inline Polygon2D<num_type> toPolygon(const Box2D<num_type> & box);

///@brief converts polyline2d with width's contour to polygon2d
template <typename num_type>
inline Polygon2D<num_type> toPolygon(const Polyline2D<num_type> & polyline, num_type width);

/**
 * @brief gets inscribed polygon of circle
 * @param[in] c input circle
 * @param[in] div number of sides of the inscribed polygon
 * @return inscribed polygon of circle
 */
template <typename num_type>
inline Polygon2D<num_type> InscribedPolygon(const Circle<num_type> & c, size_t div);

/**
 * @brief gets circumscribed polygon of circle
 * @param[in] c input circle 
 * @param[in] div number of sides of the circumscribed polygon
 * @return circumscribed polygon of circle
 */
template <typename num_type>
inline Polygon2D<num_type> CircumscribedPolygon(const Circle<num_type> & c, size_t div);

///@brief gets diametral circle of a segment2d
template <typename num_type>
inline Circle<float_type<num_type> > DiametralCircle(const Segment2D<num_type> & s) { return DiametralCircle(s[0], s[1]); }

///@brief gets diametral circle of two points that form a segment2d
template <typename num_type>
inline Circle<float_type<num_type> > DiametralCircle(const Point2D<num_type> & p1, const Point2D<num_type> & p2);

///@brief gets circumcircle of a triangle2d
template <typename num_type>
inline Circle<float_type<num_type> > CircumCircle(const Triangle2D<num_type> & tri) { return CircumCircle(tri[0], tri[1], tri[2]); }

///@brief gets circumcircle of three points that form a triangle2d
template <typename num_type>
inline Circle<float_type<num_type> > CircumCircle(const Point2D<num_type> & p1, const Point2D<num_type> & p2, const Point2D<num_type> & p3);

/**
 * @brief gets circumcircle of three points that form a triangle2d
 * @param[in] p1 one of three points that form a triangle2d
 * @param[in] p2 one of three points that form a triangle2d
 * @param[in] p3 one of three points that form a triangle2d
 * @param[out] radius2 squared radius of the circumcircle 
 * @return floating points type 2d coordinate of the circumcircle origin
 */
template <typename num_type>
inline Point2D<float_type<num_type> > CircumCircle(const Point2D<num_type> & p1, const Point2D<num_type> & p2, const Point2D<num_type> & p3, float_type<num_type> & radius2);

/**
 * @brief tests if a point is inside the circumcirle form by three points
 * @param[in] p1 one of three points that form a triangle2d
 * @param[in] p2 one of three points that form a triangle2d
 * @param[in] p3 one of three points that form a triangle2d
 * @param[in] p4 testing point
 * @param[in] considerTouch treats point that on circle edge as inside if considerTouch=true  
 * @return the result whether the point is inside the circumcircle 
 */
template <typename num_type>
inline bool isInCircumCircle(const Point2D<num_type> & p1, const Point2D<num_type> & p2, const Point2D<num_type> & p3, const Point2D<num_type> & p4, bool considerTouch = true);

/**
 * @brief circumradius to shortest edge ratio of a triangle
 * @tparam point_t point type, could be Point2D
 * @param[in] p1 one of three points that form a triangle2d
 * @param[in] p2 one of three points that form a triangle2d
 * @param[in] p3 one of three points that form a triangle2d
 * @return squared ratio of circum radius to shortest edge
 * @note a triangle's circumradius-to-shortest edge ratio r/l is related to its smallest angle theta_min by the formula r/l = 1/(2*sin(theta_min))
 */
template <typename point_t, typename std::enable_if<traits::is_2d_point_t<point_t>::value, bool>::type = true>//todo point3d
inline coor_f<point_t> CircumRadius2ShortestEdgeRatioSq(const point_t & p1, const point_t & p2, const point_t & p3);

/**
 * @brief circumradius to shortest edge ratio of a triangle
 * @tparam point_t point type, could be Point2D
 * @param[in] p1 one of three points that form a triangle2d
 * @param[in] p2 one of three points that form a triangle2d
 * @param[in] p3 one of three points that form a triangle2d
 * @return ratio of circumradius to shortest edge
 * /
template <typename point_t, typename std::enable_if<traits::is_2d_point_t<point_t>::value, bool>::type = true>//todo point3d
inline coor_f<point_t> CircumRadius2ShortestEdgeRatio(const point_t & p1, const point_t & p2, const point_t & p3);

/**
 * @brief gets angle formed by two vectors
 * @tparam vector_t vector type, could be Vector2D
 * @param[in] a one of the input vectors
 * @param[in] b one of the input vectors
 * @return angle formed by vector a and b, unit: rad, range[0, 2pi)
 */
template <typename vector_t, typename std::enable_if<traits::is_2d_point_t<vector_t>::value, bool>::type = true>
inline coor_f<vector_t> Angle(const vector_t & a, const vector_t & b);

/**
 * @brief gets angle formed by three points
 * @tparam point_t point type, could be Point2D
 * @param[in] s start point 
 * @param[in] p mid point that form the angle 
 * @param[in] e end point
 * @return angle of vector(p->s) and vector(p->e), unit: rad, range[0, 2pi)
 */
template <typename point_t, typename std::enable_if<traits::is_point_t<point_t>::value, bool>::type = true>//angle of vec(p, s) and vec(p, e)
inline coor_f<point_t> Angle(const point_t & s, const point_t & p, const point_t & e) { return Angle(s - p, e - p); }

/**
 * @brief gets inner angle formed by two vectors
 * @tparam vector_t vector type, could be Vector2D or Vector3D
 * @param[in] a one of the input vectors
 * @param[in] b one of the input vectors
 * @return inner angle formed vector a and b, unit: rad, range[0, pi]
 */
template <typename vector_t, typename std::enable_if<traits::is_point_t<vector_t>::value, bool>::type = true>
inline coor_f<vector_t> InnerAngle(const vector_t & a, const vector_t & b);

/**
 * @brief gets inner angle formed by three points
 * @tparam point_t point type, could be Point2D or Point3D
 * @param[in] s start point 
 * @param[in] p mid point that form the angle 
 * @param[in] e end point
 * @return angle of vector(p->s) and vector(p->e), unit: rad, range[0, pi]
 */
template <typename point_t, typename std::enable_if<traits::is_point_t<point_t>::value, bool>::type = true>//angle of vec(p, s) and vec(p, e)
inline coor_f<point_t> InnerAngle(const point_t & s, const point_t & p, const point_t & e) { return InnerAngle(s - p, e - p); }

/**
 * @brief gets inner angle of vertex in a triangle
 * @tparam vertex vertex index of a triangle, range[0, 2]
 * @tparam triangle_t triangle type, could be Triangle2D or Triangle3D
 * @param[in] t input triangle 
 * @return inner angle of the triangle vertex, unit: rad
 */
template <size_t vertex, typename triangle_t, typename std::enable_if<traits::is_triangle_t<triangle_t>::value, bool>::type = true>
inline coor_f<triangle_t> InteriorAngle(const triangle_t & t);

/**
 * @brief gets inner angle of vertex in a triangle
 * @tparam triangle_t triangle type, could be Triangle2D or Triangle3D
 * @param t input triangle 
 * @param vertex vertex index of a triangle, range[0, 2]
 * @return inner angle of the triangle vertex, unit: rad
 */
template <typename triangle_t, typename std::enable_if<traits::is_triangle_t<triangle_t>::value, bool>::type = true>
inline coor_f<triangle_t> InteriorAngle(const triangle_t & t, const size_t vertex);

/**
 * @brief gets inner angles of three vertex in a trangle formed by three points
 * @tparam point_t point type, could be Point2D and Point3D
 * @param[in] p1 one of three points
 * @param[in] p2 one of three points
 * @param[in] p3 one of three points
 * @return three inner angles from vertex0 to vertex 2, unit: rad
 */
template <typename point_t, typename std::enable_if<traits::is_point_t<point_t>::value, bool>::type = true>
inline std::array<coor_f<point_t>, 3> InteriorAngles(const point_t & p1, const point_t & p2, const point_t & p3);

/**
 * @brief gets a point in plane that has minimum distance to the given point
 * @param[in] p the given point
 * @param[in] plane input plane
 * @return point in plane that has minimum distance to the point p
 */
template <typename num_type>
inline Point3D<float_type<num_type> > ClosestPointOnPlane(const Point3D<num_type> & p, const Plane<num_type> & plane);

/**
 * @brief gets distance of a point to the given plane
 * @param[in] p the given point 
 * @param[in] plane input plane
 * @return perpendicular distance from point p to the plane
 */
template <typename num_type>
inline float_type<num_type> PointPlaneDistance(const Point3D<num_type> & p, const Plane<num_type> & plane);

/**
 * @brief gets distance of point to a infinite line formed by two points
 * @tparam point_t point type, could be Point2D or Point3D
 * @param[in] p the given point 
 * @param[in] a one of the point that form the infinite line
 * @param[in] b one of the point that form the infinite line
 * @return squared projection distance of point to line
 */
template <typename point_t, typename std::enable_if<traits::is_point_t<point_t>::value, bool>::type = true>
inline coor_f<point_t> PointLineDistanceSq(const point_t & p, const point_t & a, const point_t & b);

/**
 * @brief gets point in a segment that has minimum distance to the given point
 * @tparam point_t point type, could be Point2D or Point3D, should have same dimension with the given segment
 * @tparam segment_t segment type, could be Segment2D or Segment3D, should have same dimension with the given point
 * @param[in] p the given point
 * @param[in] s the given segment
 * @return point in segment that has minium distance to the given point 
 */
template <typename point_t, typename segment_t,
          typename std::enable_if<traits::is_same_dim_t<point_t, segment_t>::value &&
                                  traits::is_same_coor_t<point_t, segment_t>::value, bool>::type = true>
inline point_t ClosestPointOnSegment(const point_t & p, const segment_t & s);

/**
 * @brief gets distance of given point and segment
 * @tparam point_t point type, could be Point2D or Point3D, should have same dimension with the given segment
 * @tparam segment_t segment type, could be Segment2D or Segment3D, should have same dimension with the given point
 * @param[in] p the given point
 * @param[in] s the given segment
 * @return squared distance of the given point and segment
 */
template <typename point_t, typename segment_t,
          typename std::enable_if<traits::is_same_dim_t<point_t, segment_t>::value &&
                                  traits::is_same_coor_t<point_t, segment_t>::value, bool>::type = true>
inline coor_f<point_t> PointSegmentDistanceSq(const point_t & p, const segment_t & s);

/**
 * @brief gets the point in box that has minimum distance with the given point
 * @tparam point_t point type, could be Point2D or Point3D, should have same dimension with the given box
 * @tparam box_t box type, could be Box2D or Box3D, should have same dimension with the given point
 * @param[in] p the given point 
 * @param[in] b the given box
 * @return point in box that has minimum distance with point p
 */
template <typename point_t, typename box_t,
          typename std::enable_if<traits::is_same_dim_t<point_t, box_t>::value &&
                                  traits::is_same_coor_t<point_t, box_t>::value, bool>::type = true>
inline point_t ClosestPointInBox(const point_t & p, const box_t & b);

/**
 * @brief gets distance of given point and box
 * @tparam point_t point type, could be Point2D or Point3D, should have same dimension with the given box
 * @tparam box_t box type, could be Box2D or Box3D, should have same dimension with the given point
 * @param[in] p the given point 
 * @param[in] b the given box
 * @return squared distance of given point and box
 */
template <typename point_t, typename box_t,
          typename std::enable_if<traits::is_same_dim_t<point_t, box_t>::value &&
                                  traits::is_same_coor_t<point_t, box_t>::value, bool>::type = true>
inline typename point_t::coor_t PointBoxDistanceSq(const point_t & p, const box_t & b);

/**
 * @brief checks if the point is on the left/right or on the segment
 * @param[in] p the given point
 * @param[in] seg the given segment
 * @return location result
 */
template <typename point_t, typename segment_t, 
          typename std::enable_if<traits::is_same_dim_t<point_t, segment_t>::value, bool>::type = true>
inline PointLineLocation GetPointSegmentLocation(const point_t & p, const segment_t & seg);

/**
 * @brief checks if the point is on the left/right or on the given line
 * @param[in] p the given point
 * @param[in] line the given line
 * @return location result
 */
template <typename point_t, typename line_t, 
          typename std::enable_if<traits::is_same_dim_t<point_t, line_t>::value, bool>::type = true>
inline PointLineLocation GetPointLineLocation(const point_t & p, const line_t & line);

/**
 * @brief checks if the point is on the left/right or on the given line formed by two points
 * @param[in] p the given point
 * @param[in] v1 one of the point that formed the line
 * @param[in] v2 one of the point that formed the line
 * @return location result
 */
template <typename point_t1, typename point_t2,
          typename std::enable_if<traits::is_same_dim_t<point_t1, point_t2>::value, bool>::type = true>
inline PointLineLocation GetPointLineLocation(const point_t1 & p, const point_t2 & v1, const point_t2 & v2);

/**
 * @brief checks if the point is inside/outside or on the edge of a triangle formed by three points
 * @param[in] p the given point 
 * @param[in] v1 one of three points that formed the trangle
 * @param[in] v2 one of three points that formed the trangle
 * @param[in] v3 one of three points that formed the trangle
 * @return location result 
 */
template <typename point_t1, typename point_t2,
          typename std::enable_if<traits::is_same_dim_t<point_t1, point_t2>::value, bool>::type = true>
inline PointTriangleLocation GetPointTriangleLocation(const point_t1 & p, const point_t2 & v1, const point_t2 & v2, const point_t2 & v3);

/**
 * @brief checks if the point is inside/outside or on the edge of a triangle
 * @param[in] p the given point 
 * @param[in] tri the given triangle 
 * @return location result 
 */
template <typename point_t, typename triangle_t,
          typename std::enable_if<traits::is_same_dim_t<point_t, triangle_t>::value, bool>::type = true>
inline PointTriangleLocation GetPointTriangleLocation(const point_t & p, const triangle_t & tri);

/**
 * @brief checks if two segments intersected
 * @param[in] s1 one of two segments
 * @param[in] s2 one of two segments
 * @param considerTouch treat touched segments as intersected if considerTouch=true
 * @return whether two segments intersected 
 */
template <typename segment_t1, typename segment_t2,
          typename std::enable_if<traits::is_2d_geometry_t<segment_t1>::value && 
                                  traits::is_2d_geometry_t<segment_t2>::value, bool>::type = true>
inline bool Intersects(const segment_t1 & s1, const segment_t2 & s2, bool considerTouch = true);

/**
 * @brief gets intersection points of two segments
 * @param[in] s1 one of two segments
 * @param[in] s2 one of two segments
 * @param[out] points intersection points
 * @return true if segment `s1` and `s2` intersected
 */
template <typename num_type>
inline bool Intersection(const Segment2D<num_type> & s1, const Segment2D<num_type> & s2, std::vector<Point2D<num_type> > & points);

/**
 * @brief gets intersection points of segment and infinite line
 * @param[in] s input segment 
 * @param[in] line input line 
 * @param[out] points intersect points 
 * @return true if segment intersect with line 
 */
template <typename num_type>
inline bool Intersection(const Segment2D<num_type> & s, const Line2D<num_type> & line, std::vector<Point2D<num_type> > & points);

/**
 * @brief checks if one geometry contains another one
 * @tparam geometry_t1 geometry type, could be Segment2D, Triangle2D, Box2D, Polygon2D, PolygonWithHoles2D, Box3D 
 * @tparam geometry_t2 geometry type, could be Point2D, Segment2D, Triangle2D, Box2D, Polygon2D, PolygonWithHoles2D, Point3D, Segment3D, Box3D
 * @param[in] g1 the outer geometry 
 * @param[in] g2 the inner geometry
 * @param[in] considerTouch treats `g1` contains `g2` even `g1` touched with `g2` if considerTouch=true
 * @return true whether geometry `g1` contains `g2`
 */
template <typename geometry_t1, typename geometry_t2,
          typename std::enable_if<traits::is_geometry_t<geometry_t1>::value &&
                                  traits::is_geometry_t<geometry_t2>::value, bool>::type = true>
inline bool Contains(const geometry_t1 & g1, const geometry_t2 & g2, bool considerTouch = true);

///@brief a help function to make consistance with Extent API, input could be Point2D or Point3D
template <typename point_t, typename std::enable_if<traits::is_point_t<point_t>::value, bool>::type = true>
inline point_t Extent(const point_t & p) { return p; }

///@brief gets envelope box of input segment, could be Segment2D or Segment3D
template <typename segment_t, typename std::enable_if<traits::is_segment_t<segment_t>::value, bool>::type = true>
inline box_type<segment_t> Extent(const segment_t & segment);

///@brief gets envelope box of input triangle, could be Triangle2D or Triangle3D
template <typename triangle_t, typename std::enable_if<traits::is_triangle_t<triangle_t>::value, bool>::type = true>
inline box_type<triangle_t> Extent(const triangle_t & triangle);

///@brief gets envelope box of input box, could be Box2D or Box3D
template <typename box_t, typename std::enable_if<traits::is_box_t<box_t>::value, bool>::type = true>
inline box_t Extent(const box_t & box) { return box; }

///@brief gets envelope box of input polygon, could be Polygon2D or PolygonWithHoles2D
template <typename polygon_t, typename std::enable_if<traits::is_polygon_t<polygon_t>::value ||
                                                       traits::is_polygon_with_holes_t<polygon_t>::value, bool>::type = true>
inline box_type<polygon_t> Extent(const polygon_t & polygon);

/**
 * @brief gets envelope box of a sequence with geometries
 * @tparam iterator collection iterator, the iterator value type could be one of 
 * Point2D, Segment2D, Triangle2D, Box2D, Polygon2D, PolygonWithHoles2D, Point3D, Segment3D, Triangle3D, Box3D
 * @param[in] begin  iterator to the beginning of the sequence
 * @param[in] end iterator to the ending of the sequence
 * @return Box2D with geometry2d input or Box3D with geometry3d input
 */
template <typename iterator,
          typename std::enable_if<traits::is_geometry_t<
          typename std::iterator_traits<iterator>::value_type>::value, bool>::type = true>
inline box_type<typename std::iterator_traits<iterator>::value_type> Extent(iterator begin, iterator end);

///@brief gets envelope box of input Polyline2D
template <typename num_type>
inline Box2D<num_type> Extent(const Polyline2D<num_type> & polyline);

#if GENERIC_CURRENT_BOOST_LIBRARY_VER >= 165
/**
 * @brief gets convex hull of a sequence of polygons
 * @tparam polygon_t polygon type, should be Polygon2D
 * @param[in] polygons input polygon collections
 * @return the convex hull of input polygons 
 */
template <typename polygon_t, template <typename, typename> class container, template <typename> class allocator = std::allocator>
inline polygon_t ConvexHull(const container<polygon_t, allocator<polygon_t> > & polygons);
#endif

}//namespace geometry
}//namespace generic

#include "Utility.ipp"
#endif //GENERIC_GEOMETRY_UTILITY_HPP
