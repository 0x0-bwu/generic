#ifndef GENERIC_GEOMETRY_GEOMETRYTRAITS_HPP
#define GENERIC_GEOMETRY_GEOMETRYTRAITS_HPP
#include "generic/common/Traits.hpp"
#include "Geometries.hpp"
#include "Vector.hpp"
namespace generic  {
namespace geometry {

namespace traits{

template <typename geometry_t>
struct is_point_t : public std::false_type {};

template <typename num_type>
struct is_point_t<Point2D<num_type> > : public std::true_type {};

template <typename num_type>
struct is_point_t<Point3D<num_type> > : public std::true_type {};

template <typename geometry_t>
struct is_segment_t : public std::false_type {};

template <typename num_type>
struct is_segment_t<Segment2D<num_type> > : public std::true_type {};

template <typename num_type>
struct is_segment_t<Segment3D<num_type> > : public std::true_type {};

template <typename geometry_t>
struct is_triangle_t : public std::false_type {};

template <typename num_type>
struct is_triangle_t<Triangle2D<num_type> > : public std::true_type {};

template <typename num_type>
struct is_triangle_t<Triangle3D<num_type> > : public std::true_type {};

template <typename geometry_t>
struct is_box_t : public std::false_type {};

template <typename num_type>
struct is_box_t<Box2D<num_type> > : public std::true_type {};

template <typename num_type>
struct is_box_t<Box3D<num_type> > : public std::true_type {};

template <typename geometry_t>
struct is_polyline_t : public std::false_type {};

template <typename num_type>
struct is_polyline_t<Polyline2D<num_type> > : public std::true_type {};

template <typename geometry_t>
struct is_polygon_t : public std::false_type {};

template <typename num_type>
struct is_polygon_t<Polygon2D<num_type> > : public std::true_type {};

template <typename geometry_t>
struct is_polygon_with_holes_t : public std::false_type {};

template <typename num_type>
struct is_polygon_with_holes_t<PolygonWithHoles2D<num_type> > : public std::true_type {};

template <typename geometry_t>
struct is_vector_t : public std::false_type {};

template <typename num_type>
struct is_vector_t<Vector2D<num_type> > : public std::true_type {};

template <typename num_type>
struct is_vector_t<Vector3D<num_type> > : public std::true_type {};

template <typename num_type, size_t N>
struct is_vector_t<VectorN<num_type, N> > : public std::true_type {};

template <typename point_t>
using is_2d_point_t = typename std::conditional<is_point_t<point_t>::value && point_t::dim == 2, std::true_type, std::false_type>::type;

template <typename point_t>
using is_3d_point_t = typename std::conditional<is_point_t<point_t>::value && point_t::dim == 3, std::true_type, std::false_type>::type;

template <typename vector_t>
using is_2d_vector_t = is_2d_point_t<vector_t>;

template <typename vector_t>
using is_3d_vector_t = is_3d_point_t<vector_t>;

template <typename geometry_t>
using is_geometry_t = typename std::conditional<is_point_t<geometry_t>::value ||
                                                is_segment_t<geometry_t>::value ||
                                                is_triangle_t<geometry_t>::value ||
                                                is_box_t<geometry_t>::value ||
                                                is_polygon_t<geometry_t>::value ||
                                                is_polygon_with_holes_t<geometry_t>::value, std::true_type, std::false_type>::type;
                                                
template <typename geometry_t>
using is_2d_geometry_t = typename std::conditional<is_geometry_t<geometry_t>::value && geometry_t::dim == 2, std::true_type, std::false_type>::type;

template <typename geometry_t>
using is_3d_geometry_t = typename std::conditional<is_geometry_t<geometry_t>::value && geometry_t::dim == 3, std::true_type, std::false_type>::type;

template <typename geometry_t>
using is_2d_surf_geom_t = typename std::conditional<(is_triangle_t<geometry_t>::value ||
                                                     is_box_t<geometry_t>::value ||
                                                     is_polygon_t<geometry_t>::value ||
                                                     is_polygon_with_holes_t<geometry_t>::value) && is_2d_geometry_t<geometry_t>::value,
                                                     std::true_type, std::false_type>::type;

template <typename geometry_t1, typename geometry_t2>
using is_same_dim_t = typename std::conditional<geometry_t1::dim == geometry_t2::dim, std::true_type, std::false_type>::type;

using same_coor_t = std::true_type;
using different_coor_t = std::false_type;

template <typename geometry_t1, typename geometry_t2>
using is_same_coor_t = typename std::conditional<std::is_same<typename geometry_t1::coor_t, typename geometry_t2::coor_t>::value, same_coor_t, different_coor_t>::type;

template <typename geometry_t>
using is_floating_t = typename std::conditional<std::is_floating_point<typename geometry_t::coor_t>::value, std::true_type, std::false_type>::type;

template <typename coor_t1, typename coor_t2>
using narrow_coor_t = typename std::conditional<sizeof(coor_t1) < sizeof(coor_t2), coor_t1, coor_t2>::type;

template <typename coor_t1, typename coor_t2>
using wide_coor_t = typename std::conditional<sizeof(coor_t1) < sizeof(coor_t2), coor_t2, coor_t1>::type;

//if all integrals, return integral type with wider size
//if all floatings, return floating type with narrower size
//if integrals and floatings, return floating
template <typename... args>
struct common_coor_t;

template <typename coor_t>
struct common_coor_t<coor_t>
{
    using type = coor_t;
};

template <typename coor_t1, typename coor_t2>
struct common_coor_t<coor_t1, coor_t2>
{
    using type = typename std::conditional<
                 common::integral_type_check<coor_t1, coor_t2>::value, wide_coor_t<coor_t1, coor_t2>,
                 typename std::conditional<
                 common::floating_type_check<coor_t1, coor_t2>::value, narrow_coor_t<coor_t1, coor_t2>,
                 typename std::conditional<
                 std::is_floating_point<coor_t1>::value, coor_t1, coor_t2>::type>::type>::type;
};

template <typename coor_t, typename... args>
struct common_coor_t<coor_t, args...>
{
    using type = typename common_coor_t<coor_t, typename common_coor_t<args...>::type>::type;
};

}//namespace traits

template <typename geometry_t, typename std::enable_if<traits::is_geometry_t<geometry_t>::value || 
                                                       traits::is_vector_t<geometry_t>::value, bool>::type = true>
using coor_f = common::float_type<typename geometry_t::coor_t>;

template <typename geometry_t, typename std::enable_if<traits::is_geometry_t<geometry_t>::value, bool>::type = true>
using point_f = typename std::conditional<traits::is_2d_geometry_t<geometry_t>::value, 
                                            Point2D<coor_f<geometry_t> >,
                                            Point3D<coor_f<geometry_t> > >::type;

template <typename vector_t, typename std::enable_if<traits::is_vector_t<vector_t>::value, bool>::type = true>
using vector_f = typename std::conditional<traits::is_2d_point_t<vector_t>::value, Vector2D<coor_f<vector_t> >,
                 typename std::conditional<traits::is_3d_point_t<vector_t>::value, Vector3D<coor_f<vector_t> >,
                                           VectorN<common::float_type<typename vector_t::coor_t>, vector_t::dim> >::type>::type;

template <typename geometry_t, typename std::enable_if<traits::is_geometry_t<geometry_t>::value, bool>::type = true>
using segment_type = typename std::conditional<traits::is_2d_geometry_t<geometry_t>::value, 
                                            Segment2D<typename geometry_t::coor_t>,
                                            Segment3D<typename geometry_t::coor_t> >::type;

template <typename geometry_t, typename std::enable_if<traits::is_geometry_t<geometry_t>::value, bool>::type = true>
using triangle_type = typename std::conditional<traits::is_2d_geometry_t<geometry_t>::value, 
                                            Triangle2D<typename geometry_t::coor_t>,
                                            Triangle3D<typename geometry_t::coor_t> >::type;

template <typename geometry_t, typename std::enable_if<traits::is_geometry_t<geometry_t>::value, bool>::type = true>
using box_type = typename std::conditional<traits::is_2d_geometry_t<geometry_t>::value, 
                                            Box2D<typename geometry_t::coor_t>,
                                            Box3D<typename geometry_t::coor_t> >::type;

template <typename geometry_t, typename std::enable_if<traits::is_geometry_t<geometry_t>::value, bool>::type = true>
using polygon_type = Polygon2D<typename geometry_t::coor_t>;

template <typename geometry_t, typename std::enable_if<traits::is_geometry_t<geometry_t>::value, bool>::type = true>
using polygon_with_holes_type = PolygonWithHoles2D<typename geometry_t::coor_t>;

struct geometry_2d_point_tag {};
struct geometry_2d_segment_tag {};
struct geometry_2d_triangle_tag {};
struct geometry_2d_box_tag {};
struct geometry_2d_polygon_tag {};
struct geometry_2d_polygon_with_holes_tag{};
struct geometry_3d_point_tag {};
struct geometry_3d_segment_tag {};
struct geometry_3d_triangle_tag {};
struct geometry_3d_box_tag {};
struct geometry_unkonwn_tag {};

template <typename geometry_t>
struct geometry_tag {};

template <typename num_type>
struct geometry_tag<Point2D<num_type> > { using tag = geometry_2d_point_tag; };

template <typename num_type>
struct geometry_tag<Segment2D<num_type> > { using tag = geometry_2d_segment_tag; };

template <typename num_type>
struct geometry_tag<Triangle2D<num_type> > { using tag = geometry_2d_triangle_tag; };

template <typename num_type>
struct geometry_tag<Box2D<num_type> > { using tag = geometry_2d_box_tag; };

template <typename num_type>
struct geometry_tag<Polygon2D<num_type> > { using tag = geometry_2d_polygon_tag; };

template <typename num_type>
struct geometry_tag<PolygonWithHoles2D<num_type> > { using tag = geometry_2d_polygon_with_holes_tag; };

template <typename num_type>
struct geometry_tag<Point3D<num_type> > { using tag = geometry_3d_point_tag; };

template <typename num_type>
struct geometry_tag<Segment3D<num_type> > { using tag = geometry_3d_segment_tag; };

template <typename num_type>
struct geometry_tag<Triangle3D<num_type> > { using tag = geometry_3d_triangle_tag; };

template <typename num_type>
struct geometry_tag<Box3D<num_type> > { using tag = geometry_3d_box_tag; };

}//namespace geometry
}//namespace generic
#endif//GENERIC_GEOMETRY_GEOMETRYTRAITS_HPP