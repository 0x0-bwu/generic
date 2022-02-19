/**
 * @file BoostGeometryRegister.hpp
 * @author bwu
 * @brief Adaption of the generic geometry models to boost geometry concept
 * @version 0.1
 * @date 2022-02-22
 */
#ifndef GENERIC_GEOMETRY_BOOSTGEOMETRYREGISTER_HPP
#define GENERIC_GEOMETRY_BOOSTGEOMETRYREGISTER_HPP
#include <boost/geometry/core/coordinate_dimension.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <boost/geometry/core/coordinate_type.hpp>
#include <boost/geometry/core/mutable_range.hpp>
#include <boost/geometry/core/point_order.hpp>
#include <boost/range/mutable_iterator.hpp>
#include <boost/geometry/core/closure.hpp>
#include <boost/geometry/core/access.hpp>
#include <boost/range/const_iterator.hpp>
#include <boost/geometry/core/tags.hpp>
#include <boost/geometry/geometry.hpp>
#include <boost/geometry/core/cs.hpp>
#include "Geometries.hpp"

namespace generic { 
namespace geometry{
    template <typename num_type>
    class PolygonWithHoles2DRange;
}
}

namespace boost {
namespace geometry {
namespace traits {

using namespace generic::geometry;
//point
template <typename num_type>
struct tag<Point2D<num_type> >
{ using type = point_tag; };

template <typename num_type>
struct dimension<Point2D<num_type> > : boost::mpl::int_<2> {};

template <typename num_type>
struct coordinate_type<Point2D<num_type> >
{ using type = num_type; };

template <typename num_type>
struct coordinate_system<Point2D<num_type> >
{ using type = cs::cartesian; };

template <typename num_type>
struct access<Point2D<num_type>, 0>
{
    static inline num_type get(const Point2D<num_type> & p) { return p[0]; }
    static inline void set(Point2D<num_type> & p, const num_type & n) { p[0] = n; }
};

template <typename num_type>
struct access<Point2D<num_type>, 1>
{
    static inline num_type get(const Point2D<num_type> & p) { return p[1]; }
    static inline void set(Point2D<num_type> & p, const num_type & n) { p[1] = n; }
};

//segment
template <typename num_type>
struct tag<Segment2D<num_type> >
{ using type = segment_tag; };

template <typename num_type>
struct point_type<Segment2D<num_type> >
{ using type = Point2D<num_type>; };

template <typename num_type, size_t dimension>
struct indexed_access<Segment2D<num_type>, 0, dimension>
{
    static inline num_type get(const Segment2D<num_type> & s) { return geometry::get<dimension>(s[0]); }
    static inline void set(Segment2D<num_type> & s, const num_type & n) { geometry::set<dimension>(s[0], n); }
};

template <typename num_type, size_t dimension>
struct indexed_access<Segment2D<num_type>, 1, dimension>
{
    static inline num_type get(const Segment2D<num_type> & s) { return geometry::get<dimension>(s[1]); }
    static inline void set(Segment2D<num_type> & s, const num_type & n) { geometry::set<dimension>(s[1], n); }
};

//box
template <typename num_type>
struct tag<Box2D<num_type> >
{ using type = box_tag; };

template <typename num_type>
struct point_type<Box2D<num_type> >
{ using type = Point2D<num_type>; };

template <typename num_type, size_t dimension>
struct indexed_access<Box2D<num_type>, 0, dimension>
{
    static inline num_type get(const Box2D<num_type> & b) { return geometry::get<dimension>(b[0]); }
    static inline void set(Box2D<num_type> & b, const num_type & n) { geometry::set<dimension>(b[0], n); }
};

template <typename num_type, size_t dimension>
struct indexed_access<Box2D<num_type>, 1, dimension>
{
    static inline num_type get(const Box2D<num_type> & b) { return geometry::get<dimension>(b[1]); }
    static inline void set(Box2D<num_type> & b, const num_type & n) { geometry::set<dimension>(b[1], n); }
};

//triangle
template <typename num_type>
struct tag<Triangle2D<num_type> >
{ using type = ring_tag; };

template <typename num_type>
struct closure<Triangle2D<num_type> >
{ static const closure_selector value = open; };

//linestring
template <typename num_type>
struct tag<Polyline2D<num_type> >
{ using type = linestring_tag; };

//polygon
template <typename num_type>
struct tag<Polygon2D<num_type> >
{ using type = ring_tag; };

template <typename num_type>
struct closure<Polygon2D<num_type> >
{ static const closure_selector value = open; };

template <typename num_type>
struct clear<Polygon2D<num_type> >
{
    static inline void apply(Polygon2D<num_type> & polygon)
    {
        polygon.Clear();
    }
};

template <typename num_type>
struct push_back<Polygon2D<num_type> >
{
    static inline void apply(Polygon2D<num_type> & polygon, const Point2D<num_type> & point)
    {
        polygon.Insert(polygon.End(), point);
    }
};

template <typename num_type>
struct resize<Polygon2D<num_type> >
{
    using point_type = typename Polygon2D<num_type>::point_t;
    static inline void apply(Polygon2D<num_type> & polygon, size_t size)
    {
        polygon.Resize(size);
    }
};

//polygon with holes
template <typename num_type>
struct tag<PolygonWithHoles2D<num_type> >
{ using type = polygon_tag; };

template <typename num_type>
struct ring_const_type<PolygonWithHoles2D<num_type> >
{ using type = const typename PolygonWithHoles2D<num_type>::outline_type; };

template <typename num_type>
struct ring_mutable_type<PolygonWithHoles2D<num_type> >
{ using type = typename PolygonWithHoles2D<num_type>::outline_type; };

template <typename num_type>
struct interior_const_type<PolygonWithHoles2D<num_type> >
{ using type = const PolygonWithHoles2DRange<num_type>; };

template <typename num_type>
struct interior_mutable_type<PolygonWithHoles2D<num_type> >
{ using type = PolygonWithHoles2DRange<num_type>; };

template <typename num_type>
struct exterior_ring<PolygonWithHoles2D<num_type> >
{
    using outline_type = typename PolygonWithHoles2D<num_type>::outline_type;
    using const_outline_type = typename PolygonWithHoles2D<num_type>::outline_type const;
    static inline outline_type & get(PolygonWithHoles2D<num_type> & pwh)
    {
        return pwh.outline;
    }

    static inline const_outline_type & get(const PolygonWithHoles2D<num_type> & pwh)
    {
        return pwh.outline;
    }
};

template <typename num_type>
struct interior_rings<PolygonWithHoles2D<num_type> >
{
    static inline PolygonWithHoles2DRange<num_type> get(PolygonWithHoles2D<num_type> & pwh)
    {
        return PolygonWithHoles2DRange<num_type>(pwh.holes.begin(), pwh.holes.end());
    }

    static inline const PolygonWithHoles2DRange<num_type> get(const PolygonWithHoles2D<num_type> & pwh)
    {
        return PolygonWithHoles2DRange<num_type>(pwh.holes.begin(), pwh.holes.end());
    }
};

}//namespace traits
}//namespace geometry
}//namespace boost

namespace boost {
    
using namespace generic::geometry;
//triangle
template <typename num_type>
struct range_mutable_iterator<Triangle2D<num_type> >
{
    using type = typename Triangle2D<num_type>::iterator;
};

template <typename num_type>
struct range_const_iterator<Triangle2D<num_type> >
{
    using type = typename Triangle2D<num_type>::const_iterator;
};

template<typename num_type>
struct range_size<Triangle2D<num_type> >
{
    using type = std::size_t;
};

//polygon
template <typename num_type>
struct range_mutable_iterator<Polygon2D<num_type> >
{
    using type = typename Polygon2D<num_type>::point_iterator;
};

template <typename num_type>
struct range_const_iterator<Polygon2D<num_type> >
{
    using type = typename Polygon2D<num_type>::const_point_iterator;
};

template<typename num_type>
struct range_size<Polygon2D<num_type> >
{
    using type = std::size_t;
};

//polygon with hole
template <typename num_type>
struct range_mutable_iterator<PolygonWithHoles2DRange<num_type> >
{
    using type = typename PolygonWithHoles2DRange<num_type>::hole_iterator;
};

template <typename num_type>
struct range_const_iterator<PolygonWithHoles2DRange<num_type> >
{
    using type = typename PolygonWithHoles2DRange<num_type>::const_hole_iterator;
};

template<typename num_type>
struct range_size<PolygonWithHoles2DRange<num_type> >
{
    using type = std::size_t;
};

}//namesapce boost

namespace generic {
namespace geometry{
//triangle
template<typename num_type>
inline typename Triangle2D<num_type>::iterator
range_begin(Triangle2D<num_type> & triangle)
{
    return triangle.Begin();
}

template<typename num_type>
inline typename Triangle2D<num_type>::const_iterator
range_begin(const Triangle2D<num_type> & triangle)
{
    return triangle.ConstBegin();
}

template<typename num_type>
inline typename Triangle2D<num_type>::iterator
range_end(Triangle2D<num_type> & triangle)
{
    return triangle.End();
}

template<typename num_type>
inline typename Triangle2D<num_type>::const_iterator
range_end(const Triangle2D<num_type> & triangle)
{
    return triangle.ConstEnd();
}

template<typename num_type>
inline std::size_t range_calculate_size(const Triangle2D<num_type> & triangle)
{
    return 3;
}

//polygon
template<typename num_type>
inline typename Polygon2D<num_type>::point_iterator
range_begin(Polygon2D<num_type> & polygon)
{
    return polygon.Begin();
}

template<typename num_type>
inline typename Polygon2D<num_type>::const_point_iterator
range_begin(const Polygon2D<num_type> & polygon)
{
    return polygon.ConstBegin();
}

template<typename num_type>
inline typename Polygon2D<num_type>::point_iterator
range_end(Polygon2D<num_type> & polygon)
{
    return polygon.End();
}

template<typename num_type>
inline typename Polygon2D<num_type>::const_point_iterator
range_end(const Polygon2D<num_type> & polygon)
{
    return polygon.ConstEnd();
}

template<typename num_type>
inline std::size_t range_calculate_size(const Polygon2D<num_type> & polygon)
{
    return polygon.Size();
}

//polygon with holes
template <typename num_type>
class PolygonWithHoles2DRange
{
public:
    using hole_iterator = typename PolygonWithHoles2D<num_type>::hole_iterator;
    using const_hole_iterator = typename PolygonWithHoles2D<num_type>::const_hole_iterator;
    PolygonWithHoles2DRange() = delete;
    PolygonWithHoles2DRange(hole_iterator begin, hole_iterator end)
        : m_begin(begin), m_end(end), m_bConst(false) {}
    
    PolygonWithHoles2DRange(const_hole_iterator begin, const_hole_iterator end)
        : m_constBegin(begin), m_constEnd(end), m_bConst(true) {}

    hole_iterator Begin() { assert(!m_bConst); return m_begin; }
    hole_iterator End() { assert(!m_bConst); return m_end; }
    const_hole_iterator ConstBegin() const { assert(m_bConst); return m_constBegin; }
    const_hole_iterator ConstEnd() const { assert(m_bConst); return m_constEnd; }

    std::size_t Size() const;

private:
    hole_iterator m_begin;
    hole_iterator m_end;
    const_hole_iterator m_constBegin;
    const_hole_iterator m_constEnd;
    bool m_bConst = false;
};

template <typename num_type>
inline std::size_t PolygonWithHoles2DRange<num_type>::Size() const
{
    if(m_bConst) return std::distance(m_constBegin, m_constEnd);
    else return std::distance(m_begin, m_end);
}

template <typename num_type>
inline typename PolygonWithHoles2DRange<num_type>::hole_iterator
range_begin(PolygonWithHoles2DRange<num_type> & range) { return range.Begin(); }

template <typename num_type>
inline typename PolygonWithHoles2DRange<num_type>::const_hole_iterator
range_begin(const PolygonWithHoles2DRange<num_type> & range) { return range.ConstBegin(); }

template <typename num_type>
inline typename PolygonWithHoles2DRange<num_type>::hole_iterator
range_end(PolygonWithHoles2DRange<num_type> & range) { return range.End(); }

template <typename num_type>
inline typename PolygonWithHoles2DRange<num_type>::const_hole_iterator
range_end(const PolygonWithHoles2DRange<num_type> & range) { return range.ConstEnd(); }

template<typename num_type>
inline std::size_t range_calculate_size(const PolygonWithHoles2DRange<num_type> & range)
{ return range.Size(); }

}//namespace geometry
}//namespace generic

//adaptor
namespace generic  {
namespace geometry {

namespace bg = boost::geometry;

template <typename num_type>
using boost_point2d_t = bg::model::point<num_type, 2, bg::cs::cartesian>;

template <typename num_type, bool clockwise = false, bool closed = false,
          template <typename, typename> class container = std::vector,
          template <typename> class allocator = std::allocator>
using boost_ring_t = bg::model::ring<boost_point2d_t<num_type>, clockwise, closed, container, allocator>;

template <typename num_type, bool clockwise = false, bool closed = false,
          template <typename, typename> class point_list = std::vector,
          template <typename, typename> class ring_list = std::vector,
          template <typename> class point_alloc = std::allocator,
          template <typename> class ring_alloc = std::allocator>
using boost_polygon_t = bg::model::polygon<boost_point2d_t<num_type>, clockwise, closed, point_list, ring_list, point_alloc, ring_alloc>;

template <typename num_type>
inline Point2D<num_type> toPoint2D(const boost_point2d_t<num_type> & point)
{
    return Point2D<num_type>(bg::get<0>(point), bg::get<1>(point));
}

template <typename num_type, bool clockwise = false, bool closed = false,
          template <typename, typename> class container = std::vector,
          template <typename> class allocator = std::allocator>
inline Polygon2D<num_type> toPolygon2D(const boost_ring_t<num_type, clockwise, closed, container, allocator> & ring)
{
    Polygon2D<num_type> polygon;
    for(const auto & point : ring)
        polygon << toPoint2D(point);
    if(polygon.Size() && (polygon.Back() == polygon.Front()))
        polygon.PopBack();
    return polygon;
}

template <typename num_type, bool clockwise = false, bool closed = false,
          template <typename, typename> class point_list = std::vector,
          template <typename, typename> class ring_list = std::vector,
          template <typename> class point_alloc = std::allocator,
          template <typename> class ring_alloc = std::allocator>
inline PolygonWithHoles2D<num_type> toPolygonWithHoles2D(const boost_polygon_t<num_type, clockwise, closed, point_list, ring_list, point_alloc, ring_alloc> & polygon)
{
    PolygonWithHoles2D<num_type> pwh;
    pwh.outline = toPolygon2D(polygon.outer());
    for(const auto & hole : polygon.inners())
        pwh.holes.emplace_back(toPolygon2D(hole));
    return pwh;
}
}//namespace geometry
}//namespace generic

#endif//GENERIC_GEOMETRY_BOOSTGEOMETRYREGISTER_HPP