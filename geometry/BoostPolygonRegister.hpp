#ifndef GENERIC_GEOMETRY_BOOSTPOLYGONREGISTER_HPP
#define GENERIC_GEOMETRY_BOOSTPOLYGONREGISTER_HPP
#include "Geometries.hpp"
#include <boost/polygon/polygon.hpp>
#include <type_traits>
namespace boost {
namespace polygon {

using namespace generic::geometry;
//point
template <typename num_type>
struct geometry_concept<Point2D<num_type> >
{
    using type = point_concept;
};

template <typename num_type>
struct point_traits<Point2D<num_type> >
{
    using coordinate_type = num_type ;
    static inline coordinate_type get(const Point2D<num_type> & point, orientation_2d orient)
    {
        if(orient == HORIZONTAL) return point[0];
        else return point[1];
    }
};

template <typename num_type>
struct point_mutable_traits<Point2D<num_type> >
{
    using coordinate_type = num_type ;
    static inline void set(Point2D<num_type> & point, orientation_2d orient, num_type value)
    {
        if(orient == HORIZONTAL) point[0] = value;
        else point[1] = value;
    }

    static inline Point2D<num_type> construct(num_type x, num_type y)
    {
        return Point2D<num_type>(x, y);
    }

    static inline Point2D<num_type> construct(const point_data<num_type> & p)
    {
        return Point2D<num_type>(p.get(HORIZONTAL), p.get(VERTICAL));
    }
};

//segment
template <typename num_type>
struct geometry_concept<Segment2D<num_type> >
{
    using type = segment_concept;
};

template <typename num_type>
struct segment_traits<Segment2D<num_type> >
{
    using coordinate_type = num_type ;
    using point_type = Point2D<num_type> ;

    static inline point_type get(const Segment2D<num_type> & segment, direction_1d dir)
    {
        if(dir == LOW) return segment[0];
        return segment[1];
    }
};

template <typename num_type>
struct segment_mutable_traits<Segment2D<num_type> > {
    using coordinate_type = num_type;
    using point_type = Point2D<num_type>;

    static inline void set(Segment2D<num_type> & segment, direction_1d dir, const point_type & point)
    {
        if(dir == LOW) segment[0] = point;
        else segment[1] = point;
    }

    static inline Segment2D<num_type> construct(const point_type & low, const point_type & high)
    {
        return Segment2D<num_type>(low, high);
    }

    static inline Segment2D<num_type> construct(const point_data<num_type> & low, const point_data<num_type> & high)
    {
        Point2D<num_type> p1(low.get(HORIZONTAL), low.get(VERTICAL));
        Point2D<num_type> p2(high.get(HORIZONTAL), high.get(VERTICAL));
        return Segment2D<num_type>(p1, p2);
    }
};

//triangle
template <typename num_type>
struct geometry_concept<Triangle2D<num_type> >
{
    using type = polygon_concept;
};

template <typename num_type>
struct polygon_traits_general<Triangle2D<num_type> >
{
    using coordinate_type = num_type;
    using point_type = Point2D<num_type> ;
    using iterator_type = typename Triangle2D<num_type>::const_iterator;

    static inline iterator_type begin_points(const Triangle2D<num_type> & triangle)
    {
        return triangle.ConstBegin();
    }

    static inline iterator_type end_points(const Triangle2D<num_type> & triangle)
    {
        return triangle.ConstEnd();
    }

    static inline std::size_t size(const Triangle2D<num_type> & triangle)
    {
        return 3;
    }

    static inline winding_direction winding(const Triangle2D<num_type> & triangle)
    {
        if(triangle.isCCW()) return counterclockwise_winding;
        else return clockwise_winding;
    }
};

//box
template <typename num_type>
struct geometry_concept<Box2D<num_type> >
{
    using type = rectangle_concept;
};

template <typename num_type>
struct rectangle_traits<Box2D<num_type> >
{
    using coordinate_type = num_type;
    using interval_type = interval_data<num_type>;
    static inline interval_type get(const Box2D<num_type> & box, orientation_2d orient)
    {
        if(orient == HORIZONTAL) return interval_type(box[0][0], box[1][0]);
        else return interval_type(box[0][1], box[1][1]);
    }
};

template <typename num_type>
struct rectangle_mutable_traits<Box2D<num_type> >
{
    using interval_type = interval_data<num_type>;
    static inline void set(Box2D<num_type> & box, orientation_2d orient, const interval_type & interval)
    {
        if (orient == HORIZONTAL) {
            box[0][0] = interval.low();
            box[1][0] = interval.high();
        }
        else {
            box[0][1] = interval.low();
            box[1][1] = interval.high();
        }
    }

    static inline Box2D<num_type> construct(const interval_type & interval_horizontal, const interval_type & interval_vertical)
    {
        return Box2D<num_type>(interval_horizontal.low(), interval_vertical.low(), interval_horizontal.high(), interval_vertical.high());
    }
};

//polygon
template <typename num_type>
struct geometry_concept<Polygon2D<num_type> >
{
    using type = polygon_concept ;
};

template <typename num_type>
struct polygon_traits_general<Polygon2D<num_type> >
{
    using coordinate_type = typename Polygon2D<num_type>::coor_t;
    using point_type = typename Polygon2D<num_type>::point_t;
    using iterator_type = typename Polygon2D<num_type>::const_point_iterator;
    static inline iterator_type begin_points(const Polygon2D<num_type> & polygon)
    {
        return polygon.ConstBegin();
    }

    static inline iterator_type end_points(const Polygon2D<num_type> & polygon)
    {
        return polygon.ConstEnd();
    }
    
    static inline std::size_t size(const Polygon2D<num_type> & polygon)
    {
        return polygon.Size();
    }

    static inline winding_direction winding(const Polygon2D<num_type> &)
    {
        return unknown_winding;
    }
};

template <typename num_type>
struct polygon_mutable_traits<Polygon2D<num_type> >
{
    template <typename iterator, typename std::enable_if<std::is_same<
              typename std::iterator_traits<iterator>::value_type, point_data<num_type> >::value, bool>::type = true>
    static inline Polygon2D<num_type> & set_points(Polygon2D<num_type> & polygon,
                                                    iterator input_begin, iterator input_end)
    {
        polygon.Clear();
        polygon.GetPoints().reserve(std::distance(input_begin, input_end));
        for(auto iter = input_begin; iter != input_end; ++iter)
            polygon << point_mutable_traits<Point2D<num_type> >::construct(*iter);
        return polygon;
    }

    template <typename iterator, typename std::enable_if<std::is_same<
              typename std::iterator_traits<iterator>::value_type, Point2D<num_type> >::value, bool>::type = true>
    static inline Polygon2D<num_type> & set_points(Polygon2D<num_type> & polygon,
                                                    iterator input_begin, iterator input_end)
    {
        polygon.Clear();
        polygon.Insert(polygon.End(), input_begin, input_end);
        return polygon;
    }
};

//polygon with holes
template <typename num_type>
struct geometry_concept<PolygonWithHoles2D<num_type> >
{
    using type = polygon_with_holes_concept ;
};

template <typename num_type>
struct polygon_traits_general<PolygonWithHoles2D<num_type> >
{
    using coordinate_type = typename PolygonWithHoles2D<num_type>::coor_t;
    using point_type = typename PolygonWithHoles2D<num_type>::point_t;
    using iterator_type = typename PolygonWithHoles2D<num_type>::const_point_iterator;
    static inline iterator_type begin_points(const PolygonWithHoles2D<num_type> & pwh)
    {
        return pwh.outline.ConstBegin();
    }

    static inline iterator_type end_points(const PolygonWithHoles2D<num_type> & pwh)
    {
        return pwh.outline.ConstEnd();
    }
    
    static inline std::size_t size(const PolygonWithHoles2D<num_type> & pwh)
    {
        return pwh.outline.Size();
    }

    static inline winding_direction winding(const PolygonWithHoles2D<num_type> & pwh)
    {
        return unknown_winding;
    }
};

template <typename num_type>
struct polygon_mutable_traits<PolygonWithHoles2D<num_type> >
{
    template <typename iterator>
    static inline PolygonWithHoles2D<num_type> & set_points(PolygonWithHoles2D<num_type> & pwh,
                                                    iterator input_begin, iterator input_end)
    {
        pwh.outline.Clear();
        pwh.outline.Insert(pwh.outline.End(), input_begin, input_end);
        return pwh;
    }
};

template <typename num_type>
struct polygon_with_holes_traits<PolygonWithHoles2D<num_type> >
{
    using hole_type = typename PolygonWithHoles2D<num_type>::hole_type;
    using iterator_holes_type = typename PolygonWithHoles2D<num_type>::const_hole_iterator;

    static inline iterator_holes_type begin_holes(const PolygonWithHoles2D<num_type> & pwh)
    {
        return pwh.ConstBeginHoles();
    }

    static inline iterator_holes_type end_holes(const PolygonWithHoles2D<num_type> & pwh)
    {
        return pwh.ConstEndHoles();
    }

    static inline std::size_t size_holes(const PolygonWithHoles2D<num_type> & pwh)
    {
        return pwh.HolesSize();
    }
};

template <typename num_type>
struct polygon_with_holes_mutable_traits<PolygonWithHoles2D<num_type> >
{
    template <typename iterator>
    static inline PolygonWithHoles2D<num_type> & set_holes(PolygonWithHoles2D<num_type> & pwh,
                                                                iterator input_begin, iterator input_end)
    {
        pwh.holes.clear();
        pwh.holes.insert(pwh.holes.end(), input_begin, input_end);
        return pwh;
    }
};

}//namespace polygon
}//namespace boost

#endif//GENERIC_GEOMETRY_BOOSTPOLYGONREGISTER_HPP
