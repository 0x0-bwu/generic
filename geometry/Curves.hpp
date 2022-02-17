#ifndef GENERIC_GEOMETRY_CURVES_HPP
#define GENERIC_GEOMETRY_CURVES_HPP
#include "generic/math/MathUtility.hpp"
#include "Predicates.hpp"
#include "Sphere.hpp"
namespace generic {
namespace geometry{

template <typename num_type>
class Arc
{
    using float_t = float_type<num_type>;
public:
    static const size_t dim = 2;
    Point2D<num_type> start, origin, end;
    bool isCCW;
    Arc(const Point2D<num_type> & start, const Point2D<num_type> & origin, const Point2D<num_type> & end, bool ccw = true)
     : start(start), origin(origin), end(end), isCCW(ccw) {}
};

template <typename num_type>
class Arc3
{
    using float_t = float_type<num_type>;
public:
    static const size_t dim = 2;
    Point2D<num_type> start, mid, end;
    Arc3(const Point2D<num_type> & start, const Point2D<num_type> & mid, const Point2D<num_type> & end)
     : start(start), mid(mid), end(end) {}

    Arc<num_type> toArc() const
    {   
        float_t r2;
        auto origin = Circle<num_type>::CircumCircle(start, mid, end, r2).template Cast<num_type>();
        auto isCCW = Point2D<num_type>::isCCW(start, mid, end);
        return Arc<num_type>(start, origin, end, isCCW);
    }
};

}//namespace geometry
}//namespace generic
#endif//GENERIC_GEOMETRY_CURVES_HPP
