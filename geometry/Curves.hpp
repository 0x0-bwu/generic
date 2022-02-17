#ifndef GENERIC_GEOMETRY_CURVES_HPP
#define GENERIC_GEOMETRY_CURVES_HPP
#include "generic/math/MathUtility.hpp"
#include "Point.hpp"
#include <complex>
namespace generic {
namespace geometry{

template <typename num_type>
class Arc
{
    using float_t = float_type<num_type>;
public:
    static const size_t dim = 2;
    float_t radian;//+ccw, -cw
    Point2D<num_type> origin, start;
    Arc(const Point2D<num_type> & origin, const Point2D<num_type> & start, float_t radian)
     : radian(radian), origin(origin), start(start) {}

    void GetStartAlphaMag(float_t & alpha, float_t & mag) const;
    Point2D<num_type> EndPoint() const;
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
};

template <typename num_type>
void Arc<num_type>::GetStartAlphaMag(float_t & alpha, float_t & mag) const
{
    auto v = start - origin;
    auto c = std::complex<float_t>(v[0], v[1]);
    mag = std::abs(c); alpha = std::arg(c);
    if(alpha < 0) alpha += math::pi_2;
}

template <typename num_type>
Point2D<num_type> Arc<num_type>::EndPoint() const
{
    float_t alpha, mag;
    GetStartAlphaMag(alpha, mag);
    alpha += radian;
    return Point2D<num_type>(mag * std::cos(alpha), mag * std::sin(alpha)) + origin;
}

}//namespace geometry
}//namespace generic
#endif//GENERIC_GEOMETRY_CURVES_HPP
