/**
 * @file Curves.hpp
 * @author bwu
 * @brief Model of some curve concepts
 * @version 0.1
 * @date 2022-02-14 
 */
#ifndef GENERIC_GEOMETRY_CURVES_HPP
#define GENERIC_GEOMETRY_CURVES_HPP
#include "generic/math/MathUtility.hpp"
#include "Point.hpp"
#include <complex>
namespace generic {
namespace geometry{

///@brief Arc model construct by origin point, start point and the arc radian
template <typename num_type>
class Arc
{
    using float_t = float_type<num_type>;
public:
    static const size_t dim = 2;
    float_t radian;//+ccw, -cw
    Point2D<num_type> origin, start;
    /**
     * @brief constructs an Arc object
     * @param origin the origin of the arc
     * @param start the start point of the arc
     * @param radian the arc radian, the arc is CCW if radian is positive, or CW otherwise
     */
    Arc(const Point2D<num_type> & origin, const Point2D<num_type> & start, float_t radian)
     : radian(radian), origin(origin), start(start) {}

    /**
     * @brief gets the polar alpha and magnitude from origin to start point
     * @param[out] alpha the alpha of vector from origin to start point, range[0, 360)
     * @param[out] mag the magnitude of vector from origin to start point
     */
    void GetStartAlphaMag(float_t & alpha, float_t & mag) const;
    ///@brief gets the end point location
    Point2D<num_type> EndPoint() const;
};

///@brief Arc model construct by three points not on the same straight line
template <typename num_type>
class Arc3
{
    using float_t = float_type<num_type>;
public:
    static const size_t dim = 2;
    Point2D<num_type> start, mid, end;
    /**
     * @brief constructs an Arc3 object
     * @param start start point of Arc3
     * @param mid mid point of Arc3
     * @param end end point of Arc3
     */
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
