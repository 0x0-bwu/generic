/**
 * @file Line.hpp
 * @author bwu
 * @brief Model of line concept
 * @version 0.1
 * @date 2022-02-14
 */
#ifndef GENERIC_GEOMETRY_LINE_HPP
#define GENERIC_GEOMETRY_LINE_HPP
#include "Segment.hpp"
#include "Vector.hpp"
namespace generic  {
namespace geometry {
using generic::common::float_type;

template <typename num_type>
class Line2D
{
    using float_t = float_type<num_type>;
public:
    float_t a = 0;
    float_t b = 0;
    float_t c = 0;
    bool normalized = false;
    static const size_t dim = 2;
    ///@brief constructs an invalid line2d
    Line2D(){}
    ///@brief constructs a line with coefficient a, b and c
    Line2D(float_t _a, float_t _b, float_t _c);

    ///@brief checks if this line equals to line l
    bool operator== (const Line2D<num_type> & l) const;
    ///@brief checks if this line not equals to line l
    bool operator!= (const Line2D<num_type> & l) const;

    ///@brief normalizes the coefficient a, b, and c, make a^2 + b^2 + c^2 = 1
    void Normalize();
    ///@brief checks if the line is valid(a and b not equal to zero at same time)
    bool isValid() const;

    /**
     * @brief checks if two lines intersected and get the intersect point
     * @param[in] line1 one of the input line
     * @param[in] line2 one of the input line
     * @param[in] point the intersect location if two line intersected.
     * @return whether two lines intersected 
     */
    static bool Intersects(const Line2D<num_type> & line1, const Line2D<num_type> & line2, Point2D<float_t> & point);

    /**
     * @brief gets side value of line and point, usually use the sign to check the point line location, quick than Distance(line, point)
     * @param[in] line the input line
     * @param[in] point the input line
     * @return side value of point to line
     */
    template <typename point_t, typename std::enable_if<point_t::dim == 2, bool>::type = true>
    static float_t SideValue(const Line2D<num_type> & line, const point_t & point);

    /**
     * @brief gets distance of line and point
     * @param[in] line the input line
     * @param[in] point the input point
     * @return distance of point to line
     */
    template <typename point_t, typename std::enable_if<point_t::dim == 2, bool>::type = true>
    static float_t Distance(const Line2D<num_type> & line, const point_t & point);
};

template <typename num_type>
inline Line2D<num_type>::Line2D(float_t _a, float_t _b, float_t _c) : a(_a), b(_b), c(_c) {}

template <typename num_type>
inline bool Line2D<num_type>::operator== (const Line2D<num_type> & l) const
{
    auto norm = Vector3D<float_t>(a, b, c).CrossProduct(Vector3D<float_t>(l.a, l.b, l.c));
    return norm == Vector3D<float_t>(0, 0, 0);
}

template <typename num_type>
inline bool Line2D<num_type>::operator!= (const Line2D<num_type> & l) const
{
    return !(*this == l);
}

template <typename num_type>
inline void Line2D<num_type>::Normalize()
{
    const float_t norm_inv = 1 / std::sqrt(a * a + b * b);
    a *= norm_inv; b *= norm_inv; c *= norm_inv;
    normalized = true;
}

template <typename num_type>
inline bool Line2D<num_type>::isValid() const
{
    return !(math::EQ(a, float_t(0)) && math::EQ(b, float_t(0)));
}

template <typename num_type>
inline bool Line2D<num_type>::Intersects(const Line2D<num_type> & line1, const Line2D<num_type> & line2, Point2D<float_t> & point)
{
    const float_t denominator = line1.b * line2.a - line1.a * line2.b;
    if(math::EQ(denominator, float_t(0))) return false;

    point[0] = (line1.c * line2.b - line1.b * line2.c) / denominator;
    point[1] = (line1.a * line2.c - line1.c * line2.a) / denominator;
    return true;
}

template <typename num_type>
template <typename point_t, typename std::enable_if<point_t::dim == 2, bool>::type>
inline float_type<num_type> Line2D<num_type>::SideValue(const Line2D<num_type> & line, const point_t & point)
{
    return line.a * point[0] + line.b * point[1] + line.c;
}

template <typename num_type>
template <typename point_t, typename std::enable_if<point_t::dim == 2, bool>::type>
inline float_type<num_type> Line2D<num_type>::Distance(const Line2D<num_type> & line, const point_t & point)
{
    return SideValue(line, point) / std::sqrt(line.a * line.a + line.b * line.b);
}

template <typename num_type>
inline Line2D<num_type> makeLineByPointAndNormal(const Point2D<num_type> & p, const Vector2D<num_type> & n)
{
    float_type<num_type> a = n[0];
    float_type<num_type> b = n[1];
    float_type<num_type> c = -n[0] * p[0] - n[1] * p[1];
    return Line2D<num_type>(a, b, c);
}

template <typename num_type>
inline Line2D<num_type> makeLineByTwoPoints(const Point2D<num_type> & p1, const Point2D<num_type> & p2)
{
    Vector2D<num_type> vec = p2 - p1;
    Vector2D<num_type> norm(-vec[1], vec[0]);
    return makeLineByPointAndNormal(p1, norm);
}

template <typename num_type>
inline Line2D<num_type> makeLineBySegment(const Segment2D<num_type> & s)
{
    return makeLineByTwoPoints(s[0], s[1]);
}

}//namespace geometry
}//namespace generic
#endif//GENERIC_GEOMETRY_LINE_HPP