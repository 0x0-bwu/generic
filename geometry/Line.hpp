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
    Line2D(){}
    Line2D(float_t _a, float_t _b, float_t _c);

    bool operator== (const Line2D<num_type> & l) const;
    bool operator!= (const Line2D<num_type> & l) const;

    void Normalize();
    bool isValid() const;
    static bool Intersects(const Line2D<num_type> & line1, const Line2D<num_type> & line2, Point2D<float_t> & point);

    template <typename point_t, typename std::enable_if<point_t::dim == 2, bool>::type = true>
    static float_t SideValue(const Line2D<num_type> & line, const point_t & point);

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