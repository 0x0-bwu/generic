#ifndef GENERIC_GEOMETRY_SPHERE_HPP
#define GENERIC_GEOMETRY_SPHERE_HPP
#include "generic/math/MathUtility.hpp"
#include "Common.hpp"
#include "Point.hpp"
namespace generic {
namespace geometry{
template <typename num_type>
class Circle
{
    using float_t = float_type<num_type>;
public:
    num_type r;
    Point2D<num_type> o;
    using coor_t = num_type;
    static const size_t dim = 2;

public:
    Circle(){ r = 0; }
    Circle(const Point2D<num_type> & origin, num_type radius);

    template <typename point_t>
    bool Contains(const point_t & p) const;

    ///@brief get circumcircle of three points that not in common line
    static Circle<float_t> CircumCircle(const Point2D<num_type> & p1, const Point2D<num_type> & p2, const Point2D<num_type> & p3);
    static Point2D<float_t> CircumCircle(const Point2D<num_type> & p1, const Point2D<num_type> & p2, const Point2D<num_type> & p3, float_t & radius2);
};

template <typename num_type>
class Sphere
{
public:
    num_type r;
    Point3D<num_type> o;
    using coor_t = num_type;
    static const size_t dim = 3;

public:
    Sphere(){ r = 0; }
    Sphere(const Point3D<num_type> & o, num_type r);
};

template <typename num_type>
Circle<num_type>::Circle(const Point2D<num_type> & origin, num_type radius){ o = origin, r = radius; }

template <typename num_type>
template <typename point_t>
bool Circle<num_type>::Contains(const point_t & p) const
{
    auto distSq = (p[0] - o[0]) * (p[0] - o[0]) + (p[1] - o[1]) * (p[1] - o[1]);
    return !math::GT(distSq, r * r);
}

template <typename num_type>
inline Circle<float_type<num_type> > Circle<num_type>::CircumCircle(const Point2D<num_type> & p1, const Point2D<num_type> & p2, const Point2D<num_type> & p3)
{
    float_type<num_type> r2(0);
    Circle<float_type<num_type> > circle;
    circle.o = CircumCircle(p1, p2, p3, r2);
    circle.r = std::sqrt(r2);
    return circle;
}

template <typename num_type>//return coor of circle
inline Point2D<float_type<num_type> > Circle<num_type>::CircumCircle(const Point2D<num_type> & p1, const Point2D<num_type> & p2, const Point2D<num_type> & p3, float_type<num_type> & radius2)
{   
    using float_t = float_type<num_type>;
    float_t epsilon = std::numeric_limits<float_t>::epsilon();
    float_t a1 = p2[1] - p1[1];
    float_t a2 = p3[1] * p3[1] - p1[1] * p1[1] + p3[0] * p3[0] - p1[0] * p1[0];
    float_t b1 = p3[1] - p1[1];
    float_t b2 = p2[1] * p2[1] - p1[1] * p1[1] + p2[0] * p2[0] - p1[0] * p1[0];
    float_t c = 2.0 * ((p3[0] - p1[0]) * (p2[1] - p1[1]) - (p2[0] - p1[0]) * (p3[1] - p1[1]));
    if(std::fabs(c) < epsilon) c = std::copysign(epsilon, c);
    float_t x = (a1 * a2 - b1 * b2) / c;

    float_t d1 = p2[0] - p1[0];
    float_t d2 = p3[0] * p3[0] - p1[0] * p1[0] + p3[1] * p3[1] - p1[1] * p1[1];
    float_t e1 = p3[0] - p1[0];
    float_t e2 = p2[0] * p2[0] - p1[0] * p1[0] + p2[1] * p2[1] - p1[1] * p1[1];
    float_t f = 2.0 * ((p3[1] - p1[1]) * (p2[0] - p1[0]) - (p2[1] - p1[1]) * (p3[0] - p1[0]));
    if(std::fabs(f) < epsilon) f = std::copysign(epsilon, f);
    float_t y = (d1 * d2 - e1 * e2) / f;

    radius2 = float_t((x - p1[0]) * (x - p1[0]) + (y - p1[1]) * (y - p1[1]));
    return Point2D<float_t>(x, y);
}

template <typename num_type>
Sphere<num_type>::Sphere(const Point3D<num_type> & origin, num_type radius){ o = origin, r = radius; }
    
}//namespace geometry
}//namespace generic
#endif//GENERIC_GEOMETRY_SPHERE_HPP
