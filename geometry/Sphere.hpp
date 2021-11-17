#ifndef GENERIC_GEOMETRY_SPHERE_HPP
#define GENERIC_GEOMETRY_SPHERE_HPP
#include "math/MathUtility.hpp"
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
Sphere<num_type>::Sphere(const Point3D<num_type> & origin, num_type radius){ o = origin, r = radius; }
    
}//namespace geometry
}//namespace generic
#endif//GENERIC_GEOMETRY_SPHERE_HPP
