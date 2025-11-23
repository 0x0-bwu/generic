/**
 * @file Sphere.hpp
 * @author bwu
 * @brief Model of circle and sphere concept
 * @version 0.1
 * @date 2022-02-22
 */
#pragma once
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
    ///@brief constructs a circle at point (0, 0) with zero radius
    Circle(){ r = 0; }
    ///@brief constructs a circle at point `origin` with radius `raiuds`
    Circle(const Point2D<num_type> & origin, num_type radius);

    ///@brief checks if circle contains point `p`, not robust
    template <typename point_t>
    bool Contains(const point_t & p) const;
};

template <typename num_type>
class Ellipse
{
    using float_t = float_type<num_type>;
public:
    num_type a, b;
    Point2D<num_type> o;
    using coor_t = num_type;
    static const size_t dim = 2;
public:
    ///@brief constructs an ellipse at point (0, 0) with semi-major axis `a` and semi-minor axis `b`
    Ellipse(){ a = b = 0; }
    ///@brief constructs an ellipse at point `origin` with semi-major axis `a` and semi-minor axis `b`
    Ellipse(const Point2D<num_type> & origin, num_type a, num_type b);
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
    ///@brief constructs a sphere at point (0, 0, 0) with zero radius
    Sphere(){ r = 0; }
    ///@brief constructs a sphere at point `o` with radius `r`
    Sphere(const Point3D<num_type> & o, num_type r);
};

template <typename num_type>
Circle<num_type>::Circle(const Point2D<num_type> & origin, num_type radius){ o = origin, r = radius; }

template <typename num_type>
Ellipse<num_type>::Ellipse(const Point2D<num_type> & origin, num_type a, num_type b){ o = origin, this->a = a, this->b = b; }

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