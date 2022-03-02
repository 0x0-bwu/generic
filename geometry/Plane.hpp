/**
 * @file Plane.hpp
 * @author bwu
 * @brief Model of plane concept
 * @version 0.1
 * @date 2022-02-22
 */
#ifndef GENERIC_GEOMETRY_PLANE_HPP
#define GENERIC_GEOMETRY_PLANE_HPP
#include "generic/common/Traits.hpp"
#include "Point.hpp"
#include "Vector.hpp"
namespace generic  {
namespace geometry {
using generic::common::float_type;
template <typename num_type>
class Plane
{
    using float_t = float_type<num_type>;
public:
    ///@brief constructs a plane with the normal vector and through point
    Plane(const Vector3D<num_type> & normal, const Point3D<num_type> & p);
    ///@brief constructs a plane with three points that not on the same straight line
    Plane(const Point3D<num_type> & p1, const Point3D<num_type> & p2, const Point3D<num_type> & p3);

    ///@brief gets plane normal, points x on the plane satisfy Dot(n, x) = d, d is the result of D()
    Vector3D<float_t> Normal() const { return m_normal; }
    ///@brief the distance of this plane from the origin
    float_t D() const { return m_dot; }

private:
    Vector3D<float_t> m_normal;
    float_t m_dot;
};

template <typename num_type
inline Plane<num_type>::Plane(const Vector3D<num_type> & normal, const Point3D<num_type> & p)
{
    m_normal = Normalize(normal);
    m_dot = DotProduct(m_normal, p.template Cast<float_t>());
}

template <typename num_type>
inline Plane<num_type>::Plane(const Point3D<num_type> & p1, const Point3D<num_type> & p2, const Point3D<num_type> & p3)
{
    m_normal = Normalize((p2 - p1).CrossProduct(p3 - p1));
    m_dot = m_normal.Dot(p1.template Cast<float_t>());
}
}//namespace geometry
}//namespace generic
#endif//GENERIC_GEOMETRY_PLANE_HPP