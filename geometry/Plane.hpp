#ifndef GENERIC_GEOMETRY_PLANE_HPP
#define GENERIC_GEOMETRY_PLANE_HPP
#include "common/Traits.hpp"
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
    Plane(const Point3D<num_type> & p1, const Point3D<num_type> & p2, const Point3D<num_type> & p3);

    Vector3D<float_t> Normal() const { return m_normal; }
    float_t D() const { return m_dot; }

private:
    Vector3D<float_t> m_normal;
    float_t m_dot;
};

template <typename num_type>
inline Plane<num_type>::Plane(const Point3D<num_type> & p1, const Point3D<num_type> & p2, const Point3D<num_type> & p3)
{
    m_normal = Normalize((p2 - p1).CrossProduct(p3 - p1));
    m_dot = m_normal.Dot(p1.template Cast<float_t>());
}
}//namespace geometry
}//namespace generic
#endif//GENERIC_GEOMETRY_PLANE_HPP