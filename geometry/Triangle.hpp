#ifndef GENERIC_GEOMETRY_TRIANGLE_HPP
#define GENERIC_GEOMETRY_TRIANGLE_HPP
#include "common/Traits.hpp"
#include "Common.hpp"
#include "Point.hpp"
#include "Box.hpp"
#include <array>
#include <cmath>
namespace generic  {
namespace geometry {
using generic::common::float_type;
template <typename num_type>
class Triangle2D
{
    using float_t = float_type<num_type>;
public:
    static const size_t dim = 2;
    using coor_t = num_type;
    using point_t = Point2D<coor_t>;
    using iterator = typename std::array<Point2D<num_type>, 3>::iterator;
    using const_iterator = typename std::array<Point2D<num_type>, 3>::const_iterator;

    Triangle2D(){}
    Triangle2D(const Point2D<num_type> & p0, const Point2D<num_type> & p1, const Point2D<num_type> & p2);

    Point2D<num_type> & operator[](size_t v);
    const Point2D<num_type> & operator[](size_t v) const;

    Box2D<num_type> BoundingBox() const;
    Point2D<float_t> Center() const;
    float_t Area() const;

    bool isCCW() const;
    void Reverse();
    WindingDirection & GetWindingDirection();
    const WindingDirection & GetWindingDirection() const;

    template<typename other_num_type>
    Triangle2D<other_num_type> Cast() const;

    iterator Begin() { return m_vertices.begin(); }
    iterator End() { return m_vertices.end(); }
    const_iterator ConstBegin() const { return m_vertices.begin(); }
    const_iterator ConstEnd() const { return m_vertices.end(); }

    static bool isCCW(const Point2D<num_type> & p1, const Point2D<num_type> & p2, const Point2D<num_type> & p3);

private:
    std::array<Point2D<num_type>, 3> m_vertices;
    mutable WindingDirection m_direction = WindingDirection::Unknown;
};

template <typename num_type>
class Triangle3D
{
    using float_t = float_type<num_type>;
public:
    static const size_t dim = 3;
    using coor_t = num_type;
    using point_t = Point3D<coor_t>;

    Triangle3D(){}
    Triangle3D(const Point3D<num_type> & p0, const Point3D<num_type> & p1, const Point3D<num_type> & p2);

    Point3D<num_type> & operator[](int v);
    const Point3D<num_type> & operator[](int v) const;

    Box3D<num_type> BoundingBox() const;
    Point3D<float_t> Center() const;
    float_t Area() const;
 
    template<typename other_num_type>
    Triangle3D<other_num_type> Cast() const;

private:
    std::array<Point3D<num_type>, 3> m_vertices;
};

template <typename num_type>
Triangle2D<num_type>::Triangle2D(const Point2D<num_type> & p0, const Point2D<num_type> & p1, const Point2D<num_type> & p2)
{ m_vertices = { p0, p1, p2 }; }

template <typename num_type>
inline Point2D<num_type> & Triangle2D<num_type>::operator[](size_t v)
{ 
    m_direction = WindingDirection::Unknown;
    return m_vertices[v]; 
}

template <typename num_type>
inline const Point2D<num_type> & Triangle2D<num_type>::operator[](size_t v) const { return m_vertices[v]; }

template <typename num_type>
inline Box2D<num_type> Triangle2D<num_type>::BoundingBox() const
{
    Box2D<num_type> bbox;
    for(size_t i = 0; i < 3; ++i)
        bbox |= m_vertices[i];
    return bbox;
}

template <typename num_type>
inline Point2D<float_type<num_type> > Triangle2D<num_type>::Center() const
{
    const float_type<num_type> inv_3 = 1 / 3.0;
    Point2D<float_type<num_type> > ct;
    for(size_t i = 0; i < 3; ++i)
        ct += m_vertices[i].template Cast<float_type<num_type> >() * inv_3;
    return ct;
}

template <typename num_type>
inline float_type<num_type> Triangle2D<num_type>::Area() const
{
    auto v1 = m_vertices[1] - m_vertices[0];
    auto v2 = m_vertices[2] - m_vertices[0];
    return std::fabs(float_type<num_type>(0.5) * v1.CrossProduct(v2));
}

template <typename num_type>
inline bool Triangle2D<num_type>::isCCW() const
{
    if(m_direction != WindingDirection::Unknown)
        return m_direction == WindingDirection::CounterClockwise;
    
    if(isCCW(m_vertices[0], m_vertices[1], m_vertices[2]))
        m_direction = WindingDirection::CounterClockwise;
    else m_direction = WindingDirection::Clockwise;
    return m_direction == WindingDirection::CounterClockwise;
}

template <typename num_type>
inline void Triangle2D<num_type>::Reverse()
{
    if(m_direction == WindingDirection::CounterClockwise) m_direction = WindingDirection::Clockwise;
    if(m_direction == WindingDirection::Clockwise) m_direction = WindingDirection::CounterClockwise;
    std::swap(m_vertices[2], m_vertices[1]);
}

template <typename num_type>
inline WindingDirection & Triangle2D<num_type>::GetWindingDirection()
{
    return m_direction;
}

template <typename num_type>
inline const WindingDirection & Triangle2D<num_type>::GetWindingDirection() const
{
    return m_direction;
}

template <typename num_type>
template<typename other_num_type>
inline Triangle2D<other_num_type> Triangle2D<num_type>::Cast() const
{
    return Triangle2D<other_num_type>(m_vertices[0].template Cast<other_num_type>(),
                                      m_vertices[1].template Cast<other_num_type>(),
                                      m_vertices[2].template Cast<other_num_type>());
}

template <typename num_type>
inline bool Triangle2D<num_type>::isCCW(const Point2D<num_type> & p1, const Point2D<num_type> & p2, const Point2D<num_type> & p3)
{
    auto v1 = p2 - p1, v2 = p3 - p1;
    return math::GT(v1.CrossProduct(v2), num_type(0));
}

template <typename num_type>
Triangle3D<num_type>::Triangle3D(const Point3D<num_type> & p0, const Point3D<num_type> & p1, const Point3D<num_type> & p2)
{ m_vertices = { p0, p1, p2 }; }

template <typename num_type>
inline Point3D<num_type> & Triangle3D<num_type>::operator[](int v) { return m_vertices[v]; }

template <typename num_type>
inline const Point3D<num_type> & Triangle3D<num_type>::operator[](int v) const { return m_vertices[v]; }

template <typename num_type>
inline Box3D<num_type> Triangle3D<num_type>::BoundingBox() const
{
    Box3D<num_type> bbox;
    for(size_t i = 0; i < 3; ++i)
        bbox |= m_vertices[i];
    return bbox;
}

template <typename num_type>
inline Point3D<float_type<num_type> > Triangle3D<num_type>::Center() const
{
    const float_type<num_type> inv_3 = 1 / 3.0;
    Point3D<float_type<num_type> > ct;
    for(size_t i = 0; i < 3; ++i)
        ct += m_vertices[i].template Cast<float_type<num_type> >() * inv_3;
    return ct;
}

template <typename num_type>
inline float_type<num_type> Triangle3D<num_type>::Area() const
{
    auto v1 = m_vertices[1] - m_vertices[0];
    auto v2 = m_vertices[2] - m_vertices[0];
    auto norm = std::sqrt(float_type<num_type>(v1.CrossProduct(v2).NormSquare()));
    return float_type<num_type>(0.5) * norm;
}

template <typename num_type>
template<typename other_num_type>
inline Triangle3D<other_num_type> Triangle3D<num_type>::Cast() const
{
    return Triangle3D<other_num_type>(m_vertices[0].template Cast<other_num_type>(),
                                      m_vertices[1].template Cast<other_num_type>(),
                                      m_vertices[2].template Cast<other_num_type>());
}

}//namespace geometry
}//namespace generic
#endif //GENERIC_GEOMETRY_TRIANGLE_HPP
