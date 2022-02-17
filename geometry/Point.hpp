#ifndef GENERIC_GEOMETRY_POINT_HPP
#define GENERIC_GEOMETRY_POINT_HPP
#include "generic/math/MathUtility.hpp"
namespace generic  {
namespace geometry {
using std::size_t;
using generic::common::float_type;
template <typename num_type>
class Point2D
{
    using float_t = float_type<num_type>;
public:
    static const size_t dim = 2;
    using coor_t = num_type;

    Point2D();
    Point2D(num_type x, num_type y);

    bool operator== (const Point2D<num_type> & p) const;
    bool operator!= (const Point2D<num_type> & p) const;

    num_type & operator[] (size_t dim);
    const num_type & operator[] (size_t dim) const;
    Point2D operator+ (const Point2D<num_type> & coor) const;
    Point2D operator- () const;
    Point2D operator- (const Point2D<num_type> & coor) const;
    Point2D operator* (float_t scale) const;
    Point2D operator/ (float_t scale) const;
    void operator+= (const Point2D<num_type> & coor);
    void operator-= (const Point2D<num_type> & coor);
    void operator*= (float_t scale);
    void operator/= (float_t scale);
    Point2D operator* (const Point2D<num_type> & a) const;

    num_type Dot(const Point2D<num_type> & a) const;
    num_type NormSquare() const;
    float_t Norm2() const;
    num_type CrossProduct(const Point2D<num_type> & a) const;

    template<typename other_num_type>
    Point2D<other_num_type> Cast() const;

    static bool isCCW(const Point2D<num_type> & p1, const Point2D<num_type> & p2, const Point2D<num_type> & p3);

private:
    num_type m_coor[2];
};

template <typename num_type>
class Point3D
{
    using float_t = float_type<num_type>;
public:
    static const size_t dim = 3;
    using coor_t = num_type;

    Point3D();
    Point3D(num_type x, num_type y, num_type z);

    bool operator== (const Point3D<num_type> & p) const;
    bool operator!= (const Point3D<num_type> & p) const;

    num_type & operator[] (size_t dim);
    const num_type & operator[] (size_t dim) const;

    Point3D operator+ (const Point3D<num_type> & coor) const;
    Point3D operator- () const;
    Point3D operator- (const Point3D<num_type> & coor) const;
    Point3D operator* (float_t scale) const;
    Point3D operator/ (float_t scale) const;
    void operator+= (const Point3D<num_type> & coor);
    void operator-= (const Point3D<num_type> & coor);
    void operator*= (float_t scale);
    void operator/= (float_t scale);
    Point3D operator* (const Point3D<num_type> & a) const;

    num_type Dot(const Point3D<num_type> & a) const;
    num_type NormSquare() const;
    float_t Norm2() const;
    Point3D CrossProduct(const Point3D<num_type> & a) const;

    template<typename other_num_type>
    Point3D<other_num_type> Cast() const;
    
private:
    num_type m_coor[3];
};

template <typename num_type>
inline Point2D<num_type>::Point2D() { m_coor[0] = 0; m_coor[1] = 0; }

template <typename num_type>
inline Point2D<num_type>::Point2D(num_type x, num_type y) { m_coor[0] = x; m_coor[1] = y; }

template <typename num_type>
inline bool Point2D<num_type>::operator== (const Point2D<num_type> & p) const
{
    return math::EQ(m_coor[0], p[0]) && math::EQ(m_coor[1], p[1]);
}

template <typename num_type>
inline bool Point2D<num_type>::operator!= (const Point2D<num_type> & p) const
{
    return !(*this == p);
}

template <typename num_type>
inline num_type & Point2D<num_type>::operator[] (size_t dim) { return m_coor[dim]; }

template <typename num_type>
inline const num_type & Point2D<num_type>::operator[] (size_t dim) const { return m_coor[dim]; }

template <typename num_type>
inline Point2D<num_type> Point2D<num_type>::operator+ (const Point2D<num_type> & coor) const 
{
    return Point2D(m_coor[0] + coor[0], m_coor[1] + coor[1]);
}

template <typename num_type>
inline Point2D<num_type> Point2D<num_type>::operator- () const
{
    return Point2D(-m_coor[0], -m_coor[1]);
}

template <typename num_type>
inline Point2D<num_type> Point2D<num_type>::operator- (const Point2D<num_type> & coor) const 
{
    return Point2D(m_coor[0] - coor[0], m_coor[1] - coor[1]);
}

template <typename num_type>
inline Point2D<num_type> Point2D<num_type>::operator* (float_t scale) const
{
    return Point2D<num_type>(scale * m_coor[0], scale * m_coor[1]);
}

template <typename num_type>
inline Point2D<num_type> Point2D<num_type>::operator/ (float_t scale) const
{
    scale = math::SafeInv(scale);
    return Point2D<num_type>(scale * m_coor[0], scale * m_coor[1]);
}

template <typename num_type>
inline void Point2D<num_type>::operator+= (const Point2D<num_type> & coor)
{
    m_coor[0] += coor[0]; m_coor[1] += coor[1];
}

template <typename num_type>
inline void Point2D<num_type>::operator-= (const Point2D<num_type> & coor)
{
    m_coor[0] -= coor[0]; m_coor[1] -= coor[1];
}

template <typename num_type>
inline void Point2D<num_type>::operator*= (float_t scale)
{
    m_coor[0] *= scale; m_coor[1] *= scale;
}

template <typename num_type>
inline void Point2D<num_type>::operator/= (float_t scale)
{
    *this *= math::SafeInv(scale);
}

template <typename num_type>
inline Point2D<num_type> Point2D<num_type>::operator* (const Point2D<num_type> & a) const
{
    return Point2D(a[0] * m_coor[0], a[1] * m_coor[1]);
}

template <typename num_type>
inline num_type Point2D<num_type>::Dot(const Point2D<num_type> & a) const
{
    return m_coor[0] * a[0] + m_coor[1] * a[1];
}

template <typename num_type>
inline num_type Point2D<num_type>::NormSquare() const
{
    return m_coor[0] * m_coor[0] + m_coor[1] * m_coor[1];
}

template <typename num_type>
inline float_type<num_type> Point2D<num_type>::Norm2() const
{
    return std::sqrt(NormSquare());
}

template <typename num_type>
inline num_type Point2D<num_type>::CrossProduct(const Point2D<num_type> & a) const
{
    return m_coor[0] * a[1] - m_coor[1] * a[0];
}

template <typename num_type>
template <typename other_num_type>
inline Point2D<other_num_type> Point2D<num_type>::Cast() const
{
    return Point2D<other_num_type>(other_num_type(m_coor[0]), other_num_type(m_coor[1]));
}

template <typename num_type>
inline bool Point2D<num_type>::isCCW(const Point2D<num_type> & p1, const Point2D<num_type> & p2, const Point2D<num_type> & p3)
{
    auto v1 = p2 - p1, v2 = p3 - p1;
    return math::GT(v1.CrossProduct(v2), num_type(0));
}

template <typename num_type>
inline Point3D<num_type>::Point3D() { m_coor[0] = 0; m_coor[1] = 0; m_coor[2] = 0; }

template <typename num_type>
inline Point3D<num_type>::Point3D(num_type x, num_type y, num_type z)
{ 
    m_coor[0] = x; m_coor[1] = y; m_coor[2] = z;
}

template <typename num_type>
inline bool Point3D<num_type>::operator== (const Point3D<num_type> & p) const
{
    return math::EQ(m_coor[0], p[0]) && math::EQ(m_coor[1], p[1]) && math::EQ(m_coor[2], p[2]);
}

template <typename num_type>
inline bool Point3D<num_type>::operator!= (const Point3D<num_type> & p) const
{
    return !(*this == p);
}

template <typename num_type>
inline num_type & Point3D<num_type>::operator[] (size_t dim) { return m_coor[dim]; }

template <typename num_type>
inline const num_type & Point3D<num_type>::operator[] (size_t dim) const { return m_coor[dim]; }

template <typename num_type>
inline Point3D<num_type> Point3D<num_type>::operator+ (const Point3D<num_type> & coor) const 
{
    return Point3D(m_coor[0] + coor[0], m_coor[1] + coor[1], m_coor[2] + coor[2]);
}

template <typename num_type>
inline Point3D<num_type> Point3D<num_type>::operator- () const
{
    return Point3D(-m_coor[0], -m_coor[1], -m_coor[2]);
}

template <typename num_type>
inline Point3D<num_type> Point3D<num_type>::operator- (const Point3D<num_type> & coor) const 
{
    return Point3D(m_coor[0] - coor[0], m_coor[1] - coor[1], m_coor[2] - coor[2]);
}

template <typename num_type>
inline Point3D<num_type> Point3D<num_type>::operator* (float_t scale) const
{
    return Point3D<num_type>(scale * m_coor[0], scale * m_coor[1], scale * m_coor[2]);
}

template <typename num_type>
inline Point3D<num_type> Point3D<num_type>::operator/ (float_t scale) const
{
    scale = math::SafeInv(scale);
    return Point3D<num_type>(scale * m_coor[0], scale * m_coor[1], scale * m_coor[2]);
}

template <typename num_type>
inline void Point3D<num_type>::operator+= (const Point3D<num_type> & coor)
{
    m_coor[0] += coor[0]; m_coor[1] += coor[1]; m_coor[2] += coor[2];
}

template <typename num_type>
inline void Point3D<num_type>::operator-= (const Point3D<num_type> & coor)
{
    m_coor[0] -= coor[0]; m_coor[1] -= coor[1]; m_coor[2] -= coor[2];
}

template <typename num_type>
inline void Point3D<num_type>::operator*= (float_t scale)
{
    m_coor[0] *= scale; m_coor[1] *= scale; m_coor[2] *= scale;
}

template <typename num_type>
inline void Point3D<num_type>::operator/= (float_t scale)
{
    *this *= math::SafeInv(scale); 
}

template <typename num_type>
inline Point3D<num_type> Point3D<num_type>::operator* (const Point3D<num_type> & a) const
{
    return Point3D(a[0] * m_coor[0], a[1] * m_coor[1], a[2] * m_coor[2]);
}

template <typename num_type>
inline num_type Point3D<num_type>::Dot(const Point3D<num_type> & a) const
{
    return m_coor[0] * a[0] + m_coor[1] * a[1] + m_coor[2] * a[2];
}

template <typename num_type>
inline num_type Point3D<num_type>::NormSquare() const
{
    return m_coor[0] * m_coor[0] + m_coor[1] * m_coor[1] + m_coor[2] * m_coor[2];
}

template <typename num_type>
inline float_type<num_type> Point3D<num_type>::Norm2() const
{
    return std::sqrt(NormSquare());
}

template <typename num_type>
inline Point3D<num_type> Point3D<num_type>::CrossProduct(const Point3D<num_type> & a) const
{
    return Point3D<num_type>(m_coor[1] * a[2] - m_coor[2] * a[1],
                             m_coor[2] * a[0] - m_coor[0] * a[2],
                             m_coor[0] * a[1] - m_coor[1] * a[0]);
}

template <typename num_type>
template <typename other_num_type>
inline Point3D<other_num_type> Point3D<num_type>::Cast() const
{
    return Point3D<other_num_type>(other_num_type(m_coor[0]), other_num_type(m_coor[1]), other_num_type(m_coor[2]));
}

}//namespace geometry
}//namespace generic
#endif //GENERIC_GEOMETRY_POINT_HPP