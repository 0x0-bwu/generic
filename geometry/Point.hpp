/**
 * @file Point.hpp
 * @author bwu
 * @brief Model of point2d and point3d concept
 * @version 0.1
 * @date 2022-02-22
 */
#pragma once
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

    ///@brief constructs a point2d(0, 0)
    Point2D();
    ///@brief constructs a point2d(x, y)
    Point2D(num_type x, num_type y);

    ///@brief check if this point equals to p
    bool operator== (const Point2D<num_type> & p) const;
    ///@brief check if this point not equals to p
    bool operator!= (const Point2D<num_type> & p) const;

    ///@brief access the x, y coordinate by index 0-1
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

    ///@brief dot pruduct of this point with a
    num_type Dot(const Point2D<num_type> & a) const;
    ///@brief squared norm2 value
    num_type NormSquare() const;
    ///@brief norm2 vaule
    float_t Norm2() const;
    ///@brief cross product of this point with a
    num_type CrossProduct(const Point2D<num_type> & a) const;

    ///@brief converts to point with other number type explicitly
    template<typename other_num_type>
    Point2D<other_num_type> Cast() const;

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

    ///@brief constructs a point3d(0, 0, 0)
    Point3D();
    ///@brief constructs a point3d(x, y, z)
    Point3D(num_type x, num_type y, num_type z);

    ///@brief check if this point equals to p
    bool operator== (const Point3D<num_type> & p) const;
    ///@brief check if this point not equals to p
    bool operator!= (const Point3D<num_type> & p) const;

    ///@brief access the x, y, zcoordinate by index 0-2
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

    ///@brief dot pruduct of this point with a
    num_type Dot(const Point3D<num_type> & a) const;
    ///@brief squared norm2 value
    num_type NormSquare() const;
    ///@brief norm2 vaule
    float_t Norm2() const;
    ///@brief cross product of this point with a
    Point3D CrossProduct(const Point3D<num_type> & a) const;

    ///@brief converts to point with other number type explicitly
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
    if constexpr (std::is_integral_v<num_type>)
        return Point2D<num_type>(std::round(scale * m_coor[0]), std::round(scale * m_coor[1]));
    else return Point2D<num_type>(scale * m_coor[0], scale * m_coor[1]);
}

template <typename num_type>
inline Point2D<num_type> Point2D<num_type>::operator/ (float_t scale) const
{
    return *this * math::SafeInv(scale);
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
    *this = *this * scale;
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
    if constexpr (std::is_floating_point_v<num_type> && std::is_integral_v<other_num_type>)
        return Point2D<other_num_type>(std::round(m_coor[0]), std::round(m_coor[1]));
    else return Point2D<other_num_type>(other_num_type(m_coor[0]), other_num_type(m_coor[1]));
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
    if constexpr (std::is_integral_v<num_type>)
        return Point3D<num_type>(std::round(scale * m_coor[0]), std::round(scale * m_coor[1]), std::round(scale * m_coor[2]));
    else return Point3D<num_type>(scale * m_coor[0], scale * m_coor[1], scale * m_coor[2]);
}

template <typename num_type>
inline Point3D<num_type> Point3D<num_type>::operator/ (float_t scale) const
{
    return *this * math::SafeInv(scale);
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
    *this = *this * scale;
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
    if constexpr (std::is_floating_point_v<num_type> && std::is_integral_v<other_num_type>)
        return Point3D<other_num_type>(std::round(m_coor[0]), std::round(m_coor[1]), std::round(m_coor[2]));
    else return Point3D<other_num_type>(other_num_type(m_coor[0]), other_num_type(m_coor[1]), other_num_type(m_coor[2]));
}

}//namespace geometry
}//namespace generic
