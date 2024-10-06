/**
 * @file Box.hpp
 * @author bwu
 * @brief Model of axis-aligned bounding box2d and box3d concept
 * @version 0.1
 * @date 2022-02-22
 */
#pragma once
#include "Vector.hpp"
namespace generic  {
namespace geometry {
using generic::common::float_type;

template <typename num_type>
class Box2D
{
    using float_t = float_type<num_type>;
public:
    const static size_t dim = 2;
    using coor_t = num_type;
    using point_t = Point2D<coor_t>;
    
    ///@brief constructs an invalid box2d
    Box2D();
    ///@brief constructs a box2d by LL point and UR point
    Box2D(num_type llx, num_type lly, num_type urx, num_type ury);
    ///@brief constructs a box2d formed by two points
    Box2D(const Point2D<num_type> & p1, const Point2D<num_type> & p2);

    ///@brief checks if this box equal to box b
    bool operator== (const Box2D & b) const;
    ///@brief checks if this box not equal to box b
    bool operator!= (const Box2D & b) const;

    ///@brief accesses LL and UR point by index 0 and 1
    Point2D<num_type> & operator[] (size_t i);
    const Point2D<num_type> & operator[] (size_t i) const;

    ///@brief gets scaled box by factor `scale`
    Box2D<num_type> operator* (float_t scale) const;
    ///@brief gets scaled box by factor 1 / `scale`
    Box2D<num_type> operator/ (float_t scale) const;

    ///@brief scalew box by factor `scale`
    void operator*= (float_t scale);
    ///@brief scalew box by factor 1 / `scale`
    void operator/= (float_t scale);

    ///@brief checks if this box is inside box b
    bool operator<  (const Box2D<num_type> & b) const;
    ///@brief checks if box b is inside this box
    bool operator>  (const Box2D<num_type> & b) const;
    ///@brief checks if point p is inside this box
    bool operator>  (const Point2D<num_type> & p) const;
    ///@brief checks if this box is inside box b, consider touch
    bool operator<= (const Box2D<num_type> & b) const;
    ///@brief checks if box b is inside this box, consider touch
    bool operator>= (const Box2D<num_type> & b) const;
    ///@brief checks if point p is inside this box, consider touch
    bool operator>= (const Point2D<num_type> & p) const;
    ///@brief gets union box of this box with box b
    Box2D<num_type> operator+ (const Box2D<num_type> & b) const;
    ///@brief self union with box b
    void operator+= (const Box2D<num_type> & b);
    ///@brief same as operator+ box b
    void operator|= (const Box2D<num_type> & b);
    ///@brief gets union box of this box with point p
    Box2D<num_type> operator| (const Point2D<num_type> & p) const;
    ///@brief same as operator+
    void operator|= (const Point2D<num_type> & p);
    ///@brief gets intersect box of this box with box b
    Box2D<num_type> operator& (const Box2D<num_type> & b) const;
    ///@brief self intersection with box b
    void operator&= (const Box2D<num_type> & b);

    ///@brief shift self by vector v
    void operator+= (const Vector2D<num_type> & v);
    ///#brief shift self by vector -v
    void operator-= (const Vector2D<num_type> & v);
    ///@brief shift box by vector v
    Box2D<num_type> operator+ (const Vector2D<num_type> & v) const;
    ///@brief shift box by vector -v
    Box2D<num_type> operator- (const Vector2D<num_type> & v) const;

    ///@brief gets box center
    Point2D<float_t> Center() const;
    ///@brief gets box length in axis-x
    num_type Length() const;
    ///@brief gets box width in axis-y
    num_type Width() const;
    ///@brief gets box area
    num_type Area() const;
    ///@brief checks if LL point is under UR point
    bool isValid() const;
    ///@brief reorders the LL point and UR point location
    void Normalize();
    ///@brief reverses LL and UR point location of a valid box
    void SetInvalid();
    /// @brief scale the box from center
    void Scale(float_t scale);

    ///@brief converts this box to box with other number type explicitly
    template<typename other_num_type>
    Box2D<other_num_type> Cast() const;

    ///@brief collision test of two boxes
    static bool Collision(const Box2D<num_type> & a, const Box2D<num_type> & b, bool considerTouch = false);

private:
    static bool Collision(const Box2D<num_type> & a, const Box2D<num_type> & b, std::true_type);//consider touch
    static bool Collision(const Box2D<num_type> & a, const Box2D<num_type> & b, std::false_type);//not consider touch

private:
    Point2D<num_type> m_corner[2];
};

template <typename num_type>
class Box3D
{
    using float_t = common::float_type<num_type>;
public:
    const static size_t dim = 3;
    using coor_t = num_type;
    using point_t = Point3D<coor_t>;

    ///@brief constructs an invalid box3d
    Box3D();
    ///@brief constructs a box3d by LL point and UR point
    Box3D(num_type llx, num_type lly, num_type llz, num_type urx, num_type ury, num_type urz);
    ///@brief constructs a box3d formed by two points
    Box3D(const Point3D<num_type> & p1, const Point3D<num_type> & p2);

    ///@brief checks if this box equal to box b
    bool operator== (const Box3D & b) const;
    ///@brief checks if this box not equal to box b
    bool operator!= (const Box3D & b) const;

    ///@brief accesses LL and UR point by index 0 and 1
    Point3D<num_type> & operator[] (size_t i);
    const Point3D<num_type> & operator[] (size_t i) const;

    ///@brief gets scaled box by factor `scale`
    Box3D<num_type> operator* (float_t scale) const;
    ///@brief gets scaled box by factor 1 / `scale`
    Box3D<num_type> operator/ (float_t scale) const;

    ///@brief scalew box by factor `scale`
    void operator*= (float_t scale);
    ///@brief scalew box by factor 1 / `scale`
    void operator/= (float_t scale);

    ///@brief checks if this box is inside box b
    bool operator<  (const Box3D<num_type> & b) const;
    ///@brief checks if box b is inside this box
    bool operator>  (const Box3D<num_type> & b) const;
    ///@brief checks if point p is inside this box
    bool operator>  (const Point3D<num_type> & p) const;
    ///@brief checks if this box is inside box b, consider touch
    bool operator<= (const Box3D<num_type> & b) const;
    ///@brief checks if box b is inside this box, consider touch
    bool operator>= (const Box3D<num_type> & b) const;
    ///@brief checks if point p is inside this box, consider touch
    bool operator>= (const Point3D<num_type> & p) const;
    ///@brief gets union box of this box with box b
    Box3D<num_type> operator+ (const Box3D<num_type> & b) const;
    ///@brief self union with box b
    void operator+= (const Box3D<num_type> & b);
    ///@brief same as operator+
    void operator|= (const Box3D<num_type> & b);
    ///@brief self union with point p
    void operator|= (const Point3D<num_type> & p);
    ///@brief gets intersect box of this box with box b
    Box3D<num_type> operator& (const Box3D<num_type> & b) const;
    ///@brief self intersection with box b
    void operator&= (const Box3D<num_type> & b);

    ///@brief shift self by vector v
    void operator+= (const Vector3D<num_type> & v);
    ///#brief shift self by vector -v
    void operator-= (const Vector3D<num_type> & v);
    ///@brief shift box by vector v
    Box3D<num_type> operator+ (const Vector3D<num_type> & v) const;
    ///@brief shift box by vector -v
    Box3D<num_type> operator- (const Vector3D<num_type> & v) const;

    ///@brief gets box center
    Point3D<float_t> Center() const;
    ///@brief gets box length in axis-x
    num_type Length() const;
    ///@brief gets box width in axis-y
    num_type Width() const;
    ///@brief gets box height in axis-z
    num_type Height() const;
    ///@brief gets diagonal vector from LL to UR
    Point3D<num_type> Diagonal() const;
    ///@brief gets half of the surface area of this box
    num_type HalfArea() const;
    ///@brief gets surface area of this box
    num_type SurfArea() const;
    ///@brief gets volume of this box
    num_type Volume() const;
    ///@brief gets axis index with largest range
    size_t LargestAxis() const;
    ///@brief checks if LL point is under UR point
    bool isValid() const;
    ///@brief reorders the LL point and UR point location
    void Normalize();
    ///@brief reverses LL and UR point location of a valid box
    void SetInvalid();
    /// @brief scale the box from center
    void Scale(float_t scale);

    ///@brief converts this box to box with other number type explicitly
    template<typename other_num_type>
    Box3D<other_num_type> Cast() const;

    ///@brief collision test of two boxes
    static bool Collision(const Box3D<num_type> & a, const Box3D<num_type> & b, bool considerTouch = false);

private:
    static bool Collision(const Box3D<num_type> & a, const Box3D<num_type> & b, std::true_type);//consider touch
    static bool Collision(const Box3D<num_type> & a, const Box3D<num_type> & b, std::false_type);//not consider touch

private:
    Point3D<num_type> m_corner[2];
};

template <typename num_type>
inline Box2D<num_type>::Box2D()
{
    num_type max = std::numeric_limits<num_type>::max();
    m_corner[0] = Point2D<num_type>(max, max);
    m_corner[1] = Point2D<num_type>(-max, -max);
}

template <typename num_type>
inline Box2D<num_type>::Box2D(num_type llx, num_type lly, num_type urx, num_type ury)
{
    m_corner[0] = Point2D<num_type>(llx, lly);
    m_corner[1] = Point2D<num_type>(urx, ury);
    Normalize();
}

template <typename num_type>
inline Box2D<num_type>::Box2D(const Point2D<num_type> & p1, const Point2D<num_type> & p2)
{
    m_corner[0] = p1;
    m_corner[1] = p2;
    Normalize();
}

template <typename num_type>
inline bool Box2D<num_type>::operator== (const Box2D & b) const
{
    return m_corner[0] == b[0] && m_corner[1] == b[1];
}

template <typename num_type>
inline bool Box2D<num_type>::operator!= (const Box2D & b) const
{
    return !(*this == b);
}

template <typename num_type>
inline Point2D<num_type> & Box2D<num_type>::operator[] (size_t i) { return m_corner[i]; }

template <typename num_type>
inline const Point2D<num_type> & Box2D<num_type>::operator[] (size_t i) const { return m_corner[i]; }

template <typename num_type>
inline Box2D<num_type> Box2D<num_type>::operator* (float_t scale) const
{
    auto box = *this;
    box *= scale;
    return box;
}

template <typename num_type>
inline Box2D<num_type> Box2D<num_type>::operator/ (float_t scale) const
{
    auto box = *this;
    box /= scale;
    return box;
}

template <typename num_type>
inline void Box2D<num_type>::operator*= (float_t scale)
{
    m_corner[0] *= scale;
    m_corner[1] *= scale;
}

template <typename num_type>
inline void Box2D<num_type>::operator/= (float_t scale)
{
    m_corner[0] /= scale;
    m_corner[1] /= scale;   
}

template <typename num_type>
inline bool Box2D<num_type>::operator< (const Box2D<num_type> & b) const
{
    return math::GT(m_corner[0][0], b[0][0]) && math::GT(m_corner[0][1], b[0][1]) &&
           math::LT(m_corner[1][0], b[1][0]) && math::LT(m_corner[1][1], b[1][1]);
}

template <typename num_type>
inline bool Box2D<num_type>::operator> (const Box2D<num_type> & b) const
{
    return math::LT(m_corner[0][0], b[0][0]) && math::LT(m_corner[0][1], b[0][1]) &&
           math::GT(m_corner[1][0], b[1][0]) && math::GT(m_corner[1][1], b[1][1]);
}

template <typename num_type>
inline bool Box2D<num_type>::operator> (const Point2D<num_type> & p) const
{
    return math::LT(m_corner[0][0], p[0]) && math::LT(m_corner[0][1], p[1]) &&
           math::GT(m_corner[1][0], p[0]) && math::GT(m_corner[1][1], p[1]);
}

template <typename num_type>
inline bool Box2D<num_type>::operator<= (const Box2D<num_type> & b) const
{
    return math::GE(m_corner[0][0], b[0][0]) && math::GE(m_corner[0][1], b[0][1]) &&
           math::LE(m_corner[1][0], b[1][0]) && math::LE(m_corner[1][1], b[1][1]);
}

template <typename num_type>
inline bool Box2D<num_type>::operator>= (const Box2D<num_type> & b) const
{
    return math::LE(m_corner[0][0], b[0][0]) && math::LE(m_corner[0][1], b[0][1]) &&
           math::GE(m_corner[1][0], b[1][0]) && math::GE(m_corner[1][1], b[1][1]);
}

template <typename num_type>
inline bool Box2D<num_type>::operator>= (const Point2D<num_type> & p) const
{
    return math::LE(m_corner[0][0], p[0]) && math::LE(m_corner[0][1], p[1]) &&
           math::GE(m_corner[1][0], p[0]) && math::GE(m_corner[1][1], p[1]);
}

template <typename num_type>
inline Box2D<num_type> Box2D<num_type>::operator+ (const Box2D<num_type> & b) const
{
    Box2D<num_type> box;
    box[0][0] = std::min(m_corner[0][0], b[0][0]);
    box[0][1] = std::min(m_corner[0][1], b[0][1]);
    box[1][0] = std::max(m_corner[1][0], b[1][0]);
    box[1][1] = std::max(m_corner[1][1], b[1][1]);
    return box;
}

template <typename num_type>
inline void Box2D<num_type>::operator+= (const Box2D<num_type> & b)
{
    (*this) = (*this) + b;
}

template <typename num_type>
inline void Box2D<num_type>::operator|= (const Box2D<num_type> &b)
{
    *this += b;
}

template <typename num_type>
inline Box2D<num_type> Box2D<num_type>::operator| (const Point2D<num_type> & p) const
{
    Box2D<num_type> box;
    box[0][0] = std::min(m_corner[0][0], p[0]);
    box[0][1] = std::min(m_corner[0][1], p[1]);
    box[1][0] = std::max(m_corner[1][0], p[0]);
    box[1][1] = std::max(m_corner[1][1], p[1]);
    return box;
}

template <typename num_type>
inline void Box2D<num_type>::operator|= (const Point2D<num_type> & p)
{
    m_corner[0][0] = std::min(m_corner[0][0], p[0]);
    m_corner[0][1] = std::min(m_corner[0][1], p[1]);
    m_corner[1][0] = std::max(m_corner[1][0], p[0]);
    m_corner[1][1] = std::max(m_corner[1][1], p[1]);
}

template <typename num_type>
inline Box2D<num_type> Box2D<num_type>::operator& (const Box2D<num_type> & b) const
{
    Box2D<num_type> box;
    box[0][0] = std::max(m_corner[0][0], b[0][0]);
    box[0][1] = std::max(m_corner[0][1], b[0][1]);
    box[1][0] = std::min(m_corner[1][0], b[1][0]);
    box[1][1] = std::min(m_corner[1][1], b[1][1]);
    return box;
}

template <typename num_type>
inline void Box2D<num_type>::operator&= (const Box2D<num_type> & b)
{
    (*this) = (*this) & b;
}

template <typename num_type>
inline void Box2D<num_type>::operator+= (const Vector2D<num_type> & v)
{
    m_corner[0] += v;
    m_corner[1] += v;
}

template <typename num_type>
inline void Box2D<num_type>::operator-= (const Vector2D<num_type> & v)
{
    m_corner[0] -= v;
    m_corner[1] -= v;
}

template <typename num_type>
inline Box2D<num_type> Box2D<num_type>::operator+ (const Vector2D<num_type> & v) const
{
    auto box = *this;
    box += v;
    return box;
}

template <typename num_type>
inline Box2D<num_type> Box2D<num_type>::operator- (const Vector2D<num_type> & v) const
{
    auto box = *this;
    box -= v;
    return box;
}

template <typename num_type>
inline Point2D<float_type<num_type> > Box2D<num_type>::Center() const
{
    float_t x = 0.5 * m_corner[0][0] + 0.5 * m_corner[1][0];
    float_t y = 0.5 * m_corner[0][1] + 0.5 * m_corner[1][1];
    return Point2D<float_t>(x, y);
}

template <typename num_type>
inline num_type Box2D<num_type>::Length() const
{
    return m_corner[1][0] - m_corner[0][0];
}

template <typename num_type>
inline num_type Box2D<num_type>::Width() const
{
    return m_corner[1][1] - m_corner[0][1];
}

template <typename num_type>
inline num_type Box2D<num_type>::Area() const
{
    return Width() * Length();
}

template <typename num_type>
inline bool Box2D<num_type>::isValid() const
{
    return math::LE(m_corner[0][0], m_corner[1][0]) && math::LE(m_corner[0][1], m_corner[1][1]);
}

template <typename num_type>
inline void Box2D<num_type>::Normalize()
{
    for(size_t i = 0; i < 2; ++i){
        num_type min = std::min(m_corner[0][i], m_corner[1][i]);
        num_type max = std::max(m_corner[0][i], m_corner[1][i]);
        m_corner[0][i] = min; m_corner[1][i] = max;
    }
}

template <typename num_type>
inline void Box2D<num_type>::SetInvalid()
{
    *this = Box2D();
}

template <typename num_type>
inline void Box2D<num_type>::Scale(float_t scale)
{
    auto ct = Center().template Cast<num_type>();
    *this -= ct;
    *this *= scale;
    *this += ct;
}

template <typename num_type>
template<typename other_num_type>
inline Box2D<other_num_type> Box2D<num_type>::Cast() const
{
    Box2D<other_num_type> box;
    box[0] = m_corner[0].template Cast<other_num_type>();
    box[1] = m_corner[1].template Cast<other_num_type>();
    return box;
}

template <typename num_type>
inline bool Box2D<num_type>::Collision(const Box2D<num_type> & a, const Box2D<num_type> & b, bool considerTouch)
{
    if(considerTouch)
        return Collision(a, b, std::true_type{});
    else return Collision(a, b, std::false_type{});
}

template <typename num_type>
inline bool Box2D<num_type>::Collision(const Box2D<num_type> & a, const Box2D<num_type> & b, std::true_type)
{
    if(math::GE(a[1][0], b[0][0]) && math::LE(a[0][0], b[1][0]) &&
       math::GE(a[1][1], b[0][1]) && math::LE(a[0][1], b[1][1]))
        return true;
    return false;
}

template <typename num_type>
inline bool Box2D<num_type>::Collision(const Box2D<num_type> & a, const Box2D<num_type> & b, std::false_type)
{
    if(math::GT(a[1][0], b[0][0]) && math::LT(a[0][0], b[1][0]) &&
       math::GT(a[1][1], b[0][1]) && math::LT(a[0][1], b[1][1]))
        return true;
    return false;
}

template <typename num_type>
inline Box3D<num_type>::Box3D()
{
    num_type max = std::numeric_limits<num_type>::max();
    m_corner[0] = Point3D<num_type>( max,  max,  max);
    m_corner[1] = Point3D<num_type>(-max, -max, -max);
}

template <typename num_type>
inline Box3D<num_type>::Box3D(num_type llx, num_type lly, num_type llz, num_type urx, num_type ury, num_type urz)
{
    m_corner[0] = Point3D<num_type>(llx, lly, llz);
    m_corner[1] = Point3D<num_type>(urx, ury, urz);
    Normalize();
}

template <typename num_type>
inline Box3D<num_type>::Box3D(const Point3D<num_type> & p1, const Point3D<num_type> & p2)
{
    m_corner[0] = p1;
    m_corner[1] = p2;
    Normalize();
}

template <typename num_type>
inline bool Box3D<num_type>::operator== (const Box3D & b) const
{
    return m_corner[0] == b[0] && m_corner[1] == b[1];
}

template <typename num_type>
inline bool Box3D<num_type>::operator!= (const Box3D & b) const
{
    return !(*this == b);
}

template <typename num_type>
inline Point3D<num_type> & Box3D<num_type>::operator[] (size_t i) { return m_corner[i]; }

template <typename num_type>
inline const Point3D<num_type> & Box3D<num_type>::operator[] (size_t i) const { return m_corner[i]; }

template <typename num_type>
inline Box3D<num_type> Box3D<num_type>::operator* (float_t scale) const
{
    auto box = *this;
    box *= scale;
    return box;
}

template <typename num_type>
inline Box3D<num_type> Box3D<num_type>::operator/ (float_t scale) const
{
    auto box = *this;
    box /= scale;
    return box;
}

template <typename num_type>
inline void Box3D<num_type>::operator*= (float_t scale)
{
    m_corner[0] *= scale;
    m_corner[1] *= scale;
}

template <typename num_type>
inline void Box3D<num_type>::operator/= (float_t scale)
{
    m_corner[0] /= scale;
    m_corner[1] /= scale;   
}

template <typename num_type>
inline bool Box3D<num_type>::operator< (const Box3D<num_type> & b) const
{
    return math::GT(m_corner[0][0], b[0][0]) && math::GT(m_corner[0][1], b[0][1]) && math::GT(m_corner[0][2], b[0][2]) &&
           math::LT(m_corner[1][0], b[1][0]) && math::LT(m_corner[1][1], b[1][1]) && math::LT(m_corner[1][2], b[1][2]);
}

template <typename num_type>
inline bool Box3D<num_type>::operator> (const Box3D<num_type> & b) const
{
    return math::LT(m_corner[0][0], b[0][0]) && math::LT(m_corner[0][1], b[0][1]) && math::LT(m_corner[0][2], b[0][2]) &&
           math::GT(m_corner[1][0], b[1][0]) && math::GT(m_corner[1][1], b[1][1]) && math::GT(m_corner[1][2], b[1][2]);
}

template <typename num_type>
inline bool Box3D<num_type>::operator> (const Point3D<num_type> & p) const
{
    return math::LT(m_corner[0][0], p[0]) && math::LT(m_corner[0][1], p[1]) && math::LT(m_corner[0][2], p[2]) &&
           math::GT(m_corner[1][0], p[0]) && math::GT(m_corner[1][1], p[1]) && math::GT(m_corner[1][2], p[2]);
}

template <typename num_type>
inline bool Box3D<num_type>::operator<= (const Box3D<num_type> & b) const
{
    return math::GE(m_corner[0][0], b[0][0]) && math::GE(m_corner[0][1], b[0][1]) && math::GE(m_corner[0][2], b[0][2]) &&
           math::LE(m_corner[1][0], b[1][0]) && math::LE(m_corner[1][1], b[1][1]) && math::LE(m_corner[1][2], b[1][2]);
}

template <typename num_type>
inline bool Box3D<num_type>::operator>= (const Box3D<num_type> & b) const
{
    return math::LE(m_corner[0][0], b[0][0]) && math::LE(m_corner[0][1], b[0][1]) && math::LE(m_corner[0][2], b[0][2]) &&
           math::GE(m_corner[1][0], b[1][0]) && math::GE(m_corner[1][1], b[1][1]) && math::GE(m_corner[1][2], b[1][2]);
}

template <typename num_type>
inline bool Box3D<num_type>::operator>= (const Point3D<num_type> & p) const
{
    return math::LE(m_corner[0][0], p[0]) && math::LE(m_corner[0][1], p[1]) && math::LE(m_corner[0][2], p[2]) &&
           math::GE(m_corner[1][0], p[0]) && math::GE(m_corner[1][1], p[1]) && math::GE(m_corner[1][2], p[2]);
}

template <typename num_type>
inline Box3D<num_type> Box3D<num_type>::operator+ (const Box3D<num_type> & b) const
{
    Box3D<num_type> box;
    for(size_t i = 0; i < 3; ++i){
        box[0][i] = std::min(m_corner[0][i], b[0][i]);
        box[1][i] = std::max(m_corner[1][i], b[1][i]);
    }
    return box;
}

template <typename num_type>
inline void Box3D<num_type>::operator+= (const Box3D<num_type> & b)
{
    (*this) = (*this) + b;
}

template <typename num_type>
inline void Box3D<num_type>::operator|= (const Box3D<num_type> & b)
{
    *this += b;
}

template <typename num_type>
inline void Box3D<num_type>::operator|= (const Point3D<num_type> & p)
{
    for(size_t i = 0; i < 3; ++i){
        m_corner[0][i] = std::min(m_corner[0][i], p[i]);
        m_corner[1][i] = std::max(m_corner[1][i], p[i]);
    }
}

template <typename num_type>
inline Box3D<num_type> Box3D<num_type>::operator& (const Box3D<num_type> & b) const
{
    Box3D<num_type> box;
    for(size_t i = 0; i < 3; ++i){
        box[0][i] = std::max(m_corner[0][i], b[0][i]);
        box[1][i] = std::min(m_corner[1][i], b[1][i]);
    }
    return box;
}

template <typename num_type>
inline void Box3D<num_type>::operator&= (const Box3D<num_type> & b)
{
    (*this) = (*this) & b;
}

template <typename num_type>
inline void Box3D<num_type>::operator+= (const Vector3D<num_type> & v)
{
    m_corner[0] += v;
    m_corner[1] += v;
}

template <typename num_type>
inline void Box3D<num_type>::operator-= (const Vector3D<num_type> & v)
{
    m_corner[0] -= v;
    m_corner[1] -= v;
}

template <typename num_type>
inline Box3D<num_type> Box3D<num_type>::operator+ (const Vector3D<num_type> & v) const
{
    auto box = *this;
    box += v;
    return box;
}

template <typename num_type>
inline Box3D<num_type> Box3D<num_type>::operator- (const Vector3D<num_type> & v) const
{
    auto box = *this;
    box -= v;
    return box;
}

template <typename num_type>
inline Point3D<float_type<num_type> > Box3D<num_type>::Center() const
{
    float_t x = 0.5 * m_corner[0][0] + 0.5 * m_corner[1][0];
    float_t y = 0.5 * m_corner[0][1] + 0.5 * m_corner[1][1];
    float_t z = 0.5 * m_corner[0][2] + 0.5 * m_corner[1][2];
    return Point3D<float_t>(x, y, z);
}

template <typename num_type>
inline num_type Box3D<num_type>::Length() const { return m_corner[1][0] - m_corner[0][0]; }

template <typename num_type>
inline num_type Box3D<num_type>::Width() const { return m_corner[1][1] - m_corner[0][1]; }

template <typename num_type>
inline num_type Box3D<num_type>::Height() const { return m_corner[1][2] - m_corner[0][2]; }

template <typename num_type>
inline Point3D<num_type> Box3D<num_type>::Diagonal() const
{
    return m_corner[1] - m_corner[0];
}

template <typename num_type>
inline num_type Box3D<num_type>::HalfArea() const
{
    Point3D<num_type> dia = Diagonal();
    return (dia[0] + dia[1]) * dia[2] + dia[0] * dia[1];
}

template <typename num_type>
inline num_type Box3D<num_type>::SurfArea() const
{
    return 2 * HalfArea();
}

template <typename num_type>
inline num_type Box3D<num_type>::Volume() const
{
    Point3D<num_type> dia = Diagonal();
    return dia[0] * dia[1] * dia[2];
}

template <typename num_type>
inline size_t Box3D<num_type>::LargestAxis() const
{
    auto dia = Diagonal();
    size_t axis = 0;
    if(dia[0] < dia[1]) axis = 1;
    if(dia[axis] < dia[2]) axis = 2;
    return axis;
}

template <typename num_type>
inline bool Box3D<num_type>::isValid() const
{
    return math::LE(m_corner[0][0], m_corner[1][0]) && math::LE(m_corner[0][1], m_corner[1][1]) && math::LE(m_corner[0][2], m_corner[1][2]);
}

template <typename num_type>
inline void Box3D<num_type>::Normalize()
{
    for(int i = 0; i < 3; ++i){
        num_type min = std::min(m_corner[0][i], m_corner[1][i]);
        num_type max = std::max(m_corner[0][i], m_corner[1][i]);
        m_corner[0][i] = min; m_corner[1][i] = max;
    }
}

template <typename num_type>
inline void Box3D<num_type>::SetInvalid()
{
    *this = Box3D();
}

template <typename num_type>
inline void Box3D<num_type>::Scale(float_t scale)
{
    auto ct = Center().template Cast<num_type>();
    *this -= ct;
    *this *= scale;
    *this += ct;
}

template <typename num_type>
template<typename other_num_type>
inline Box3D<other_num_type> Box3D<num_type>::Cast() const
{
    Box3D<other_num_type> box;
    box[0] = m_corner[0].template Cast<other_num_type>();
    box[1] = m_corner[1].template Cast<other_num_type>();
    return box;
}

template <typename num_type>
inline bool Box3D<num_type>::Collision(const Box3D<num_type> & a, const Box3D<num_type> & b, bool considerTouch)
{
    if(considerTouch)
        return Collision(a, b, std::true_type{});
    else return Collision(a, b, std::false_type{});
}

template <typename num_type>
inline bool Box3D<num_type>::Collision(const Box3D<num_type> & a, const Box3D<num_type> & b, std::true_type)
{
    if(math::GE(a[1][0], b[0][0]) && math::LE(a[0][0], b[1][0]) &&
       math::GE(a[1][1], b[0][1]) && math::LE(a[0][1], b[1][1]) &&
       math::GE(a[1][2], b[0][2]) && math::LE(a[0][2], b[1][2]))
       return true;
    return false;
}

template <typename num_type>
inline bool Box3D<num_type>::Collision(const Box3D<num_type> & a, const Box3D<num_type> & b, std::false_type)
{
    if(math::GT(a[1][0], b[0][0]) && math::LT(a[0][0], b[1][0]) &&
       math::GT(a[1][1], b[0][1]) && math::LT(a[0][1], b[1][1]) &&
       math::GT(a[1][2], b[0][2]) && math::LT(a[0][2], b[1][2]))
       return true;
    return false;
}

}//namespace geometry
}//namespace generic