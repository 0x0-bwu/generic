/**
 * @file Segment.hpp
 * @author bwu
 * @brief Model of segment concept
 * @version 0.1
 * @date 2022-02-22
 */
#pragma once
#include "generic/common/Traits.hpp"
#include "Point.hpp"
#include <cmath>
namespace generic  {
namespace geometry {
using generic::common::float_type;
template <typename num_type>
class Segment2D
{
public:
    const static size_t dim = 2;
    using coor_t = num_type;
    using point_t = Point2D<coor_t>;

    ///@brief constructs a zero length segment from point(0, 0)
    Segment2D(){}
    ///@brief constructs a segment from point `s` to `e`
    Segment2D(const Point2D<num_type> & s, const Point2D<num_type> & e);

    ///@brief checks if this segment equals to `seg`
    bool operator== (const Segment2D & seg) const;
    ///@brief checks if this segment not equals to `seg`
    bool operator!= (const Segment2D & seg) const;

    ///@brief accesses segment point by index 0-1
    Point2D<num_type> & operator[](size_t i);
    const Point2D<num_type> & operator[](size_t i) const;

    ///@brief converts to segment with other number type explicitly
    template<typename other_num_type>
    Segment2D<other_num_type> Cast() const;

    ///@brief gets shortest distance of two segments
    static float_type<num_type> Distance(const Segment2D<num_type> & p, const Segment2D<num_type> & q);

private:
    Point2D<num_type> m_points[2];
};

template <typename num_type>
class Segment3D
{
public:
    const static size_t dim = 3;
    using coor_t = num_type;
    using point_t = Point3D<coor_t>;

    ///@brief constructs a zero length segment from point(0, 0, 0)
    Segment3D(){}
    ///@brief constructs a segment from point `s` to `e`
    Segment3D(const Point3D<num_type> & s, const Point3D<num_type> & e);

    ///@brief checks if this segment equals to `seg`
    bool operator== (const Segment3D & seg) const;
    ///@brief checks if this segment not equals to `seg`
    bool operator!= (const Segment3D & seg) const;

    ///@brief accesses segment point by index 0-1
    Point3D<num_type> & operator[](int i);
    const Point3D<num_type> & operator[](int i) const;

    ///@brief converts to segment with other number type explicitly
    template<typename other_num_type>
    Segment3D<other_num_type> Cast() const;

    ///@brief gets shortest distance of two segments
    static float_type<num_type> Distance(const Segment3D<num_type> & p, const Segment3D<num_type> & q);

private:
    Point3D<num_type> m_points[2];
};

template <typename num_type>
inline Segment2D<num_type>::Segment2D(const Point2D<num_type> & s, const Point2D<num_type> & e)
{
    m_points[0] = s; m_points[1] = e;
}

template <typename num_type>
inline bool Segment2D<num_type>::operator==(const Segment2D<num_type> & seg) const
{
    return m_points[0] == seg.m_points[0] && m_points[1] == seg.m_points[1];
}

template <typename num_type>
inline bool Segment2D<num_type>::operator!=(const Segment2D<num_type> & seg) const
{
    return !(*this == seg);
}

template <typename num_type>
inline Point2D<num_type> & Segment2D<num_type>::operator[](size_t i) { return m_points[i]; }

template <typename num_type>
inline const Point2D<num_type> & Segment2D<num_type>::operator[](size_t i) const { return m_points[i]; }

template <typename num_type>
template<typename other_num_type>
inline Segment2D<other_num_type>  Segment2D<num_type>::Cast() const
{
    return Segment2D<other_num_type>(m_points[0].template Cast<other_num_type>(),
                                     m_points[1].template Cast<other_num_type>());
}

template <typename num_type>
float_type<num_type> Segment2D<num_type>::Distance(const Segment2D<num_type> & p, const Segment2D<num_type> & q)
{
    Segment3D<num_type> u(Point3D<num_type>(p[0][0], p[0][1], 0), Point3D<num_type>(p[1][0], p[1][1], 0));
    Segment3D<num_type> v(Point3D<num_type>(q[0][0], q[0][1], 0), Point3D<num_type>(q[1][0], q[1][1], 0));
    return Segment3D<num_type>::Distance(u, v);
}

template <typename num_type>
inline Segment3D<num_type>::Segment3D(const Point3D<num_type> & s, const Point3D<num_type> & e)
{
    m_points[0] = s; m_points[1] = e;
}

template <typename num_type>
inline bool Segment3D<num_type>::operator==(const Segment3D<num_type> & seg) const
{
    return m_points[0] == seg.m_points[0] && m_points[1] == seg.m_points[1];
}

template <typename num_type>
inline bool Segment3D<num_type>::operator!=(const Segment3D<num_type> & seg) const
{
    return !(*this == seg);
}

template <typename num_type>
inline Point3D<num_type> & Segment3D<num_type>::operator[](int i) { return m_points[i]; }

template <typename num_type>
inline const Point3D<num_type> & Segment3D<num_type>::operator[](int i) const { return m_points[i]; }

template <typename num_type>
template<typename other_num_type>
inline Segment3D<other_num_type>  Segment3D<num_type>::Cast() const
{
    return Segment3D<other_num_type>(m_points[0].template Cast<other_num_type>(),
                                     m_points[1].template Cast<other_num_type>());
}

template <typename num_type>
float_type<num_type> Segment3D<num_type>::Distance(const Segment3D<num_type> & p, const Segment3D<num_type> & q)
{
    //Implement based on [Robust Computation of Distance Between Line Segments] by David Eberly
    using float_type = common::float_type<num_type>;
    float_type a = (p[1] - p[0]).Dot(p[1] - p[0]);
    float_type b = (p[1] - p[0]).Dot(q[1] - q[0]);
    float_type c = (q[1] - q[0]).Dot(q[1] - q[0]);
    float_type d = (p[1] - p[0]).Dot(p[0] - q[0]);
    float_type e = (q[1] - q[0]).Dot(p[0] - q[0]);
    float_type det = a * c - b * b;
    float_type bte, ctd, ate, btd, s, t;
    if(det > 0) {
        bte = b * e; ctd = c * d;
        if(bte <= ctd){
            if( e <= 0) { s = (-d >= a ? 1 : (-d > 0 ? -d / a : 0)); t = 0; }
            else if( e < c) { s = 0; t = e / c; }
            else { s = (b - d >= a ? 1 : (b - d > 0 ? (b - d) / a : 0)); t = 1; }
        }
        else{
            s = bte - ctd;
            if(s >= det){
                if(b + e <= 0) { s = (-d <= 0 ? 0 : (-d < a ? -d / a : 1)); t = 0; }
                else if(b + e < c) { s = 1; t = (b + e) / c; }
                else { s = (b - d <= 0 ? 0 : (b - d < a ? (b - d) / a : 1)); t = 1; }
            }
            else{
                ate = a * e; btd = b * d;
                if(ate <= btd) { s = (-d <= 0 ? 0 : (-d >= a ? 1 : -d / a)); t = 0; }
                else{
                    t = ate - btd;
                    if(t >= det) { s = (b - d <= 0 ? 0 : (b - d >= a ? 1 : (b - d) / a)); t = 1; }
                    else { s /= det; t /= det; }
                }
            }
        }
    }
    else{
        if(e <= 0) { s = (-d <= 0 ? 0 : (-d >= a ? 1 : -d / a)); t = 0; }
        else if(e >= c) { s = (b - d <= 0 ? 0 : (b - d >= a ? 1 : (b - d) / a)); t = 1; }
        else{ s = 0; t = e / c; }
    }
    
    Point3D<float_type> fp0(p[0][0], p[0][1], p[0][2]);
    Point3D<float_type> fp1(p[1][0], p[1][1], p[1][2]);
    Point3D<float_type> fq0(q[0][0], q[0][1], q[0][2]);
    Point3D<float_type> fq1(q[1][0], q[1][1], q[1][2]);
    Point3D<float_type> temp = (fp0 * (1 - s) + fp1 * s) - (fq0 * (1 - t) + fq1 * t);
    return std::sqrt(temp.NormSquare());
}
}//namespace geometry
}//namespace generic