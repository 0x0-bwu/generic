#ifndef GENERIC_GEOMETRY_TRANSFORM_HPP
#define GENERIC_GEOMETRY_TRANSFORM_HPP
#include <boost/qvm/quat_vec_operations.hpp>
#include <boost/qvm/vec_mat_operations.hpp>
#include <boost/qvm/quat_operations.hpp>
#include <boost/qvm/mat_operations.hpp>
#include <boost/qvm/vec_operations.hpp>
#include <boost/qvm/vec_traits.hpp>
#include <boost/qvm/quat.hpp>
#include <boost/qvm/mat.hpp>
#include <boost/qvm/vec.hpp>
#include "generic/math/Numbers.hpp"
#include "Geometries.hpp"
#include "Vector.hpp"
namespace boost {
namespace qvm {

using namespace generic::geometry;
template <typename num_type>
struct vec_traits<Vector3D<num_type> >
{
    static const int dim = 3;
    using scalar_type = num_type;

    template <int I>
    static inline scalar_type & write_element(Vector3D<num_type> & vec) { return vec[I]; }

    template <int I>
    static inline scalar_type read_element(const Vector3D<num_type> & vec) { return vec[I]; }

    static inline scalar_type & write_element_idx(int i, Vector3D<num_type> & vec) { return vec[i]; }

    static inline scalar_type read_element_idx(int i, const Vector3D<num_type> & vec) { return vec[i]; }
};

template <typename num_type, size_t N>
struct vec_traits<VectorN<num_type, N> >
{
    static const int dim = static_cast<int>(N);
    using scalar_type = num_type;

    template <int I>
    static inline scalar_type & write_element(VectorN<num_type, N> & vec) { return vec[I]; }

    template <int I>
    static inline scalar_type read_element(const VectorN<num_type, N> & vec) { return vec[I]; }

    static inline scalar_type & write_element_idx(int i, VectorN<num_type, N> & vec) { return vec[i]; }

    static inline scalar_type read_element_idx(int i, const VectorN<num_type, N> & vec) { return vec[i]; }
};
}//namespace qvm
}//namespace boost

namespace generic {
namespace geometry{

template <typename num_type>
using Matrix2x2 = boost::qvm::mat<num_type, 2, 2>;

template <typename num_type>
using Matrix3x3 = boost::qvm::mat<num_type, 3, 3>;

template <typename num_type>
using Matrix4x4 = boost::qvm::mat<num_type, 4, 4>;

template <typename float_t> class Transform2D;
template <typename float_t> class Transform3D;
template <typename float_t> class Quaternion;

template <typename float_t, typename num_type>
inline Transform2D<float_t> makeShiftTransform2D(const Vector2D<num_type> & shift)
{
    Transform2D<float_t> trans;
    trans(0, 2) = float_t(shift[0]);
    trans(1, 2) = float_t(shift[1]);
    return trans;
}

template <typename float_t>
inline Transform2D<float_t> makeRotateTransform2D(float_t rot)
{
    float_t c = std::cos(rot);
    float_t s = std::sin(rot);
    Transform2D<float_t> trans;
    trans(0, 0) = c; trans(0, 1) = -s;
    trans(1, 0) = s; trans(1, 1) =  c;
    return trans;
}

template <typename float_t>
inline Transform2D<float_t> makeScaleTransform2D(float_t scale)
{
    Transform2D<float_t> trans;
    trans(0, 0) = scale; trans(1,1) = scale;
    return trans;
}

template <typename float_t>
inline Transform2D<float_t> makeMirroredTransform2D(Axis axis)
{
    Transform2D<float_t> trans;
    if(axis == Axis::X) trans(1, 1) *= -1;
    else if(axis == Axis::Y) trans(0, 0) *= -1;
    return trans;
}

template <typename float_t, typename num_type>
inline Transform3D<float_t> makeShiftTransform3D(const Vector3D<num_type> & shift)
{
    Transform3D<float_t> trans;
    trans(0, 3) = float_t(shift[0]);
    trans(1, 3) = float_t(shift[1]);
    trans(2, 3) = float_t(shift[2]);
    return trans;
}

template <typename float_t>
inline Transform3D<float_t> makeScaleTransform3D(float_t scale)
{
    Transform3D<float_t> trans;
    trans(0, 0) = scale; trans(1,1) = scale; trans(2, 2) = scale;
    return trans;
}

template <typename float_t>
class Transform2D
{
    static_assert(std::is_floating_point<float_t>::value, "only floating point type support in transform construct!");
public:
    static const size_t dim = 2;
    Transform2D();
    explicit Transform2D(const Matrix3x3<float_t> & m);

    float_t & operator() (size_t row, size_t col);
    const float_t & operator() (size_t row, size_t col) const;

    bool Inverse();

    Transform2D<float_t> operator * (const Transform2D<float_t> & trans) const;
    Transform2D<float_t> & Prod(const Transform2D<float_t> & tran);
    
    template <typename num_type>
    Point2D<num_type> operator * (const Point2D<num_type> & point) const;

    template <typename num_type>
    Segment2D<num_type> operator * (const Segment2D<num_type> & segment) const;

    template <typename num_type>
    Triangle2D<num_type> operator * (const Triangle2D<num_type> & triangle) const;

    template <typename num_type>
    Polygon2D<num_type> operator * (const Box2D<num_type> & box) const;

    template <typename num_type>
    Polyline2D<num_type> operator * (const Polyline2D<num_type> & polyline) const;

    template <typename num_type>
    Polygon2D<num_type> operator * (const Polygon2D<num_type> & polygon) const;

    template <typename num_type>
    PolygonWithHoles2D<num_type> operator * (const PolygonWithHoles2D<num_type> & pwh) const;

    Matrix3x3<float_t> & GetMatrix() { return m_matrix; }
    const Matrix3x3<float_t> & GetMatrix() const { return m_matrix; }

    void GetCoeffs(float_t coeffs[], bool rowMajor = true);

private:
    Matrix3x3<float_t> m_matrix;
};

template <typename float_t>
class Transform3D
{
    static_assert(std::is_floating_point<float_t>::value, "only floating point type support in transform construct!");
public:
    static const size_t dim = 3;
    Transform3D();
    explicit Transform3D(const Transform2D<float_t> & trans2d);
    explicit Transform3D(const Matrix4x4<float_t> & m);

    float_t & operator() (size_t row, size_t col);
    const float_t & operator() (size_t row, size_t col) const;

    bool Inverse();
    Transform2D<float_t> GetTransfrom2D() const;

    Transform3D<float_t> operator * (const Transform3D<float_t> & trans) const;
    Transform3D<float_t> & Prod(const Transform3D<float_t> & tran);

    template <typename num_type>
    Point3D<num_type> operator * (const Point3D<num_type> & point) const;

    template <typename num_type>
    Segment3D<num_type> operator * (const Segment3D<num_type> & segment) const;

    template <typename num_type>
    Triangle3D<num_type> operator * (const Triangle3D<num_type> & triangle) const;

    Matrix4x4<float_t> & GetMatrix() { return m_matrix; }
    const Matrix4x4<float_t> & GetMatrix() const { return m_matrix; }

    void GetCoeffs(float_t coeffs[], bool rowMajor = true);

private:
    Matrix4x4<float_t> m_matrix;
};

template <typename float_t>
class Quaternion
{
    static_assert(std::is_floating_point<float_t>::value, "only floating point type support in transform construct!");
public:
    static const size_t dim = 3;
    Quaternion();
    Quaternion(float_t q0, float_t q1, float_t q2, float_t q3);
    Quaternion(const Vector3D<float_t> & axis, float_t angle);
    Quaternion(const boost::qvm::quat<float_t> & q);

    float_t & operator[](size_t dim);
    const float_t & operator[](size_t dim) const;
    Quaternion<float_t> & operator *= (const Quaternion<float_t> & q);
    Quaternion<float_t> operator * (const Quaternion<float_t> & q) const;

    void SetAxisAndAngle(const Vector3D<float_t> & axis, float angle);
    Vector3D<float_t> Axis() const;
    float_t Angle() const;

    void Invert();
    void Negate();
    void Normalize();           
    float_t Mag() const;
    float_t MagSqrt() const;
    Quaternion<float_t> Log() const;
    Quaternion<float_t> Exp() const;
    Quaternion<float_t> Inverse() const;
    Vector3D<float_t> Rotate(const Vector3D<float_t> & vec) const;
    Vector3D<float_t> InverseRotate(const Vector3D<float_t> & vec) const;

    float_t Dot(const Quaternion<float_t> & q) const;

    boost::qvm::quat<float_t> & q() { return m_q; }
    const boost::qvm::quat<float_t> & q() const { return m_q; }

    static Quaternion<float_t> LnDif(const Quaternion<float_t> &a, const Quaternion<float_t> & b);
    static Quaternion<float_t> SquadTangent(const Quaternion<float_t> & before,
                                            const Quaternion<float_t> & center,
                                            const Quaternion<float_t> & after);
    static Quaternion<float_t> Squad(const Quaternion<float_t> & a, const Quaternion<float_t> & tgA,
                                     const Quaternion<float_t> & b, const Quaternion<float_t> & tgB, float_t t);

    static Quaternion<float_t> Slerp(const Quaternion<float_t> & a, const Quaternion<float_t> & b, float_t t);

private:
    boost::qvm::quat<float_t> m_q;
};

template <typename float_t>
inline Transform2D<float_t>::Transform2D()
{
    boost::qvm::set_identity(m_matrix);
}

template <typename float_t>
inline Transform2D<float_t>::Transform2D(const Matrix3x3<float_t> & m)
{
    boost::qvm::assign(m_matrix, m);
}

template <typename float_t>
inline float_t & Transform2D<float_t>::operator() (size_t row, size_t col)
{
    return m_matrix.a[row][col];
}

template <typename float_t>
inline const float_t & Transform2D<float_t>::operator() (size_t row, size_t col) const
{
    return m_matrix.a[row][col];
}

template <typename float_t>
inline bool Transform2D<float_t>::Inverse()
{
    float_t det = boost::qvm::determinant(m_matrix);
    if(math::EQ(det, float_t(0))) return false;
    m_matrix = boost::qvm::inverse(m_matrix, det);
    return true;
}

template <typename float_t>
inline Transform2D<float_t> Transform2D<float_t>::operator * (const Transform2D<float_t> & trans) const
{
    using namespace boost::qvm;
    return Transform2D<float_t>(m_matrix * trans.GetMatrix());
}

template <typename float_t>
inline Transform2D<float_t> & Transform2D<float_t>::Prod(const Transform2D<float_t> & trans)
{
    m_matrix = trans.GetMatrix() * m_matrix;
    return *this;
}

template <typename float_t>
template <typename num_type>
inline Point2D<num_type> Transform2D<float_t>::operator * (const Point2D<num_type> & point) const
{
    Point2D<num_type> res(point);
    Transform(res, *this);
    return res;
}

template <typename float_t>
template <typename num_type>
inline Segment2D<num_type> Transform2D<float_t>::operator * (const Segment2D<num_type> & segment) const
{
    Segment2D<num_type> res(segment);
    Transform(res, *this);
    return res;
}

template <typename float_t>
template <typename num_type>
inline Triangle2D<num_type> Transform2D<float_t>::operator * (const Triangle2D<num_type> & triangle) const
{
    Triangle2D<num_type> res(triangle);
    Transform(res, *this);
    return res;
}

template <typename float_t>
template <typename num_type>
inline Polygon2D<num_type> Transform2D<float_t>::operator * (const Box2D<num_type> & box) const
{
    Point2D<num_type> p[4] = { box[0], {box[1][0], box[0][1]}, box[1], {box[0][0], box[1][1]} };
    Polygon2D<num_type> res;
    for(size_t i = 0; i < 4; ++i){
        Transform(p[i], *this);
        res << p[i];
    }
    return res;
}

template <typename float_t>
template <typename num_type>
inline Polyline2D<num_type> Transform2D<float_t>::operator * (const Polyline2D<num_type> & polyline) const
{
    Polyline2D<num_type> res = polyline;
    Transform(res, *this);
    return res;  
}

template <typename float_t>
template <typename num_type>
inline Polygon2D<num_type>  Transform2D<float_t>::operator * (const Polygon2D<num_type> & polygon) const
{
    Polygon2D<num_type> res(polygon);
    Transform(res, *this);
    return res;
}

template <typename float_t>
template <typename num_type>
inline PolygonWithHoles2D<num_type> Transform2D<float_t>::operator * (const PolygonWithHoles2D<num_type> & pwh) const
{
    PolygonWithHoles2D<num_type> res(pwh);
    Transform(res, *this);
    return res;
}

template <typename float_t>
inline void Transform2D<float_t>::GetCoeffs(float_t coeffs[], bool rowMajor)
{
    if(rowMajor){
        for(size_t i = 0; i <= dim; ++i)
            for(size_t j = 0; j <= dim; ++j)
                coeffs[(dim + 1) * i + j] = m_matrix.a[i][j];
    }
    else{
        for(size_t i = 0; i <= dim; ++i)
            for(size_t j = 0; j <= dim; ++j)
                coeffs[(dim + 1) * i + j] = m_matrix.a[j][i];
    }
}


template <typename float_t>
inline Transform3D<float_t>::Transform3D()
{
    boost::qvm::set_identity(m_matrix);
}

template <typename float_t>
inline Transform3D<float_t>::Transform3D(const Transform2D<float_t> & trans2d)
{
    boost::qvm::set_identity(m_matrix);
    for(size_t i = 0; i < 2; ++i)
        for(size_t j = 0; j < 2; ++j)
            m_matrix.a[i][j] = trans2d(i, j);
    m_matrix.a[0][3] = trans2d(0, 2);
    m_matrix.a[1][3] = trans2d(1, 2);
}

template <typename float_t>
inline Transform3D<float_t>::Transform3D(const Matrix4x4<float_t> & m)
{
    boost::qvm::assign(m_matrix, m);
}

template <typename float_t>
inline float_t & Transform3D<float_t>::operator() (size_t row, size_t col)
{
    return m_matrix.a[row][col];
}

template <typename float_t>
inline const float_t & Transform3D<float_t>::operator() (size_t row, size_t col) const
{
    return m_matrix.a[row][col];
}

template <typename float_t>
inline bool Transform3D<float_t>::Inverse()
{
    float_t det = boost::qvm::determinant(m_matrix);
    if(math::EQ(det, float_t(0))) return false;
    m_matrix = boost::qvm::inverse(m_matrix, det);
    return true;
}

template <typename float_t>
inline Transform2D<float_t> Transform3D<float_t>::GetTransfrom2D() const
{
    Transform2D<float_t> trans;
    for(size_t i = 0; i < 2; ++i)
        for(size_t j = 0; j < 2; ++j)
            trans(i, j) = m_matrix.a[i][j];
    trans(0, 2) = m_matrix.a[0][3];
    trans(1, 2) = m_matrix.a[1][3];
    return trans;
}

template <typename float_t>
inline Transform3D<float_t> Transform3D<float_t>::operator * (const Transform3D<float_t> & trans) const
{
    using namespace boost::qvm;
    return Transform3D<float_t>(m_matrix * trans.GetMatrix());
}

template <typename float_t>
inline Transform3D<float_t> & Transform3D<float_t>::Prod(const Transform3D<float_t> & trans)
{
    m_matrix = trans.GetMatrix() * m_matrix;
    return *this;
}

template <typename float_t>
template <typename num_type>
inline Point3D<num_type> Transform3D<float_t>::operator * (const Point3D<num_type> & point) const
{
    Point3D<num_type> res(point);
    Transform(res, *this);
    return res;
}

template <typename float_t>
template <typename num_type>
inline Segment3D<num_type> Transform3D<float_t>::operator * (const Segment3D<num_type> & segment) const
{
    Segment3D<num_type> res(segment);
    Transform(res, *this);
    return res;
}

template <typename float_t>
template <typename num_type>
inline Triangle3D<num_type> Transform3D<float_t>::operator * (const Triangle3D<num_type> & triangle) const
{
    Triangle3D<num_type> res(triangle);
    Transform(res, *this);
    return res;
}

template <typename float_t>
inline void Transform3D<float_t>::GetCoeffs(float_t coeffs[], bool rowMajor)
{
    if(rowMajor){
        for(size_t i = 0; i <= dim; ++i)
            for(size_t j = 0; j <= dim; ++j)
                coeffs[(dim + 1) * i + j] = m_matrix.a[i][j];
    }
    else{
        for(size_t i = 0; i <= dim; ++i)
            for(size_t j = 0; j <= dim; ++j)
                coeffs[(dim + 1) * i + j] = m_matrix.a[j][i];    }
}

template <typename float_t>
inline Quaternion<float_t>::Quaternion()
{
    m_q.a[0] = 1;
    m_q.a[1] = m_q.a[2] = m_q.a[3] = 0;
}

template <typename float_t>
inline Quaternion<float_t>::Quaternion(float_t q0, float_t q1, float_t q2, float_t q3)
{
    m_q.a[0] = q0;
    m_q.a[1] = q1; m_q.a[2] = q2; m_q.a[3] = q3;
}

template <typename float_t>
inline Quaternion<float_t>::Quaternion(const Vector3D<float_t> & axis, float_t angle)
{
    SetAxisAndAngle(axis, angle);
}

template <typename float_t>
inline Quaternion<float_t>::Quaternion(const boost::qvm::quat<float_t> & q)
{
    boost::qvm::assign(m_q, q);
}

template <typename float_t>
inline float_t & Quaternion<float_t>::operator[](size_t dim)
{
    return m_q.a[dim];
}

template <typename float_t>
inline const float_t & Quaternion<float_t>::operator[](size_t dim) const
{
    return m_q.a[dim];
}

template <typename float_t>
inline Quaternion<float_t> & Quaternion<float_t>::operator *= (const Quaternion<float_t> & q)
{
    *this = (*this) * q;
    return *this;
}

template <typename float_t>
inline Quaternion<float_t> Quaternion<float_t>::operator * (const Quaternion<float_t> & q) const
{
    return Quaternion<float_t>(m_q * q.q());
}

template <typename float_t>
inline void Quaternion<float_t>::SetAxisAndAngle(const Vector3D<float_t> & axis, float angle)
{
    boost::qvm::set_rot(m_q, axis, angle);
}

template <typename float_t>
inline Vector3D<float_t> Quaternion<float_t>::Axis() const
{
    Vector3D<float_t> res(m_q.a[1], m_q.a[2], m_q.a[3]);
    float_t sinus = res.Norm2();
    if(sinus > std::numeric_limits<float_t>::epsilon()) res /= sinus;
    return math::LE(std::acos(m_q.a[0]), float_t(math::pi_half)) ? res : -res;
}

template <typename float_t>
float_t Quaternion<float_t>::Angle() const
{
    float_t angle = 2.0 * std::acos(m_q.a[0]);
    return math::LE(angle, float_t(math::pi)) ? angle : math::pi_2 - angle;
}

template <typename float_t>
inline void Quaternion<float_t>::Invert()
{
    boost::qvm::assign(m_q, boost::qvm::inverse(m_q));
}

template <typename float_t>
inline void Quaternion<float_t>::Negate()
{
    Invert();
    m_q.a[0] = -m_q.a[0];
}

template <typename float_t>
inline void Quaternion<float_t>::Normalize()
{
    boost::qvm::normalize(m_q);
}

template <typename float_t>
inline float_t Quaternion<float_t>::Mag() const
{
    return boost::qvm::mag(m_q);
}

template <typename float_t>
inline float_t Quaternion<float_t>::MagSqrt() const
{
    return boost::qvm::mag_sqr(m_q);
}

template <typename float_t>
inline Quaternion<float_t> Quaternion<float_t>::Log() const
{
    float_t len = std::sqrt(m_q.a[1] * m_q.a[1] + m_q.a[2] * m_q.a[2] + m_q.a[3] * m_q.a[3]);
    if(len < std::numeric_limits<float_t>::epsilon())
        return Quaternion<float_t>(0, m_q.a[1], m_q.a[2], m_q.a[3]);
    else{
        float_t coef = std::acos(m_q.a[0]) / len;
        return Quaternion<float_t>(0, coef * m_q.a[1], coef * m_q.a[2], coef * m_q.a[3]);
    }
}

template <typename float_t>
inline Quaternion<float_t> Quaternion<float_t>::Exp() const
{
    float_t theta = std::sqrt(m_q.a[1] * m_q.a[1] + m_q.a[2] * m_q.a[2] + m_q.a[3] * m_q.a[3]);
    if(theta < std::numeric_limits<float_t>::epsilon())
        return Quaternion<float_t>(std::cos(theta), m_q.a[1], m_q.a[2], m_q.a[3]);
    else{
        float_t coef = std::sin(theta) / theta;
        return Quaternion<float_t>(std::cos(theta), coef * m_q.a[1], coef * m_q.a[2], coef * m_q.a[3]);
    }
}

template <typename float_t>
inline Quaternion<float_t> Quaternion<float_t>::Inverse() const
{
    return Quaternion<float_t>(boost::qvm::inverse(m_q));
}

template <typename float_t>
inline Vector3D<float_t> Quaternion<float_t>::Rotate(const Vector3D<float_t> & vec) const
{
   using namespace boost::qvm;
   return m_q * vec;
}

template <typename float_t>
inline Vector3D<float_t> Quaternion<float_t>::InverseRotate(const Vector3D<float_t> & vec) const
{
    return Inverse().Rotate(vec);
}

template <typename float_t>
inline float_t Quaternion<float_t>::Dot(const Quaternion<float_t> & q) const
{
    return boost::qvm::dot(m_q, q.q());
}

template <typename float_t>
inline Quaternion<float_t> Quaternion<float_t>::LnDif(const Quaternion<float_t> &a, const Quaternion<float_t> & b)
{
    auto dif = a.Inverse() * b;
    dif.Normalize();
    return dif.Log();
}

template <typename float_t>
inline Quaternion<float_t> Quaternion<float_t>::SquadTangent(const Quaternion<float_t> & before,
                                                                const Quaternion<float_t> & center,
                                                                    const Quaternion<float_t> & after)
{
    auto l1 = LnDif(center, before);
    auto l2 = LnDif(center, after);
    Quaternion<float_t> e;
    for(size_t i = 0; i < 4; ++i)
        e[i] = -0.25 * (l1[i] + l2[2]);
    e = center * e.Exp();
    return e;
}

template <typename float_t>
inline Quaternion<float_t> Quaternion<float_t>::Squad(const Quaternion<float_t> & a, const Quaternion<float_t> & tgA,
                                                        const Quaternion<float_t> & b, const Quaternion<float_t> & tgB, float_t t)
{
    auto ab = Slerp(a, b, t);
    auto tg = Slerp(tgA, tgB, t);
    return Slerp(ab, tg, 2.0 * t * (1.0 - t));
}

template <typename float_t>
inline Quaternion<float_t> Quaternion<float_t>::Slerp(const Quaternion<float_t> & a, const Quaternion<float_t> & b, float_t t)
{
    using namespace boost::qvm;
    return Quaternion<float_t>(slerp(a.q(), b.q(), t));
}

template <typename num_type, typename float_t>
inline void Transform(Point2D<num_type> & point, const Transform2D<float_t> & trans)
{
    VectorN<float_t, 3> vec = trans.GetMatrix() * VectorN<float_t, 3>(point[0], point[1], float_t(1));
    point = Point2D<num_type>(vec[0], vec[1]);
}

template <typename num_type, typename float_t>
inline void Transform(Segment2D<num_type> & segment, const Transform2D<float_t> & trans)
{
    Transform(segment[0], trans);
    Transform(segment[1], trans);
}

template <typename num_type, typename float_t>
inline void Transform(Triangle2D<num_type> & triangle, const Transform2D<float_t> & trans)
{
    Transform(triangle[0], trans);
    Transform(triangle[1], trans);
    Transform(triangle[2], trans);
}

template <typename num_type, typename float_t>
inline void Transform(Polyline2D<num_type> & polyline, const Transform2D<float_t> & trans)
{
    for(auto & point : polyline)
        Transform(point, trans);
}

template <typename num_type, typename float_t>
inline void Transform(Polygon2D<num_type> & polygon, const Transform2D<float_t> & trans)
{
    for(auto iter = polygon.Begin(); iter != polygon.End(); ++iter)
        Transform(*iter, trans);
}

template <typename num_type, typename float_t>
inline void Transform(PolygonWithHoles2D<num_type> & pwh, const Transform2D<float_t> & trans)
{
    Transform(pwh.outline, trans);
    auto iter = pwh.BeginHoles();
    for(; iter != pwh.EndHoles(); ++iter)
        Transform(*iter, trans);
}

template <typename num_type, typename float_t>
inline void Transform(Point3D<num_type> & point, const Transform3D<float_t> & trans)
{
    VectorN<float_t, 4> vec(trans.GetMatrix() * VectorN<float_t, 4>(point[0], point[1], point[2], float_t(1)));
    point = Point3D<num_type>(vec[0], vec[1], vec[2]);
}

template <typename num_type, typename float_t>
inline void Transform(Segment3D<num_type> & segment, const Transform3D<float_t> & trans)
{
    Transform(segment[0], trans);
    Transform(segment[1], trans);
}

template <typename num_type, typename float_t>
inline void Transform(Triangle3D<num_type> & triangle, const Transform3D<float_t> & trans)
{
    Transform(triangle[0], trans);
    Transform(triangle[1], trans);
    Transform(triangle[2], trans);
}

template <typename geometry_t, typename iterator, typename transform_t, 
          typename std::enable_if<geometry_t::dim == transform_t::dim && std::is_same<geometry_t,
          typename std::iterator_traits<iterator>::value_type>::value, bool>::type = true>
inline void Transform(iterator begin, iterator end, const transform_t & trans)
{
    for(auto iter = begin; iter != end; ++iter){
        geometry_t & geom = *iter;
        Transform(geom, trans);
    }
}
}//namespace geometry
}//namespace generic
#endif//GENERIC_GEOMETRY_TRANSFORM_HPP