/**
 * @file Transform.hpp
 * @author bwu
 * @brief Model of transform vector, matrix and quaternion concept
 * @version 0.1
 * @date 2022-02-22
 */
#pragma once
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
#include "generic/math/MathIO.hpp"
#include "generic/tools/Hash.hpp"
#include "Geometries.hpp"
#include "Vector.hpp"

namespace boost::qvm {

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
}//namespace boost::qvm

namespace generic::geometry {

template <typename num_type>
using Matrix2x2 = boost::qvm::mat<num_type, 2, 2>;

template <typename num_type>
using Matrix3x3 = boost::qvm::mat<num_type, 3, 3>;

template <typename num_type>
using Matrix4x4 = boost::qvm::mat<num_type, 4, 4>;

template <typename float_t> class Transform2D;
template <typename float_t> class Transform3D;
template <typename float_t> class Quaternion;

///@brief makes a 2d shift transform matrix by a vector2d 
template <typename float_t, typename num_type>
inline Transform2D<float_t> makeShiftTransform2D(const Vector2D<num_type> & shift)
{
    Transform2D<float_t> trans;
    trans(0, 2) = float_t(shift[0]);
    trans(1, 2) = float_t(shift[1]);
    return trans;
}

///@brief makes a 2d rotation transform matrix by radian, unit: rad
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

///@brief makes a 2d scale transform matrix by scale factor
template <typename float_t>
inline Transform2D<float_t> makeScaleTransform2D(float_t scale)
{
    Transform2D<float_t> trans;
    trans(0, 0) = scale; trans(1,1) = scale;
    return trans;
}

///@brief makes a 2d mirror transform matrix by given axis
template <typename float_t>
inline Transform2D<float_t> makeMirroredTransform2D(Axis axis)
{
    Transform2D<float_t> trans;
    if(axis == Axis::X) trans(1, 1) *= -1;
    else if(axis == Axis::Y) trans(0, 0) *= -1;
    else if(axis == Axis::Z) { trans(0, 0) *= -1; trans(1, 1) *= -1; }
    return trans;
}

///@brief makes a 3d shift transform matrix by a vector2d 
template <typename float_t, typename num_type>
inline Transform3D<float_t> makeShiftTransform3D(const Vector3D<num_type> & shift)
{
    Transform3D<float_t> trans;
    trans(0, 3) = float_t(shift[0]);
    trans(1, 3) = float_t(shift[1]);
    trans(2, 3) = float_t(shift[2]);
    return trans;
}

///@brief makes a 3d scale transform matrix by scale factor
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
    ///@brief constructs a transform2d model by matrix 3x3
    explicit Transform2D(const Matrix3x3<float_t> & m);

    ///@brief accesses matrix element by row and col index 0-2
    float_t & operator() (size_t row, size_t col);
    const float_t & operator() (size_t row, size_t col) const;

    ///@brief inverses the transform matrix, return false if the matrix det = 0
    bool Inverse();

    ///@brief gets result transform of this * trans
    Transform2D<float_t> operator * (const Transform2D<float_t> & trans) const;
    ///@brief gets result transform of trans * this 
    Transform2D<float_t> & Prod(const Transform2D<float_t> & trans);
    
    ///@brief gets transformed point of input
    template <typename num_type>
    Point2D<num_type> operator * (const Point2D<num_type> & point) const;

    ///@brief gets transformed segment of input
    template <typename num_type>
    Segment2D<num_type> operator * (const Segment2D<num_type> & segment) const;

    ///@brief gets transformed triangle of input
    template <typename num_type>
    Triangle2D<num_type> operator * (const Triangle2D<num_type> & triangle) const;

    ///@brief gets transformed box of input
    template <typename num_type>
    Polygon2D<num_type> operator * (const Box2D<num_type> & box) const;

    ///@brief gets transformed polyline of input
    template <typename num_type>
    Polyline2D<num_type> operator * (const Polyline2D<num_type> & polyline) const;

    ///@brief gets transformed polygon of input
    template <typename num_type>
    Polygon2D<num_type> operator * (const Polygon2D<num_type> & polygon) const;

    ///@brief gets transformed polygon with holes of input
    template <typename num_type>
    PolygonWithHoles2D<num_type> operator * (const PolygonWithHoles2D<num_type> & pwh) const;

    ///@brief accesses internal matrix data
    Matrix3x3<float_t> & GetMatrix() { return m_matrix; }
    const Matrix3x3<float_t> & GetMatrix() const { return m_matrix; }

    ///@brief gets matrix elements in array, adapt for OpenGL API
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
    ///@brief constructs a transform3d model from transform2d
    explicit Transform3D(const Transform2D<float_t> & trans2d);
    ///@brief constructs a transform3d model by matrix 4x4
    explicit Transform3D(const Matrix4x4<float_t> & m);

    ///@brief accesses matrix element by row and col index 0-2
    float_t & operator() (size_t row, size_t col);
    const float_t & operator() (size_t row, size_t col) const;

    ///@brief inverses the transform matrix, return false if the matrix det = 0
    bool Inverse();
    ///@brief gets transform2d(x, y) from transform3d(x, y, z)
    Transform2D<float_t> GetTransform2D() const;

    ///@brief gets result transform of this * trans
    Transform3D<float_t> operator * (const Transform3D<float_t> & trans) const;
    ///@brief gets result transform of trans * this 
    Transform3D<float_t> & Prod(const Transform3D<float_t> & tran);

    ///@brief gets transformed point of input
    template <typename num_type>
    Point3D<num_type> operator * (const Point3D<num_type> & point) const;

    ///@brief gets transformed segment of input
    template <typename num_type>
    Segment3D<num_type> operator * (const Segment3D<num_type> & segment) const;

    ///@brief gets transformed triangle of input
    template <typename num_type>
    Triangle3D<num_type> operator * (const Triangle3D<num_type> & triangle) const;

    ///@brief accesses internal matrix data
    Matrix4x4<float_t> & GetMatrix() { return m_matrix; }
    const Matrix4x4<float_t> & GetMatrix() const { return m_matrix; }

    ///@brief gets matrix elements in array, adapt for OpenGL API
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

    ///@brief accesses quaternion element by index 0-3
    float_t & operator[](size_t dim);
    const float_t & operator[](size_t dim) const;
    ///@brief multiples this with q
    Quaternion<float_t> & operator *= (const Quaternion<float_t> & q);
    ///@brief gets multiple result of this * q
    Quaternion<float_t> operator * (const Quaternion<float_t> & q) const;

    ///@brief performs a rotation around the axis at angle radians
    void SetAxisAndAngle(const Vector3D<float_t> & axis, float angle);
    ///@brief gets axis represented by this quaternion
    Vector3D<float_t> Axis() const;
    ///@brief gets angle represented by this quaternion
    float_t Angle() const;

    ///@brief multiplicative inverse of this quaternion
    void Invert();
    ///@brief negates all the coefficients of this quaternion
    void Negate();
    ///@brief normalizes the quaternion coefficients with unit quaternions
    void Normalize();
    ///@brief returns the magnitude of this           
    float_t Mag() const;
    ///@brief returns the squared magnitude of this
    float_t MagSqrt() const;
    Quaternion<float_t> Log() const;
    Quaternion<float_t> Exp() const;
    Quaternion<float_t> Inverse() const;
    Vector3D<float_t> Rotate(const Vector3D<float_t> & vec) const;
    Vector3D<float_t> InverseRotate(const Vector3D<float_t> & vec) const;

    ///@brief returns dot product of this and q
    float_t Dot(const Quaternion<float_t> & q) const;

    ///@brief accesses internal quat data
    boost::qvm::quat<float_t> & q() { return m_q; }
    const boost::qvm::quat<float_t> & q() const { return m_q; }

    static Quaternion<float_t> LnDif(const Quaternion<float_t> &a, const Quaternion<float_t> & b);
    static Quaternion<float_t> SquadTangent(const Quaternion<float_t> & before,
                                            const Quaternion<float_t> & center,
                                            const Quaternion<float_t> & after);
    static Quaternion<float_t> Squad(const Quaternion<float_t> & a, const Quaternion<float_t> & tgA,
                                     const Quaternion<float_t> & b, const Quaternion<float_t> & tgB, float_t t);

    ///@brief returns the result of spherical linear interpolation of the input quat `a`, `b` and interpolation parameter `t`
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
inline Transform2D<float_t> Transform3D<float_t>::GetTransform2D() const
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
        float_t coeff = std::acos(m_q.a[0]) / len;
        return Quaternion<float_t>(0, coeff * m_q.a[1], coeff * m_q.a[2], coeff * m_q.a[3]);
    }
}

template <typename float_t>
inline Quaternion<float_t> Quaternion<float_t>::Exp() const
{
    float_t theta = std::sqrt(m_q.a[1] * m_q.a[1] + m_q.a[2] * m_q.a[2] + m_q.a[3] * m_q.a[3]);
    if(theta < std::numeric_limits<float_t>::epsilon())
        return Quaternion<float_t>(std::cos(theta), m_q.a[1], m_q.a[2], m_q.a[3]);
    else{
        float_t coeff = std::sin(theta) / theta;
        return Quaternion<float_t>(std::cos(theta), coeff * m_q.a[1], coeff * m_q.a[2], coeff * m_q.a[3]);
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

///@brief transforms the point by the transform matrix
template <typename num_type, typename float_t>
inline void Transform(Point2D<num_type> & point, const Transform2D<float_t> & trans)
{
    VectorN<float_t, 3> vec = trans.GetMatrix() * VectorN<float_t, 3>(point[0], point[1], float_t(1));
    point = Point2D<num_type>(vec[0], vec[1]);
}

///@brief transforms the segment by the transform matrix
template <typename num_type, typename float_t>
inline void Transform(Segment2D<num_type> & segment, const Transform2D<float_t> & trans)
{
    Transform(segment[0], trans);
    Transform(segment[1], trans);
}

///@brief transforms the triangle by the transform matrix
template <typename num_type, typename float_t>
inline void Transform(Triangle2D<num_type> & triangle, const Transform2D<float_t> & trans)
{
    Transform(triangle[0], trans);
    Transform(triangle[1], trans);
    Transform(triangle[2], trans);
}

///@brief transforms the polyline by the transform matrix
template <typename num_type, typename float_t>
inline void Transform(Polyline2D<num_type> & polyline, const Transform2D<float_t> & trans)
{
    for(auto & point : polyline)
        Transform(point, trans);
}

///@brief transforms the polygon by the transform matrix
template <typename num_type, typename float_t>
inline void Transform(Polygon2D<num_type> & polygon, const Transform2D<float_t> & trans)
{
    for(auto iter = polygon.Begin(); iter != polygon.End(); ++iter)
        Transform(*iter, trans);
}

///@brief transforms the polygon with holes by the transform matrix
template <typename num_type, typename float_t>
inline void Transform(PolygonWithHoles2D<num_type> & pwh, const Transform2D<float_t> & trans)
{
    Transform(pwh.outline, trans);
    auto iter = pwh.BeginHoles();
    for(; iter != pwh.EndHoles(); ++iter)
        Transform(*iter, trans);
}

///@brief transforms the point by the transform matrix
template <typename num_type, typename float_t>
inline void Transform(Point3D<num_type> & point, const Transform3D<float_t> & trans)
{
    VectorN<float_t, 4> vec(trans.GetMatrix() * VectorN<float_t, 4>(point[0], point[1], point[2], float_t(1)));
    point = Point3D<num_type>(vec[0], vec[1], vec[2]);
}

///@brief transforms the segment by the transform matrix
template <typename num_type, typename float_t>
inline void Transform(Segment3D<num_type> & segment, const Transform3D<float_t> & trans)
{
    Transform(segment[0], trans);
    Transform(segment[1], trans);
}

///@brief transforms the triangle by the transform matrix
template <typename num_type, typename float_t>
inline void Transform(Triangle3D<num_type> & triangle, const Transform3D<float_t> & trans)
{
    Transform(triangle[0], trans);
    Transform(triangle[1], trans);
    Transform(triangle[2], trans);
}

///@brief transforms a collection of geometries by the transform matrix
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
} //namespace generic::geometry

namespace {

template <typename float_t>
inline std::ostream & operator<< (std::ostream & os, const generic::geometry::Transform2D<float_t> & trans)
{
    os << '[';
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            os << trans(i, j) << ',';
        }
        os << GENERIC_DEFAULT_EOL;
    }
    os << ']';
    return os;
}

template <typename float_t>
inline std::ostream & operator<< (std::ostream & os, const generic::geometry::Transform3D<float_t> & trans)
{
    os << '[';
    for (size_t i = 0; i < 4; ++i) {
        for (size_t j = 0; j < 4; ++j) {
            os << trans(i, j) << ',';
        }
        os << GENERIC_DEFAULT_EOL;
    }
    os << ']';
    return os;
}

} // namespace

namespace std {

template <typename float_t>
struct hash<generic::geometry::Transform2D<float_t>>
{
    size_t operator() (const generic::geometry::Transform2D<float_t> & t) const
    {
        return generic::hash::HashCombine(0, t(0, 0), t(0, 1), t(0, 2), 
                                             t(1, 0), t(1, 1), t(1, 2));
    }
};

template <typename float_t>
struct hash<generic::geometry::Transform3D<float_t>>
{
    size_t operator() (const generic::geometry::Transform3D<float_t> & t) const
    {
        return generic::hash::HashCombine(0, t(0, 0), t(0, 1), t(0, 2), t(0, 3), 
                                             t(1, 0), t(1, 1), t(1, 2), t(1, 3),
                                             t(2, 0), t(2, 1), t(2, 2), t(2, 3));
    }
};

} // namespace std