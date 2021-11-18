#ifndef GENERIC_MATH_LINEARALGEBRA_HPP
#define GENERIC_MATH_LINEARALGEBRA_HPP
#include "generic/common/Traits.hpp"
#include "generic/common/Version.hpp"
#include "MathUtility.hpp"
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/hermitian.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <initializer_list>
#include <cmath>
namespace generic{
namespace math{
namespace la {
using namespace generic::common;
using namespace generic::math;
using namespace boost::numeric;

template <typename num_type>
inline ublas::unbounded_array<num_type> makeUbArray(std::initializer_list<num_type> l) {
    ublas::unbounded_array<num_type> res(l.size());
    for(size_t i = 0; i < l.size(); ++i)
        res[i] = *(l.begin() + i);
    return res;
}

template <typename num_type, typename iterator>
inline ublas::unbounded_array<num_type> makeUbArray(iterator begin, iterator end) {
    size_t size = std::distance(begin, end);
    ublas::unbounded_array<num_type> res(size);
    for(size_t i = 0; i < size; ++i)
        res[i] = *(begin + i);
    return res;
}

template <typename num_type, size_t M, size_t N> class Matrix;
template <typename num_type, size_t N>
class Vector
{
public:
    using data_type = ublas::vector<num_type>;
    Vector() : m_data(N) { std::fill(m_data.begin(), m_data.end(), 0); }
    Vector(num_type s) : m_data(N) { std::fill(m_data.begin(), m_data.end(), s); }
    Vector(data_type data) : m_data(data) { m_data.resize(N, true); }
    Vector(std::initializer_list<num_type> l) : m_data(makeUbArray(l)) { m_data.resize(N, true); }

    template <typename... Args>
    Vector(num_type first, num_type second, Args... args) : m_data(N) { Set(first, second, args...); }

    bool operator== (const Vector<num_type, N> & v) { return  Same(*this, v); }
    bool operator!= (const Vector<num_type, N> & v) { return !Same(*this, v); }

    num_type & operator() (size_t dim) { return m_data(dim); }
    const num_type & operator() (size_t dim) const { return m_data(dim); }

    num_type & operator[] (size_t dim) { return m_data(dim); }
    const num_type & operator[] (size_t dim) const { return m_data(dim); }

    Vector operator+ (const Vector<num_type, N> & v) const;
    Vector operator- () const;
    Vector operator- (const Vector<num_type, N> & v) const;
    Vector operator* (float_type<num_type> scale) const;
    Vector operator/ (float_type<num_type> scale) const;   
    Vector operator* (const Matrix<num_type, N, N> & m) const; 

    Vector & operator += (const Vector<num_type, N> & v);
    Vector & operator -= (const Vector<num_type, N> & v);
    Vector & operator *= (float_type<num_type> scale);
    Vector & operator /= (float_type<num_type> scale);

    void Swap(Vector<num_type, N> & v) { m_data.swap(v.Data()); }
    Matrix<num_type, N, 1> T() const;
    num_type Norm1() const;
    float_type<num_type> Norm2() const;
    num_type NormSquare() const;

    data_type & Data() { return m_data; }
    const data_type & Data() const { return m_data; }

    template <typename... Args>
    void Set(Args... args);

    static bool Same(const Vector<num_type, N> & v1, const Vector<num_type, N> & v2);
private:
    data_type m_data;
};

struct matrix_general_tag {};
struct matrix_identity_tag {};
struct matrix_diagonal_tag {};
struct matrix_upper_triangle_tag {};
struct matrix_lower_triangle_tag {};
struct matrix_symmetric_tag {};
struct matrix_hermitian_tag {};
struct matrix_transform2d_tag {};
struct matrix_transform3d_tag {};
struct matrix_unknown_tag {};

template <typename num_type, size_t M = 1, size_t N = 1>
class Matrix
{
public:
    using data_type = ublas::matrix<num_type>;
    Matrix();

    Matrix(num_type s);

    template<typename matrix_tag, 
             typename std::enable_if<
                      std::is_same<matrix_tag, matrix_diagonal_tag>::value  && M == N, bool>::type = true>
    Matrix(num_type s, matrix_tag);

    Matrix(data_type data) { data.resize(M, N, true); m_data = data; }

    Matrix(const std::initializer_list<std::initializer_list<num_type> > & ll)
    {
        m_data.resize(M, N);
        auto iter1 = ll.begin();
        for(size_t i = 0; i < M && iter1 != ll.end(); ++i, ++iter1){
            const std::initializer_list<num_type> & l = *iter1;
            auto iter2 = l.begin();
            for(size_t j = 0; j < N && iter2 != l.end(); ++j, ++iter2)
                m_data(i, j) = *iter2;
        }
    }

    template <typename... Args>
    Matrix(Vector<num_type, N> first, Vector<num_type, N> second, Args... args);

    bool operator== (const Matrix<num_type, M, N> & m) { return  Same(*this, m); }
    bool operator!= (const Matrix<num_type, M, N> & m) { return !Same(*this, m); }

    num_type & operator() (size_t row, size_t col);
    const num_type & operator() (size_t row, size_t col) const;

    Matrix operator+ (const Matrix<num_type, M, N> & m) const;
    Matrix operator- () const;
    Matrix operator- (const Matrix<num_type, M, N> & m) const;
    Matrix operator* (float_type<num_type> scale) const;
    Matrix operator/ (float_type<num_type> scale) const;
    
    Vector<num_type, N> operator* (const Vector<num_type, N> & v) const;

    template <size_t O>
    Matrix<num_type, M, O> operator* (const Matrix<num_type, N, O> & m) const;

    Matrix & operator+= (const Matrix<num_type, M, N> & m);
    Matrix & operator-= (const Matrix<num_type, M, N> & m);
    Matrix & operator*= (float_type<num_type> scale);
    Matrix & operator/= (float_type<num_type> scale);
    Matrix & operator*= (const Matrix<num_type, N, N> & m);

    Vector<num_type, N> Row(size_t row) const;
    Vector<num_type, M> Col(size_t col) const;

    void Swap(Matrix<num_type, M, N> & m) { m_data.swap(m.Data()); }
    Matrix<num_type, N, M> T() const;
    float_type<num_type> Norm2() const { return std::sqrt(NormSquare()); }
    num_type NormSquare() const;

    data_type & Data() { return m_data; }
    const data_type & Data() const { return m_data; }

    template <typename... Args>
    void Set(Args... args);

    static bool Same(const Matrix<num_type, M, N> & m1, const Matrix<num_type, M, N> & m2);
private:
    data_type m_data;
};

template <typename num_type, size_t I, size_t N>
struct VectorSetter
{
    template <typename... Args>
    static void Set(Vector<num_type, N> & v, num_type s, Args... args)
    {
        v[I] = s;
        VectorSetter<num_type, I + 1, N>::Set(v, args...);
    }
};

template <typename num_type, size_t N>
struct VectorSetter<num_type, N, N>
{
    static void Set(Vector<num_type, N> &) {}
};

#if GENERIC_CURRENT_BOOST_LIBRARY_VER >= 172
template <typename num_type, size_t N>
inline num_type Vector<num_type, N>::Norm1() const
{
    return norm_1(m_data);
}

template <typename num_type, size_t N>
inline float_type<num_type> Vector<num_type, N>::Norm2() const
{
    return norm_2(m_data);
}

template <typename num_type, size_t N>
inline num_type Vector<num_type, N>::NormSquare() const
{
    return norm_2_square(m_data);
}
#else
template <typename num_type, size_t N>
inline num_type Vector<num_type, N>::Norm1() const
{
    num_type norm1(0);
    for(size_t i = 0; i < N; ++i)
        norm1 += std::abs(m_data[i]);
    return norm1;
}

template <typename num_type, size_t N>
inline float_type<num_type> Vector<num_type, N>::Norm2() const
{
    return std::sqrt(float_type<num_type>(NormSquare()));
}

template <typename num_type, size_t N>
inline num_type Vector<num_type, N>::NormSquare() const
{
    num_type normSquare(0);
    for(size_t i = 0; i < N; ++i)
        normSquare += m_data[i] * m_data[i];
    return normSquare;
}
#endif

template <typename num_type, size_t N>
template <typename... Args>
inline void Vector<num_type, N>::Set(Args... args)
{
    VectorSetter<num_type, 0, N>::Set(*this, num_type(args)...);
}

template <typename num_type, size_t N>
inline Vector<num_type, N> Vector<num_type, N>::operator+ (const Vector<num_type, N> & v) const
{
    return Vector<num_type, N>(m_data + v.Data());
}

template <typename num_type, size_t N>
inline Vector<num_type, N> Vector<num_type, N>::operator- () const
{
    return Vector<num_type, N>(0) - *this;
}

template <typename num_type, size_t N>
inline Vector<num_type, N> Vector<num_type, N>::operator- (const Vector<num_type, N> & v) const
{
    return Vector<num_type, N>(m_data - v.Data());
}

template <typename num_type, size_t N>
inline Vector<num_type, N> Vector<num_type, N>::operator* (float_type<num_type> scale) const
{
    return Vector<num_type, N>(m_data * scale);
}

template <typename num_type, size_t N>
inline Vector<num_type, N> Vector<num_type, N>::operator/ (float_type<num_type> scale) const
{
    return Vector<num_type, N>(m_data / scale);
}

template <typename num_type, size_t N>
inline Vector<num_type, N> Vector<num_type, N>::operator* (const Matrix<num_type, N, N> & m) const
{
    return Vector<num_type, N>(prod(m_data, m.Data()));
}

template <typename num_type, size_t N>
inline Vector<num_type, N> & Vector<num_type, N>::operator += (const Vector<num_type, N> & v)
{
    m_data += v.Data();
    return *this;
}

template <typename num_type, size_t N>
inline Vector<num_type, N> & Vector<num_type, N>::operator -= (const Vector<num_type, N> & v)
{
    m_data -= v.Data();
    return *this;
}

template <typename num_type, size_t N>
inline Vector<num_type, N> & Vector<num_type, N>::operator *= (float_type<num_type> scale)
{
    m_data *= scale;
    return *this;
}

template <typename num_type, size_t N>
inline Vector<num_type, N> & Vector<num_type, N>::operator /= (float_type<num_type> scale)
{
    m_data /= scale;
    return *this;
}

template <typename num_type, size_t N>
inline Matrix<num_type, N, 1> Vector<num_type, N>::T() const
{
    Matrix<num_type, N, 1> m;
    for(size_t i = 0; i < N; ++i)
        m(i, 0) = m_data(i);
    return m;
}

template <typename num_type, size_t N>
inline bool Vector<num_type, N>::Same(const Vector<num_type, N> & v1, const Vector<num_type, N> & v2)
{
    return math::EQ((v1 - v2).NormSquare(), num_type(0));
}

template <typename num_type, size_t I, size_t M, size_t N>
struct MatrixSetter
{
    template <typename... Args>
    static void Set(Matrix<num_type, M, N> & m, Vector<num_type, N> v, Args... args)
    {
        for(size_t j = 0; j < N; ++j)
            m(I, j) = v[j];
        MatrixSetter<num_type, I + 1, M, N>::Set(m, args...);
    }
};

template <typename num_type, size_t M, size_t N>
struct MatrixSetter<num_type, M, M, N>
{
    static void Set(Matrix<num_type, M, N> &) {}
};

template <typename num_type, size_t M, size_t N>
inline Matrix<num_type, M, N>::Matrix() : m_data(M, N)
{
    for(size_t i = 0; i < M; ++i)
        for(size_t j = 0; j < N; ++j)
            m_data(i, j) = 0;
}

template <typename num_type, size_t M, size_t N>
inline Matrix<num_type, M, N>::Matrix(num_type s) : Matrix()
{
    for(size_t i = 0; i < M; ++i)
        for(size_t j = 0; j < N; ++j)
            m_data(i, j) = s;
}

template <typename num_type, size_t M, size_t N>
template <typename matrix_tag, 
          typename std::enable_if<
                   std::is_same<matrix_tag, matrix_diagonal_tag>::value && M == N, bool>::type>
inline Matrix<num_type, M, N>::Matrix(num_type s, matrix_tag) : Matrix()
{
    for(size_t i = 0; i < M; ++i)
        m_data(i, i) = s;
}

template <typename num_type, size_t M, size_t N>
template <typename... Args>
inline Matrix<num_type, M, N>::Matrix(Vector<num_type, N> first, Vector<num_type, N> second, Args... args)
 : m_data(M, N)
{
    Set(first, second, args...);
}

template <typename num_type, size_t M, size_t N>
inline num_type & Matrix<num_type, M, N>::operator() (size_t row, size_t col) { return m_data(row, col); }

template <typename num_type, size_t M, size_t N>
inline const num_type & Matrix<num_type, M, N>::operator() (size_t row, size_t col) const { return m_data(row, col); }

template <typename num_type, size_t M, size_t N>
inline Matrix<num_type, M, N> Matrix<num_type, M, N>::operator+ (const Matrix<num_type, M, N> & m) const
{
    return Matrix<num_type, M, N>(m_data + m.Data());
}

template <typename num_type, size_t M, size_t N>
inline Matrix<num_type, M, N> Matrix<num_type, M, N>::operator- () const
{
    return Matrix<num_type, M, N>(0) - *this;
}

template <typename num_type, size_t M, size_t N>
inline Matrix<num_type, M, N> Matrix<num_type, M, N>::operator- (const Matrix<num_type, M, N> & m) const
{
    return Matrix<num_type, M, N>(m_data - m.Data());
}

template <typename num_type, size_t M, size_t N>
inline Matrix<num_type, M, N> Matrix<num_type, M, N>::operator* (float_type<num_type> scale) const
{
    return Matrix<num_type, M, N>(m_data * scale);
}

template <typename num_type, size_t M, size_t N>
inline Matrix<num_type, M, N> Matrix<num_type, M, N>::operator/ (float_type<num_type> scale) const
{
    return Matrix<num_type, M, N>(m_data / scale);
}
template <typename num_type, size_t M, size_t N>
inline Vector<num_type, N> Matrix<num_type, M, N>:: operator* (const Vector<num_type, N> & v) const
{
    return Vector<num_type, N>(prod(m_data, v.Data()));
}

template <typename num_type, size_t M, size_t N>
template <size_t O>
inline Matrix<num_type, M, O> Matrix<num_type, M, N>::operator* (const Matrix<num_type, N, O> & m) const
{
    return Matrix<num_type, M, O>(prod(m_data, m.Data()));
}  

template <typename num_type, size_t M, size_t N>
inline Matrix<num_type, M, N> & Matrix<num_type, M, N>::operator += (const Matrix<num_type, M, N> & m)
{
    m_data += m.Data();
    return *this;
}

template <typename num_type, size_t M, size_t N>
inline Matrix<num_type, M, N> & Matrix<num_type, M, N>::operator -= (const Matrix<num_type, M, N> & m)
{
    m_data -= m.Data();
    return *this;
}

template <typename num_type, size_t M, size_t N>
inline Matrix<num_type, M, N> & Matrix<num_type, M, N>::operator *= (float_type<num_type> scale)
{
    m_data *= scale;
    return *this;
}

template <typename num_type, size_t M, size_t N>
inline Matrix<num_type, M, N> & Matrix<num_type, M, N>::operator /= (float_type<num_type> scale)
{
    m_data /= scale;
    return *this;
}

template <typename num_type, size_t M, size_t N>
inline Matrix<num_type, M, N> & Matrix<num_type, M, N>::operator *= (const Matrix<num_type, N, N> & m)
{
    m_data = prod(m_data, m.Data());
    return *this;
}

template <typename num_type, size_t M, size_t N>
inline Vector<num_type, N> Matrix<num_type, M, N>::Row(size_t row) const
{
    Vector<num_type, N> v;
    for(size_t i = 0; i < N; ++i)
        v[i] = m_data(row, i);
    return v;
}

template <typename num_type, size_t M, size_t N>
inline Vector<num_type, M> Matrix<num_type, M, N>::Col(size_t col) const
{
    Vector<num_type, M> v;
    for(size_t j = 0; j < M; ++j)
        v[j] = m_data(j, col);
    return v;
}

template <typename num_type, size_t M, size_t N>
Matrix<num_type, N, M> Matrix<num_type, M, N>::T() const
{
    Matrix<num_type, N, M> mT;
    for(size_t i = 0; i < M; ++i)
        for(size_t j = 0; j < N; ++j)
            mT(j, i) = m_data(i, j);
    return mT;
}

template <typename num_type, size_t M, size_t N>
template <typename... Args>
inline void Matrix<num_type, M, N>::Set(Args... args)
{
    MatrixSetter<num_type, 0, M, N>::Set(*this, Vector<num_type, N>(args)...);
}

template <typename num_type, size_t M, size_t N>
inline num_type Matrix<num_type, M, N>::NormSquare() const
{
    num_type res(0);
    for(size_t i = 0; i < M; ++i)
        for(size_t j = 0; j < N; ++j)
            res += m_data(i, j) * m_data(i, j);
    return res;
}

template <typename num_type, size_t M, size_t N>
inline bool Matrix<num_type, M, N>::Same(const Matrix<num_type, M, N> & m1, const Matrix<num_type, M, N> & m2)
{
    return math::EQ((m1 - m2).NormSquare(), num_type(0));
}
}//namespace la
}//namespace math
}//namespace generic
#endif//GENERIC_MATH_LINEARALGEBRA_HPP