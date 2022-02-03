#ifndef GENERIC_GEOMETRY_VECTOR_HPP
#define GENERIC_GEOMETRY_VECTOR_HPP
#include "generic/math/MathUtility.hpp"
#include "Common.hpp"
#include "Point.hpp"
#include <array>
namespace generic  {
namespace geometry {

template <typename num_type>
using Vector2D = Point2D<num_type>;

template <typename num_type>
inline Vector2D<float_type<num_type> > Normalize(const Vector2D<num_type> & vec)
{
    float_type<num_type> norm = std::sqrt(vec.NormSquare());
    float_type<num_type> epsilon = std::numeric_limits<float_type<num_type> >::epsilon();
    if(norm < epsilon) return Vector2D<float_type<num_type> >();
    return vec.template Cast<float_type<num_type> >() / norm;
}

template <typename num_type>
using Vector3D = Point3D<num_type>;

template <typename num_type>
inline Vector3D<float_type<num_type> > Normalize(const Vector3D<num_type> & vec)
{
    float_type<num_type> norm = std::sqrt(vec.NormSquare());
    float_type<num_type> epsilon = std::numeric_limits<float_type<num_type> >::epsilon();
    if(norm < epsilon) return Vector3D<float_type<num_type> >();
    return vec.template Cast<float_type<num_type> >() / norm;
}

template <typename num_type, size_t N>
class VectorN
{
public:
    using coor_t = num_type;
    static const size_t dim = N;
    VectorN() { std::fill(m_data.begin(), m_data.end(), 0); }
    VectorN(num_type s) { std::fill(m_data.begin(), m_data.end(), s); }

    template <typename... Args>
    VectorN(num_type first, num_type second, Args... args) { Set(first, second, args...); }

    num_type& operator[](int dim);
    const num_type& operator[](int dim) const;

    bool operator== (const VectorN<num_type, N> & v) const;
    bool operator!= (const VectorN<num_type, N> & v) const;

    template <typename... Args>
    void Set(Args... args);

    num_type Dot(const VectorN<num_type, N> & v);

    static num_type NormSquare(const VectorN<num_type, N> & v1, const VectorN<num_type, N> & v2);

private:
    std::array<num_type, N> m_data;
};

template <typename num_type, size_t I, size_t N>
struct VectorNSetter
{
    template <typename... Args>
    static void Set(VectorN<num_type, N> & v, num_type s, Args... args)
    {
        v[I] = s;
        VectorNSetter<num_type, I + 1, N>::Set(v, args...);
    }
};

template <typename num_type, size_t N>
struct VectorNSetter<num_type, N, N>
{
    static void Set(VectorN<num_type, N> &) {}
};

template <typename num_type, size_t N>
inline num_type & VectorN<num_type, N>::operator[](int dim) { return m_data[dim]; }

template <typename num_type, size_t N>
inline const num_type & VectorN<num_type, N>::operator[](int dim) const { return m_data[dim]; }

template <typename num_type, size_t N>
inline bool VectorN<num_type, N>::operator== (const VectorN<num_type, N> & v) const
{
    return math::EQ(NormSquare(*this, v), num_type(0));
}

template <typename num_type, size_t N>
inline bool VectorN<num_type, N>::operator!= (const VectorN<num_type, N> & v) const
{
    return !(*this == v);
}

template <typename num_type, size_t N>
template <typename... Args>
inline void VectorN<num_type, N>::Set(Args... args)
{
    VectorNSetter<num_type, 0, N>::Set(*this, num_type(args)...);
}

template <typename num_type, size_t N>
inline num_type VectorN<num_type, N>::Dot(const VectorN<num_type, N> & v)
{
    num_type res(0);
    for(size_t i = 0; i < N; ++i)
        res += (*this)[i] * v[i];
    return res;
}

template <typename num_type, size_t N>
inline num_type VectorN<num_type, N>::NormSquare(const VectorN<num_type, N> & v1, const VectorN<num_type, N> & v2)
{
    num_type normSquare(0);
    for(size_t i = 0; i < N; ++i)
        normSquare += (v1[i] - v2[i]) * (v1[i] - v2[i]);
    return normSquare;
}

}// namespace geometry
}// namespace generic
#endif//GENERIC_GEOMETRY_VECTOR_HPP