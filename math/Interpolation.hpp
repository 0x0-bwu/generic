/**
 * @file Interpolation.hpp
 * @author bwu
 * @brief Interpolation method implementation
 * @version 0.1
 * @date 2022-02-22
 */
#pragma once
#include "generic/math/la/Common.hpp"
#include "MathUtility.hpp"
#include <vector>
namespace generic::math {

template <typename RandomAccessContainer>
class Interpolation
{
public:
    enum class Method { LINEAR = 1, CUBIC = 2};
    enum class BCType { FIRST_DERIV = 1, SECOND_DERIV = 2, NOT_A_KNOT = 3};
    using Scalar = typename RandomAccessContainer::value_type;
    static_assert(std::is_floating_point_v<Scalar>, "only floating point type supported!");

    ///@brief constructs an interpolation object with method
    Interpolation(Method method = Method::LINEAR);
    /**
     * @brief constructs an interpolation object with parameters
     * @param x input samples of x
     * @param y input samples of y
     * @param method interpolation method
     * @note if Cubic method failed with input x, y samples, it will decay to Linear method
     * @param monotonic set whether curve monotonic when interpolation
     * @param left left boundary condition type
     * @param lValue left boundary condition initial value
     * @param right right boundary condition type
     * @param rValue right boundary condition inital value
     */
    Interpolation(RandomAccessContainer x, RandomAccessContainer y, Method method, bool monotonic = false,
                    BCType left  = BCType::SECOND_DERIV, Scalar lValue = 0,
                    BCType right = BCType::SECOND_DERIV, Scalar rValue = 0);

    ///@brief sets boundary conditions and inital values
    void SetBoundary(BCType left, BCType right, Scalar lVal, Scalar rVal);

    ///@brief sets input samples of x and y
    void SetSamples(RandomAccessContainer x, RandomAccessContainer y);

    ///@brief gets interpolated value of input x
    Scalar operator() (Scalar x) const { return Interpolate(x); }

private:
    void UpdateCoefficients();
    bool makeMonotonic();
    size_t FindClosest(Scalar x) const;
    Scalar Interpolate(Scalar x) const;

private:
    Method m_method;
    bool m_monotonic{false};
    BCType m_left{BCType::SECOND_DERIV};
    BCType m_right{BCType::SECOND_DERIV};
    Scalar m_lVal{0}, m_rVal{0}, m_coef0{0};
    std::vector<Scalar> m_coef1, m_coef2, m_coef3;
    RandomAccessContainer m_x, m_y;
};

template <typename RandomAccessContainer>
inline Interpolation<RandomAccessContainer>::Interpolation(Method method)
 : m_method(method)
{
}

template <typename RandomAccessContainer>
inline Interpolation<RandomAccessContainer>::Interpolation(RandomAccessContainer x, RandomAccessContainer y, Method method, bool monotonic,
                                                           BCType left, Scalar lValue, BCType right, Scalar rValue)
 : m_method(method), m_monotonic(monotonic), m_left(left), m_right(right), m_lVal(lValue), m_rVal(rValue)
{
    SetSamples(std::move(x), std::move(y));
}

template <typename RandomAccessContainer>
inline void Interpolation<RandomAccessContainer>::SetBoundary(BCType left, BCType right, Scalar lVal, Scalar rVal)
{
    GENERIC_ASSERT(m_x.size() == 0/*should set boundary type before setting samples*/);
    m_left = left; m_right = right;
    m_lVal = lVal; m_rVal = rVal;
}

template <typename RandomAccessContainer>
inline void Interpolation<RandomAccessContainer>::SetSamples(RandomAccessContainer x, RandomAccessContainer y)
{
    const size_t size = x.size();
    GENERIC_ASSERT(size > 0 && size == y.size());
    std::swap(x, m_x);
    std::swap(y, m_y);
    m_coef1.assign(size, 0);
    m_coef2.assign(size, 0);
    m_coef3.assign(size, 0);
    if (size < 2) return;
    if (size < 3) m_method = Method::LINEAR;
    if (Method::LINEAR == m_method) {
        for(size_t i = 0, j = 1; i < size - 1; ++i, ++j)
            m_coef1[i] = Scalar(m_y[j] - m_y[i]) / Scalar(m_x[j] - m_x[i]);
        m_coef1[size - 1] = m_coef1[size - 2];
    }
    else if (Method::CUBIC == m_method){
        la::DenseMatrix<Scalar> m(size, size);
        la::DenseVector<Scalar> rhs(size);
        for(size_t i = 0, j = 1, k = 2; k < size; ++i, ++j, ++k){
            m(j, i) = 1.0 / 3.0 * (m_x[j] - m_x[i]);
            m(j, j) = 2.0 / 3.0 * (m_x[k] - m_x[i]);
            m(j, k) = 1.0 / 3.0 * (m_x[k] - m_x[j]);
            rhs[j] = (m_y[k] - m_y[j]) / (m_x[k] - m_x[j]) - (m_y[j] - m_y[i]) / (m_x[j] - m_x[i]);
        }
        //BC
        if (m_left == BCType::FIRST_DERIV) {
            m(0, 0) = 2.0 * (m_x[1] - m_x[0]);
            m(0, 1) = 1.0 * (m_x[1] - m_x[0]);
            rhs[0] = 3.0 * ((m_y[1] - m_y[0]) / (m_x[1] - m_x[0]) - m_lVal);
        }
        else if (m_left == BCType::SECOND_DERIV) {
            m(0, 0) = 2.0;
            m(0, 1) = 0.0;
            rhs[0] = m_lVal;
        }
        else if (m_left == BCType::NOT_A_KNOT) {
            m(0, 0) = - m_x[2] + m_x[1];
            m(0, 1) =   m_x[2] - m_x[0];
            m(0, 2) = - m_x[1] + m_x[0];
            rhs[0] = 0.0;
        }
        else { GENERIC_ASSERT(false); }
    
        if (m_right == BCType::FIRST_DERIV) {
            m(size - 1, size - 1) = 2.0 * (m_x[size - 1] - m_x[size - 2]);
            m(size - 1, size - 2) = 1.0 * (m_x[size - 1] - m_x[size - 2]);
            rhs[size - 1] = 3.0 * (m_rVal - (m_y[size - 1] - m_y[size - 2]) / (m_x[size - 1] - m_x[size - 2]));
        }
        else if (m_right == BCType::SECOND_DERIV) {
            m(size - 1, size - 1) = 2.0;
            m(size - 1, size - 2) = 0.0;
            rhs[size - 1] = m_rVal;
        }
        else if (m_right == BCType::NOT_A_KNOT) {
            m(size - 1, size - 3) = - m_x[size - 1] + m_x[size - 2];
            m(size - 1, size - 2) =   m_x[size - 1] - m_x[size - 3];
            m(size - 1, size - 1) = - m_x[size - 2] + m_x[size - 3]; 
        }
        else { GENERIC_ASSERT(false); }

        la::VectorView<la::DenseVector<Scalar>> b(m_coef2.data(), m_coef2.size());
        b = m.fullPivLu().solve(rhs);
        
        for (size_t i = 0, j = 1; j < size; ++i, ++j) {
            m_coef1[i] = (m_y[j] - m_y[i]) / (m_x[j] - m_x[i]) - 1.0 / 3.0 * (2.0 * m_coef2[i] + m_coef2[j]) * (m_x[j] - m_x[i]);
            m_coef3[i] = 1.0 / 3.0 * (m_coef2[j] - m_coef2[i]) / (m_x[j] - m_x[i]);
        }

        Scalar h = m_x[size - 1] - m_x[size - 2];
        m_coef3[size - 1] = 0.0;
        m_coef1[size - 1] = 3.0 * m_coef3[size - 2] * h * h + 2.0 * m_coef2[size - 2] * h + m_coef1[size - 2];
        if (m_right == BCType::FIRST_DERIV) m_coef2[size - 1] = 0.0;
    }
    m_coef0 = (m_left == BCType::FIRST_DERIV) ? 0 : m_coef3[0];
    
    if (not m_monotonic && size > 2) makeMonotonic();
}

template <typename RandomAccessContainer>
inline void Interpolation<RandomAccessContainer>::UpdateCoefficients()
{
    size_t size = m_coef1.size();
    for (size_t i = 0; i < size - 1; ++i) {
        const float_t h = m_x[i + 1] - m_x[i];
        m_coef2[i] = (3.0 * (m_y[i + 1] - m_y[i]) / h - (2.0 * m_coef1[i] + m_coef1[i + 1])) / h;
        m_coef3[i] = ((m_coef1[i + 1] - m_coef1[i]) / (3.0 * h) - 2.0 / 3.0 * m_coef2[i]) / h;
    }
    m_coef0 = (m_left == BCType::FIRST_DERIV) ? 0 : m_coef3[0];
}

template <typename RandomAccessContainer>
inline bool Interpolation<RandomAccessContainer>::makeMonotonic()
{
    GENERIC_ASSERT(m_x.size() > 2);
    bool modified = false;
    for (size_t i = 0; i < m_x.size(); ++i) {
        size_t im1 = std::max(i - 1, size_t(0));
        size_t ip1 = std::min(i + 1, m_x.size() - 1);
        if ((m_y[im1] <= m_y[i] && m_y[i] <= m_y[ip1] && m_coef1[i] < 0) ||
            (m_y[im1] >= m_y[i] && m_y[i] >= m_y[ip1] && m_coef1[i] > 0)) {
            modified = true;
            m_coef1[i] = 0;
        }
    }

    for (size_t i = 0; i < m_x.size() - 1; ++i) {
        auto avg = (m_y[i + 1] - m_y[i]) / (m_x[i + 1] - m_x[i]);
        if (EQ<Scalar>(avg, 0) && (not EQ<Scalar>(m_coef1[i], 0) or not EQ<Scalar>(m_coef1[i], 0))) {
            modified = true;
            m_coef1[i] = 0;
            m_coef1[i + 1] = 0;
        }
        else if((GE<Scalar>(m_coef1[i], 0) && GE<Scalar>(m_coef1[i + 1], 0) && GE<Scalar>(avg, 0)) ||
                (LE<Scalar>(m_coef1[i], 0) && LE<Scalar>(m_coef1[i + 1], 0) && LE<Scalar>(avg, 0))) {
            Scalar r = std::sqrt(m_coef1[i] * m_coef1[i] + m_coef1[i + 1] * m_coef1[i + 1]) / std::fabs(avg);
            if (r > 3.0) {
                modified = true;
                m_coef1[i] *= (3.0 / r);
                m_coef1[i + 1] *= (3.0 / r);
            }
        }
    }

    if(modified == true){
        UpdateCoefficients();
        m_monotonic = true;
    }
    return modified;
}

template <typename RandomAccessContainer>
inline size_t Interpolation<RandomAccessContainer>::FindClosest(Scalar x) const
{
    auto iter = std::upper_bound(m_x.begin(), m_x.end(), x);
    size_t index = std::max(int(iter - m_x.begin() - 1), 0);
    return index;
}

template <typename RandomAccessContainer>
inline typename Interpolation<RandomAccessContainer>::Scalar
Interpolation<RandomAccessContainer>::Interpolate(Scalar x) const
{
    size_t n = m_x.size();
    size_t idx = FindClosest(x);
    Scalar h = x - m_x[idx];
    if (x < m_x[0]) return (m_coef0 * h + m_coef1[0]) * h + m_y[0];
    else if (x > m_x[n - 1]) return (m_coef2[n-1] * h + m_coef1[n - 1]) * h + m_y[n - 1];
    else return ((m_coef3[idx] * h + m_coef2[idx]) * h + m_coef1[idx]) * h + m_y[idx];        
}

}//namespace generic::math
