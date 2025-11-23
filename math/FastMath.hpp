/**
 * @file FastMath.hpp
 * @author bwu
 * @brief fast math operations from Paul Mineiro's FastFloat
 * @version 0.1
 * @date 2024-10-12
 */
#pragma once
#include "generic/common/Exception.hpp"
#include <type_traits>
#include <algorithm>
#include <iterator>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <vector>
#include <cmath>

namespace generic::math {

///@brief "Fast" sine approximation, valid for x in [-PI, PI], max abs error about 4e-05
inline float FastSin(float x)
{
    static constexpr float q = 0.78444488374548933f;
    union { float f; uint32_t i; } p = { 0.20363937680730309f };
    union { float f; uint32_t i; } r = { 0.015124940802184233f };
    union { float f; uint32_t i; } s = { -0.0032225901625579573f };

    union { float f; uint32_t i; } vx = { x };
    uint32_t sign = vx.i & 0x80000000;
    vx.i = vx.i & 0x7FFFFFFF;

    float qpprox = 1.2732395447351627f * x - 0.40528473456935109f * x * vx.f;
    float qpproxsq = qpprox * qpprox;

    p.i |= sign;
    r.i |= sign;
    s.i ^= sign;

    return q * qpprox + qpproxsq * (p.f + qpproxsq * (r.f + qpproxsq * s.f));
}

///@brief "Fast" cosine approximation, valid for x in [-PI, PI], max abs error about 4e-05
inline float FastCos(float x)
{
    float offset = (x > 1.5707963267948966f) ? -4.7123889803846899f : 1.5707963267948966f;
    return FastSin(x + offset);
}

///@brief "Fast" log base 2 approximation, valid for positive x as precision allows, max abs error about 2e-4
inline float FastLog2(float x)
{
    union { float f; uint32_t i; } vx = { x };
    union { uint32_t i; float f; } mx = { (vx.i & 0x007FFFFF) | 0x3f000000 };
    float y = vx.i;
    y *= 1.1920928955078125e-7f;
    return y - 124.22551499f - 1.498030302f * mx.f - 1.72587999f / (0.3520887068f + mx.f);
}

///@brief "Faster" log base 2 approximation, valid for positive x as precision allows, max abs error about 6e-2
inline float FasterLog2(float x)
{
    union { float f; uint32_t i; } vx = { x };
    float y = (float)(vx.i);
    y *= 1.1920928955078125e-7f;
    return y - 126.94269504f;
}

///@brief "Fast" natural logarithm approximation, valid for positive x as precision allows, max abs error about 2e-4
inline float FastLog(float x)
{
    return 0.6931471805599453094f * FastLog2(x);
}

///@brief "Fast" natural logarithm approximation, valid for positive x as precision allows, max abs error about 4e-2
inline float FasterLog(float x)
{
  return 0.6931471805599453094f * FasterLog2(x);
}

///@brief "Fast" power of 2 approximation, valid for x in [ -126, 0 ] as precision allows, max abs error about 4e-5
inline float FastPow2(float p)
{
    float clipp = (p < -126) ? -126.0f : p;
    int w = clipp;
    float z = clipp - w + 1.f;
    union { uint32_t i; float f; } v = { (uint32_t) ( (1 << 23) * 
        (clipp + 121.2740575f + 27.7280233f / (4.84252568f - z) - 1.49012907f * z)
        ) };
    return v.f;
}

///@brief "Faster" power of 2 approximation, valid for x in [ -126, 0 ] as precision allows, max abs error about 3e-2
inline float FasterPow2(float p)
{
    float clipp = (p < -126) ? -126.0f : p;
    union { uint32_t i; float f; } v = { (uint32_t)( (1 << 23) * (clipp + 126.94269504f) ) };
    return v.f;
}

///@brief "Fast" exponential approximation, valid for x in [ ~ -87, 0 ] as precision allows, max abs error about 4e-5
inline float FastExp(float p)
{
    return FastPow2(1.442695040f * p);
}

///@brief "Faster" exponential approximation, valid for x in [ ~ -87, 0 ] as precision allows, max abs error about 3e-2
inline float FasterExp(float p)
{
    return FasterPow2(1.442695040f * p);
}

template <typename Float>
class FastPchip
{
public:
    FastPchip() = default;
    // x and y must be same size >= 2, and x strictly increasing
    FastPchip(std::vector<Float> && x, std::vector<Float> && y)
    {
        Reset(std::move(x), std::move(y));
    }

    void Reset(std::vector<Float> && x, std::vector<Float> && y)
    {
        GENERIC_ASSERT(x.size() == y.size());
        std::swap(m_x, x);
        std::swap(m_y, y);
        m_size = m_x.size();
        GENERIC_ASSERT(m_size >= 2);

        m_h.resize(m_size - 1);
        m_delta.resize(m_size - 1);
        // compute h and secant slopes
        for (size_t i = 0; i + 1 < m_size; ++i) {
            auto hi = m_x[i+1] - m_x[i];
            GENERIC_ASSERT(hi > 0);
            m_h[i] = hi;
            m_delta[i] = (m_y[i+1] - m_y[i]) / hi;
        }
        ComputeNodeSlopes(); // fills m
        PrecomputeCoeffs(); // fills a,b,c,d per interval

        // cache data pointers for faster hot-path access in Evaluate / LocateInterval
        m_xptr = m_x.data();
        m_aptr = m_a.data();
        m_bptr = m_b.data();
        m_cptr = m_c.data();
        m_dptr = m_d.data();
    }

    inline Float operator() (Float xq) const noexcept
    {
        return Evaluate(xq);
    }

    // Single evaluation: O(log n) search
    inline Float Evaluate(Float xq) const noexcept
    {
        auto idx = LocateInterval(xq);
        return EvalInterval(idx, xq);
    }

    // Batch for random queries: each query does binary search -> O(m log n)
    // x and y must have same length
    void EvaluateBatchRandom(const std::vector<Float>& x, std::vector<Float>& y) const
    {
        y.resize(x.size());
        for (size_t k = 0; k < x.size(); ++k)
            y[k] = Evaluate(x[k]); // hot path uses cached pointers
    }

    // Batch for monotonic-increasing queries: single linear scan -> O(n + m)
    // x must be non-decreasing (monotonic increasing). y must be sized to x.
    void EvaluateBatchMonotonic(const std::vector<Float> & x, std::vector<Float> & y) const
    {
        size_t seg = 0; // current interval index (in [0, m_size-2])
        y.resize(x.size());
        for (size_t k = 0; k < x.size(); ++k) {
            Float xq = x[k];
            // advance seg while xq > m_x[seg+1] and not last interval
            while ((seg + 1) < (m_size - 1) && xq > m_xptr[seg + 1]) ++seg;
            y[k] = EvalInterval(seg, xq);
        }
    }

private:
    // locate interval index i such that x in [x[i], x[i+1]]
    inline size_t LocateInterval(Float xq) const noexcept
    {
        // clamp outside domain to endpoints
        if (xq <= m_x.front()) return 0;
        if (xq >= m_x.back()) return m_size - 2;
        // use pointer-based upper_bound to find first x > xq, then -1 gives index
        const Float * begin = m_xptr;
        const Float * end = m_xptr + m_size;
        auto it = std::upper_bound(begin, end, xq);
        size_t pos = static_cast<size_t>(it - begin);
         size_t idx = (pos == 0) ? 0 : pos - 1;
         // ensure idx in [0, n-2]
         if (idx >= m_size - 1) idx = m_size - 2;
         return idx;
     }

     inline Float EvalInterval(size_t idx, Float xq) const noexcept
     {
         // s = xq - x[idx]
        Float s = xq - m_xptr[idx];
        // Horner for cubic a*s^3 + b*s^2 + c*s + d, using cached pointers
        return ((m_aptr[idx] * s + m_bptr[idx]) * s + m_cptr[idx]) * s + m_dptr[idx];
     }

     void ComputeNodeSlopes()
     {
        m_m.assign(m_size, 0.0);
        if (m_size == 2) {
            m_m[0] = m_m[1] = m_delta[0];
            return;
        }
        // left endpoint: weighted 3-point estimate with monotonicity clamps
        {
            Float h0 = m_h[0];
            Float h1 = m_h[1];
            Float d0 = m_delta[0];
            Float d1 = m_delta[1];
            Float m0 = ((2 * h0 + h1) * d0 - h0 * d1) / (h0 + h1);
            if ((m0 > 0) != (d0 > 0)) m0 = 0;
            if (std::fabs(m0) > 3 * std::fabs(d0)) m0 = 3 * d0;
            m_m[0] = m0;
        }
        // right endpoint: symmetric
        {
            size_t i = m_size - 1;
            Float hnm2 = m_h[m_size - 2];     // last h
            Float hnm3 = m_h[m_size - 3];     // prev h
            Float dnm2 = m_delta[m_size - 2];
            Float dnm3 = m_delta[m_size - 3];
            Float mn = ((2 * hnm2 + hnm3) * dnm2 - hnm2 * dnm3) / (hnm2 + hnm3);
            if ((mn > 0) != (dnm2 > 0)) mn = 0;
            if (std::fabs(mn) > 3 * std::fabs(dnm2)) mn = 3 * dnm2;
            m_m[i] = mn;
        }
 
        // interior slopes using Fritsch-Carlson (weighted harmonic mean) to preserve monotonicity
        for (size_t i = 1; i + 1 < m_size; ++i) {
            Float dl = m_delta[i-1];
            Float dr = m_delta[i];
            if (dl == 0.0 or dr == 0.0 or (dl > 0.0) != (dr > 0.0))
                m_m[i] = 0.0;
            else {
                Float hl = m_h[i-1];
                Float hr = m_h[i];
                Float w1 = 2 * hr + hl;
                Float w2 = hr + 2 * hl;
                // harmonic mean weighted:
                m_m[i] = (w1 + w2) / (w1 / dl + w2 / dr);
            }
        }
    }

    void PrecomputeCoeffs()
    {
        auto intervals = m_size - 1;
        m_a.resize(intervals);
        m_b.resize(intervals);
        m_c.resize(intervals);
        m_d.resize(intervals);
        for (size_t i = 0; i < intervals; ++i) {
            auto hi = m_h[i];
            auto yi = m_y[i];
            auto yi1 = m_y[i+1];
            auto mi = m_m[i];
            auto mi1 = m_m[i+1];
            m_d[i] = yi;
            m_c[i] = mi;
            // coefficients derived from Hermite basis:
            // P(s) = yi*(1 - 3(t^2) + 2(t^3)) + yi1*(3(t^2)-2(t^3)) + mi*h*(t^3-2t^2+t) + mi1*h*(t^3 - t^2)
            // with t = s/h.
            // After algebra:
            m_a[i] = (2 * yi - 2 * yi1 + mi * hi + mi1 * hi) / (hi * hi * hi);
            m_b[i] = (-3 * yi + 3 * yi1 - 2 * mi * hi - mi1 * hi) / (hi * hi);
            // c[i] and d[i] already set
        }
    }

    size_t m_size{0};
    // cached raw pointers to contiguous storage for hot path
    const Float *m_xptr{nullptr};
    const Float *m_aptr{nullptr};
    const Float *m_bptr{nullptr};
    const Float *m_cptr{nullptr};
    const Float *m_dptr{nullptr};
    std::vector<Float> m_x, m_y;
    std::vector<Float> m_h, m_delta, m_m; // node slopes
    // per-interval cubic: p(s) = a*s^3 + b*s^2 + c*s + d, s = x - x_i
    std::vector<Float> m_a, m_b, m_c, m_d;
};

} // namespace generic::math
