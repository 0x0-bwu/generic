/**
 * @file FastMath.hpp
 * @author bwu
 * @brief fast math operations from Paul Mineiro's FastFloat
 * @version 0.1
 * @date 2024-10-12
 */
#pragma once
#include <cstdint>
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

/**
 * @brief "Fast" x to the power of p approximation
 * @note Warning: Seems to have divergent segments with discontinuities for some base/exponent combinations
 */
inline float FastPow(float x, float p)
{
    return FastPow2(p * FastLog2(x));
}

///@brief "Faster" x to the power of p approximation
inline float FasterPow(float x, float p)
{
    return FasterPow2(p * FasterLog2(x));
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

} // namespace generic::math