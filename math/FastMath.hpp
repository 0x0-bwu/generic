/**
 * @file FastMath.hpp
 * @author bwu
 * @brief fast math operations
 * @version 0.1
 * @date 2024-10-12
 */
#pragma once
#include <cstdint>
#include <cmath>
namespace generic::math {

/** "Fast" sine approximation, valid for x in [-pi, pi], max abs error: 3.89917e-05
 * @note Adapted from Paul Mineiro's FastFloat
 */
inline float FastSinUniCycle(float x)
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

/** "Fast" sine approximation, valid on full x domain
 * @note Adapted from Paul Mineiro's FastFloat
 */
float FastSin(float x)
{
    const int32_t k = (int32_t)(x * 0.15915494309189534f);
    return FastSinUniCycle(((x < 0 ? -0.5f : 0.5f) + k) * 6.283185307179586f - x);
}

} // namespace generic::math