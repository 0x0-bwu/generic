#ifndef GENERIC_MATH_MATHUTILITY_HPP
#define GENERIC_MATH_MATHUTILITY_HPP
#include "common/Traits.hpp"
#include "Numbers.hpp"
#include <type_traits>
#include <algorithm>
#include <random>
#include <cmath>
namespace generic{
namespace math{
using generic::common::float_type;
static std::default_random_engine randGen;

template <typename num_type>
inline typename std::enable_if<std::is_integral<num_type>::value, num_type>::type
Random(num_type min, num_type max)//random in [min, max]
{
    std::uniform_int_distribution<num_type> u(min, max);
    return u(randGen);
}

template <typename num_type>
inline typename std::enable_if<std::is_floating_point<num_type>::value, num_type>::type
Random(num_type min, num_type max)//random in [min, max]
{
    std::uniform_real_distribution<num_type> u(min, max);
    return u(randGen);
}

template <typename num_type>
inline float_type<num_type> Rad(num_type deg) { return deg * pi / 180; }

template <typename num_type>
inline float_type<num_type> Deg(num_type rad) { return rad * 180 * pi_inv; }

template <typename num_type>
inline bool isNegative(num_type num) { return std::signbit(num); }

template <typename num_type>
inline bool isPositive(num_type num) { return num > 0; }

template <typename num_type, typename std::enable_if<std::is_integral<num_type>::value, bool>::type = true>
inline bool EQ(num_type num1, num_type num2) {
    return num1 == num2;
}

template <typename num_type, typename std::enable_if<std::is_floating_point<num_type>::value, bool>::type = true>
inline bool EQ(num_type num1, num_type num2, num_type tolerance = std::numeric_limits<num_type>::epsilon())
{
    return std::fabs(num1 - num2) <= tolerance;
}

template <typename num_type, typename std::enable_if<std::is_integral<num_type>::value, bool>::type = true>
inline bool NE(num_type num1, num_type num2)
{
    return !EQ(num1, num2);
}

template <typename num_type, typename std::enable_if<std::is_floating_point<num_type>::value, bool>::type = true>
inline bool NE(num_type num1, num_type num2, num_type tolerance = std::numeric_limits<num_type>::epsilon())
{ 
    return !EQ(num1, num2, tolerance);
}

template <typename num_type, typename std::enable_if<std::is_integral<num_type>::value, bool>::type = true>
inline bool GE(num_type num1, num_type num2) {
    return num1 >= num2;
}

template <typename num_type, typename std::enable_if<std::is_floating_point<num_type>::value, bool>::type = true>
inline bool GE(num_type num1, num_type num2, num_type tolerance = std::numeric_limits<num_type>::epsilon())
{
    return num1 > num2 || EQ(num1, num2, tolerance);
}

template <typename num_type, typename std::enable_if<std::is_integral<num_type>::value, bool>::type = true>
inline bool LE(num_type num1, num_type num2) {
    return num1 <= num2;
}

template <typename num_type, typename std::enable_if<std::is_floating_point<num_type>::value, bool>::type = true>
inline bool LE(num_type num1, num_type num2, num_type tolerance = std::numeric_limits<num_type>::epsilon())
{
    return num1 < num2 || EQ(num1, num2, tolerance);
}

template <typename num_type, typename std::enable_if<std::is_integral<num_type>::value, bool>::type = true>
inline bool GT(num_type num1, num_type num2) {
    return num1 > num2;
}

template <typename num_type, typename std::enable_if<std::is_floating_point<num_type>::value, bool>::type = true>
inline bool GT(num_type num1, num_type num2, num_type tolerance = std::numeric_limits<num_type>::epsilon())
{
    return num1 - num2 > tolerance;
}

template <typename num_type, typename std::enable_if<std::is_integral<num_type>::value, bool>::type = true>
inline bool LT(num_type num1, num_type num2) {
    return num1 < num2;
}

template <typename num_type, typename std::enable_if<std::is_floating_point<num_type>::value, bool>::type = true>
inline bool LT(num_type num1, num_type num2, num_type tolerance = std::numeric_limits<num_type>::epsilon())
{
    return num2 - num1 > tolerance;
}

template <typename num_type>
inline float_type<num_type> SafeInv(num_type scalar)
{
    if(EQ(scalar, num_type(0))){
        float_type<num_type> epsilon = std::numeric_limits<float_type<num_type> >::epsilon();
        return float_type<num_type>(1) / std::copysign(epsilon, scalar);
    }
    return float_type<num_type>(1) / scalar;
}

}//namespace math
}//namespace generic
#endif//GENERIC_MATH_MATHUTILITY_HPP