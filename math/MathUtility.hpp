/**
 * @file MathUtility.hpp
 * @author bwu
 * @brief Utility functions for math
 * @version 0.1
 * @date 2022-02-22
 */
#ifndef GENERIC_MATH_MATHUTILITY_HPP
#define GENERIC_MATH_MATHUTILITY_HPP
#include "generic/common/Traits.hpp"
#include "Numbers.hpp"
#include <type_traits>
#include <algorithm>
#include <random>
#include <cmath>
namespace generic{
namespace math{
using generic::common::float_type;
static std::default_random_engine randGen;

///@brief generates an integral number randomly in range [min, max]
template <typename num_type>
inline typename std::enable_if<std::is_integral<num_type>::value, num_type>::type
Random(num_type min, num_type max)//random in [min, max]
{
    std::uniform_int_distribution<num_type> u(min, max);
    return u(randGen);
}

///@brief generates a floating point number randomly in range [min, max]
template <typename num_type>
inline typename std::enable_if<std::is_floating_point<num_type>::value, num_type>::type
Random(num_type min, num_type max)//random in [min, max]
{
    std::uniform_real_distribution<num_type> u(min, max);
    return u(randGen);
}

///@brief converts angle in degree to radians
template <typename num_type>
inline float_type<num_type> Rad(num_type deg) { return deg * pi / 180; }

///@brief converts angle in radians to degree
template <typename num_type>
inline float_type<num_type> Deg(num_type rad) { return rad * 180 * pi_inv; }

///@brief checks if input number is negative
template <typename num_type>
inline bool isNegative(num_type num) { return std::signbit(num); }

///@brief checks if input number is posivite
template <typename num_type>
inline bool isPositive(num_type num) { return num > 0; }

///@brief checks if two integral numbers equal
template <typename num_type, typename std::enable_if<std::is_integral<num_type>::value, bool>::type = true>
inline bool EQ(num_type num1, num_type num2) {
    return num1 == num2;
}

///@brief checks if two floating point numbers equal in given tolerance
template <typename num_type, typename std::enable_if<std::is_floating_point<num_type>::value, bool>::type = true>
inline bool EQ(num_type num1, num_type num2, num_type tolerance = std::numeric_limits<num_type>::epsilon())
{
    return std::fabs(num1 - num2) <= tolerance;
}

///@brief checks if two integral numbers not equal
template <typename num_type, typename std::enable_if<std::is_integral<num_type>::value, bool>::type = true>
inline bool NE(num_type num1, num_type num2)
{
    return !EQ(num1, num2);
}

///@brief checks if two floating point numbers not equal in given tolerance
template <typename num_type, typename std::enable_if<std::is_floating_point<num_type>::value, bool>::type = true>
inline bool NE(num_type num1, num_type num2, num_type tolerance = std::numeric_limits<num_type>::epsilon())
{ 
    return !EQ(num1, num2, tolerance);
}

///@brief checks if integral number `num1` is greater than or equal to integral number `num2`
template <typename num_type, typename std::enable_if<std::is_integral<num_type>::value, bool>::type = true>
inline bool GE(num_type num1, num_type num2) {
    return num1 >= num2;
}

///@brief checks if floating point number `num1` is greater than or equal to floaing point number `num2` in given tolerance
template <typename num_type, typename std::enable_if<std::is_floating_point<num_type>::value, bool>::type = true>
inline bool GE(num_type num1, num_type num2, num_type tolerance = std::numeric_limits<num_type>::epsilon())
{
    return num1 > num2 || EQ(num1, num2, tolerance);
}

///@brief checks if integral number `num1` is less than or equal to integral number `num2`
template <typename num_type, typename std::enable_if<std::is_integral<num_type>::value, bool>::type = true>
inline bool LE(num_type num1, num_type num2) {
    return num1 <= num2;
}

///@brief checks if floating point number `num1` is less than or equal to floaing point number `num2` in given tolerance
template <typename num_type, typename std::enable_if<std::is_floating_point<num_type>::value, bool>::type = true>
inline bool LE(num_type num1, num_type num2, num_type tolerance = std::numeric_limits<num_type>::epsilon())
{
    return num1 < num2 || EQ(num1, num2, tolerance);
}

///@brief checks if integral number `num1` is greater than integral number `num2`
template <typename num_type, typename std::enable_if<std::is_integral<num_type>::value, bool>::type = true>
inline bool GT(num_type num1, num_type num2) {
    return num1 > num2;
}

///@brief checks if floating point number `num1` is greater than floaing point number `num2` in given tolerance
template <typename num_type, typename std::enable_if<std::is_floating_point<num_type>::value, bool>::type = true>
inline bool GT(num_type num1, num_type num2, num_type tolerance = std::numeric_limits<num_type>::epsilon())
{
    return num1 - num2 > tolerance;
}

///@brief checks if integral number `num1` is less than integral number `num2`
template <typename num_type, typename std::enable_if<std::is_integral<num_type>::value, bool>::type = true>
inline bool LT(num_type num1, num_type num2) {
    return num1 < num2;
}

///@brief checks if floating point number `num1` is less than floaing point number `num2` in given tolerance
template <typename num_type, typename std::enable_if<std::is_floating_point<num_type>::value, bool>::type = true>
inline bool LT(num_type num1, num_type num2, num_type tolerance = std::numeric_limits<num_type>::epsilon())
{
    return num2 - num1 > tolerance;
}

struct OpenInterval {};
struct ClosedInterval {};
struct LeftOpenRightClosed {};
struct LeftClosedRightOpen {};

namespace detail {

template <typename num_type, typename std::enable_if<std::is_integral<num_type>::value, bool>::type = true>
inline bool Within(num_type num, num_type min, num_type max, OpenInterval)
{
    return GT<num_type>(num, min) && LT<num_type>(num, max);
}

template <typename num_type, typename std::enable_if<std::is_floating_point<num_type>::value, bool>::type = true>
inline bool Within(num_type num, num_type min, num_type max, num_type tolerance, OpenInterval)
{
    return GT<num_type>(num, min, tolerance) && LT<num_type>(num, max, tolerance);
}

template <typename num_type, typename std::enable_if<std::is_integral<num_type>::value, bool>::type = true>
inline bool Within(num_type num, num_type min, num_type max, ClosedInterval)
{
    return GE<num_type>(num, min) && LE<num_type>(num, max);
}

template <typename num_type, typename std::enable_if<std::is_floating_point<num_type>::value, bool>::type = true>
inline bool Within(num_type num, num_type min, num_type max, num_type tolerance, ClosedInterval)
{
    return GE<num_type>(num, min, tolerance) && LE<num_type>(num, max, tolerance);
}

template <typename num_type, typename std::enable_if<std::is_integral<num_type>::value, bool>::type = true>
inline bool Within(num_type num, num_type min, num_type max, LeftOpenRightClosed)
{
    return GT<num_type>(num, min) && LE<num_type>(num, max);
}

template <typename num_type, typename std::enable_if<std::is_floating_point<num_type>::value, bool>::type = true>
inline bool Within(num_type num, num_type min, num_type max, num_type tolerance, LeftOpenRightClosed)
{
    return GT<num_type>(num, min, tolerance) && LE<num_type>(num, max, tolerance);
}

template <typename num_type, typename std::enable_if<std::is_integral<num_type>::value, bool>::type = true>
inline bool Within(num_type num, num_type min, num_type max, LeftClosedRightOpen)
{
    return GE<num_type>(num, min) && LT<num_type>(num, max);
}

template <typename num_type, typename std::enable_if<std::is_floating_point<num_type>::value, bool>::type = true>
inline bool Within(num_type num, num_type min, num_type max, num_type tolerance, LeftClosedRightOpen)
{
    return GE<num_type>(num, min, tolerance) && LT<num_type>(num, max, tolerance);
}

}//namespace detail

template <typename interval_type, typename num_type, 
          typename std::enable_if<std::is_integral<num_type>::value, bool>::type = true>
inline bool Within(num_type num, num_type min, num_type max)
{
    return detail::Within<num_type>(num, min, max, interval_type{});
}

template <typename interval_type, typename num_type,
          typename std::enable_if<std::is_floating_point<num_type>::value, bool>::type = true>
inline bool Within(num_type num, num_type min, num_type max, num_type tolerance = std::numeric_limits<num_type>::epsilon())
{
    return detail::Within<num_type>(num, min, max, tolerance, interval_type{});
}

///@brief returns inverse of a scalar, a huge number but not INF if input is closed or equal to zero
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