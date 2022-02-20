/**
 * @file Numbers.hpp
 * @author bwu
 * @brief Define of const numbers in math
 * @version 0.1
 * @date 2022-02-22
 */
#ifndef GENERIC_MATH_NUMBERS_HPP
#define GENERIC_MATH_NUMBERS_HPP
namespace generic {
namespace math {
    inline static constexpr double pi = 3.14159265358979323846;
    inline static constexpr double pi_2 = 2 * pi;
    inline static constexpr double pi_half = 0.5 * pi;
    inline static constexpr double pi_quarter = 0.25 * pi;
    inline static constexpr double pi_inv = 1.0 / pi;

    inline static constexpr double sqrt_2 = 1.4142135623730951;
}//namespace math
}//namespace generic
#endif//GENERIC_MATH_NUMBERS_HPP