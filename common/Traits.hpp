#ifndef GENERIC_COMMON_TRAITS_HPP
#define GENERIC_COMMON_TRAITS_HPP
#include <type_traits>
namespace generic {
namespace common {
template <typename num_type>
using float_type = typename std::conditional<
                            std::is_floating_point<num_type>::value, num_type, typename
                            std::conditional<std::is_same<long long, num_type>::value, long double, double>::type
                            >::type;

struct num_floating_tag {};
struct num_integer_tag {};

template <typename num_type>
struct num_traits_float_or_int
{
    using tag = typename std::conditional<
                            std::is_floating_point<num_type>::value, 
                            num_floating_tag, num_integer_tag>::type;
};

template <typename... args>
struct floating_type_check;

template <>
struct floating_type_check<>
{
    using type = std::true_type;
    static constexpr bool value = type{};
};

template <typename num_type, typename... args>
struct floating_type_check<num_type, args...>
{
    using type = typename std::conditional<
                          std::is_floating_point<num_type>::value &&
                          floating_type_check<args...>::value, std::true_type, std::false_type>::type;
    static constexpr bool value = type{};
};

template <typename... args>
struct integral_type_check;

template <>
struct integral_type_check<>
{
    using type = std::true_type;
    static constexpr bool value = type{};
};

template <typename num_type, typename... args>
struct integral_type_check<num_type, args...>
{
    using type = typename std::conditional<
                          std::is_integral<num_type>::value &&
                          integral_type_check<args...>::value, std::true_type, std::false_type>::type;
    static constexpr bool value = type{};
};

}//namespace common
}//namespace generic
#endif//GENERIC_COMMON_TRAITS_HPP