/**
 * @file StringHelper.hpp
 * @author bwu
 * @brief String related functions 
 * @version 0.1
 * @date 2022-02-22
 */
#ifndef GENERIC_STR_STRINGHELPER_HPP
#define GENERIC_STR_STRINGHELPER_HPP
#include <boost/algorithm/string.hpp>
#include <string>
namespace generic {
///@brief string handle functions
namespace str {

using CaseSensitive = std::true_type;
using CaseInsensitive = std::false_type;

///@brief checks if a string `s` starts with `prefix` case sensitive or case insensitive
template <typename C = CaseSensitive> 
inline bool StartsWith(const std::string & s, const std::string & prefix)
{
    if constexpr (C::value)
        return boost::starts_with(s, prefix);
    else return boost::istarts_with(s, prefix);
}

///@brief checks if a string `s` ends with `suffix` case sensitive or case insensitive
template <typename C = CaseSensitive> 
inline bool EndsWith(const std::string & s, const std::string & suffix)
{
    if constexpr (C::value)
        return boost::ends_with(s, suffix);
    else return boost::iends_with(s, suffix);
}

///@brief checks if string `s1` equals `s2` case sensitive or case insensitive
template <typename C = CaseSensitive>
inline bool Equals(const std::string & s1, const std::string & s2)
{
    if constexpr (C::value)
        return boost::equals(s1, s2);
    else return boost::iequals(s1, s2);
}

}//str
}//generic
#endif//GENERIC_STR_STRINGHELPER_HPP