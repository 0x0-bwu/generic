#ifndef GENERIC_STR_STRINGHELPER_HPP
#define GENERIC_STR_STRINGHELPER_HPP
#include <boost/algorithm/string.hpp>
#include <string>
namespace generic {
namespace str {

using CaseSensitive = std::true_type;
using CaseInsensitive = std::false_type;

template <typename C = CaseSensitive> 
inline bool StartsWith(const std::string & s, const std::string & prefix)
{
    if constexpr (C::value)
        return boost::starts_with(s, prefix);
    else return boost::istarts_with(s, prefix);
}

template <typename C = CaseSensitive> 
inline bool EndsWith(const std::string & s, const std::string & suffix)
{
    if constexpr (C::value)
        return boost::ends_with(s, suffix);
    else return boost::iends_with(s, suffix);
}

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