/**
 * @file StringHelper.hpp
 * @author bwu
 * @brief String related functions 
 * @version 0.1
 * @date 2022-02-22
 */
#pragma once
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
#include <string_view>
#include <string>
#include <regex>

///@brief string handle functions
namespace generic::str {

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

///@brief checks if string `str` contains substring `sub` case sensitive or case insensitive
template <typename C = CaseSensitive>
inline bool Contains(const std::string & str, const std::string & sub)
{
    if constexpr (C::value)
        return boost::find(str, sub);
    else return boost::ifind_first(str, sub); 
}

///@brief trims the space in `str`
inline std::string Trim(const std::string & str)
{
    return boost::trim_copy(str);
}

///@brief overload version of `Split`
inline void Split(const std::string & str, std::vector<std::string> & items, const std::string & seperators = " ")
{
    items.clear();
    using Tokenizer = boost::tokenizer<boost::char_separator<char> >;
    boost::char_separator<char> sep(seperators.c_str());
    Tokenizer tokens(str, sep);
    for(auto iter = tokens.begin(); iter != tokens.end(); ++iter)
        items.push_back(*iter);
}

///@brief split sting `str` by `seperators`
inline std::vector<std::string> Split(const std::string & str, const std::string & seperators = " ")
{
    std::vector<std::string> items;
    Split(str, items, seperators);
    return items;
}

///@brief split sting `str` by `seperator`
inline std::vector<std::string> Split(const std::string & str, char c)
{
    std::vector<std::string> items;
    Split(str, items, std::string(1, c));
    return items;
}

///@brief checks if a string `str` matches a wildcard pattern `pattern`
inline bool WildcardMatch(std::string_view str, std::string_view pattern)
{
    std::regex starReplace("\\*");
    std::regex questionmakrReplace("\\?");
    auto wildcardPattern = std::regex_replace(
        std::regex_replace(pattern.data(), starReplace, ".*"), questionmakrReplace, ".");
    std::regex re(wildcardPattern);
    return std::regex_match(str.begin(), str.end(), re);
}

} // namespace generic::str
