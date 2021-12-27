#ifndef GENERIC_STR_STRINGHELPER_HPP
#define GENERIC_STR_STRINGHELPER_HPP
#include <string>
namespace generic {
namespace str {

inline bool StartsWith(const std::string & s, const std::string & prefix)
{
    if(prefix.size() > s.size()) return false;
    return s.rfind(prefix, 0) == 0;
}

inline bool EndsWith(const std::string & s, const std::string & suffix)
{
    if (suffix.size() > s.size()) return false;
    return std::equal(suffix.rbegin(), suffix.rend(), s.rbegin());
}

}//str
}//generic
#endif//GENERIC_STR_STRINGHELPER_HPP