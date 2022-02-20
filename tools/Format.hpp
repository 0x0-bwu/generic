/**
 * @file Format.hpp
 * @author bwu
 * @brief String formatting
 * @version 0.1
 * @date 2022-02-22
 */
#ifndef GENERIC_FORMAT_HPP
#define GENERIC_FORMAT_HPP
#include <boost/format.hpp>
#include <iomanip>
namespace generic {
namespace format {

using Format = boost::format;

inline Format & FormatArgs(Format & format) { return format; }

///@brief makes a format object by args
template <typename Arg, typename ... Args>
inline Format & FormatArgs(Format & format, Arg && arg, Args &&... args)
{
    format % arg;
    return FormatArgs(format, args...);
}

///@brief returns formated string by input args
template <typename ... Args>
inline std::string Format2String(const std::string & fmt, Args &&... args)
{
    Format format(fmt);
    format = FormatArgs(format, args...);
    return format.str();
}

class FormatHelper
{
public:
    /**
     * @brief append integral number by specified format
     * @tparam width integral minumum width when convert to string
     * @param dest string object to append integral after
     * @param i input integral
     * @param pad padding the rest of the width with char `pad`, default is space
     */
    template <size_t width>
    static void AppendInt(std::string & dest, int i, char pad = char(32))
    {   
        using namespace boost::io;
        auto fmt = Format("%1%") % group(std::setfill(pad), std::dec, std::setw(width), i);
        dest.append(fmt.str().substr(0, width));
    }
};

}//namespace format
}//namespace generic
#endif//GENERIC_FORMAT_HPP