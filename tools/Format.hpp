/**
 * @file Format.hpp
 * @author bwu
 * @brief String formatting
 * @version 0.1
 * @date 2022-02-22
 */
#pragma once
#include "generic/common/Traits.hpp"
#include <boost/format.hpp>
#include <iterator>
#include <iomanip>
#include <sstream>

namespace generic {
///@brief string format
namespace fmt {

using Format = boost::format;

/**
 * @brief Brief description of FormatArgs.
 * @param format
 * @return inline Format &
 */
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
inline std::string Fmt2Str(const std::string & fmt, Args &&... args)
{
    Format format(fmt);
    format = FormatArgs(format, args...);
    return format.str();
}

///@brief returns formated container contents
template <typename T, std::enable_if_t<not std::is_same_v<T, std::string> and common::iterable<T>, bool> = true>
inline std::string Fmt2Str(const T & container, std::string_view splitter)
{
    std::stringstream ss;
    std::copy(container.cbegin(), container.cend(), std::ostream_iterator<typename T::value_type>(ss, splitter.data()));
    return ss.str();
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

}//namespace fmt
}//namespace generic
