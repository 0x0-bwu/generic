#ifndef GENERIC_FORMAT_HPP
#define GENERIC_FORMAT_HPP
#include <boost/format.hpp>
#include <iomanip>
namespace generic {
namespace format {

using Format = boost::format;

inline Format & FormatArgs(Format & format) { return format; }

template <typename Arg, typename ... Args>
inline Format & FormatArgs(Format & format, Arg && arg, Args &&... args)
{
    format % arg;
    return FormatArgs(format, args...);
}

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
    template <size_t width>
    static void AppendInt(std::string & dest, int i, char pad = char(32))
    {   
        using namespace boost::io;
        auto fmt = Format("%1%") % group(std::setfill(pad), std::dec, std::setw(width), i);
        dest.append(fmt.str());
    }
};

}//namespace format
}//namespace generic
#endif//GENERIC_FORMAT_HPP