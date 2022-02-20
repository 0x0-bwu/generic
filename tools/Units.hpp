/**
 * @file Units.hpp
 * @author bwu
 * @brief Unit definition and functions
 * @version 0.1
 * @date 2022-02-22
 */
#ifndef GENERIC_UNITS_HPP
#define GENERIC_UNITS_HPP
#include <unistd.h>
#include <string>
#include <cmath>
namespace generic{
///@brief units definition and functions
namespace unit{

enum class Time { Picosecond = -4, Nanosecond = -3, Microsecond = -2, Millisecond = -1, Second = 0, Minute = 60, Hour = 3600 };

enum class Length { Nanometer = -3, Micrometer = -2, Millimeter  = -1, Meter = 0 , Inch = 100 };

///@brief gets scale from input `unit` to second
inline double Scale2Second(Time unit)
{
    if(int(unit) < 0) return std::pow(10, int(unit) * 3);
    else if(int(unit) == 0) return 1.0;
    else return double(unit);
}

///@brief gets scale from input `unit` to millisecond
inline double Scale2Millisecond(Time unit)
{
    return Scale2Second(unit) * 1e3;
}

///@brief gets scale from input `unit` to meter
inline double Scale2Meter(Length unit)
{
    switch (unit) {
        case Length::Inch : return 0.0254;
        default : return std::pow(10, int(unit) * 3);
    }
}

}//namespace unit
}//namespace generic

namespace {
using namespace generic;
inline std::string toString(unit::Time unit)
{
    using namespace unit;
    switch(unit){
        case Time::Picosecond  : return "ps";
        case Time::Nanosecond  : return "ns";
        case Time::Microsecond : return "us";
        case Time::Millisecond : return "ms";
        case Time::Second : return "s";
        case Time::Minute : return "min";
        case Time::Hour   : return "h";
        default : return "";
    }       
}
}

#endif//GENERIC_UNITS_HPP