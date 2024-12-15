/**
 * @file Units.hpp
 * @author bwu
 * @brief Unit definition and functions
 * @version 0.1
 * @date 2022-02-22
 */
#pragma once
#include <unistd.h>
#include <string>
#include <cmath>
namespace generic{
///@brief units definition and functions
namespace unit{

enum class Time { Picosecond = -4, Nanosecond = -3, Microsecond = -2, Millisecond = -1, Second = 0, Minute = 60, Hour = 3600 };

enum class Length { Nanometer = -3, Micrometer = -2, Millimeter  = -1, Meter = 0 , Inch = 100 };

enum class Temperature { Celsius, Kelvins };

enum class Capacitance { F = 0, MF = -1, UF = -2, NF = -3, PF = -4, FF = -5 };

enum class Resistance { KOHM = 1, OHM = 0 };

enum class Current { A = 0, MA = -1, UA = -2 };

enum class Voltage { V = 0, MV = -1, UV = -2 };

enum class Power { W = 0, MW = -1, UW = -2, NW = -3, PW = -4 };

///@brief gets scale from input `unit` to second
inline double Scale2Second(Time unit)
{
    if(int(unit) <= 0) return std::pow(10, int(unit) * 3);
    return double(unit);
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

///@brief transfer celsius degree to kelvins
inline double Celsius2Kelvins(double t)
{
    return t + 273.15;
}

///@brief transfer kelvins degree to celsius
inline double Kelvins2Celsius(double t)
{
    return t - 273.15;
}

///@brief gets scale from input `unit` to SI unit
template <typename Unit, typename Scalar = float>
inline Scalar Scale2SI(Unit unit)
{
    return std::pow(Scalar(10), (int(unit) * 3));
}

template <typename Scalar = float>
inline Scalar Scale2SI(Time t)
{
    return Scale2Second(t);
}

template <typename Scalar = float>
inline Scalar Scale2SI(Length l)
{
    return Scale2Meter(l);
}

template <typename Scalar = float>
inline Scalar Scalar2SI(Temperature t) = delete;

}//namespace unit

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

}//namespace generic
