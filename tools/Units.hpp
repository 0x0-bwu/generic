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

enum class Time { PICOSECOND = -12, NANOSECOND = -9, MICROSECOND = -6, MILLISECOND = -3, SECOND = 0, MINUTE = 60, HOUR = 3600 };

enum class Length { ANGSTROM = -10, NANOMETER = -9, MICROMETER = -6, MILLIMETER = -3, METER = 0, INCH = 100 };

enum class Temperature { CELSIUS, KELVINS};

enum class Capacitance { F = 0, MF = -3, UF = -6, NF = -9, PF = -12, FF = -15 };

enum class Resistance { KOHM = 1, OHM = 0 };

enum class Current { A = 0, MA = -1, UA = -2 };

enum class Voltage { V = 0, MV = -1, UV = -2 };

enum class Power { W = 0, MW = -1, UW = -2, NW = -3, PW = -4 };

///@brief gets scale from input `unit` to second
inline float Scale2Second(Time unit)
{
    if (auto i = int(unit); i <= 0)
        return std::pow(10, i);
    return float(unit);
}

///@brief gets scale from input `unit` to millisecond
inline float Scale2Millisecond(Time unit)
{
    return Scale2Second(unit) * 1e3f;
}

///@brief gets scale from input `unit` to meter
inline float Scale2Meter(Length unit)
{
    switch (unit) {
        case Length::INCH : return 0.0254;
        default : return std::pow(10, int(unit));
    }
}

///@brief transfer celsius degree to kelvins
inline float Celsius2Kelvins(float t)
{
    return t + 273.15f;
}

///@brief transfer kelvins degree to celsius
inline float Kelvins2Celsius(float t)
{
    return t - 273.15f;
}

///@brief gets scale from input `unit` to SI unit
template <typename Unit, typename Scalar = float>
inline Scalar Scale2SI(Unit unit)
{
    return std::pow(Scalar(10), int(unit));
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
        case Time::PICOSECOND  : return "ps";
        case Time::NANOSECOND  : return "ns";
        case Time::MICROSECOND : return "us";
        case Time::MILLISECOND : return "ms";
        case Time::SECOND : return "sec";
        case Time::MINUTE : return "min";
        case Time::HOUR   : return "hour";
        default : return "";
    }       
}

}//namespace generic
