/**
 * @file Color.hpp
 * @author bwu
 * @brief Color related functions
 * @version 0.1
 * @date 2022-02-22
 */
#ifndef GENERIC_COLOR_HPP
#define GENERIC_COLOR_HPP
#include "generic/math/MathUtility.hpp"
namespace generic{
///@brief color related definitions, classes and functions
namespace color {

//a|R|G|B
inline static constexpr int32_t black = 0x00000000;
inline static constexpr int32_t white = 0xFFFFFFFF;
inline static constexpr int32_t red   = 0xFFFF0000;
inline static constexpr int32_t green = 0xFF00FF00;
inline static constexpr int32_t blue  = 0xFF0000FF;

enum class ColorMap { Grey = 0, Jet };

///@brief converts an 8 bit rgb color to an integral number
inline int RGBToInt(int r, int g, int b)
{
    int32_t c = 0xFF;
    c = (c << 8) | r;
    c = (c << 8) | g;
    c = (c << 8) | b;
    return c;
}

///@brief converts an 8 bit rgb color with alpha channel to an integral number
inline int RGBaToInt(int r, int g, int b, int a)
{
    int32_t c = RGBToInt(r, g, b);
    c &= 0xFFFFFF;
    c |= a << 24;
    return c;
}

///@brief gets 8 bit rgb values from an integral number
inline void RGBFromInt(int c, int & r, int & g, int & b)
{
    b = 0xFF & c; c >>= 8;
    g = 0xFF & c; c >>= 8;
    r = 0xFF & c; c >>= 8;
}

///@brief gets 8 bit rgb and alpha values from an integral number
inline void RGBaFromInt(int c, int & r, int & g, int & b, int & a)
{
    RGBFromInt(c, r, g, b);
    a = 0xFF & (c >> 24);
}

///@brief gets random 8 bit rgb color
inline void RandomRGB(int & r, int & g, int & b)
{
    r = 0xFF & math::Random(0, 255);
    g = 0xFF & math::Random(0, 255);
    b = 0xFF & math::Random(0, 255);
}

///@brief mapping a scalar in range [0, 1] to jet colors rgb value
inline void toJet(double scalar, int & r, int & g, int & b)
{
    scalar *= 8;
    if(scalar < 0) { r = 0; g = 0; b = 127; }
    else if(scalar < 1){ r = 0; g = 0; b = int(scalar * 128) + 128; }
    else if(scalar < 3){ r = 0; g = int(scalar * 128) - 128; b = 255; }
    else if(scalar < 5){ r = int(scalar * 128) - 384; g = 255; b = 639 - int(scalar * 128); }
    else if(scalar < 7){ r = 255; g = 895 - int(scalar * 128); b = 0; }
    else if(scalar < 8) { r= 1151 - int(scalar * 128); g = 0; b = 0; }
    else { r = 127; g = 0; b = 0; }
}

/**
 * @brief Mapping scalar value to color by color map
 * @param[in] scalar input scalar in range [0, 1]
 * @param[out] r 8 bit red channal value 
 * @param[out] g 8 bit green channal value
 * @param[out] b 8 bit blue channla value
 * @param[in] colorMap color map used to map scaler
 */
inline void RGBFromScalar(double scalar, int & r, int & g, int & b, ColorMap colorMap = ColorMap::Jet)
{
    switch(colorMap) {
        case ColorMap::Jet : { toJet(scalar, r, g, b); break; }
        default : { r = 0; g = 0; b = 0; break; }
    }
}

}//namespace color
}//namespace generic
#endif//GENERIC_COLOR_HPP