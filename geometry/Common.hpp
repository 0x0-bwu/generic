/**
 * @file Common.hpp
 * @author bwu
 * @brief Common geometry definition
 * @version 0.1
 * @date 2022-02-22
 */
#pragma once
#include "generic/common/Exception.hpp"

namespace generic  {
///@brief geometry classes and functions
namespace geometry {

    enum class Axis { X = 0, Y = 1, Z = 2 };

    enum class Orientation2D { Horizontal = 0, Vertical = 1 }; 

    enum class PointLineLocation { OnLine = 0, Left, Right = 2};

    enum class PointTriangleLocation { Outside = 0, OnEdge1, OnEdge2, OnEdge3, Inside = 4};

    enum class WindingDirection { Clockwise = 0, CounterClockwise = 1, Unknown = 2 };
} // namespace geometry
} // namespace generic
