/**
 * @file Trapezoid.hpp
 * @author bwu
 * @brief Model of trapezoid concept
 * @version 0.1
 * @date 2022-02-22
 */
#ifndef GENERIC_GEOMETRY_TRAPEZOID_HPP
#define GENERIC_GEOMETRY_TRAPEZOID_HPP
#include "generic/common/Exception.hpp"
#include "generic/common/Traits.hpp"
#include "Common.hpp"
#include "Point.hpp"
namespace generic {
namespace geometry{
/*

Orientation: Horizental
   corners[1]_____length[1]
            /     \
           /       \
corners[0]/_________\length[0]

Orientation: Vertical

length[1]
|\
| \
|  \
|   |length[0]
|   |corners[0]
|  /
| /
|/
corners[1]

*/
using generic::common::float_type;

template <typename num_type>
class Trapezoid
{
    using float_t = float_type<num_type>;
public:
    const static size_t dim = 2;
    using coor_t = num_type;
    using point_t = Point2D<coor_t>;

    Orientation2D direction;
    num_type length[2];
    Point2D<num_type> corner[2];

    Trapezoid() = default;
    Trapezoid(const Point2D<num_type> & c1, const Point2D<num_type> & c2, num_type len1, num_type len2, Orientation2D d)
     : direction(d)
    {
        length[0] = len1; length[1] = len2;
        corner[0] = c1;   corner[1] = c2;
        Normalize();
    }

    ///@brief gets corner point by index 0-3
    Point2D<num_type> operator[] (size_t index) const;

    ///@brief checks whether this trapezoid is valid
    bool isValid() const;
    ///@brief makes this trapezoid valid
    void Normalize();
};

template <typename num_type>
inline bool Trapezoid<num_type>::isValid() const
{
    if(std::signbit(length[0]) || std::signbit(length[1])) return false;
    if(direction == Orientation2D::Horizontal && math::GT(corner[0][1], corner[1][1])) return false;
    if(direction == Orientation2D::Vertical   && math::GT(corner[1][0], corner[0][0])) return false;
    return true;
}

template <typename num_type>
inline void Trapezoid<num_type>::Normalize()
{
    auto swap = [this] () mutable { std::swap(corner[0], corner[1]); std::swap(length[0], length[1]); };
    if(direction == Orientation2D::Horizontal && math::GT(corner[0][1], corner[1][1])) swap();
    if(direction == Orientation2D::Vertical   && math::GT(corner[1][0], corner[0][0])) swap();
    auto i = (direction == Orientation2D::Horizontal) ? 0 : 1;
    for(auto n = 0; n < 2; ++n){
        if(std::signbit(length[n])){
            corner[n][i] += length[n];
            length[n] = -length[n];
        }
    }
}

template <typename num_type>
inline Point2D<num_type> Trapezoid<num_type>::operator[] (size_t index) const
{
    bool h = (direction == Orientation2D::Horizontal);
    if(0 == index) return corner[0];
    else if(1 == index) return corner[0] + (h ? Point2D<num_type>(length[0], 0) : Point2D<num_type>(0, length[0]));
    else if(2 == index) return corner[1] + (h ? Point2D<num_type>(length[1], 0) : Point2D<num_type>(0, length[1]));
    else if(3 == index) return corner[1];
    else GENERIC_THROW(std::out_of_range("index out out range"))
}

/**
 * @brief tries constructing a trapezoid with four points and given direction
 * @param[in] points input points
 * @param[in] direction trapezoid direction
 * @param[out] res whether construct trapezoid successfully
 * @return constructed trapezoid
 */
template <typename num_type>
inline Trapezoid<num_type> toTrapezoid(const std::array<Point2D<num_type>, 4> & points, Orientation2D direction, bool & res)
{
    res = true;
    std::vector<size_t> indices{0, 1, 2, 3};
    auto i = (direction == Orientation2D::Horizontal) ? 1 : 0;
    auto j = (direction == Orientation2D::Horizontal) ? 0 : 1;
    auto cmpX = [&points](size_t i1, size_t i2) { return points[i1][0] < points[i2][0]; };
    auto cmpY = [&points](size_t i1, size_t i2) { return points[i1][1] > points[i2][1]; };
    if(direction == Orientation2D::Horizontal)
        std::sort(indices.begin(), indices.end(), cmpY);
    else std::sort(indices.begin(), indices.end(), cmpX);
    if(math::NE(points[indices[0]][i], points[indices[1]][i])) res = false;
    if(math::NE(points[indices[2]][i], points[indices[3]][i])) res = false;
    if(!res) return Trapezoid<num_type>();

    return Trapezoid<num_type>{points[indices[0]],
                               points[indices[2]],
                               points[indices[1]][j] - points[indices[0]][j],
                               points[indices[3]][j] - points[indices[2]][j], 
                               direction};
}

/**
 * @brief tries constructing a trapezoid with four points
 * @param[in] points input points
 * @param[out] res whether construct trapezoid successfully
 * @return constructed trapezoid
 */
template <typename num_type>
inline Trapezoid<num_type> toTrapezoid(const std::array<Point2D<num_type>, 4> & points, bool & res)
{
    auto trapezoid = toTrapezoid(points, Orientation2D::Horizontal, res);
    if(!res) trapezoid = toTrapezoid(points, Orientation2D::Vertical, res);
    return trapezoid;
}

}//namespace geometry
}//namespace generic

#endif//GENERIC_GEOMETRY_TRAPEZOID_HPP