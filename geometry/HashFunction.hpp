/**
 * @file HashFunction.hpp
 * @author bwu
 * @brief Hush functions of geometries
 * @version 0.1
 * @date 2022-02-14
 */
#ifndef GENERIC_GEOMETRY_HASHFUNCTION_HPP
#define GENERIC_GEOMETRY_HASHFUNCTION_HPP
#include "Point.hpp"
#include <boost/functional/hash.hpp>
#include <functional>
namespace generic  {
namespace geometry {

template <typename num_type>
struct PointHash
{
};

template <>
struct PointHash<int32_t>
{
    size_t operator()(const Point2D<int32_t> & point) const noexcept
    {
        size_t seed(0);
        boost::hash_combine(seed, point[0]);
        boost::hash_combine(seed, point[1]);
        return seed;
    }
};

template <>
struct PointHash<int64_t>
{
    size_t operator()(const Point2D<int64_t> & point) const noexcept
    {
        size_t seed(0);
        boost::hash_combine(seed, point[0]);
        boost::hash_combine(seed, point[1]);
        return seed;
    }
};

}//namespace geometry
}//namespace generic
#endif//GENERIC_GEOMETRY_HASHFUNCTION_HPP