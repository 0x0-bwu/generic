/**
 * @file HashFunction.hpp
 * @author bwu
 * @brief Hush functions of geometries
 * @version 0.1
 * @date 2022-02-22
 */
#pragma once
#include "generic/tools/Hash.hpp"
#include "Geometries.hpp"

namespace std {

template <typename num_type>
struct hash<generic::geometry::Point2D<num_type>>
{
    size_t operator() (const generic::geometry::Point2D<num_type> & p) const
    {
        size_t seed{0};
        generic::hash::HashCombine(seed, p[0], p[1]);
        return seed;
    }
};

template <typename num_type>
struct hash<generic::geometry::Point3D<num_type>>
{
    size_t operator() (const generic::geometry::Point3D<num_type> & p) const
    {
        return generic::hash::HashCombine(0, p[0], p[1], p[2]);
    }
};

template <typename num_type>
struct hash<generic::geometry::Segment2D<num_type>>
{
    size_t operator() (const generic::geometry::Segment2D<num_type> & s) const
    {
        auto hasher = hash<generic::geometry::Point2D<num_type>>();
        return generic::hash::HashCombine(0, hasher(s[0]), hasher(s[1]));
    }
};

template <typename num_type>
struct hash<generic::geometry::Segment3D<num_type>>
{
    size_t operator() (const generic::geometry::Segment3D<num_type> & s) const
    {
        auto hasher = hash<generic::geometry::Point3D<num_type>>();
        return generic::hash::HashCombine(0, hasher(s[0]), hasher(s[1]));
    }
};

template <typename num_type>
struct hash<generic::geometry::Triangle2D<num_type>>
{
    size_t operator() (const generic::geometry::Triangle2D<num_type> & t) const
    {
        auto hasher = hash<generic::geometry::Point2D<num_type>>();
        return generic::hash::HashCombine(0, hasher(t[0]), hasher(t[1]), hasher(t[2]));
    }
};

template <typename num_type>
struct hash<generic::geometry::Triangle3D<num_type>>
{
    size_t operator() (const generic::geometry::Triangle3D<num_type> & t) const
    {
        auto hasher = hash<generic::geometry::Point3D<num_type>>();
        return generic::hash::HashCombine(0, hasher(t[0]), hasher(t[1]), hasher(t[2]));
    }
};

template <typename num_type>
struct hash<generic::geometry::Box2D<num_type>>
{
    size_t operator() (const generic::geometry::Box2D<num_type> & b) const
    {
        auto hasher = hash<generic::geometry::Point2D<num_type>>();
        return generic::hash::HashCombine(0, hasher(b[0]), hasher(b[1]));
    }
};

template <typename num_type>
struct hash<generic::geometry::Box3D<num_type>>
{
    size_t operator() (const generic::geometry::Box3D<num_type> & b) const
    {
        auto hasher = hash<generic::geometry::Point3D<num_type>>();
        return generic::hash::HashCombine(0, hasher(b[0]), hasher(b[1]));
    }
};

template <typename num_type>
struct hash<generic::geometry::Polyline2D<num_type>>
{
    size_t operator() (const generic::geometry::Polyline2D<num_type> & p) const
    {
        auto hasher = hash<generic::geometry::Point2D<num_type>>();
        return generic::hash::OrderedHash(p, hasher);
    }
};

template <typename num_type>
struct hash<generic::geometry::Polygon2D<num_type>>
{
    size_t operator() (const generic::geometry::Polygon2D<num_type> & p) const
    {
        auto hasher = hash<generic::geometry::Point2D<num_type>>();
        return generic::hash::OrderedHash(p.GetPoints(), hasher);
    }
};

template <typename num_type>
struct hash<generic::geometry::PolygonWithHoles2D<num_type>>
{
    size_t operator() (const generic::geometry::PolygonWithHoles2D<num_type> & p) const
    {
        auto hasher = hash<generic::geometry::Polygon2D<num_type>>();
        auto holeHash = generic::hash::OrderedHash(p.holes, hasher);
        return generic::hash::HashCombine(hasher(p.outline), holeHash);
    }
};

} // namespace std