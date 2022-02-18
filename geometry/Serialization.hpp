/**
 * @file Serialization.hpp
 * @author bwu
 * @brief Boost serialization support for geometry classes
 * @version 0.1
 * @date 2022-02-14
 */
#ifndef GENERIC_GEOMETRY_SERIALIZATION_HPP
#define GENERIC_GEOMETRY_SERIALIZATION_HPP
#define GENERIC_GEOMETRY_ARCHIVE_VERSION 1
#include "generic/common/Archive.hpp"
#include "Geometries.hpp"
#include "Common.hpp"
namespace boost {
namespace serialization {

using namespace generic::geometry;
template <typename Archive, typename num_type>
void serialize(Archive & ar, Point2D<num_type> & p, const unsigned int)
{
    ar & boost::serialization::make_nvp("x", p[0]);
    ar & boost::serialization::make_nvp("y", p[1]);
}

template <typename Archive, typename num_type>
void serialize(Archive & ar, Point3D<num_type> & p, const unsigned int)
{
    ar & boost::serialization::make_nvp("x", p[0]);
    ar & boost::serialization::make_nvp("y", p[1]);
    ar & boost::serialization::make_nvp("z", p[2]);
}

template <typename Archive, typename num_type>
void serialize(Archive & ar, Segment2D<num_type> & s, const unsigned int)
{
    ar & boost::serialization::make_nvp("low", s[0]);
    ar & boost::serialization::make_nvp("high", s[1]);
}

template <typename Archive, typename num_type>
void serialize(Archive & ar, Segment3D<num_type> & s, const unsigned int)
{
    ar & boost::serialization::make_nvp("low", s[0]);
    ar & boost::serialization::make_nvp("high", s[1]);
}

template <typename Archive, typename num_type>
void serialize(Archive & ar, Triangle2D<num_type> & t, const unsigned int)
{
    ar & boost::serialization::make_nvp("v1", t[0]);
    ar & boost::serialization::make_nvp("v2", t[1]);
    ar & boost::serialization::make_nvp("v3", t[2]);
    ar & boost::serialization::make_nvp("direction", t.GetWindingDirection());
}

template <typename Archive, typename num_type>
void serialize(Archive & ar, Triangle3D<num_type> & t, const unsigned int)
{
    ar & boost::serialization::make_nvp("v1", t[0]);
    ar & boost::serialization::make_nvp("v2", t[1]);
    ar & boost::serialization::make_nvp("v3", t[2]);
}

template <typename Archive, typename num_type>
void serialize(Archive & ar, Box2D<num_type> & b, const unsigned int)
{
    ar & boost::serialization::make_nvp("ll", b[0]);
    ar & boost::serialization::make_nvp("ur", b[1]);
}

template <typename Archive, typename num_type>
void serialize(Archive & ar, Box3D<num_type> & b, const unsigned int)
{
    ar & boost::serialization::make_nvp("ll", b[0]);
    ar & boost::serialization::make_nvp("ur", b[1]);
}

template <typename Archive, typename num_type>
void serialize(Archive & ar, Polygon2D<num_type> & p, const unsigned int)
{
    ar & boost::serialization::make_nvp("points", p.GetPoints());
}

template <typename Archive, typename num_type>
void serialize(Archive & ar, PolygonWithHoles2D<num_type> & pwh, const unsigned int)
{
    ar & boost::serialization::make_nvp("outline", pwh.outline);
    ar & boost::serialization::make_nvp("holes", pwh.holes);
}

}//namespace serialization
}//namespace boost
#endif//GENERIC_GEOMETRY_SERIALIZATION_HPP