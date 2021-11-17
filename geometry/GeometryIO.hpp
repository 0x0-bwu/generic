#ifndef GENERIC_GEOMETRY_GEOMETRYIO_HPP
#define GENERIC_GEOMETRY_GEOMETRYIO_HPP
#include <boost/geometry/io/wkt/write.hpp>
#include <boost/geometry/io/wkt/read.hpp>
#include "BoostGeometryRegister.hpp"
#include "GeometryTraits.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

#ifdef BOOST_GIL_IO_PNG_SUPPORT
#include "Rasterization.hpp"
#include "Utility.hpp"
#include "tools/Color.hpp"
#include <boost/gil/extension/io/png.hpp>
#include <boost/gil.hpp>
#include <png.h>
#endif

namespace {

using namespace generic::geometry;

template <typename num_type>
inline std::ostream & operator<< (std::ostream & os, const Point2D<num_type> & p)
{
    return os << boost::geometry::wkt(p);
}

template <typename num_type>
inline std::istream & operator >> (std::istream & is, Point2D<num_type> & p)
{
    std::string s;
    std::getline(is, s);
    try { boost::geometry::read_wkt(s, p); }
    catch (...) { is.setstate(std::ios::failbit); }
    return is;
}

template <typename num_type>
inline std::ostream & operator<< (std::ostream & os, const Point3D<num_type> & p)
{
    char sp(32);
    return os << "POINT(" << p[0] << sp << p[1] << sp << p[2] << ")";
}

template <typename num_type>
inline std::ostream & operator<< (std::ostream & os, const Segment2D<num_type> & s)
{
    char sp(32);
    os << "SEGMENT(";
    os << s[0][0] << sp << s[0][1] << ",";
    os << s[1][0] << sp << s[1][1] << ")";
    return os;
}

template <typename num_type>
inline std::ostream & operator<< (std::ostream & os, const Segment3D<num_type> & s)
{   
    char sp(32);
    os << "SEGMENT(";
    os << s[0][0] << sp << s[0][1] << sp << s[0][2] << ",";
    os << s[1][0] << sp << s[1][1] << sp << s[1][2] << ")";
    return os;
}

template <typename num_type>
inline std::ostream & operator<< (std::ostream & os, const Triangle2D<num_type> & t)
{
    char sp(32);
    os << "TRIANGLE(";
    os << t[0][0] << sp << t[0][1] << ",";
    os << t[1][0] << sp << t[1][1] << ",";
    os << t[2][0] << sp << t[2][1] << ")";
    return os;
}

template <typename num_type>
inline std::ostream & operator<< (std::ostream & os, const Triangle3D<num_type> & t)
{
    char sp(32);
    os << "TRIANGLE(";
    os << t[0][0] << sp << t[0][1] << sp << t[0][2] << ",";
    os << t[1][0] << sp << t[1][1] << sp << t[1][2] << ",";
    os << t[2][0] << sp << t[2][1] << sp << t[2][2] << ")";
    return os;
}

template <typename num_type>
inline std::ostream & operator<< (std::ostream & os, const Box2D<num_type> & b)
{
    char sp(32);
    os << "BOX(";
    os << b[0][0] << sp << b[0][1] << ",";
    os << b[1][0] << sp << b[1][1] << ")";
    return os;
}

template <typename num_type>
inline std::ostream & operator<< (std::ostream & os, const Box3D<num_type> & b)
{
    char sp(32);
    os << "BOX(";
    os << b[0][0] << sp << b[0][1] << sp << b[0][2] << ",";
    os << b[1][0] << sp << b[1][1] << sp << b[1][2] << ")";
    return os;
}

template <typename num_type>
inline std::ostream & operator << (std::ostream & os, const Polygon2D<num_type> & p)
{
    return os << boost::geometry::wkt(p);
}

template <typename num_type>
inline std::istream & operator >> (std::istream & is, Polygon2D<num_type> & p)
{
    std::string s;
    std::getline(is, s);
    try { boost::geometry::read_wkt(s, p); }
    catch (...) { is.setstate(std::ios::failbit); }
    if(p.Size() && (p.ConstFront() == p.ConstBack())) p.PopBack();
    return is;
}

template <typename num_type>
inline std::ostream & operator << (std::ostream & os, const PolygonWithHoles2D<num_type> & pwh)
{
    return os << boost::geometry::wkt(pwh);
}

template <typename num_type>
inline std::istream & operator >> (std::istream & is, PolygonWithHoles2D<num_type> & pwh)
{   
    using point_t = boost::geometry::model::point<num_type, 2, boost::geometry::cs::cartesian>;
    using polygon_t = boost::geometry::model::polygon<point_t>;
    polygon_t p;
    std::string s;
    std::getline(is, s);
    try { boost::geometry::read_wkt(s, p); }
    catch (...) { is.setstate(std::ios::failbit); }
    pwh = toPolygonWithHoles2D(p);
    return is;
}

template <typename geometry_t, typename std::enable_if<traits::is_geometry_t<geometry_t>::value, bool>::type = true>
inline std::string toString(const geometry_t & geometry)
{
    std::stringstream ss;
    ss << geometry;
    return ss.str();
}

}

namespace generic  {
namespace geometry {

class GeometryIO
{
public:
    template <typename geometry, typename iterator,
              typename std::enable_if<std::is_same<geometry,
              typename std::iterator_traits<iterator>::value_type>::value, bool>::type = true>
    static bool Write(const std::string & file, iterator begin, iterator end)
    {
        std::ofstream out(file);
        if(!out.is_open()) return false;
        for(auto iter = begin; iter != end; ++iter){
            typename std::iterator_traits<iterator>::reference r = *iter;
            if(iter != begin) out << std::endl;
            out << r;
        }
        out.close();
        return true;
    }

    template <typename geometry, typename iterator>
    static bool Read(const std::string & file, iterator result, std::string * err = nullptr)
    {
        std::ifstream in(file);
        if(!in.is_open()) {
            if(err) *err = "Error: fail to open: " + file;
            return false;
        }

        size_t line = 0;
        while(!in.eof()){
            line++;
            if(in.peek() == '#'){
                std::string tmp;
                std::getline(in, tmp);
                continue;
            }
            geometry g; in >> g;
            if(in.fail()){
                if(err) *err = "Error: unrecognized syntax in line: " + std::to_string(line);
                return false;
            }
            *result++ = std::move(g);
        }
        in.close();
        return true;
    }

#ifdef BOOST_GIL_IO_PNG_SUPPORT

    template <typename geometry, typename iterator,
              typename std::enable_if<std::is_same<geometry,
              typename std::iterator_traits<iterator>::value_type>::value, bool>::type = true>
    static void WriteImgView(boost::gil::rgb8_image_t::view_t & view, iterator begin, iterator end,
                             const typename geometry::coor_t stride,
                             const Point2D<typename geometry::coor_t> & ref, int rgb)
    {
        using namespace boost::gil;
        using coor_t = typename geometry::coor_t;

        auto width = view.width();
        auto height = view.height();

        int r, g, b;
        generic::color::RGBFromInt(rgb, r, g, b);
        auto color = rgb8_pixel_t(r, g, b);
        std::vector<std::array<int, 2> > grids;
        for(auto iter = begin;iter != end; ++iter){
            grids.clear();
            Rasterization::Rasterize(*iter, stride, grids, ref);
            for(const auto & grid : grids){
                if(grid[0] < 0 || grid[0] >= width) continue;
                if(grid[1] < 0 || grid[1] >= height) continue;
                view(grid[0], height - grid[1] - 1) = color;
            }
        }
    }

    template <typename geometry, typename iterator,
              typename std::enable_if<std::is_same<geometry,
              typename std::iterator_traits<iterator>::value_type>::value, bool>::type = true>
    static bool WritePNG(const std::string & filename, iterator begin, iterator end, size_t width = 512,
                         int color = generic::color::black, int bgColor = generic::color::white)
    {
        using namespace boost::gil;
        using coor_t = typename geometry::coor_t;

        auto bbox = Extent(begin, end);
        auto stride = bbox.Length() / coor_t(width);
        size_t height = static_cast<size_t>(double(width) / bbox.Length() * bbox.Width());

        rgb8_image_t img(width, height);
        rgb8_image_t::view_t v = boost::gil::view(img);

        int r, g, b;
        generic::color::RGBFromInt(bgColor, r, g, b);
        fill_pixels(v, rgb8_pixel_t(r, g, b));

        WriteImgView<geometry>(v, begin, end, stride, bbox[0], color);

        write_view(filename, boost::gil::view(img), png_tag());
        return true;
    }

#endif//BOOST_GIL_IO_PNG_SUPPORT
};

}//namespace geometry
}//namespace generic

#endif//GENERIC_GEOMETRY_GEOMETRYIO_HPP
