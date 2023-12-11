/**
 * @file GeometryIO.hpp
 * @author bwu
 * @brief I/O functions of geometries
 * @version 0.1
 * @date 2022-02-22
 */
#pragma once
#include <boost/geometry/io/wkt/write.hpp>
#include <boost/geometry/io/wkt/read.hpp>
#include "BoostGeometryRegister.hpp"
#include "GeometryTraits.hpp"
#include <filesystem>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

#if BOOST_GIL_IO_PNG_SUPPORT
#include "Rasterization.hpp"
#include "Utility.hpp"
#include "generic/tools/Color.hpp"
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
    return os << "POINT(" << p[0] << sp << p[1] << sp << p[2] << ')';
}

template <typename num_type>
inline std::ostream & operator<< (std::ostream & os, const Segment2D<num_type> & s)
{
    char sp(32);
    os << "SEGMENT(";
    os << s[0][0] << sp << s[0][1] << ',';
    os << s[1][0] << sp << s[1][1] << ')';
    return os;
}

template <typename num_type>
inline std::ostream & operator<< (std::ostream & os, const Segment3D<num_type> & s)
{   
    char sp(32);
    os << "SEGMENT(";
    os << s[0][0] << sp << s[0][1] << sp << s[0][2] << ',';
    os << s[1][0] << sp << s[1][1] << sp << s[1][2] << ')';
    return os;
}

template <typename num_type>
inline std::ostream & operator<< (std::ostream & os, const Triangle2D<num_type> & t)
{
    char sp(32);
    os << "TRIANGLE(";
    os << t[0][0] << sp << t[0][1] << ',';
    os << t[1][0] << sp << t[1][1] << ',';
    os << t[2][0] << sp << t[2][1] << ')';
    return os;
}

template <typename num_type>
inline std::ostream & operator<< (std::ostream & os, const Triangle3D<num_type> & t)
{
    char sp(32);
    os << "TRIANGLE(";
    os << t[0][0] << sp << t[0][1] << sp << t[0][2] << ',';
    os << t[1][0] << sp << t[1][1] << sp << t[1][2] << ',';
    os << t[2][0] << sp << t[2][1] << sp << t[2][2] << ')';
    return os;
}

template <typename num_type>
inline std::ostream & operator<< (std::ostream & os, const Box2D<num_type> & b)
{
    char sp(32);
    os << "BOX(";
    os << b[0][0] << sp << b[0][1] << ',';
    os << b[1][0] << sp << b[1][1] << ')';
    return os;
}

template <typename num_type>
inline std::ostream & operator<< (std::ostream & os, const Box3D<num_type> & b)
{
    char sp(32);
    os << "BOX(";
    os << b[0][0] << sp << b[0][1] << sp << b[0][2] << ',';
    os << b[1][0] << sp << b[1][1] << sp << b[1][2] << ')';
    return os;
}

template <typename num_type>
inline std::ostream & operator<< (std::ostream & os, const Polyline2D<num_type> & l)
{
    char sp(32);
    os << "LINE(";
    for (size_t i = 0; i < l.size() - 1; ++i) {
        os << l[i][0] << sp << l[i][1] << ',';
    }
    os << l.back()[0] << sp << l.back()[1] << ')';
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

///@brief class with geometrie I/O functions
class GeometryIO
{
public:
    /**
     * @brief writes a collection of geometry to VTK(e Visualization Toolkit) file
     * @param[in] file the output *.vtk filename
     * @param[in] begin iterator to the beginning of the geometry collection
     * @param[in] end iterator to the ending of the geometry collection
     * @return whether write file successfully
     */
    template <typename geometry, typename iterator,
              typename std::enable_if<traits::is_polygon_t<geometry>::value && std::is_same<geometry,
              typename std::iterator_traits<iterator>::value_type>::value, bool>::type = true>
    static bool WriteVTK(std::string_view filename, iterator begin, iterator end, typename geometry::coor_t zRef = 0)
    {
        std::ofstream out(filename.data());
        if (not out.is_open()) return false;

        using Point = typename geometry::point_t;
        using IndexPolygon = std::vector<size_t>;

        std::vector<Point> points;
        std::vector<IndexPolygon> polygons;
        polygons.reserve(std::distance(begin, end));

        auto addPoint = [&points](const Point & p) mutable { points.push_back(p); return points.size() - 1; };
        for(auto iter = begin; iter != end; ++iter){
            typename std::iterator_traits<iterator>::reference r = *iter;
            IndexPolygon polygon(r.Size());
            for(size_t i = 0; i < r.Size(); ++i)
                polygon[i] = addPoint(r[i]);
            polygons.emplace_back(std::move(polygon));
        }

        char sp(32);
        out << "# vtk DataFile Version 2.0" << GENERIC_DEFAULT_EOL;
        out << "Unstructured Grid" << GENERIC_DEFAULT_EOL;
        out << "ASCII" << GENERIC_DEFAULT_EOL;
        out << "DATASET UNSTRUCTURED_GRID" << GENERIC_DEFAULT_EOL;
        out << "POINTS" << sp << points.size() << sp << common::toString<typename geometry::coor_t>() << GENERIC_DEFAULT_EOL;
        for(const auto & point : points) {
            out << point[0] << sp << point[1] << sp << zRef << GENERIC_DEFAULT_EOL;
        }
        out << GENERIC_DEFAULT_EOL;
        out << "CELLS" << sp << polygons.size() << sp << points.size() + polygons.size() << GENERIC_DEFAULT_EOL;
        for(const auto & polygon : polygons) {
            out << polygon.size();
            for(const auto & index : polygon)
                out << sp << index;
            out << GENERIC_DEFAULT_EOL;
        }

        out << "CELL_TYPES" << sp << polygons.size() << GENERIC_DEFAULT_EOL;
        for([[maybe_unused]] const auto & polygon : polygons)
            out << '7' << GENERIC_DEFAULT_EOL;

        out.close();
        
        return true;
    }

    /**
     * @brief writes a collection of geometry to WKT(Well-Know Text) file
     * @param[in] file the output *.wkt filename
     * @param[in] begin iterator to the beginning of the geometry collection
     * @param[in] end iterator to the ending of the geometry collection
     * @return whether write file successfully
     */
    template <typename geometry, typename iterator,
              typename std::enable_if<std::is_same<geometry,
              typename std::iterator_traits<iterator>::value_type>::value, bool>::type = true>
    static bool WriteWKT(std::string_view filename, iterator begin, iterator end)
    {
        std::ofstream out(filename);
        if (not out.is_open()) return false;
        for (auto iter = begin; iter != end; ++iter){
            typename std::iterator_traits<iterator>::reference r = *iter;
            if(iter != begin) out << GENERIC_DEFAULT_EOL;
            out << r;
        }
        out.close();
        return true;
    }

    /**
     * @brief reads a collection of geometry from WKT(Well-Know Text) file
     * @param[in] file the input *.wkt filename
     * @param[in] result back insert iterator of the geometry collection
     * @param[out] err error message if failed to read
     * @return whether read file successfully 
     */
    template <typename geometry, typename iterator>
    static bool ReadWKT(std::string_view filename, iterator result, std::string * err = nullptr)
    {
        std::ifstream in(filename);
        if(!in.is_open()) {
            if(err) *err = "Error: fail to open: " + std::string(filename);
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

#if BOOST_GIL_IO_PNG_SUPPORT

    /**
     * @brief draws a collection of geometry to an image file with png format
     * @param[in] filename the output *.png file name
     * @param[in] begin iterator to the beginning of the geometry collection
     * @param[in] end iterator to the beginning of the geometry collection
     * @param[in] width pixel width of the image
     * @param[in] color geometry outline color
     * @param[in] bgColor back ground color
     * @return whether read file successfully  
     */
    template <typename geometry, typename iterator, std::enable_if_t<std::is_same<geometry,
              typename std::iterator_traits<iterator>::value_type>::value, bool> = true>
    static bool WritePNG(std::string_view filename, iterator begin, iterator end, size_t width = 512,
                         int color = generic::color::black, int bgColor = generic::color::white)
    {
        using namespace boost::gil;
        using coor_t = typename geometry::coor_t;

        auto bbox = Extent(begin, end);
        if(!bbox.isValid()) return false;

        //expand the canvas by 10%
        auto expX = bbox.Length() * 0.1;
        auto expY = bbox.Width() * 0.1;
        bbox[0][0] -= expX; bbox[0][1] -= expY;
        bbox[1][0] += expX; bbox[1][1] += expY;

        auto stride = bbox.Length() / coor_t(width);
        size_t height = static_cast<size_t>(double(width) / bbox.Length() * bbox.Width());
        if(width == 0 || height == 0) return false;

        rgb8_image_t img(width, height);
        rgb8_image_t::view_t v = boost::gil::view(img);

        int r, g, b;
        generic::color::RGBFromInt(bgColor, r, g, b);
        fill_pixels(v, rgb8_pixel_t(r, g, b));

        WriteImgView<geometry>(v, begin, end, Vector2D<coor_t>(stride, stride), bbox[0], color);

        auto dirPath = std::filesystem::path(filename).parent_path();
        std::filesystem::create_directories(dirPath);
        if (not std::filesystem::exists(dirPath)) return false;

        write_view(filename.data(), boost::gil::view(img), png_tag());
        return true;
    }

private:
    template <typename geometry, typename iterator,
              typename std::enable_if<std::is_same<geometry,
              typename std::iterator_traits<iterator>::value_type>::value, bool>::type = true>
    static void WriteImgView(boost::gil::rgb8_image_t::view_t & view, iterator begin, iterator end,
                             const Vector2D<typename geometry::coor_t> & stride,
                             const Point2D<typename geometry::coor_t> & ref, int rgb)
    {
        using namespace boost::gil;
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

#endif//BOOST_GIL_IO_PNG_SUPPORT
};

}//namespace geometry
}//namespace generic