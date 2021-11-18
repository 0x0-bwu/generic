#ifndef GENERIC_GEOMETRY_RASTERIZATION_HPP
#define GENERIC_GEOMETRY_RASTERIZATION_HPP
#include "generic/common/Exception.hpp"
#include "Geometries.hpp"
#include <vector>
#include <array>
namespace generic {
namespace geometry {
using common::float_type;
class Rasterization
{
public:
    template <typename num_type>
    static std::array<int, 2> Rasterize(const Point2D<num_type> & p, num_type stride, const Point2D<num_type> & ref = Point2D<num_type>{0, 0});

    template <typename num_type>
    static void Rasterize(const Point2D<num_type> & p, num_type stride, std::vector<std::array<int, 2> > & grids, const Point2D<num_type> & ref = Point2D<num_type>{0, 0});
    
    template <typename num_type>
    static void Rasterize(const Point2D<num_type> & s, const Point2D<num_type> & e, num_type stride, std::vector<std::array<int, 2> > & grids, const Point2D<num_type> & ref = Point2D<num_type>{0, 0});
    
    template <typename num_type>
    static void Rasterize(const Segment2D<num_type> & seg, num_type stride, std::vector<std::array<int, 2> > & grids, const Point2D<num_type> & ref = Point2D<num_type>{0, 0});

    template <typename num_type>
    static void Rasterize(const Point2D<num_type> & p1, const Point2D<num_type> & p2, const Point2D<num_type> & p3, num_type stride, std::vector<std::array<int, 2> > & grids, const Point2D<num_type> & ref = Point2D<num_type>{0, 0});
    
    template <typename num_type>
    static void Rasterize(const Triangle2D<num_type> & tri, num_type stride, std::vector<std::array<int, 2> > & grids, const Point2D<num_type> & ref = Point2D<num_type>{0, 0});
    
    template <typename num_type>
    static void Rasterize(const Box2D<num_type> & box, num_type stride, std::vector<std::array<int, 2> > & grids, const Point2D<num_type> & ref = Point2D<num_type>{0, 0});
    
    template <typename num_type>
    static void Rasterize(const Polygon2D<num_type> & polygon, num_type stride, std::vector<std::array<int, 2> > & grids, const Point2D<num_type> & ref = Point2D<num_type>{0, 0});

    template <typename num_type>
    static void Rasterize(const PolygonWithHoles2D<num_type> & pwh, num_type stride, std::vector<std::array<int, 2> > & grids, const Point2D<num_type> & ref = Point2D<num_type>{0, 0});

    template <typename num_type>
    static Polygon2D<num_type> Rasterize(const Polygon2D<num_type> & polygon, num_type stride, const Point2D<num_type> & ref = Point2D<num_type>{0, 0});
};

template <typename num_type>
inline std::array<int, 2> Rasterization::Rasterize(const Point2D<num_type> & p, num_type stride, const Point2D<num_type> & ref)
{
    GENERIC_ASSERT(stride > 0)
    using float_t = float_type<num_type>;
    float_t invRes = 1.0 / float_t(stride);
    num_type dx = p[0] - ref[0];
    num_type dy = p[1] - ref[1];
    auto ix = static_cast<int>(dx * invRes);
    auto iy = static_cast<int>(dy * invRes);
    if(math::isNegative(dx)) ix--;
    if(math::isNegative(dy)) iy--;
    return std::array<int, 2>{ix, iy};
}

template <typename num_type>
inline void Rasterization::Rasterize(const Point2D<num_type> & p, num_type stride, std::vector<std::array<int, 2> > & grids, const Point2D<num_type> & ref)
{
    GENERIC_ASSERT(stride > 0)
    grids.push_back(Rasterize(p, stride, ref));
}

template <typename num_type>
inline void Rasterization::Rasterize(const Point2D<num_type> & s, const Point2D<num_type> & e, num_type stride, std::vector<std::array<int, 2> > & grids, const Point2D<num_type> & ref)
{
    GENERIC_ASSERT(stride > 0)
    using float_t = float_type<num_type>;
    auto sign = [](int x) { return x > 0 ? 1 : (x < 0 ? -1 : 0); };
    auto head = std::array<int, 2>{int((s[0] - ref[0]) / stride), int((s[1] - ref[1]) / stride)};
    auto tail = std::array<int, 2>{int((e[0] - ref[0]) / stride), int((e[1] - ref[1]) / stride)};

    auto dx = std::abs(head[0] - tail[0]);
    auto dy = std::abs(head[1] - tail[1]);
    auto sx = sign(tail[0] - head[0]);
    auto sy = sign(tail[1] - head[1]);

    bool swapped = false;
    if(dy > dx){
        std::swap(dx, dy);
        swapped = true;
    }

    auto d = 2 * dy - dx;
    auto x = head[0];
    auto y = head[1];

    for(auto i = 1; i <= dx; ++i){
        grids.push_back({x, y});
        while(d >= 0){
            if(swapped) x += sx;
            else y += sy;
            d = d - 2 * dx;
        }
        if(swapped) y += sy;
        else x += sx;
        d += 2 * dy;
    }
}

// template <typename num_type>
// inline void Rasterization::Rasterize(const Point2D<num_type> & s, const Point2D<num_type> & e, num_type stride, std::vector<std::array<int, 2> > & grids, const Point2D<num_type> & ref)
// {
//     using float_t = float_type<num_type>;
//     float_t invRes = 1.0 / float_t(stride);    
//     auto dx = e[0] - s[0];
//     auto dy = e[1] - s[1];
//     int sx = math::isPositive(dx) ? 1 : -1;
//     int sy = math::isPositive(dy) ? 1 : -1;

//     dx = std::fabs(dx * invRes);
//     dy = std::fabs(dy * invRes);

//     bool swapped = false;
//     if(dy > dx){
//         std::swap(dx, dy);
//         swapped = true;
//     }

//     auto d = 2.0 * dx - dy;
//     int step = std::ceil(dx);
//     auto [x, y] = Rasterize(s, stride, ref);
//     grids.push_back({x, y});
//     for(int i = 0; i < step; ++i){
//         if(d < 0){
//             if(swapped){
//                 y += sy;
//                 grids.push_back({x, y});
//             }
//             else{
//                 x += sx;
//                 grids.push_back({x, y});
//             }
//             d += 2.0 * dy;
//         }
//         else {
//             x += sx;
//             y += sy;
//             grids.push_back({x, y});
//             d += 2.0 * dy - 2.0 * dx;
//         }
//     }   
// }

template <typename num_type>
inline void Rasterization::Rasterize(const Segment2D<num_type> & seg, num_type stride, std::vector<std::array<int, 2> > & grids, const Point2D<num_type> & ref)
{
    Rasterize<num_type>(seg[0], seg[1], stride, grids, ref);
}

template <typename num_type>
inline void Rasterization::Rasterize(const Point2D<num_type> & p1, const Point2D<num_type> & p2, const Point2D<num_type> & p3, num_type stride, std::vector<std::array<int, 2> > & grids, const Point2D<num_type> & ref)
{
    Rasterize<num_type>(p1, p2, stride, grids, ref);
    Rasterize<num_type>(p2, p3, stride, grids, ref);
    Rasterize<num_type>(p3, p1, stride, grids, ref);
}

template <typename num_type>
inline void Rasterization::Rasterize(const Triangle2D<num_type> & triangle, num_type stride, std::vector<std::array<int, 2> > & grids, const Point2D<num_type> & ref)
{
    Rasterize(triangle[0], triangle[1], triangle[2], stride, grids, ref);
}

template <typename num_type>
inline void Rasterization::Rasterize(const Box2D<num_type> & box, num_type stride, std::vector<std::array<int, 2> > & grids, const Point2D<num_type> & ref)
{
    Rasterize(box[0], Point2D<num_type>(box[1][0], box[0][1]), stride, grids, ref);
    Rasterize(Point2D<num_type>(box[1][0], box[0][1]), box[1], stride, grids, ref);
    Rasterize(box[1], Point2D<num_type>(box[0][0], box[1][1]), stride, grids, ref);
    Rasterize(Point2D<num_type>(box[0][0], box[1][1]), box[0], stride, grids, ref);
}
    
template <typename num_type>
inline void Rasterization::Rasterize(const Polygon2D<num_type> & polygon, num_type stride, std::vector<std::array<int, 2> > & grids, const Point2D<num_type> & ref)
{
    size_t size = polygon.Size();
    for(auto i = 0; i < size; ++i){
        const auto & s = polygon[i];
        const auto & e = polygon[(i + 1) % size];
        std::vector<std::array<int, 2> > tmp;
        Rasterize(s, e, stride, tmp, ref);
        grids.insert(grids.end(), tmp.begin(), tmp.end());
    }
}

template <typename num_type>
inline void Rasterization::Rasterize(const PolygonWithHoles2D<num_type> & pwh, num_type stride, std::vector<std::array<int, 2> > & grids, const Point2D<num_type> & ref)
{
    Rasterize(pwh.outline, stride, grids, ref);
    for(const auto & hole : pwh.holes)
        Rasterize(hole, stride, grids, ref);
}

template <typename num_type>
inline Polygon2D<num_type> Rasterization::Rasterize(const Polygon2D<num_type> & polygon, num_type stride, const Point2D<num_type> & ref)
{
    std::vector<std::array<int, 2> > grids;
    Rasterize<num_type>(polygon, stride, grids, ref);

    auto makePoint = [&stride, &ref](const std::array<int, 2> & grid)
    {
        return Point2D<num_type>(ref[0] + grid[0] * stride, ref[1] + grid[1] * stride);
    };

    auto size = grids.size();
    Polygon2D<num_type> out;
    out.GetPoints().reserve(size);
    for(size_t i = 0; i < size; ++i){
        const auto & curr = grids[i];
        const auto & next = grids[(i + 1) % size];
        out << makePoint(curr);

        auto dx = next[0] - curr[0];
        auto dy = next[1] - curr[1];
        if(dx != 0 && dy != 0 && dx != dy){
            if(std::fabs(dx) > std::fabs(dy))
                out << makePoint({next[0], curr[1]});
            else out << makePoint({curr[0], next[1]});
        }          
    }
    return out;
}

}//namespace geometry
}//namespace generic
#endif//GENERIC_GEOMETRY_RASTERIZATION_HPP