/**
 * @file Rasterization.hpp
 * @author bwu
 * @brief Rasterize geometry outline to pixel indices
 * @version 0.1
 * @date 2022-02-22
 */
#pragma once
#include "generic/common/Exception.hpp"
#include "Geometries.hpp"
#include "Vector.hpp"
#include <vector>
#include <array>
namespace generic {
namespace geometry {
using common::float_type;
class Rasterization
{
public:
    ///@brief gets rasterized pixel index of a point `p`, from reference point `ref` with x, y stride length `stride`
    template <typename num_type>
    static std::array<int, 2> Rasterize(const Point2D<num_type> & p, const Vector2D<num_type> & stride, const Point2D<num_type> & ref = Point2D<num_type>{0, 0});

    ///@brief gets rasterized pixel indices of a point, from reference point `ref` with x, y stride length `stride`
    template <typename num_type>
    static void Rasterize(const Point2D<num_type> & p, const Vector2D<num_type> & stride, std::vector<std::array<int, 2> > & grids, const Point2D<num_type> & ref = Point2D<num_type>{0, 0});
    
    ///@brief gets rasterized pixel indices of segment formed by point `s` to `e`, from reference point `ref` with x, y stride length `stride`
    template <typename num_type>
    static void Rasterize(const Point2D<num_type> & s, const Point2D<num_type> & e, const Vector2D<num_type> & stride, std::vector<std::array<int, 2> > & grids, const Point2D<num_type> & ref = Point2D<num_type>{0, 0});
    
    ///@brief gets rasterized pixel indices of segment `seg`, from reference point `ref` with x, y stride length `stride`
    template <typename num_type>
    static void Rasterize(const Segment2D<num_type> & seg, const Vector2D<num_type> & stride, std::vector<std::array<int, 2> > & grids, const Point2D<num_type> & ref = Point2D<num_type>{0, 0});

    ///@brief gets rasterized pixel indices of triangle formed by points `p1`, `p2` and `p3`, from reference point `ref` with x, y stride length `stride`
    template <typename num_type>
    static void Rasterize(const Point2D<num_type> & p1, const Point2D<num_type> & p2, const Point2D<num_type> & p3, const Vector2D<num_type> & stride, std::vector<std::array<int, 2> > & grids, const Point2D<num_type> & ref = Point2D<num_type>{0, 0});
    
    ///@brief gets rasterized pixel indices of triangle `tri`, from reference point `ref` with x, y stride length `stride`
    template <typename num_type>
    static void Rasterize(const Triangle2D<num_type> & tri, const Vector2D<num_type> & stride, std::vector<std::array<int, 2> > & grids, const Point2D<num_type> & ref = Point2D<num_type>{0, 0});
    
    ///@brief gets rasterized pixel indices of box `box`, from reference point `ref` with x, y stride length `stride`
    template <typename num_type>
    static void Rasterize(const Box2D<num_type> & box, const Vector2D<num_type> & stride, std::vector<std::array<int, 2> > & grids, const Point2D<num_type> & ref = Point2D<num_type>{0, 0});
    
    ///@brief gets rasterized pixel indices of polygon `polygon`, from reference point `ref` with x, y stride length `stride`
    template <typename num_type>
    static void Rasterize(const Polygon2D<num_type> & polygon, const Vector2D<num_type> & stride, std::vector<std::array<int, 2> > & grids, const Point2D<num_type> & ref = Point2D<num_type>{0, 0});

    ///@brief gets rasterized pixel indices of polygon with holes `pwh`, from reference point `ref` with x, y stride length `stride`
    template <typename num_type>
    static void Rasterize(const PolygonWithHoles2D<num_type> & pwh, const Vector2D<num_type> & stride, std::vector<std::array<int, 2> > & grids, const Point2D<num_type> & ref = Point2D<num_type>{0, 0});

    ///@brief gets rasterized polygon from a polygon from reference point `ref` with x, y stride length `stride`, which each point on the rasterization grid
    template <typename num_type>
    static Polygon2D<num_type> Rasterize(const Polygon2D<num_type> & polygon, const Vector2D<num_type> & stride, const Point2D<num_type> & ref = Point2D<num_type>{0, 0});
};

template <typename num_type>
inline std::array<int, 2> Rasterization::Rasterize(const Point2D<num_type> & p, const Vector2D<num_type> & stride, const Point2D<num_type> & ref)
{
    GENERIC_ASSERT(stride[0] > 0 && stride[1] > 0)
    using float_t = float_type<num_type>;
    float_t invResX = 1.0 / float_t(stride[0]);
    float_t invResY = 1.0 / float_t(stride[1]);
    num_type dx = p[0] - ref[0];
    num_type dy = p[1] - ref[1];
    auto ix = static_cast<int>(dx * invResX);
    auto iy = static_cast<int>(dy * invResY);
    if(math::isNegative(dx)) ix--;
    if(math::isNegative(dy)) iy--;
    return std::array<int, 2>{ix, iy};
}

template <typename num_type>
inline void Rasterization::Rasterize(const Point2D<num_type> & p, const Vector2D<num_type> & stride, std::vector<std::array<int, 2> > & grids, const Point2D<num_type> & ref)
{
    GENERIC_ASSERT(stride[0] > 0 && stride[1] > 0)
    grids.push_back(Rasterize(p, stride, ref));
}

template <typename num_type>
inline void Rasterization::Rasterize(const Point2D<num_type> & s, const Point2D<num_type> & e, const Vector2D<num_type> & stride, std::vector<std::array<int, 2> > & grids, const Point2D<num_type> & ref)
{
    GENERIC_ASSERT(stride[0] > 0 && stride[1] > 0)
    auto sign = [](int x) { return x > 0 ? 1 : (x < 0 ? -1 : 0); };
    auto head = std::array<int, 2>{int((s[0] - ref[0]) / stride[0]), int((s[1] - ref[1]) / stride[1])};
    auto tail = std::array<int, 2>{int((e[0] - ref[0]) / stride[0]), int((e[1] - ref[1]) / stride[1])};

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

template <typename num_type>
inline void Rasterization::Rasterize(const Segment2D<num_type> & seg, const Vector2D<num_type> & stride, std::vector<std::array<int, 2> > & grids, const Point2D<num_type> & ref)
{
    Rasterize<num_type>(seg[0], seg[1], stride, grids, ref);
}

template <typename num_type>
inline void Rasterization::Rasterize(const Point2D<num_type> & p1, const Point2D<num_type> & p2, const Point2D<num_type> & p3, const Vector2D<num_type> & stride, std::vector<std::array<int, 2> > & grids, const Point2D<num_type> & ref)
{
    Rasterize<num_type>(p1, p2, stride, grids, ref);
    Rasterize<num_type>(p2, p3, stride, grids, ref);
    Rasterize<num_type>(p3, p1, stride, grids, ref);
}

template <typename num_type>
inline void Rasterization::Rasterize(const Triangle2D<num_type> & triangle, const Vector2D<num_type> & stride, std::vector<std::array<int, 2> > & grids, const Point2D<num_type> & ref)
{
    Rasterize(triangle[0], triangle[1], triangle[2], stride, grids, ref);
}

template <typename num_type>
inline void Rasterization::Rasterize(const Box2D<num_type> & box, const Vector2D<num_type> & stride, std::vector<std::array<int, 2> > & grids, const Point2D<num_type> & ref)
{
    Rasterize(box[0], Point2D<num_type>(box[1][0], box[0][1]), stride, grids, ref);
    Rasterize(Point2D<num_type>(box[1][0], box[0][1]), box[1], stride, grids, ref);
    Rasterize(box[1], Point2D<num_type>(box[0][0], box[1][1]), stride, grids, ref);
    Rasterize(Point2D<num_type>(box[0][0], box[1][1]), box[0], stride, grids, ref);
}
    
template <typename num_type>
inline void Rasterization::Rasterize(const Polygon2D<num_type> & polygon, const Vector2D<num_type> & stride, std::vector<std::array<int, 2> > & grids, const Point2D<num_type> & ref)
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
inline void Rasterization::Rasterize(const PolygonWithHoles2D<num_type> & pwh, const Vector2D<num_type> & stride, std::vector<std::array<int, 2> > & grids, const Point2D<num_type> & ref)
{
    Rasterize(pwh.outline, stride, grids, ref);
    for(const auto & hole : pwh.holes)
        Rasterize(hole, stride, grids, ref);
}

template <typename num_type>
inline Polygon2D<num_type> Rasterization::Rasterize(const Polygon2D<num_type> & polygon, const Vector2D<num_type> & stride, const Point2D<num_type> & ref)
{
    std::vector<std::array<int, 2> > grids;
    Rasterize<num_type>(polygon, stride, grids, ref);

    auto makePoint = [&stride, &ref](const std::array<int, 2> & grid)
    {
        return Point2D<num_type>(ref[0] + grid[0] * stride[0], ref[1] + grid[1] * stride[1]);
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