/**
 * @file Utility.ipp
 * @author bwu
 * @brief Implementation of funtions in header 'Utility.hpp'
 * @version 0.1
 * @date 2022-02-22
 */
#pragma once
#include "Utility.hpp"
namespace generic {
namespace geometry {
using generic::common::float_type;

namespace detail {

template <typename num_type, typename = void>
class Point2DIndexMap {};

template <typename num_type>
class Point2DIndexMap<num_type, typename std::enable_if<std::is_integral<num_type>::value>::type>
{
public:
    using Point = Point2D<num_type>;

    void Clear() { m_ptMap.clear(); }
    void Insert(const Point & p, size_t index) { m_ptMap.insert(std::make_pair(p, index)); }
    bool Find(const Point & p, size_t & index) const
    {
        for(int i : {0, -1, 1}) {
            for(int j : {0, -1, 1}) {
                auto iter = m_ptMap.find(Point(p[0] + i, p[1] + j));
                if(iter != m_ptMap.end()) {
                    index = iter->second;
                    return true;
                }
            }
        }
        return false;
    }
private:
    std::unordered_map<Point, size_t> m_ptMap;
};

template <typename num_type>
class Point2DIndexMap<num_type, typename std::enable_if<std::is_floating_point<num_type>::value>::type>
{
public:
    using Point = Point2D<num_type>;
    using IndexPoint = std::pair<Point, size_t>;

    void Clear() { m_rtree.clear(); }
    void Insert(const Point & p, size_t index) { m_rtree.insert(std::make_pair(p, index)); }
    bool Find(const Point & p, size_t & index) const
    {
        std::vector<IndexPoint> query;
        m_rtree.query(boost::geometry::index::nearest(p, 1), std::back_inserter(query));
        if(query.empty()) return false;
        if(query.front().first == p) {
            index = query.front().second;
            return true;
        }
        return false;
    }
private:
    boost::geometry::index::rtree<IndexPoint, boost::geometry::index::linear<16, 4> > m_rtree;
};
}//namespace detail

template<typename vector_t, typename std::enable_if<traits::is_vector_t<vector_t>::value, bool>::type>
inline auto Inverse(const vector_t & vec) -> vector_f<vector_t>
{
    vector_f<vector_t> inv;
    for(size_t i = 0; i < vector_f<vector_t>::dim; ++i)
        inv[i] = coor_f<vector_t>(1) / vec[i];
    return inv;
}

template<typename vector_t, typename std::enable_if<traits::is_vector_t<vector_t>::value, bool>::type>
inline auto SafeInverse(const vector_t & vec) -> vector_f<vector_t>
{
    vector_f<vector_t> inv;
    auto epsilon = std::numeric_limits<float_type<coor_f<vector_t> > >::epsilon();
    for(size_t i = 0; i < vector_f<vector_t>::dim; ++i)
        inv[i] = coor_f<vector_t>(1) / (std::fabs(vec[i]) < epsilon ? std::copysign(epsilon, vec[i]) : vec[i]);
    return inv;
}

template<typename num_type>
inline Arc<num_type> toArc(const Arc3<num_type> & arc3)
{   
    float_type<num_type> r2;
    auto origin = CircumCircle(arc3.start, arc3.mid, arc3.end, r2).template Cast<num_type>();
    auto radian = Angle(arc3.start, origin, arc3.end);
    auto isCCW = Triangle2D<num_type>::isCCW(arc3.start, arc3.mid, arc3.end);
    return Arc<num_type>(origin, arc3.start, isCCW ? radian : -radian);
}

template <typename num_type>
inline Polyline2D<num_type> toPolyline(const Arc<num_type> & arc, size_t div)
{
    GENERIC_ASSERT(div != 0);
    bool ccw = arc.radian > 0;
    auto step = std::fabs(arc.radian) / div;

    Polyline2D<num_type> polyline;
    polyline.reserve(div + 1);

    float_type<num_type> alpha, mag;
    arc.GetStartAlphaMag(alpha, mag);
    polyline.emplace_back(arc.start);
    for(size_t i = 1; i < div; ++i){
        auto theta = ccw ? alpha + step * i : alpha - step * i;
        auto point = Point2D<num_type>(mag * std::cos(theta), mag * std::sin(theta)) + arc.origin;
        polyline.emplace_back(std::move(point));
    }
    polyline.push_back(arc.EndPoint());
    return polyline;
}

template <typename num_type>
inline Polyline2D<num_type> toPolyline(const Arc3<num_type> & arc3, size_t div)
{
    return toPolyline(toArc(arc3), div);
}

template <typename num_type>
inline Polygon2D<num_type> toPolygon(const Triangle2D<num_type> & tri)
{
    Polygon2D<num_type> p;
    if(tri.isCCW()) { p << tri[0] << tri[1] << tri[2];}
    else { p << tri[0] << tri[2] << tri[1]; }
    return p;
}

template <typename num_type>
inline Polygon2D<num_type> toPolygon(const Box2D<num_type> & box)
{
    Polygon2D<num_type> p;
    p << box[0] << Point2D<num_type>(box[1][0], box[0][1]);
    p << box[1] << Point2D<num_type>(box[0][0], box[1][1]);
    return p;
}

template <typename num_type>
inline Polygon2D<num_type> toPolygon(const Ellipse<num_type> & ellipse, size_t div)
{
    GENERIC_ASSERT(div != 0);
    Polygon2D<num_type> p;
    auto step = math::pi * 2 / div;
    for(size_t i = 0; i < div; ++i){
        auto theta = step * i;
        auto point = Point2D<num_type>(ellipse.a * std::cos(theta), ellipse.b * std::sin(theta)) + ellipse.o;
        p << point;
    }
    return p;
}

template <typename num_type>
inline Polygon2D<num_type> toPolygon(const Polyline2D<num_type> & polyline, num_type width, size_t roundDiv)
{
    using float_t = float_type<num_type>;
    GENERIC_ASSERT(polyline.size() > 0);
    
    auto roundOrForward = [] (float_t value) {
        if constexpr (std::is_integral<num_type>::value)
            return std::round(value);
        return value;
    };
    
    float_t theta[2] = {0, 0};
    std::deque<Point2D<num_type> > points;
    auto iter = polyline.begin();
    for(; iter != polyline.end(); ++iter){
        auto prev = iter, next = iter;
        if(prev != polyline.begin()) -- prev;
        ++next;
        if(next == polyline.end()) --next;

        Point2D<num_type> point = *iter;
        Point2D<num_type> neighbors[2] = { *prev, *next };
        if(point == neighbors[0] && point == neighbors[1]) continue;

        std::fill(theta, theta+2, std::numeric_limits<float_t>::max());
        for(size_t i = 0; i < 2; ++i){
            auto dx = neighbors[i][0] - point[0];
            auto dy = neighbors[i][1] - point[1];

            if(math::NE<num_type>(dx, 0) || math::NE<num_type>(dy, 0)){
                float_t sinTheta = dy / std::sqrt(dx * dx + dy * dy);
                float_t cosTheta = dx / std::sqrt(dx * dx + dy * dy);
                theta[i] = math::GE<float_t>(sinTheta, 0) ? std::acos(cosTheta) : math::pi_2 - std::acos(cosTheta);
            }
        }
        if(theta[0] == std::numeric_limits<float_t>::max()) theta[0] = theta[1] + math::pi;
        else if(theta[1] == std::numeric_limits<float_t>::max()) theta[1] = theta[0] + math::pi;
        
        if (iter == prev){
            point[0] = roundOrForward(point[0] + width / 2.0 * std::cos(theta[0]));
            point[1] = roundOrForward(point[1] + width / 2.0 * std::sin(theta[0]));
        }
        if (iter == next){
            point[0] = roundOrForward(point[0] + width / 2.0 * std::cos(theta[1]));
            point[1] = roundOrForward(point[1] + width / 2.0 * std::sin(theta[1]));
        }

        float_t targetTheta = 0.5 * (theta[0] + theta[1]);
        Point2D<num_type> pps[2];//polygon points
        if(points.empty()){
            float_t targetDx = width / 2.0 * std::cos(targetTheta);
            float_t targetDy = width / 2.0 * std::sin(targetTheta);
            pps[0] = Point2D<num_type>(roundOrForward(point[0] + targetDx), roundOrForward(point[1] + targetDy));
            pps[1] = Point2D<num_type>(roundOrForward(point[0] - targetDx), roundOrForward(point[1] - targetDy));
        }
        else {
            for(size_t i = 0; i < 2; ++i){
                auto pp = (i == 0) ? points.front() : points.back();//prev point
                pps[i][0] = roundOrForward((pp[1] * std::cos(theta[0]) * std::cos(targetTheta) -
                                        point[1] * std::cos(theta[0]) * std::cos(targetTheta) - 
                                        pp[0] * std::cos(targetTheta) * std::sin(theta[0]) + 
                                        point[0] * std::cos(theta[0]) * std::sin(targetTheta)) / 
                                        (std::cos(theta[0]) * std::sin(targetTheta) - std::cos(targetTheta) * std::sin(theta[0])));

                pps[i][1] = roundOrForward((pp[1] * std::cos(theta[0]) * std::sin(targetTheta) - 
                                        point[1] * std::cos(targetTheta) * std::sin(theta[0]) - 
                                        pp[0] * std::sin(theta[0]) * std::sin(targetTheta) + 
                                        point[0] * std::sin(theta[0]) * std::sin(targetTheta)) /
                                        (std::cos(theta[0]) * std::sin(targetTheta) - std::cos(targetTheta) * std::sin(theta[0])));
            }
        }
        points.push_front(pps[0]);
        points.push_back(pps[1]);
    }

    Polygon2D<num_type> polygon;
    polygon.Insert(polygon.End(), points.begin(), points.end());

    if(polygon.Size() == 0) {
        Point2D<num_type> half(roundOrForward(0.5 * width), roundOrForward(0.5 * width));
        Box2D<num_type> box(polyline.front() - half, polyline.front() + half);
        polygon = toPolygon(box);
    }

    if (roundDiv > 0) polygon = RoundCorners(polygon, num_type(0.5 * width), roundDiv);
    return polygon;
}

template <typename num_type>
inline Point2D<float_type<num_type>> InscribedCircle(const Point2D<num_type> & s, const Point2D<num_type> & p, const Point2D<num_type> & e, num_type r, Point2D<float_type<num_type>> & fps, Point2D<float_type<num_type>> & fpe)
{
    auto a = 0.5 * InnerAngle(s, p, e);;
    auto l = r / std::tan(a);
    auto spLen = Distance(s, p);
    auto epLen = Distance(e, p);
    for (size_t i = 0; i < 2; ++i) {
        fps[i] = (s[i] - p[i]) * l / spLen + p[i];
        fpe[i] = (e[i] - p[i]) * l / epLen + p[i]; 
    }
    auto scale = 1 / (2 * std::cos(a) * std::cos(a));
    Point2D<float_type<num_type>> o;
    for (size_t i = 0; i < 2; ++i)
        o[i] = scale * (fps[i] + fpe[i] - 2 * p[i]) + p[i];
    return o;
}

template <typename num_type>
inline Polygon2D<num_type> InscribedPolygon(const Circle<num_type> & c, size_t div)
{
    assert(div >= 3);
    Polygon2D<num_type> polygon;
    polygon.GetPoints().reserve(div);
    
    float_t ang = math::pi_2 / div;
    for(size_t i = 0; i < div; ++i)
        polygon << (c.o + Point2D<num_type>(std::cos(ang * i) * c.r, std::sin(ang * i) * c.r));
    
    return polygon;
}

template <typename num_type>
inline Polygon2D<num_type> CircumscribedPolygon(const Circle<num_type> & c, size_t div)
{
    assert(div >= 3);
    float_t ang = math::pi / div;
    return InscribedPolygon(Circle<num_type>(c.o, c.r / std::cos(ang)), div);
}
template <typename num_type>
inline Circle<float_type<num_type> > DiametralCircle(const Point2D<num_type> & p1, const Point2D<num_type> & p2)
{
    Circle<float_type<num_type> > c;
    c.o = (p1 + p2).template Cast<float_type<num_type> >() * 0.5;
    c.r = Distance(p1, p2) * 0.5;
    return c;
}

template <typename num_type>
inline Circle<float_type<num_type> > CircumCircle(const Point2D<num_type> & p1, const Point2D<num_type> & p2, const Point2D<num_type> & p3)
{
    float_type<num_type> r2(0);
    Circle<float_type<num_type> > circle;
    circle.o = CircumCircle(p1, p2, p3, r2);
    circle.r = std::sqrt(r2);
    return circle;
}

template <typename num_type>//return coor of circle
inline Point2D<float_type<num_type> > CircumCircle(const Point2D<num_type> & p1, const Point2D<num_type> & p2, const Point2D<num_type> & p3, float_type<num_type> & radius2)
{   
    using float_t = float_type<num_type>;
    float_t epsilon = std::numeric_limits<float_t>::epsilon();
    float_t a1 = p2[1] - p1[1];
    float_t a2 = p3[1] * p3[1] - p1[1] * p1[1] + p3[0] * p3[0] - p1[0] * p1[0];
    float_t b1 = p3[1] - p1[1];
    float_t b2 = p2[1] * p2[1] - p1[1] * p1[1] + p2[0] * p2[0] - p1[0] * p1[0];
    float_t c = 2.0 * ((p3[0] - p1[0]) * (p2[1] - p1[1]) - (p2[0] - p1[0]) * (p3[1] - p1[1]));
    if(std::fabs(c) < epsilon) c = std::copysign(epsilon, c);
    float_t x = (a1 * a2 - b1 * b2) / c;

    float_t d1 = p2[0] - p1[0];
    float_t d2 = p3[0] * p3[0] - p1[0] * p1[0] + p3[1] * p3[1] - p1[1] * p1[1];
    float_t e1 = p3[0] - p1[0];
    float_t e2 = p2[0] * p2[0] - p1[0] * p1[0] + p2[1] * p2[1] - p1[1] * p1[1];
    float_t f = 2.0 * ((p3[1] - p1[1]) * (p2[0] - p1[0]) - (p2[1] - p1[1]) * (p3[0] - p1[0]));
    if(std::fabs(f) < epsilon) f = std::copysign(epsilon, f);
    float_t y = (d1 * d2 - e1 * e2) / f;

    radius2 = float_t((x - p1[0]) * (x - p1[0]) + (y - p1[1]) * (y - p1[1]));
    return Point2D<float_t>(x, y);
}

template <typename num_type>
inline bool isInCircumCircle(const Point2D<num_type> & p1, const Point2D<num_type> & p2, const Point2D<num_type> & p3, const Point2D<num_type> & p4, bool considerTouch)
{
    using float_t = float_type<num_type>;
    auto res = predicates::adaptive::inCircle<float_t>(p1[0], p1[1], p2[0], p2[1], p3[0], p3[1], p4[0], p4[1]);
    if(math::isNegative(res)) return false;
    if(!considerTouch && res == float_t(0)) return false;
    return true;
}

template <typename point_t, typename std::enable_if<traits::is_2d_point_t<point_t>::value, bool>::type>//todo point3d
inline coor_f<point_t> CircumRadius2ShortestEdgeRatioSq(const point_t & p1, const point_t & p2, const point_t & p3)
{
    using float_t = float_type<typename point_t::coor_t>;
    auto minEdgeSq = std::min(DistanceSq(p1, p2), DistanceSq(p2, p3));
    minEdgeSq = std::min(minEdgeSq, DistanceSq(p3, p1));

    float_t radiusSq(0);
    CircumCircle(p1, p2, p3, radiusSq);
    return radiusSq / minEdgeSq;
}

template <typename point_t, typename std::enable_if<traits::is_2d_point_t<point_t>::value, bool>::type>//todo point3d
inline coor_f<point_t> CircumRadius2ShortestEdgeRatio(const point_t & p1, const point_t & p2, const point_t & p3)
{
    return std::sqrt(CircumRadius2ShortestEdgeRatioSq(p1, p2, p3));
}

template <typename vector_t, typename std::enable_if<traits::is_2d_point_t<vector_t>::value, bool>::type>
inline coor_f<vector_t> Angle(const vector_t & v)
{
    auto result = std::atan2(v[1], v[0]);
    if(math::isNegative(result)) result += math::pi_2;
    return result; 
}

template <typename vector_t, typename std::enable_if<traits::is_2d_point_t<vector_t>::value, bool>::type>
inline coor_f<vector_t> Angle(const vector_t & a, const vector_t & b)
{
    auto result = std::atan2(CrossProduct(a, b), DotProduct(a, b));
    if(math::isNegative(result)) result += math::pi_2;
    return result;
}

template <typename vector_t, typename std::enable_if<traits::is_point_t<vector_t>::value, bool>::type>
inline coor_f<vector_t> InnerAngle(const vector_t & a, const vector_t & b)
{
    using float_t = float_type<typename vector_t::coor_t>;
    auto dot = DotProduct(a, b);
    float_t na = a.Norm2(), nb = b.Norm2();
    if(math::EQ(na, float_t(0)) || math::EQ(nb, float_t(0))) return 0;
    return std::acos(dot / (na * nb));
}

template <size_t vertex, typename triangle_t, typename std::enable_if<traits::is_triangle_t<triangle_t>::value, bool>::type>
inline coor_f<triangle_t> InteriorAngle(const triangle_t & t)
{
    const auto & p1 = t[(vertex + 2) % 3];
    const auto & p2 = t[(vertex + 0) % 3];
    const auto & p3 = t[(vertex + 1) % 3];
    return InnerAngle(p1, p2, p3);
}

template <typename triangle_t, typename std::enable_if<traits::is_triangle_t<triangle_t>::value, bool>::type>
inline coor_f<triangle_t> InteriorAngle(const triangle_t & t, const size_t vertex)
{
    const auto & p1 = t[(vertex + 2) % 3];
    const auto & p2 = t[(vertex + 0) % 3];
    const auto & p3 = t[(vertex + 1) % 3];
    return InnerAngle(p1, p2, p3);
}

template <typename point_t, typename std::enable_if<traits::is_point_t<point_t>::value, bool>::type>
inline std::array<coor_f<point_t>, 3> InteriorAngles(const point_t & p1, const point_t & p2, const point_t & p3)
{
    std::array<float_type<typename point_t::coor_t>, 3> angles;
    angles[0] = InnerAngle(p3, p1, p2);
    angles[1] = InnerAngle(p1, p2, p3);
    angles[2] = math::pi - angles[0] - angles[1];
    return angles;
}

template <typename num_type>
inline Point3D<float_type<num_type> > ClosestPointOnPlane(const Point3D<num_type> & p, const Plane<num_type> & plane)
{
    float_type<num_type> t = DotProduct(plane.Normal(), p.template Cast<float_type<num_type> >()) - plane.D();
    return p.template Cast<float_type<num_type> >() - plane.Normal() * t;
}

template <typename num_type>
inline float_type<num_type> PointPlaneDistance(const Point3D<num_type> & p, const Plane<num_type> & plane)
{
    return DotProduct(p.template Cast<float_type<num_type> >(), plane.Normal()) - plane.D();
}

template <typename point_t, typename std::enable_if<traits::is_point_t<point_t>::value, bool>::type>//return projection dist^2
inline coor_f<point_t> PointLineDistanceSq(const point_t & p, const point_t & a, const point_t & b)
{
    using num_t = typename point_t::coor_t;
    using float_t = float_type<num_t>;

    auto ab = b - a, ap = p - a;
    num_t e = DotProduct(ap, ab);
    num_t f = DotProduct(ab, ab);
    return std::max<float_t>(0, DotProduct(ap, ap) - float_t(e * e) / f);
}

template <typename point_t, typename segment_t,
          typename std::enable_if<traits::is_same_dim_t<point_t, segment_t>::value &&
                                  traits::is_same_coor_t<point_t, segment_t>::value, bool>::type>
inline point_t ClosestPointOnSegment(const point_t & p, const segment_t & s)
{
    using num_t = typename point_t::coor_t;
    using float_t = float_type<num_t>;

    const auto & a = s[0];
    const auto & b = s[1];
    auto ab = b - a;
    num_t t = DotProduct(p - a, ab);

    point_t d; float_t tf;
    if(math::LE(t, num_t(0))){ t = num_t(0); d = a; }
    else {
        num_t denom = DotProduct(ab, ab);
        if(math::GE(t, denom)) { t = num_t(1); d = b; }
        else { tf = float_t(t) / denom; d = a + ab * tf; }
    }
    return d;
}

template <typename point_t, typename segment_t,
          typename std::enable_if<traits::is_same_dim_t<point_t, segment_t>::value &&
                                  traits::is_same_coor_t<point_t, segment_t>::value, bool>::type>//return dist^2
inline coor_f<point_t> PointSegmentDistanceSq(const point_t & p, const segment_t & s)
{
    using num_t = typename point_t::coor_t;
    using float_t = float_type<num_t>;

    const auto & a = s[0];
    const auto & b = s[1];
    auto ab = b - a;
    auto ap = p - a;
    auto bp = p - b;
    num_t e = DotProduct(ap, ab);
    if(math::LE(e, num_t(0))) return DotProduct(ap, ap);
    num_t f = DotProduct(ab, ab);
    if(math::GE(e, f)) return DotProduct(bp, bp);
    return DotProduct(ap, ap) - float_t(e) * float_t(e) / f;
}

template <typename point_t>
inline coor_f<point_t> PointSegmentDistanceSq(const point_t & p, const point_t & a, const point_t & b)
{
    using segment_type = segment_type<point_t>;
    return PointSegmentDistanceSq(p, segment_type(a, b));
}

template <typename point_t, typename box_t,
          typename std::enable_if<traits::is_same_dim_t<point_t, box_t>::value &&
                                  traits::is_same_coor_t<point_t, box_t>::value, bool>::type>
inline point_t ClosestPointInBox(const point_t & p, const box_t & b)
{
    using num_t = typename point_t::coor_t;
    point_t c;
    for(size_t i = 0; i < point_t::dim; ++i){
        num_t t = p[i];
        t = std::max(t, b[0][i]);
        t = std::min(t, b[1][i]);
        c[i] = t;
    }
    return c;
}

template <typename point_t, typename box_t,
          typename std::enable_if<traits::is_same_dim_t<point_t, box_t>::value &&
                                  traits::is_same_coor_t<point_t, box_t>::value, bool>::type>//return dist^2
inline typename point_t::coor_t PointBoxDistanceSq(const point_t & p, const box_t & b)
{
    typename point_t::coor_t distSq(0);
    for(size_t i = 0; i < point_t::dim; ++i){
        auto t = p[i];
        if(t < b[0][i]) distSq += (b[0][i] - t) * (b[0][i] - t);
        if(t > b[1][i]) distSq += (t - b[1][i]) * (t - b[1][i]);
    }
    return distSq;
}

template <typename point_t, typename segment_t, 
          typename std::enable_if<traits::is_same_dim_t<point_t, segment_t>::value, bool>::type>
inline PointLineLocation GetPointSegmentLocation(const point_t & p, const segment_t & seg)
{
    return GetPointLineLocation(p, seg[0], seg[1]);
}

template <typename point_t, typename line_t, 
          typename std::enable_if<traits::is_same_dim_t<point_t, line_t>::value, bool>::type>
inline PointLineLocation GetPointLineLocation(const point_t & p, const line_t & line)
{
    using float_t = float_type<typename point_t::coor_t>;
    float_t res = line_t::SideValue(line, p);
    if(math::EQ(res, float_t(0))) return PointLineLocation::OnLine;
    else if(0 < res) return PointLineLocation::Left;
    else return PointLineLocation::Right;
} 

template <typename point_t1, typename point_t2,
          typename std::enable_if<traits::is_same_dim_t<point_t1, point_t2>::value, bool>::type>
inline PointLineLocation GetPointLineLocation(const point_t1 & p, const point_t2 & v1, const point_t2 & v2)
{
    static_assert(point_t1::dim == 2, "point line location only applied for 2d");
    using float_t = float_type<typename traits::common_coor_t<typename point_t1::coor_t, typename point_t2::coor_t>::type>;
    auto res = predicates::adaptive::Orient2D<float_t>(v1[0], v1[1], v2[0], v2[1], p[0], p[1]);
    if(float_t(0) == res) return PointLineLocation::OnLine;
    else if(0 < res) return PointLineLocation::Left;
    else return PointLineLocation::Right;
}

template <typename point_t1, typename point_t2,
          typename std::enable_if<traits::is_same_dim_t<point_t1, point_t2>::value, bool>::type>
inline PointTriangleLocation GetPointTriangleLocation(const point_t1 & p, const point_t2 & v1, const point_t2 & v2, const point_t2 & v3)
{
    return GetPointTriangleLocation(p, triangle_type<point_t2>(v1, v2, v3));
}

template <typename point_t, typename triangle_t,
          typename std::enable_if<traits::is_same_dim_t<point_t, triangle_t>::value, bool>::type>
inline PointTriangleLocation GetPointTriangleLocation(const point_t & p, const triangle_t & tri)
{
    std::array<bool, 3> sign{true, true, true};
    for(int i = 0; i < 3; ++i){
        auto loc = GetPointLineLocation(p, tri[i], tri[(i + 1) % 3]);
        if(PointLineLocation::OnLine == loc)
            return PointTriangleLocation(i + 1);
        
        sign[i] = PointLineLocation::Right == loc ? true : false;
    }
    if(sign[0] == sign[1] && sign[1] == sign[2]) return PointTriangleLocation::Inside;
    return PointTriangleLocation::Outside;
}

template <typename segment_t1, typename segment_t2,
          typename std::enable_if<traits::is_2d_geometry_t<segment_t1>::value && 
                                  traits::is_2d_geometry_t<segment_t2>::value, bool>::type = true>
inline bool IntersectsImp(const segment_t1 & s1, const segment_t2 & s2, bool considerTouch, std::false_type)
{
    using common_t = typename traits::common_coor_t<typename segment_t1::coor_t, typename segment_t2::coor_t>::type;
    return IntersectsImp(s1.template Cast<common_t>(), s2.template Cast<common_t>(), considerTouch, std::true_type{});
}

template <typename segment_t1, typename segment_t2,
          typename std::enable_if<traits::is_2d_geometry_t<segment_t1>::value && 
                                  traits::is_2d_geometry_t<segment_t2>::value, bool>::type = true>
inline bool IntersectsImp(const segment_t1 & s1, const segment_t2 & s2, bool considerTouch, std::true_type)
{
    return boost::polygon::intersects(s1, s2, considerTouch);
}

template <typename segment_t1, typename segment_t2,
          typename std::enable_if<traits::is_2d_geometry_t<segment_t1>::value && 
                                  traits::is_2d_geometry_t<segment_t2>::value, bool>::type>
inline bool Intersects(const segment_t1 & s1, const segment_t2 & s2, bool considerTouch)
{
    return IntersectsImp(s1, s2, considerTouch, typename traits::is_same_coor_t<segment_t1, segment_t2>{});
}

template <typename num_type>
inline bool Intersection(const Segment2D<num_type> & s1, const Segment2D<num_type> & s2, std::vector<Point2D<num_type> > & points)
{
    return boost::geometry::intersection(s1, s2, points);
}

template <typename num_type>
inline bool Intersection(const Segment2D<num_type> & s, const Line2D<num_type> & line, std::vector<Point2D<num_type> > & points)
{
    points.clear();
    using float_t = float_type<num_type>;
    float_t s1 = Line2D<num_type>::SideValue(line, s[0]);
    float_t s2 = Line2D<num_type>::SideValue(line, s[1]);
    if(math::EQ(s1, float_t(0))) points.push_back(s[0]);
    if(math::EQ(s2, float_t(0))) points.push_back(s[1]);
    if(points.size()) return true;

    if(std::signbit(s1) == std::signbit(s2)) return false;
    auto p = s[0] + (s[1] - s[0]) * (s1 / (s1 - s2));
    points.push_back(p);
    return true;
}

template <typename num_type>
inline bool Intersection(const Segment3D<num_type> & s, const Plane<num_type> & plane, Point3D<float_type<num_type> > & point)
{
    using float_t = float_type<num_type>;
    auto n = plane.Normal();
    auto vec = (s[1] - s[0]).template Cast<float_t>();
    auto t = (plane.D() - DotProduct(n, s[0].template Cast<float_t>())) / DotProduct(n, vec);
    if(math::GE<float_t>(t, 0) && math::LE<float_t>(t, 1)) {
        point = s[0].template Cast<float_t>() + vec.template Cast<float_t>() * t;
        return true;
    }
    return false;
}

template <typename geometry_t1, typename geometry_t2,
          typename std::enable_if<(traits::is_segment_t<geometry_t1>::value  && traits::is_point_t<geometry_t2>::value)   ||
                                  (traits::is_segment_t<geometry_t1>::value  && traits::is_segment_t<geometry_t2>::value) ||
                                  (traits::is_triangle_t<geometry_t1>::value && traits::is_point_t<geometry_t2>::value)   ||
                                  (traits::is_box_t<geometry_t1>::value      && traits::is_point_t<geometry_t2>::value)   ||
                                  (traits::is_box_t<geometry_t1>::value      && traits::is_box_t<geometry_t2>::value)     ||
                                  (traits::is_polygon_t<geometry_t1>::value  && traits::is_point_t<geometry_t2>::value)    ||
                                  (traits::is_polygon_with_holes_t<geometry_t1>::value && traits::is_point_t<geometry_t2>::value), bool>::type = true>
inline bool ContainsImp2D(const geometry_t1 & g1, const geometry_t2 & g2, bool considerTouch)
{
    return boost::polygon::contains(g1, g2, considerTouch);
}

template <typename geometry_t1, typename geometry_t2, 
          typename std::enable_if<(traits::is_triangle_t<geometry_t1>::value ||
                                   traits::is_box_t<geometry_t1>::value ||
                                   traits::is_polygon_t<geometry_t1>::value) && traits::is_segment_t<geometry_t2>::value, bool>::type = true>
inline bool ContainsImp2D(const geometry_t1 & g1, const geometry_t2 & g2, bool considerTouch)
{
    return ContainsImp2D(g1, g2[0], considerTouch) &&
           ContainsImp2D(g1, g2[1], considerTouch);
}

template <typename geometry_t1, typename geometry_t2,
          typename std::enable_if<(traits::is_triangle_t<geometry_t1>::value ||
                                   traits::is_box_t<geometry_t1>::value ||
                                   traits::is_polygon_t<geometry_t1>::value) && traits::is_triangle_t<geometry_t2>::value, bool>::type = true>
inline bool ContainsImp2D(const geometry_t1 & g1, const geometry_t2 & g2, bool considerTouch)
{
    return ContainsImp2D(g1, g2[0], considerTouch) &&
           ContainsImp2D(g1, g2[1], considerTouch) &&
           ContainsImp2D(g1, g2[2], considerTouch);  
}

template <typename geometry_t1, typename geometry_t2,
          typename std::enable_if<(traits::is_triangle_t<geometry_t1>::value ||
                                   traits::is_polygon_t<geometry_t1>::value) && traits::is_box_t<geometry_t2>::value, bool>::type = true>
inline bool ContainsImp2D(const geometry_t1 & g1, const geometry_t2 & g2, bool considerTouch)
{
    using Point = typename geometry_t2::point_t;
    return ContainsImp2D(g1, g2[0], considerTouch) &&
           ContainsImp2D(g1, g2[1], considerTouch) &&
           ContainsImp2D(g1, Point(g2[1][0], g2[0][1]), considerTouch) &&
           ContainsImp2D(g1, Point(g2[0][0], g2[1][1]), considerTouch);
}

template <typename geometry_t1, typename geometry_t2,
          typename std::enable_if<(traits::is_polygon_t<geometry_t1>::value || traits::is_polygon_with_holes_t<geometry_t1>::value) &&
                                  (traits::is_polygon_t<geometry_t2>::value || traits::is_polygon_with_holes_t<geometry_t2>::value), bool>::type = true>
inline bool ContainsImp2D(const geometry_t1 & g1, const geometry_t2 & g2, bool considerTouch)
{
    //not tested yet!
    if(considerTouch)
        return boost::geometry::covered_by(g2, g1);
    else return boost::geometry::within(g2, g1);
}

template <typename geometry_t1, typename geometry_t2,
          typename std::enable_if<traits::is_box_t<geometry_t1>::value && traits::is_point_t<geometry_t2>::value, bool>::type = true>
inline bool ContainsImp3D(const geometry_t1 & g1, const geometry_t2 & g2, bool considerTouch)
{
    if(considerTouch){
        return math::LE(g1[0][0], g2[0]) && math::GE(g1[1][0], g2[0]) &&
               math::LE(g1[0][1], g2[1]) && math::GE(g1[1][1], g2[1]) && 
               math::LE(g1[0][2], g2[2]) && math::GE(g1[1][2], g2[2]);
    }
    else{
        return math::LE(g1[0][0], g2[0]) && math::GE(g1[1][0], g2[0]) &&
               math::LE(g1[0][1], g2[1]) && math::GE(g1[1][1], g2[1]) && 
               math::LE(g1[0][2], g2[2]) && math::GE(g1[1][2], g2[2]);
    }
}

template <typename geometry_t1, typename geometry_t2,
          typename std::enable_if<traits::is_box_t<geometry_t1>::value &&
                                 (traits::is_segment_t<geometry_t2>::value ||
                                  traits::is_box_t<geometry_t2>::value), bool>::type = true>
inline bool ContainsImp3D(const geometry_t1 & g1, const geometry_t2 & g2, bool considerTouch)
{
    return ContainsImp3D(g1, g2[0], considerTouch) &&
           ContainsImp3D(g1, g2[1], considerTouch);
}

template <typename geometry_t1, typename geometry_t2,
          typename std::enable_if<traits::is_geometry_t<geometry_t1>::value &&
                                  traits::is_geometry_t<geometry_t2>::value, bool>::type = true>
inline bool ContainsDetail(const geometry_t1 & g1, const geometry_t2 & g2, bool considerTouch, traits::different_coor_t)
{
    using common_t = typename traits::common_coor_t<typename geometry_t1::coor_t, typename geometry_t2::coor_t>::type;
    if constexpr (std::is_same<common_t, typename geometry_t1::coor_t>::value)
        return ContainsDetail(g1, g2.template Cast<common_t>(), considerTouch, traits::same_coor_t{});
    else return ContainsDetail(g1.template Cast<common_t>(), g2, considerTouch, traits::same_coor_t{});
}

template <typename geometry_t1, typename geometry_t2,
          typename std::enable_if<traits::is_2d_geometry_t<geometry_t1>::value &&
                                  traits::is_2d_geometry_t<geometry_t2>::value, bool>::type = true>
inline bool ContainsDetail(const geometry_t1 & g1, const geometry_t2 & g2, bool considerTouch, traits::same_coor_t)
{
    return ContainsImp2D(g1, g2, considerTouch);
}

template <typename geometry_t1, typename geometry_t2,
          typename std::enable_if<traits::is_3d_geometry_t<geometry_t1>::value &&
                                  traits::is_3d_geometry_t<geometry_t2>::value, bool>::type = true>
inline bool ContainsDetail(const geometry_t1 & g1, const geometry_t2 & g2, bool considerTouch, traits::same_coor_t)
{
    return ContainsImp3D(g1, g2, considerTouch);
}

template <typename geometry_t1, typename geometry_t2,
          typename std::enable_if<traits::is_geometry_t<geometry_t1>::value &&
                                  traits::is_geometry_t<geometry_t2>::value, bool>::type>
inline bool Contains(const geometry_t1 & g1, const geometry_t2 & g2, bool considerTouch)
{
    return ContainsDetail(g1, g2, considerTouch, typename traits::is_same_coor_t<geometry_t1, geometry_t2>::type{});
}

template <typename segment_t, typename std::enable_if<traits::is_segment_t<segment_t>::value, bool>::type>
inline box_type<segment_t> Extent(const segment_t & segment)
{
    box_type<segment_t> bbox;
    bbox |= segment[0];
    bbox |= segment[1];
    return bbox;
}

template <typename triangle_t, typename std::enable_if<traits::is_triangle_t<triangle_t>::value, bool>::type>
inline box_type<triangle_t> Extent(const triangle_t & triangle)
{
    box_type<triangle_t> bbox;
    bbox |= triangle[0];
    bbox |= triangle[1];
    bbox |= triangle[2];
    return bbox;
}

template <typename polygon_t, typename std::enable_if<traits::is_polygon_t<polygon_t>::value ||
                                                       traits::is_polygon_with_holes_t<polygon_t>::value, bool>::type>
inline box_type<polygon_t> Extent(const polygon_t & polygon)
{
    box_type<polygon_t> bbox;
    boost::polygon::extents(bbox, polygon);
    return bbox;
}

template <typename iterator,
          typename std::enable_if<traits::is_geometry_t<
          typename std::iterator_traits<iterator>::value_type>::value, bool>::type>
inline box_type<typename std::iterator_traits<iterator>::value_type> Extent(iterator begin, iterator end)
{
    box_type<typename std::iterator_traits<iterator>::value_type> bbox;
    for(auto iter = begin; iter != end; ++iter){
        typename std::iterator_traits<iterator>::reference r = *iter;
        bbox |= Extent(r);
    }
    return bbox;
}

template <typename num_type>
inline Box2D<num_type> Extent(const Polyline2D<num_type> & polyline)
{
    return boost::geometry::return_envelope<Box2D<num_type> >(polyline);
}

///@brief scales a box by factor
template <typename box_type>
inline void Scale(box_type & box, coor_f<box_type> factor)
{
    auto center = box.Center();
    auto fbox = box.template Cast<coor_f<box_type>>();
    fbox -= center;
    fbox *= factor;
    fbox += center;
    box = fbox.template Cast<typename box_type::coor_t>();
}

template <typename num_type>
inline void Simplify(Polygon2D<num_type> & polygon, std::vector<Polygon2D<num_type> > & holes)
{
    holes.clear();
    if(polygon.Front() == polygon.Back()) polygon.PopBack();
    if(polygon.Size() <= 3) return;

    size_t size = polygon.Size();
    struct PtNode { size_t prev; size_t next; };
    std::vector<PtNode> nodeList(size);
    for(size_t i = 0; i < size; ++i) {
        nodeList[i].prev = (i + size - 1) % size;
        nodeList[i].next = (i + 1) % size;
    }

    size_t res;
    std::vector<Polygon2D<num_type>> tmp;
    detail::Point2DIndexMap<num_type> ptMap;
    for(size_t i = 0; i < polygon.Size(); ++i) {
        auto found = ptMap.Find(polygon[i], res);
        if(found) {
            size_t prev = res;
            size_t curr = i;
            size_t next = nodeList[curr].next;

            nodeList[curr].next = prev;
            Polygon2D<num_type> hole;
            size_t start = prev;
            size_t index = start;
            while(start != nodeList[index].next) {
                hole << polygon[index];
                index = nodeList[index].next;
            }
            Simplify(hole, tmp);
            if(hole.Size() >= 3)
                holes.emplace_back(std::move(hole));
            prev = nodeList[prev].prev;
            nodeList[prev].next = curr;
            nodeList[curr].prev = prev;
            nodeList[curr].next = next;
        }
        ptMap.Insert(polygon[i], i);
    }

    size_t start = size - 1;
    size_t index = start;
    Polygon2D<num_type> simplified;
    while(start != nodeList[index].next) {
        simplified << polygon[index];
        index = nodeList[index].next;
    }
    simplified << polygon[index];
    std::swap(polygon, simplified);
}

template <typename num_type>
inline Polygon2D<num_type> RoundCorners(const Polygon2D<num_type> & polygon, num_type radius, size_t circleDiv)
{
    auto result = polygon;
    auto deltaAngle = math::pi_2 / circleDiv;
    if (not result.isCCW()) result.Reverse();
    auto & points = result.GetPoints();
    auto iter = points.begin();
    using float_t = float_type<num_type>;
    Point2D<float_t> fp1, fp2;
    bool stop = false;
    while (not stop) {
        auto it1 = iter;
        auto it2 = it1; it2++;
        if (it2 == points.end()) {
            it2 = points.begin();
            stop = true;
        }
        auto it3 = it2; it3++;
        if (it3 == points.end()) it3 = points.begin();
        auto o = InscribedCircle(*it1, *it2, *it3, radius, fp1, fp2);
        if (DistanceSq(*it1, *it2) < DistanceSq(fp1.template Cast<num_type>(), *it2) ||
            DistanceSq(*it3, *it2) < DistanceSq(fp2.template Cast<num_type>(), *it2)) {
            iter++; continue;
        }
        else {
            auto startAngle = Angle(fp1 - o);
            auto theta = InnerAngle(fp1, o, fp2); 
            size_t divide = theta / deltaAngle;
            std::vector<Point2D<num_type>> corners(divide);
            auto orient = CrossProduct(
                it2->template Cast<float_t>() - fp1,
                fp2 - it2->template Cast<float_t>()) > 0 ? 1.0 : -1.0;
            for (size_t i = 0; i < divide; ++i) {
                corners[i][0] = std::cos(startAngle + i * orient * deltaAngle) * radius + o[0];
                corners[i][1] = std::sin(startAngle + i * orient * deltaAngle) * radius + o[1];
            }
            if (not math::EQ(divide * deltaAngle, theta))
                corners.emplace_back(fp2.template Cast<num_type>());
            auto it = points.erase(it2);
            iter = points.insert(it, corners.begin(), corners.end());
            std::advance(iter, corners.size() - 1);
        }
    }
    return result;
}

template <typename polygon_t, template <typename, typename> class container, template <typename> class allocator>
inline polygon_t ConvexHull(const container<polygon_t, allocator<polygon_t> > & polygons)
{
    using multi_point_t = boost::geometry::model::multi_point<typename polygon_t::point_t>;
    multi_point_t points;
    for(const auto & polygon : polygons){
        const auto & pts = polygon.GetPoints();
        for(const auto & p : pts) points.push_back(p);
    }

    polygon_t hull;
    boost::geometry::convex_hull<multi_point_t, polygon_t>(points, hull);
    return hull;
}

}//namespace geometry
}//namespace generic