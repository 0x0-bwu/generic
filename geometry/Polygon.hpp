#ifndef GENERIC_GEOMETRY_POLYGON_HPP
#define GENERIC_GEOMETRY_POLYGON_HPP
#include "Point.hpp"
#include <vector>
namespace generic  {
namespace geometry {
using generic::common::float_type;

template <typename num_type>
using Polyline2D = std::vector<Point2D<num_type> >;

template <typename num_type>
class Polygon2D
{
    using float_t = float_type<num_type>;
public:
    const static size_t dim = 2;
    using coor_t = num_type;
    using point_t = Point2D<coor_t>;
    using value_type = point_t;
    using point_container = std::vector<point_t>;
    using point_iterator = typename point_container::iterator;
    using const_point_iterator = typename point_container::const_iterator;

    Polygon2D(){}
    Polygon2D & operator<< (point_t point) { m_points.push_back(point); return  *this; }

    Point2D<num_type> & operator[] (size_t i) { return m_points[i]; }
    const Point2D<num_type> & operator[] (size_t i) const { return m_points[i]; }
    
    void Clear() { m_points.clear(); }
    size_t Size() const { return m_points.size(); }
    float_t Area() const;
    bool isCCW() const;

    template <typename input_iterator>
    void Insert(point_iterator position, input_iterator begin, input_iterator end);
    void Insert(point_iterator position, const point_t & point);

    point_iterator Begin() { return m_points.begin(); }
    point_iterator End() { return m_points.end(); }
    const_point_iterator ConstBegin() const { return m_points.begin(); }
    const_point_iterator ConstEnd() const { return m_points.end(); }
    point_t & Front() { return m_points.front(); }
    point_t & Back() { return m_points.back(); }
    const point_t & ConstFront() const { return m_points.front(); }
    const point_t & ConstBack() const { return m_points.back(); }
    void Resize(size_t size) { m_points.resize(size); }
    void PopBack() { if(m_points.size()) m_points.pop_back(); }
    point_container & GetPoints() { return m_points; }
    const point_container & GetPoints() const { return m_points; }
    void Reverse() { std::reverse(m_points.begin(), m_points.end()); }
    void Set(point_container points) { m_points = std::move(points); }

    template <typename other_num_type>
    Polygon2D<other_num_type> Cast() const;

    static void Clean(Polygon2D<num_type> & polygon);
    
private:
    point_container m_points;
};

template <typename num_type>
inline float_type<num_type> Polygon2D<num_type>::Area() const
{
    float_t area = 0;
    size_t size = Size();
    if(size < 3) return area;
    for(size_t i = 0, j = size - 1; i < size; ++i){
        area += (m_points[j][0] + m_points[i][0]) * (m_points[j][1] - m_points[i][1]);
        j = i;
    }
    return -area * 0.5;
}

template <typename num_type>
inline bool Polygon2D<num_type>::isCCW() const
{
    return !math::isNegative(Area());
}

template <typename num_type>
template <typename input_iterator>
inline void Polygon2D<num_type>::Insert(point_iterator position, input_iterator begin, input_iterator end)
{
    m_points.insert(position, begin, end);
}

template <typename num_type>
inline void Polygon2D<num_type>::Insert(point_iterator position, const point_t & point)
{
    m_points.insert(position, point);
}

template <typename num_type>
template <typename other_num_type>
inline Polygon2D<other_num_type> Polygon2D<num_type>::Cast() const
{
    Polygon2D<other_num_type> polygon;
    for(const auto & p : m_points)
        polygon << p.template Cast<other_num_type>();
    return polygon;
}

template <typename num_type>
inline void Polygon2D<num_type>::Clean(Polygon2D<num_type> & polygon)
{
    auto size = polygon.Size();
    point_container points;
    for(size_t i = 0; i < size; ++i){
        const auto & prev = polygon[(i + size - 1) % size];
        const auto & curr = polygon[i];
        const auto & next = polygon[(i + 1) % size];
        if(curr == next) continue;
        if(math::EQ(prev[0], curr[0]) && math::EQ(curr[0], next[0])) continue;
        if(math::EQ(prev[1], curr[1]) && math::EQ(curr[1], next[1])) continue;
        points.push_back(curr);   
    }
    polygon.Set(std::move(points));
    if(polygon.Size() != size) Clean(polygon);
}

} //namespace geometry
} //namespace generic
#endif//GENERIC_GEOMETRY_POLYGON_HPP