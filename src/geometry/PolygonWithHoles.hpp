/**
 * @file PolygonWithHoles.hpp
 * @author bwu
 * @brief Model of polygon with holes concept
 * @version 0.1
 * @date 2022-02-22
 */
#pragma once
#include "Polygon.hpp"
#include <list>
namespace generic  {
namespace geometry {

template <typename num_type>
class PolygonWithHoles2D
{
public:
    const static size_t dim = 2;
    using coor_t = num_type;
    using point_t = Point2D<coor_t>;
    using outline_type = Polygon2D<num_type>;
    using point_iterator = typename outline_type::point_iterator;
    using const_point_iterator = typename outline_type::const_point_iterator;
    using hole_type = Polygon2D<coor_t>;
    using hole_container = std::vector<hole_type>;
    using hole_iterator = typename hole_container::iterator;
    using const_hole_iterator = typename hole_container::const_iterator;

    outline_type outline;
    hole_container holes;

    PolygonWithHoles2D(){}
    
    ///@brief clears outline points and holes data
    void Clear() { outline.clear(); holes.clear(); }
    ///@brief checks whether has hole
    bool hasHole() const { return HolesSize() > 0; }
    ///@brief gets total hole numbers
    size_t HolesSize() const { return holes.size(); }

    hole_iterator BeginHoles() { return holes.begin(); }
    hole_iterator EndHoles() { return holes.end(); }
    const_hole_iterator ConstBeginHoles() const { return holes.begin(); }
    const_hole_iterator ConstEndHoles() const { return holes.end(); }  

    ///@brief converts to polygon with holes with other number type explicitly
    template <typename other_num_type>
    PolygonWithHoles2D<other_num_type> Cast() const;  
};

template <typename num_type>
template <typename other_num_type>
inline PolygonWithHoles2D<other_num_type> PolygonWithHoles2D<num_type>::Cast() const
{
    PolygonWithHoles2D<other_num_type> pwh;
    pwh.outline = outline.template Cast<other_num_type>();
    for(const auto & hole : holes)
        pwh.holes.emplace_back(hole.template Cast<other_num_type>());
    return pwh;
}

}//namespace geometry
}//namespace generic