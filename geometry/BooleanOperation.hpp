#ifndef GENERIC_GEOMETRY_BOOLEANOPERATION_HPP
#define GENERIC_GEOMETRY_BOOLEANOPERATION_HPP
#include "BoostPolygonRegister.hpp"
#include "GeometryTraits.hpp"
#include <boost/polygon/polygon_set_traits.hpp>
#include <boost/polygon/polygon_set_data.hpp>
namespace generic  {
namespace geometry {
namespace boolean  {

template <typename num_type>
using PolygonSet2D = boost::polygon::polygon_set_data<num_type>;

template <typename geometry_t1, typename geometry_t2, 
          template <typename, typename> class container, 
          template <typename> class allocator = std::allocator,
          typename std::enable_if<traits::is_2d_surf_geom_t<geometry_t1>::value &&
                                  traits::is_2d_surf_geom_t<geometry_t2>::value, bool>::type = true>
inline void Unite(const geometry_t1 & g1, const geometry_t2 & g2, container<polygon_type<geometry_t1>, allocator<polygon_type<geometry_t1> > > & results)
{
    using namespace boost::polygon::operators;
    PolygonSet2D<typename geometry_t1::coor_t> polygonSet;    
    polygonSet += g1;
    polygonSet += g2;
    results.clear();
    polygonSet.get(results);
}

template <typename iterator, 
          template <typename, typename> class container,
          template <typename> class allocator = std::allocator,
          typename std::enable_if<traits::is_2d_surf_geom_t<
          typename std::iterator_traits<iterator>::value_type>::value, bool>::type = true>
inline void Unite(iterator begin, iterator end,
                  container<polygon_type<typename std::iterator_traits<iterator>::value_type>,
                  allocator<polygon_type<typename std::iterator_traits<iterator>::value_type> > > & results)
{
    using namespace boost::polygon::operators;
    using coor_t = typename std::iterator_traits<iterator>::value_type::coor_t;
    PolygonSet2D<coor_t> polygonSet;
    for(auto iter = begin; iter != end; ++iter) 
        polygonSet += *iter;
    results.clear();
    polygonSet.get(results);
}

template <typename geometry_t1, typename geometry_t2, 
          template <typename, typename> class container, 
          template <typename> class allocator = std::allocator,
          typename std::enable_if<traits::is_2d_surf_geom_t<geometry_t1>::value &&
                                  traits::is_2d_surf_geom_t<geometry_t2>::value, bool>::type = true>
inline void Intersect(const geometry_t1 & g1, const geometry_t2 & g2, container<polygon_type<geometry_t1>, allocator<polygon_type<geometry_t1> > > & results)
{
    using namespace boost::polygon::operators;
    PolygonSet2D<typename geometry_t1::coor_t> polygonSet;    
    polygonSet += g1;
    polygonSet &= g2;
    results.clear();
    polygonSet.get(results);
}

template <typename iterator1, typename iterator2,
          template <typename, typename> class container,
          template <typename> class allocator = std::allocator,
          typename std::enable_if<traits::is_2d_surf_geom_t<typename std::iterator_traits<iterator1>::value_type>::value &&
                                  traits::is_2d_surf_geom_t<typename std::iterator_traits<iterator2>::value_type>::value, bool>::type = true>
inline void Intersect(iterator1 begin1, iterator1 end1, iterator2 begin2, iterator2 end2,
                      container<polygon_type<typename std::iterator_traits<iterator1>::value_type>,
                      allocator<polygon_type<typename std::iterator_traits<iterator1>::value_type> > > & results)
{
    using namespace boost::polygon::operators;
    using coor_t = typename std::iterator_traits<iterator1>::value_type::coor_t;
    PolygonSet2D<coor_t> polygonSet1, polygonSet2;
    for(auto iter1 = begin1; iter1 != end1; ++iter1) 
        polygonSet1 += *iter1;
    for(auto iter2 = begin2; iter2 != end2; ++iter2)
        polygonSet2 += *iter2;
    polygonSet1 *= polygonSet2;
    results.clear();
    polygonSet1.get(results);
}

template <typename geometry_t1, typename geometry_t2, 
          template <typename, typename> class container, 
          template <typename> class allocator = std::allocator,
          typename std::enable_if<traits::is_2d_surf_geom_t<geometry_t1>::value &&
                                  traits::is_2d_surf_geom_t<geometry_t2>::value, bool>::type = true>
inline void Subtract(const geometry_t1 & g1, const geometry_t2 & g2, container<polygon_type<geometry_t1>, allocator<polygon_type<geometry_t1> > > & results)
{
    using namespace boost::polygon::operators;
    PolygonSet2D<typename geometry_t1::coor_t> polygonSet;    
    polygonSet += g1;
    polygonSet -= g2;
    results.clear();
    polygonSet.get(results);
}

template <typename iterator1, typename iterator2,
          template <typename, typename> class container,
          template <typename> class allocator = std::allocator,
          typename std::enable_if<traits::is_2d_surf_geom_t<typename std::iterator_traits<iterator1>::value_type>::value &&
                                  traits::is_2d_surf_geom_t<typename std::iterator_traits<iterator2>::value_type>::value, bool>::type = true>
inline void Subtract(iterator1 begin1, iterator1 end1, iterator2 begin2, iterator2 end2,
                     container<polygon_type<typename std::iterator_traits<iterator1>::value_type>,
                     allocator<polygon_type<typename std::iterator_traits<iterator1>::value_type> > > & results)
{
    using namespace boost::polygon::operators;
    using coor_t = typename std::iterator_traits<iterator1>::value_type::coor_t;
    PolygonSet2D<coor_t> polygonSet;
    for(auto iter = begin1; iter != end1; ++iter) 
        polygonSet += *iter;
    for(auto iter = begin2; iter != end2; ++iter)
        polygonSet -= *iter;
    results.clear();
    polygonSet.get(results);
}

template <typename geometry_t1, typename geometry_t2, 
          template <typename, typename> class container, 
          template <typename> class allocator = std::allocator,
          typename std::enable_if<traits::is_2d_surf_geom_t<geometry_t1>::value &&
                                  traits::is_2d_surf_geom_t<geometry_t2>::value, bool>::type = true>
inline void Xor(const geometry_t1 & g1, const geometry_t2 & g2, container<polygon_type<geometry_t1>, allocator<polygon_type<geometry_t1> > > & results)
{
    using namespace boost::polygon::operators;
    PolygonSet2D<typename geometry_t1::coor_t> polygonSet;    
    polygonSet += g1;
    polygonSet ^= g2;
    results.clear();
    polygonSet.get(results);
}

template <typename iterator1, typename iterator2,
          template <typename, typename> class container,
          template <typename> class allocator = std::allocator,
          typename std::enable_if<traits::is_2d_surf_geom_t<typename std::iterator_traits<iterator1>::value_type>::value &&
                                  traits::is_2d_surf_geom_t<typename std::iterator_traits<iterator2>::value_type>::value, bool>::type = true>
inline void Xor(iterator1 begin1, iterator1 end1, iterator2 begin2, iterator2 end2,
                container<polygon_type<typename std::iterator_traits<iterator1>::value_type>,
                allocator<polygon_type<typename std::iterator_traits<iterator1>::value_type> > > & results)
{
    using namespace boost::polygon::operators;
    using coor_t = typename std::iterator_traits<iterator1>::value_type::coor_t;
    PolygonSet2D<coor_t> polygonSet1, polygonSet2;
    for(auto iter1 = begin1; iter1 != end1; ++iter1) 
        polygonSet1 += *iter1;
    for(auto iter2 = begin2; iter2 != end2; ++iter2)
        polygonSet2 += *iter2;
    polygonSet1 ^= polygonSet2;
    results.clear();
    polygonSet1.get(results);
}

} //namespace boolean
} //namespace geometry
} //namespace generic
#endif//GENERIC_GEOMETRY_BOOLEANOPERATION_HPP