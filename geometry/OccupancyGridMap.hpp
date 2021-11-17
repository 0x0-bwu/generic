#ifndef GENERIC_GEOMETRY_OCCUPANCYGRIDMAP_HPP
#define GENERIC_GEOMETRY_OCCUPANCYGRIDMAP_HPP
#include "thread/ThreadPool.hpp"
#include "math/MathUtility.hpp"
#include "common/Exception.hpp"
#include "tools/Tools.hpp"
#include "BooleanOperation.hpp"
#include "Triangulation.hpp"
#include "Rasterization.hpp"
#include "Triangulator.hpp"
#include "Geometries.hpp"
#include "GeometryIO.hpp"
#include "Trapezoid.hpp"
#include "Utility.hpp"
#include <boost/lockfree/queue.hpp>
#include <vector>
#include <atomic>

#ifdef BOOST_GIL_IO_PNG_SUPPORT
#include "tools/Color.hpp"
#include <boost/gil/extension/io/png.hpp>
#include <boost/gil.hpp>
#include <png.h>
#endif

namespace generic {
namespace geometry {

template <typename Occupancy>
class OccupancyGridMap
{
    using Container = std::vector<std::vector<Occupancy> >;
public:
    using ResultType = Occupancy;
    OccupancyGridMap(size_t width, size_t height)
     : m_width(width),
       m_height(height),
       m_grids(width, std::vector<Occupancy>(height, Occupancy{}))
    {}
    ~OccupancyGridMap() = default;

    Occupancy & operator[] (size_t i) { return m_grids[i / m_height][i % m_height]; }
    const Occupancy & operator[](size_t i) const { return m_grids[i / m_height][i % m_height]; }

    Occupancy & operator() (size_t x, size_t y) { return m_grids[x][y]; }
    const Occupancy & operator() (size_t x, size_t y) const { return m_grids[x][y]; }

    void Resize(size_t width, size_t height)
    {
        m_width = width;
        m_height = height;
        m_grids.resize(width);
        for(auto & col : m_grids)
            col.resize(height);
    }

    void Reset()
    {
        for(auto & col : m_grids){
            std::fill(col.begin(), col.end(), 0.0f);
        }
    }

    void Clear() { m_grids.clear(); }

    size_t Width() const { return m_width; }
    size_t Height() const { return m_height; }
    size_t Size() const { return Width() * Height(); }

    template <typename Compare>
    const Occupancy & MaxOccupancy(Compare && cmp) const
    {
        auto [i, j] = MaxElement(std::forward<Compare>(cmp));
        return m_grids[i][j];
    }

    template <typename Compare>
    std::pair<size_t, size_t> MaxElement(Compare && cmp) const
    {
        size_t xMax(0), yMax(0);
        for(size_t i = 0; i < m_width; ++i){
            for(size_t j = 0; j < m_height; ++j){
                if(cmp(m_grids[i][j], m_grids[xMax][yMax])){
                    xMax = i; yMax = j;
                }
            }
        }
        return std::make_pair(xMax, yMax);
    }

#ifdef BOOST_GIL_IO_PNG_SUPPORT
    template <typename RGBaFunc>
    bool WriteImgProfile(const std::string & filename, RGBaFunc && rgbaFunc)
    {
        using namespace boost::gil;
        rgb8_image_t img(m_width, m_height);
        rgb8_image_t::view_t v = view(img);
        for(auto i = 0; i < m_width; ++i){
            for(auto j = 0; j < m_height; ++j){
                auto [r, g, b, a] = rgbaFunc(m_grids[i][j]);
                v(i, m_height - j - 1) = rgb8_pixel_t(r, g, b);
            }
        }
        write_view(filename, view(img), png_tag());
        return true;
    }
#endif

private:
    size_t m_width;
    size_t m_height;
    Container m_grids;
};

class OccupancyGridMappingFactory
{  
public:
    template <typename Occupancy>
    using GridMap = OccupancyGridMap<Occupancy>;

    template <typename property_type>
    struct Product { int x, y; property_type property; double ratio; };

    template <typename property_type>
    using ProductPipe = boost::lockfree::queue<Product<property_type>, boost::lockfree::fixed_sized<false> >;

    template <typename num_type>
    struct GridCtrl
    {
        size_t threads = 1;
        num_type stride;
        Box2D<num_type> bbox;
        Point2D<num_type> ref;

        GridCtrl(const Box2D<num_type> & _bbox, num_type _stride, size_t _threads = 1)
         : GridCtrl(_bbox, _bbox[0], _stride, _threads) {}
        
        GridCtrl(const Box2D<num_type> & _bbox, const Point2D<num_type> & _ref, num_type _stride, size_t _threads = 1)
         : threads(_threads), stride(_stride), bbox(_bbox), ref(_ref) {}
    };

    template <typename num_type>
    struct GridWorkItem
    {
        size_t begin, end;
        GridCtrl<num_type> ctrl;
    };

    template <typename num_type>
    static std::pair<size_t, size_t> GetGridMapSize(const Box2D<num_type> & bbox, num_type stride)
    {
        using float_t = common::float_type<num_type>;
        size_t width = static_cast<size_t>(std::ceil(float_t(bbox.Length()) / stride));
        size_t height = static_cast<size_t>(std::ceil(float_t(bbox.Width()) / stride));
        return std::make_pair(width, height);
    }

    template <typename property_type, typename Object, typename GeomGetter, typename Occupancy, typename BlendFunc,
              typename std::enable_if<traits::is_2d_geometry_t<typename std::result_of<GeomGetter(const Object&)>::type>::value &&
                                     (traits::is_polygon_t<typename std::result_of<GeomGetter(const Object&)>::type>::value ||
                                      traits::is_box_t<typename std::result_of<GeomGetter(const Object&)>::type>::value), bool>::type = true>
    static void Map2Grid(const std::vector<property_type> & properties, const std::vector<Object> & objects, GeomGetter && getter,
                         const GridCtrl<typename std::result_of<GeomGetter(const Object&)>::type::coor_t> & ctrl, GridMap<Occupancy> & gridMap, BlendFunc && blend)
    {
        static_assert(std::is_trivial<Product<property_type> >::value, "only trivial property supported!");
        GENERIC_ASSERT(properties.size() == objects.size())
        
        using geom_t = typename std::result_of<GeomGetter(const Object&)>::type;
        using coor_t = typename geom_t::coor_t;
        std::vector<geom_t> geometries;
        geometries.reserve(objects.size());
        for(const auto & object : objects){
            geometries.push_back(getter(object));
        }

        std::vector<const geom_t *> geomPtrs;
        geomPtrs.reserve(geometries.size());
        for(const auto & geom : geometries) geomPtrs.push_back(&geom);
        Map2Grid<property_type, geom_t, Occupancy, BlendFunc>(properties, geomPtrs, ctrl, gridMap, std::forward<BlendFunc>(blend));
    }

    template <typename property_type, typename geometry_type, typename Occupancy, typename BlendFunc>
    static void Map2Grid(const std::vector<property_type> & properties, const std::vector<const geometry_type * > & geometries, const GridCtrl<typename geometry_type::coor_t> & ctrl, GridMap<Occupancy> & gridMap, BlendFunc && blend)
    {
        static_assert(std::is_trivial<Product<property_type> >::value, "only trivial property supported!");
        GENERIC_ASSERT(properties.size() == geometries.size())

        std::atomic_bool done{false};
        auto pipeline = std::make_unique<ProductPipe<property_type> >(32767);

        // For single-thread
        // Gridding<property_type, geometry_type>(properties, geometries, ctrl, *pipeline);
        // done.store(true);
        // Mapping<property_type, Occupancy, BlendFunc>(*pipeline, done, gridMap, std::forward<BlendFunc>(blend));

        std::thread gridding(&OccupancyGridMappingFactory::Gridding<property_type, geometry_type>, std::ref(properties), std::ref(geometries), std::ref(ctrl), std::ref(*pipeline));
        std::thread mapping(&OccupancyGridMappingFactory::Mapping<property_type, Occupancy, BlendFunc>, std::ref(*pipeline), std::ref(done), std::ref(gridMap), std::ref(blend));
        
        gridding.join();
        done.store(true);
        mapping.join();
    }

private:
    template <typename property_type, typename Occupancy, typename BlendFunc>
    static void Mapping(ProductPipe<property_type> & pipeline, std::atomic_bool & done, GridMap<Occupancy> & gridMap, BlendFunc && blend)
    {
        using Result = Product<property_type>;
        auto width = static_cast<int>(gridMap.Width());
        auto height = static_cast<int>(gridMap.Height());
        auto mapping = [&](const Result & res)
        {
            if(res.x < 0 || res.x >= width) return;
            if(res.y < 0 || res.y >= height) return;
            auto & origin = gridMap(res.x, res.y);
            blend(origin, res);
        };

        while(!done.load()){
            while(pipeline.consume_one(mapping)){}
        }
        pipeline.consume_all(mapping);
    }

    template <typename property_type, typename geometry_type>
    static void Gridding(const std::vector<property_type> & properties, const std::vector<const geometry_type * > & geometries, const GridCtrl<typename geometry_type::coor_t> & ctrl, ProductPipe<property_type> & pipeline)
    {
        thread::ThreadPool pool(ctrl.threads);
        size_t size = geometries.size();
        size_t blocks = pool.Threads();
        size_t blockSize = size / blocks;

        size_t begin = 0;
        using num_type = typename geometry_type::coor_t;
        for(size_t i = 0; i < blocks && blockSize > 0; ++i){
            size_t end = begin + blockSize;
            pool.Submit(std::bind(&OccupancyGridMappingFactory::GriddingImp<property_type, geometry_type>, std::ref(properties), std::ref(geometries), GridWorkItem<num_type>{begin, end, ctrl}, std::ref(pipeline)));
            begin = end;
        }
        size_t end = size;
        if(begin != end)
            pool.Submit(std::bind(&OccupancyGridMappingFactory::GriddingImp<property_type, geometry_type>, std::ref(properties), std::ref(geometries), GridWorkItem<num_type>{begin, end, ctrl}, std::ref(pipeline)));
        
        pool.Wait();
    }

    template <typename property_type, typename geometry_type,
              typename std::enable_if<traits::is_2d_geometry_t<geometry_type>::value && traits::is_polygon_t<geometry_type>::value, bool>::type = true>
    static void GriddingImp(const std::vector<property_type> & properties, const std::vector<const geometry_type * > & polygons, GridWorkItem<typename geometry_type::coor_t> workItem, ProductPipe<property_type> & pipeline)
    {
        //TriangulationGriddingImp<property_type, geometry_type>(properties, polygons, workItem, pipeline);
        TrapezoidationGriddingImp<property_type, typename geometry_type::coor_t>(properties, polygons, workItem, pipeline);
    }

    template <typename property_type, typename geometry_type,
              typename std::enable_if<traits::is_2d_geometry_t<geometry_type>::value && traits::is_box_t<geometry_type>::value, bool>::type = true>
    static void GriddingImp(const std::vector<property_type> & properties, const std::vector<const geometry_type * > & boxes, GridWorkItem<typename geometry_type::coor_t> workItem, ProductPipe<property_type> & pipeline)
    {
        using num_type = typename geometry_type::coor_t; 
        for(size_t i = workItem.begin; i < workItem.end; ++i){
            const auto & property = properties[i];
            const auto & box = boxes[i];
            if(nullptr == box) continue;

            GriddingRectangle<property_type, num_type>(property, *box, workItem.ctrl, pipeline);
        }
    }

    template <typename property_type, typename num_type>
    static void TriangulationGriddingImp(const std::vector<property_type> & properties, const std::vector<const Polygon2D<num_type> * > & polygons, GridWorkItem<num_type> workItem, ProductPipe<property_type> & pipeline)
    {
        for(size_t i = workItem.begin; i < workItem.end; ++i){
            const auto & property = properties[i];
            const auto & polygon = polygons[i];
            if(nullptr == polygon) continue;

            tri::Triangulation<Point2D<num_type> > triangulation;
            GENERIC_ASSERT(TriangulationPolygon(polygon, triangulation))

            GriddingTriangulation<property_type, num_type>(property, triangulation, workItem.ctrl, pipeline);
        }
    }

    template <typename property_type, typename num_type>
    static void TrapezoidationGriddingImp(const std::vector<property_type> & properties, const std::vector<const Polygon2D<num_type> * > & polygons, GridWorkItem<num_type> workItem, ProductPipe<property_type> & pipeline,  const Orientation2D o = Orientation2D::Vertical)
    {
        for(size_t i = workItem.begin; i < workItem.end; ++i){
            const auto & property = properties[i];
            const auto & polygon = polygons[i];
            if(nullptr == polygon) continue;

            std::vector<Polygon2D<num_type> > trapezoids;
            TrapezoidationPolygon(polygon, trapezoids, o);

            GriddingTrapezoidation<property_type, num_type>(property, trapezoids, o, workItem.ctrl, pipeline);
        }
    }

    template <typename num_type>
    static bool TrapezoidationPolygon(const Polygon2D<num_type> * polygon, std::vector<Polygon2D<num_type> > & trapezoids, const Orientation2D o)
    {
        using namespace boost::polygon;
        trapezoids.clear();
        boolean::PolygonSet2D<num_type> ps;
        ps.insert(*polygon);
        ps.get_trapezoids(trapezoids, (o == Orientation2D::Horizontal) ? orientation_2d_enum::HORIZONTAL : orientation_2d_enum::VERTICAL);
    } 

    template <typename num_type>
    static bool TriangulationPolygon(const Polygon2D<num_type> * polygon, tri::Triangulation<Point2D<num_type> > & triangulation)
    {
        using Point = Point2D<num_type>;
        using Edge = tri::IndexEdge;
        using Triangulator = tri::Triangulator2D<num_type>;
        if(nullptr == polygon) return false;

        std::list<Edge> edges;
        size_t size = polygon->Size();
        for(size_t i = 0; i < size; ++i)
            edges.push_back(Edge{i, (i + 1) % size});

        auto points = polygon->GetPoints();
        if constexpr (std::is_integral<num_type>::value){
            tri::RemoveDuplicatesAndRemapEdges(points, edges, num_type(2));
        }

        triangulation.Clear();
        Triangulator triangulator(triangulation);
        try {
            triangulator.InsertVertices(points.begin(), points.end(), [](const Point & p){ return p[0]; }, [](const Point & p){ return p[1]; });
            triangulator.InsertEdges(edges.begin(), edges.end(), [](const Edge & e){ return e.v1(); }, [](const Edge & e){ return e.v2(); });
            triangulator.EraseOuterTriangles();
        }
        catch ( ... ){
            return false;
        }

        return true;
    }

    template <typename property_type, typename num_type>
    static void GriddingTriangulation(const property_type & property, const tri::Triangulation<Point2D<num_type> > & triangulation, const GridCtrl<num_type> & ctrl, ProductPipe<property_type> & pipeline)
    {
        using Utility = tri::TriangulationUtility<Point2D<num_type> >;
        for(auto it = 0; it < triangulation.triangles.size(); ++it){
            auto triangle = Utility::GetTriangle(triangulation, it);
            GriddingTriangle<property_type, num_type>(property, triangle, ctrl, pipeline);
        }
    }

    template <typename property_type, typename num_type>
    static void GriddingTrapezoidation(const property_type & property, std::vector<Polygon2D<num_type> > & trapezoidation, const Orientation2D o, const GridCtrl<num_type> & ctrl, ProductPipe<property_type> & pipeline)
    {
        for(auto & trapezoid : trapezoidation){
            if(trapezoid.Front() == trapezoid.Back()) trapezoid.PopBack();
            GriddingTrapezoid<property_type, num_type>(property, trapezoid, o, ctrl, pipeline);
        }
    }

    template <typename property_type, typename num_type>
    static void GriddingTrapezoid(const property_type & property, const Polygon2D<num_type> & trapezoid, const Orientation2D o, const GridCtrl<num_type> & ctrl, ProductPipe<property_type> & pipeline)
    {
        std::list<Box2D<num_type> > rects;
        std::list<Triangle2D<num_type> > triangles;
        DecomposeTrapezoid(trapezoid, o, rects, triangles);

        for(const auto & rect : rects)
            GriddingRectangle<property_type, num_type>(property, rect, ctrl, pipeline);
        
        for(const auto & triangle : triangles)
            GriddingTriangle<property_type, num_type>(property, triangle, ctrl, pipeline);
    }

    template <typename num_type>
    static void DecomposeTrapezoid(const Polygon2D<num_type> & trapezoid, const Orientation2D o, std::list<Box2D<num_type> > & rects, std::list<Triangle2D<num_type> > & triangles)
    {
        rects.clear();
        triangles.clear();
        auto size = trapezoid.Size();
        if(3 == size){
            triangles.push_back({trapezoid[0], trapezoid[1], trapezoid[2]});
            return;
        }
        else if(4 == size){
            bool res;
            auto points = std::array<Point2D<num_type>, 4>{trapezoid[0], trapezoid[1], trapezoid[2], trapezoid[3]};
            auto t = toTrapezoid(points, o, res);
            GENERIC_ASSERT(res);

            DecomposeTrapezoid(t, rects, triangles);
            return;
        }
        else {
            auto polygon = trapezoid;
            Polygon2D<num_type>::Clean(polygon);
            GENERIC_ASSERT(polygon.Size() != trapezoid.Size())
            DecomposeTrapezoid(polygon, o, rects, triangles);
        }
    }

    template <typename num_type>
    static void DecomposeTrapezoid(const Trapezoid<num_type> & trapezoid, std::list<Box2D<num_type> > & rects, std::list<Triangle2D<num_type> > & triangles)
    {
        GENERIC_ASSERT(trapezoid.isValid())
        Point2D<num_type> p[4] = { trapezoid[0], trapezoid[1], trapezoid[2], trapezoid[3] };
        if(math::EQ<num_type>(trapezoid.length[0], 0)){
            triangles.push_back({p[0], p[2], p[3]});
            return;
        }
        if(math::EQ<num_type>(trapezoid.length[1], 0)){
            triangles.push_back({p[0], p[1], p[3]});
            return;
        }

        auto i = trapezoid.direction == Orientation2D::Horizontal ? 0 : 1;
        auto h = trapezoid.direction == Orientation2D::Horizontal ? true : false;
        auto d03 = p[3][i] - p[0][i];
        auto d12 = p[2][i] - p[1][i];
        if(math::isNegative(d03) && math::isNegative(d03 + trapezoid.length[1])){
            triangles.push_back({p[0], p[1], p[2]});
            triangles.push_back({p[0], p[2], p[3]});
            return;
        }
        if(math::isPositive(d03) && math::isPositive(d03 - trapezoid.length[0])){
            triangles.push_back({p[0], p[1], p[2]});
            triangles.push_back({p[0], p[2], p[3]});
            return; 
        }
        auto s03 = h ? Point2D<num_type>(d03, 0) : Point2D<num_type>(0, d03);
        auto s12 = h ? Point2D<num_type>(d12, 0) : Point2D<num_type>(0, d12);
        auto p03 = math::isPositive(d03) ? p[0] + s03 : p[3] - s03;
        auto p12 = math::isPositive(d12) ? p[2] - s12 : p[1] + s12;

        if(math::NE<num_type>(d03, 0)) triangles.push_back({p[0], p03, p[3]});
        if(math::NE<num_type>(d12, 0)) triangles.push_back({p[1], p[2], p12});
        Box2D<num_type> rect;
        rect |= p03;
        rect |= p12;
        rect |= math::isNegative(d03) ? p[0] : p[3];
        rect |= math::isNegative(d12) ? p[2] : p[1];
        rects.push_back(std::move(rect));
    }

    template <typename property_type, typename num_type>
    static void GriddingRectangle(const property_type & property, const Box2D<num_type> & rect, const GridCtrl<num_type> & ctrl, ProductPipe<property_type> & pipeline)
    {
        auto [sx, sy] = Rasterization::Rasterize(rect[0], ctrl.stride, ctrl.ref);
        auto [ex, ey] = Rasterization::Rasterize(rect[1], ctrl.stride, ctrl.ref);

        bool swapped = false;
        if(math::GT(ey - sy, ex - sx)){
            std::swap(sx, sy);
            std::swap(ex, ey);
            swapped = true;
        }

        auto boundLLStride = [&rect, &ctrl](int index, size_t coor) { return ctrl.ref[coor] + (index + 1) * ctrl.stride - rect[0][coor]; };
        auto boundURStride = [&rect, &ctrl](int index, size_t coor) { return rect[1][coor] - index * ctrl.stride - ctrl.ref[coor]; };
        auto middleStride  = [&rect](size_t coor) {return rect[1][coor] - rect[0][coor]; };
        auto gridLength = [&](int i, int lb, int ub, size_t coor)
        {
            auto len = ctrl.stride;
            if(lb == ub) len = middleStride(coor);
            else if(i == lb) len = boundLLStride(i, coor);
            else if(i == ub) len = boundURStride(i, coor);
            return len;
        };

        auto area = ctrl.stride * ctrl.stride;
        using Result = Product<property_type>;
        for(auto i = sx; i <= ex; ++i){
            auto len1 = gridLength(i, sx, ex, swapped ? 1 : 0);
            GENERIC_ASSERT(len1 >= 0)
            for(auto j = sy; j <= ey; ++j){
                auto len2 = gridLength(j, sy, ey, swapped ? 0 : 1);
                GENERIC_ASSERT(len2 >= 0)
                auto ratio = double(len1 * len2) / area;
                auto res = swapped ? Result{j, i, property, ratio} : Result{i, j, property, ratio};
                while(!pipeline.push(res));
            }
        }
    }

    template <typename property_type, typename num_type>
    static void GriddingTriangle(const property_type & property, const Triangle2D<num_type> & triangle, const GridCtrl<num_type> & ctrl, ProductPipe<property_type> & pipeline)
    {
        std::vector<std::array<int, 2> > grids;
        std::map<int, std::set<int> > orderedGrids;
        Rasterization::Rasterize(triangle, ctrl.stride, grids, ctrl.ref);
        for(const auto & grid : grids){
            if(!orderedGrids.count(grid[0]))
                orderedGrids.insert(std::make_pair(grid[0], std::set<int>{}));
            orderedGrids[grid[0]].insert(grid[1]);
        }

        auto getBox = [&ctrl](int x, int y)
        {
            Point2D<num_type> ll = ctrl.ref + Point2D<num_type>(x * ctrl.stride, y * ctrl.stride);
            Point2D<num_type> ur = ll + Point2D<num_type>(ctrl.stride, ctrl.stride);
            return Box2D<num_type>(ll, ur);
        };
        
        auto area = ctrl.stride * ctrl.stride;
        using Result = Product<property_type>;
        auto iter_x = orderedGrids.begin();
        for(; iter_x != orderedGrids.end(); ++iter_x){
            const auto & x = iter_x->first;
            const auto & orderedYs = iter_x->second;
            auto sy = *orderedYs.begin();
            auto ey = *orderedYs.rbegin();
            for(auto y = sy; y <= ey; ++y){
                double ratio(0);
                auto box = getBox(x, y);
                if(Contains(box, triangle)) ratio = double(triangle.Area()) / area;
                else if(Contains(triangle, box)) ratio = 1.0;
                else {
                    std::list<Polygon2D<num_type> > intersects;
                    boolean::Intersect(box, triangle, intersects);
                    for(const auto & intersect : intersects)
                        ratio += double(boost::polygon::area(intersect)) / area;
                }  
                while(!pipeline.push(Result{x, y, property, ratio}));
            }
        }
    }
};

//demo
template <typename num_type>
class DensityGridMapCalculator
{
    using Occupancy = float;
    using Property = const Polygon2D<num_type> *;
public:
    using DensityGridMap = OccupancyGridMap<Occupancy>;
    using Factory = OccupancyGridMappingFactory;
    using Product = Factory::Product<Property>;

    void Insert(const Polygon2D<num_type> & polygon, bool isHole = false)
    {
        if(!isHole){
            m_solids.push_back(&polygon);
            m_solidPropties.push_back(&polygon);        
        }
        else{
            m_holes.push_back(&polygon);
            m_holePropties.push_back(&polygon);
        }
    }

    void Insert(const PolygonWithHoles2D<num_type> & pwh)
    {
        Insert(pwh.outline);
        for(const auto & hole : pwh.holes)
            Insert(hole, true);
    }

    std::unique_ptr<DensityGridMap> CalculateGridMap(const Box2D<num_type> & bbox, num_type stride, size_t threads = 1)
    {
        auto [width, height] = Factory::GetGridMapSize(bbox, stride);
        auto gridMap = std::make_unique<DensityGridMap>(width, height);

        auto ctrl = typename Factory::GridCtrl<num_type>(bbox, stride, threads);
        
        auto solidBlend = [](typename DensityGridMap::ResultType & res, const Product & p) { res += p.ratio;};
        auto holeBlend =  [](typename DensityGridMap::ResultType & res, const Product & p) { res -= p.ratio;};

        Factory::Map2Grid<Property>(m_solidPropties, m_solids, ctrl, *gridMap, solidBlend);
        Factory::Map2Grid<Property>(m_holePropties, m_holes, ctrl, *gridMap, holeBlend);

        return gridMap;
    }
private:
    std::vector<Property> m_solidPropties;
    std::vector<Property> m_holePropties;
    std::vector<const Polygon2D<num_type> * > m_solids;
    std::vector<const Polygon2D<num_type> * > m_holes;
};

}//namespace geometry
}//namespace generic
#endif//GENERIC_GEOMETRY_OCCUPANCYGRIDMAP_HPP