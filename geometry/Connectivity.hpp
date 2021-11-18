#ifndef GENERIC_GEOMETRY_CONNECTIVITY_HPP
#define GENERIC_GEOMETRY_CONNECTIVITY_HPP
#include "BoostPolygonRegister.hpp"
#include "generic/topology/IndexGraph.hpp"
#include "generic/thread/ThreadPool.hpp"
#include "GeometryTraits.hpp"
#include <boost/variant.hpp>
#include <unordered_map>
#include <vector>
#include <list>
#include <set>
namespace generic  {
namespace geometry {

using namespace topology;
using GeomConnGraph = SparseIndexGraph;

/**
 * @brief extract objects's physical connection
 * 
 * @tparam T object type, could be object or pointer.
 * @tparam GeomGetter a functor type to convert T to internal geometry type
 * @note the functor should have operator() (const & T) function
 * @note the functor return type should be object or reference in one of the Triangle2D, Box2D, Polygon2D, PolygonWithHoles2D
 * @param[in] objects a vector of objects for connectivity extraction
 * @param[in] geomGetter the functor object to convert T to internal geometry type
 * @param[out] connection return the connection index graph in adjacent list, the index is same as objects's vector index
 */
template <typename T, typename GeomGetter,
          typename std::enable_if<traits::is_2d_surf_geom_t<
          typename std::remove_const<typename std::remove_reference<
          typename std::result_of<GeomGetter(const T&)>::type>::type>::type>::value, bool>::type = true>
inline void ConnectivityExtraction(const std::vector<T> & objects,  GeomGetter && geomGetter, std::vector<std::set<int> > & connection)
{
    using coor_t = typename std::decay<typename std::result_of<GeomGetter(const T&)>::type>::type::coor_t;
    
    connection.clear();
    size_t size = objects.size();
    if(0 == size) return;

    connection.resize(size);
    boost::polygon::connectivity_extraction<coor_t> ce;
    for(size_t i = 0; i < size; ++i)
        GENERIC_ASSERT(i == ce.insert(geomGetter(objects[i])))

    ce.extract(connection);
}

/// A function class for layer geometry connectivity extraction
template <typename num_type>
class ConnectivityExtractor
{
    using Geometry = boost::variant<Triangle2D<num_type>, Box2D<num_type>, Polygon2D<num_type>, PolygonWithHoles2D<num_type> >;
    using GeometryContainer = std::vector<Geometry>;
    using GeomIndexContainer = std::vector<index_t>;
    using LayerGeometriesMap = std::unordered_map<index_t, GeomIndexContainer>;
    using LayerConnections = UndirectedIndexEdgeSet;
    using GeomConnections = std::unordered_map<index_t, std::unordered_set<index_t> >;
public:
    virtual ~ConnectivityExtractor() = default;

    /**
     * @brief add object to connectivity extractor
     * 
     * @tparam T object type, could be object, reference or pointer
     * @tparam GeomGetter a functor type to convert T to internal geometry type
     * @note the functor should have operator() (const & T) function, the return type should be one of the Triangle2D, Box2D, Polygon2D, PolygonWithHoles2D
     * @note the coor of the return type should be consitant with class template num_type
     * @param layer the index of the layer the obj belongs to
     * @param obj the geometry object for connectivity extraction
     * @param getter the functor object to convert T to internal geometry type
     * @return index_t global index of obj in the extraction result graph
     */
    template <typename T, typename GeomGetter,
              typename std::enable_if<traits::is_2d_surf_geom_t<
              typename std::result_of<GeomGetter(const T&)>::type>::value, bool>::type = true>
    index_t AddObject(index_t layer, const T & obj, GeomGetter && getter)
    {
        using coor_t = typename std::result_of<GeomGetter(const T&)>::type::coor_t;
        static_assert(std::is_same<coor_t, num_type>::value, "different coordinate type of GeomGetter return type and ConnectivityExtractor!");

        index_t index = m_geometries.size();
        m_geometries.push_back(getter(obj));
        if(!m_layerGeoms.count(layer))
            m_layerGeoms.insert(std::make_pair(layer, GeomIndexContainer()));
        m_layerGeoms[layer].push_back(index);
        return index;
    }

    /**
     * @brief add objects to connectivity extractor by iterator
     * 
     * @tparam Iterator objects container iterator
     * @param begin the begin iterator of objects container
     * @param end the end iterator of objects container
     * @return std::pair<index_t, index_t> global indices range from iterator begin to end in the extraction result graph
     */
    template <typename Iterator, typename GeomGetter,
              typename std::enable_if<traits::is_2d_surf_geom_t<
              typename std::result_of<GeomGetter(const typename Iterator::value_type&)>::type>::value, bool>::type = true>
    std::pair<index_t, index_t> AddObjects(index_t layer, Iterator begin, Iterator end, GeomGetter && getter)
    {
        using coor_t = typename std::result_of<GeomGetter(const typename Iterator::value_type&)>::type::coor_t;
        static_assert(std::is_same<coor_t, num_type>::value, "different coordinate type of GeomGetter return type and ConnectivityExtractor!");

        if(!m_layerGeoms.count(layer))
            m_layerGeoms.insert(std::make_pair(layer, GeomIndexContainer()));

        index_t index = m_geometries.size();
        auto range = std::make_pair(index, index);
        for(auto iter = begin; iter != end; ++iter){
            m_geometries.push_back(getter(*iter));
            m_layerGeoms[layer].push_back(range.second++);
        }
        range.second--;
        return range;
    }

    /**
     * @brief add connection b/w two layers
     * 
     * @param layer1 index of layer 1
     * @param layer2 index of layer 2
     */
    void AddLayerConnection(index_t layer1, index_t layer2)
    {
        m_layerConns.insert({layer1, layer2});
    }

    /**
     * @brief extract connection result of input objects.
     * 
     * @param threads threads number allowed in multi-threading extraction.
     * @return std::unique_ptr<GeomConnGraph> the extraction result, represented by a spares index graph
     * @note function return nullptr if no geometry exists.
     */
    std::unique_ptr<GeomConnGraph> Extract(size_t threads = 1)
    {
        if(0 == m_geometries.size()) return nullptr;
        auto graph = std::make_unique<GeomConnGraph>(m_geometries.size());

        std::unordered_set<index_t> unconnectedLayers;
        for(const auto & layerGeom : m_layerGeoms)
            unconnectedLayers.insert(layerGeom.first);

        std::list<std::future<std::unique_ptr<GeomConnections> > > futures;
        thread::ThreadPool pool(threads);

        //external
        for(const auto & layerConn : m_layerConns){
            const auto & layer1 = layerConn.v1();
            const auto & layer2 = layerConn.v2();
            if(m_layerGeoms.count(layer1) && m_layerGeoms.count(layer2)){
                auto connection = pool.Submit(std::bind(&ConnectivityExtractor<num_type>::ExtractLayersConnection, this, layer1, layer2));
                futures.emplace_back(std::move(connection));
                unconnectedLayers.erase(layer1);
                unconnectedLayers.erase(layer2);
            }
        }

        //internal
        for(const auto & layer : unconnectedLayers){
            auto connection = pool.Submit(std::bind(&ConnectivityExtractor<num_type>::ExtractLayerConnection, this, layer));
            futures.emplace_back(std::move(connection));
        }
        
        for(auto & future : futures){
            auto connection = future.get();
            AddConnection(graph.get(), std::move(connection));
        }
        
        return graph;
    }

    /// @brief clear former geometries and layer connections
    void Clear()
    {
        m_layerGeoms.clear();
        m_geometries.clear();
        m_layerConns.clear();
    }

private:
    void FillLayerGeometries(index_t layer, std::list<index_t> & geomIndices)
    {
        if(m_layerGeoms.count(layer)){
            const auto & geoms = m_layerGeoms[layer];
            geomIndices.insert(geomIndices.end(), geoms.begin(), geoms.end());
        }
    }
    
    void AddConnection(GeomConnGraph * graph, std::unique_ptr<GeomConnections> connection)
    {
        if(nullptr == connection) return;
        for(const auto & conn : *connection){
            index_t i = conn.first;
            for(index_t j : conn.second)
                AddEdge(i, j, *graph);
        }
    }

    std::unique_ptr<GeomConnections> ExtractLayerConnection(index_t layer)
    {
        std::list<index_t> geomIndices;
        FillLayerGeometries(layer, geomIndices);
        return ExtractGeomConnection(geomIndices);
    }

    std::unique_ptr<GeomConnections> ExtractLayersConnection(index_t layer1, index_t layer2)
    {
        std::list<index_t> geomIndices;
        FillLayerGeometries(layer1, geomIndices);
        FillLayerGeometries(layer2, geomIndices);
        return ExtractGeomConnection(geomIndices);
    }

    std::unique_ptr<GeomConnections> ExtractGeomConnection(const std::list<index_t> & geomIndices)
    {        
        index_t size = geomIndices.size();
        if(0 == size) return nullptr;

        std::vector<std::set<int> > tmp(size);
        std::unordered_map<int, index_t> idxMap;
        boost::polygon::connectivity_extraction<num_type> ce;
        auto addOneGeom = [this, &ce](index_t index)
        {
            int id = -1;
            if(auto * tri = boost::get<Triangle2D<num_type> >(&m_geometries[index]))
                id = ce.template insert<Triangle2D<num_type> >(*tri);
            else if(auto * box = boost::get<Box2D<num_type> >(&m_geometries[index]))
                id = ce.template insert<Box2D<num_type> >(*box);
            else if(auto * poly = boost::get<Polygon2D<num_type> >(&m_geometries[index]))
                id = ce.template insert<Polygon2D<num_type> >(*poly);
            else if(auto * pwh = boost::get<PolygonWithHoles2D<num_type> >(&m_geometries[index]))
                id = ce.template insert<PolygonWithHoles2D<num_type> >(*pwh);
            else GENERIC_THROW(std::runtime_error("unknow geometry type!"))

            return id;
        };

        for(const auto & index : geomIndices){
            idxMap.insert(std::make_pair(addOneGeom(index), index));
        }
        ce.extract(tmp);

        auto connection = std::make_unique<GeomConnections>();
        for(int i = 0; i < static_cast<int>(tmp.size()); ++i){
            if(tmp.at(i).empty()) continue;
            if(!(connection->count(idxMap.at(i))))
                connection->insert(std::make_pair(idxMap.at(i), std::unordered_set<index_t>()));
            auto & connectObj = connection->at(idxMap.at(i));
            for(auto j : tmp.at(i)){
                connectObj.insert(idxMap.at(j));
            }
        }
        return connection;
    }
    
private:
    LayerGeometriesMap m_layerGeoms;
    GeometryContainer  m_geometries;
    LayerConnections   m_layerConns;
};

}//namespace geometry
}//namespace generic

#endif//GENERIC_GEOMETRY_CONNECTIVITY_HPP