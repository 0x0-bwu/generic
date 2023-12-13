/**
 * @file Topology.hpp
 * @author bwu
 * @brief Geometry topology relationship
 * @version 0.1
 * @date 2022-02-22
 */
#pragma once
#include "generic/topology/IndexGraph.hpp"
#include "Geometries.hpp"
#include <fstream>
#include <string>
namespace generic  {
namespace geometry {
using namespace generic::topology;

///@brief represents a class that extract vertex connection from geometries
template <typename num_type>
class GeoTopology2D
{
    using Point = Point2D<num_type>;
    using Connection = SparseIndexGraph;

public:
    ///@brief adds a polygon to the topology network
    void AddGeometry(const Polygon2D<num_type> & polygon);
    ///@brief adds a polygon with holes to the topology network
    void AddGeometry(const PolygonWithHoles2D<num_type> & pwh);

    ///@brief gets all points in the topology network
    const std::vector<Point> & GetAllPoints() const { return m_points; }
    ///@brief gets all edges in the topology network
    void GetAllEdges(std::list<std::pair<size_t, size_t> > & edges) const;

    ///@brief imports vertex/edge topology network from file
    bool Read(std::string_view filename, std::string * err = nullptr);
    ///@brief exports vertex/edge topology network to file
    bool Write(std::string_view filename, std::string * err = nullptr) const;

    ///@brief clear all points and connection in this topology network
    void Clear();
private:
    std::vector<Point> m_points;
    Connection m_connection;
};

template <typename num_type>
inline void GeoTopology2D<num_type>::AddGeometry(const Polygon2D<num_type> & polygon)
{
    size_t size = polygon.Size();
    size_t offset = m_points.size();
    m_points.resize(size + offset);
    auto iter = polygon.ConstBegin();
    for(size_t i = 0; i < size; ++i){
        size_t j = (i + 1) % size;
        m_points[i + offset] = polygon[i];
        AddEdge(i + offset, j + offset, m_connection);
    }
}

template <typename num_type>
inline void GeoTopology2D<num_type>::AddGeometry(const PolygonWithHoles2D<num_type> & pwh)
{
    AddGeometry(pwh.outline);
    for(const auto & hole : pwh.holes)
        AddGeometry(hole);
}

template <typename num_type>
inline void GeoTopology2D<num_type>::GetAllEdges(std::list<std::pair<size_t, size_t> > & edges) const
{
    edges.clear();
    auto iters = topology::Edges(m_connection);
    for(auto iter = iters.first; iter != iters.second; ++iter){
        auto edge = std::make_pair(Source(*iter, m_connection), Target(*iter, m_connection));
        edges.emplace_back(std::move(edge));
    }
}

template <typename num_type>
inline bool GeoTopology2D<num_type>::Read(std::string_view filename, std::string * err)
{
    std::ifstream in(filename.data());
    if(!in.is_open()) {
        if(err) *err = "Error: fail to open: " + std::string(filename);
        return false;
    }

    Clear();
    char temp;
    size_t size;
    while(!in.eof()){
        if(in.peek() == 'P'){
            in >> temp >> size;
            m_points.resize(size);
            for(size_t i = 0; i < size; ++i){
                in >> m_points[i][0] >> m_points[i][1];
            }
        }
        else if(in.peek() == 'E'){
            in >> temp >> size;
            size_t u, v;
            for(size_t i = 0; i < size; ++i){
                in >> u >> v;
                AddEdge(u, v, m_connection);
            }
        }
        else{
            std::string tmp;
            std::getline(in, tmp);
        }
    }

    in.close();
    return true;
}

template <typename num_type>
inline bool GeoTopology2D<num_type>::Write(std::string_view filename, std::string * err) const
{
    std::ofstream out(filename.data());
    if (not out.is_open()){
        if (err) *err = "Error: fail to open: " + std::string(filename);
        return false;
    }

    char sp(32);
    out << "P" << sp << m_points.size() << GENERIC_DEFAULT_EOL;
    for(const auto & p : m_points){
        out << p[0] << sp << p[1] << GENERIC_DEFAULT_EOL;
    }
    out << GENERIC_DEFAULT_EOL;

    std::list<std::pair<size_t, size_t> > edges;
    GetAllEdges(edges);
    
    out << "E" << sp << edges.size() << GENERIC_DEFAULT_EOL;
    for(const auto & e : edges){
        out << e.first << sp << e.second << GENERIC_DEFAULT_EOL;
    }
    out.close();
    return true;
}

template <typename num_type>
inline void GeoTopology2D<num_type>::Clear()
{
    m_points.clear();
    m_connection.clear();
}
}//namespace geometry
}//namespace generic