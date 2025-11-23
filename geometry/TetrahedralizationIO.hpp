/**
 * @file TetrahedralizationIO.hpp
 * @author bwu
 * @brief I/O functions of geometries
 * @version 0.1
 * @date 2022-02-22
 */
#pragma once
#include "generic/tools/StringHelper.hpp"
#include "generic/common/Traits.hpp"
#include "Tetrahedralization.hpp"
namespace generic {
namespace geometry {
namespace tet {

template <typename element_t, size_t start_index>
bool LoadElesFromFile(std::ifstream & in, size_t & line, std::vector<element_t> & elements, std::string * err = nullptr)
{
    elements.clear();
    
    char sp(32);
    std::string tmp;
    size_t dim, elSize;
    //header
    while(!in.eof()){
        line++;
        std::getline(in, tmp);
        if(tmp.empty()) continue;

        auto items = str::Split(tmp, sp);
        if(items.size() < 2){
            if(err) *err = "Error: column mismatch in line " + std::to_string(line);
            return false;  
        }
        elSize = std::stol(items[0]);
        dim = std::stol(items[1]);
        if(dim != std::tuple_size<element_t>::value){
            if(err) *err = "Error: dimension mismatch";
            return false;  
        }
        break;
    }

    size_t elCount(0);
    elements.resize(elSize);
    while(!in.eof() && elCount != elSize){
        line++;
        std::getline(in, tmp);
        if(tmp.empty()) continue;

        auto items = str::Split(tmp, sp);
        if(items.size() < dim + 1){
            if(err) *err = "Error: column mismatch in line " + std::to_string(line);
            return false;  
        }

        for(size_t i = 0; i < dim; ++i){
            elements[elCount][i] = std::stol(items[i + 1]) - start_index;
        }
        elCount++;
    }
    return elCount == elSize;    
}

template <size_t start_index>
void WriteEdgesToFile(std::ofstream & out, const std::vector<IndexEdge> & edges)
{
    char sp(32);
    out << edges.size() << sp << 1 << GENERIC_DEFAULT_EOL;
    size_t index = start_index;
    for(const auto & edge : edges){
        out << index++ << sp;
        out << edge.v1() + start_index << sp;
        out << edge.v2() + start_index << sp;
        out << 1 << GENERIC_DEFAULT_EOL;
    }
}

template <size_t start_index>
bool LoadEdgesFromFile(std::ifstream & in, size_t & line, std::vector<IndexEdge> & edges, std::string * err = nullptr)
{
    edges.clear();
    
    char sp(32);
    std::string tmp;
    size_t egSize;
    //header
    while(!in.eof()){
        line++;
        std::getline(in, tmp);
        if(tmp.empty()) continue;

        auto items = str::Split(tmp, sp);
        if(items.size() != 2){
            if(err) *err = "Error: column mismatch in line " + std::to_string(line);
            return false;  
        }
        egSize = std::stol(items[0]);
        break;
    }

    size_t egCount(0);
    size_t idx1, idx2;
    edges.resize(egSize);
    while(!in.eof() && egCount != egSize){
        line++;
        std::getline(in, tmp);
        if(tmp.empty()) continue;

        auto items = str::Split(tmp, sp);
        if(items.size() < 3){
            if(err) *err = "Error: column mismatch in line " + std::to_string(line);
            return false;  
        }

        idx1 = std::stol(items[1]) - start_index;
        idx2 = std::stol(items[2]) - start_index;
        edges[egCount] = IndexEdge(idx1, idx2);
        egCount++;
    }
    return egCount == egSize;
}

template <typename point_t, size_t start_index>
void WritePointsToFile(std::ofstream & out, const std::vector<point_t> & points)
{
    char sp(32);
    out << points.size() << sp << point_t::dim << sp << 0 << sp << 0 << GENERIC_DEFAULT_EOL;

    size_t index = start_index;
    for(const auto & point : points){
        out << index++;
        for(size_t i = 0; i < point_t::dim; ++i){
            out << sp << point[i];
        }
        out << GENERIC_DEFAULT_EOL;
    }
}

template <typename point_t>
bool LoadPointsFromFile(std::ifstream & in, size_t & line, std::vector<point_t> & points, std::string * err = nullptr)
{   
    points.clear();

    char sp(32);
    std::string tmp;
    size_t dim, ptSize;
    //header
    while(!in.eof()){
        line++;
        std::getline(in, tmp);
        if(tmp.empty()) continue;

        auto items = str::Split(tmp, sp);
        if(items.size() != 4){
            if(err) *err = "Error: column mismatch in line " + std::to_string(line);
            return false;  
        }
        ptSize = std::stol(items[0]);
        dim = std::stol(items[1]);
        if(dim != point_t::dim){
            if(err) *err = "Error: dimension mismatch";
            return false;  
        }
        break;
    }

    //points
    auto toCoord = [](const std::string & num)
    {
        using coor_t = typename point_t::coor_t;
/**
 * @brief Brief description of constexpr.
 * @param std::is_integral<coor_t>::value
 * @return if
 */
        if constexpr (std::is_integral<coor_t>::value){
            return static_cast<coor_t>(std::stoi(num));
        }
        else return static_cast<coor_t>(std::stod(num));
    };
    size_t ptCount(0);
    points.resize(ptSize);
    while(!in.eof() && ptCount != ptSize){
        line++;
        std::getline(in, tmp);
        if(tmp.empty()) continue;

        auto items = str::Split(tmp, sp);
        if(items.size() < (dim + 1)){
            if(err) *err = "Error: column mismatch in line " + std::to_string(line);
            return false;  
        }
        points[ptCount][0] = toCoord(items[1]);
        points[ptCount][1] = toCoord(items[2]);
        if constexpr (3 == point_t::dim)
            points[ptCount][2] = toCoord(items[3]);
        
        ptCount++;
    }
    return ptCount == ptSize;
}

template <typename point_t, size_t start_index>
inline void WriteSurfacesToFile(std::ofstream & out, const std::vector<typename PiecewiseLinearComplex<point_t>::Surface> & surfaces)
{
    char sp(32);
    out << surfaces.size() << sp << 0 << GENERIC_DEFAULT_EOL;

    size_t holeIdx = start_index;
    for(const auto & surface : surfaces){
        out << surface.faces.size() << sp;
        out << surface.holes.size() << sp;
        out << 0 << GENERIC_DEFAULT_EOL;

        for(const auto & face : surface.faces){
            out << face.size();
            for(const auto & v : face){
                out << sp << (v + start_index);
            }
            out << GENERIC_DEFAULT_EOL;
        }
        
        for(const auto & hole : surface.holes){
            out << holeIdx++;
            for(size_t i = 0; i < point_t::dim; ++i){
                out << sp << hole[i];
            }
            out << GENERIC_DEFAULT_EOL;
        }
    }
}

template <typename point_t, size_t start_index>
inline bool LoadSurfacesFromFile(std::ifstream & in, size_t & line,  std::vector<typename PiecewiseLinearComplex<point_t>::Surface> & surfaces, std::string * err = nullptr)
{
    char sp(32);
    size_t sfSize;
    std::string tmp;
    //header
    while(!in.eof()){
        line++;
        std::getline(in, tmp);
        if(tmp.empty()) continue;

        auto items = str::Split(tmp, sp);
        if(items.size() != 2){
            if(err) *err = "Error: column mismatch in line " + std::to_string(line);
            return false;  
        }
        sfSize = std::stol(items[0]);
        break;
    }

    auto toCoord = [](const std::string & num)
    {
        using coor_t = typename point_t::coor_t;
/**
 * @brief Brief description of constexpr.
 * @param std::is_integral<coor_t>::value
 * @return if
 */
        if constexpr (std::is_integral<coor_t>::value){
            return static_cast<coor_t>(std::stoi(num));
        }
        else return static_cast<coor_t>(std::stod(num));
    };
    using PolyLine = typename PiecewiseLinearComplex<point_t>::IdxPolyLine;
    //surface
    size_t sfCount(0);
    size_t dim = point_t::dim;
    [[maybe_unused]] size_t currMark, currFaces, currHoles;
    surfaces.resize(sfSize);
    while(!in.eof() && sfCount != sfSize){
        line++;
        std::getline(in, tmp);
        if(tmp.empty()) continue;

        auto items = str::Split(tmp, sp);
        if(items.size() == 0){
            if(err) *err = "Error: column mismatch in line " + std::to_string(line);
            return false;  
        }
        currMark  = 0;
        currFaces = 0;
        currHoles = 0;
        currFaces = std::stol(items[0]);
        if(items.size() > 1) currHoles = std::stol(items[1]);
        if(items.size() > 2) currMark = std::stol(items[2]);

        //faces
        size_t vSize, fCount(0);
        surfaces[sfCount].faces.resize(currFaces);
        while(!in.eof() && fCount != currFaces){
            line++;
            std::getline(in, tmp);
            if(tmp.empty()) continue;

            auto items = str::Split(tmp, sp);
            if(items.size() == 0){
                if(err) *err = "Error: column mismatch in line " + std::to_string(line);
                return false;  
            }

            PolyLine polyLine;
            vSize = std::stol(items[0]);
            if(items.size() != (vSize + 1)){
                if(err) *err = "Error: column mismatch in line " + std::to_string(line);
                return false;  
        
            }
            polyLine.resize(vSize);
            for(size_t i = 0; i < vSize; ++i)
                polyLine[i] = std::stol(items[i + 1]) - start_index;
            surfaces[sfCount].faces[fCount] = std::move(polyLine);
            fCount++;
        }

        //holes
        size_t hCount(0);
        surfaces[sfCount].holes.resize(currHoles);
        while(!in.eof() && hCount != currHoles){
            line++;
            std::getline(in, tmp);
            if(tmp.empty()) continue;

            auto items = str::Split(tmp, sp);
            if(items.size() != (dim + 1)){
                if(err) *err = "Error: column mismatch in line " + std::to_string(line);
                return false;  
            }
            surfaces[sfCount].holes[hCount][0] = toCoord(items[1]);
            surfaces[sfCount].holes[hCount][1] = toCoord(items[2]);
            if constexpr (3 == point_t::dim)
                surfaces[sfCount].holes[hCount][2] = toCoord(items[3]);

            hCount++;
        }
        sfCount++;
    }
    return sfCount == sfSize;
}

template <typename point_t, size_t start_index = 1>
inline bool WritePlcToNodeAndEdgeFiles(std::string_view filename, const PiecewiseLinearComplex<point_t> & plc, std::string * err = nullptr)
{
    auto basename = std::string(filename);
    std::ofstream nodeOut(basename + ".node");
    std::ofstream edgeOut(basename + ".edge");
    if (not nodeOut.is_open()) {
        if (err) *err = "Error: fail to open: " + basename + ".node";
        return false;
    }

    if (not edgeOut.is_open()){
        if (err) *err = "Error: fail to open: " + basename + ".edge";
        return false;
    }

    WritePointsToFile<point_t, start_index>(nodeOut, plc.points);

    std::vector<IndexEdge> edges;
    edges.reserve(plc.surfaces.size());
    for(const auto & surface : plc.surfaces){
        for(const auto & face : surface.faces){
            auto size = face.size();
            for(size_t i = 0; i < size; ++i){
                size_t j = (i + 1) % size;
                edges.emplace_back(IndexEdge(face[i], face[j]));
            }
        }
    }
    WriteEdgesToFile<start_index>(edgeOut, edges);

    nodeOut.close();
    edgeOut.close();
    return true;
}

template <typename point_t, size_t start_index = 1>
inline bool WritePlcToPolyFile(std::string_view filename, const PiecewiseLinearComplex<point_t> & plc, std::string * err = nullptr)
{
    std::ofstream out(filename);
    if(!out.is_open()){
        if(err) *err = "Error: fail to open: " + std::string(filename);
        return false;
    }

    WritePointsToFile<point_t, start_index>(out, plc.points);
    WriteSurfacesToFile<point_t, start_index>(out, plc.surfaces);

    out.close();
    return true;
}

template <typename point_t, size_t start_index = 1>
inline bool LoadPlcFromPolyFile(std::string_view filename, PiecewiseLinearComplex<point_t> & plc, std::string * err = nullptr)
{
    std::ifstream in(filename);
    if(!in.is_open()) {
        if(err) *err = "Error: fail to open: " + std::string(filename);
        return false;
    }

    plc.Clear();
    size_t line = 0;
    if(!LoadPointsFromFile(in, line, plc.points, err)) return false;
    if(!LoadSurfacesFromFile<point_t, start_index>(in, line, plc.surfaces, err)) return false;

    in.close();
    return true;
}

template <typename point_t>
bool LoadPointsFromNodeFile(std::string_view filename, std::vector<point_t> & points, std::string * err = nullptr)
{
    std::ifstream in(filename.data());
    if (not in.is_open()){
        if(err) *err = "Error: fail to open: " + std::string(filename);
        return false;
    }

    size_t line = 0;
    return LoadPointsFromFile(in, line, points, err);
} 

template <size_t start_index = 1>
bool LoadEdgesFromEdgeFile(std::string_view filename, std::vector<IndexEdge> & edges, std::string * err = nullptr)
{
    std::ifstream in(filename.data());
    if (not in.is_open()){
        if(err) *err = "Error: fail to open: " + std::string(filename);
        return false;
    }

    size_t line = 0;
    return LoadEdgesFromFile<start_index>(in, line, edges, err);
}

template <typename element_t, size_t start_index = 1>
bool LoadElementsFromEleFile(std::string_view filename, std::vector<element_t> & elements, std::string * err = nullptr)
{
    std::ifstream in(filename.data());
    if (not in.is_open()){
        if(err) *err = "Error: fail to open: " + std::string(filename);
        return false;
    }

    size_t line = 0;
    return LoadElesFromFile<element_t, start_index>(in, line, elements, err);
}

template <typename point_t>
bool WriteVtkFile(std::string_view filename, const Tetrahedralization<point_t> & t, std::string * err = nullptr)
{
    std::ofstream out(filename.data());
    if (not out.is_open()){
        if (err) *err = "Error: fail to open: " + std::string(filename);
        return false;
    }

    char sp(32);
    out << "# vtk DataFile Version 2.0" << GENERIC_DEFAULT_EOL;
    out << "Unstructured Grid" << GENERIC_DEFAULT_EOL;
    out << "ASCII" << GENERIC_DEFAULT_EOL;
    out << "DATASET UNSTRUCTURED_GRID" << GENERIC_DEFAULT_EOL;
    out << "POINTS" << sp << t.points.size() << sp << toString<typename point_t::coor_t>() << GENERIC_DEFAULT_EOL;
    for(const auto & point : t.points){
        out << point[0] << sp << point[1] << sp << point[2] << GENERIC_DEFAULT_EOL;
    }

    out << GENERIC_DEFAULT_EOL; 
    out << "CELLS" << sp << t.tetrahedrons.size() << sp << t.tetrahedrons.size() * 5 << GENERIC_DEFAULT_EOL;
    for(const auto & tetrahedron : t.tetrahedrons){
        out << '4';
        for(size_t i = 0; i < 4; ++i){
            out << sp << tetrahedron.vertices[i];
        }
        out << GENERIC_DEFAULT_EOL; 
    }
    out << GENERIC_DEFAULT_EOL;

    out << "CELL_TYPES" << sp << t.tetrahedrons.size() << GENERIC_DEFAULT_EOL;
    for([[maybe_unused]] const auto & tetrahedron : t.tetrahedrons){
        out << "10" << GENERIC_DEFAULT_EOL;
    }

    out.close();
    return true;
}

}//namespace tet
}//namespace geometry
}//namespace generic
