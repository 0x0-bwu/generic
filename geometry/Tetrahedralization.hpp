#ifndef GENERIC_GEOMETRY_TET_TETRAHEDRALIZATION_HPP
#define GENERIC_GEOMETRY_TET_TETRAHEDRALIZATION_HPP
#include "generic/topology/IndexGraph.hpp"
#include "Point.hpp"
#include <vector>
#include <array>
#include <list>
namespace generic {
namespace geometry {
namespace tet {

using IndexEdge = topology::UndirectedIndexEdge;
using IndexEle4 = std::array<size_t,  4>;
using IndexEle10 = std::array<size_t, 10>;

template <typename point_t>
class PiecewiseLinearComplex
{
public:
    using IdxPolyLine = std::vector<size_t>;
    struct Surface
    {
        std::vector<IdxPolyLine> faces;
        std::vector<point_t> holes;
    };

    std::vector<point_t> points;
    std::vector<Surface> surfaces;
    
    void Clear();
};

template <typename point_t>
void PiecewiseLinearComplex<point_t>::Clear()
{
    points.clear();
    surfaces.clear();
}

struct IndexVertex;
struct IndexTetrahedron;
using index_t = topology::index_t;
using PosIdx = index_t;
using VerIdx = index_t;
using TetIdx = index_t;

}//namespace tet
}//namespace geometry
}//namespace generic
#endif//GENERIC_GEOMETRY_TET_TETRAHEDRALIZATION_HPP
