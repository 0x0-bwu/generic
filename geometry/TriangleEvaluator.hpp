/**
 * @file TriangleEvaluator.hpp
 * @author bwu
 * @brief Utility class for triangulation quality evaluation
 * @version 0.1
 * @date 2022-02-22
 */
#ifndef GENERIC_GEOMETRY_TRI_TRIANGLEEVALUATOR_HPP
#define GENERIC_GEOMETRY_TRI_TRIANGLEEVALUATOR_HPP
#include "generic/common/Traits.hpp"
#include "generic/math/Numbers.hpp"
#include "Triangulation.hpp"
#include "Utility.hpp"
namespace generic  {
namespace geometry {
namespace tri {

using generic::common::float_type;

template <size_t bins>
struct TriangleEvaluation
{
    size_t nodes;
    size_t elements;
    double minAngle;
    double maxAngle;
    double minEdgeLen;
    double maxEdgeLen;
    std::array<double, bins> triAngleHistogram;
    std::array<double, bins> triEdgeLenHistogram;
};

/**
 * @brief represents a class that take triangulation evaluation
 * @tparam point_t triangulation point type
 * @tparam bins histogram numbers of triangulation angle and edge distribution
 */
template <typename point_t, size_t bins = 10>
class TriangleEvaluator
{
public:
    using coor_t = typename point_t::coor_t;
    using float_t = float_type<coor_t>;
    using Utility = TriangulationUtility<point_t>;

    // static const size_t bins = 10;
    using Evaluation = TriangleEvaluation<bins>;

    explicit TriangleEvaluator(const Triangulation<point_t> & tri,
                               const TriIdxSet & skipT = {},
                               const VerIdxSet & skipV = {})
        : m_tri(tri), m_skipT(skipT), m_skipV(skipV) {}

    ///@brief generates the evaluation result
    Evaluation Report() const;

    ///@brief gets total edge size of the triangulation
    index_t EdgeSize() const
    {
        size_t size = 0;
        for(size_t it = 0; it < m_tri.triangles.size(); ++it){
            if(m_skipT.count(it)) continue;
            const auto & triangle = m_tri.triangles[it];
            for(auto neighbor : triangle.neighbors){
                neighbor == noNeighbor ? size += 2 : size += 1;
            }
        }
        return size / 2;
    }
    ///@brief gets total vertex size of the triangulation
    index_t VertexSize()    const { return m_tri.vertices.size()  - m_skipV.size(); }
    ///@brief gets total triangle size of the triangulation
    index_t TriangleSize()  const { return m_tri.triangles.size() - m_skipT.size(); }
    ///@brief gets minimum interior angle in triangulation
    float_t MinimumAngle()  const { return MinimumAngle(m_tri,  0, TriangleSize(), m_skipT); }
    ///@brief gets maximum interior angle in triangulation
    float_t MaximumAngle()  const { return MaximumAngle(m_tri,  0, TriangleSize(), m_skipT); }
    ///@brief gets minimum edge length in triangulation
    float_t MinEdgeLength() const { return MinEdgeLength(m_tri, 0, TriangleSize(), m_skipT); }
    ///@brief gets maximum edge length in triangulation
    float_t MaxEdgeLength() const { return MaxEdgeLength(m_tri, 0, TriangleSize(), m_skipT); }
    ///@brief gets average edge length of the triangulation
    float_t AveEdgeLength() const { return AveEdgeLength(m_tri, 0, TriangleSize(), m_skipT); }
    ///@brief gets average minimum interior angle of triangulation
    float_t AveMinimumAngle() const { return AveMinimumAngle(m_tri, 0, TriangleSize(), m_skipT); }
    
    ///@brief gets minimum angle distribution counts histogram with range [0, 180] in bins size
    std::array<size_t, bins> MinimumAngleHistogram() const { return MinimumAngleHistogram<bins>(m_tri, 0, TriangleSize(), m_skipT); }
    ///@brief gets edge length distribution counts histogram with length range [min, max] in bins size
    std::array<size_t, bins> EdgeLengthHistogram(float_t min, float_t max) const;

    static float_t MinimumAngle(const Triangulation<point_t> & tri, TriIdx begin, TriIdx end, const TriIdxSet & skipT = {});
    static float_t MaximumAngle(const Triangulation<point_t> & tri, TriIdx begin, TriIdx end, const TriIdxSet & skipT = {});
    static float_t MaxEdgeLength(const Triangulation<point_t> & tri, TriIdx begin, TriIdx end, const TriIdxSet & skipT = {});
    static float_t MinEdgeLength(const Triangulation<point_t> & tri, TriIdx begin, TriIdx end, const TriIdxSet & skipT = {});
    static float_t AveEdgeLength(const Triangulation<point_t> & tri, TriIdx begin, TriIdx end, const TriIdxSet & skipT = {});
    static float_t AveMinimumAngle(const Triangulation<point_t> & tri, TriIdx begin, TriIdx end, const TriIdxSet & skipT = {});

    template <size_t colums>
    static std::array<size_t, colums> MinimumAngleHistogram(const Triangulation<point_t> & tri, TriIdx begin, TriIdx end, const TriIdxSet & skipT = {});

    template <size_t colums>
    static std::array<size_t, colums> EdgeLengthHistogramDoubleCount(const Triangulation<point_t> & tri, TriIdx begin, TriIdx end, float_t min, float_t max, const TriIdxSet & skipT = {});

    static float_t MinimumAngle(const Triangulation<point_t> & tri, TriIdx it);
    static float_t MaximumAngle(const Triangulation<point_t> & tri, TriIdx it);
    static float_t MaxEdgeLength(const Triangulation<point_t> & tri, TriIdx it);
    static float_t MinEdgeLength(const Triangulation<point_t> & tri, TriIdx it);
    static float_t AveEdgeLength(const Triangulation<point_t> & tri, TriIdx it);
    static float_t EdgeLenSq(const Triangulation<point_t> & tri, TriIdx it, index_t ie);
private:
    const Triangulation<point_t> & m_tri;
    const TriIdxSet & m_skipT;
    const VerIdxSet & m_skipV;
};

template <typename point_t>
inline typename TriangleEvaluator<point_t>::Evaluation TriangleEvaluator<point_t>::Report() const
{
    std::array<double, bins> angles, edges;
    auto minEdge = MinEdgeLength();
    auto maxEdge = MaxEdgeLength();
    auto angleCounts = MinimumAngleHistogram();
    auto edgeCounts = EdgeLengthHistogram(minEdge, maxEdge);
    auto angleSize = std::accumulate(angleCounts.begin(), angleCounts.end(), 0);
    auto edgeSize = std::accumulate(edgeCounts.begin(), edgeCounts.end(), 0);
    for(size_t i = 0; i < bins; ++i){
        angles[i] = float_t(angleCounts[i]) / angleSize;
        edges[i] = float_t(edgeCounts[i]) / edgeSize;
    }

    return Evaluation{ VertexSize(), TriangleSize(), static_cast<double>(MinimumAngle()), static_cast<double>(MaximumAngle()), static_cast<double>(minEdge), static_cast<double>(maxEdge), angles, edges };
}

template <typename point_t>
inline typename TriangleEvaluator<point_t>::float_t
TriangleEvaluator<point_t>::MinimumAngle(const Triangulation<point_t> & tri, TriIdx begin, TriIdx end, const TriIdxSet & skipT)
{
    float_t min(math::pi);
    for(TriIdx i = begin; i < end; ++i){
        if(skipT.count(i)) continue;
        min = std::min(min, MinimumAngle(tri, i));
    }
    return min;
}

template <typename point_t>
inline typename TriangleEvaluator<point_t>::float_t
TriangleEvaluator<point_t>::MaximumAngle(const Triangulation<point_t> & tri, TriIdx begin, TriIdx end, const TriIdxSet & skipT)
{
    float_t max(0);
    for(TriIdx i = begin; i < end; ++i){
        if(skipT.count(i)) continue;
        max = std::max(max, MaximumAngle(tri, i));
    }
    return max; 
}

template <typename point_t>
inline typename TriangleEvaluator<point_t>::float_t
TriangleEvaluator<point_t>::MaxEdgeLength(const Triangulation<point_t> & tri, TriIdx begin, TriIdx end, const TriIdxSet & skipT)
{
    float_t max(0);
    for(TriIdx i = begin; i < end; ++i){
        if(skipT.count(i)) continue;
        max = std::max(max, MaxEdgeLength(tri, i));
    }
    return max;
}

template <typename point_t>
inline typename TriangleEvaluator<point_t>::float_t
TriangleEvaluator<point_t>::MinEdgeLength(const Triangulation<point_t> & tri, TriIdx begin, TriIdx end, const TriIdxSet & skipT)
{
    float_t min = std::numeric_limits<float_t>::max();
    for(TriIdx i = begin; i < end; ++i){
        if(skipT.count(i)) continue;
        min = std::min(min, MinEdgeLength(tri, i));
    }
    return min;
}

template <typename point_t>
inline typename TriangleEvaluator<point_t>::float_t
TriangleEvaluator<point_t>::AveEdgeLength(const Triangulation<point_t> & tri, TriIdx begin, TriIdx end, const TriIdxSet & skipT)
{
    size_t  num = 0;
    float_t sum = 0;
    for(TriIdx i = begin; i < end; ++i){
        if(skipT.count(i)) continue;
        sum += AveEdgeLength(tri, i);
        num += 1;
    }
    return sum / num;
}

template <typename point_t>
inline typename TriangleEvaluator<point_t>::float_t
TriangleEvaluator<point_t>::AveMinimumAngle(const Triangulation<point_t> & tri, TriIdx begin, TriIdx end, const TriIdxSet & skipT)
{
    size_t  num = 0;
    float_t sum = 0;
    for(TriIdx i = begin; i < end; ++i){
        if(skipT.count(i)) continue;
        sum += MinimumAngle(tri, i);
        num += 1;
    }
    return sum / num;
}

template <typename point_t>
inline std::array<size_t, TriangleEvaluator<point_t>::bins>
TriangleEvaluator<point_t>::EdgeLengthHistogram(float_t min, float_t max) const
{
   auto histogram = EdgeLengthHistogramDoubleCount<bins>(m_tri, 0, TriangleSize(), min, max, m_skipT);
   for(auto & count : histogram) count /= 2;
   return histogram;
}

template <typename point_t>
template <size_t colums>
inline  std::array<size_t, colums> TriangleEvaluator<point_t>::MinimumAngleHistogram(const Triangulation<point_t> & tri, TriIdx begin, TriIdx end, const TriIdxSet & skipT)
{
    std::array<size_t, colums> histogram {0};
    float_t coef = 3 * colums * math::pi_inv;
    auto index = [=](float_t angle) { return std::min(colums - 1, std::max(size_t(0), size_t(coef * angle))); };
    for(TriIdx i = begin; i < end; ++i){
        if(skipT.count(i)) continue;
        histogram[index(MinimumAngle(tri, i))]++;
    }
    return histogram;
}

template <typename point_t>
template <size_t colums>
inline std::array<size_t, colums> TriangleEvaluator<point_t>::EdgeLengthHistogramDoubleCount(const Triangulation<point_t> & tri, TriIdx begin, TriIdx end, float_t min, float_t max, const TriIdxSet & skipT)
{
    std::array<size_t, colums> histogram {0};
    float_t coef = float_t(colums) / (max - min);
    auto index = [=](float_t len) { return std::min(colums - 1, std::max(size_t(0), size_t(std::fma(len, coef, -min * coef)))); };
    for(TriIdx it = begin; it < end; ++it){
        if(skipT.count(it)) continue;
        const auto & triangle = tri.triangles[it];
        for(size_t i = 0; i < 3; ++i){
            float len = Utility::GetEdgeLength(tri, it, i);
            triangle.neighbors[i] == noNeighbor ? histogram[index(len)] += 2 : histogram[index(len)] += 1;
        }
    }
    return histogram;
}

template <typename point_t>
inline typename TriangleEvaluator<point_t>::float_t
TriangleEvaluator<point_t>::MinimumAngle(const Triangulation<point_t> & tri, TriIdx it)
{
    return Utility::GetMinimumAngle(tri, it);
}

template <typename point_t>
inline typename TriangleEvaluator<point_t>::float_t
TriangleEvaluator<point_t>::MaximumAngle(const Triangulation<point_t> & tri, TriIdx it)
{
    return Utility::GetMaximumAngle(tri, it);
}

template <typename point_t>
inline typename TriangleEvaluator<point_t>::float_t
TriangleEvaluator<point_t>::MaxEdgeLength(const Triangulation<point_t> & tri, TriIdx it)
{
    auto maxE = Utility::GetMaxLenEdge(tri, it);
    return Utility::GetEdgeLength(tri,maxE);
}

template <typename point_t>
inline typename TriangleEvaluator<point_t>::float_t
TriangleEvaluator<point_t>::MinEdgeLength(const Triangulation<point_t> & tri, TriIdx it)
{
    auto minE = Utility::GetMinLenEdge(tri, it);
    return Utility::GetEdgeLength(tri,minE);
}

template <typename point_t>
inline typename TriangleEvaluator<point_t>::float_t
TriangleEvaluator<point_t>::AveEdgeLength(const Triangulation<point_t> & tri, TriIdx it)
{
    float_t ave = 0;
    for(size_t i = 0; i < 3; ++i)
        ave += Utility::GetEdgeLength(tri, it, i) / 3;
    return ave;
}

}//namespace tri
}//namespace geometry
}//namespace generic
#endif//GENERIC_GEOMETRY_TRI_TRIANGLEEVALUATOR_HPP
