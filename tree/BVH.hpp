/**
 * @file BVH.hpp
 * @author bwu
 * @brief Model of bounding volume hierarchy tree concept and related algorithms
 * @version 0.1
 * @date 2022-02-22
 */
#pragma once
#include "generic/common/Traits.hpp"
#include "generic/geometry/Box.hpp"
#include "generic/geometry/Vector.hpp"
#include "generic/geometry/Triangle.hpp"
#include "generic/geometry/Utility.hpp"
#include "Builder.hpp"
#include <algorithm>
#include <cassert>
#include <atomic>
#include <memory>
#include <vector>
#include <stack>
#include <array>
#include <cmath>
#include <list>
namespace generic{
namespace tree{
///@brief bvh tree related
namespace bvh
{

using generic::common::float_type;
using generic::geometry::Point3D;
using generic::geometry::Vector3D;
using generic::geometry::Box3D;
using generic::geometry::Triangle3D;
using TopDownTaskSpawner = tree::TopDownTaskSpawner;
///@brief model of bounding volume hierarchy tree concept
template <typename num_type >
struct BVH
{
    struct BVHNode
    {
        Box3D<num_type > boundary;
        ///@brief leaf size in a node, 0 indicates this node is not a leaf node
        size_t primitiveCount;
        /**
         * @brief if is leaf node, get primitve by primives[BVH.primIndices[BVHNode.firstChildOrPrim]]
         * primive indices in a leaf node is from firstChildOrPrim to firstChilOrPrim + primitiveCount - 1
         * if is internal node, get node by BVH.nodes[BVHNode.firstChildOrPrim]
         */
        size_t firstChildOrPrim;

        bool operator== (const BVHNode & n) const
        {
            return boundary == n.boundary &&
                    primitiveCount == n.primitiveCount &&
                    firstChildOrPrim == n.firstChildOrPrim;
        }

        bool operator!= (const BVHNode & n) const
        {
            return !(*this == n);
        }

        bool isLeaf() const { return primitiveCount != 0; }
    };

    static size_t Sibling(size_t index) {
        assert(index != 0);
        return index % 2 == 1 ? index + 1 : index - 1;
    }

    static bool isLeftSibling(size_t index) {
        assert(index != 0);
        return index % 2 == 1;
    }

    size_t AddSubNodes()
    {
        return nodeCount.fetch_add(2);
    }
    
    std::vector<BVHNode> nodes;
    std::vector<size_t> primIndices;    
    std::atomic<size_t> nodeCount = { 0 };
};

///@brief bvh tree utility
class BVHUtility
{
public:
/**
 * @brief calculate bbox and center of input primitives that used for bvh tree building
 * 
 * @tparam primitive input primitive type 
 * @tparam num_type coordinate type of bvh
 * @tparam extent functor to calculate bbox of primitive object
 * @tparam centroid functor to calculate center of primitive object
 * @param[in] primitives input primitives
 * @param[out] bboxes calculated boxes
 * @param[out] centers calculated center points
 * @return bounding box that cover all input primitives
 */
template< typename primitive, typename num_type, typename extent, typename centroid >
static Box3D<num_type> CalculateBBoxAndCenter(const std::vector<primitive * > & primitives, 
                                        std::vector<Box3D<num_type > > & bboxes, 
                                        std::vector<Point3D<float_type<num_type> > > & centers)
{
    using primIter = typename std::vector<primitive * >::const_iterator;
    using boxIter = typename std::vector<Box3D<num_type > >::iterator;
    using ctIter = typename std::vector<Point3D<float_type<num_type> > >::iterator;

    bboxes.resize(primitives.size());
    centers.resize(primitives.size());
    return CalculateBlockBBoxAndCenter
            <num_type, primIter, boxIter, ctIter, extent, centroid>
            (primitives.begin(), primitives.end(), bboxes.begin(), centers.begin());
}

template<typename num_type, 
         typename prim_iterator, typename box_iterator, typename ct_iterator,
         typename extent, typename centroid >
static Box3D<num_type> CalculateBlockBBoxAndCenter(prim_iterator begin, prim_iterator end, 
                                                box_iterator boxBegin, ct_iterator ctBegin)
{
    extent ext;
    centroid ctd;
    Box3D<num_type > boundary;
    for(auto iter = begin; iter != end; ++iter, ++boxBegin,++ctBegin){
        *boxBegin = ext(*(*iter));
        *ctBegin = ctd(*(*iter));
        boundary |= (*boxBegin);
    }
    return boundary;
}
};

///@brief SAH bvh tree build algorithm
class SahBasedAlgorithm
{
protected:
    ~SahBasedAlgorithm() {}

    template <typename num_type>
    num_type ComputeCost(const BVH<num_type> & bvh, num_type traversalCost) const
    {
        num_type cost = 0;
        size_t nodeCount = bvh.nodeCount;
        for(size_t i = 0; i < nodeCount; ++i){
            if(bvh.nodes[i].isLeaf())
                cost += bvh.nodes[i].boundary.HalfArea() * bvh.nodes[i].primitiveCount;
            else
                cost += traversalCost * bvh.nodes[i].boundary.HalfArea();
        }
        return cost / bvh.nodes[0].boundary.HalfArea();
    }
};

///@brief bvh tree builder based on SAH algorithm
template <typename, size_t, size_t, typename> class BinnedSahBuildTask;
template <typename num_type, size_t bin_count, size_t max_leaf, typename TaskSpawner = TopDownTaskSpawner>
class BinnedSahBuilder
{
    using BuildTask = BinnedSahBuildTask<num_type, bin_count, max_leaf, TaskSpawner>;
    friend BuildTask;

    BVH<num_type> & m_bvh;
    num_type m_traversalCost = 1;
public:
    size_t maxDepth = 1024;

    ///@brief construct a builder with empty bvh data
    BinnedSahBuilder(BVH<num_type> & bvh) : m_bvh(bvh) {}

    /**
     * @brief builds the bvh tree based on inputs
     * @param[in] boundary whole boundary of the bvh
     * @param[in] bboxes input bbboxes got from primitives 
     * @param[in] centers input center points got from primitives
     */
    void Build(const Box3D<num_type> & boundary,
                const std::vector<Box3D<num_type> > & bboxes,
                const std::vector<Point3D<float_type<num_type> > > & centers)
    {
        size_t size = bboxes.size();
        assert(size > 0);
        assert(size == centers.size());

        m_bvh.nodes.resize(2 * size + 1);
        m_bvh.primIndices.resize(size);

        m_bvh.nodeCount = 1;
        m_bvh.nodes[0].boundary = boundary;

        for(size_t i = 0; i < size; ++i)
            m_bvh.primIndices[i] = i;
        
        BuildTask firstTask(*this, bboxes, centers);
        WorkItem item {0, 0, size, 0};

        {//ensure the run finished in multi-thread;
            TaskSpawner spawner;
            spawner.RunTask(firstTask, item);
        }

        m_bvh.nodes.resize(m_bvh.nodeCount);
    }
};

template <typename num_type, size_t bin_count, size_t max_leaf, typename TaskSpawner>
class BinnedSahBuildTask
{
    using float_t = float_type<num_type>;
    using Builder = BinnedSahBuilder<num_type, bin_count, max_leaf, TaskSpawner>;
    struct Bin
    {
        Box3D<num_type> bbox;
        size_t primitiveCount;
        num_type rightCost;
    };

    std::array<Bin, bin_count> m_binsPerAxis[3];

    Builder & m_builder;
    const std::vector<Box3D<num_type> > & m_boxes;
    const std::vector<Point3D<float_t> > & m_centers;

    std::pair<num_type, size_t> FindSplit(size_t axis)
    {
        auto & bins = m_binsPerAxis[axis];

        Box3D<num_type> currentBBox;
        size_t currentCount = 0;
        for(size_t i = bin_count - 1; i > 0; --i){
            currentBBox |= bins[i].bbox;
            currentCount += bins[i].primitiveCount;
            bins[i].rightCost = currentBBox.HalfArea() * currentCount;
        }

        currentBBox = Box3D<num_type>();
        currentCount = 0;

        auto bestSplit = std::pair<num_type, size_t>(std::numeric_limits<num_type>::max(), bin_count);
        for(size_t i = 0; i < bin_count -1; ++i){
            currentBBox |= bins[i].bbox;
            currentCount += bins[i].primitiveCount;
            num_type cost = currentBBox.HalfArea() * currentCount + bins[i + 1].rightCost;
            if(cost < bestSplit.first)
                bestSplit = std::make_pair(cost, i + 1);
        }
        return bestSplit;
    }

public:
    BinnedSahBuildTask(Builder & builder,
                        const std::vector<Box3D<num_type > > & boxes,
                        const std::vector<Point3D<float_t > > & centers)
    : m_builder(builder), m_boxes(boxes), m_centers(centers)
    {}

    std::list<WorkItem> Build(const WorkItem & item)
    {
        std::list<WorkItem> items;
        BVH<num_type> & bvh = m_builder.m_bvh;
        auto & node = bvh.nodes[item.nodeIndex];

        auto makeLeaf = [](typename BVH<num_type>::BVHNode & _node, size_t _begin, size_t _end)
        {
            _node.firstChildOrPrim = _begin;
            _node.primitiveCount = _end - _begin;
        };

        if(item.WorkSize() <= 1 || item.depth >= m_builder.maxDepth){
            makeLeaf(node, item.begin, item.end);
            return items;
        }

        std::vector<size_t > & primIndices = bvh.primIndices;
        std::pair<num_type, size_t > bestSplits[3];
        
        Box3D<num_type > bbox = node.boundary;
        auto center2Bin = geometry::SafeInverse(bbox.Diagonal()) * bin_count;
        auto binOffset = 
        Point3D<float_t >(-bbox[0][0], -bbox[0][1], -bbox[0][2]) * center2Bin;
        auto computeBinIndex = [=](const Point3D<float_t > & center, size_t axis)
        {
            int binIndex = std::round(std::fma(center[axis], center2Bin[axis], binOffset[axis]));
            return std::min(bin_count - 1, size_t(std::max(0, binIndex)));
        };

        for(size_t axis = 0; axis < 3; ++axis){
            for(Bin & bin : m_binsPerAxis[axis]){
                bin.bbox = Box3D<num_type>();
                bin.primitiveCount = 0; 
            }
        }

        for(size_t i = item.begin; i < item.end; ++i){
            size_t primitiveIndex = bvh.primIndices[i];
            for(size_t axis = 0; axis < 3; ++axis){
                Bin & bin = m_binsPerAxis[axis][computeBinIndex(m_centers[primitiveIndex], axis)];
                bin.primitiveCount ++;
                bin.bbox |= m_boxes[primitiveIndex];
            }
        }

        for(size_t axis = 0; axis < 3; ++axis){
            bestSplits[axis] = FindSplit(axis);
        }

        size_t bestAxis = 0;
        if(bestSplits[bestAxis].first > bestSplits[1].first) bestAxis = 1;
        if(bestSplits[bestAxis].first > bestSplits[2].first) bestAxis = 2;
        
        size_t splitIndex = bestSplits[bestAxis].second;
        
        num_type maxSplitCost = node.boundary.HalfArea() * (item.WorkSize() - m_builder.m_traversalCost);
        if(bestSplits[bestAxis].second == bin_count || bestSplits[bestAxis].first >= maxSplitCost){
            if(item.WorkSize() > max_leaf){
                bestAxis = node.boundary.LargestAxis();

                for(size_t i = 0, count = 0; i < bin_count - 1; ++i){
                    count += m_binsPerAxis[bestAxis][i].primitiveCount;
                    
                    if(count >= item.WorkSize() * 2 / 5 + 1){
                        splitIndex = i + 1;
                        break;
                    }
                }
            }
            else{
                makeLeaf(node, item.begin, item.end);
                return items;
            }
        }
        auto pred = [&](size_t i){ return computeBinIndex(m_centers[i], bestAxis) < splitIndex; };
        size_t beginRight = std::partition(primIndices.begin() + item.begin, primIndices.begin() + item.end, pred) - primIndices.begin();

        if(beginRight > item.begin && beginRight < item.end){
            size_t firstChild = bvh.AddSubNodes();

            auto & left = bvh.nodes[firstChild + 0];
            auto & right = bvh.nodes[firstChild + 1];
            node.firstChildOrPrim = firstChild;
            node.primitiveCount = 0;

            auto & bins = m_binsPerAxis[bestAxis];
            Box3D<num_type> leftBox, rightBox;
            for(size_t i = 0; i < bestSplits[bestAxis].second; ++i) leftBox |= bins[i].bbox;
            for(size_t i = splitIndex; i < bin_count; ++i) rightBox |= bins[i].bbox;

            left.boundary = leftBox;
            right.boundary = rightBox;

            items.emplace_back(WorkItem(firstChild + 0, item.begin, beginRight, item.depth + 1));
            items.emplace_back(WorkItem(firstChild + 1, beginRight, item.end, item.depth + 1));
            return items;
        }

        makeLeaf(node, item.begin, item.end);
        return items;
    }
};

///@brief utility class for collision detection of input bvh
template <typename num_type>
class BVHCollisionDetector
{
    typedef typename BVH<num_type>::BVHNode Node;
public:
    ///@brief constructs a dector with input bvh, will treat touched boxes as collide if considerTouch=true
    explicit BVHCollisionDetector(const BVH<num_type> & bvh, bool considerTouch = false)
     : m_bvh(bvh)
     , m_bConsiderTouch(considerTouch){}

    /**
     * @brief gets collision dectect results
     * 
     * @param[out] result list of pairs of collide boxes indices
     * @return wheher collison detected in this bvh 
     */
    bool CollisionDectect(std::list<std::pair<size_t, size_t> > & result)
    {
        result.clear();
        const Node & root = m_bvh.nodes[0];
        CollisionDectect_(result, root);
        return result.size();
    }

private:
    void CollisionDectect_(std::list<std::pair<size_t, size_t> > & result, const Node & node)
    {
        if(node.isLeaf()) return;
        CollisionDectect_(result, m_bvh.nodes[node.firstChildOrPrim], m_bvh.nodes[node.firstChildOrPrim + 1]);
        CollisionDectect_(result, m_bvh.nodes[node.firstChildOrPrim]);//left
        CollisionDectect_(result, m_bvh.nodes[node.firstChildOrPrim + 1]);//right
    }

    void CollisionDectect_(std::list<std::pair<size_t, size_t> > & result, const Node & a, const Node & b)
    {
        if(!Box3D<num_type>::Collision(a.boundary, b.boundary, m_bConsiderTouch)) return;

        if(a.isLeaf()){
            if(b.isLeaf()){
                size_t pa = a.firstChildOrPrim, pb = b.firstChildOrPrim;
                for(size_t i = 0; i < a.primitiveCount; ++i)
                    for(size_t j = 0;j < b.primitiveCount; ++j)
                        result.push_back(std::make_pair(m_bvh.primIndices[pa + i],
                                                        m_bvh.primIndices[pb + j]));//board-phase
            }
            else{
                CollisionDectect_(result, a, m_bvh.nodes[b.firstChildOrPrim]);//b.left
                CollisionDectect_(result, a, m_bvh.nodes[b.firstChildOrPrim + 1]);//b.right
            }
        }
        else {
            if(b.isLeaf()){
                CollisionDectect_(result, m_bvh.nodes[a.firstChildOrPrim], b);//a.left
                CollisionDectect_(result, m_bvh.nodes[a.firstChildOrPrim + 1], b);//a.right
            }
            else{
                CollisionDectect_(result, m_bvh.nodes[a.firstChildOrPrim], m_bvh.nodes[b.firstChildOrPrim]);//a.left, b.left
                CollisionDectect_(result, m_bvh.nodes[a.firstChildOrPrim], m_bvh.nodes[b.firstChildOrPrim + 1]);//a.left, b.right
                CollisionDectect_(result, m_bvh.nodes[a.firstChildOrPrim + 1], m_bvh.nodes[b.firstChildOrPrim]);//a.right, b.left
                CollisionDectect_(result, m_bvh.nodes[a.firstChildOrPrim + 1], m_bvh.nodes[b.firstChildOrPrim + 1]);//a.right, b.right
            }
        }
    }

private:
    const BVH<num_type> & m_bvh;
    bool m_bConsiderTouch;
};

///@brief a model of ray concept in 3d world
template <typename num_type>
struct Ray {
    Point3D<num_type> origin;
    Vector3D<num_type> direction;
    num_type tmin = 0;
    num_type tmax = 0;

    Ray() = default;
    ///@brief constructs a ray object with origin point and ray direction vector
    Ray(const Point3D<num_type> & _origin,
        const Vector3D<num_type>& _direction,
        num_type _tmin = num_type(0),
        num_type _tmax = std::numeric_limits<num_type>::max())
        : origin(_origin), direction(_direction), tmin(_tmin), tmax(_tmax)
    {}
};

template <typename num_type>
struct Intersection
{
    num_type t, u, v;
    num_type Distance() const { return t; }
    bool Hit(const Ray<num_type> & ray) const
    { return (ray.tmax >= t && t >= ray.tmin && u >= 0 && v >= 0 && (1 - u - v) >= 0); }
    Intersection(num_type _t, num_type _u, num_type _v)
    { t = _t; u = _u; v = _v; }
};

template<typename num_type, bool any_hit = false>
class TriangleIntersector
{
public:
    struct Result
    {
        int primIndex = -1;
        Intersection<num_type> intersection = Intersection<num_type>(std::numeric_limits<num_type>::max(), 0, 0);
        bool Hit(const Ray<num_type> & ray) const { return primIndex >= 0 && intersection.Hit(ray); }
        num_type Distance() const { return intersection.Distance(); }
    };
    
    TriangleIntersector(const BVH<num_type> & bvh, const std::vector<Triangle3D<num_type> * > & primitives)
     : m_bvh(bvh), m_primitives(primitives)
    {}

    static constexpr bool anyHit = any_hit;

    Result Intersect(size_t index, const Ray<num_type> & ray) const
    {
        if(index > m_bvh.primIndices.size()) return Result();
        else if(m_bvh.primIndices[index] > m_primitives.size()) return Result();

        auto intersection = Intersect_(*m_primitives[m_bvh.primIndices[index]], ray);
        if(intersection.Hit(ray)){
            Result res;
            res.primIndex = m_bvh.primIndices[index];
            res.intersection = std::move(intersection);
            return res;
        }
        return Result();
    }

private:
    //Implement based on [Moller Trumbore Algorithm]
    Intersection<num_type> Intersect_(const Triangle3D<num_type> & tri, const Ray<num_type> & ray) const
    {
        auto e1 = tri[1] - tri[0];
        auto e2 = tri[2] - tri[0];
        auto s = ray.origin - tri[0];
        auto s1 = geometry::CrossProduct(ray.direction, e2);
        auto s2 = geometry::CrossProduct(s, e1);
        auto coeff = 1.0 / geometry::DotProduct(s1, e1);
        num_type t = coeff * geometry::DotProduct(s2, e2);
        num_type u = coeff * geometry::DotProduct(s1, s);
        num_type v = coeff *  geometry::DotProduct(s2, ray.direction);
        return Intersection<num_type>(t, u, v);
    }

private:
    const BVH<num_type> & m_bvh;
    const std::vector<Triangle3D<num_type> * > & m_primitives;
};

template <typename num_type>
struct NodeIntersector
{
    using float_type = common::float_type<num_type>;
    std::array<int, 3> octant;
    Vector3D<float_type> scaledOrigin;
    Vector3D<float_type> inverseDirection;
    typedef typename BVH<num_type>::BVHNode Node;

    NodeIntersector(const Ray<num_type> & ray)
     : octant { std::signbit(ray.direction[0]),
                std::signbit(ray.direction[1]),
                std::signbit(ray.direction[2])}
    {
        inverseDirection = geometry::SafeInverse(ray.direction);
        scaledOrigin = - ray.origin.template Cast<float_type>() * inverseDirection;
    }

    num_type IntersectAxis(int axis, const num_type & p) const
    {
        return std::fma(p, inverseDirection[axis], scaledOrigin[axis]);
    }

    std::pair<num_type, num_type> Intersect(const Node & node, const Ray<num_type> & ray) const
    {
        Vector3D<num_type> entry, exit;
        for(int i = 0; i < 3; ++i){
            num_type coor1 = node.boundary[0][i];
            num_type coor2 = node.boundary[1][i];
            entry[i] = IntersectAxis(i, coor1);
            exit[i] = IntersectAxis(i, coor2);
            if(octant[i]) std::swap(entry[i], exit[i]);
        }
        return std::make_pair(
            std::max(std::max(entry[0], entry[1]), std::max(entry[2], ray.tmin)),
            std::min(std::min(exit[0], exit[1]), std::min(exit[2], ray.tmax)));
    }
};

template <typename num_type, typename node_intersector = NodeIntersector<num_type> >
class SingleRayTraverser
{
    using Node = typename BVH<num_type>::BVHNode;
public:
    struct Statistics{
        size_t traversalSteps = 0;
        size_t intersections = 0;
    };

    SingleRayTraverser(const BVH<num_type> & bvh) : m_bvh(bvh) {}

    template <typename primitive_intersector>
    typename primitive_intersector::Result
    Traverse(const Ray<num_type> & ray, primitive_intersector & intersector, Statistics * stat = nullptr) const
    {
        return Intersect_(ray, intersector, stat);
    }

private:
    template <typename primitive_intersector>
    typename primitive_intersector::Result 
    IntersectLeaf_(const Node & node, Ray<num_type> & ray, primitive_intersector & intersector,
                    typename primitive_intersector::Result & bestHit, Statistics * stat) const
    {
        assert(node.isLeaf());
        size_t begin = node.firstChildOrPrim;
        size_t end = begin + node.primitiveCount;
        if(stat) stat->intersections += end - begin;
        for(size_t i = begin; i < end; ++i){
            auto hit = intersector.Intersect(i, ray);
            if(hit.Hit(ray) && hit.Distance() < bestHit.Distance()){
                bestHit = hit;
                if(intersector.anyHit)
                    return bestHit;
                ray.tmax = hit.Distance();
            }
        }
        return bestHit;
    }

    template <typename primitive_intersector>
    typename primitive_intersector::Result 
    Intersect_(Ray<num_type> ray, primitive_intersector & intersector, Statistics * stat) const
    {
        auto bestHit = typename primitive_intersector::Result();
        if(m_bvh.nodes.size() && m_bvh.nodes.front().isLeaf())
            return IntersectLeaf_(m_bvh.nodes[0], ray, intersector, bestHit, stat);
        
        node_intersector nodeIntersector(ray);

        std::stack<size_t> nodeStack;
        size_t idxL = m_bvh.nodes[0].firstChildOrPrim;
        const auto * childL = &m_bvh.nodes[idxL];
        while(true) {
            if(stat) stat->traversalSteps++;

            auto * childR = &m_bvh.nodes[idxL + 1];
            auto timeL = nodeIntersector.Intersect(*childL, ray);
            auto timeR = nodeIntersector.Intersect(*childR, ray);

            if(timeL.first <= timeL.second && timeL.second > 0){
                if(childL->isLeaf()){
                    if(IntersectLeaf_(*childL, ray, intersector, bestHit, stat).Hit(ray) && intersector.anyHit)
                       break;
                    childL = nullptr;
                }
            }
            else childL = nullptr;

            if(timeR.first <= timeR.second && timeR.second > 0){
                if(childR->isLeaf()){
                    if(IntersectLeaf_(*childR, ray, intersector, bestHit, stat).Hit(ray) && intersector.anyHit)
                        break;
                    childR = nullptr;
                } 
            }
            else childR = nullptr;

            if(childL){
                if(childR){
                    if(timeL.first > timeR.first)
                        std::swap(childL, childR);
                    nodeStack.push(childR->firstChildOrPrim);
                }
                idxL = childL->firstChildOrPrim;
                childL = &m_bvh.nodes[idxL];
            }
            else if (childR) {
                idxL = childR->firstChildOrPrim;
                childL = &m_bvh.nodes[idxL];
            }
            else {
                if(nodeStack.empty()) break;
                idxL = nodeStack.top();
                childL = &m_bvh.nodes[idxL];
                nodeStack.pop();
            }
        }
            
        return bestHit;
    }   

private:
    const BVH<num_type> & m_bvh;
};
}//namespace bvh
}//namespace tree
}//namespace generic