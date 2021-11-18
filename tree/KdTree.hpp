#ifndef GENERIC_TREE_KDTREE_HPP
#define GENERIC_TREE_KDTREE_HPP
#include "generic/math/MathUtility.hpp"
#include "generic/geometry/Vector.hpp"
#include "generic/common/Traits.hpp"
#include "Builder.hpp"
#include <unordered_set>
#include <algorithm>
#include <atomic>
#include <vector>
#include <queue>
#include <stack>
#include <list>
namespace generic{
namespace tree {

using generic::common::float_type;
using generic::geometry::VectorN;
template <typename num_type, size_t K>
struct KdTree
{
    struct KdNode
    {
        size_t split;
        size_t left;
        size_t right;
        size_t firstPrim;
        size_t primCount;
        static constexpr size_t noLeaf = std::numeric_limits<size_t>::max();

        void makeLeaf() { left = right = noLeaf; }
        bool isLeaf() const { return left == noLeaf && right == noLeaf; }
        bool hasChildL() const { return left != noLeaf; }
        bool hasChildR() const { return right != noLeaf; }
    };

    size_t AddNode()
    {
        return nodeCount.fetch_add(1);
    }

    std::vector<KdNode> nodes;
    std::vector<size_t> primIndices;
    std::atomic<size_t> nodeCount = { 0 };
    std::vector<VectorN<num_type, K> > vectors;
};

class KdTreeUtility
{
public:
    template <typename primitive, size_t dimension, typename num_type, typename vectorizer>
    static void CalculateVector(const std::vector<primitive * > & primitives, std::vector<VectorN<num_type, dimension> > & vectors)
    {
        using primIter = typename std::vector<primitive * >::const_iterator;
        using vecIter = typename std::vector<VectorN<num_type, dimension> >::iterator;
        vectors.resize(primitives.size());
        CalculateVector<num_type, dimension, primIter, vecIter, vectorizer>(primitives.begin(), primitives.end(),vectors.begin());
    }

    template<typename num_type, size_t dimension,
             typename prim_iterator, typename vec_iterator, typename vectorizer >
    static void CalculateVector(prim_iterator begin, prim_iterator end, vec_iterator vecBegin)
    {
        vectorizer maker;
        for(auto iter = begin; iter != end; ++iter, ++vecBegin){
            *vecBegin = maker(*(*iter));
        }
    }

    template <typename num_type, size_t K>//return [prim index, distance square]
    static std::list<std::pair<size_t, num_type> > FindKNearestPrims(const KdTree<num_type, K> & tree, const VectorN<num_type, K> & vec, size_t k)
    {
        std::list<std::pair<size_t, num_type> > items;
        if(tree.nodeCount == 0) return items;

        struct HeapComp{
            bool operator() (const std::pair<size_t, num_type> & i, const std::pair<size_t, num_type> & j) const
            { return i.second < j.second; }
        };

        using NearestHeap = std::priority_queue<
                            std::pair<size_t, num_type>, 
                            std::vector<std::pair<size_t, num_type> >, HeapComp >;
        
        NearestHeap nearest;
        std::stack<size_t> path;
        std::unordered_set<size_t> visited;

        auto updateNearestHeap = [&](size_t nodeIndex) mutable {
            if(visited.count(nodeIndex)) return;
            visited.insert(nodeIndex);
            path.push(nodeIndex);
            const typename KdTree<num_type, K>::KdNode & node = tree.nodes[nodeIndex];
            for(size_t i = 0; i < node.primCount; ++i){
                size_t index = node.firstPrim + i;
                const VectorN<num_type, K> & v = tree.vectors[index];
                num_type normSquare = VectorN<num_type, K>::NormSquare(vec, v);
                auto pr = std::make_pair(tree.primIndices[index], normSquare);

                if(nearest.size() < k)
                    nearest.push(pr);
                else if(normSquare < nearest.top().second){
                    nearest.pop();
                    nearest.push(pr);
                }
            }
        };
        
        size_t nodeIndex = 0;
        while(nodeIndex != KdTree<num_type, K>::KdNode::noLeaf){
            updateNearestHeap(nodeIndex);
            const typename KdTree<num_type, K>::KdNode & node = tree.nodes[nodeIndex];
            const VectorN<num_type, K> & v = tree.vectors[node.firstPrim];
            nodeIndex = vec[node.split] < v[node.split] ? node.left : node.right;
        }

        while(!path.empty()){
            size_t nodeIndex = path.top();
            path.pop();
            
            const typename KdTree<num_type, K>::KdNode & node = tree.nodes[nodeIndex];
            if(node.isLeaf()) continue;

            if(nearest.size() < k){
                if(node.hasChildL())
                    updateNearestHeap(node.left);
                if(node.hasChildR())
                    updateNearestHeap(node.right);
            }
            else{
                const auto & v = tree.vectors[node.firstPrim];
                num_type nodeSplitVal = v[node.split];
                num_type coorSplitVal = vec[node.split];
                num_type heapTopDist = nearest.top().second;
                if(coorSplitVal > nodeSplitVal){
                    if(node.hasChildR())
                        updateNearestHeap(node.right);
                    if((coorSplitVal - nodeSplitVal) < heapTopDist && node.hasChildL())
                        updateNearestHeap(node.left);
                }
                else{
                    if(node.hasChildL())
                        updateNearestHeap(node.left);
                    if((nodeSplitVal - coorSplitVal) < heapTopDist && node.hasChildR())
                        updateNearestHeap(node.right);
                }
            }
        }

        while(!nearest.empty()){
            items.emplace_back(nearest.top());
            nearest.pop();
        }

        return items;
    }

    template <typename num_type, size_t K>//return [prim index, distance square]
    static std::list<std::pair<size_t, num_type> > FindRNearestPrims(const KdTree<num_type, K> & tree, const VectorN<num_type, K> & vec, float_type<num_type> r)
    {
        using NearestList = std::list<std::pair<size_t, num_type> >;
        NearestList nearest;
        if(tree.nodeCount == 0) return nearest;

        float_type<num_type> r2 = r * r;
        std::stack<size_t> path;
        std::unordered_set<size_t> visited;

        auto updateNearestList = [&](size_t nodeIndex) mutable {
            if(visited.count(nodeIndex)) return;
            visited.insert(nodeIndex);
            path.push(nodeIndex);
            const typename KdTree<num_type, K>::KdNode & node = tree.nodes[nodeIndex];
            for(size_t i = 0; i < node.primCount; ++i){
                size_t index = node.firstPrim + i;
                const VectorN<num_type, K> & v = tree.vectors[index];
                num_type normSquare = VectorN<num_type, K>::NormSquare(vec, v);
                if(math::LE(float_type<num_type>(normSquare), r2))
                    nearest.emplace_back(std::make_pair(tree.primIndices[index], normSquare));
            }
        };
        
        size_t nodeIndex = 0;
        while(nodeIndex != KdTree<num_type, K>::KdNode::noLeaf){
            updateNearestList(nodeIndex);
            const typename KdTree<num_type, K>::KdNode & node = tree.nodes[nodeIndex];
            const VectorN<num_type, K> & v = tree.vectors[node.firstPrim];
            nodeIndex = vec[node.split] < v[node.split] ? node.left : node.right;
        }

        while(!path.empty()){
            size_t nodeIndex = path.top();
            path.pop();
            
            const typename KdTree<num_type, K>::KdNode & node = tree.nodes[nodeIndex];
            if(node.isLeaf()) continue;

            const auto & v = tree.vectors[node.firstPrim];
            num_type nodeSplitVal = v[node.split];
            num_type coorSplitVal = vec[node.split];
            float_type<num_type> dist = std::fabs(nodeSplitVal - coorSplitVal);
            if(coorSplitVal > nodeSplitVal){
                if(node.hasChildR())
                    updateNearestList(node.right);
                if(math::LE(dist, r) && node.hasChildL())
                    updateNearestList(node.left);
            }
            else{
                if(node.hasChildL())
                    updateNearestList(node.left);
                if(math::LE(dist, r) && node.hasChildR())
                    updateNearestList(node.right);
            }
        }

        return nearest;
    }

};

namespace kdtree{
using generic::geometry::VectorN;

enum class PlaneSplitMethod { Sequential = 0, MaxRange = 1, MaxVariance = 2 };
enum class ValueSplitMethod { Random = 0, Median = 1 };

template <typename, size_t, typename> class BuildTask;
template <typename num_type, size_t K, typename TaskSpawner = TopDownTaskSpawner>
class TreeBuilder
{
    using Task = BuildTask<num_type, K, TaskSpawner>;
    friend Task;
public:
    TreeBuilder(KdTree<num_type, K> & tree,
                  PlaneSplitMethod planeMethod = PlaneSplitMethod::MaxRange,
                  ValueSplitMethod valueMethod = ValueSplitMethod::Median)
     : m_kdTree(tree), m_planeMethod(planeMethod), m_valueMethod(valueMethod) {}

    void Build(const std::vector<VectorN<num_type, K> > & vectors)
    {
        size_t size = vectors.size();
        assert(size > 0);

        m_kdTree.nodes.resize(size);
        m_kdTree.primIndices.resize(size);
        m_kdTree.vectors.resize(size);

        for(size_t i = 0; i < size; ++i)
            m_kdTree.primIndices[i] = i;

        size_t split = 0;
        if(PlaneSplitMethod::MaxRange == m_planeMethod)
            split = FindMaxRangeSplitPlane(m_kdTree.primIndices, 0, size, vectors);
        else if(PlaneSplitMethod::MaxVariance == m_planeMethod)
            split = FindMaxVarianceSplitPlane(m_kdTree.primIndices, 0, size, vectors);

        size_t index = m_kdTree.AddNode();
        m_kdTree.nodes[index].split = split;
        
        Task firstTask(*this, vectors);
        WorkItem item {0, 0, size, 0};

        {//ensure the run finished in multi-thread;
            TaskSpawner spawner;
            spawner.RunTask(firstTask, item);
        }
        m_kdTree.nodes.resize(m_kdTree.nodeCount);
    }

private:
    static size_t FindMaxRangeSplitPlane(const std::vector<size_t> & primIndices, size_t begin, size_t end,
                          const std::vector<VectorN<num_type, K> > & vectors)
    {
        if((end - begin) <= 1) return 0;

        size_t currBestDim = 0;
        num_type currMaxRange = 0;
        num_type currMinVal, currMaxVal;
        for(size_t dim = 0; dim < K; ++dim){
            num_type val = vectors[primIndices[begin]][dim];
            currMinVal = currMaxVal = val;
            for(size_t i = begin; i < end; ++i){
                val = vectors[primIndices[begin]][dim];
                currMinVal = std::min(currMinVal, val);
                currMaxVal = std::max(currMaxVal, val);
            }
            if((currMaxVal - currMinVal) > currMaxRange){
                currMaxRange = currMaxVal - currMinVal;
                currBestDim = dim;
            }
        }
        return currBestDim;
    }

   static size_t FindMaxVarianceSplitPlane(const std::vector<size_t> & primIndices, size_t begin, size_t end,
                          const std::vector<VectorN<num_type, K> > & vectors)
    {
        if((end - begin) <= 1) return 0;

        size_t currBestDim = 0;
        float_type<num_type> currMaxVariance = 0;
        for(size_t dim = 0; dim < K; ++dim){
            float_type<num_type> s(0), ss(0), inv(1.0 / (end - begin));
            for(size_t i = begin; i < end; ++i){
                num_type val = vectors[primIndices[begin]][dim];                                                                                                                                                                                                                 
                s += val * inv;
                ss += val * val * inv;
            }
            float_type<num_type> var = ss - s * s;
            if(var > currMaxVariance){
                currMaxVariance = var;
                currBestDim = dim;
            }
        }
        return currBestDim;
    }

private:
    KdTree<num_type, K> & m_kdTree;
    PlaneSplitMethod m_planeMethod;
    ValueSplitMethod m_valueMethod;
};

template <typename num_type, size_t K, typename TaskSpawner>
class BuildTask
{
    using Tree = KdTree<num_type, K>;
    using Builder = TreeBuilder<num_type, K, TaskSpawner>;

    Builder & m_builder;
    const std::vector<VectorN<num_type, K> > & m_vectors;
public:
    BuildTask(Builder & builder, const std::vector<VectorN<num_type, K> > & vectors)
     : m_builder(builder), m_vectors(vectors)
    {}

    std::list<WorkItem> Build(const WorkItem & item)
    {
        std::list<WorkItem> items;
        if(item.WorkSize() == 0) return items;

        auto & kdTree = m_builder.m_kdTree;
        auto & vectors = kdTree.vectors;
        auto & primIndices = kdTree.primIndices;
        auto & node = kdTree.nodes[item.nodeIndex];

        num_type splitVal;
        size_t splitIndex = node.split;
        if(ValueSplitMethod::Median == m_builder.m_valueMethod)
            splitVal = FindMedianSplitValue(splitIndex, item.begin, item.end);
        else 
            splitVal = FindRandomSplitValue(splitIndex, item.begin, item.end);

        auto lt = [&](size_t i){ return math::LT(m_vectors[i][splitIndex], splitVal); };
        size_t beginRight = std::partition(primIndices.begin() + item.begin, primIndices.begin() + item.end, lt) - primIndices.begin();
        auto eq = [&](size_t i){ return math::EQ(m_vectors[i][splitIndex], splitVal); };
        size_t endRight = std::partition(primIndices.begin() + beginRight, primIndices.begin() + item.end, eq) - primIndices.begin();
        node.firstPrim = beginRight;
        node.primCount = endRight - beginRight;
        for(size_t i = beginRight; i < endRight; ++i)
            vectors[i] = m_vectors[primIndices[i]];

        size_t nextSplit = (splitIndex + 1) % K;
        if(item.begin < beginRight){
            size_t left = kdTree.AddNode();
            node.left = left;

            if(PlaneSplitMethod::MaxRange == m_builder.m_planeMethod)
                nextSplit = Builder::FindMaxRangeSplitPlane(primIndices, item.begin, beginRight, m_vectors);
            else if(PlaneSplitMethod::MaxVariance == m_builder.m_planeMethod)
                nextSplit = Builder::FindMaxVarianceSplitPlane(primIndices, item.begin, beginRight, m_vectors);
            auto & leftNode = kdTree.nodes[left];
            leftNode.split = nextSplit;

            items.emplace_back(WorkItem(left, item.begin, beginRight, item.depth + 1));
        }
        else { node.left = Tree::KdNode::noLeaf; }

        if(endRight < item.end){
            size_t right = kdTree.AddNode();
            node.right = right;

            if(PlaneSplitMethod::MaxRange == m_builder.m_planeMethod)
                nextSplit = Builder::FindMaxRangeSplitPlane(primIndices, endRight, item.end, m_vectors);
            else if(PlaneSplitMethod::MaxVariance == m_builder.m_planeMethod)
                nextSplit = Builder::FindMaxVarianceSplitPlane(primIndices, endRight, item.end, m_vectors);
            auto & rightNode = kdTree.nodes[right];
            rightNode.split = nextSplit;

            items.emplace_back(WorkItem(right, endRight, item.end, item.depth + 1));
        }
        else { node.right = Tree::KdNode::noLeaf; }

        return items;
    }

private:
    num_type FindRandomSplitValue(size_t splitIndex, size_t begin, size_t end)
    {
        auto & kdTree = m_builder.m_kdTree;
        auto & primIndices = kdTree.primIndices;
        auto splitVal = m_vectors[primIndices[math::Random(begin, end - 1)]][splitIndex];
        return splitVal;
    }

    num_type FindMedianSplitValue(size_t splitIndex, size_t begin, size_t end)
    {
        auto & kdTree = m_builder.m_kdTree;
        auto & primIndices = kdTree.primIndices;
        auto comp = [&](size_t i, size_t j){ return m_vectors[i][splitIndex] < m_vectors[j][splitIndex];};

        std::nth_element(primIndices.begin() + begin,
                         primIndices.begin() + (begin + end) / 2,
                         primIndices.begin() + end, comp);

        auto splitVal = m_vectors[primIndices[(begin + end) / 2]][splitIndex];
        return splitVal;
    }
};
}//namespace kdtree
}//namespace tree
}//namespace generic
#endif//GENERIC_TREE_KDTREE_HPP