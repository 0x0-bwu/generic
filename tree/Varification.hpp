/**
 * @file Varification.hpp
 * @author bwu
 * @brief Varification for tree build corectness and quality check
 * @version 0.1
 * @date 2022-02-22
 */
#ifndef GENERIC_TREE_VARIFICATION_HPP
#define GENERIC_TREE_VARIFICATION_HPP
#include "KdTree.hpp"
#include "BVH.hpp"
#include <string>
namespace generic{
namespace tree{

template <typename num_type, size_t K, typename primitive, typename vectorizer>
class KdTreeVarification
{
    size_t m_maxDepth = 0;
    vectorizer m_vectorizer;
    const KdTree<num_type, K> & m_tree;
    const std::vector<primitive * > & m_primitives;

    struct HeapComp{
        bool operator() (const std::pair<size_t, num_type> & i, const std::pair<size_t, num_type> & j) const
        { return i.second < j.second; }
    };

    using NearestHeap = std::priority_queue<
                        std::pair<size_t, num_type>, 
                        std::vector<std::pair<size_t, num_type> >, HeapComp >;
                            
public:
    KdTreeVarification(const KdTree<num_type, K> & tree, const std::vector<primitive * > & prims)
     : m_tree(tree), m_primitives(prims) {}

    bool VarifyTreeStructure(size_t & maxDepth)
    {
        if(m_tree.nodeCount == 0) return false;
        const typename KdTree<num_type, K>::KdNode & node = m_tree.nodes[0];
        size_t count = 0;
        m_maxDepth = 0;
        if(!VarifyNodeStructure_(node, count, 0)) return false;
        if(count != m_primitives.size()) return false;
        maxDepth = m_maxDepth;
        return true;
    }

    bool VarifyKNearest(const VectorN<num_type, K> & vec, size_t k)
    {
        NearestHeap varifyHeap;
        for(size_t i = 0; i < m_primitives.size(); ++i){
            VectorN<num_type, K> v = m_vectorizer(*m_primitives[i]);
            num_type normSquare = VectorN<num_type, K>::NormSquare(v, vec);
            auto pr = std::make_pair(i, normSquare);
            if(varifyHeap.size() < k) varifyHeap.emplace(std::move(pr));
            else {
                if(normSquare < varifyHeap.top().second){
                    varifyHeap.pop();
                    varifyHeap.emplace(std::move(pr));
                }
            }
        }

        auto kNearest = KdTreeUtility::FindKNearestPrims(m_tree, vec, k);
        
        std::vector<size_t> varify, test;
        while(!varifyHeap.empty()){
            varify.push_back(varifyHeap.top().first);
            varifyHeap.pop();
        }

        for(const auto & nearest : kNearest)
            test.push_back(nearest.first);
        
        std::sort(varify.begin(), varify.end());
        std::sort(test.begin(), test.end());

        if(varify != test) return false;
        return true;
    }

    bool VarifyRNearest(const VectorN<num_type, K> & vec, float_type<num_type> r)
    {
        std::list<std::pair<size_t, num_type> > rNearest;
        
        while(rNearest.size() == 0){
            rNearest = KdTreeUtility::FindRNearestPrims(m_tree, vec, r);
            r *= 2;
        }

        r /= 2;
        std::vector<size_t> varify, test;
        for(size_t i = 0; i < m_primitives.size(); ++i){
            VectorN<num_type, K> v = m_vectorizer(*m_primitives[i]);
            num_type normSquare = VectorN<num_type, K>::NormSquare(v, vec);
            if(normSquare <= r * r) varify.push_back(i);
        }

        for(const auto & nearest : rNearest)
            test.push_back(nearest.first);
        
        std::sort(varify.begin(), varify.end());
        std::sort(test.begin(), test.end());

        if(varify != test) return false;
        return true;  
    }

private:
    bool VarifyNodeStructure_(const typename KdTree<num_type, K>::KdNode & node, size_t & count, size_t currDepth)
    {
        count += node.primCount;
        if(currDepth > m_maxDepth) m_maxDepth = currDepth;

        size_t split = node.split;
        num_type splitVal = m_tree.vectors[node.firstPrim][split];
        for(size_t i = 0; i < node.primCount; ++i){
            VectorN<num_type, K> varifyVec = m_vectorizer(*m_primitives[m_tree.primIndices[node.firstPrim + i]]);
            VectorN<num_type, K> testVec = m_tree.vectors[node.firstPrim + i];
            if(testVec != varifyVec) {
                return false;
            }
            if(!math::EQ(varifyVec[split], splitVal)) {
                return false;
            }
        }

        if(node.hasChildL()){
            const typename KdTree<num_type, K>::KdNode & left = m_tree.nodes[node.left];
            for(size_t i = 0; i < left.primCount; ++i){
                const auto & vec = m_tree.vectors[left.firstPrim + i];
                if(vec[split] >= splitVal) {
                    return false;
                }
            }
            if(!VarifyNodeStructure_(left, count, currDepth + 1)) {
                return false;
            }
        }

        if(node.hasChildR()){
            const typename KdTree<num_type, K>::KdNode & right = m_tree.nodes[node.right];
            for(size_t i = 0; i < right.primCount; ++i){
                const auto & vec = m_tree.vectors[right.firstPrim + i];
                if(vec[split] <= splitVal) {
                    return false;
                }
            }
            if(!VarifyNodeStructure_(right, count, currDepth + 1)) {
                return false;
            }
        }

        return true;
    }
};

template<typename num_type, typename primitive, typename extent>
class BVHVarification
{
    extent m_extent;
    size_t m_maxDepth = 0;
    size_t m_boxFault = 0;
    const bvh::BVH<num_type> & m_bvh;
    const std::vector<primitive * > & m_primitives;

    public:
    BVHVarification(const bvh::BVH<num_type> & bvh, const std::vector<primitive * > & prims)
     : m_bvh(bvh), m_primitives(prims) {}

    bool VarifyTreeStructure(double & boxTolerance, size_t & maxDepth)
    {
        if(m_bvh.nodeCount == 0) return false;
        const typename bvh::BVH<num_type>::BVHNode & node = m_bvh.nodes[0];
        size_t count = 0;
        m_maxDepth = 0;
        if(!VarifyNodeStructure_(node, count, 0)) return false;
        if(count != m_primitives.size()) return false;

        boxTolerance = double(m_boxFault) / m_bvh.nodeCount;
        maxDepth = m_maxDepth;
        return true;
    }

private:
    bool VarifyNodeStructure_(const typename bvh::BVH<num_type>::BVHNode & node, size_t & count, size_t currDepth)
    {
        count += node.primitiveCount;
        if(currDepth > m_maxDepth) m_maxDepth = currDepth;

        if(!node.isLeaf()){
            const auto & left  = m_bvh.nodes[node.firstChildOrPrim];
            const auto & right = m_bvh.nodes[node.firstChildOrPrim + 1];
            if(!(node.boundary >= left.boundary)) m_boxFault++;
            if(!(node.boundary >= right.boundary)) m_boxFault++;
            if(!VarifyNodeStructure_(left, count, currDepth + 1)) return false;
            if(!VarifyNodeStructure_(right,count, currDepth + 1)) return false; 
        }
        return true;
    }
};

}//namespace tree
}//namespace generic
#endif//GENERIC_TREE_VARIFICATION_HPP