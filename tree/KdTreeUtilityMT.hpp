/**
 * @file KdTreeUtilityMT.hpp
 * @author bwu
 * @brief kd tree utility with multi-threads support
 * @version 0.1
 * @date 2022-02-22
 */
#ifndef GENERIC_TREE_KDTREEUTILITYMT_HPP
#define GENERIC_TREE_KDTREEUTILITYMT_HPP
#include "BuilderMT.hpp"
#include "KdTree.hpp"
namespace generic {
namespace tree {
using generic::thread::ThreadPool;
///@brief kd tree utility with multi-threads support
class KdTreeUtilityMT
{
public:
    template <typename primitive, size_t dimension, typename num_type, typename vectorizer>
    static void CalculateVector(const std::vector<primitive * > & primitives, std::vector<VectorN<num_type, dimension> > & vectors)
    {
        using primIter = typename std::vector<primitive * >::const_iterator;
        using vecIter = typename std::vector<VectorN<num_type, dimension> >::iterator;

        size_t size = primitives.size();
        if(0 == size) return;
        vectors.resize(size);

        ThreadPool pool;
        size_t blocks = pool.Threads();
        if(0 == blocks) blocks = 1;
        size_t blockSize = size / blocks;

        primIter begin = primitives.begin();
        vecIter vecBegin = vectors.begin();
        for(size_t i = 0; i < blocks && blockSize > 0; ++i){
            primIter end = begin;
            std::advance(end, blockSize);
            pool.Submit([=]{ return KdTreeUtility::CalculateVector<num_type, dimension, primIter, vecIter, vectorizer>(begin, end, vecBegin); });
            std::advance(vecBegin, blockSize);
            begin = end;
        }
        primIter end = primitives.end();
        if(begin != end)
            pool.Submit([=]{ return KdTreeUtility::CalculateVector<num_type, dimension, primIter, vecIter, vectorizer>(begin, end, vecBegin); });
        
        pool.Wait();
    }
};
}//namespace tree
}//namespace generic
#endif//GENERIC_TREE_KDTREEUTILITYMT_HPP