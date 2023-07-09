/**
 * @file BVHUtilityMT.hpp
 * @author bwu
 * @brief BVH tree utility with multi-threads support
 * @version 0.1
 * @date 2022-02-22
 */
#pragma once
#include "BuilderMT.hpp"
#include "BVH.hpp"
namespace generic{
namespace tree{
namespace bvh
{
using generic::thread::ThreadPool;
///@brief bvh tree utility with multi-threads support
class BVHUtilityMT
{
public:
template< typename primitive, typename num_type, typename extent, typename centroid >
static Box3D<num_type> CalculateBBoxAndCenter(const std::vector<primitive * > & primitives, 
                                        std::vector<Box3D<num_type > > & bboxes, 
                                        std::vector<Point3D<float_type<num_type> > > & centers)
{
    using primIter = typename std::vector<primitive * >::const_iterator;
    using boxIter = typename std::vector<Box3D<num_type > >::iterator;
    using ctIter = typename std::vector<Point3D<float_type<num_type> > >::iterator;

    size_t size = primitives.size();
    if(0 == size) return Box3D<num_type>();

    bboxes.resize(size);
    centers.resize(size);

    ThreadPool pool;
    size_t blocks = pool.Threads();
    if(0 == blocks) blocks = 1;
    size_t blockSize = size / blocks;

    primIter begin = primitives.begin();
    boxIter boxBegin = bboxes.begin();
    ctIter ctBegin = centers.begin();
    std::vector<std::future<Box3D<num_type> > > futures(blocks);
    for(size_t i = 0; i < blocks && blockSize > 0; ++i){
        primIter end = begin;
        std::advance(end, blockSize);
        futures[i] = pool.Submit([=]{
            return BVHUtility::CalculateBlockBBoxAndCenter<num_type, primIter, boxIter, ctIter, extent, centroid>(begin, end, boxBegin, ctBegin); });
        
        std::advance(boxBegin, blockSize);
        std::advance(ctBegin, blockSize);
        begin = end;
    }
    primIter end = primitives.end();
    Box3D<num_type> boundary = BVHUtility::CalculateBlockBBoxAndCenter<num_type, primIter, boxIter, ctIter, extent, centroid>(begin, end, boxBegin, ctBegin);
    for(size_t i = 0; i < blocks; ++i){
        boundary |= futures[i].get();
    }

    return boundary;
}
};
}//namespace bvh
}//namespace tree
}//namespace generic