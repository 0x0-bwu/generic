/**
 * @file IO.hpp
 * @author bwu
 * @brief I/O function for trees
 * @version 0.1
 * @date 2022-02-22
 */
#ifndef GENERIC_TREE_IO_HPP
#define GENERIC_TREE_IO_HPP
#include "KdTree.hpp"
#include <iostream>
#include <string>
namespace generic{
namespace tree {
inline std::string toString(const kdtree::PlaneSplitMethod & method)
{
    std::string str;
    using namespace kdtree;
    if(PlaneSplitMethod::Sequential == method) str = "Sequential";
    else if(PlaneSplitMethod::MaxRange == method) str = "MaxRange";
    else if(PlaneSplitMethod::MaxVariance == method) str = "MaxVariance";
    return str;
}

inline std::string toString(const kdtree::ValueSplitMethod & method)
{
    std::string str;
    using namespace kdtree;
    if(ValueSplitMethod::Random == method) str = "Random";
    else if(ValueSplitMethod::Median == method) str = "Median";
    return str;
}
}//namespace tree
}//namespace generic
#endif//GENERIC_TREE_IO_HPP