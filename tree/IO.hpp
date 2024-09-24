/**
 * @file IO.hpp
 * @author bwu
 * @brief I/O function for trees
 * @version 0.1
 * @date 2022-02-22
 */
#pragma once
#include "generic/tools/FileSystem.hpp"
#include "generic/tools/Format.hpp"
#include "RandomTree.hpp"
#include "KdTree.hpp"
#include <iostream>
#include <fstream>
#include <string>

template <typename Stream>
inline Stream & operator<< (Stream & os, const generic::tree::RandomTree & tree)
{
    os << "nid=[" << fmt::Fmt2Str(tree.nid, ",") << ']' << GENERIC_DEFAULT_EOL;
    os << "p=[" << fmt::Fmt2Str(tree.p, ",") << ']' << GENERIC_DEFAULT_EOL;
    os << "leaves=[" << fmt::Fmt2Str(tree.leaves, ",") << ']' << GENERIC_DEFAULT_EOL;
    return os;
}

namespace generic{

inline std::string toString(const tree::kdtree::PlaneSplitMethod & method)
{
    std::string str;
    using namespace tree::kdtree;
    if(PlaneSplitMethod::Sequential == method) str = "Sequential";
    else if(PlaneSplitMethod::MaxRange == method) str = "MaxRange";
    else if(PlaneSplitMethod::MaxVariance == method) str = "MaxVariance";
    return str;
}

inline std::string toString(const tree::kdtree::ValueSplitMethod & method)
{
    std::string str;
    using namespace tree::kdtree;
    if(ValueSplitMethod::Random == method) str = "Random";
    else if(ValueSplitMethod::Median == method) str = "Median";
    return str;
}
namespace tree {
namespace io {
inline bool WriteDot(const RandomTree & tree, std::string_view filename)
{
    fs::CreateDir(fs::DirName(filename));
    std::ofstream out(filename.data());
    if (not out.is_open()) return false;

    using namespace fmt;
    const char sp(32), eol('\n');
    out << Fmt2Str("strict digraph \"Tree %1%\" {", tree.Size()) << eol;
    out << sp << "node [label=\"\\N\"];" << eol;
    for (size_t i = 0; i < tree.Size(); ++i)
        out << sp << Fmt2Str("%1% [label=\"%2%\"];", tree.nid.at(i), tree.nid.at(i)) << eol;
    for (size_t i = 1; i < tree.Size(); ++i)
        out << sp << Fmt2Str("%1% -> %2% [label=\"%3%\"];", i, tree.p.at(i), "") << eol;

    out << '}' << eol;
    out.close();
    return true;
} 
} // namespace io
} // namespace tree
} // namespace generic