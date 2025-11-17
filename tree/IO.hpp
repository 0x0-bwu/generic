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
    using namespace generic;
    os << "nid=[" << fmt::Fmt2Str(tree.nid, ",") << ']' << GENERIC_DEFAULT_EOL;
    os << "p=[" << fmt::Fmt2Str(tree.p, ",") << ']' << GENERIC_DEFAULT_EOL;
    os << "leaves=[" << fmt::Fmt2Str(tree.leaves, ",") << ']' << GENERIC_DEFAULT_EOL;
    return os;
}

namespace generic{

inline std::string toString(const tree::kdtree::PlaneSplitMethod & method)
{
    using namespace tree::kdtree;
    if(PlaneSplitMethod::Sequential == method) return "Sequential";
    else if(PlaneSplitMethod::MaxRange == method) return "MaxRange";
    else if(PlaneSplitMethod::MaxVariance == method) return "MaxVariance";
    return "";
}

inline std::string toString(const tree::kdtree::ValueSplitMethod & method)
{
    using namespace tree::kdtree;
    if(ValueSplitMethod::Random == method) return "Random";
    else if(ValueSplitMethod::Median == method) return "Median";
    return "";
}
namespace tree {
namespace io {
inline bool WriteDOT(const RandomTree & tree, std::string_view filename)
{
    fs::CreateDir(fs::DirName(filename));
    std::ofstream out(filename.data());
    if (not out.is_open()) return false;

    using namespace fmt;
    constexpr char sp(32);
    out << Fmt2Str("strict digraph \"Tree %1%\" {", tree.Size()) << GENERIC_DEFAULT_EOL;
    out << sp << "node [label=\"\\N\"];" << GENERIC_DEFAULT_EOL;
    for (size_t i = 0; i < tree.Size(); ++i)
        out << sp << Fmt2Str("%1% [label=\"%2%\"];", tree.nid.at(i), tree.nid.at(i)) << GENERIC_DEFAULT_EOL;
    for (size_t i = 1; i < tree.Size(); ++i)
        out << sp << Fmt2Str("%1% -> %2% [label=\"%3%\"];", i, tree.p.at(i), "") << GENERIC_DEFAULT_EOL;

    out << '}' << GENERIC_DEFAULT_EOL;
    out.close();
    return true;
} 

template <typename num_type>
inline bool WriteVTK(const bvh::BVH<num_type> & tree, std::string_view filename)
{
    fs::CreateDir(fs::DirName(filename));
    std::ofstream out(filename.data());
    if (not out.is_open()) return false;

    constexpr char sp(32);
    out << "# vtk DataFile Version 2.0" << GENERIC_DEFAULT_EOL;
    out << "Unstructured Grid" << GENERIC_DEFAULT_EOL;
    out << "ASCII" << GENERIC_DEFAULT_EOL;
    out << "DATASET UNSTRUCTURED_GRID" << GENERIC_DEFAULT_EOL;
    out << "POINTS" << sp << tree.nodes.size() * 8 << sp << toString<num_type>() << GENERIC_DEFAULT_EOL;
    for (const auto & node : tree.nodes) {
        const auto & ll = node.boundary[0];
        const auto & ur = node.boundary[1];
        auto streamOut = [&](const auto & x, const auto & y, const auto & z) {
            out << x[0] << sp << y[1] << sp << z[2] << GENERIC_DEFAULT_EOL;
        };
        streamOut(ll, ur, ll);
        streamOut(ur, ur, ll);
        streamOut(ll, ur, ur);
        streamOut(ur, ur, ur);
        streamOut(ll, ll, ll);
        streamOut(ur, ll, ll);
        streamOut(ll, ll, ur);
        streamOut(ur, ll, ur);
    }
    out << GENERIC_DEFAULT_EOL;
    out << "CELLS" << sp << tree.nodes.size() << sp << tree.nodes.size() * 9 << GENERIC_DEFAULT_EOL;
    size_t index = 0;
    for ([[maybe_unused]] auto & node : tree.nodes) {
        out << 8;
        for (size_t i = 0; i < 8; ++i)
            out << sp << index++;
        out << GENERIC_DEFAULT_EOL;
    }

    out << GENERIC_DEFAULT_EOL;
    out << "CELL_TYPES" << sp << tree.nodes.size() << GENERIC_DEFAULT_EOL;
    for ([[maybe_unused]] auto & node : tree.nodes)
        out << "11" << GENERIC_DEFAULT_EOL;

    out << GENERIC_DEFAULT_EOL;
    out << "CELL_DATA" << sp << tree.nodes.size() << GENERIC_DEFAULT_EOL;
    out << "SCALARS SCALARS FLOAT 1" << GENERIC_DEFAULT_EOL;
    out << "LOOKUP_TABLE BBOX" << GENERIC_DEFAULT_EOL;

    out << 1 << GENERIC_DEFAULT_EOL;
    for (size_t i = 1; i < tree.nodes.size(); ++i)
        out << (bvh::BVH<num_type>::isLeftSibling(i) ? 0 : 1) << GENERIC_DEFAULT_EOL;

    out << GENERIC_DEFAULT_EOL;
    out << "LOOKUP_TABLE BBOX 3" << GENERIC_DEFAULT_EOL;
    out << "1 0 0 1" << GENERIC_DEFAULT_EOL;
    out << "0 1 0 1" << GENERIC_DEFAULT_EOL;
    out << "0 0 1 1" << GENERIC_DEFAULT_EOL;

    out.close();
    return true;
}

} // namespace io
} // namespace tree
} // namespace generic
