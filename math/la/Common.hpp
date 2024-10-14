/**
 * @file Common.hpp
 * @author bwu
 * @brief common define
 * @version 0.1
 * @date 2024-05-10
 */
#pragma once
#include "generic/math/MathUtility.hpp"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <cmath>

#ifdef GENERIC_BOOST_SERIALIZATION_SUPPORT
#include "generic/common/Archive.hpp"
#endif //GENERIC_BOOST_SERIALIZATION_SUPPORT
namespace generic::math {

///@brief linear algebra related classes and functions
namespace la {
using namespace generic::common;
using namespace generic::math;

template <typename Scalar, int Rows = Eigen::Dynamic>
using DenseVector = Eigen::Matrix<Scalar, Rows, 1>;

template <typename Scalar, int Rows = Eigen::Dynamic, int Cols = Eigen::Dynamic>
using DenseMatrix = Eigen::Matrix<Scalar, Rows, Cols>;

// Convenience template for using Eigen's special allocator with vectors
template<typename Scalar, int Rows = Eigen::Dynamic, int Cols = Eigen::Dynamic>
using MatrixVector = std::vector<DenseMatrix<Scalar, Rows, Cols>, Eigen::aligned_allocator<DenseMatrix<Scalar, Rows, Cols> > >;

template <typename Scalar>
using Triplets = std::vector<Eigen::Triplet<Scalar> >;

template <typename Scalar>
using SparseMatrix = Eigen::SparseMatrix<Scalar>;

template<typename PlainObjectType, int MapOptions = Eigen::AlignmentType::Aligned, typename StrideType = Eigen::Stride<0,0>>
using VectorView = Eigen::Map<PlainObjectType, MapOptions, StrideType>;

template<typename PlainObjectType, int MapOptions = Eigen::AlignmentType::Aligned, typename StrideType = Eigen::Stride<0,0>>
using MatrixView = Eigen::Map<PlainObjectType, MapOptions, StrideType>;

using PermutMatrix = Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, size_t>;

}//namespace la
}//namespace generic::math

#ifdef GENERIC_BOOST_SERIALIZATION_SUPPORT

namespace boost::serialization {

template<class Archive, typename Scalar, int Rows, int Cols, int Ops, int MaxRows, int MaxCols>
inline void save(Archive & ar,const Eigen::Matrix<Scalar, Rows, Cols, Ops, MaxRows, MaxCols> & g, const unsigned int)
{
    int rows = g.rows();
    int cols = g.cols();

    ar & make_nvp("rows", rows);
    ar & make_nvp("cols", cols);
    ar & make_nvp("data", make_array(g.data(), rows * cols));
}

template<class Archive, typename Scalar, int Rows, int Cols, int Ops, int MaxRows, int MaxCols>
inline void load(Archive & ar, Eigen::Matrix<Scalar, Rows, Cols, Ops, MaxRows, MaxCols> & g, const unsigned int)
{
    int rows, cols;
    ar & make_nvp("rows", rows);
    ar & make_nvp("cols", cols);
    g.resize(rows, cols);
    ar & make_nvp("data", make_array(g.data(), rows * cols));
}

template<class Archive, typename Scalar, int Rows, int Cols, int Ops, int MaxRows, int MaxCols>
inline void serialize(Archive & ar, Eigen::Matrix<Scalar, Rows, Cols, Ops, MaxRows, MaxCols> & g, const unsigned int version)
{
    split_free(ar, g, version);
}

template <class Archive, typename Scalar>
inline void save(Archive & ar, const Eigen::Triplet<Scalar> & m, const unsigned int)
{
    ar & make_nvp("row", m.row());
    ar & make_nvp("col", m.col());
    ar & make_nvp("value", m.value());
}

template <class Archive, typename Scalar>
inline void load(Archive & ar, Eigen::Triplet<Scalar> & m, const unsigned int)
{
    int row, col;
    Scalar value;
    ar & make_nvp("row", row);
    ar & make_nvp("col", col);
    ar & make_nvp("value", value);
    m = Eigen::Triplet<Scalar>(row, col, value);
}

template <class Archive, typename Scalar>
inline void serialize(Archive & ar, Eigen::Triplet<Scalar> & m, const unsigned int version)
{
    split_free(ar, m, version);
}

template <class Archive, typename Scalar, int Ops,typename Index>
inline void save(Archive & ar, const Eigen::SparseMatrix<Scalar, Ops, Index> & m, const unsigned int)
{
    int innerSize = m.innerSize();
    int outerSize = m.outerSize();
    std::vector<Eigen::Triplet<Scalar>> triplets;
    triplets.reserve(m.nonZeros());
    for (int i = 0; i < outerSize; ++i) {
        using Iter = typename Eigen::SparseMatrix<Scalar, Ops, Index>::InnerIterator;
        for (Iter it(m,i); it; ++it) {
            triplets.emplace_back(it.row(), it.col(), it.value());
        }
    }
    ar & make_nvp("inner_size", innerSize);
    ar & make_nvp("outer_size", outerSize);
    ar & make_nvp("triplets", triplets);
}

template <class Archive, typename Scalar, int Ops, typename Index>
inline void load(Archive & ar, Eigen::SparseMatrix<Scalar, Ops, Index> & m, const unsigned int)
{
    int innerSize;
    int outerSize;
    ar & make_nvp("inner_size", innerSize);
    ar & make_nvp("outer_size", outerSize);
    int rows = m.IsRowMajor ? outerSize : innerSize;
    int cols = m.IsRowMajor ? innerSize : outerSize;
    m.resize(rows,cols);
    std::vector<Eigen::Triplet<Scalar>> triplets;
    ar & make_nvp("triplets", triplets);
    m.setFromTriplets(triplets.begin(), triplets.end());
}

template <class Archive, typename Scalar, int Ops, typename Index>
inline void serialize(Archive & ar, Eigen::SparseMatrix<Scalar, Ops, Index> & m, const unsigned int version)
{
    split_free(ar,m,version);
}

} // namespace boost::serialization

#endif //GENERIC_BOOST_SERIALIZATION_SUPPORT