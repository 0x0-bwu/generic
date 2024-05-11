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
namespace generic{
namespace math{
///@brief linear algebra related classes and functions
namespace la {
using namespace generic::common;
using namespace generic::math;

template <typename num_type, int rows = Eigen::Dynamic>
using DenseVector = Eigen::Matrix<num_type, rows, 1>;

template <typename num_type, int rows = Eigen::Dynamic, int cols = Eigen::Dynamic>
using DenseMatrix = Eigen::Matrix<num_type, rows, cols>;

// Convenience template for using Eigen's special allocator with vectors
template<typename num_type, int rows = Eigen::Dynamic, int cols = Eigen::Dynamic>
using MatrixVector = std::vector<DenseMatrix<num_type, rows, cols>, Eigen::aligned_allocator<DenseMatrix<num_type, rows, cols> > >;

template <typename num_type>
using Triplets = std::vector<Eigen::Triplet<num_type> >;

template <typename num_type>
using SparseMatrix = Eigen::SparseMatrix<num_type>;

using PermutMatrix = Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, size_t>;

}//namespace la
}//namespace math
}//namespace generic