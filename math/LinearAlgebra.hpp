/**
 * @file LinearAlgebra.hpp
 * @author bwu
 * @brief Model of vector and matrix concept in static dimension size, implement based on boost ublas
 * @version 0.1
 * @date 2022-02-22
 */
#pragma once
#include "generic/common/Version.hpp"
#include "generic/common/Traits.hpp"
#include "MathUtility.hpp"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <cmath>
namespace generic{
namespace math{
///@brief linear algebra related classes and functions
namespace la {
using namespace generic::common;
using namespace generic::math;

template <typename num_type>
using DenseVector = Eigen::Matrix<num_type, Eigen::Dynamic, 1>;

template <typename num_type>
using DenseMatrix = Eigen::Matrix<num_type, Eigen::Dynamic, Eigen::Dynamic>;

template <typename num_type>
using Triplets = std::vector<Eigen::Triplet<num_type> >;

template <typename num_type>
using SparseMatrix = Eigen::SparseMatrix<num_type>;

using PermutMatrix = Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, size_t>;

}//namespace la
}//namespace math
}//namespace generic