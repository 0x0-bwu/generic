/**
 * @file MathIO.hpp
 * @author bwu
 * @brief I/O functions of math
 * @version 0.1
 * @date 2022-02-22 
 */
#pragma once
#include "generic/tools/FileSystem.hpp"
#include "LinearAlgebra.hpp"
#include <iostream>
#include <fstream>
#include <complex>

#if BOOST_GIL_IO_PNG_SUPPORT
#include "generic/tools/Color.hpp"
#include <boost/gil/extension/io/png.hpp>
#include <boost/gil/extension/numeric/sampler.hpp>
#include <boost/gil/extension/numeric/resample.hpp>
#include <boost/gil.hpp>
#include <png.h>
#endif

namespace generic::math {

namespace io {

#if BOOST_GIL_IO_PNG_SUPPORT
template <typename num_type>
inline static bool PatternView(const la::DenseMatrix<num_type> & m, const std::string & filename, size_t width = 512)
{
    auto dir = generic::fs::DirName(filename);
    if (not generic::fs::CreateDir(dir)) return false;

    num_type min = m.minCoeff();
    num_type range = m.maxCoeff() - min;
    using namespace boost::gil;
    gray8_image_t img(m.cols(), m.rows());
    gray8_image_t::view_t v = view(img);
    for (Eigen::Index row = 0; row < m.rows(); ++row){
        for (Eigen::Index col = 0; col < m.cols(); ++col){
            v(col, row) = 255 - 255 * (m(row, col) - min) / range;
        }
    }
    size_t height = double(width) / m.cols() * m.rows();
    gray8_image_t out(width, height);
    resize_view(v, view(out), bilinear_sampler());
    write_view(filename, v, png_tag());
    return true;
}

template <typename num_type>
inline static bool PatternView(const la::SparseMatrix<num_type> & m, const std::string & filename, size_t width = 512, bool maxMode = false)
{
    auto dir = generic::fs::DirName(filename);
    if (not generic::fs::CreateDir(dir)) return false;

    num_type min = m.coeffs().minCoeff();
    num_type range = m.coeffs().maxCoeff() - min;
    
    auto rows = m.rows();
    auto cols = m.cols();
    double padding = cols / double(width);
    size_t height = double(width) / cols * rows;
    
    using namespace boost::gil;
    gray8_image_t img(width, height);
    gray8_image_t::view_t v = view(img);
    for (Eigen::Index k = 0; k < m.outerSize(); ++k) {
        for (typename Eigen::SparseMatrix<num_type>::InnerIterator it(m,k); it; ++it) {
            auto val = 255 - 255 * (it.value() - min) / range;
            auto x = std::min<size_t>(it.col() / padding, width - 1);
            auto y = std::min<size_t>(it.row() / padding, height - 1);
            v(x, y) = maxMode ? std::max<uint8_t>(v(x, y), val) : (v(x, y) > 0 ? 0.5 * (v(x, y) + val) : val);
        }
    }
    write_view(filename, v, png_tag());
    return true;
}
#endif

} // namespace io
} // namespace generic::math