/**
 * @file MathIO.hpp
 * @author bwu
 * @brief I/O functions of math
 * @version 0.1
 * @date 2022-02-22 
 */
#ifndef GENERIC_MATH_MATHIO_HPP
#define GENERIC_MATH_MATHIO_HPP
#include "LinearAlgebra.hpp"
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <iostream>
#include <fstream>
#include <complex>
namespace {
using namespace generic::math;
using namespace generic::math::la;

///@brief out stream a vector
template <typename num_type, size_t N>
inline std::ostream & operator<< (std::ostream & os, const Vector<num_type, N> & v)
{
    return os << v.Data();
}

///@brief out stream a matrix
template <typename num_type, size_t M, size_t N>
inline std::ostream & operator<< (std::ostream & os, const Matrix<num_type, M, N> & m)
{
    return os << m.Data();
}
}

namespace generic {
namespace math {
namespace la {
using namespace boost::numeric::ublas;

class MatrixIO
{
public:
    template <typename num_type, bool zero_based = true>
    static bool ReadSparseMatrixComplex(const std::string & dia, const std::string & offDia, mapped_matrix<std::complex<num_type> > & m, std::string * err = nullptr)
    {
        std::ifstream in(dia);
        if(!in.is_open()) {
            if(err) *err = "Error: fail to open: " + dia;
            return false;
        }

        size_t i, j;
        num_type real, imag;
        size_t row = m.size1();
        size_t col = m.size2();
        while(!in.eof()){
            in >> i >> real >> imag;
            if(zero_based) m(i, i) = std::complex<num_type>(real, imag);
            else m(i - 1, i - 1) = std::complex<num_type>(real, imag);
        }

        in.close();
        in.open(offDia);
        if(!in.is_open()) {
            if(err) *err = "Error: fail to open: " + dia;
            return false;
        }
        
        while(!in.eof()){
            in >> i >> j >> real >> imag;
            if(zero_based) m(i, j) = std::complex<num_type>(real, imag);
            else m(i - 1, j - 1) = std::complex<num_type>(real, imag);
        }
        in.close();
        return true;
    }

    template <typename num_type, bool zero_based = true>
    static bool ReadSparseMatrix(const std::string & dia, const std::string & offDia, coordinate_matrix<num_type> & m, std::string * err = nullptr)
    {
        std::ifstream in(dia);
        if(!in.is_open()) {
            if(err) *err = "Error: fail to open: " + dia;
            return false;
        }

        size_t i, j;
        num_type real, imag;
        size_t row = m.size1();
        size_t col = m.size2();
        while(!in.eof()){
            in >> i >> real >> imag;
            if(zero_based) m(i, i) = real;
            else m(i - 1, i - 1) = real;
        }

        in.close();
        in.open(offDia);
        if(!in.is_open()) {
            if(err) *err = "Error: fail to open: " + dia;
            return false;
        }
        
        while(!in.eof()){
            in >> i >> j >> real >> imag;
            if(zero_based) m(i, j) = real;
            else m(i - 1, j - 1) = real;
        }
        in.close();
        return true;
    }
};
}//namespace la
}//namespace math
}//namespace generic
#endif//GENERIC_MATH_MATHIO_HPP