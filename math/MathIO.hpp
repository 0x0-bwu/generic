/**
 * @file MathIO.hpp
 * @author bwu
 * @brief I/O functions of math
 * @version 0.1
 * @date 2022-02-22 
 */
#pragma once
#include "LinearAlgebra.hpp"
#include <iostream>
#include <fstream>
#include <complex>

#if BOOST_GIL_IO_PNG_SUPPORT
#include "generic/tools/Color.hpp"
#include <boost/gil/extension/io/png.hpp>
#include <boost/gil.hpp>
#include <png.h>
#endif
