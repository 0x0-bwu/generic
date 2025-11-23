/**
 * @file TestCommon.hpp
 * @author bwu
 * @brief Unit test common header files
 * @version 0.1
 * @date 2022-02-22
 */
#pragma once
#define BOOST_TEST_INCLUDED
// #define BOOST_TEST_DYN_LINK
#include "generic/tools/FileSystem.hpp"
#include <boost/test/included/unit_test.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/test_tools.hpp>
#include <boost/mpl/list.hpp>
#include <functional>
#include <iostream>
#include <numeric>
#include <vector>
#include <list>

inline std::string GetTestOutDataPath()
{
    return generic::fs::DirName(__FILE__).string() + "/data/out";
}

inline std::string GetTestInDataPath()
{
    return generic::fs::DirName(__FILE__).string() + "/data/in";
}