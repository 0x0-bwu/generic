/**
 * @file TestGraph.hpp
 * @author bwu
 * @brief Unit test of classes and functions in namespace graph
 * @version 0.1
 * @date 2024-11-29
 */
#pragma once

#include "generic/test/TestCommon.hpp"
#include "generic/utils/LinearMap.hpp"
#include "generic/utils/Index.hpp"

void t_utils()
{
    using namespace generic::utils;

    // Index
    {
        struct Tag {};
        auto a = Index<Tag>(2);
        BOOST_CHECK(a == Index<Tag>(2));
        BOOST_CHECK(a < Index<Tag>(3));
        BOOST_CHECK(size_t(a) == 2);
        BOOST_CHECK(bool(a));
        a.makeInvalid();
        BOOST_CHECK(not bool(a));
    }

    // Linear Map
    {
        LinearMap<size_t, std::string> lm(3);
        BOOST_CHECK(lm.Size() == 3);

        lm.Insert(5, "test");
        BOOST_CHECK(lm.Size() == 6);

        auto s = lm[5];
        BOOST_CHECK(s == "test");
    }
}

test_suite * create_utils_test_suite()
{
    using namespace boost::unit_test;
    test_suite * utils_suite = BOOST_TEST_SUITE("s_utils");
    //
    utils_suite->add(BOOST_TEST_CASE(&t_utils));
    //
    return utils_suite;
}