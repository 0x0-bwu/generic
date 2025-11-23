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
#include "generic/utils/Version.hpp"
#include "generic/utils/ZipView.hpp"
#include "generic/utils/Index.hpp"

void t_utils_index()
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
}

void t_utils_linear_map()
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
}

void t_utils_version()
{
    using namespace generic::utils;
    // Version
    {
        BOOST_CHECK(Version(1, 2, 3) == Version(102003));
        BOOST_CHECK(Version(99, 98, 999) < Version(9999999));
        BOOST_CHECK(Version(911126).toString() == "9.11.126");
        BOOST_CHECK(Version(911126).toInt() == 911126);
        BOOST_CHECK(std::hash<Version>()(Version(42)) != std::hash<Version>()(Version(41)));
    }
}

void t_utils_zip_view()
{
    using namespace generic::utils;
    // ZipView
    {
        std::vector<int> test1 = {1, 2, 3}, verify1;
        std::vector<char> test2 = {'a', 'b', 'c'}, verify2;
        std::vector<float> test3 = {10.1f, 20.2f, 30.3f}, verify3;

        auto zipped = Zip(test1, test2, test3);
        for (auto [t1, t2, t3] : zipped) {
            static_assert(std::is_same_v<decltype(t1), int&>);
            static_assert(std::is_same_v<decltype(t2), char&>);
            static_assert(std::is_same_v<decltype(t3), float&>);
            verify1.push_back(t1);
            verify2.push_back(t2);
            verify3.push_back(t3);
        }
        BOOST_TEST(test1 == verify1, boost::test_tools::per_element());
        BOOST_TEST(test2 == verify2, boost::test_tools::per_element());
        BOOST_TEST(test3 == verify3, boost::test_tools::per_element());
    }
}

test_suite * create_utils_test_suite()
{
    using namespace boost::unit_test;
    test_suite * utils_suite = BOOST_TEST_SUITE("s_utils");
    //
    utils_suite->add(BOOST_TEST_CASE(&t_utils_index));
    utils_suite->add(BOOST_TEST_CASE(&t_utils_linear_map));
    utils_suite->add(BOOST_TEST_CASE(&t_utils_version));
    utils_suite->add(BOOST_TEST_CASE(&t_utils_zip_view));

    //
    return utils_suite;
}