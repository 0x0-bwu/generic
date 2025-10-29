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
    struct Tag {};
    using Index = generic::utils::Index<Tag>;
    
    // Test constructors
    LinearMap<Index, int> map1;
    BOOST_CHECK(map1.Empty());
    BOOST_CHECK_EQUAL(map1.Size(), 0);
    
    LinearMap<Index, int> map2(5);
    BOOST_CHECK_EQUAL(map2.Size(), 5);
    
    LinearMap<Index, int> map3(3, 42);
    BOOST_CHECK_EQUAL(map3.Size(), 3);
    BOOST_CHECK_EQUAL(map3[0], 42);
    BOOST_CHECK_EQUAL(map3[1], 42);
    BOOST_CHECK_EQUAL(map3[2], 42);
    
    // Test Append
    LinearMap<Index, int> map4;
    auto idx1 = map4.Append(10);
    auto idx2 = map4.Append(20);
    auto idx3 = map4.Append(30);
    BOOST_CHECK_EQUAL(map4.Size(), 3);
    BOOST_CHECK_EQUAL(map4[idx1], 10);
    BOOST_CHECK_EQUAL(map4[idx2], 20);
    BOOST_CHECK_EQUAL(map4[idx3], 30);
    
    // Test operator[]
    map4[idx1] = 100;
    BOOST_CHECK_EQUAL(map4[idx1], 100);
    
    // Test Insert
    LinearMap<Index, int> map5;
    map5.Insert(Index(0), 5);
    map5.Insert(Index(2), 15);
    BOOST_CHECK_EQUAL(map5.Size(), 3);
    BOOST_CHECK_EQUAL(map5[Index(0)], 5);
    BOOST_CHECK_EQUAL(map5[Index(2)], 15);
    
    // Test Resize
    map5.Resize(5, 99);
    BOOST_CHECK_EQUAL(map5.Size(), 5);
    BOOST_CHECK_EQUAL(map5[Index(3)], 99);
    BOOST_CHECK_EQUAL(map5[Index(4)], 99);
    
    // Test Contain
    BOOST_CHECK(map5.Contain(Index(0)));
    BOOST_CHECK(map5.Contain(Index(4)));
    BOOST_CHECK(!map5.Contain(Index(5)));
    
    // Test Clear
    map5.Clear();
    BOOST_CHECK(map5.Empty());
    BOOST_CHECK_EQUAL(map5.Size(), 0);
    
    // Test iterators
    LinearMap<Index, int> map6;
    map6.Append(1);
    map6.Append(2);
    map6.Append(3);
    
    int sum = 0;
    for (auto it = map6.Begin(); it != map6.End(); ++it) {
        sum += *it;
    }
    BOOST_CHECK_EQUAL(sum, 6);
    
    // Test range-based for loop
    sum = 0;
    for (int val : map6) {
        sum += val;
    }
    BOOST_CHECK_EQUAL(sum, 6);
    
    // Test Swap
    LinearMap<Index, int> map7;
    map7.Append(100);
    map6.Swap(map7);
    BOOST_CHECK_EQUAL(map6.Size(), 1);
    BOOST_CHECK_EQUAL(map6[Index(0)], 100);
    BOOST_CHECK_EQUAL(map7.Size(), 3);
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