/**
 * @file TestGraph.hpp
 * @author bwu
 * @brief Unit test of classes and functions in namespace graph
 * @version 0.1
 * @date 2024-11-29
 */
#pragma once

#include "generic/test/TestCommon.hpp"
#include "generic/graph/model/DirectedGraph.hpp"
#include "generic/graph/traits/BoostGraphTraits.hpp"

test_suite * create_graph_test_suite()
{
    using namespace boost::unit_test;
    test_suite * graph_suite = BOOST_TEST_SUITE("s_graph");
    //
    // graph_suite->add(BOOST_TEST_CASE(&));
    //
    return graph_suite;
}