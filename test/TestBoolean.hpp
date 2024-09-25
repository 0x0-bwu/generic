/**
 * @file TestBoolean.hpp
 * @author bwu
 * @brief Unit test of classes and functions in namespace boolean
 * @version 0.1
 * @date 2024-09-26
 */
#pragma once
#include "generic/test/TestCommon.hpp"
#include "generic/boolean/Operation.hpp"

using namespace boost::unit_test;
using namespace generic::boolean;

void t_operation()
{
    BOOST_CHECK(All(true, false, true, false) == false);
    BOOST_CHECK(Any(true, false, true, false) == true);
}

test_suite * create_boolean_test_suite()
{
    test_suite * boolean_suite = BOOST_TEST_SUITE("s_boolean");
    //
    boolean_suite->add(BOOST_TEST_CASE(&t_operation));
    //
    return boolean_suite;
}

