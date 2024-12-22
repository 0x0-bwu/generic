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
#include "generic/boolean/Expression.hpp"

using namespace boost::unit_test;

void t_expression()
{
    // using namespace generic::boolean::expr;
    // Expression e;
    // std::stringstream ss;
    // {
    //     std::string s = "(a & bbb[0]) ^ ((c & d) | (a & b)";
    //     auto res = ParseExpression(s.c_str(), e);
    //     BOOST_CHECK(res);
    //     ss << e;
    //     std::cout << ss.str() << std::endl;
    //     BOOST_CHECK(ss.str() == "((a & b[0]) ^ ((c & d) | (a & b)))");
    // }
    
    // {
    //     std::string s = "a_[999] & b ^ (c & d | a & b)";
    //     auto res =ParseExpression(s.c_str(), e);
    //     BOOST_CHECK(res);
    //     ss.str(""); ss << e;
    //     std::cout << ss.str() << std::endl;
    //     BOOST_CHECK(ss.str() == "((a_[999] & b) ^ ((c & d) | (a & b)))");
    // }
}

void t_operation()
{
    using namespace generic::boolean;
    BOOST_CHECK(All(true, false, true, false) == false);
    BOOST_CHECK(Any(true, false, true, false) == true);
}

test_suite * create_boolean_test_suite()
{
    test_suite * boolean_suite = BOOST_TEST_SUITE("s_boolean");
    //
    boolean_suite->add(BOOST_TEST_CASE(&t_operation));
    boolean_suite->add(BOOST_TEST_CASE(&t_expression));
    //
    return boolean_suite;
}

