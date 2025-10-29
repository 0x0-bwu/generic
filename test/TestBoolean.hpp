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
    using namespace generic::boolean::expr;
    Expression e;
    std::stringstream ss;
    {
        std::string err;
        std::string s = "(a & bbb[0]) ^ ((c & d) | (a & b)";
        auto res = ParseExpression(s.c_str(), e, &err);
        BOOST_CHECK(not res);
        BOOST_CHECK(err == "Unexpected end of content. expecting \")\" line 1\n");
    }
    {
        std::string s = "a_\\[999] & b ^ (c & d | a & b)";
        auto res = ParseExpression(s.c_str(), e);
        BOOST_CHECK(res);
        ss.str(""); ss << e;
        BOOST_CHECK(ss.str() == "((a_[999] & b) ^ ((c & d) | (a & b)))");
    }
    {
        // Test simple expression
        std::string s = "a & b";
        auto res = ParseExpression(s.c_str(), e);
        BOOST_CHECK(res);
        ss.str(""); ss << e;
        BOOST_CHECK(ss.str() == "(a & b)");
    }
    {
        // Test OR expression
        std::string s = "x | y | z";
        auto res = ParseExpression(s.c_str(), e);
        BOOST_CHECK(res);
        ss.str(""); ss << e;
        BOOST_CHECK(ss.str() == "((x | y) | z)");
    }
    {
        // Test XOR expression
        std::string s = "p ^ q";
        auto res = ParseExpression(s.c_str(), e);
        BOOST_CHECK(res);
    }
    {
        // Test nested parentheses
        std::string s = "((a & b) | (c & d))";
        auto res = ParseExpression(s.c_str(), e);
        BOOST_CHECK(res);
    }
}

void t_operation()
{
    using namespace generic::boolean;
    BOOST_CHECK(All(true, false, true, false) == false);
    BOOST_CHECK(Any(true, false, true, false) == true);
    
    // Test All with all true
    BOOST_CHECK(All(true, true, true) == true);
    
    // Test All with all false
    BOOST_CHECK(All(false, false, false) == false);
    
    // Test Any with all false
    BOOST_CHECK(Any(false, false, false) == false);
    
    // Test Any with all true
    BOOST_CHECK(Any(true, true, true) == true);
    
    // Test single argument
    BOOST_CHECK(All(true) == true);
    BOOST_CHECK(All(false) == false);
    BOOST_CHECK(Any(true) == true);
    BOOST_CHECK(Any(false) == false);
    
    // Test two arguments
    BOOST_CHECK(All(true, true) == true);
    BOOST_CHECK(All(true, false) == false);
    BOOST_CHECK(All(false, true) == false);
    BOOST_CHECK(All(false, false) == false);
    
    BOOST_CHECK(Any(true, true) == true);
    BOOST_CHECK(Any(true, false) == true);
    BOOST_CHECK(Any(false, true) == true);
    BOOST_CHECK(Any(false, false) == false);
    
    // Test many arguments
    BOOST_CHECK(All(true, true, true, true, true, true) == true);
    BOOST_CHECK(All(true, true, true, false, true, true) == false);
    BOOST_CHECK(Any(false, false, false, true, false, false) == true);
    BOOST_CHECK(Any(false, false, false, false, false, false) == false);
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
