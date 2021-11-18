#ifndef TEST_TESTMATH_HPP
#define TEST_TESTMATH_HPP
#define BOOST_TEST_INCLUDED
#include <boost/test/unit_test.hpp>
#include <boost/test/test_tools.hpp>
#include <boost/mpl/list.hpp>
#include "generic/math/MathUtility.hpp"
#include "generic/math/LinearAlgebra.hpp"
#include "generic/math/MathIO.hpp"
using namespace boost::unit_test;
using namespace generic;
using namespace generic::math;
using t_math_num_types = boost::mpl::list<int, float, double, long int, long double>;

BOOST_TEST_CASE_TEMPLATE_FUNCTION(t_math_utility_t, num_type)
{
    //Random
    for(auto i = 0; i < 10; ++i){
        num_type max(999999999);
        num_type min = -max;
        num_type rand = Random(min, max);
        BOOST_CHECK_LE(rand, max);
        BOOST_CHECK_GE(rand, min);
    }
}

void t_math_utility()
{
    //Compare
    BOOST_CHECK(EQ(1, 1) == true);
    BOOST_CHECK(EQ(1e-4, 1.1e-4) == false);
    BOOST_CHECK(EQ(1e-4, 1.1e-4, 1e-5) == true);
    BOOST_CHECK(LE(1.1e-4, 1e-4) == false);
    BOOST_CHECK(LE(1.1e-4, 1e-4, 1e-5) == true);
    BOOST_CHECK(GE(1e-4, 1.1e-4) == false);
    BOOST_CHECK(GE(1e-4, 1.1e-4, 1e-5) == true);
}

BOOST_TEST_CASE_TEMPLATE_FUNCTION(t_math_linear_algebra_t, math_num_types)
{
    using namespace la;
    using num = math_num_types;
    //Vector
    Vector<num, 3> v1(1), v2(1, 1, 1), v3{1, 2, 3};
    BOOST_CHECK(v1 == v2);
    v2.Swap(v3);
    BOOST_CHECK(v1 == v3);
    BOOST_CHECK((-v1 == Vector<num, 3>(-1)));
    BOOST_CHECK(((v1 + v3) == Vector<num, 3>(2)));
    BOOST_CHECK(((v1 - v3) == Vector<num, 3>()));
    BOOST_CHECK(((v2 * 2 / 2) == v2));
    BOOST_CHECK((v1.T() == Matrix<num, 3, 1>(1)));

    //Matrix
    Matrix<num, 2, 3> m1(1), m2(v2, v3), m3{ {8, 12, 16}, {8, 8, 8} };
    BOOST_CHECK((m1 + m2) * 2 == m3 / 2);
    BOOST_CHECK(((m3 - m3) == Matrix<num, 2, 3>(0)));
    BOOST_CHECK((m2.Row(0) == v2));
    BOOST_CHECK((m2.T().Col(1) == v3));
    BOOST_CHECK(((m3 / 4 * v2.T()).T() * Vector<num, 2>(12, -20).T()) == 0);
    Matrix<num, 3, 3> m4{ {5, -4, 2}, {-1, 2, 3}, {-2, 1, 0}};
    BOOST_CHECK(m4.NormSquare() == 64);

    Vector<num, 3> v4{0, 1, 2};
    Matrix<num, 3, 3> m5{{0, 1, 2}, {3, 4, 5}, {6, 7, 8}};
    BOOST_CHECK((m5 * v4 == Vector<num, 3>{5, 14, 23}));
    BOOST_CHECK((v4 * m5 == Vector<num, 3>{15, 18, 21}));
}

test_suite * create_math_test_suite()
{
    test_suite * math_suite = BOOST_TEST_SUITE("s_math");
    //
    math_suite->add(BOOST_TEST_CASE_TEMPLATE(t_math_utility_t, t_math_num_types));
    math_suite->add(BOOST_TEST_CASE(&t_math_utility));
    math_suite->add(BOOST_TEST_CASE_TEMPLATE(t_math_linear_algebra_t, t_math_num_types));
    //
    return math_suite;
}
#endif//TEST_TESTMATH_HPP
