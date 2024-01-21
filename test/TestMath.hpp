/**
 * @file TestMath.hpp
 * @author bwu
 * @brief Unit test of classes and functions in namespace math
 * @version 0.1
 * @date 2022-02-22
 */
#pragma once
#include "generic/test/TestCommon.hpp"
#include "generic/tools/FileSystem.hpp"
#include "generic/math/MathUtility.hpp"
#include "generic/math/LinearAlgebra.hpp"
#include "generic/math/PolynomialFit.hpp"
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

void t_math_io()
{
    using namespace la;
    using namespace io;
#if BOOST_GIL_IO_PNG_SUPPORT
    std::string outDir = fs::DirName(__FILE__).string() + "/data/out";
    DenseMatrix<float> dense(5, 5);
    for (Eigen::Index i = 0; i < dense.rows(); ++i) {
        dense(i, i) = i;
    }
    dense(0, 4) = 4;
    BOOST_CHECK(PatternView(dense, outDir + "/dense.png"));

    Triplets<double> t;
    SparseMatrix<double> sparse(1e6, 1e6);
    for (Eigen::Index i = 0; i < 1e3; ++i) {
        auto row = Random<int>(0, sparse.rows() - 1);
        auto col = Random<int>(row, sparse.rows() - 1);
        t.emplace_back(row, col, Random<double>(0, 1));
    }
    sparse.setFromTriplets(t.begin(), t.end());
    BOOST_CHECK(PatternView(sparse, outDir + "/sparse.png"));
#endif
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
    BOOST_CHECK(Within<OpenInterval>(1e-4, 1.0e-4, 1.1e-4) == false);
    BOOST_CHECK(Within<LeftClosedRightOpen>(1e-4, 1.0e-4, 1.1e-4) == true);

    //Equation Solve
    double rc{10}, tr{20};
    auto func  = [&] (auto x) { return x - rc + rc * std::exp(-(0.5 * tr + x) / rc); };
    auto dfunc = [&] (auto x) { return 1 - std::exp(-(0.5 * tr + x) / rc); };
    BOOST_CHECK_CLOSE(Bisection<double>(func, 0, 10, 1e-10), 8.41406, 1e-1);
    BOOST_CHECK_CLOSE(NewtonRaphson<double>(func, dfunc, 10, 1e-10), 8.41406, 1e-1);
}

void t_math_polynomial_fit()
{
    // y = 2 - x;
    auto coeffs = PolyFit(std::vector<double>{0, 1, 2}, std::vector<double>{2, 1, 0}, 1);
    BOOST_CHECK(coeffs[0] ==  2);
    BOOST_CHECK(coeffs[1] == -1);

    auto yValues = Polyval(coeffs, {3, 4});
    BOOST_CHECK(yValues[0] == -1);
    BOOST_CHECK(yValues[1] == -2);
}

BOOST_TEST_CASE_TEMPLATE_FUNCTION(t_math_linear_algebra_t, math_num_types)
{
    using namespace la;
    using num = math_num_types;
}

test_suite * create_math_test_suite()
{
    test_suite * math_suite = BOOST_TEST_SUITE("s_math");
    //
    math_suite->add(BOOST_TEST_CASE_TEMPLATE(t_math_utility_t, t_math_num_types));
    math_suite->add(BOOST_TEST_CASE(&t_math_io));
    math_suite->add(BOOST_TEST_CASE(&t_math_utility));
    math_suite->add(BOOST_TEST_CASE(&t_math_polynomial_fit));
    math_suite->add(BOOST_TEST_CASE_TEMPLATE(t_math_linear_algebra_t, t_math_num_types));
    //
    return math_suite;
}