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
#include "generic/math/la/Common.hpp"
#include "generic/math/Interpolation.hpp"
#include "generic/math/PolynomialFit.hpp"
#include "generic/math/LookupTable.hpp"
#include "generic/math/FastMath.hpp"
#include "generic/math/MathIO.hpp"
#include "generic/math/Filter.hpp"

#include <boost/math/interpolators/pchip.hpp>

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
#ifdef GENERIC_BOOST_GIL_IO_PNG_SUPPORT
    DenseMatrix<float> dense(5, 5);
    for (Eigen::Index i = 0; i < dense.rows(); ++i) {
        dense(i, i) = i;
    }
    dense(0, 4) = 4;
    BOOST_CHECK(PatternView(dense, GetTestOutDataPath() + "/dense.png"));

    Triplets<double> t;
    SparseMatrix<double> sparse(1e6, 1e6);
    for (Eigen::Index i = 0; i < 1e3; ++i) {
        auto row = Random<int>(0, sparse.rows() - 1);
        auto col = Random<int>(row, sparse.rows() - 1);
        t.emplace_back(row, col, Random<double>(0, 1));
    }
    sparse.setFromTriplets(t.begin(), t.end());
    BOOST_CHECK(PatternView(sparse, GetTestOutDataPath() + "/sparse.png"));
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
    
    //toBits
    {
        BOOST_CHECK(toBits(0.0003153f).to_string() == "00111001101001010100111011011010");
    }

    //Equation Solve
    double rc{10}, tr{20};
    auto func  = [&] (auto x) { return x - rc + rc * std::exp(-(0.5 * tr + x) / rc); };
    auto dfunc = [&] (auto x) { return 1 - std::exp(-(0.5 * tr + x) / rc); };
    BOOST_CHECK_CLOSE(Bisection<double>(func, 0, 10, 1e-10), 8.41406, 1e-1);
    BOOST_CHECK_CLOSE(NewtonRaphson<double>(func, dfunc, 10, 1e-10), 8.41406, 1e-1);
    
    //Mean and Variance
    {
        std::list<double> input{83.1192,85.2815,83.3486,84.1985,85.8013,83.8569,85.4694,86.4357,82.5886,86.266,87.5939,83.908};
        auto [mean, variance] = MeanAndVariance(input.begin(), input.end());
        BOOST_CHECK_CLOSE(mean, 84.8230, 0.1);
        BOOST_CHECK_CLOSE(variance, 2.37963, 0.1);
    }
    
    //Random Series
    {
        auto series = generic::math::RandomSeries<double>(1e5, 30, 3);
        auto result = generic::math::MeanAndVariance(series.begin(), series.end());
        BOOST_CHECK_CLOSE(result.first, 30, 1);
        BOOST_CHECK_CLOSE(result.second, 9, 1);
    }

    //Next/Prev Pwo Tow
    {
        BOOST_CHECK(2 == PrevPowTwo(math::pi));
        BOOST_CHECK(4 == NextPowTwo(math::pi));
    }
}

void t_math_polynomial_fit()
{
    double t = 1e-6;
    // y = 2 - x;
    auto coeffs = PolyFit(std::vector<double>{0, 1, 2}, std::vector<double>{2, 1, 0}, 1);
    BOOST_CHECK_CLOSE(coeffs[0],  2, t);
    BOOST_CHECK_CLOSE(coeffs[1], -1, t);

    auto yValues = PolyVal(coeffs, std::vector<double>{3, 4});
    BOOST_CHECK_CLOSE(yValues[0], -1, t);
    BOOST_CHECK_CLOSE(yValues[1], -2, t);
}

void t_math_filter()
{
    // simple moving average
    {
        auto before = std::vector<double>{0.004713,  0.005984,  0.007222,  0.008698,  0.010947,  0.013219,  0.014986,  0.016344,  0.017281,  0.017872,  0.018200,  0.018428,  0.018614,  0.018798,  0.018978,  0.019107,  0.019163,  0.019146,  0.019153,  0.019239,  0.019361,  0.019523,  0.019816,  0.020364,  0.020660,  0.021032,  0.021221,  0.023790,  0.025440,  0.026777,  0.027209,  0.026668,  0.026257,  0.026438,  0.027470,  0.027989,  0.027361,  0.028195,  0.028478,  0.029014,  0.029014,  0.029014,  0.029014,  0.029014};
        auto after = SimpleMovingAverage(before.data(), before.size(), 2);
        BOOST_CHECK_CLOSE(after[ 0], 0.0054690, 1e-3);
        BOOST_CHECK_CLOSE(after[43], 0.0289874, 1e-3);
        BOOST_CHECK_CLOSE(after[10], 0.0178824, 1e-3);
    }
}

void t_math_fast_math()
{
    float maxErr{0};
    constexpr size_t total = 1e5;
    for (size_t i = 0; i < total; ++i) {
        float num = Random<float>(-pi, pi);
        maxErr = std::max(maxErr, std::abs(std::sin(num) - FastSin(num)));
    }
    BOOST_CHECK(maxErr < 4e-5);

    maxErr = 0;
    for (size_t i = 0; i < total; ++i) {
        float num = Random<float>(-pi, pi);
        maxErr = std::max(maxErr, std::abs(std::cos(num) - FastCos(num)));
    }
    BOOST_CHECK(maxErr < 4e-5);

    maxErr = 0;
    for (size_t i = 0; i < total; ++i) {
        float num = Random<float>(0 + std::numeric_limits<float>::epsilon(), std::numeric_limits<float>::max());
        maxErr = std::max(maxErr, std::abs(std::log2(num) - FastLog2(num)));
    }
    BOOST_CHECK(maxErr < 2e-4);

    maxErr = 0;
    for (size_t i = 0; i < total; ++i) {
        float num = Random<float>(0 + std::numeric_limits<float>::epsilon(), std::numeric_limits<float>::max());
        maxErr = std::max(maxErr, std::abs(std::log2(num) - FasterLog2(num)));
    }
    BOOST_CHECK(maxErr < 6e-2);

    maxErr = 0;
    for (size_t i = 0; i < total; ++i) {
        float num = Random<float>(0 + std::numeric_limits<float>::epsilon(), std::numeric_limits<float>::max());
        maxErr = std::max(maxErr, std::abs(std::log(num) - FastLog(num)));
    }
    BOOST_CHECK(maxErr < 2e-4);

    maxErr = 0;
    for (size_t i = 0; i < total; ++i) {
        float num = Random<float>(0 + std::numeric_limits<float>::epsilon(), std::numeric_limits<float>::max());
        maxErr = std::max(maxErr, std::abs(std::log(num) - FasterLog(num)));
    }
    BOOST_CHECK(maxErr < 4e-2);

    maxErr = 0;
    for (size_t i = 0; i < total; ++i) {
        float num = Random<float>(-126, 0);
        maxErr = std::max(maxErr, std::abs(std::pow(2.0f, num) - FastPow2(num)));
    }
    BOOST_CHECK(maxErr < 4.5e-5);

    maxErr = 0;
    for (size_t i = 0; i < total; ++i) {
        float num = Random<float>(-126, 0);
        maxErr = std::max(maxErr, std::abs(std::pow(2.0f, num) - FasterPow2(num)));
    }
    BOOST_CHECK(maxErr < 3e-2);

    maxErr = 0;
    for (size_t i = 0; i < total; ++i) {
        float num = Random<float>(-87, 0);
        maxErr = std::max(maxErr, std::abs(std::exp(num) - FastExp(num)));
    }
    BOOST_CHECK(maxErr < 4.5e-5);

    maxErr = 0;
    for (size_t i = 0; i < total; ++i) {
        float num = Random<float>(-87, 0);
        maxErr = std::max(maxErr, std::abs(std::exp(num) - FasterExp(num)));
    }
    BOOST_CHECK(maxErr < 3e-2);

    // FastPchip correctness & batch tests
    {
        std::vector<float> xs{0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f};
        std::vector<float> ys(xs.size());
        for (size_t i = 0; i < xs.size(); ++i)
            ys[i] = std::sin(2.0f * 3.1415926f * xs[i]) + 1.0f;

        auto xs2 = xs;
        auto ys2 = ys;
        FastPchip<float> fastp(std::move(xs), std::move(ys));
        boost::math::interpolators::pchip<std::vector<float>> boostp(std::move(xs2), std::move(ys2));

        // random-query comparison
        float maxErr = 0.0f;
        const size_t Q = 2000;
        for (size_t k = 0; k < Q; ++k) {
            float q = Random<float>(0.1f, 0.9f);
            float a = fastp.Evaluate(q);
            float b = boostp(q);
            maxErr = std::max(maxErr, std::abs(a - b));
        }
        BOOST_CHECK_LT(maxErr, 1e-4f);

        // monotonic batch evaluation consistency
        std::vector<float> qxs{0.15f, 0.25f, 0.35f, 0.45f, 0.55f, 0.65f, 0.75f, 0.85f};
        std::vector<float> out;
        fastp.EvaluateBatchMonotonic(qxs, out);
        BOOST_CHECK_EQUAL(out.size(), qxs.size());
        for (size_t i = 0; i < qxs.size(); ++i) {
            float b = boostp(qxs[i]);
            BOOST_CHECK_SMALL(std::abs(out[i] - b), 1e-4f);
        }
    }
}

void t_math_interpolation()
{   
    // cubic
    {
        std::vector<double> x{0, 1, 3, 5, 9}, y, yy;
        std::transform(x.begin(), x.end(), std::back_inserter(y), [](auto x){ return  x * x * x + 2 * x * x + 3 * x + 4; });
        using Interp = math::Interpolation<std::vector<double>>;
        auto interp = Interp(x, y, Interp::Method::CUBIC, true);
        std::transform(x.begin(), x.end(), std::back_inserter(yy), [&](auto v){ return interp(v); });
        BOOST_TEST(y == yy, boost::test_tools::per_element());
    }
}

void t_math_lookup_table()
{
    using namespace math;
    // 1D
    {
        using Lut = LookupTable<int, 1>;
        Lut::Values values{2, 4};
        Lut::Indices indices{Lut::Values{1, 3}};
        Lut lut(indices, values);
        for (size_t i = 0; i < indices.size(); ++i) {
            BOOST_CHECK(lut[0][i] == indices[0][i]);
            BOOST_CHECK(lut(i) == values[i]);
        }
    }   
    // 2D
    {  
        using Lut = LookupTable<int, 2>;
        Lut::Values values{2, 4, 6, 6, 12, 18};
        Lut::Indices indices{Lut::Values{1, 3}, Lut::Values{2, 4, 6}};
        Lut lut(indices, values);
        for (size_t i = 0; i < indices.size(); ++i) {
            for (size_t j = 0; j < indices.at(i).size(); ++j) {
                BOOST_CHECK(lut[i][j] == indices[i][j]);
            }
        }

        size_t count = 0;
        for (size_t i = 0; i < indices.at(0).size(); ++i) {
            for (size_t j = 0; j < indices.at(1).size(); ++j) {
                BOOST_CHECK(lut(i, j) == values[count]);
                ++count;
            }
        }
    }
    //3D
    {
        using Lut = LookupTable<int, 3>;
        Lut::Values values{6, 12, 18, 24, 12, 24, 36, 48, 18, 36, 54, 72, 18, 36, 54, 72, 36, 72, 108, 144, 54, 108, 162, 216};
        Lut::Indices indices{Lut::Values{1, 3}, Lut::Values{2, 4, 6}, Lut::Values{3, 6, 9, 12}};
        Lut lut(indices, values);
        for (size_t i = 0; i < indices.size(); ++i) {
            for (size_t j = 0; j < indices.at(i).size(); ++j) {
                BOOST_CHECK(lut[i][j] == indices[i][j]);
            }
        }
        size_t count = 0;
        for (size_t i = 0; i < indices.at(0).size(); ++i) {
            for (size_t j = 0; j < indices.at(1).size(); ++j) {
                for (size_t k = 0; k < indices.at(2).size(); ++k) {
                    BOOST_CHECK(lut(i, j, k) == values[count]);
                    ++count;
                }
            }
        }
    }
}

BOOST_TEST_CASE_TEMPLATE_FUNCTION(t_math_linear_algebra_t, math_num_types)
{
    // using namespace la;
    // using num = math_num_types;
}

test_suite * create_math_test_suite()
{
    test_suite * math_suite = BOOST_TEST_SUITE("s_math");
    //
    math_suite->add(BOOST_TEST_CASE_TEMPLATE(t_math_utility_t, t_math_num_types));
    math_suite->add(BOOST_TEST_CASE(&t_math_io));
    math_suite->add(BOOST_TEST_CASE(&t_math_filter));
    math_suite->add(BOOST_TEST_CASE(&t_math_utility));
    math_suite->add(BOOST_TEST_CASE(&t_math_fast_math));
    math_suite->add(BOOST_TEST_CASE(&t_math_lookup_table));
    math_suite->add(BOOST_TEST_CASE(&t_math_interpolation));
    math_suite->add(BOOST_TEST_CASE(&t_math_polynomial_fit));
    math_suite->add(BOOST_TEST_CASE_TEMPLATE(t_math_linear_algebra_t, t_math_num_types));
    //
    return math_suite;
}