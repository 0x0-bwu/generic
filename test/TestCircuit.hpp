/**
 * @file TestCircuit.hpp
 * @author bwu
 * @brief Unit test of classes and functions in namespace ckt
 * @version 0.1
 * @date 2023-07-23
 */
#pragma once
#include "TestCommon.hpp"
#include "generic/circuit/Simulator.hpp"
#include "generic/circuit/MOR.hpp"

using namespace boost::unit_test;
using namespace generic;
using namespace generic::ckt;
using t_ckt_num_types = boost::mpl::list<float, double>;

void t_dense_ckt_simulator()
{
    using float_t = double; 
    float_t vdd = 0.88;
    float_t ts  = 200e-12;
    auto steps = 1000;
    auto vfun = [vdd, ts](size_t, auto t) -> float_t { return t > ts ? vdd : t * vdd / ts; };

    auto ckt = DenseCircuit<float_t>(3, {0}, {1, 2});
    ckt.SetR(0, 1, 0.01);
    ckt.SetR(1, 2, 1e2);
    ckt.SetC(1, 1e-12);
    ckt.SetC(2, 1e-12);

    auto im = Intermidiate<float_t>(ckt.Build(), false);
    using namespace boost::numeric;
    using StateType = typename Intermidiate<float_t>::StateType;
	using ErrorStepperType = odeint::runge_kutta_cash_karp54<StateType>;
    
    StateType initState(im.StateSize(), 0);
	odeint::integrate_adaptive(
        odeint::make_controlled(float_t{1e-12}, float_t{1e-10}, ErrorStepperType{}),
        Simulator(im, vfun), initState, float_t{0}, float_t{ts * 2}, float_t{ts * 2 / steps});

    auto outs = im.State2Output(initState);
    BOOST_CHECK(outs.size() == 2);
    BOOST_CHECK_CLOSE(outs.front(), 0.8800, 1e-2);
    BOOST_CHECK_CLOSE(outs.back(), 0.82852, 1e-2);
}

void t_dense_ckt_cross_talk()
{
    using float_t = double;
    float_t vdd = 0.88;
    float_t ts  = 200e-12;
    auto steps = 1000;
    auto vfun = [vdd, ts](size_t i, auto t) -> float_t
    {
        if (i == 0) return t > ts ? vdd : t * vdd / ts;
        else return t < ts ? vdd : vdd - (t - ts) * vdd / ts;
    };

    auto ckt = DenseCircuit<float_t>(12, {0, 6}, {3, 9});
    ckt.SetR(0, 1, 0.01);
    ckt.SetR(6, 7, 0.01);

    ckt.SetC(1, 1e-12);
    ckt.SetC(7, 1e-12);
    ckt.SetR(1, 2, 10);
    ckt.SetR(7, 8, 10);

    ckt.SetC(2, 1e-12);
    ckt.SetC(8, 1e-12);
    ckt.SetR(2, 3, 10);
    ckt.SetR(8, 9, 10);

    ckt.SetC(3, 1e-12);
    ckt.SetC(9, 1e-12);
    ckt.SetR(3, 4, 10);
    ckt.SetR(9, 10, 10);
    ckt.SetC(3, 9, 1e-12);
    
    ckt.SetC(4, 1e-12);
    ckt.SetC(10, 1e-12);
    ckt.SetR(4, 5, 10);
    ckt.SetR(10, 11, 10);

    ckt.SetC(5, 1e-12);
    ckt.SetC(11, 1e-12);
    const auto & m = ckt.Build();
    auto im = Intermidiate<float_t>(m, false);
    using namespace boost::numeric;
    using StateType = typename Intermidiate<float_t>::StateType;
	using ErrorStepperType = odeint::runge_kutta_cash_karp54<StateType>;

    struct Observer
    {
        float_t & max3;
        float_t & max9;
        const Intermidiate<float_t> & im;
        Observer(const Intermidiate<float_t> & im, float_t & max3, float_t & max9) : max3(max3), max9(max9), im(im) {}

        void operator() (const StateType & x, float_t t) {
            auto out = im.State2Output(x);
            max3 = std::max(max3, out.front());
            max9 = std::max(max9, out.back());
        }
    };
    
    StateType out;
    float_t max3{0}, max9{0};
    StateType state(im.StateSize(), 0);
	odeint::integrate_adaptive(
        odeint::make_controlled(float_t{1e-12}, float_t{1e-10}, ErrorStepperType{}),
        Simulator(im, vfun), state, float_t{0}, float_t{ts * 2}, float_t{ts * 2 / steps}, Observer(im, max3, max9));
    
    BOOST_CHECK_CLOSE(max3, 0.784533, 1e-2);
    BOOST_CHECK_CLOSE(max9, 0.842692, 1e-2);
}

void t_mor_ckt_cross_talk()
{
    using float_t = double;
    float_t vdd = 0.88;
    float_t ts  = 200e-12;
    auto steps = 1000;
    auto vfun = [vdd, ts](size_t i, auto t) -> float_t
    {
        if (i == 0) return t > ts ? vdd : t * vdd / ts;
        else return t < ts ? vdd : vdd - (t - ts) * vdd / ts;
    };
    auto ckt = SparseCircuit<float_t>(12, {0, 6}, {3, 9});
    ckt.SetR(0, 1, 0.01);
    ckt.SetR(6, 7, 0.01);

    ckt.SetC(1, 1e-12);
    ckt.SetC(7, 1e-12);
    ckt.SetR(1, 2, 10);
    ckt.SetR(7, 8, 10);

    ckt.SetC(2, 1e-12);
    ckt.SetC(8, 1e-12);
    ckt.SetR(2, 3, 10);
    ckt.SetR(8, 9, 10);

    ckt.SetC(3, 1e-12);
    ckt.SetC(9, 1e-12);
    ckt.SetR(3, 4, 10);
    ckt.SetR(9, 10, 10);
    ckt.SetC(3, 9, 1e-12);
    
    ckt.SetC(4, 1e-12);
    ckt.SetC(10, 1e-12);
    ckt.SetR(4, 5, 10);
    ckt.SetR(10, 11, 10);

    ckt.SetC(5, 1e-12);
    ckt.SetC(11, 1e-12);   
    const auto & m = ckt.Build();
    auto rm = Reduce(m, 5);    
    auto im = Intermidiate<float_t>(rm.m, false);
    using namespace boost::numeric;
    using StateType = typename Intermidiate<float_t>::StateType;
	using ErrorStepperType = odeint::runge_kutta_cash_karp54<StateType>;

    struct Observer
    {
        float_t & max3;
        float_t & max9;
        const Intermidiate<float_t> & im;
        Observer(const Intermidiate<float_t> & im, float_t & max3, float_t & max9) : max3(max3), max9(max9), im(im) {}

        void operator() (const StateType & x, float_t t) {
            auto out = im.State2Output(x);
            max3 = std::max(max3, out.front());
            max9 = std::max(max9, out.back());
        }
    };
    
    StateType out;
    float_t max3{0}, max9{0};
    StateType state(im.StateSize(), 0);
	odeint::integrate_adaptive(
        odeint::make_controlled(float_t{1e-12}, float_t{1e-10}, ErrorStepperType{}),
        Simulator(im, vfun), state, float_t{0}, float_t{ts * 2}, float_t{ts * 2 / steps}, Observer(im, max3, max9));
    
    BOOST_CHECK_CLOSE(max3, 0.784533, 2.0);
    BOOST_CHECK_CLOSE(max9, 0.842692, 2.5);   
}

BOOST_TEST_CASE_TEMPLATE_FUNCTION(t_mor_retrieve_pi_model, float_t)
{
    auto ckt = SparseCircuit<float_t, false>(2, {0}, {0, 1}, 0.984375);
    ckt.SetR(0, 1, 364.250);
    ckt.SetC(0, 0.00345);
    ckt.SetC(1, 0.00345 + 0.000433);
    const auto & m = ckt.Build();
    auto rm = Reduce(m, 2); 
    size_t mm = 4;
    auto momentsOrigin = mna::Moments(DenseMatrix<float_t>(m.C),
                                      DenseMatrix<float_t>(m.G),
                                      DenseMatrix<float_t>(m.B), 
                                      DenseMatrix<float_t>(m.L), mm);
    auto momentsReduce = Moments(rm, mm);
                                    
    for (size_t i = 0; i < mm; ++i) {
        BOOST_CHECK_CLOSE(momentsOrigin.at(i)(0, 0), momentsReduce.at(i)(0, 0), 1e-2);
        BOOST_CHECK_CLOSE(momentsOrigin.at(i)(1, 0), momentsReduce.at(i)(1, 0), 1e-2);
    }
    auto pi = RetrievePiModel(rm.m, 0, ckt.rd);
    BOOST_CHECK_CLOSE(pi.c1, 0.00345 , 1e-2);
    BOOST_CHECK_CLOSE(pi.c2, 0.003883, 1e-2);
    BOOST_CHECK_CLOSE(pi.r , 364.25  , 1e-2);
}

test_suite * create_circuit_test_suite()
{
    test_suite * circuit_suite = BOOST_TEST_SUITE("s_circuit");
    //
    circuit_suite->add(BOOST_TEST_CASE(&t_dense_ckt_simulator));
    circuit_suite->add(BOOST_TEST_CASE(&t_dense_ckt_cross_talk));
    circuit_suite->add(BOOST_TEST_CASE(&t_mor_ckt_cross_talk));
    circuit_suite->add(BOOST_TEST_CASE_TEMPLATE(t_mor_retrieve_pi_model, t_ckt_num_types));
    //
    return circuit_suite;
}