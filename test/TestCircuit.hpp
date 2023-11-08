/**
 * @file TestCircuit.hpp
 * @author bwu
 * @brief Unit test of classes and functions in namespace ckt
 * @version 0.1
 * @date 2023-07-23
 */
#pragma once
#include "generic/test/TestCommon.hpp"
#include "generic/circuit/Simulator.hpp"

using namespace boost::unit_test;
using namespace generic;
using namespace generic::ckt;
using t_ckt_num_types = boost::mpl::list<double>;

BOOST_TEST_CASE_TEMPLATE_FUNCTION(t_simulator_t, float_t)
{
    float_t vdd = 0.88;
    float_t ts  = 200e-12;
    auto steps = 1000;
    auto vfun = [vdd, ts](size_t, auto t) -> float_t { return t > ts ? vdd : t * vdd / ts; };

    auto ckt = DenseCircuit<float_t>(3, {0}, {1, 2});
    ckt.SetR(0, 1, 0.01);
    ckt.SetR(1, 2, 1e2);
    ckt.SetC(1, 1e-12);
    ckt.SetC(2, 1e-12);

    auto im = Intermidiate(ckt, true);
    using namespace boost::numeric;
    using StateType = typename Intermidiate<float_t>::StateType;
	using ErrorStepperType = odeint::runge_kutta_cash_karp54<StateType>;
    
    StateType initState(im.StateSize(), 0);
	odeint::integrate_adaptive(
        odeint::make_controlled(float_t{1e-12}, float_t{1e-10}, ErrorStepperType{}),
        Simulator(im, vfun), initState, float_t{0}, float_t{ts * 2}, float_t{ts * 2 / steps});

    auto outs = im.State2Output(initState);
    for (auto out : outs)
        std::cout << out << ',';
    BOOST_CHECK(outs.size() == 2);
    BOOST_CHECK_CLOSE(outs.front(), 0.8800, 1e-2);
    BOOST_CHECK_CLOSE(outs.back(), 0.82852, 1e-2);
}

test_suite * create_circuit_test_suite()
{
    test_suite * circuit_suite = BOOST_TEST_SUITE("s_circuit");
    //
    circuit_suite->add(BOOST_TEST_CASE_TEMPLATE(t_simulator_t, t_ckt_num_types));
    //
    return circuit_suite;
}