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

#if EIGEN_LIBRARY_SUPPORT
using namespace generic::ckt;

void t_simulator()
{
    double vdd = 0.88;
    double ts  = 200e-12;
    auto steps = 1000;
    auto vfun = [vdd, ts](size_t, auto t) { return t > ts ? vdd : t * vdd / ts; };

    auto ckt = CircuitV(2, {0}, {1});
    ckt.SetR(0, 1, 1e2);
    ckt.SetC(1, 1e-12);
    ckt.SetC(2, 1e-12);

    Intermidiate im(ckt, false);
    StateType initState(im.StateSize(), 0);

    using namespace boost::numeric;
	using ErrorStepperType = odeint::runge_kutta_cash_karp54<StateType>;
	odeint::integrate_adaptive(
        odeint::make_controlled(1e-12, 1e-10,ErrorStepperType{}),
        Simulator(im, vfun), initState, 0.0, ts * 2, ts * 2 / steps);

    std::vector<double> outs;
    outs = im.State2Output(initState);
    BOOST_CHECK(outs.size() == 2);
    BOOST_CHECK_CLOSE(outs.back(), 0.82852, 1e-2);
}

#endif //EIGEN_LIBRARY_SUPPORT

test_suite * create_circuit_test_suite()
{
    test_suite * circuit_suite = BOOST_TEST_SUITE("s_circuit");
    //
#if EIGEN_LIBRARY_SUPPORT
    circuit_suite->add(BOOST_TEST_CASE(&t_simulator));
#endif //EIGEN_LIBRARY_SUPPORT
    //
    return circuit_suite;
}