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

template <typename State2OutFunc>
struct DelayCalculator
{
    double & delay;
    bool recorded{false};
    State2OutFunc state2out;
    DelayCalculator(State2OutFunc && state2out, double & delay)
     : state2out(std::move(state2out)), delay(delay) {}

    void operator() (const StateType & state, double t)
    {
        if (not recorded) {
            auto out = state2out(state);
            if (out.back() > 0.44) {
                delay = t;
                recorded = true;
            }
        }
    }
};

void t_simulator()
{
    double vdd = 0.88;
    double ts  = 200e-12;
    auto steps = 1000;
    auto vfun = [vdd, ts](size_t, auto t) { return t > ts ? vdd : t * vdd / ts; };

    auto ckt = Circuit(3, {0});
    ckt.SetR(0, 1, 1);
    ckt.SetR(1, 2, 1e4);
    ckt.SetC(1, 1e-14);
    ckt.SetC(2, 1e-14);

    auto simulator = Simulator(ckt, std::move(vfun), false);
    StateType initState(simulator.StateSize(), 0);

    double delay{0};
    auto func = [&](const StateType & x){ return simulator.State2Output(x); };
    auto calc = DelayCalculator(std::move(func), delay);

    using namespace boost::numeric;
	using ErrorStepperType = odeint::runge_kutta_cash_karp54<StateType>;
	odeint::integrate_adaptive(
        odeint::make_controlled(1e-12, 1e-10,ErrorStepperType{}),
        simulator, initState, 0.0, ts * 2, ts * 2 / steps, calc);

    BOOST_CHECK_CLOSE(delay, 1.84194e-10, 1e-3);
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