#include "generic/circuit/Simulator.hpp"
#include "generic/circuit/MOR.hpp"
#include "generic/tools/Format.hpp"
#include <boost/math/interpolators/pchip.hpp>
#include <iostream>

using namespace generic;
using namespace generic::ckt;

void test()
{
    using float_t = double;
    constexpr double minR = 0.01; 
    auto ckt = DenseCircuit<float_t>(11, {10}, {1, 2});
    ckt.SetR(0, 9, minR);
    ckt.SetR(9, 8, 0.001);
    ckt.SetR(9, 7, 108.647);
    ckt.SetR(7, 6, 74.6107);
    ckt.SetR(7, 4, 93.7641);
    ckt.SetR(6, 5, 0.001);
    ckt.SetR(6, 2, minR);
    ckt.SetR(4, 3, 0.001);
    ckt.SetR(4, 1, minR);
    ckt.SetR(10, 0, minR);
    
    // ckt.SetC(0, minC);
    ckt.SetC(1, 0.001912);
    ckt.SetC(2, 0.003199);


    auto im = Intermidiate<float_t>(ckt.Build(), true);
    using namespace boost::numeric;
    using StateType = typename Intermidiate<float_t>::StateType;
	using ErrorStepperType = odeint::runge_kutta_cash_karp54<StateType>;
    StateType initState(im.StateSize(), 0);

    std::vector<float_t> ts{0.000000,  3.235554,  5.883425,  8.024378,  9.789495, 11.402321, 12.936249, 14.435669, 16.252773, 18.872063};
    std::vector<float_t> vs{0,         0.1,       0.2,       0.3,       0.4,      0.5,       0.6,       0.7,       0.8,       0.9,     };
    auto interp = boost::math::interpolators::pchip<std::vector<float_t>>(std::vector<float_t>(ts), std::vector<float_t>(vs));
    for (size_t i = 1; i < 10; ++i) {
        float_t t0 = ts[i - 1];
        float_t t1 = ts[i];
        float_t dt = t1 - t0;
        auto vfun = [&](size_t, float_t t)-> double { 
            return interp(t);
        };
        odeint::integrate_adaptive(
            odeint::make_controlled(float_t{1e-6}, float_t{1e-4}, ErrorStepperType{}),
            Simulator(im, std::move(vfun)), initState, t0, t1, dt / 10);
        auto outs = im.State2Output(initState);
        std::cout << "t1: " << t1 << " out: " << fmt::Fmt2Str(outs, ",") << std::endl;
    }
}

int main()
{
    test();
    return 0;
}