#include "generic/circuit/Simulator.hpp"
#include "generic/circuit/MOR.hpp"

using namespace generic::ckt;
int main()
{
    using float_t = double;
    float_t vdd = 0.88;
    float_t ts  = 200e-12;
    auto steps = 1000;
    auto vfun = [vdd, ts](size_t i, auto t) -> float_t
    {
        if (i == 0)
            return t > ts ? vdd : t * vdd / ts;
        else return t > ts ? 0 : vdd - t * vdd / ts;
    };

    // auto ckt = SparseCircuit<float_t>(11, {0, 6, 11}, {0, 6, 11});
    auto ckt = DenseCircuit<float_t>(11, {0, 6}, {0, 6, 11});
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

    const auto & m = ckt.Build();
    // const auto x = Prima(m, 5);

    // auto xt = x.transpose();
    // auto rC = xt * m.C * x;
    // auto rG = xt * m.G * x;
    // auto rB = xt * m.B;
    // auto rL = xt * m.L;
    
    // std::cout <<  "x:\n" <<  x << std::endl;
    // std::cout << "rC:\n" << rC << std::endl;
    // std::cout << "rG:\n" << rG << std::endl;
    // std::cout << "rB:\n" << rB << std::endl;
    // std::cout << "rL:\n" << rL << std::endl;

    auto im = Intermidiate<float_t>(m, false);
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
    std::cout << std::endl;
}