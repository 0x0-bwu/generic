#include "generic/circuit/Simulator.hpp"
#include "generic/circuit/MOR.hpp"
#include "generic/math/Numbers.hpp"
#include <fstream>
using namespace generic;
using namespace generic::ckt;
int main()
{
    using float_t = double;
    float_t vdd = 88;
    float_t ts  = 200e-12;
    auto steps = 1000;
    auto vfun = [vdd, ts](size_t i, auto t) -> float_t
    {
        // if (i == 0) return t > ts ? vdd : t * vdd / ts;
        // else return t < ts ? vdd : vdd - (t - ts) * vdd / ts;
        if (i == 0) return 44 * std::sin(2 * math::pi / ts * t) + 44;
        else return 44 * std::cos(2 * math::pi / ts * t) + 44;
    };
    auto ckt = DenseCircuit<float_t>(12, {0, 6}, {1, 3, 7, 9});
    // auto ckt = SparseCircuit<float_t>(12, {0, 6}, {1, 3, 7, 9});
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
    // std::cout << "G:\n" << m.G << std::endl;
    // std::cout << "C:\n" << m.C << std::endl;
    // std::cout << "B:\n" << m.B << std::endl;
    // std::cout << "L:\n" << m.L << std::endl;
    // auto rm = Reduce(m, 6);    
    // std::cout <<  "x:\n" << rm.x   << std::endl;
    // std::cout << "rC:\n" << rm.m.C << std::endl;
    // std::cout << "rG:\n" << rm.m.G << std::endl;
    // std::cout << "rB:\n" << rm.m.B << std::endl;
    // std::cout << "rL:\n" << rm.m.L << std::endl;

    auto im = Intermidiate<float_t>(m, true);
    using namespace boost::numeric;
    using StateType = typename Intermidiate<float_t>::StateType;
	using ErrorStepperType = odeint::runge_kutta_cash_karp54<StateType>;

    struct Observer
    {
        StateType & out;
        std::ostream & os;
        const Intermidiate<float_t> & im;
        Observer(const Intermidiate<float_t> & im, StateType & out, std::ostream & os) : out(out), os(os), im(im) {}

        void operator() (const StateType & x, float_t t) {
            os << t;
            im.State2Output(x, out);
            for (auto & o : out)
                os << char(32) << o;
            os << std::endl;
        }

    };
    
    StateType out;
    StateType state(im.StateSize(), 0);
    std::ofstream os("./out.txt");
	odeint::integrate_adaptive(
        odeint::make_controlled(float_t{1e-12}, float_t{1e-10}, ErrorStepperType{}),
        Simulator(im, vfun), state, float_t{0}, float_t{ts * 100}, ts / 10, Observer(im, out, os));

    auto outs = im.State2Output(state);
    for (auto out : outs)
        std::cout << out << ',';
    std::cout << std::endl;
}