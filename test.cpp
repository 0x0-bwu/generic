
#include <boost/numeric/odeint.hpp>
#include <fstream>

#include "circuit/Simulator.hpp"
using namespace generic::ckt;
struct PushBackStateAndTime
{
    std::vector<double> & times;
    std::vector<StateType> & states;

    PushBackStateAndTime(std::vector<double> & times, std::vector<StateType> & states)
    : times(times), states(states) {}

    void operator() (const StateType & x , double t)
    {
        states.push_back(x);
        times.push_back(t);
    }
};

void test1()
{
    double vdd = 0.88;
    double ts = 200e-12;
    double dt = ts / 1000;
    auto vfun = [vdd, ts](size_t, auto t) { return t > ts ? vdd : t * vdd / ts; };

    auto ckt = Circuit(3, {0}, {1, 2});
    ckt.SetR(0, 1, 1);
    ckt.SetR(1, 2, 1e4);
    ckt.SetC(1, 1e-14);
    ckt.SetC(2, 1e-14);

    auto simulator = Simulator(ckt, std::move(vfun));
    StateType x(simulator.StateSize(), 0);

	std::vector<double> times;
    std::vector<StateType> stateHistory;
    
    using namespace boost::numeric;
	using ErrorStepperType = odeint::runge_kutta_cash_karp54<StateType>;
	odeint::integrate_adaptive(
        odeint::make_controlled(1e-10, 1e-6,ErrorStepperType{}),
        simulator, x, 0.0, ts * 5, dt, PushBackStateAndTime(times, stateHistory));

	std::ofstream out("./out1.txt");
	for (size_t i = 0; i < times.size(); ++i) {
		auto output = simulator.State2Output(stateHistory.at(i));
		out << times.at(i) << ' ' << output[0] << ' ' << output[1] << '\n';
	}
	out.close();
}

void test2()
{
    double vdd = 0.88;
    double ts = 200e-12;
    double dt = ts / 1000;
    auto vFun = [vdd, ts](size_t, auto t) { return t > ts ? vdd : t * vdd / ts; };

    auto ckt = Circuit(6, {0}, {});
    ckt.SetR(0, 5, 55.935);
    ckt.SetR(5, 2, 4.22393);
    ckt.SetR(5, 4, 55.935);
    ckt.SetR(4, 3, 107.263);
    ckt.SetR(3, 1, 55.935);

    ckt.SetC(0, 1e-6);
    ckt.SetC(5, 1.596e-5);
    ckt.SetC(2, 2.02e-5);
    ckt.SetC(4, 0.000369);
    ckt.SetC(3, 0.000369);
    ckt.SetC(1, 1e-6);

    auto simulator = Simulator(ckt, std::move(vFun));
    StateType x(simulator.StateSize(), 0);

	std::vector<double> times;
    std::vector<StateType> stateHistory;
    
    using namespace boost::numeric;
	using ErrorStepperType = odeint::runge_kutta_cash_karp54<StateType>;
	odeint::integrate_adaptive(
        odeint::make_controlled(1e-10, 1e-6,ErrorStepperType{}),
        simulator, x, 0.0, ts * 5, dt, PushBackStateAndTime(times, stateHistory));

	std::ofstream out("./out2.txt");
	for (size_t i = 0; i < times.size(); ++i) {
		auto output = simulator.State2Output(stateHistory.at(i));
		out << times.at(i) << ' ' << output[0] << ' ' << output[1] << '\n';
	}
	out.close();
}

int main()
{
    test1();
    test2();
}