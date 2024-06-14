#include "generic/circuit/Simulator.hpp"
#include "generic/circuit/MOR.hpp"
using namespace generic;
using namespace generic::ckt;

template <class circuit_t>
inline circuit_t makeCircuit(double rd)
{
    auto ckt = circuit_t(2, {0}, {0, 1}, rd);
    ckt.SetR(0, 1, 364.250);

    ckt.SetC(0, 0.00345);
    ckt.SetC(1, 0.00345 + 0.000433);
    return ckt;
}

int main()
{
    using float_t = double;
    auto ckt = makeCircuit<SparseCircuit<float_t, false>>(0.984375);    
    const auto & m = ckt.Build();
    std::cout << "mna: " << m << std::endl;
    auto rm = Reduce(m, 2); 

    std::cout << "origin: " << m << std::endl;
    std::cout << "reduce: " << rm.m << std::endl;
    size_t mm = 4;
    auto momentsOrigin = mna::Moments(DenseMatrix<float_t>(m.C),
                                      DenseMatrix<float_t>(m.G),
                                      DenseMatrix<float_t>(m.B), 
                                      DenseMatrix<float_t>(m.L), mm);
    auto momentsReduce = Moments(rm, mm);
                                    
    for (size_t i = 0; i < mm; ++i) {
        std::cout << "origin moments " << i << ": " << std::endl;
        std::cout << momentsOrigin.at(i) << std::endl;
        std::cout << "reduce moments " << i << ": " << std::endl;
        std::cout << momentsReduce.at(i) << std::endl;
    }

    // auto dckt = makeCircuit<DenseCircuit<float_t, false>>(0.984375);    
    // auto dm = dckt.Build();

    auto pi = RetrievePiModel(rm.m, 0, ckt.rd);
    std::cout << "pi model: " << pi << std::endl;
}