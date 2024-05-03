#include "generic/circuit/Simulator.hpp"
#include "generic/circuit/MOR.hpp"
using namespace generic;
using namespace generic::ckt;

template <typename float_t, template <typename> class circuit_t>
inline circuit_t<float_t> makeCircuit()
{
    auto ckt = circuit_t<float_t>(12, {0, 6}, {3, 9});
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

    return ckt;
}

int main()
{
    using float_t = double;
    auto ckt = makeCircuit<float_t, SparseCircuit>();    
    const auto & m = ckt.Build();
    auto rm = Reduce(m, 6); 

    std::cout << "origin: " << m << std::endl;
    std::cout << "reduce: " << rm.m << std::endl;
    auto momentsOrigin = mna::Moments(DenseMatrix<float_t>(m.C),
                                      DenseMatrix<float_t>(m.G),
                                      DenseMatrix<float_t>(m.B), 
                                      DenseMatrix<float_t>(m.L), 4);
    auto momentsReduce = Moments(rm, 4);
                                    
    for (size_t i = 0; i < 3; ++i) {
        std::cout << "origin moments " << i << ": " << std::endl;
        std::cout << momentsOrigin.at(i) << std::endl;
        std::cout << "reduce moments " << i << ": " << std::endl;
        std::cout << momentsReduce.at(i) << std::endl;
    }
    // auto [nm, p] = mna::RegularizeSu(m);
    // std::cout << "origin: " << m << std::endl;
    // std::cout << "regularize: " << nm << std::endl;
    // std::cout << "permutation: " << p.indices() << std::endl;

    // auto pi = RetrievePiModel(rm, 0, 0.01);
    // std::cout << "pi model: " << pi << std::endl;
}