/**
 * @file Simulator.hpp
 * @author bwu
 * @brief static/transient circuit simulator
 * @version 0.1
 * @date 2023-07-16
 */

#include "MNA.hpp"

#ifdef EIGEN_LIBRARY_SUPPORT
#include <memory>
#include <set>
namespace generic::ckt {

using namespace Eigen;
class Circuit
{
public:
    Circuit(size_t nodes, std::set<size_t> sources, std::set<size_t> probs = {})
     : m_nodes(nodes)
     , m_probs(std::move(probs))
     , m_sources(std::move(sources))
    {
        m_G.reset(new MatrixXd); *m_G = MatrixXd::Zero(Nodes() + SNodes(), Nodes() + SNodes());
        m_C.reset(new MatrixXd); *m_C = MatrixXd::Zero(Nodes() + SNodes(), Nodes() + SNodes());
        m_B.reset(new MatrixXd); *m_B = MatrixXd::Zero(Nodes() + SNodes(), SNodes());
        m_L.reset(new MatrixXd); *m_L = MatrixXd::Zero(Nodes() + SNodes(), PNodes());
        
        size_t i = 0;
        for (auto n : m_sources) {
            B()(i + Nodes(), i) = -1;
            mna::StampI(G(), n, i + Nodes());  ++i;
        }
        i = 0;
        if (not Probs().empty()) { for (auto n : Probs()) L()(n, i++) = 1;}
        else { for (; i < Nodes(); ++i) L()(i, i) = 1; }
    }

    virtual ~Circuit() = default;

    void SetR(size_t i, size_t j, double r)
    {
        mna::Stamp(G(), i, j, 1 / r);
    }

    void SetC(size_t i, size_t j, double c)
    {
        mna::Stamp(C(), i, j, c);
    }

    void SetC(size_t i, double c)
    {
        mna::Stamp(C(), i, c);
    }
    
    MatrixXd & G() { return *m_G; } 
    MatrixXd & C() { return *m_C; } 
    MatrixXd & B() { return *m_B; } 
    MatrixXd & L() { return *m_L; } 

    const MatrixXd & G() const { return *m_G; } 
    const MatrixXd & C() const { return *m_C; } 
    const MatrixXd & B() const { return *m_B; } 
    const MatrixXd & L() const { return *m_L; } 

    size_t Nodes() const { return m_nodes; }
    size_t SNodes() const { return Sources().size(); }
    size_t PNodes() const { return Probs().empty() ? Nodes() : Probs().size(); }

    const std::set<size_t> & Probs() const { return m_probs; }
    const std::set<size_t> & Sources() const { return m_sources; }

private:
    size_t m_nodes;
    std::set<size_t> m_probs;
    std::set<size_t> m_sources;
    std::unique_ptr<MatrixXd> m_G{nullptr};
    std::unique_ptr<MatrixXd> m_C{nullptr};
    std::unique_ptr<MatrixXd> m_B{nullptr};
    std::unique_ptr<MatrixXd> m_L{nullptr};
};

using StateType = std::vector<double>;
template <typename VoltageFunc>
struct Simulator
{
    VectorXd u;
    MatrixXd coeff;
    MatrixXd input;
    MatrixXd Lred;
    const Circuit & ckt;
    VoltageFunc voltFun;
    bool verbose{false};
    Simulator(const Circuit & ckt, VoltageFunc && voltFun, bool verbose = false)
     : ckt(ckt), voltFun(std::move(voltFun)), verbose(verbose)
    {
        MatrixXd Gred, Cred, Bred;
        if (verbose) {
            std::cout << "G:\n" << ckt.G() << std::endl;
            std::cout << "C:\n" << ckt.C() << std::endl;
            std::cout << "B:\n" << ckt.B() << std::endl;
            std::cout << "L:\n" << ckt.L() << std::endl;
        }

        std::tie(Gred, Cred, Bred, Lred) = mna::RegularizeSuDynamic(ckt.G(), ckt.C(), ckt.B(), ckt.L());
        
        if (verbose) {
            std::cout << "Gred:\n" << Gred << std::endl;
            std::cout << "Cred:\n" << Cred << std::endl;
            std::cout << "Bred:\n" << Bred << std::endl;
            std::cout << "Lred:\n" << Lred << std::endl;
        }

        u.resize(ckt.Sources().size());
        coeff = Cred.ldlt().solve(-1.0 * Gred);
        input = Cred.ldlt().solve(Bred);
    }

    void operator() (const StateType x, StateType & dxdt, double t)
    {
        for (size_t i = 0; i < ckt.SNodes(); ++i) u(i) = voltFun(i, t);

        Map<const VectorXd> xvec(x.data(), x.size());
        Map<VectorXd> result(dxdt.data(), dxdt.size());
        result = coeff * xvec + input * u;  
    }

    size_t StateSize() const { return coeff.cols(); }

    std::vector<double> State2Output(const StateType & x)
    {
        std::vector<double> result(ckt.PNodes());
        Map<const VectorXd> xvec(x.data(), x.size());
        Map<VectorXd> ovec(result.data(), result.size());
        ovec = Lred.transpose() * xvec;
        return result;
    }
};

} // namespace generic::ckt

#endif // EIGEN_LIBRARY_SUPPORT