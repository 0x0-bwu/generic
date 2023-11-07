/**
 * @file Simulator.hpp
 * @author bwu
 * @brief static/transient circuit simulator
 * @version 0.1
 * @date 2023-07-16
 */
#pragma once

#include "MNA.hpp"

#if EIGEN_LIBRARY_SUPPORT 
#include <boost/numeric/odeint.hpp>
#include <memory>
#include <set>
namespace generic::ckt {

using namespace Eigen;
class CircuitV
{
public:
    CircuitV(size_t nodes, std::set<size_t> sources, std::set<size_t> probs)
     : m_nodes(nodes)
     , m_probs(std::move(probs))
     , m_sources(std::move(sources))
    {
        m_G = MatrixXd::Zero(Nodes() + SNodes(), Nodes() + SNodes());
        m_C = MatrixXd::Zero(Nodes() + SNodes(), Nodes() + SNodes());
        m_B = MatrixXd::Zero(Nodes() + SNodes(), SNodes());
        m_L = MatrixXd::Zero(Nodes() + SNodes(), PNodes());
        
        size_t i = 0, j = 0;
        for (auto n : m_sources) {
            B()(i + Nodes(), i) = -1;
            mna::StampI(G(), n, i + Nodes());  ++i;
        }
        for (auto p : m_probs) {
            L()(p, j++) = 1;
        }
    }

    virtual ~CircuitV() = default;

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
    
    MatrixXd & G() { return m_G; } 
    MatrixXd & C() { return m_C; } 
    MatrixXd & B() { return m_B; } 
    MatrixXd & L() { return m_L; } 

    const MatrixXd & G() const { return m_G; } 
    const MatrixXd & C() const { return m_C; } 
    const MatrixXd & B() const { return m_B; } 
    const MatrixXd & L() const { return m_L; } 

    size_t Nodes() const { return m_nodes; }
    size_t SNodes() const { return Sources().size(); }
    size_t PNodes() const { return Probs().size(); }
    const std::set<size_t> & Probs() const { return m_probs; }
    const std::set<size_t> & Sources() const { return m_sources; }

private:
    size_t m_nodes;
    std::set<size_t> m_probs;
    std::set<size_t> m_sources;
    MatrixXd m_G, m_C, m_B, m_L;
};

using StateType = std::vector<double>;
struct Intermidiate
{
    VectorXd u;
    MatrixXd coeff;
    MatrixXd input;
    const CircuitV & ckt;
    MatrixXd iG, iC, iB, iL;
    Eigen::PermutationMatrix<Dynamic, Dynamic, size_t> permut;
    Intermidiate(const CircuitV & ckt, bool verbose = false)
     : ckt(ckt)
    {
        u.resize(ckt.Sources().size());
        std::tie(iG, iC, iB, iL, permut) = mna::RegularizeSuDynamic(ckt.G(), ckt.C(), ckt.B(), ckt.L(), verbose);

        auto dcomp = iC.ldlt();
        coeff = dcomp.solve(-1.0 * iG);
        input = dcomp.solve(iB);
        u.resize(input.cols());

        if (verbose) {
            std::cout << "G:\n" << ckt.G() << std::endl;
            std::cout << "C:\n" << ckt.C() << std::endl;
            std::cout << "B:\n" << ckt.B() << std::endl;
            std::cout << "L:\n" << ckt.L() << std::endl;

            std::cout << "iG:\n" << iG << std::endl;
            std::cout << "iC:\n" << iC << std::endl;
            std::cout << "iB:\n" << iB << std::endl;
            std::cout << "iL:\n" << iL << std::endl;

            std::cout << "coeff:\n" << coeff << std::endl;
            std::cout << "input:\n" << input << std::endl;

            std::cout << "p:\n" << permut.indices() << std::endl;
        }
    }
    void State2Output(const StateType & x, StateType & result) const
    {
        result.resize(iL.cols());
        Map<const VectorXd> xvec(x.data(), x.size());
        Map<VectorXd> ovec(result.data(), result.size());
        ovec = iL.transpose() * xvec;
    }

    StateType State2Output(const StateType & x) const
    {
        StateType results;
        State2Output(x, results);
        return results;
    }

    size_t StateSize() const
    {
        return iG.cols();
    }
};

template <typename VoltFunc>
struct Simulator
{
    Intermidiate & im;
    const VoltFunc & vf;
    Simulator(Intermidiate & im, const VoltFunc & vf) : im(im), vf(vf) {}

    void operator() (const StateType & x, StateType & dxdt, double t)
    {   
        for (size_t i = 0; i < im.ckt.SNodes(); ++i) im.u(i) = vf(i, t);
        Map<const VectorXd> xvec(x.data(), x.size());
        Map<VectorXd> result(dxdt.data(), dxdt.size());
        result = im.coeff * xvec + im.input * im.u;  
    }
};

} // namespace generic::ckt

#endif // GENERIC_EIGEN_LIBRARY_SUPPORT