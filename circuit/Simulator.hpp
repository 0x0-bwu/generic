/**
 * @file Simulator.hpp
 * @author bwu
 * @brief static/transient circuit simulator
 * @version 0.1
 * @date 2023-07-16
 */
#pragma once

#include "MNA.hpp"
#include <boost/numeric/odeint.hpp>
#include <set>
namespace generic::ckt {

template <typename Float>
struct DenseCircuit
{
    using VectorType = Eigen::Matrix<Float, Eigen::Dynamic, 1>;
    using MatrixType = Eigen::Matrix<Float, Eigen::Dynamic, Eigen::Dynamic>;
    size_t nodes;
    MatrixType G, C, B, L;
    std::set<size_t> prob;
    std::set<size_t> source;
    DenseCircuit(size_t nodes, std::set<size_t> source, std::set<size_t> prob)
     : nodes(nodes), prob(prob), source(source)
    {
        G = MatrixType::Zero(nodes + source.size(), nodes + source.size());
        C = MatrixType::Zero(nodes + source.size(), nodes + source.size());
        B = MatrixType::Zero(nodes + source.size(), source.size());
        L = MatrixType::Zero(nodes + source.size(), prob.size());
    
        size_t i = 0, j = 0;
        for (auto n : source) {
            B(i + nodes, i) = -1;
            mna::StampI(G, n, i + nodes);  ++i;
        }
        for (auto p : prob) L(p, j++) = 1;
    }

    void SetR(size_t i, size_t j, Float r) { mna::Stamp(G, i, j, 1 / r); }
    void SetC(size_t i, size_t j, Float c) { mna::Stamp(C, i, j, c); }
    void SetC(size_t i, Float c) { mna::Stamp(C, i, c); }
};

template <typename Float>
struct Intermidiate
{
    using StateType = std::vector<Float>;
    using VectorType = typename DenseCircuit<Float>::VectorType;
    using MatrixType = typename DenseCircuit<Float>::MatrixType;
    using PermutationMatrixType = Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, size_t>;
    VectorType u;
    MatrixType coeff, input;
    PermutationMatrixType permut;
    MatrixType iG, iC, iB, iL, iLT;
    const DenseCircuit<Float> & ckt;
    Intermidiate(const DenseCircuit<Float> & ckt, bool verbose = false)
     : ckt(ckt)
    {
        u.resize(ckt.source.size());
        std::tie(iG, iC, iB, iL, permut) = mna::RegularizeSuDynamic(ckt.G, ckt.C, ckt.B, ckt.L, verbose);
        iLT = iL.transpose();
        auto dcomp = iC.ldlt();
        coeff = dcomp.solve(-1 * iG);
        input = dcomp.solve(iB);
        u.resize(input.cols());
        if (verbose) {
            std::cout << "G:\n" << ckt.G << std::endl;
            std::cout << "C:\n" << ckt.C << std::endl;
            std::cout << "B:\n" << ckt.B << std::endl;
            std::cout << "L:\n" << ckt.L << std::endl;
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
        Eigen::Map<const VectorType> xvec(x.data(), x.size());
        Eigen::Map<VectorType> ovec(result.data(), result.size());
        ovec = iLT * xvec;
    }

    StateType State2Output(const StateType & x) const
    {
        StateType results;
        State2Output(x, results);
        return results;
    }

    size_t StateSize() const { return iG.cols(); }
};

template <typename Float, typename VoltFunc>
struct Simulator
{
    const VoltFunc & vf;
    Intermidiate<Float> & im;
    using StateType = typename Intermidiate<Float>::StateType;
    Simulator(Intermidiate<Float> & im, const VoltFunc & vf) : vf(vf), im(im) {}

    void operator() (const StateType & x, StateType & dxdt, double t)
    {   
        using VectorType = typename  Intermidiate<Float>::VectorType;
        for (size_t i = 0; i < im.ckt.source.size(); ++i) im.u(i) = vf(i, t);
        Eigen::Map<VectorType> result(dxdt.data(), dxdt.size());
        Eigen::Map<const VectorType> xvec(x.data(), x.size());
        result = im.coeff * xvec + im.input * im.u;  
    }
};

} // namespace generic::ckt
