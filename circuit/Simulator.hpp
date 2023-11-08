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
    MNA<MatrixType> m;
    std::set<size_t> prob;
    std::set<size_t> source;
    DenseCircuit(size_t nodes, std::set<size_t> source, std::set<size_t> prob)
     : nodes(nodes), prob(prob), source(source)
    {
        m.G = MatrixType::Zero(nodes + source.size(), nodes + source.size());
        m.C = MatrixType::Zero(nodes + source.size(), nodes + source.size());
        m.B = MatrixType::Zero(nodes + source.size(), source.size());
        m.L = MatrixType::Zero(nodes + source.size(), prob.size());
    
        size_t i = 0, j = 0;
        for (auto n : source) {
            m.B(i + nodes, i) = -1;
            mna::StampI(m.G, n, i + nodes);  ++i;
        }
        for (auto p : prob) m.L(p, j++) = 1;
    }

    void SetR(size_t i, size_t j, Float r) { mna::Stamp(m.G, i, j, 1 / r); }
    void SetC(size_t i, size_t j, Float c) { mna::Stamp(m.C, i, j, c); }
    void SetC(size_t i, Float c) { mna::Stamp(m.C, i, c); }
    const MNA<MatrixType> & Build() const { return m; }
};

template <typename Float>
struct SparseCircuit
{
    using MatrixType = Eigen::SparseMatrix<Float>;
    using Triplets = std::vector<Eigen::Triplet<Float> >;
    size_t nodes;
    std::set<size_t> prob;
    std::set<size_t> source;
    mutable MNA<MatrixType> m;
    mutable Triplets tC, tG, tB, tL;
    SparseCircuit(size_t nodes, std::set<size_t> source, std::set<size_t> prob)
     : nodes(nodes), prob(prob), source(source)
    {
    }

    void SetR(size_t i, size_t j, Float r) { mna::Stamp(tG, i, j, 1 / r); }
    void SetC(size_t i, size_t j, Float c) { mna::Stamp(tC, i, j, c); }
    void SetC(size_t i, Float c) { mna::Stamp(tC, i, c); }

    const MNA<MatrixType> & Build() const
    {
        m.G = MatrixType(nodes + source.size(), nodes + source.size());
        m.C = MatrixType(nodes + source.size(), nodes + source.size());
        m.B = MatrixType(nodes + source.size(), source.size());
        m.L = MatrixType(nodes + source.size(), prob.size());
    
        size_t i = 0, j = 0;
        for (auto n : source) {
            tB.emplace_back(i + nodes, i, -1);
            mna::StampI(tG, n, i + nodes);  ++i;
        }
        for (auto p : prob)
            tL.emplace_back(p, j++, 1);

        m.C.setFromTriplets(tC.begin(), tC.end());
        m.G.setFromTriplets(tG.begin(), tG.end());
        m.B.setFromTriplets(tB.begin(), tB.end());
        m.L.setFromTriplets(tL.begin(), tL.end());
        tC.clear();
        tG.clear();
        tB.clear();
        tL.clear();
        return m;
    }
};

template <typename Float>
struct Intermidiate
{
    using StateType = std::vector<Float>;
    using VectorType = typename DenseCircuit<Float>::VectorType;
    using MatrixType = typename DenseCircuit<Float>::MatrixType;
    using PermutationMatrixType = Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, size_t>;
    
    VectorType u;
    MatrixType iLT;
    size_t stateSize{0};
    MatrixType coeff, input;
    const MNA<MatrixType> & m;
    Intermidiate(const MNA<MatrixType> & m, bool verbose = false)
     : m(m)
    {
        u.resize(m.B.cols());
        auto [iG, iC, iB, iL] = mna::RegularizeSuDynamic(m, verbose);
        auto dcomp = iC.ldlt();
        coeff = dcomp.solve(-1 * iG);
        input = dcomp.solve(iB);
        u.resize(input.cols());
        stateSize = iG.cols();
        iLT = iL.transpose();
        if (verbose) {
            std::cout << "G:\n" << m.G << std::endl;
            std::cout << "C:\n" << m.C << std::endl;
            std::cout << "B:\n" << m.B << std::endl;
            std::cout << "L:\n" << m.L << std::endl;
            std::cout << "iG:\n" << iG << std::endl;
            std::cout << "iC:\n" << iC << std::endl;
            std::cout << "iB:\n" << iB << std::endl;
            std::cout << "iL:\n" << iL << std::endl;
            std::cout << "coeff:\n" << coeff << std::endl;
            std::cout << "input:\n" << input << std::endl;
        }
    }

    void State2Output(const StateType & x, StateType & result) const
    {
        result.resize(iLT.rows());
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

    size_t StateSize() const { return stateSize; }
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
        for (size_t i = 0; i < im.m.B.cols(); ++i) im.u(i) = vf(i, t);
        Eigen::Map<VectorType> result(dxdt.data(), dxdt.size());
        Eigen::Map<const VectorType> xvec(x.data(), x.size());
        result = im.coeff * xvec + im.input * im.u;  
    }
};

} // namespace generic::ckt
