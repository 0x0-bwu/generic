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
    MatrixType rLT;
    PermutMatrix p;
    size_t stateSize{0};
    MatrixType coeff, input;
    const MNA<MatrixType> & m;
    Intermidiate(const MNA<MatrixType> & m, bool verbose = false)
     : m(m)
    {
        u.resize(m.B.cols());
        auto [rG, rC, rB, rL] = mna::RegularizeSuDynamic(m, p);
        auto dcomp = rC.ldlt();
        coeff = dcomp.solve(-1 * rG);
        input = dcomp.solve(rB);
        u.resize(input.cols());
        stateSize = rG.cols();
        rLT = rL.transpose();
        if (verbose) {
            std::cout << "G:\n" << m.G << std::endl;
            std::cout << "C:\n" << m.C << std::endl;
            std::cout << "B:\n" << m.B << std::endl;
            std::cout << "L:\n" << m.L << std::endl;
            std::cout << "rG:\n" << rG << std::endl;
            std::cout << "rC:\n" << rC << std::endl;
            std::cout << "rB:\n" << rB << std::endl;
            std::cout << "rL:\n" << rL << std::endl;
            std::cout << "coeff:\n" << coeff << std::endl;
            std::cout << "input:\n" << input << std::endl;
        }
    }

    StateType InitState(const StateType & input) const
    {
        GENERIC_ASSERT(input.size() == static_cast<size_t>(m.G.cols()));
        StateType state(m.G.cols());
        Eigen::Map<const VectorType> in(input.data(), input.size());
        Eigen::Map<VectorType> out(state.data(), state.size());
        out = p * in;
        state.resize(stateSize);
        return state;
    }

    void State2Output(const StateType & x, StateType & out) const
    {
        out.resize(rLT.rows());
        Eigen::Map<VectorType> ovec(out.data(), out.size());
        Eigen::Map<const VectorType> xvec(x.data(), x.size());
        ovec = rLT * xvec;
    }

    StateType State2Output(const StateType & x) const
    {
        StateType out;
        State2Output(x, out);
        return out;
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

    void operator() (const StateType & x, StateType & dxdt, Float t)
    {   
        using VectorType = typename  Intermidiate<Float>::VectorType;
        for (int i = 0; i < im.m.B.cols(); ++i) im.u(i) = vf(i, t);
        Eigen::Map<VectorType> result(dxdt.data(), dxdt.size());
        Eigen::Map<const VectorType> xvec(x.data(), x.size());
        result = im.coeff * xvec + im.input * im.u;  
    }
};

} // namespace generic::ckt
