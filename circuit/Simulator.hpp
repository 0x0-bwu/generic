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

template <typename Float, bool Augmented = true>
struct DenseCircuit
{
    using VectorType = Eigen::Matrix<Float, Eigen::Dynamic, 1>;
    using MatrixType = Eigen::Matrix<Float, Eigen::Dynamic, Eigen::Dynamic>;
    Float rd{1}; 
    size_t nodes;
    MNA<MatrixType> m;
    std::vector<size_t> prob;
    std::vector<size_t> source;
    DenseCircuit(size_t nodes, std::vector<size_t> source, std::vector<size_t> prob, Float rd = 1)
     : rd(rd), nodes(nodes), prob(prob), source(source)
    {
        size_t augment = Augmented ? source.size() : 0;
        m.G = MatrixType::Zero(nodes + augment, nodes + augment);
        m.C = MatrixType::Zero(nodes + augment, nodes + augment);
        m.B = MatrixType::Zero(nodes + augment, source.size());
        m.L = MatrixType::Zero(nodes + augment, prob.size());
    
        size_t i = 0, j = 0;
        for (auto n : source) {
            if (Augmented) {
                m.B(i + nodes, i) = -1;
                mna::StampI(m.G, n, i + nodes);  ++i;
            }
            else {
                m.B(n, i) = Float(1 / rd);
                mna::Stamp(m.G, n, Float(1 / rd));
            }
        }
        for (auto p : prob) m.L(p, j++) = 1;
    }

/**
 * @brief Brief description of SetR.
 * @param i
 * @param j
 * @param r
 * @return void
 */
    void SetR(size_t i, size_t j, Float r) { mna::Stamp(m.G, i, j, 1 / r); }
/**
 * @brief Brief description of SetR.
 * @param i
 * @param r
 * @return void
 */
    void SetR(size_t i, Float r) { mna::Stamp(m.G, i, 1 / r); }

/**
 * @brief Brief description of SetC.
 * @param i
 * @param j
 * @param c
 * @return void
 */
    void SetC(size_t i, size_t j, Float c) { mna::Stamp(m.C, i, j, c); }
/**
 * @brief Brief description of SetC.
 * @param i
 * @param c
 * @return void
 */
    void SetC(size_t i, Float c) { mna::Stamp(m.C, i, c); }
    const MNA<MatrixType> & Build() const { return m; }
};

template <typename Float, bool Augmented = true>
struct SparseCircuit
{
    using MatrixType = Eigen::SparseMatrix<Float>;
    using Triplets = std::vector<Eigen::Triplet<Float> >;
    Float rd{1}; 
    size_t nodes;
    std::vector<size_t> prob;
    std::vector<size_t> source;
    mutable MNA<MatrixType> m;
    mutable Triplets tC, tG, tB, tL;
    SparseCircuit(size_t nodes, std::vector<size_t> source, std::vector<size_t> prob, Float rd = 1)
     : rd(rd), nodes(nodes), prob(prob), source(source)
    {
    }

/**
 * @brief Brief description of SetR.
 * @param i
 * @param j
 * @param r
 * @return void
 */
    void SetR(size_t i, size_t j, Float r) { mna::Stamp(tG, i, j, 1 / r); }
/**
 * @brief Brief description of SetC.
 * @param i
 * @param j
 * @param c
 * @return void
 */
    void SetC(size_t i, size_t j, Float c) { mna::Stamp(tC, i, j, c); }
/**
 * @brief Brief description of SetC.
 * @param i
 * @param c
 * @return void
 */
    void SetC(size_t i, Float c) { mna::Stamp(tC, i, c); }

    const MNA<MatrixType> & Build() const
    {
        size_t augment = Augmented ? source.size() : 0;
        m.G = MatrixType(nodes + augment, nodes + augment);
        m.C = MatrixType(nodes + augment, nodes + augment);
        m.B = MatrixType(nodes + augment, source.size());
        m.L = MatrixType(nodes + augment, prob.size());
    
        size_t i = 0, j = 0;
        for (auto n : source) {
            if (Augmented) {
                tB.emplace_back(i + nodes, i, -1);
                mna::StampI(tG, n, i + nodes);  ++i;
            }
            else {
                tB.emplace_back(n, i, Float(1 / rd));
                mna::Stamp(tG, n, Float(1 / rd));
            }
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
        auto [rG, rC, rB, rL] = mna::RegularizeSu(m, p);
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
            std::cout << "P:\n"  << p.indices() << std::endl; 
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
