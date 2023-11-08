/**
 * @file MNA.hpp
 * @author bwu
 * @brief fork MNA implementation from jefftrull
 * @version 0.1
 * @date 2023-07-16
 */
#pragma once

#include "generic/common/Exception.hpp"
#include "generic/common/Macros.hpp"

#include <Eigen/Sparse>
#include <Eigen/Dense>

#include <iostream>
#include <numeric>
namespace generic {

namespace ckt::mna {

// Functions for implementing MNA with an Eigen matrix
template<typename M, typename Float>
inline void Stamp(M & matrix, size_t i, size_t j, Float g)
{
   // General stamp: conductance at [i,i] and [j,j],
   // -conductance at [i,j] and [j,i], summed with existing values
   matrix(i, i) += g;
   matrix(j, j) += g;
   matrix(i, j) -= g;
   matrix(j, i) -= g;

}

// for when the other end of the device is at GND
template<typename M, typename Float>
inline void Stamp(M & matrix, size_t i, Float g)
{
   matrix(i, i) += g;
}

// for voltage sources (inputs)
template<typename M>
inline void StampI(M & matrix, size_t vnodeno, size_t istateno)
{
   // just basically marks the connection between the inductor (or voltage source)
   // and the voltage, because unlike capacitance, both are state variables.

   matrix(vnodeno, istateno) = 1;   // current is taken *into* inductor or vsource
   matrix(istateno, vnodeno) = -1;
}

// Overloads for when the goal is to create an Eigen SparseMatrix
// In this case we build up lists of (row, col, value) triplets to be set all at once
template<typename Float>
inline void Stamp(typename std::vector<Eigen::Triplet<Float> >& tlist, size_t i, size_t j, Float g)
{
   // General stamp: conductance at [i,i] and [j,j],
   // -conductance at [i,j] and [j,i]
   // Eigen takes care of summing these for us
   tlist.emplace_back(i, i, g);
   tlist.emplace_back(j, j, g);
   tlist.emplace_back(i, j, -g);
   tlist.emplace_back(j, i, -g);

}

template<typename Float>
inline void Stamp(typename std::vector<Eigen::Triplet<Float> >& tlist, size_t i, Float g)
{
   tlist.emplace_back(i, i, g);
}

template<typename Float>
inline void StampI(typename std::vector<Eigen::Triplet<Float> >& tlist, size_t vnodeno, size_t istateno)
{
   tlist.emplace_back(vnodeno, istateno, 1);
   tlist.emplace_back(istateno, vnodeno, -1);
}


template<typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
inline bool IsSingular(const Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols> & m) {
   // A singular matrix has at least one zero eigenvalue -
   // in theory, at least... but due to machine precision we can have "nearly singular"
   // matrices that misbehave.  Comparing rank instead is safer because it uses thresholds
   // for near-zero values.

   GENERIC_ASSERT(m.rows() == m.cols());   // singularity has no meaning for a non-square matrix
   return (m.fullPivLu().rank() != m.rows());
}

template<class Derived, unsigned Mode>
inline bool IsSingular(const Eigen::TriangularView<Derived, Mode> & m) {
   // a "triangular view" is singular if any diagonal element is 0
   // we could use "diagonal" to do this reduction if it were a proper Matrix
   for (int i=0; i < m.rows(); ++i) {
      if (m.coeff(i, i) == 0.0) {
         return true;
      }
   }
   return false;
}

// Convenience template for using Eigen's special allocator with vectors
template<typename Float, int nrows, int ncols>
using MatrixVector = std::vector<Eigen::Matrix<Float, nrows, ncols>, Eigen::aligned_allocator<Eigen::Matrix<Float, nrows, ncols> > >;

// Calculate moments of given system in MNA form
template<int icount, int ocount, int scount, typename Float>
inline MatrixVector<Float, ocount, icount>
Moments(Eigen::Matrix<Float, scount, scount> const & G,
        Eigen::Matrix<Float, scount, scount> const & C,
        Eigen::Matrix<Float, scount, icount> const & B,
        Eigen::Matrix<Float, scount, ocount> const & L,
        Eigen::Matrix<Float, ocount, icount> const & E,
        size_t count) {
    MatrixVector<Float, ocount, icount> result;

    GENERIC_ASSERT(!IsSingular(G));
    auto G_QR = G.fullPivHouseholderQr();
    Eigen::Matrix<Float, scount, scount> A = -G_QR.solve(C);
    Eigen::Matrix<Float, scount, icount> R = G_QR.solve(B);

    result.push_back(L.transpose() * R + E);   // incorporate feedthrough into first moment
    Eigen::Matrix<Float, scount, scount> AtotheI = A;
    for (size_t i = 1; i < count; ++i) {
        result.push_back(L.transpose() * AtotheI * R);
        AtotheI = A * AtotheI;
    }

    return result;
}

// take a circuit's linear system description in G, C, B, L form and compress it so
// the resulting C array is non-singular.  Operation depends on runtime data, so
// output array dimensions are "Dynamic"
template<int icount, int ocount, int scount, typename Float>
inline std::tuple<
        Eigen::Matrix<Float, Eigen::Dynamic, Eigen::Dynamic>,   // G result
        Eigen::Matrix<Float, Eigen::Dynamic, Eigen::Dynamic>,   // C result
        Eigen::Matrix<Float, Eigen::Dynamic, icount>,    // B result
        Eigen::Matrix<Float, Eigen::Dynamic, ocount> >   // L result
RegularizeSu(
        Eigen::Matrix<Float, scount, scount> const & G,
        Eigen::Matrix<Float, scount, scount> const & C,
        Eigen::Matrix<Float, scount, icount> const & B,
        Eigen::Matrix<Float, scount, ocount> const & L) {

    // Use the techniques described in Su (Proc 15th ASP-DAC, 2002) to reduce
    // this set of equations so the state variable derivatives have coefficients
    // Otherwise we cannot integrate to get the time domain result...

    // Use Eigen reductions to find zero rows
    auto zero_rows = (C.array() == 0.0).rowwise().all();   // per row "all zeros"
    size_t zero_count = zero_rows.count();
    size_t state_count = static_cast<size_t>(C.rows());
    size_t nonzero_count = state_count - zero_count;

    // 1. Generate permutation matrix to move zero rows to the bottom
    Eigen::PermutationMatrix<scount, scount, size_t> permut;
    permut.setIdentity(state_count);      // start with null permutation
    size_t i, j;
    for (i = 0, j=(state_count-1); i < j;) {
        // loop invariant: rows > j are all zero; rows < i are not
        while ((i < state_count) && !zero_rows(i)) ++i;
        while ((j > 0) && zero_rows(j)) --j;
        if (i < j) {
            // exchange rows i and j via the permutation vector
            permut.applyTranspositionOnTheRight(i, j);
            ++i; --j;
        }
    }

    // 2. Apply permutation to MNA matrices
    using EqnMatrixType = Eigen::Matrix<Float, scount, scount>;
    EqnMatrixType CP = permut * C * permut;          // permute rows and columns
    EqnMatrixType GP = permut * G * permut;
    Eigen::Matrix<Float, Eigen::Dynamic, icount> BP = permut * B; // permute only rows
    Eigen::Matrix<Float, Eigen::Dynamic, ocount> LP = permut * L;
    
    // 3. Produce reduced equations following Su (Proc. 15th ASP-DAC, 2002)s
    using MatrixD = Eigen::Matrix<Float, Eigen::Dynamic, Eigen::Dynamic>;
    auto G11 = GP.topLeftCorner(nonzero_count, nonzero_count);
    auto G12 = GP.topRightCorner(nonzero_count, zero_count);
    MatrixD G21 = GP.bottomLeftCorner(zero_count, nonzero_count);
    MatrixD G22 = GP.bottomRightCorner(zero_count, zero_count);

    auto L1 = LP.topRows(nonzero_count);
    auto L2 = LP.bottomRows(zero_count);

    auto B1 = BP.topRows(nonzero_count);
    auto B2 = BP.bottomRows(zero_count);

    MatrixD Cred = CP.topLeftCorner(nonzero_count, nonzero_count);

    GENERIC_ASSERT(!IsSingular(G22));
    auto G22QR = G22.fullPivLu();

    MatrixD G22invG21 = G22QR.solve(G21);
    auto G22invB2 = G22QR.solve(B2);
    MatrixD Gred = G11 - G12 * G22invG21;

    Eigen::Matrix<Float, Eigen::Dynamic, ocount> Lred = (L1.transpose() - L2.transpose() * G22invG21).transpose();
    Eigen::Matrix<Float, Eigen::Dynamic, icount> Bred = B1 - G12 * G22invB2;

    // This approach presumes no feedthrough (input-to-output) term
    // MatrixD D = L2.transpose() * G22invB2;
    // GENERIC_ASSERT(D.isZero());

    return std::make_tuple(Gred, Cred, Bred, Lred);
}

template <typename Float>
inline std::tuple<
Eigen::Matrix<Float, Eigen::Dynamic, Eigen::Dynamic>,//G
Eigen::Matrix<Float, Eigen::Dynamic, Eigen::Dynamic>,//C
Eigen::Matrix<Float, Eigen::Dynamic, Eigen::Dynamic>,//B
Eigen::Matrix<Float, Eigen::Dynamic, Eigen::Dynamic>,//L
Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, size_t> >//P
RegularizeSuDynamic(const Eigen::Matrix<Float, Eigen::Dynamic, Eigen::Dynamic> & G, 
                    const Eigen::Matrix<Float, Eigen::Dynamic, Eigen::Dynamic> & C,
                    const Eigen::Matrix<Float, Eigen::Dynamic, Eigen::Dynamic> & B,
                    const Eigen::Matrix<Float, Eigen::Dynamic, Eigen::Dynamic> & L, bool verbose = false)
{
    // Use the techniques described in Su (Proc 15th ASP-DAC, 2002) to reduce
    // this set of equations so the state variable derivatives have coefficients
    // Otherwise we cannot integrate to get the time domain result...

    const size_t scount = G.rows();
    [[maybe_unused]] const size_t icount = B.cols();
    [[maybe_unused]] const size_t ocount = L.cols();

    // Use Eigen reductions to find zero rows
    auto zero_rows = (C.array() == 0.0).rowwise().all();   // per row "all zeros"
    size_t zero_count = zero_rows.count();
    size_t state_count = static_cast<size_t>(C.rows());
    size_t nonzero_count = state_count - zero_count;

    // 1. Generate permutation matrix to move zero rows to the bottom
    Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, size_t> permut;
    permut.resize(scount);
    permut.setIdentity(state_count);      // start with null permutation
    size_t i, j;
    for (i = 0, j=(state_count-1); i < j;) {
        // loop invariant: rows > j are all zero; rows < i are not
        while ((i < state_count) && !zero_rows(i)) ++i;
        while ((j > 0) && zero_rows(j)) --j;
        if (i < j) {
            // exchange rows i and j via the permutation vector
            permut.applyTranspositionOnTheRight(i, j);
            ++i; --j;
        }
    }

    // 2. Apply permutation to MNA matrices
    auto CP = permut * C * permut;          // permute rows and columns
    auto GP = permut * G * permut;
    auto BP = permut * B; // permute only rows
    auto LP = permut * L;

    if (verbose) {
        std::cout << "CP:\n " << CP << std::endl;
        std::cout << "GP:\n " << GP << std::endl;
        std::cout << "BP:\n " << BP << std::endl;
        std::cout << "LP:\n " << LP << std::endl;
    }
    
    // 3. Produce reduced equations following Su (Proc. 15th ASP-DAC, 2002)
    auto G11 = GP.topLeftCorner(nonzero_count, nonzero_count);
    auto G12 = GP.topRightCorner(nonzero_count, zero_count);
    auto G21 = GP.bottomLeftCorner(zero_count, nonzero_count);
    auto G22 = GP.bottomRightCorner(zero_count, zero_count);

    if (verbose) {
        std::cout << "G11:\n " << G11 << std::endl;
        std::cout << "G12:\n " << G12 << std::endl;
        std::cout << "G21:\n " << G21 << std::endl;
        std::cout << "G22:\n " << G22 << std::endl;
    }

    auto L1 = LP.topRows(nonzero_count);
    auto L2 = LP.bottomRows(zero_count);    
    auto B1 = BP.topRows(nonzero_count);
    auto B2 = BP.bottomRows(zero_count);

    if (verbose) {
        std::cout << "L1:\n " << L1 << std::endl;
        std::cout << "L2:\n " << L2 << std::endl;
        std::cout << "B1:\n " << B1 << std::endl;
        std::cout << "B2:\n " << B2 << std::endl;  
    }

    auto Cred = CP.topLeftCorner(nonzero_count, nonzero_count);

    GENERIC_ASSERT(not IsSingular(G22));
    auto G22QR = G22.fullPivLu();
    auto G22invG21 = G22QR.solve(G21);
    auto G22invB2 = G22QR.solve(B2);
    auto Gred = G11 - G12 * G22invG21;

    auto Lred = (L1.transpose() - L2.transpose() * G22invG21).transpose();
    auto Bred = B1 - G12 * G22invB2;

    if (verbose) {
        std::cout << "Cred:\n " << Cred << std::endl;
        std::cout << "Gred:\n " << Gred << std::endl;
        std::cout << "Bred:\n " << Bred << std::endl;
        std::cout << "Lred:\n " << Lred << std::endl;
    }

    // This approach presumes no feedthrough (input-to-output) term
    // auto D = L2.transpose() * G22invB2;
    // GENERIC_ASSERT(D.isZero());

    return std::make_tuple(Gred, Cred, Bred, Lred, permut);
}

namespace detail {

// Implementation of Natarajan regularization
// Each iteration of this process can produce an input derivative term that we subsequently
// absorb into the state variable once the process is complete.  This means potentially
// a series of input derivative coefficients (B's).  We hide that from the users by delegating here:

template<int icount, int ocount, int scount, typename Float = double>
inline std::tuple<
        Eigen::Matrix<Float, Eigen::Dynamic, Eigen::Dynamic>,   // G result
        Eigen::Matrix<Float, Eigen::Dynamic, Eigen::Dynamic>,   // C result
        Eigen::Matrix<Float, Eigen::Dynamic, icount>,    // B result
        Eigen::Matrix<Float, Eigen::Dynamic, ocount>,    // D result
        Eigen::Matrix<Float, ocount, icount> >    // E result (feedthrough)
RegularizeNatarajan(
        const Eigen::Matrix<Float, scount, scount> & G,
        const Eigen::Matrix<Float, scount, scount> & C,
        const MatrixVector<Float, scount, icount>  & B, // in decreasing order of derived-ness
        const Eigen::Matrix<Float, scount, ocount> & D) {

    // Implements the algorithm in [Natarajan]
    // Circuits, Devices and Systems, IEE Proceedings G, June 1991

    using MatrixD = Eigen::Matrix<Float, Eigen::Dynamic, Eigen::Dynamic>;

    // Step 1: put C into "Row Echelon" form by performing LU factorization
    auto lu = C.fullPivLu();
    auto k = lu.rank();
    if (k == C.rows()) {
        // C is already non-singular
        auto E = Eigen::Matrix<Float, ocount, icount>::Zero();
        return std::make_tuple(G, C, B.back(), D, E);
    }

    MatrixD U = lu.matrixLU().template triangularView<Eigen::Upper>();
    auto L    = lu.matrixLU().template triangularView<Eigen::UnitLower>();

    // Step 2: "perform the same elementary operations on G and B"
    // given that C = P.inverse() * L * U * Q.inverse()
    // (from source it seems that permutationP/Q is inverse)
    // then to get the new G we reverse those operations:
    auto & Cprime = U;   // note we may have small non-zero values in bottom rows, but they will be ignored
    auto P = lu.permutationP();
    auto Q = lu.permutationQ();

    MatrixD Gprime = L.solve(P * G * Q);                   // rows and columns
    MatrixVector<Float, scount, icount> Bprime;
    std::transform(B.begin(), B.end(), std::back_inserter(Bprime),
                   [&L, &P](Matrix<Float, scount, icount> const& b) -> Matrix<Float, scount, icount> {
                       return L.solve(P * b);              // rows only
                   });

    // The D input is like L in PRIMA but this algorithm uses the transpose
    auto Dprime = D.transpose() * Q;          // columns only

    // Step 3: "Convert [G21 G22] matrix into row echelon form starting from the last row"

    // decompose the bottom rows
    // The author performed a standard gaussian elimination on the bottom rows of the G
    // matrix, rotated 180 degrees, i.e. as though the matrix were upside down, to get
    // a lower triangular matrix in G22 when the original orientation is restored.
    // Christoph Hertzberg pointed out to me (on IRC) that given pivoting etc. the
    // initial orientation of the matrix is irrelevant and we can simply reorder the
    // necessary bits afterward.  Thus, my new plan:

    // Perform LU decomposition, then exchange G21 and G22, reversing the former, i.e.:
    // X X X  Y Y Y        (G2)       Q Q Q  X 0 0
    // 0 X X  Z Z Z    ============>  Z Z Z  X X 0
    // 0 0 X  Q Q Q                   Y Y Y  X X X
    // In this approach only the original G21 gets column reversed
    // In addition, for performance we should keep these blocks as lazy expressions to
    // avoid reconstructing the full-sized matrix

    auto G2_LU = Gprime.bottomRows(C.rows() - k).fullPivLu();
    // former lower left gets reversed and becomes lower right
    MatrixD G22R = G2_LU.matrixLU().leftCols(C.rows() - k).reverse();
    auto G22   = G22R.template triangularView<Eigen::Lower>();
    // former lower right just gets row permutation to match G22
    auto G21   = G2_LU.matrixLU().rightCols(k).colwise().reverse();  // reverses rows

    // top blocks of G get the column permutations but nothing else
    // both get column permutation from LU
    auto G1  = Gprime.topRows(k) * G2_LU.permutationQ();
    // upper left gets same column reversal as lower right did and becomes G12
    auto G12 = G1.leftCols(C.rows() - k).rowwise().reverse();    // reverses columns
    // upper right just becomes G11
    auto G11 = G1.rightCols(k);

    // Step 4: "Carry out the same row operations in the B matrix"
    // This means the row permutation and multiplication by L from the LU,
    // plus the row reversing permutation we did to make G22 lower triangular.
    // These operations are only applied to the lower part of B
 
    // extract and apply L operation from reversed G2
    auto G2_L = G2_LU.matrixLU().leftCols(C.rows() - k).template triangularView<Eigen::UnitLower>();
    GENERIC_ASSERT(!IsSingular(G2_L));
    MatrixVector<Float, scount, icount> Bnew;
    std::transform(Bprime.begin(), Bprime.end(), std::back_inserter(Bnew),
                   [k, &G2_L, &G2_LU]
                   (Eigen::Matrix<Float, scount, icount> bn) {
                       auto B2 = bn.bottomRows(bn.rows() - k);
                       bn.bottomRows(bn.rows() - k) =
                           G2_L.solve(G2_LU.permutationP() * B2).colwise().reverse();
                       return bn;
                   });

    // Step 5: "Interchange the columns in the G, C, and D matrices... such that G22 is non-singular"
    // Since we have done a full pivot factorization of G2 I assume G22 is already non-singular,
    // so the only thing left to do is reorder the C and D matrices according to the G2 factorization

    // same column permutations as G1 for C1 (only applies to top rows; bottom is zero)
    auto C1  = Cprime.topRows(k) * G2_LU.permutationQ();
    auto C12 = C1.leftCols(C.rows() - k).rowwise().reverse();    // reverses columns
    auto C11 = C1.rightCols(k);

    // D columns also get permuted
    MatrixD pq_view = G2_LU.permutationQ();
    auto D1  = Dprime * G2_LU.permutationQ();
    auto D02 = D1.leftCols(C.rows() - k).rowwise().reverse();
    auto D01 = D1.rightCols(k);

    // Step 6: compute reduced matrices using equations given in paper

    GENERIC_ASSERT(!IsSingular(G22));
    // G22 is a "TriangularView" so has a simple solve method based on back-substitution

    MatrixD Gfinal  = G11 - G12 * G22.solve(G21);
    MatrixD Cfinal  = C11 - C12 * G22.solve(G21);
    Eigen::Matrix<Float, ocount, Eigen::Dynamic> Dfinal
                    = D01 - D02 * G22.solve(G21);

    Eigen::Matrix<Float, Eigen::Dynamic, icount> B02 = Bnew.back().bottomRows(Bnew.back().rows() - k);
    Eigen::Matrix<Float, ocount, icount> E1 = D02 * G22.solve(B02);

    // reduce the entire series of B's to the new size
    // Performing the same substitution as in Natarajan beginning with eqn [5]
    // but with additional input derivatives present.  Adding B11/B12 multiplying a first
    // derivative of Ws demonstrates that each additional input derivative term contributes:
    // Bn1 - G12 * G22^-1 * Bn2  to its own term, and
    //     - C12 * G22^-1 * Bn2  to the derivative n+1 coefficient,
    // once reduced.
    MatrixVector<Float, Eigen::Dynamic, icount> Btrans;
    // n+1's first (equation 9d)
    std::transform(Bnew.begin(), Bnew.end(), std::back_inserter(Btrans),
                   [k, &C12, &G22](Eigen::Matrix<Float, scount, icount> const& Bn) {
                       auto Bn2 = Bn.bottomRows(Bn.rows() - k);
                       return Eigen::Matrix<Float, Eigen::Dynamic, icount>(-C12 * G22.solve(Bn2));
                   });
    Btrans.push_back(Eigen::Matrix<Float, Eigen::Dynamic, icount>::Zero(k, icount));  // contribution from n-1 is 0 (nonexistent)

    // n's next, shifted by one (equation 9c)
    std::transform(Bnew.begin(), Bnew.end(), Btrans.begin()+1, Btrans.begin()+1,
                   [k, &G12, &G22](Eigen::Matrix<Float, scount, icount> const& Bn,
                                   Eigen::Matrix<Float, Eigen::Dynamic, icount> const& Bnm1_contribution)
                   -> Eigen::Matrix<Float, Eigen::Dynamic, icount> {  // without explicitly declared return type Eigen
                                                                    // will keep references to these locals:
                       auto Bn1 = Bn.topRows(k);
                       auto Bn2 = Bn.bottomRows(Bn.rows() - k);

                       return Bn1 - G12 * G22.solve(Bn2) + Bnm1_contribution;
                   });

    // If Cfinal is singular, we need to repeat this analysis on the new matrices
    if (IsSingular(Cfinal)) {
        Eigen::Matrix<Float, Eigen::Dynamic, ocount> Dtrans = Dfinal.transpose();   // no implicit conversion on fn tmpl args
        auto recursiveResult = RegularizeNatarajan<icount, ocount, Eigen::Dynamic>(Gfinal, Cfinal, Btrans, Dtrans);
        return std::make_tuple(std::get<0>(recursiveResult),  // G
                               std::get<1>(recursiveResult),  // C
                               std::get<2>(recursiveResult),  // B
                               std::get<3>(recursiveResult),  // D
                               std::get<4>(recursiveResult) + E1);  // combine E
    }

    // We've found a non-singular Cfinal and a set of B's
    // We need to apply a transformation suggested by Chen (TCAD July 2012) to eliminate
    // all input derivative terms.  Chen gives only the simplest case, for B0 * Ws + B1 * Ws' :
    // Br = B0 - Gr * Cr^-1 * B1
    // based on a variable substitution of:
    // Xnew = X - Cr^-1 * B1 * Ws
    // and mentions the rest should be done "recursively".  I believe the general case is:
    // Br = B0 - Gr * Cr^-1 * (B1 - Gr * Cr^-1 *(B2 - ... ))
    Eigen::Matrix<Float, Eigen::Dynamic, icount> Bfinal = Eigen::Matrix<Float, Eigen::Dynamic, icount>::Zero(k, icount);
    Bfinal = std::accumulate(
        // starting with the first (most derived) coefficient, compute above expression for Br:
        Btrans.begin(), Btrans.end(), Bfinal,
        [Gfinal, Cfinal](Eigen::Matrix<Float, Eigen::Dynamic, icount> const& acc,
                         Eigen::Matrix<Float, Eigen::Dynamic, icount> const& B) -> decltype(Bfinal) {
            return B - Gfinal * Cfinal.fullPivHouseholderQr().solve(acc);
        });

    // The variable substitution for the 2nd derivative case is:
    // Xnew = X - Cr^-1 * (B2 * Ws' - (Gr * Cr^-1 * B2 - B1) * Ws)
    // Making this substitution in the output equation Y = D * X + E * Ws gives
    // Y = D * Xnew + D * Cr^-1 * (B1 - Gr * Cr^-1 * B2) * Ws + Cr^-1 * B2 * Ws'
    // however, if the Ws' term is nonzero the system is ill-formed:
    if (Btrans.size() >= 3) {
        Eigen::Matrix<Float, Eigen::Dynamic, icount> CinvB = Cfinal.fullPivLu().solve(*(Btrans.rbegin()+2));
        GENERIC_ASSERT(CinvB.isZero());
    }

    // now I can calculate the new value for E, which can only be:
    // E = E1 + D * Cr^-1 * B1
    // because, thanks to the GENERIC_ASSERTion, all other terms must be 0
    Eigen::Matrix<Float, ocount, icount> Efinal = E1 + Dfinal * Cfinal.fullPivHouseholderQr().solve(*(Btrans.rbegin()+1));

    return std::make_tuple(Gfinal, Cfinal, Bfinal,
                           Dfinal.transpose(),  // for PRIMA compatibility
                           Efinal);
}
}  // namespace detail

// user-facing function (only one "B" parameter)
template<int icount, int ocount, int scount, typename Float>
inline std::tuple<
        Eigen::Matrix<Float, Eigen::Dynamic, Eigen::Dynamic>,   // G result
        Eigen::Matrix<Float, Eigen::Dynamic, Eigen::Dynamic>,   // C result
        Eigen::Matrix<Float, Eigen::Dynamic, icount>,    // B result
        Eigen::Matrix<Float, Eigen::Dynamic, ocount>,    // D result
        Eigen::Matrix<Float, ocount, icount> >    // E result (feedthrough)
RegularizeNatarajan(
        const Eigen::Matrix<Float, scount, scount> & G,
        const Eigen::Matrix<Float, scount, scount> & C,
        const Eigen::Matrix<Float, scount, icount> & B,
        const Eigen::Matrix<Float, scount, ocount> & D) {
    return detail::RegularizeNatarajan(G, C, MatrixVector<Float, scount, icount>(1, B), D);
}

} // namespace ckt::mna
} // namespace generic