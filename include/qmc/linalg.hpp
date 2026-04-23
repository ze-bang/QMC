// SPDX-License-Identifier: MIT
//
// Minimal dense linear algebra used by the DQMC engine.
//
// The DQMC algorithm needs only a handful of operations on small to
// moderate (typically 64x64 ... 256x256) real matrices:
//
//   * matrix-matrix product           (GEMM)
//   * inverse                         (LU with partial pivoting)
//   * log |det A|, sign det A         (LU)
//   * Householder QR                  (numerical stabilization of B-products)
//   * (I + A)^{-1} for general A      (LU on (I + A))
//
// We deliberately do not depend on BLAS / LAPACK / Eigen. Performance
// is fine for all sizes typical of pedagogical DQMC studies (8x8 - 16x16
// lattices, M = 64 - 256). For production-scale studies one should
// link against optimized BLAS+LAPACK; the Matrix wrapper here uses a
// contiguous row-major buffer so this is a one-line change.

#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

#include "qmc/types.hpp"

namespace qmc::la {

// Row-major dense matrix of `Real`. Trivial value semantics, with a
// flat std::vector<Real> backing store for cache locality.
class Matrix {
public:
    Matrix() = default;
    Matrix(int rows, int cols, Real fill = 0.0)
        : rows_(rows), cols_(cols), data_(static_cast<std::size_t>(rows) * cols, fill) {}

    static Matrix identity(int n) {
        Matrix M(n, n, 0.0);
        for (int i = 0; i < n; ++i) M(i, i) = 1.0;
        return M;
    }

    int rows() const noexcept { return rows_; }
    int cols() const noexcept { return cols_; }
    Real*       data()       noexcept { return data_.data(); }
    const Real* data() const noexcept { return data_.data(); }

    Real& operator()(int i, int j) noexcept {
        return data_[static_cast<std::size_t>(i) * cols_ + j];
    }
    Real operator()(int i, int j) const noexcept {
        return data_[static_cast<std::size_t>(i) * cols_ + j];
    }

    void resize(int r, int c, Real fill = 0.0) {
        rows_ = r; cols_ = c;
        data_.assign(static_cast<std::size_t>(r) * c, fill);
    }
    void fill(Real v) { std::fill(data_.begin(), data_.end(), v); }
    void swap(Matrix& o) noexcept {
        std::swap(rows_, o.rows_); std::swap(cols_, o.cols_);
        data_.swap(o.data_);
    }

    Real trace() const {
        const int n = std::min(rows_, cols_);
        Real s = 0.0;
        for (int i = 0; i < n; ++i) s += (*this)(i, i);
        return s;
    }

    Real frobenius_norm() const {
        Real s = 0.0;
        for (Real v : data_) s += v * v;
        return std::sqrt(s);
    }

private:
    int                rows_{0};
    int                cols_{0};
    std::vector<Real>  data_;
};

// ---------------------------------------------------------------------------
// Basic BLAS-style operations
// ---------------------------------------------------------------------------

// C = alpha * A * B + beta * C. Naive triple-loop with the inner two
// loops reordered for unit-stride access to A and C. Good enough for
// the matrix sizes of interest; replace with a BLAS GEMM for speed.
inline void gemm(Real alpha, const Matrix& A, const Matrix& B,
                 Real beta, Matrix& C) {
    const int M = A.rows();
    const int K = A.cols();
    const int N = B.cols();
    if (B.rows() != K || C.rows() != M || C.cols() != N) {
        throw std::invalid_argument("la::gemm: dimension mismatch");
    }
    if (beta == 0.0) {
        std::fill(C.data(), C.data() + static_cast<std::size_t>(M) * N, 0.0);
    } else if (beta != 1.0) {
        const std::size_t total = static_cast<std::size_t>(M) * N;
        for (std::size_t k = 0; k < total; ++k) C.data()[k] *= beta;
    }
    for (int i = 0; i < M; ++i) {
        for (int k = 0; k < K; ++k) {
            const Real aik = alpha * A(i, k);
            if (aik == 0.0) continue;
            const Real* Brow = B.data() + static_cast<std::size_t>(k) * N;
            Real*       Crow = C.data() + static_cast<std::size_t>(i) * N;
            for (int j = 0; j < N; ++j) Crow[j] += aik * Brow[j];
        }
    }
}

inline Matrix matmul(const Matrix& A, const Matrix& B) {
    Matrix C(A.rows(), B.cols(), 0.0);
    gemm(1.0, A, B, 0.0, C);
    return C;
}

// In-place multiply rows by a diagonal: A = D * A.
inline void diag_mul_left(const std::vector<Real>& d, Matrix& A) {
    const int M = A.rows();
    const int N = A.cols();
    if ((int)d.size() != M) throw std::invalid_argument("diag_mul_left: dim mismatch");
    for (int i = 0; i < M; ++i) {
        const Real di = d[i];
        Real* row = A.data() + static_cast<std::size_t>(i) * N;
        for (int j = 0; j < N; ++j) row[j] *= di;
    }
}

// In-place multiply columns by a diagonal: A = A * D.
inline void diag_mul_right(Matrix& A, const std::vector<Real>& d) {
    const int M = A.rows();
    const int N = A.cols();
    if ((int)d.size() != N) throw std::invalid_argument("diag_mul_right: dim mismatch");
    for (int i = 0; i < M; ++i) {
        Real* row = A.data() + static_cast<std::size_t>(i) * N;
        for (int j = 0; j < N; ++j) row[j] *= d[j];
    }
}

// ---------------------------------------------------------------------------
// LU decomposition with partial pivoting
// ---------------------------------------------------------------------------
//
// Factors A = P * L * U in-place. Returns the row-permutation pivot
// array `piv` and the parity (+1 / -1) of the permutation. Throws on
// exact singularity.

struct LU {
    Matrix             A;       // L (below diag, unit) + U (on/above diag)
    std::vector<int>   piv;     // length n, row permutation
    int                parity;  // sign of det due to permutations
};

inline LU lu_factor(Matrix A) {
    const int n = A.rows();
    if (A.cols() != n) throw std::invalid_argument("lu_factor: matrix must be square");
    std::vector<int> piv(n);
    for (int i = 0; i < n; ++i) piv[i] = i;
    int parity = +1;

    for (int k = 0; k < n; ++k) {
        // Partial pivoting: find row index with the largest |A(i, k)|.
        int p = k;
        Real best = std::fabs(A(k, k));
        for (int i = k + 1; i < n; ++i) {
            const Real v = std::fabs(A(i, k));
            if (v > best) { best = v; p = i; }
        }
        if (best == 0.0) throw std::runtime_error("lu_factor: singular matrix");
        if (p != k) {
            // Swap rows k and p in the dense matrix and the pivot vector.
            for (int j = 0; j < n; ++j) std::swap(A(k, j), A(p, j));
            std::swap(piv[k], piv[p]);
            parity = -parity;
        }
        const Real pivot = A(k, k);
        for (int i = k + 1; i < n; ++i) {
            const Real factor = A(i, k) / pivot;
            A(i, k) = factor;
            for (int j = k + 1; j < n; ++j) {
                A(i, j) -= factor * A(k, j);
            }
        }
    }
    return {std::move(A), std::move(piv), parity};
}

// log |det A| from an LU factorization, with the sign of det.
struct SignedLogDet { Real log_abs_det; int sign; };
inline SignedLogDet logdet(const LU& lu) {
    const int n = lu.A.rows();
    Real ld = 0.0;
    int  sg = lu.parity;
    for (int i = 0; i < n; ++i) {
        const Real d = lu.A(i, i);
        ld += std::log(std::fabs(d));
        if (d < 0.0) sg = -sg;
    }
    return {ld, sg};
}

inline SignedLogDet logdet(const Matrix& A) { return logdet(lu_factor(A)); }

inline Real det(const Matrix& A) {
    const auto sld = logdet(A);
    return static_cast<Real>(sld.sign) * std::exp(sld.log_abs_det);
}

// Solve  A * X = B  via the LU factorization.
inline Matrix lu_solve(const LU& lu, Matrix B) {
    const int n = lu.A.rows();
    if (B.rows() != n) throw std::invalid_argument("lu_solve: dim mismatch");
    const int nrhs = B.cols();
    // Apply row permutation P to B (B := P * B).
    Matrix Bp(n, nrhs, 0.0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < nrhs; ++j) Bp(i, j) = B(lu.piv[i], j);
    }
    // Forward solve L * Y = P * B.
    for (int k = 0; k < n; ++k) {
        for (int i = k + 1; i < n; ++i) {
            const Real lik = lu.A(i, k);
            for (int j = 0; j < nrhs; ++j) Bp(i, j) -= lik * Bp(k, j);
        }
    }
    // Backward solve U * X = Y.
    for (int k = n - 1; k >= 0; --k) {
        const Real ukk = lu.A(k, k);
        for (int j = 0; j < nrhs; ++j) Bp(k, j) /= ukk;
        for (int i = 0; i < k; ++i) {
            const Real uik = lu.A(i, k);
            for (int j = 0; j < nrhs; ++j) Bp(i, j) -= uik * Bp(k, j);
        }
    }
    return Bp;
}

inline Matrix inverse(const Matrix& A) {
    const int n = A.rows();
    if (A.cols() != n) throw std::invalid_argument("inverse: matrix must be square");
    return lu_solve(lu_factor(A), Matrix::identity(n));
}

// ---------------------------------------------------------------------------
// Symmetric eigendecomposition (Jacobi rotations)
// ---------------------------------------------------------------------------
//
// Used to compute  exp(-dt * K)  for a symmetric kinetic matrix K via
// K = Q diag(w) Q^T  ->  exp(-dt K) = Q diag(exp(-dt w)) Q^T.
// Jacobi has poor asymptotic complexity (O(n^3) per sweep, ~log(eps^-1)
// sweeps) but is dependable and trivial to implement. For the lattice
// sizes of interest (n <= 256) it is more than fast enough since the
// kinetic exponential is computed exactly *once* at simulation start.

struct SymEig { std::vector<Real> w; Matrix V; };

inline SymEig sym_eig(Matrix A) {
    const int n = A.rows();
    if (A.cols() != n) throw std::invalid_argument("sym_eig: not square");
    Matrix V = Matrix::identity(n);
    constexpr int    max_sweeps = 100;
    constexpr Real   tol        = 1e-14;

    for (int sweep = 0; sweep < max_sweeps; ++sweep) {
        // Sum of squared off-diagonal entries; converged when small.
        Real off = 0.0;
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j) off += A(i, j) * A(i, j);
        if (off < tol * tol) break;

        for (int p = 0; p < n - 1; ++p) {
            for (int q = p + 1; q < n; ++q) {
                const Real apq = A(p, q);
                if (std::fabs(apq) < 1e-300) continue;
                const Real app = A(p, p);
                const Real aqq = A(q, q);
                const Real theta = 0.5 * (aqq - app) / apq;
                const Real t = (theta >= 0.0)
                    ?  1.0 / (theta + std::sqrt(1.0 + theta * theta))
                    : -1.0 / (-theta + std::sqrt(1.0 + theta * theta));
                const Real c = 1.0 / std::sqrt(1.0 + t * t);
                const Real s = t * c;

                A(p, p) = app - t * apq;
                A(q, q) = aqq + t * apq;
                A(p, q) = 0.0;
                A(q, p) = 0.0;
                for (int i = 0; i < n; ++i) {
                    if (i == p || i == q) continue;
                    const Real aip = A(i, p);
                    const Real aiq = A(i, q);
                    A(i, p) = c * aip - s * aiq;  A(p, i) = A(i, p);
                    A(i, q) = s * aip + c * aiq;  A(q, i) = A(i, q);
                }
                for (int i = 0; i < n; ++i) {
                    const Real vip = V(i, p);
                    const Real viq = V(i, q);
                    V(i, p) = c * vip - s * viq;
                    V(i, q) = s * vip + c * viq;
                }
            }
        }
    }
    std::vector<Real> w(n);
    for (int i = 0; i < n; ++i) w[i] = A(i, i);
    return {std::move(w), std::move(V)};
}

// Build  exp(-dt * K)  for a symmetric K via its eigendecomposition.
inline Matrix expm_sym(const Matrix& K, Real dt) {
    auto eig = sym_eig(K);
    const int n = K.rows();
    std::vector<Real> e(n);
    for (int i = 0; i < n; ++i) e[i] = std::exp(-dt * eig.w[i]);
    // exp(-dt K) = V * diag(e) * V^T.
    Matrix VD = eig.V;
    diag_mul_right(VD, e);
    Matrix R(n, n, 0.0);
    // R = VD * V^T
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            Real s = 0.0;
            for (int k = 0; k < n; ++k) s += VD(i, k) * eig.V(j, k);
            R(i, j) = s;
        }
    }
    return R;
}

// ---------------------------------------------------------------------------
// Householder QR
// ---------------------------------------------------------------------------
//
// Factors A (m x n, m >= n) into  A = Q * R  with Q (m x m) orthogonal
// and R (m x n) upper-triangular. Used in the DQMC stabilization
// scheme to keep the product of B-matrices well-conditioned.

struct QR { Matrix Q; Matrix R; };

inline QR qr_factor(Matrix A) {
    const int m = A.rows();
    const int n = A.cols();
    if (m < n) throw std::invalid_argument("qr_factor: need m >= n");
    Matrix Q = Matrix::identity(m);

    for (int k = 0; k < n; ++k) {
        // Compute the Householder vector for column k below the diag.
        Real sigma = 0.0;
        for (int i = k; i < m; ++i) sigma += A(i, k) * A(i, k);
        if (sigma == 0.0) continue;
        const Real norm_x = std::sqrt(sigma);
        const Real alpha  = (A(k, k) >= 0.0) ? -norm_x : norm_x;
        std::vector<Real> v(m - k, 0.0);
        v[0] = A(k, k) - alpha;
        for (int i = 1; i < m - k; ++i) v[i] = A(i + k, k);
        Real vnorm2 = 0.0;
        for (Real vi : v) vnorm2 += vi * vi;
        if (vnorm2 == 0.0) continue;
        const Real beta = 2.0 / vnorm2;

        // Apply H = I - beta v v^T to A from the left (rows k..m-1).
        for (int j = k; j < n; ++j) {
            Real dot = 0.0;
            for (int i = 0; i < m - k; ++i) dot += v[i] * A(i + k, j);
            const Real f = beta * dot;
            for (int i = 0; i < m - k; ++i) A(i + k, j) -= f * v[i];
        }
        // Accumulate Q := Q * H by applying H from the right to Q.
        for (int i = 0; i < m; ++i) {
            Real dot = 0.0;
            for (int j = 0; j < m - k; ++j) dot += Q(i, j + k) * v[j];
            const Real f = beta * dot;
            for (int j = 0; j < m - k; ++j) Q(i, j + k) -= f * v[j];
        }
    }
    // R := upper-triangular part of A.
    Matrix R(m, n, 0.0);
    for (int i = 0; i < n; ++i)
        for (int j = i; j < n; ++j) R(i, j) = A(i, j);
    return {std::move(Q), std::move(R)};
}

// ---------------------------------------------------------------------------
// Stable evaluation of  G = (I + A)^{-1}
// ---------------------------------------------------------------------------
//
// When A is the product of many ill-conditioned propagator matrices we
// must avoid forming I + A naively (the result loses precision). The
// approach below uses the QR factorization of A: A = Q R, then
//
//     I + A = Q (Q^T + R)        (since Q^T Q = I)
//     G     = (Q^T + R)^{-1} Q^T
//
// (Q^T + R) is upper-Hessenberg-like and is solved by a plain LU. For
// the matrix sizes typical of DQMC studies (n = 64 - 256) this gives
// machine-precision Green's functions for moderate beta.

inline Matrix stable_green(const Matrix& A) {
    const int n = A.rows();
    if (A.cols() != n) throw std::invalid_argument("stable_green: not square");
    auto qr = qr_factor(A);
    Matrix M = qr.R;
    // M += Q^T   (i.e. M(i,j) += Q(j,i))
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) M(i, j) += qr.Q(j, i);
    // RHS  Q^T  (n x n)
    Matrix QT(n, n, 0.0);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) QT(i, j) = qr.Q(j, i);
    return lu_solve(lu_factor(std::move(M)), QT);
}

} // namespace qmc::la
