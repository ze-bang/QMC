// SPDX-License-Identifier: MIT
#include <cmath>
#include <vector>

#include "qmc/linalg.hpp"
#include "qmc/rng.hpp"
#include "test_runner.hpp"

using qmc::la::Matrix;
using qmc::Pcg32;

namespace {

Matrix random_matrix(int n, Pcg32& rng) {
    Matrix M(n, n, 0.0);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            M(i, j) = 2.0 * rng.uniform() - 1.0;
    return M;
}

Matrix random_symmetric(int n, Pcg32& rng) {
    Matrix M(n, n, 0.0);
    for (int i = 0; i < n; ++i) {
        M(i, i) = 2.0 * rng.uniform() - 1.0;
        for (int j = i + 1; j < n; ++j) {
            const double v = 2.0 * rng.uniform() - 1.0;
            M(i, j) = v; M(j, i) = v;
        }
    }
    return M;
}

double max_abs_diff(const Matrix& A, const Matrix& B) {
    double d = 0.0;
    for (int i = 0; i < A.rows(); ++i)
        for (int j = 0; j < A.cols(); ++j)
            d = std::max(d, std::fabs(A(i, j) - B(i, j)));
    return d;
}

} // namespace

QMC_TEST(linalg_gemm_identity) {
    Pcg32 rng(11);
    auto A = random_matrix(8, rng);
    auto I = Matrix::identity(8);
    auto AI = qmc::la::matmul(A, I);
    auto IA = qmc::la::matmul(I, A);
    QMC_REQUIRE(max_abs_diff(A, AI) < 1e-12);
    QMC_REQUIRE(max_abs_diff(A, IA) < 1e-12);
}

QMC_TEST(linalg_inverse_round_trip) {
    Pcg32 rng(17);
    for (int n : {2, 4, 8, 16}) {
        // Random + lambda * I to ensure non-singularity.
        auto A = random_matrix(n, rng);
        for (int i = 0; i < n; ++i) A(i, i) += 5.0;
        auto Ainv = qmc::la::inverse(A);
        auto AAi  = qmc::la::matmul(A, Ainv);
        QMC_REQUIRE(max_abs_diff(AAi, Matrix::identity(n)) < 1e-9);
    }
}

QMC_TEST(linalg_logdet_matches_known_value) {
    Matrix A(3, 3, 0.0);
    // Construct A with known det = 6.
    A(0, 0) = 2; A(0, 1) = 0; A(0, 2) = 0;
    A(1, 0) = 1; A(1, 1) = 3; A(1, 2) = 0;
    A(2, 0) = 4; A(2, 1) = 5; A(2, 2) = 1;
    const double d = qmc::la::det(A);
    QMC_REQUIRE_NEAR(d, 6.0, 1e-12);
}

QMC_TEST(linalg_qr_orthogonal_and_recovers_A) {
    Pcg32 rng(3);
    for (int n : {3, 8, 16}) {
        auto A  = random_matrix(n, rng);
        auto qr = qmc::la::qr_factor(A);
        auto QtQ = qmc::la::matmul(
            // Q^T
            [&] {
                Matrix QT(n, n, 0.0);
                for (int i = 0; i < n; ++i)
                    for (int j = 0; j < n; ++j) QT(i, j) = qr.Q(j, i);
                return QT;
            }(),
            qr.Q);
        QMC_REQUIRE(max_abs_diff(QtQ, Matrix::identity(n)) < 1e-9);
        auto QR = qmc::la::matmul(qr.Q, qr.R);
        QMC_REQUIRE(max_abs_diff(QR, A) < 1e-9);
    }
}

QMC_TEST(linalg_sym_eig_diagonalizes) {
    Pcg32 rng(33);
    const int n = 6;
    auto A = random_symmetric(n, rng);
    auto e = qmc::la::sym_eig(A);
    // Reconstruct A = V diag(w) V^T and compare.
    Matrix VD = e.V;
    qmc::la::diag_mul_right(VD, e.w);
    Matrix Ar(n, n, 0.0);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            double s = 0.0;
            for (int k = 0; k < n; ++k) s += VD(i, k) * e.V(j, k);
            Ar(i, j) = s;
        }
    QMC_REQUIRE(max_abs_diff(Ar, A) < 1e-8);
}

QMC_TEST(linalg_expm_sym_check_known_2x2) {
    // K = [[0, 1], [1, 0]] has eigenvalues +/- 1, so
    //   exp(-dt K) = [[cosh(dt), -sinh(dt)], [-sinh(dt), cosh(dt)]].
    Matrix K(2, 2, 0.0);
    K(0, 1) = 1.0; K(1, 0) = 1.0;
    const double dt = 0.7;
    auto E = qmc::la::expm_sym(K, dt);
    QMC_REQUIRE_NEAR(E(0, 0), std::cosh(dt), 1e-12);
    QMC_REQUIRE_NEAR(E(0, 1), -std::sinh(dt), 1e-12);
    QMC_REQUIRE_NEAR(E(1, 0), -std::sinh(dt), 1e-12);
    QMC_REQUIRE_NEAR(E(1, 1), std::cosh(dt), 1e-12);
}

QMC_TEST(linalg_stable_green_matches_naive) {
    Pcg32 rng(7);
    for (int n : {4, 8, 12}) {
        // Random "small" A so that I + A is well-conditioned.
        auto A = random_matrix(n, rng);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) A(i, j) *= 0.3;
        auto G_stable = qmc::la::stable_green(A);
        // Naive reference: (I + A)^{-1}.
        Matrix IpA = A;
        for (int i = 0; i < n; ++i) IpA(i, i) += 1.0;
        auto G_ref = qmc::la::inverse(IpA);
        QMC_REQUIRE(max_abs_diff(G_stable, G_ref) < 1e-8);
    }
}
