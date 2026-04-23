// SPDX-License-Identifier: MIT
//
// Determinant Quantum Monte Carlo (DQMC) for the repulsive Hubbard
// model on a bipartite lattice.
//
// Algorithm:    Blankenbecler-Scalapino-Sugar (BSS, 1981) with the
//               Hirsch (1983) discrete Hubbard-Stratonovich
//               transformation that couples to S^z. Single-spin-flip
//               local updates are accepted using rank-1 Sherman-Morrison
//               formulas for the equal-time Green's function. The
//               propagator product is recomputed from scratch every
//               `n_stab` time slices via a QR-based stable formula
//               (see `la::stable_green`).
//
// Conventions:
//   * Trotter step:  beta = L_tau * dt
//   * B-matrix:      B_l_sigma = exp(-dt * K) * exp(V_l_sigma)
//   * Symmetric HS:  cosh(alpha) = exp(dt * U / 2)
//   * Equal-time G:  G_sigma(l)_ij = <c_i c_j^+>     on time slice l,
//                    so G_sigma(l) = (I + A_sigma(l))^{-1}, where
//                    A_sigma(l) = B_l B_{l-1} ... B_1 B_Ltau ... B_{l+1}.
//
// Sign problem:    on a bipartite lattice with mu = 0 (half-filling) the
//                  two spin determinants are equal and Z is sign-free.
//                  Away from half-filling there is in general a sign
//                  problem that grows exponentially with beta U; the
//                  engine still runs and reports the average sign.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include "qmc/hubbard.hpp"
#include "qmc/lattice.hpp"
#include "qmc/linalg.hpp"
#include "qmc/rng.hpp"
#include "qmc/types.hpp"

namespace qmc {

struct DqmcParams {
    Real beta   = 4.0;     // inverse temperature
    int  L_tau  = 40;      // number of Trotter slices  (dt = beta / L_tau)
    int  n_stab = 10;      // recompute G from scratch every n_stab slices
};

struct DqmcStats {
    std::uint64_t flips_proposed{0};
    std::uint64_t flips_accepted{0};
    std::uint64_t sweeps_done{0};
    Real          sign_running_sum{0.0};   // sum of sign(det^2) (always +1 at HF)
    Real          max_recompute_drift{0.0};

    Real acceptance() const {
        return flips_proposed
                   ? static_cast<Real>(flips_accepted) / flips_proposed
                   : 0.0;
    }
    Real average_sign() const {
        return sweeps_done ? sign_running_sum / sweeps_done : 1.0;
    }
};

class DqmcEngine {
public:
    DqmcEngine(const Lattice& lat,
               const HubbardModel& model,
               const DqmcParams& params,
               Pcg32 rng);

    // Run one full sweep (every site at every time slice gets a
    // single-spin-flip proposal). Returns the configuration sign for
    // this sweep (+1 except in the rare case of a sign change).
    Real sweep();

    const Lattice&      lattice() const { return lat_; }
    const HubbardModel& model()   const { return model_; }
    const DqmcParams&   params()  const { return params_; }
    const DqmcStats&    stats()   const { return stats_; }
    Real                dt()      const { return dt_; }
    Real                alpha()   const { return alpha_; }

    // Equal-time Green's function for spin sigma at the *current*
    // time slice (l_now). Always recomputed from the latest stable
    // factorisation, so this is safe to call any time after sweep().
    const la::Matrix& green_up()   const { return G_up_;   }
    const la::Matrix& green_down() const { return G_down_; }

    // Auxiliary field accessor: hs_field_(l, i) in {-1, +1}.
    int hs(int l, int i) const { return hs_[l * lat_.n_sites() + i]; }
    int& hs(int l, int i)       { return hs_[l * lat_.n_sites() + i]; }

private:
    // Build  exp(V_l_sigma)  as a *diagonal* (returned as a vector).
    std::vector<Real> exp_V_diag_(int l, int sigma) const;

    // Multiply  M  on the left by  B_l_sigma  in-place: M <- B_l M.
    void mul_B_left_ (int l, int sigma, la::Matrix& M) const;
    // Multiply  M  on the right by  B_l_sigma^{-1}: M <- M * B_l^{-1}.
    void mul_Binv_right_(int l, int sigma, la::Matrix& M) const;

    // Recompute G_sigma at slice l (forward sweep) from a stable
    // QR-based product of B-matrices. Returns max element-wise drift
    // versus the incoming `G` (used for monitoring).
    Real recompute_green_(int l, int sigma, la::Matrix& G) const;

    void wrap_green_up_(int l);   // G_up   <-  B_l_up   * G_up   * B_l_up^{-1}
    void wrap_green_dn_(int l);

    // --- members ---
    const Lattice&     lat_;
    HubbardModel       model_;
    DqmcParams         params_;
    Pcg32              rng_;
    Real               dt_{0.0};
    Real               alpha_{0.0};                // HS coupling

    la::Matrix         expK_;       // exp(-dt * K)            (Ns x Ns)
    la::Matrix         expmK_;      // exp(+dt * K)            (Ns x Ns)
    la::Matrix         G_up_;       // current equal-time G_up
    la::Matrix         G_down_;     // current equal-time G_dn

    std::vector<int>   hs_;         // L_tau * Ns auxiliary spins (+/- 1)

    DqmcStats          stats_;
};

// ---------------------------------------------------------------------------
// Implementation
// ---------------------------------------------------------------------------

inline DqmcEngine::DqmcEngine(const Lattice& lat,
                              const HubbardModel& model,
                              const DqmcParams& params,
                              Pcg32 rng)
    : lat_(lat),
      model_(model),
      params_(params),
      rng_(rng) {
    model_.validate();
    if (params_.beta <= 0.0)  throw std::invalid_argument("DQMC: beta must be > 0");
    if (params_.L_tau <= 0)   throw std::invalid_argument("DQMC: L_tau must be > 0");
    if (params_.n_stab <= 0)  params_.n_stab = 1;
    dt_    = params_.beta / params_.L_tau;
    // Hirsch HS:  exp(-dt U (n_up-1/2)(n_dn-1/2))
    //   = (e^{-dt U/4}/2) sum_{s=+/-1} exp(alpha s (n_up - n_dn))
    //   with cosh(alpha) = exp(dt U / 2).
    if (model_.U > 0.0) {
        alpha_ = std::acosh(std::exp(0.5 * dt_ * model_.U));
    } else {
        alpha_ = 0.0;     // U=0: aux field decouples entirely
    }

    const int Ns = lat_.n_sites();
    auto K  = model_.kinetic_matrix(lat_);
    expK_   = la::expm_sym(K, +dt_);   // exp(-dt * K)
    expmK_  = la::expm_sym(K, -dt_);   // exp(+dt * K)

    // Random initial Hubbard-Stratonovich field.
    hs_.assign(static_cast<std::size_t>(params_.L_tau) * Ns, +1);
    for (auto& s : hs_) s = rng_.bernoulli(0.5) ? +1 : -1;

    // Initialize Green's functions at slice l = 0 via the stable
    // (chunked-QR) recompute path used during sweeps.
    G_up_   = la::Matrix(Ns, Ns, 0.0);
    G_down_ = la::Matrix(Ns, Ns, 0.0);
    recompute_green_(0, +1, G_up_);
    recompute_green_(0, -1, G_down_);
}

inline std::vector<Real> DqmcEngine::exp_V_diag_(int l, int sigma) const {
    const int Ns = lat_.n_sites();
    std::vector<Real> d(Ns);
    const Real a = (sigma == +1) ? alpha_ : -alpha_;
    for (int i = 0; i < Ns; ++i) {
        d[i] = std::exp(a * hs_[l * Ns + i]);
    }
    return d;
}

inline void DqmcEngine::mul_B_left_(int l, int sigma, la::Matrix& M) const {
    // B_l = exp(-dt K) * exp(V_l)
    // M <- B_l * M  =  expK * (diag(d) * M)
    auto d = exp_V_diag_(l, sigma);
    la::diag_mul_left(d, M);
    M = la::matmul(expK_, M);
}

inline void DqmcEngine::mul_Binv_right_(int l, int sigma, la::Matrix& M) const {
    // B_l^{-1} = exp(-V_l) * exp(+dt K)
    // M <- M * B_l^{-1}  =  (M * expmK) * diag(1/d)
    auto d = exp_V_diag_(l, sigma);
    M = la::matmul(M, expmK_);
    std::vector<Real> dinv(d.size());
    for (std::size_t i = 0; i < d.size(); ++i) dinv[i] = 1.0 / d[i];
    la::diag_mul_right(M, dinv);
}

inline Real DqmcEngine::recompute_green_(int l, int sigma, la::Matrix& G) const {
    // Build  A_sigma(l) = B_{l-1} B_{l-2} ... B_0 B_{Ltau-1} ... B_{l+1} B_l
    // by left-multiplications, with *incremental* QR stabilization every
    // n_stab steps.  This is the standard Loh-Gubernatis-style scheme:
    // we maintain (Q, R) such that the partial product equals Q * R at
    // all times.  Each chunk multiplies a few B-matrices into Q, then
    // re-factorizes  (B_chunk * Q) = Q_new * R_chunk_local  and updates
    // R := R_chunk_local * R.  Q remains orthogonal throughout, which
    // tames the catastrophic conditioning of  exp(-dt K)^{Lt}  for
    // moderate to low temperatures.  The final  G = (I + Q R)^{-1}
    // is then obtained from the stable formula
    //
    //     G = (Q^T + R)^{-1} Q^T
    //
    // which is conditioned by the actual G itself rather than by A.
    const int Ns    = lat_.n_sites();
    const int Lt    = params_.L_tau;
    const int chunk = std::max(1, params_.n_stab);

    la::Matrix Q = la::Matrix::identity(Ns);
    la::Matrix R = la::Matrix::identity(Ns);
    la::Matrix C = la::Matrix::identity(Ns);   // accumulator for the current chunk
    int chunk_count = 0;

    auto flush_chunk = [&]() {
        // Refactorize  (C * Q) = Q_new * R_local;  R := R_local * R.
        la::Matrix M  = la::matmul(C, Q);
        auto qr = la::qr_factor(std::move(M));
        Q = std::move(qr.Q);
        R = la::matmul(qr.R, R);
        C = la::Matrix::identity(Ns);
        chunk_count = 0;
    };

    for (int step = 0; step < Lt; ++step) {
        const int idx = (l + step) % Lt;
        mul_B_left_(idx, sigma, C);
        if (++chunk_count == chunk) flush_chunk();
    }
    if (chunk_count > 0) flush_chunk();

    // Stable Green: G = (Q^T + R)^{-1} Q^T.
    la::Matrix M(Ns, Ns, 0.0);
    for (int i = 0; i < Ns; ++i)
        for (int j = 0; j < Ns; ++j) M(i, j) = R(i, j) + Q(j, i);
    la::Matrix QT(Ns, Ns, 0.0);
    for (int i = 0; i < Ns; ++i)
        for (int j = 0; j < Ns; ++j) QT(i, j) = Q(j, i);
    la::Matrix Gnew = la::lu_solve(la::lu_factor(std::move(M)), std::move(QT));
    // Compute drift versus the in-place wrapped G (for monitoring).
    Real drift = 0.0;
    for (int i = 0; i < Ns; ++i)
        for (int j = 0; j < Ns; ++j)
            drift = std::max(drift, std::fabs(Gnew(i, j) - G(i, j)));
    G = std::move(Gnew);
    return drift;
}

inline void DqmcEngine::wrap_green_up_(int l) {
    // After advancing past slice l (forward),  G(l+1) = B_l G(l) B_l^{-1}
    // with  B_l = expK * diag(d)  and  B_l^{-1} = diag(1/d) * expmK.
    // The correct factor ordering is therefore
    //
    //   G' = expK * (diag(d) * G * diag(1/d)) * expmK,
    //
    // NOT  expK * diag(d) * G * expmK * diag(1/d) -- diag(1/d) does
    // not commute with expmK.
    auto d = exp_V_diag_(l, +1);
    std::vector<Real> dinv(d.size());
    for (std::size_t i = 0; i < d.size(); ++i) dinv[i] = 1.0 / d[i];
    la::diag_mul_left (d,    G_up_);
    la::diag_mul_right(G_up_, dinv);
    G_up_ = la::matmul(expK_,  G_up_);
    G_up_ = la::matmul(G_up_,  expmK_);
}
inline void DqmcEngine::wrap_green_dn_(int l) {
    auto d = exp_V_diag_(l, -1);
    std::vector<Real> dinv(d.size());
    for (std::size_t i = 0; i < d.size(); ++i) dinv[i] = 1.0 / d[i];
    la::diag_mul_left (d,      G_down_);
    la::diag_mul_right(G_down_, dinv);
    G_down_ = la::matmul(expK_,    G_down_);
    G_down_ = la::matmul(G_down_,  expmK_);
}

inline Real DqmcEngine::sweep() {
    const int Ns = lat_.n_sites();
    const int Lt = params_.L_tau;
    Real cfg_sign = +1.0;

    for (int l = 0; l < Lt; ++l) {
        // Single-spin-flip local updates at every site of slice l.
        for (int i = 0; i < Ns; ++i) {
            ++stats_.flips_proposed;
            const int s_old = hs_[l * Ns + i];
            // Delta_sigma = exp(-2 alpha s_old * sigma) - 1
            //  sigma = +1: Delta_up = exp(-2 alpha s_old) - 1
            //  sigma = -1: Delta_dn = exp(+2 alpha s_old) - 1
            const Real e_minus = std::exp(-2.0 * alpha_ * s_old);
            const Real e_plus  = std::exp(+2.0 * alpha_ * s_old);
            const Real Du = e_minus - 1.0;
            const Real Dd = e_plus  - 1.0;
            const Real Ru = 1.0 + Du * (1.0 - G_up_  (i, i));
            const Real Rd = 1.0 + Dd * (1.0 - G_down_(i, i));
            const Real R  = Ru * Rd;
            if (R == 0.0) continue;
            const Real p_acc = std::min(1.0, std::fabs(R));
            if (rng_.uniform() < p_acc) {
                ++stats_.flips_accepted;
                if (R < 0.0) cfg_sign = -cfg_sign;
                // Sherman-Morrison rank-1 update of the equal-time Green's
                // function for a single-spin flip at site i.  In convention
                // (a) (G_sigma(l) = [I + A(l)]^{-1} with A(l) = B_{l-1} ...
                // B_{l+1} B_l), an aux-spin flip changes A(l) by a rank-1
                // matrix on the *right*:  A' = A (I + Delta E_{ii}).  The
                // resulting update is
                //
                //     G'_{jk} = G_{jk} - (Delta / R) * (delta_{ji} - G_{ji}) * G_{ik},
                //
                // i.e. an outer product of column i of (I - G) with row i of G.
                auto sherman = [&](la::Matrix& G, Real Delta, Real Rsig) {
                    if (Delta == 0.0) return;
                    const Real f = Delta / Rsig;
                    std::vector<Real> u(Ns), row(Ns);
                    for (int j = 0; j < Ns; ++j) u[j]   = ((j == i) ? 1.0 : 0.0) - G(j, i);
                    for (int k = 0; k < Ns; ++k) row[k] = G(i, k);
                    for (int j = 0; j < Ns; ++j) {
                        const Real fj = f * u[j];
                        if (fj == 0.0) continue;
                        Real* Grow = G.data() + static_cast<std::size_t>(j) * Ns;
                        for (int k = 0; k < Ns; ++k) Grow[k] -= fj * row[k];
                    }
                };
                sherman(G_up_,   Du, Ru);
                sherman(G_down_, Dd, Rd);
                hs_[l * Ns + i] = -s_old;
            }
        }

        // Advance to the next time slice.
        const int next = (l + 1) % Lt;
        wrap_green_up_(l);
        wrap_green_dn_(l);
        // Periodic stabilization: recompute G fresh from a stable
        // QR-based product. Use the maximum element-wise drift versus
        // the wrapped G as a numerical-quality monitor.
        if ((next % params_.n_stab) == 0) {
            const Real du = recompute_green_(next, +1, G_up_);
            const Real dd = recompute_green_(next, -1, G_down_);
            stats_.max_recompute_drift =
                std::max(stats_.max_recompute_drift, std::max(du, dd));
        }
    }

    ++stats_.sweeps_done;
    stats_.sign_running_sum += cfg_sign;
    return cfg_sign;
}

} // namespace qmc
