// SPDX-License-Identifier: MIT
//
// Equal-time observables for the DQMC engine on the Hubbard model.
//
// Conventions (matching `DqmcEngine`):
//   * G_sigma(i, j) = <c_i c_j^+>     (single-particle Green's function)
//   * <n_i sigma>   = 1 - G_sigma(i, i)
//   * <n_i up n_i down> = (1 - G_up(i,i)) * (1 - G_dn(i,i))    (Wick)
//
// Two-point estimators use Wick's theorem on the Hubbard-Stratonovich
// configuration: for a fixed s the system is non-interacting and
// Wick contractions are exact.

#pragma once

#include <cmath>

#include "qmc/dqmc_engine.hpp"
#include "qmc/observable.hpp"

namespace qmc {

class DqmcMeasurements {
public:
    explicit DqmcMeasurements(const DqmcEngine& eng)
        : eng_(eng),
          density_("density"),
          double_occ_("double_occupancy"),
          mz_sq_("<m_z^2>"),
          ms_sq_("<m_stag^2>"),
          struct_factor_pi_("S(pi,pi)"),
          sign_("<sign>") {}

    void clear() {
        density_.clear();
        double_occ_.clear();
        mz_sq_.clear();
        ms_sq_.clear();
        struct_factor_pi_.clear();
        sign_.clear();
    }

    void measure(Real config_sign = 1.0) {
        const auto& Gu = eng_.green_up();
        const auto& Gd = eng_.green_down();
        const auto& lat = eng_.lattice();
        const auto& sub = lat.sublattice();
        const int Ns = lat.n_sites();

        // <n_sigma>_i = 1 - G_sigma(i, i).
        Real n_total = 0.0;
        Real docc    = 0.0;
        Real mz2     = 0.0;
        Real ms2     = 0.0;
        for (int i = 0; i < Ns; ++i) {
            const Real nu = 1.0 - Gu(i, i);
            const Real nd = 1.0 - Gd(i, i);
            n_total += nu + nd;
            docc    += nu * nd;
        }
        // <m_z^2> = (1/Ns^2) sum_{ij} <S^z_i S^z_j>
        // and the staggered version  <m_s^2> with sign  (-1)^{sub_i + sub_j}.
        // S^z_i = (n_up - n_down)/2.
        // Using Wick:
        //   <n_iu n_ju> = nu_i nu_j + delta_ij nu_i (1 - nu_i)
        //                 - G_u(i,j) * G_u(j,i)        (off-diagonal i!=j)
        // For DQMC, the cleaner standard form:
        //   <S^z_i S^z_j> = (1/4)
        //       * { (1 - G_u(i,i))(1 - G_u(j,j)) + (1 - G_d(i,i))(1 - G_d(j,j))
        //           - 2 (1 - G_u(i,i))(1 - G_d(j,j))
        //           + (delta_ij - G_u(j,i)) G_u(i,j)
        //           + (delta_ij - G_d(j,i)) G_d(i,j) }.
        for (int i = 0; i < Ns; ++i) {
            for (int j = 0; j < Ns; ++j) {
                const Real nui = 1.0 - Gu(i, i);
                const Real nuj = 1.0 - Gu(j, j);
                const Real ndi = 1.0 - Gd(i, i);
                const Real ndj = 1.0 - Gd(j, j);
                const Real delta = (i == j) ? 1.0 : 0.0;
                const Real exch_u = (delta - Gu(j, i)) * Gu(i, j);
                const Real exch_d = (delta - Gd(j, i)) * Gd(i, j);
                const Real SzSz = 0.25 * (
                    nui * nuj + ndi * ndj - 2.0 * nui * ndj + exch_u + exch_d);
                mz2 += SzSz;
                const Real sign_ij = ((sub[i] ^ sub[j]) ? -1.0 : 1.0);
                ms2 += sign_ij * SzSz;
            }
        }
        const Real invN  = 1.0 / Ns;
        const Real invN2 = invN * invN;

        density_         .add(config_sign * n_total * invN);
        double_occ_      .add(config_sign * docc    * invN);
        mz_sq_           .add(config_sign * mz2     * invN2);
        ms_sq_           .add(config_sign * ms2     * invN2);
        struct_factor_pi_.add(config_sign * ms2     * invN);
        sign_            .add(config_sign);
    }

    const Observable& density()        const { return density_; }
    const Observable& double_occ()     const { return double_occ_; }
    const Observable& mz_sq()          const { return mz_sq_; }
    const Observable& ms_sq()          const { return ms_sq_; }
    const Observable& structure_factor_pi() const { return struct_factor_pi_; }
    const Observable& sign()           const { return sign_; }

    void report(std::ostream& os) const {
        auto line = [&](const Observable& o) { o.format(os); os << '\n'; };
        line(sign_);
        line(density_);
        line(double_occ_);
        line(mz_sq_);
        line(ms_sq_);
        line(struct_factor_pi_);
    }

private:
    const DqmcEngine& eng_;
    Observable        density_;
    Observable        double_occ_;
    Observable        mz_sq_;
    Observable        ms_sq_;
    Observable        struct_factor_pi_;
    Observable        sign_;
};

} // namespace qmc
