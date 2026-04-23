// SPDX-License-Identifier: MIT
//
// Standard observables for the spin-1/2 AF Heisenberg model evaluated
// from an SSE configuration.
//
// All formulas below assume:
//   * the operator string was generated for the bond Hamiltonian with
//     constant offset C = 1/4 per bond,
//   * spins are stored as ±1 (so S^z = s/2 with s ∈ {−1,+1}),
//   * the lattice is bipartite with sublattice index ∈ {0,1}.
//
// References:
//   A. W. Sandvik, "Computational studies of quantum spin systems",
//   AIP Conf. Proc. 1297, 135 (2010), Sec. 5.

#pragma once

#include <cmath>

#include "qmc/observable.hpp"
#include "qmc/sse_engine.hpp"

namespace qmc {

class Measurements {
public:
    explicit Measurements(const SseEngine& engine)
        : engine_(engine),
          energy_("energy_per_site"),
          energy_sq_("energy_per_site_sq"),
          n_op_("<n>"),
          n_op_sq_("<n^2>"),
          mz_abs_("|m_z|"),
          mz_sq_("m_z^2"),
          ms_abs_("|m_stag|"),
          ms_sq_("m_stag^2"),
          ms4_("m_stag^4"),
          chi_uniform_("chi_uniform"),
          chi_stag_("chi_stag") {}

    void clear() {
        energy_.clear();
        energy_sq_.clear();
        n_op_.clear();
        n_op_sq_.clear();
        mz_abs_.clear();
        mz_sq_.clear();
        ms_abs_.clear();
        ms_sq_.clear();
        ms4_.clear();
        chi_uniform_.clear();
        chi_stag_.clear();
        n_jack_.clear();
        n_sq_jack_.clear();
    }

    void measure() {
        const auto& lat   = engine_.lattice();
        const auto& spins = engine_.spins();
        const auto& sub   = lat.sublattice();
        const Real  beta  = engine_.beta();
        const Real  J     = engine_.model().J;
        const Real  E0    = engine_.model().energy_offset(lat);
        const int   Ns    = lat.n_sites();
        const Real  invN  = 1.0 / Ns;

        // Energy: E = -<n>/beta + E_offset, then per site.
        const Length n  = engine_.n_op();
        const Real   En = -static_cast<Real>(n) / beta + E0;
        const Real   E_per_site = En * invN;
        energy_.add(E_per_site);
        energy_sq_.add(E_per_site * E_per_site);
        n_op_.add(static_cast<Real>(n));
        n_op_sq_.add(static_cast<Real>(n) * static_cast<Real>(n));

        // Magnetizations (uniform and staggered) along z.
        // m_z = (1/N) sum_i s_i / 2
        // m_s = (1/N) sum_i (-1)^{sub(i)} s_i / 2
        Real mz = 0.0, ms = 0.0;
        for (Site i = 0; i < Ns; ++i) {
            const Real si = 0.5 * spins[i];
            mz += si;
            ms += sub[i] ? -si : si;
        }
        mz *= invN;
        ms *= invN;
        const Real mz2 = mz * mz;
        const Real ms2 = ms * ms;
        mz_abs_.add(std::abs(mz));
        mz_sq_.add(mz2);
        ms_abs_.add(std::abs(ms));
        ms_sq_.add(ms2);
        ms4_.add(ms2 * ms2);

        // Static susceptibilities. The integrated form is
        //     chi = beta * <m^2>  (for ground-state-like saturated runs:
        //                          actually beta * <(M)^2> where M is the
        //                          *total* magnetization summed once).
        // We use Sandvik's "fast" form: chi_q = (beta / N) * <(M_q)^2>
        // for a static spin component. For the operator-string sampling
        // of <m^2>, the standard estimator is the average of the
        // squared total magnetization at a single time slice (since
        // SSE preserves total S^z, this is exact for the uniform
        // susceptibility). For the staggered susceptibility the simple
        // single-slice estimator is used here as a baseline; a more
        // accurate integrated form using the operator string is left
        // for a future revision.
        const Real Mz = mz * Ns;     // total uniform magnetization
        const Real Ms = ms * Ns;     // total staggered magnetization
        chi_uniform_.add(beta * Mz * Mz * invN);
        chi_stag_.add(   beta * Ms * Ms * invN);

        // Jackknife inputs for derived quantities (specific heat).
        n_jack_.add(static_cast<Real>(n));
        n_sq_jack_.add(static_cast<Real>(n) * static_cast<Real>(n));

        (void)J;
    }

    // Specific heat per site, computed via the SSE estimator
    //     C = (<n^2> - <n>^2 - <n>) / N
    // (Sandvik, AIP 2010, eq. (113)). Returned with a jackknife error.
    std::pair<Real, Real> specific_heat() const {
        std::vector<JackknifeSamples*> sources{
            const_cast<JackknifeSamples*>(&n_jack_),
            const_cast<JackknifeSamples*>(&n_sq_jack_)};
        const Real Ns = engine_.lattice().n_sites();
        return JackknifeSamples::jackknife(sources, [Ns](const std::vector<Real>& m) {
            const Real n  = m[0];
            const Real n2 = m[1];
            return (n2 - n * n - n) / Ns;
        });
    }

    // Binder cumulant of the staggered magnetization,
    //     U = 1 - <m_s^4> / (3 <m_s^2>^2)
    // Useful for finite-size scaling at the AF transition (irrelevant
    // for d=1 but reported nonetheless).
    Real binder_cumulant() const {
        const Real m2 = ms_sq_.mean();
        const Real m4 = ms4_.mean();
        if (m2 <= 0.0) return 0.0;
        return 1.0 - m4 / (3.0 * m2 * m2);
    }

    // Direct accessors -------------------------------------------------
    const Observable& energy()      const { return energy_; }
    const Observable& n_op()        const { return n_op_; }
    const Observable& mz_abs()      const { return mz_abs_; }
    const Observable& mz_sq()       const { return mz_sq_; }
    const Observable& ms_abs()      const { return ms_abs_; }
    const Observable& ms_sq()       const { return ms_sq_; }
    const Observable& ms4()         const { return ms4_; }
    const Observable& chi_uniform() const { return chi_uniform_; }
    const Observable& chi_stag()    const { return chi_stag_; }

    void report(std::ostream& os) const {
        auto line = [&](const Observable& o) { o.format(os); os << '\n'; };
        line(energy_);
        line(n_op_);
        line(mz_abs_);
        line(mz_sq_);
        line(ms_abs_);
        line(ms_sq_);
        line(chi_uniform_);
        line(chi_stag_);
        const auto [cv, dcv] = specific_heat();
        os << std::setw(28) << std::left << "specific_heat (jack)"
           << std::setw(16) << std::scientific << std::setprecision(8) << cv
           << " +/- " << std::setw(12) << std::scientific << std::setprecision(3) << dcv << '\n';
        os << std::setw(28) << std::left << "binder_cumulant"
           << std::setw(16) << std::scientific << std::setprecision(8) << binder_cumulant()
           << '\n';
    }

private:
    const SseEngine&  engine_;
    Observable        energy_;
    Observable        energy_sq_;
    Observable        n_op_;
    Observable        n_op_sq_;
    Observable        mz_abs_;
    Observable        mz_sq_;
    Observable        ms_abs_;
    Observable        ms_sq_;
    Observable        ms4_;
    Observable        chi_uniform_;
    Observable        chi_stag_;

    JackknifeSamples  n_jack_;
    JackknifeSamples  n_sq_jack_;
};

} // namespace qmc
