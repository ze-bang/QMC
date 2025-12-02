/**
 * @file measurements.cpp
 * @brief Implementation of measurement routines
 */

#include "sse/measurements.hpp"
#include <fmt/format.h>
#include <iostream>
#include <cmath>

namespace sse {

Measurements::Measurements(const SSEConfig& config, Real beta, int n_bins)
    : beta_(beta),
      n_sites_(config.getLattice().numSites()),
      n_operators_("n_operators", n_bins),
      energy_("energy", n_bins),
      specific_heat_("specific_heat", n_bins),
      magnetization_("magnetization", n_bins),
      magnetization_sq_("magnetization_sq", n_bins),
      magnetization_4th_("magnetization_4th", n_bins),
      susceptibility_("susceptibility", n_bins),
      stag_magnetization_("stag_magnetization", n_bins),
      binder_ratio_("binder_ratio", n_bins) {
    reset();
}

void Measurements::reset() {
    accum_n_ = 0.0;
    accum_n_sq_ = 0.0;
    accum_mag_ = 0.0;
    accum_mag_sq_ = 0.0;
    accum_stag_mag_ = 0.0;
    n_samples_ = 0;
}

void Measurements::startBin() {
    n_operators_.startBin();
    energy_.startBin();
    specific_heat_.startBin();
    magnetization_.startBin();
    magnetization_sq_.startBin();
    magnetization_4th_.startBin();
    susceptibility_.startBin();
    stag_magnetization_.startBin();
    binder_ratio_.startBin();
}

void Measurements::endBin() {
    n_operators_.endBin();
    energy_.endBin();
    specific_heat_.endBin();
    magnetization_.endBin();
    magnetization_sq_.endBin();
    magnetization_4th_.endBin();
    susceptibility_.endBin();
    stag_magnetization_.endBin();
    binder_ratio_.endBin();
}

void Measurements::measure(const SSEConfig& config) {
    measureEnergy(config);
    measureMagnetization(config);
}

void Measurements::measureEnergy(const SSEConfig& config) {
    // Energy estimator in SSE:
    // ⟨H⟩ = -⟨n⟩/β + const
    // where n is the number of operators
    
    Real n = static_cast<Real>(config.numOperators());
    n_operators_.add(n);
    
    // Energy per site (including offset)
    BondIdx n_bonds = config.getLattice().numBonds();
    const VertexDataCollection& vd = config.getVertexData();
    Real offset = 0.0;
    for (int t = 0; t < vd.numBondTypes(); ++t) {
        offset += vd.get(t).getEnergyOffset();
    }
    offset *= static_cast<Real>(n_bonds) / vd.numBondTypes();
    
    Real E = -n / beta_ + offset;
    Real E_per_site = E / n_sites_;
    energy_.add(E_per_site);
    
    // Specific heat from fluctuations: C = β²(⟨n²⟩ - ⟨n⟩²)/N
    accum_n_ += n;
    accum_n_sq_ += n * n;
    ++n_samples_;
    
    if (n_samples_ > 1) {
        Real avg_n = accum_n_ / n_samples_;
        Real avg_n_sq = accum_n_sq_ / n_samples_;
        Real var_n = avg_n_sq - avg_n * avg_n;
        Real C = beta_ * beta_ * var_n / n_sites_;
        specific_heat_.add(C);
    }
}

void Measurements::measureMagnetization(const SSEConfig& config) {
    const auto& spins = config.getSpins();
    const Lattice& lat = config.getLattice();
    
    // Uniform magnetization: m = (1/N) Σ Sz_i = (1/N) Σ (s_i - 1/2)
    Real m = 0.0;
    for (SiteIdx s = 0; s < n_sites_; ++s) {
        m += spins[s] - 0.5;  // Convert 0/1 to -1/2, +1/2
    }
    m /= n_sites_;
    
    magnetization_.add(m);
    magnetization_sq_.add(m * m);
    magnetization_4th_.add(m * m * m * m);
    
    // Susceptibility: χ = β N ⟨m²⟩
    Real chi = beta_ * n_sites_ * m * m;
    susceptibility_.add(chi);
    
    // Staggered magnetization (for bipartite lattices)
    // m_s = (1/N) Σ (-1)^i Sz_i
    Real m_stag = 0.0;
    if (lat.type() == LatticeType::Square || lat.type() == LatticeType::Chain) {
        auto dims = lat.dimensions();
        if (!dims.empty()) {
            if (lat.type() == LatticeType::Chain) {
                for (SiteIdx s = 0; s < n_sites_; ++s) {
                    int sign = (s % 2 == 0) ? 1 : -1;
                    m_stag += sign * (spins[s] - 0.5);
                }
            } else {
                SiteIdx Lx = dims[0];
                for (SiteIdx s = 0; s < n_sites_; ++s) {
                    int x = s % Lx;
                    int y = s / Lx;
                    int sign = ((x + y) % 2 == 0) ? 1 : -1;
                    m_stag += sign * (spins[s] - 0.5);
                }
            }
        }
        m_stag /= n_sites_;
    }
    stag_magnetization_.add(m_stag * m_stag);
    
    // Binder ratio: U = 1 - ⟨m⁴⟩/(3⟨m²⟩²)
    accum_mag_ += m * m;
    accum_mag_sq_ += m * m * m * m;
    
    if (n_samples_ > 1 && accum_mag_ > 0) {
        Real avg_m2 = accum_mag_ / n_samples_;
        Real avg_m4 = accum_mag_sq_ / n_samples_;
        if (avg_m2 > 1e-10) {
            Real U = 1.0 - avg_m4 / (3.0 * avg_m2 * avg_m2);
            binder_ratio_.add(U);
        }
    }
}

void Measurements::print() const {
    std::cout << "\n=== Measurement Results ===\n";
    std::cout << fmt::format("Energy/site:        {:.6f} ± {:.6f}\n",
                             energy_.mean(), energy_.stdError());
    std::cout << fmt::format("Specific heat/site: {:.6f} ± {:.6f}\n",
                             specific_heat_.mean(), specific_heat_.stdError());
    std::cout << fmt::format("⟨m²⟩:               {:.6f} ± {:.6f}\n",
                             magnetization_sq_.mean(), magnetization_sq_.stdError());
    std::cout << fmt::format("Susceptibility:     {:.6f} ± {:.6f}\n",
                             susceptibility_.mean(), susceptibility_.stdError());
    std::cout << fmt::format("⟨m_s²⟩:             {:.6f} ± {:.6f}\n",
                             stag_magnetization_.mean(), stag_magnetization_.stdError());
    std::cout << fmt::format("Binder ratio:       {:.6f} ± {:.6f}\n",
                             binder_ratio_.mean(), binder_ratio_.stdError());
    std::cout << fmt::format("⟨n⟩:                {:.1f}\n", n_operators_.mean());
}

std::unordered_map<std::string, std::pair<Real, Real>> Measurements::getResults() const {
    return {
        {"energy", {energy_.mean(), energy_.stdError()}},
        {"specific_heat", {specific_heat_.mean(), specific_heat_.stdError()}},
        {"magnetization_sq", {magnetization_sq_.mean(), magnetization_sq_.stdError()}},
        {"susceptibility", {susceptibility_.mean(), susceptibility_.stdError()}},
        {"stag_magnetization", {stag_magnetization_.mean(), stag_magnetization_.stdError()}},
        {"binder_ratio", {binder_ratio_.mean(), binder_ratio_.stdError()}},
        {"n_operators", {n_operators_.mean(), n_operators_.stdError()}}
    };
}

// CorrelationMeasurements implementation

CorrelationMeasurements::CorrelationMeasurements(const SSEConfig& config, int n_bins)
    : n_sites_(config.getLattice().numSites()) {
    correlations_.reserve(n_sites_);
    for (SiteIdx r = 0; r < n_sites_; ++r) {
        correlations_.emplace_back(fmt::format("corr_{}", r), n_bins);
    }
}

void CorrelationMeasurements::measure(const SSEConfig& config) {
    const auto& spins = config.getSpins();
    
    // Reference site (site 0)
    Real s0 = spins[0] - 0.5;
    
    // Compute ⟨S_0 · S_r⟩
    for (SiteIdx r = 0; r < n_sites_; ++r) {
        Real sr = spins[r] - 0.5;
        correlations_[r].add(s0 * sr);
    }
}

Real CorrelationMeasurements::correlation(SiteIdx r) const {
    return correlations_[r].mean();
}

Real CorrelationMeasurements::correlationError(SiteIdx r) const {
    return correlations_[r].stdError();
}

Real CorrelationMeasurements::structureFactor(const std::vector<Real>& k) const {
    // S(k) = (1/N) Σ_r exp(i k·r) ⟨S_0·S_r⟩
    // For real correlations, this simplifies to Σ cos(k·r) ⟨S_0·S_r⟩
    
    Real S_real = 0.0;
    // This would need lattice geometry to compute properly
    // Simplified version for 1D:
    for (SiteIdx r = 0; r < n_sites_; ++r) {
        Real kr = k[0] * r;  // Assuming 1D
        S_real += std::cos(kr) * correlations_[r].mean();
    }
    return S_real / n_sites_;
}

void CorrelationMeasurements::reset() {
    for (auto& corr : correlations_) {
        // BinnedObservable doesn't have reset, would need to add
    }
}

} // namespace sse
