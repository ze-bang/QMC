#pragma once

/**
 * @file measurements.hpp
 * @brief Physical observables and measurements for SSE QMC
 * 
 * Implements estimators for various physical quantities:
 * - Energy
 * - Specific heat
 * - Magnetization
 * - Susceptibility
 * - Correlation functions
 * - Structure factors
 */

#include <vector>
#include <array>
#include <string>
#include <unordered_map>
#include <cmath>
#include "sse/types.hpp"
#include "sse/sse_config.hpp"

namespace sse {

/**
 * @brief Statistics accumulator for single observable
 */
class Observable {
public:
    Observable(const std::string& name = "") : name_(name) {}
    
    void add(Real value) {
        sum_ += value;
        sum_sq_ += value * value;
        ++count_;
    }
    
    void reset() {
        sum_ = 0.0;
        sum_sq_ = 0.0;
        count_ = 0;
    }
    
    Real mean() const {
        return count_ > 0 ? sum_ / count_ : 0.0;
    }
    
    Real variance() const {
        if (count_ < 2) return 0.0;
        Real m = mean();
        return (sum_sq_ / count_ - m * m) * count_ / (count_ - 1);
    }
    
    Real stdError() const {
        return count_ > 0 ? std::sqrt(variance() / count_) : 0.0;
    }
    
    int64_t count() const { return count_; }
    const std::string& name() const { return name_; }
    Real sum() const { return sum_; }
    Real sumSq() const { return sum_sq_; }
    
private:
    std::string name_;
    Real sum_ = 0.0;
    Real sum_sq_ = 0.0;
    int64_t count_ = 0;
};

/**
 * @brief Binned statistics for error estimation
 */
class BinnedObservable {
public:
    BinnedObservable(const std::string& name = "", int n_bins = 100) 
        : name_(name), n_bins_(n_bins), bins_(n_bins) {}
    
    void startBin() {
        current_bin_.reset();
    }
    
    void add(Real value) {
        current_bin_.add(value);
    }
    
    void endBin() {
        if (current_bin_.count() > 0) {
            bins_[current_bin_idx_ % n_bins_] = current_bin_.mean();
            ++current_bin_idx_;
        }
    }
    
    Real mean() const {
        int n = std::min(current_bin_idx_, n_bins_);
        if (n == 0) return 0.0;
        Real sum = 0.0;
        for (int i = 0; i < n; ++i) sum += bins_[i];
        return sum / n;
    }
    
    Real stdError() const {
        int n = std::min(current_bin_idx_, n_bins_);
        if (n < 2) return 0.0;
        Real m = mean();
        Real var = 0.0;
        for (int i = 0; i < n; ++i) {
            Real d = bins_[i] - m;
            var += d * d;
        }
        var /= (n - 1);
        return std::sqrt(var / n);
    }
    
    const std::string& name() const { return name_; }
    
private:
    std::string name_;
    int n_bins_;
    std::vector<Real> bins_;
    Observable current_bin_;
    int current_bin_idx_ = 0;
};

/**
 * @brief Measurements class for SSE simulation
 */
class Measurements {
public:
    /**
     * @brief Construct measurements for given configuration
     */
    Measurements(const SSEConfig& config, Real beta, int n_bins = 100);
    
    /**
     * @brief Perform all measurements on current configuration
     */
    void measure(const SSEConfig& config);
    
    /**
     * @brief Start new bin
     */
    void startBin();
    
    /**
     * @brief End current bin
     */
    void endBin();
    
    /**
     * @brief Reset all measurements
     */
    void reset();
    
    // Accessors for results
    Real energy() const { return energy_.mean(); }
    Real energyError() const { return energy_.stdError(); }
    
    Real specificHeat() const { return specific_heat_.mean(); }
    Real specificHeatError() const { return specific_heat_.stdError(); }
    
    Real magnetization() const { return magnetization_.mean(); }
    Real magnetizationError() const { return magnetization_.stdError(); }
    
    Real magnetizationSquared() const { return magnetization_sq_.mean(); }
    Real magnetizationSquaredError() const { return magnetization_sq_.stdError(); }
    
    Real susceptibility() const { return susceptibility_.mean(); }
    Real susceptibilityError() const { return susceptibility_.stdError(); }
    
    Real stagMagnetization() const { return stag_magnetization_.mean(); }
    Real stagMagnetizationError() const { return stag_magnetization_.stdError(); }
    
    Real binderRatio() const { return binder_ratio_.mean(); }
    Real binderRatioError() const { return binder_ratio_.stdError(); }
    
    int64_t avgOperators() const { 
        return static_cast<int64_t>(n_operators_.mean()); 
    }
    
    /**
     * @brief Print measurement summary
     */
    void print() const;
    
    /**
     * @brief Get all observables as map
     */
    std::unordered_map<std::string, std::pair<Real, Real>> getResults() const;
    
private:
    Real beta_;
    SiteIdx n_sites_;
    
    // Observables
    BinnedObservable n_operators_;
    BinnedObservable energy_;
    BinnedObservable specific_heat_;
    BinnedObservable magnetization_;
    BinnedObservable magnetization_sq_;
    BinnedObservable magnetization_4th_;
    BinnedObservable susceptibility_;
    BinnedObservable stag_magnetization_;
    BinnedObservable binder_ratio_;
    
    // For accumulated measurements within sweep
    Real accum_n_;
    Real accum_n_sq_;
    Real accum_mag_;
    Real accum_mag_sq_;
    Real accum_stag_mag_;
    int n_samples_;
    
    void measureEnergy(const SSEConfig& config);
    void measureMagnetization(const SSEConfig& config);
};

/**
 * @brief Correlation function measurements
 */
class CorrelationMeasurements {
public:
    CorrelationMeasurements(const SSEConfig& config, int n_bins = 100);
    
    /**
     * @brief Measure spin-spin correlations
     */
    void measure(const SSEConfig& config);
    
    /**
     * @brief Get correlation function ⟨S_0 · S_r⟩
     */
    Real correlation(SiteIdx r) const;
    Real correlationError(SiteIdx r) const;
    
    /**
     * @brief Get structure factor S(k)
     */
    Real structureFactor(const std::vector<Real>& k) const;
    
    void startBin() { /* Start bin for all correlations */ }
    void endBin() { /* End bin for all correlations */ }
    void reset();
    
private:
    SiteIdx n_sites_;
    std::vector<BinnedObservable> correlations_;
};

} // namespace sse
