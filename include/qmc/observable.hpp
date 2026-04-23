// SPDX-License-Identifier: MIT
//
// Observable accumulator with logarithmic binning analysis.
//
// We use the Flyvbjerg-Petersen blocking method: incoming samples are
// pushed into bin level 0. Whenever level k has two samples they are
// averaged and the result is propagated up to level k+1; the
// integrated autocorrelation time becomes visible as a *plateau* in
// the variance estimate as a function of bin level.
//
// Concretely, for level k with B_k samples, the standard error
//
//     err_k = sqrt( var_k / (B_k - 1) )
//
// rises with k as long as samples on level k are correlated, and
// flattens off once the bin size exceeds the autocorrelation time.
// We report the maximum err_k across levels with B_k >= 8 as the
// observable's error bar -- a robust estimator that does not require
// the user to specify a bin size in advance.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

#include "qmc/types.hpp"

namespace qmc {

class Observable {
public:
    explicit Observable(std::string name = "obs") : name_(std::move(name)) {}

    void clear() {
        levels_.clear();
        partials_.clear();
        n_samples_ = 0;
    }

    // Add a single sample.
    void add(Real x) {
        ++n_samples_;
        std::size_t k = 0;
        Real sample   = x;
        while (true) {
            if (k >= levels_.size()) {
                levels_.emplace_back();
                partials_.push_back({false, 0.0});
            }
            levels_[k].push(sample);
            auto& part = partials_[k];
            if (!part.has) {
                part.has   = true;
                part.value = sample;
                break;
            }
            sample     = 0.5 * (part.value + sample);
            part.has   = false;
            part.value = 0.0;
            ++k;
        }
    }

    const std::string& name() const { return name_; }
    std::uint64_t      count() const { return n_samples_; }
    std::size_t        n_levels() const { return levels_.size(); }

    Real mean() const {
        return levels_.empty() ? 0.0 : levels_[0].mean();
    }

    // Standard error using the maximum across well-populated bin levels.
    Real stderr_() const {
        Real best = 0.0;
        for (std::size_t k = 0; k < levels_.size(); ++k) {
            const auto& L = levels_[k];
            if (L.count < 8) break;
            const Real e = L.std_error();
            if (e > best) best = e;
        }
        if (best == 0.0 && !levels_.empty()) best = levels_[0].std_error();
        return best;
    }

    // Estimated integrated autocorrelation time from the ratio of the
    // converged blocked-error squared to the level-0 error squared.
    Real tau_int() const {
        if (levels_.empty() || levels_[0].count < 16) return 1.0;
        const Real e0 = levels_[0].std_error();
        const Real e  = stderr_();
        if (e0 <= 0.0) return 1.0;
        const Real r = e / e0;
        // err_k = err_0 * sqrt(2 * tau_int) -> tau = (err_k / err_0)^2 / 2
        return std::max(0.5, 0.5 * r * r);
    }

    void format(std::ostream& os) const {
        os << std::setw(28) << std::left << name_
           << std::setw(16) << std::scientific << std::setprecision(8) << mean()
           << " +/- " << std::setw(12) << std::scientific << std::setprecision(3) << stderr_()
           << "  tau_int=" << std::fixed << std::setprecision(2) << tau_int()
           << "  N=" << count();
    }

    std::string summary() const {
        std::ostringstream os;
        format(os);
        return os.str();
    }

private:
    struct Level {
        std::uint64_t count{0};
        Real          sum{0.0};
        Real          sum_sq{0.0};

        void push(Real x) {
            ++count;
            sum    += x;
            sum_sq += x * x;
        }
        Real mean() const { return count ? sum / static_cast<Real>(count) : 0.0; }
        Real variance() const {
            if (count < 2) return 0.0;
            const Real m = mean();
            const Real v = sum_sq / static_cast<Real>(count) - m * m;
            return v < 0.0 ? 0.0 : v;
        }
        Real std_error() const {
            if (count < 2) return 0.0;
            return std::sqrt(variance() / static_cast<Real>(count - 1));
        }
    };
    struct Partial {
        bool has{false};
        Real value{0.0};
    };

    std::string          name_;
    std::vector<Level>   levels_;
    std::vector<Partial> partials_;
    std::uint64_t        n_samples_{0};
};

// Compute (mean(f) - mean(g)^2) style derived quantities with
// jackknife errors. The jackknife requires raw samples, which we hold
// in a small companion class:

class JackknifeSamples {
public:
    void add(Real x) { x_.push_back(x); }
    void clear()     { x_.clear(); }
    std::size_t size() const { return x_.size(); }
    const std::vector<Real>& data() const { return x_; }

    // Generic jackknife: f is a callable that takes a (vector of jack
    // means for each component) and returns the derived quantity.
    template <class F>
    static std::pair<Real, Real> jackknife(
        const std::vector<JackknifeSamples*>& sources, F&& f) {
        const std::size_t k = sources.size();
        const std::size_t n = sources[0]->size();
        for (auto* s : sources) {
            if (s->size() != n) {
                throw std::runtime_error("jackknife: mismatched sample counts");
            }
        }
        if (n < 2) return {0.0, 0.0};
        std::vector<Real> sums(k, 0.0);
        for (std::size_t c = 0; c < k; ++c) {
            for (Real v : sources[c]->data()) sums[c] += v;
        }
        std::vector<Real> jack_means(k);
        for (std::size_t c = 0; c < k; ++c) jack_means[c] = sums[c] / static_cast<Real>(n);
        const Real full = f(jack_means);

        // Jackknife replicates.
        Real mean_jack = 0.0;
        Real var_jack  = 0.0;
        std::vector<Real> jm(k);
        std::vector<Real> reps(n);
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t c = 0; c < k; ++c) {
                jm[c] = (sums[c] - sources[c]->data()[i]) / static_cast<Real>(n - 1);
            }
            const Real r = f(jm);
            reps[i]    = r;
            mean_jack += r;
        }
        mean_jack /= static_cast<Real>(n);
        for (Real r : reps) {
            const Real d = r - mean_jack;
            var_jack += d * d;
        }
        const Real err = std::sqrt(static_cast<Real>(n - 1) / static_cast<Real>(n) * var_jack);
        // Bias-corrected jackknife mean (rarely needed but cheap):
        const Real biased_mean = full - static_cast<Real>(n - 1) * (mean_jack - full);
        return {biased_mean, err};
    }

private:
    std::vector<Real> x_;
};

} // namespace qmc
