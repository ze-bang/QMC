// SPDX-License-Identifier: MIT
#include <cmath>

#include "qmc/observable.hpp"
#include "qmc/rng.hpp"
#include "test_runner.hpp"

using qmc::Observable;
using qmc::JackknifeSamples;
using qmc::Pcg32;

QMC_TEST(observable_iid_normal_mean_and_error) {
    Pcg32 r(123);
    Observable o("x");
    constexpr int N = 1 << 16;
    // Box-Muller for N(0,1).
    for (int i = 0; i < N; i += 2) {
        const double u1 = std::max(1e-12, r.uniform());
        const double u2 = r.uniform();
        const double z1 = std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * M_PI * u2);
        const double z2 = std::sqrt(-2.0 * std::log(u1)) * std::sin(2.0 * M_PI * u2);
        o.add(z1);
        o.add(z2);
    }
    QMC_REQUIRE(std::fabs(o.mean()) < 0.05);
    // For iid N(0,1) the standard error is 1/sqrt(N) ~= 3.9e-3.
    QMC_REQUIRE(o.stderr_() < 0.02);
    // Autocorrelation should look like ~1.
    QMC_REQUIRE(o.tau_int() < 5.0);
}

QMC_TEST(observable_correlated_increases_tau_int) {
    Pcg32 r(42);
    Observable o("ar1");
    constexpr int N   = 1 << 16;
    constexpr double a = 0.9;     // AR(1) with strong correlation
    double x = 0.0;
    for (int i = 0; i < N; ++i) {
        const double u1 = std::max(1e-12, r.uniform());
        const double u2 = r.uniform();
        const double z = std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * M_PI * u2);
        x = a * x + std::sqrt(1 - a * a) * z;
        o.add(x);
    }
    // Theoretical tau_int = (1+a)/(1-a) = 19 for a=0.9.
    QMC_REQUIRE(o.tau_int() > 5.0);
}

QMC_TEST(jackknife_variance_recovers_naive_estimator) {
    JackknifeSamples a, b;
    Pcg32 r(7);
    constexpr int N = 5000;
    for (int i = 0; i < N; ++i) {
        const double u = r.uniform();
        a.add(u);
        b.add(u);
    }
    // f = <a> - <b>^2 with mean values <a>=<b>=0.5 -> f = 0.25.
    auto [val, err] = JackknifeSamples::jackknife({&a, &b}, [](const auto& m) {
        return m[0] - m[1] * m[1];
    });
    QMC_REQUIRE_NEAR(val, 0.25, 5.0e-3);
    QMC_REQUIRE(err > 0.0);
    QMC_REQUIRE(err < 0.05);
}
