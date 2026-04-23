// SPDX-License-Identifier: MIT
#include <array>
#include <cmath>

#include "qmc/rng.hpp"
#include "test_runner.hpp"

using qmc::Pcg32;

QMC_TEST(rng_reproducible_with_same_seed) {
    Pcg32 a(42), b(42);
    for (int i = 0; i < 1000; ++i) {
        QMC_REQUIRE(a() == b());
    }
}

QMC_TEST(rng_uniform_in_unit_interval) {
    Pcg32 r(123);
    for (int i = 0; i < 100000; ++i) {
        const double u = r.uniform();
        QMC_REQUIRE(u >= 0.0);
        QMC_REQUIRE(u < 1.0);
    }
}

QMC_TEST(rng_uniform_int_bounded) {
    Pcg32 r(7);
    for (int n : {2, 3, 7, 1023, 1000000}) {
        for (int i = 0; i < 10000; ++i) {
            QMC_REQUIRE(r.uniform_int(n) < static_cast<unsigned>(n));
        }
    }
}

QMC_TEST(rng_uniform_mean_close_to_half) {
    Pcg32 r(2024);
    constexpr int N = 1 << 18;
    double sum = 0.0;
    for (int i = 0; i < N; ++i) sum += r.uniform();
    const double mean = sum / N;
    QMC_REQUIRE_NEAR(mean, 0.5, 5.0e-3);
}

QMC_TEST(rng_chi_squared_uniform_int_bins) {
    Pcg32 r(99);
    constexpr int K = 16;
    std::array<long long, K> bins{};
    constexpr int N = 1 << 18;
    for (int i = 0; i < N; ++i) ++bins[r.uniform_int(K)];
    const double expected = static_cast<double>(N) / K;
    double chi2 = 0.0;
    for (auto b : bins) {
        const double d = b - expected;
        chi2 += d * d / expected;
    }
    // 99% quantile of chi^2 with 15 dof is ~30.6; allow some headroom.
    QMC_REQUIRE(chi2 < 60.0);
}
