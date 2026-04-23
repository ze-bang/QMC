// SPDX-License-Identifier: MIT
//
// Throughput benchmark. Reports millions of operator-string updates per
// second (a useful figure of merit independent of the lattice geometry).

#include <chrono>
#include <cstdio>
#include <string>

#include "qmc/qmc.hpp"

namespace {

void run(const std::string& kind, int Lx, int Ly, qmc::Real beta, int sweeps) {
    auto lat = qmc::make_lattice(kind, Lx, Ly);
    qmc::HeisenbergModel m{1.0};
    qmc::Pcg32 rng(123456789ULL);
    qmc::SseEngine eng(lat, m, beta, rng, 64);
    for (int i = 0; i < 2000; ++i) eng.thermalize_step();

    const auto t0 = std::chrono::steady_clock::now();
    std::uint64_t total_legs = 0;
    for (int i = 0; i < sweeps; ++i) {
        eng.mc_step();
        total_legs = eng.stats().legs_traversed;
    }
    const auto   t1   = std::chrono::steady_clock::now();
    const double dt   = std::chrono::duration<double>(t1 - t0).count();
    const double Mops = sweeps / dt / 1e6;
    const double Mlegs= total_legs / dt / 1e6;
    std::printf("  %-22s  N=%5d  beta=%5.2f  %6d sweeps  %.2f s   "
                "%.2f Msweeps/s   %.2f Mlegs/s   <n>=%lld\n",
                lat.name().c_str(), lat.n_sites(), beta, sweeps, dt, Mops, Mlegs,
                static_cast<long long>(eng.stats().max_n_observed));
}

} // namespace

int main() {
    std::printf("qmc_sse benchmark v%s\n\n", qmc::kVersion);
    run("chain",     32,  0, 1.0, 20000);
    run("chain",     64,  0, 2.0, 10000);
    run("square",     8,  8, 1.0,  5000);
    run("square",    16, 16, 1.0,  2000);
    run("honeycomb",  6,  6, 1.0,  5000);
    return 0;
}
