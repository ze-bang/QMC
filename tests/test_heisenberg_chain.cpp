// SPDX-License-Identifier: MIT
//
// Physics regression tests for the 1D AF Heisenberg chain.
//
// Reference values:
//
//   (a) High-temperature limit: for J*beta -> 0 the energy per site
//       vanishes (the trace of the Heisenberg Hamiltonian is zero),
//       and the operator density satisfies <n> -> beta * J * N_b / 4
//       (= beta * E_offset). We check both up to O(beta * J) corrections.
//
//   (b) Ground state of the L=4 ring (closed-form / ED):
//           E_0 = -2 J,  E_0 / N = -0.5    (J = 1)
//       At beta = 8 the SSE result should already be on top of E_0
//       to better than 1 % of |E_0| (which is comfortable for a quick
//       regression test).
//
// We deliberately keep the statistics modest so that the test suite
// finishes in seconds.

#include <cmath>
#include <cstdio>

#include "qmc/qmc.hpp"
#include "test_runner.hpp"

QMC_TEST(heisenberg_chain_high_temperature_limits) {
    using namespace qmc;
    auto lat = make_chain(8);
    HeisenbergModel m{1.0};
    Pcg32 rng(12345);
    const Real beta = 0.05;
    SseEngine eng(lat, m, beta, rng, 16);
    for (int i = 0; i < 2000; ++i) eng.thermalize_step();
    Measurements meas(eng);
    for (int i = 0; i < 20000; ++i) {
        eng.mc_step();
        meas.measure();
    }
    // (i) Energy per site -> 0 as beta -> 0; the leading correction is
    //     ~ -3/16 * z * J^2 * beta = -0.01875 here. With our finite
    //     statistics we just require |E/N| <= 0.05.
    const Real e = meas.energy().mean();
    std::printf("    high-T E/N = %g +/- %g (expected ~ -0.019, |.|<0.05)\n",
                e, meas.energy().stderr_());
    QMC_REQUIRE(std::fabs(e) < 0.05);

    // (ii) <n> -> beta * E_offset = beta * J * N_b / 4 in the high-T
    //      limit. For our parameters this is 0.05 * 8/4 = 0.1.
    const Real n  = meas.n_op().mean();
    const Real n0 = beta * m.energy_offset(lat);
    std::printf("    high-T <n> = %g (expected ~ %g, leading order)\n", n, n0);
    QMC_REQUIRE(std::fabs(n - n0) < 0.05);
}

QMC_TEST(heisenberg_chain_L4_low_temperature) {
    // The 4-site AF Heisenberg ring with PBC has E_0 = -2 J (S_total=0).
    // At beta = 8 the simulation should already be very close to the
    // ground state. We allow a 5% tolerance to keep the test fast.
    using namespace qmc;
    auto lat = make_chain(4);
    HeisenbergModel m{1.0};
    Pcg32 rng(54321);
    const Real beta = 8.0;
    SseEngine eng(lat, m, beta, rng, 32);
    for (int i = 0; i < 5000; ++i) eng.thermalize_step();
    Measurements meas(eng);
    for (int i = 0; i < 30000; ++i) {
        eng.mc_step();
        meas.measure();
    }
    const Real e   = meas.energy().mean();
    const Real err = meas.energy().stderr_();
    const Real expected = -2.0 / 4.0;     // -0.5 per site
    std::printf("    L=4 beta=8 E/N = %g +/- %g  (expected %g)\n",
                e, err, expected);
    QMC_REQUIRE(std::fabs(e - expected) < 0.025);
}

QMC_TEST(heisenberg_chain_uniform_susceptibility_positive) {
    // chi_uniform >= 0 by construction; should also be O(beta) at high T
    // (Curie law) and saturate to a finite value at low T for a finite
    // antiferromagnet -- here we only check non-negativity and finiteness.
    using namespace qmc;
    auto lat = make_chain(8);
    HeisenbergModel m{1.0};
    Pcg32 rng(7);
    SseEngine eng(lat, m, 1.0, rng, 32);
    for (int i = 0; i < 2000; ++i) eng.thermalize_step();
    Measurements meas(eng);
    for (int i = 0; i < 5000; ++i) {
        eng.mc_step();
        meas.measure();
    }
    QMC_REQUIRE(meas.chi_uniform().mean() >= 0.0);
    QMC_REQUIRE(meas.chi_stag().mean()    >  0.0);
    QMC_REQUIRE(std::isfinite(meas.chi_uniform().mean()));
    QMC_REQUIRE(std::isfinite(meas.chi_stag().mean()));
}
