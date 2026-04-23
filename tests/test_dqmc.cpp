// SPDX-License-Identifier: MIT
//
// Physics regression tests for the DQMC engine.
//
// We exploit the fact that at U = 0 the Hubbard-Stratonovich auxiliary
// field decouples (alpha = 0) and the simulation samples the exact
// non-interacting Green's function
//
//     G_sigma(i, j)  =  [ (I + e^{-beta K})^{-1} ]_{ij}.
//
// This is a *deterministic* check: no statistics are needed because
// every aux-field configuration gives identically the same Green's
// function. We therefore compare the engine's Green's function to the
// closed-form value to high precision.
//
// We also exercise a small interacting cluster (4-site Hubbard ring at
// half filling) and check that the half-filling sum rules
//     <n> = 1,    <n_up n_down> in [0, 1/2]
// hold within reasonable statistical bounds.

#include <cmath>
#include <cstdio>

#include "qmc/qmc.hpp"
#include "test_runner.hpp"

QMC_TEST(dqmc_free_fermion_green_matches_closed_form) {
    using namespace qmc;
    auto lat = make_chain(8);
    HubbardModel model{1.0, 0.0, 0.0};   // U = 0, mu = 0
    DqmcParams p; p.beta = 4.0; p.L_tau = 40; p.n_stab = 10;
    Pcg32 rng(101);
    DqmcEngine eng(lat, model, p, rng);

    // Reference  G_ref = (I + exp(-beta K))^{-1}.
    auto K = model.kinetic_matrix(lat);
    auto Bbeta = la::expm_sym(K, p.beta);
    la::Matrix IpB = Bbeta;
    for (int i = 0; i < lat.n_sites(); ++i) IpB(i, i) += 1.0;
    auto G_ref = la::inverse(IpB);

    // The engine's G_up at slice 0 should equal G_ref to numerical
    // precision (alpha = 0 -> aux field has no effect).
    const auto& Gu = eng.green_up();
    double err = 0.0;
    for (int i = 0; i < lat.n_sites(); ++i)
        for (int j = 0; j < lat.n_sites(); ++j)
            err = std::max(err, std::fabs(Gu(i, j) - G_ref(i, j)));
    std::printf("    free-fermion |G - G_ref|_max = %.3e\n", err);
    QMC_REQUIRE(err < 1e-9);

    // Run a sweep -- with U = 0 the Green's function should be
    // *invariant* (every Sherman-Morrison call has Delta = 0, but the
    // wrap+stabilize cycle still happens and must not corrupt G).
    eng.sweep();
    err = 0.0;
    for (int i = 0; i < lat.n_sites(); ++i)
        for (int j = 0; j < lat.n_sites(); ++j)
            err = std::max(err, std::fabs(eng.green_up()(i, j) - G_ref(i, j)));
    std::printf("    free-fermion drift after 1 sweep = %.3e\n", err);
    QMC_REQUIRE(err < 1e-9);
}

QMC_TEST(dqmc_free_fermion_density_at_half_filling) {
    // For U = 0, mu = 0 on a bipartite lattice the density is exactly 1.
    using namespace qmc;
    auto lat = make_square(2, 2);
    HubbardModel model{1.0, 0.0, 0.0};
    DqmcParams p; p.beta = 2.0; p.L_tau = 20; p.n_stab = 5;
    Pcg32 rng(7);
    DqmcEngine eng(lat, model, p, rng);
    DqmcMeasurements meas(eng);
    for (int i = 0; i < 30; ++i) eng.sweep();
    for (int i = 0; i < 100; ++i) { eng.sweep(); meas.measure(); }
    const Real n = meas.density().mean();
    std::printf("    free-fermion <n> = %g (expected 1.0)\n", n);
    QMC_REQUIRE_NEAR(n, 1.0, 1e-9);   // deterministic at U = 0
}

QMC_TEST(dqmc_half_filled_repulsive_density_one) {
    // Repulsive Hubbard at mu = 0 on bipartite lattice: <n> = 1 by
    // particle-hole symmetry, regardless of U. This is a powerful
    // consistency check on the entire DQMC pipeline.
    using namespace qmc;
    auto lat = make_chain(4);
    HubbardModel model{1.0, 4.0, 0.0};
    DqmcParams p; p.beta = 2.0; p.L_tau = 20; p.n_stab = 5;
    Pcg32 rng(2024);
    DqmcEngine eng(lat, model, p, rng);
    DqmcMeasurements meas(eng);
    for (int i = 0; i < 100; ++i) eng.sweep();
    for (int i = 0; i < 400; ++i) { eng.sweep(); meas.measure(); }
    const Real n   = meas.density().mean();
    const Real err = meas.density().stderr_();
    std::printf("    repulsive U=4, beta=2 chain L=4: <n> = %g +/- %g\n",
                n, err);
    QMC_REQUIRE(std::fabs(n - 1.0) < 0.02);  // PH-symmetric, robust
    // Double occupancy lies in (0, 1/2).
    const Real d = meas.double_occ().mean();
    std::printf("    double occupancy = %g\n", d);
    QMC_REQUIRE(d > 0.0);
    QMC_REQUIRE(d < 0.5);
    // Average sign should be exactly +1 at half filling on a bipartite
    // lattice (sign-problem free).
    QMC_REQUIRE_NEAR(meas.sign().mean(), 1.0, 1e-12);
}
