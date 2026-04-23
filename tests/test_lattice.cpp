// SPDX-License-Identifier: MIT
#include "qmc/lattice.hpp"
#include "test_runner.hpp"

QMC_TEST(lattice_chain_basic) {
    auto L = qmc::make_chain(8);
    QMC_REQUIRE(L.n_sites() == 8);
    QMC_REQUIRE(L.n_bonds() == 8);             // PBC -> N bonds
    QMC_REQUIRE(L.is_bipartite());
    QMC_REQUIRE_NEAR(L.avg_coordination(), 2.0, 1e-12);
}

QMC_TEST(lattice_square_basic) {
    auto L = qmc::make_square(4, 4);
    QMC_REQUIRE(L.n_sites() == 16);
    QMC_REQUIRE(L.n_bonds() == 32);
    QMC_REQUIRE(L.is_bipartite());
    QMC_REQUIRE_NEAR(L.avg_coordination(), 4.0, 1e-12);
}

QMC_TEST(lattice_honeycomb_basic) {
    auto L = qmc::make_honeycomb(3, 3);
    QMC_REQUIRE(L.n_sites() == 18);
    QMC_REQUIRE(L.n_bonds() == 27);            // 3 bonds per A site
    QMC_REQUIRE(L.is_bipartite());
    QMC_REQUIRE_NEAR(L.avg_coordination(), 3.0, 1e-12);
}

QMC_TEST(lattice_chain_odd_not_bipartite) {
    auto L = qmc::make_chain(5);
    QMC_REQUIRE(!L.is_bipartite());
}

QMC_TEST(lattice_site_bonds_consistent) {
    auto L = qmc::make_square(3, 4);
    const auto& sb = L.site_bonds();
    int total = 0;
    for (auto& v : sb) total += static_cast<int>(v.size());
    QMC_REQUIRE(total == 2 * L.n_bonds());     // each bond contributes to 2 sites
}
