// SPDX-License-Identifier: MIT
#include "qmc/lattice.hpp"
#include "qmc/operator_string.hpp"
#include "test_runner.hpp"

QMC_TEST(operator_string_pack_unpack) {
    using namespace qmc;
    for (Bond b : {0, 1, 17, 1023, 1 << 20}) {
        for (int t : {0, 1}) {
            const OpCode op = pack_op(b, t);
            QMC_REQUIRE(op_bond(op) == b);
            QMC_REQUIRE(op_type(op) == t);
            QMC_REQUIRE(!op_is_identity(op));
            QMC_REQUIRE(op_is_diagonal(op) == (t == 0));
            QMC_REQUIRE(op_is_offdiag(op)  == (t == 1));
        }
    }
    QMC_REQUIRE(op_is_identity(kIdentity));
}

QMC_TEST(operator_string_n_op_tracking) {
    using namespace qmc;
    OperatorString os;
    os.resize(8);
    QMC_REQUIRE(os.n_op() == 0);
    os.set(0, pack_op(0, 0));
    os.set(3, pack_op(1, 1));
    os.set(7, pack_op(2, 0));
    QMC_REQUIRE(os.n_op() == 3);
    os.set(3, kIdentity);
    QMC_REQUIRE(os.n_op() == 2);
    os.set(7, pack_op(2, 1));      // type change, not removal
    QMC_REQUIRE(os.n_op() == 2);
}

QMC_TEST(linked_vertices_chain_simple) {
    using namespace qmc;
    auto lat = make_chain(4);
    OperatorString os;
    os.resize(4);
    // Insert two diagonal operators on bonds (0,1) and (2,3).
    os.set(0, pack_op(0, 0));
    os.set(2, pack_op(2, 0));

    LinkedVertices lvl;
    SpinConfig spins(4, +1);
    lvl.build(os, lat, spins);
    QMC_REQUIRE(lvl.n_vertices() == 2);
    // Each vertex has 4 legs; total = 8.
    QMC_REQUIRE((int)lvl.link().size() == 8);
}
