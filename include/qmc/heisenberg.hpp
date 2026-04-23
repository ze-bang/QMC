// SPDX-License-Identifier: MIT
//
// Antiferromagnetic spin-1/2 Heisenberg model parameters.
//
//     H = J  Σ_<ij>  S_i · S_j   ,   J > 0
//
// On a bipartite lattice we apply the Marshall sign transformation and
// expand around the positive-definite bond Hamiltonian
//
//     H = -J Σ_b ( H_{1,b} + H_{2,b} )  +  J · N_b / 4 · I
//
// with
//     H_{1,b} = 1/4 - S^z_i S^z_j         (diagonal,    matrix elt 1/2)
//     H_{2,b} = (1/2)(S+_i S-_j + h.c.)   (off-diagonal,matrix elt 1/2)
//
// Both have the same matrix element 1/2 on antiparallel spins, so we
// sample bonds uniformly with a per-operator weight  J / 2.

#pragma once

#include <cassert>
#include <stdexcept>

#include "qmc/lattice.hpp"
#include "qmc/types.hpp"

namespace qmc {

struct HeisenbergModel {
    Real J = 1.0;            // antiferromagnetic exchange (J > 0)

    // Bond matrix element (same for diagonal and off-diagonal pieces).
    Real bond_matrix_element() const { return 0.5 * J; }

    // Constant energy offset coming from the C = 1/4 shift.
    Real energy_offset(const Lattice& lat) const {
        return 0.25 * J * lat.n_bonds();
    }

    void validate(const Lattice& lat) const {
        if (J <= 0.0) {
            throw std::invalid_argument(
                "HeisenbergModel: J must be positive (this implementation "
                "uses the Marshall sign transformation valid for AF on "
                "bipartite lattices).");
        }
        if (!lat.is_bipartite()) {
            throw std::invalid_argument(
                "HeisenbergModel: lattice is not bipartite -- the loop "
                "update is sign-problem free only on bipartite lattices.");
        }
    }
};

} // namespace qmc
