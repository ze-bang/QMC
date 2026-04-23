// SPDX-License-Identifier: MIT
//
// Single-band Hubbard model on a generic lattice:
//
//   H = -t  sum_{<ij>,sigma} (c_{i sigma}^+ c_{j sigma} + h.c.)
//       - mu sum_{i,sigma} n_{i sigma}
//       + U  sum_i (n_{i up} - 1/2)(n_{i down} - 1/2)
//
// The symmetric (n - 1/2) form makes the half-filling point sit
// exactly at mu = 0 on bipartite lattices, by particle-hole symmetry.
//
// The kinetic single-particle matrix K_{ij} stored here uses the
// convention
//
//      K_{ij} = -t   if i,j are nearest neighbours,
//      K_{ii} = -mu                ,
//
// so that the non-interacting one-body part of H is  c^+ K c.

#pragma once

#include <stdexcept>

#include "qmc/lattice.hpp"
#include "qmc/linalg.hpp"
#include "qmc/types.hpp"

namespace qmc {

struct HubbardModel {
    Real t  = 1.0;       // hopping
    Real U  = 4.0;       // on-site interaction (U > 0 = repulsive)
    Real mu = 0.0;       // chemical potential (mu = 0 = half-filling on bipartite)

    // Build the Ns x Ns hopping + chemical-potential matrix K.
    la::Matrix kinetic_matrix(const Lattice& lat) const {
        const int Ns = lat.n_sites();
        la::Matrix K(Ns, Ns, 0.0);
        for (int i = 0; i < Ns; ++i) K(i, i) = -mu;
        for (const auto& bd : lat.bonds()) {
            K(bd.i, bd.j) -= t;
            K(bd.j, bd.i) -= t;
        }
        return K;
    }

    void validate() const {
        if (t <= 0.0)
            throw std::invalid_argument("HubbardModel: t must be positive");
        // U == 0 is allowed (free fermions, used as a regression test).
        if (U < 0.0)
            throw std::invalid_argument(
                "HubbardModel: U < 0 (attractive Hubbard) is not supported "
                "by this implementation; the Hirsch HS field would couple "
                "to charge instead of spin.");
    }
};

} // namespace qmc
