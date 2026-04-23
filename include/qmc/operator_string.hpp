// SPDX-License-Identifier: MIT
//
// SSE operator string and linked-vertex list.
//
// The operator string is a length-M array of `OpCode`s; each non-identity
// entry is one bond operator acting at imaginary time index p in [0, M).
// `n` counts non-identity entries.
//
// The linked-vertex list is the data structure that powers the
// operator-loop update (Sandvik 1999, "Stochastic series expansion method
// with operator-loop cluster updates", PRB 59, R14157). For each non-
// identity operator at position p we have 4 legs:
//
//        leg2          leg3
//         |             |          (top of vertex,    spins after operator)
//         +------(p)----+
//         |             |          (bottom of vertex, spins before operator)
//        leg0          leg1
//
// In the flattened arrays:
//     vertex[v].leg[k]   for k in {0,1,2,3}   k = (top<<1) | which_site
// We use a global integer index `4*v + k` for legs, and store
//     X[4*v + k]  = 4*v' + k'   (the leg connected to this one)

#pragma once

#include <cassert>
#include <cstdint>
#include <vector>

#include "qmc/lattice.hpp"
#include "qmc/types.hpp"

namespace qmc {

class OperatorString {
public:
    OperatorString() = default;

    void resize(Length M) {
        ops_.assign(static_cast<std::size_t>(M), kIdentity);
    }

    Length size() const noexcept { return static_cast<Length>(ops_.size()); }
    Length n_op() const noexcept { return n_op_; }

    OpCode  operator[](Length p) const noexcept { return ops_[p]; }
    OpCode& operator[](Length p)       noexcept { return ops_[p]; }

    void set(Length p, OpCode op) noexcept {
        OpCode old = ops_[p];
        if (old == kIdentity && op != kIdentity) ++n_op_;
        else if (old != kIdentity && op == kIdentity) --n_op_;
        ops_[p] = op;
    }

    // Adaptive growth of M during thermalization. Aims for M slightly
    // larger than the largest n observed, so that the n / (M - n)
    // balance in the diagonal update is not bottlenecked by truncation.
    bool grow_if_needed(Real safety = 1.3) {
        const Length needed =
            static_cast<Length>(static_cast<Real>(n_op_) * safety + 4.0);
        if (needed <= size()) return false;
        ops_.resize(static_cast<std::size_t>(needed), kIdentity);
        return true;
    }

    const std::vector<OpCode>& data() const noexcept { return ops_; }

    void clear() {
        ops_.clear();
        n_op_ = 0;
    }

private:
    std::vector<OpCode> ops_;
    Length              n_op_{0};
};

// Linked-vertex list -------------------------------------------------------
//
// We rebuild the linked list from scratch before every loop sweep. This
// is O(M + N_sites) work and keeps the data structures cache-friendly.

class LinkedVertices {
public:
    void build(const OperatorString& ops, const Lattice& lat,
               const SpinConfig& spins);

    // Number of non-identity operators in the current build (= N_op).
    int n_vertices() const noexcept { return n_vertices_; }

    // Position p in the operator string of vertex v.
    Length vertex_position(int v) const { return v_pos_[v]; }

    // The leg array: `link_[leg]` is the leg connected to `leg`, or -1
    // for legs that loop back to themselves (free spin tracks).
    std::vector<std::int32_t>&       link()       noexcept { return link_; }
    const std::vector<std::int32_t>& link() const noexcept { return link_; }

    // First / last leg encountered for each site (used to wrap the
    // periodic time direction). -1 if the site has no operators.
    const std::vector<std::int32_t>& first_leg() const noexcept { return first_; }
    const std::vector<std::int32_t>& last_leg () const noexcept { return last_;  }

    // Map vertex -> position in the operator string.
    const std::vector<Length>& v_pos() const noexcept { return v_pos_; }

private:
    int                          n_vertices_{0};
    std::vector<Length>          v_pos_;     // v -> p
    std::vector<std::int32_t>    link_;      // 4*Nv legs
    std::vector<std::int32_t>    first_;     // per site, first leg encountered
    std::vector<std::int32_t>    last_;      // per site, last  leg encountered
};

inline void LinkedVertices::build(const OperatorString& ops,
                                  const Lattice& lat,
                                  const SpinConfig& /*spins*/) {
    const Length M  = ops.size();
    const int    Ns = lat.n_sites();

    v_pos_.clear();
    v_pos_.reserve(static_cast<std::size_t>(ops.n_op()));
    for (Length p = 0; p < M; ++p) {
        if (!op_is_identity(ops[p])) v_pos_.push_back(p);
    }
    n_vertices_ = static_cast<int>(v_pos_.size());

    link_.assign(static_cast<std::size_t>(4 * n_vertices_), -1);
    first_.assign(Ns, -1);
    last_ .assign(Ns, -1);

    // Track the most recent unmatched leg per site as we sweep upward.
    std::vector<std::int32_t> last_leg_at_site(Ns, -1);

    for (int v = 0; v < n_vertices_; ++v) {
        const Length p     = v_pos_[v];
        const OpCode op    = ops[p];
        const Bond   b     = op_bond(op);
        const auto&  bd    = lat.bond(b);
        const Site   s0    = bd.i;
        const Site   s1    = bd.j;
        const std::int32_t leg_b0 = 4 * v + 0;   // bottom-left  (site s0, in)
        const std::int32_t leg_b1 = 4 * v + 1;   // bottom-right (site s1, in)
        const std::int32_t leg_t0 = 4 * v + 2;   // top-left     (site s0, out)
        const std::int32_t leg_t1 = 4 * v + 3;   // top-right    (site s1, out)

        // Connect bottom legs to the previous unmatched top leg at each
        // site (or record as "first" if none yet).
        const std::int32_t prev0 = last_leg_at_site[s0];
        if (prev0 == -1) {
            first_[s0] = leg_b0;
        } else {
            link_[prev0] = leg_b0;
            link_[leg_b0] = prev0;
        }
        last_leg_at_site[s0] = leg_t0;

        const std::int32_t prev1 = last_leg_at_site[s1];
        if (prev1 == -1) {
            first_[s1] = leg_b1;
        } else {
            link_[prev1] = leg_b1;
            link_[leg_b1] = prev1;
        }
        last_leg_at_site[s1] = leg_t1;
    }

    // Close the cyclic boundary in imaginary time.
    for (Site s = 0; s < Ns; ++s) {
        const std::int32_t f = first_[s];
        const std::int32_t l = last_leg_at_site[s];
        if (f != -1) {
            link_[l] = f;
            link_[f] = l;
            last_[s] = l;
        }
    }
}

} // namespace qmc
