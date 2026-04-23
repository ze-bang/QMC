// SPDX-License-Identifier: MIT
//
// Stochastic Series Expansion (SSE) engine with operator-loop updates
// for the spin-1/2 antiferromagnetic Heisenberg model on bipartite
// lattices.
//
// Reference algorithm:
//   A. W. Sandvik, "Stochastic series expansion method with operator-loop
//   cluster updates", Phys. Rev. B 59, R14157 (1999).
//
// One Monte Carlo step consists of:
//   1) a *diagonal* update that sweeps the operator string and proposes
//      insertion / removal of diagonal bond operators with the standard
//      Metropolis acceptance ratios derived from the SSE weight,
//   2) building the *linked-vertex* list and performing an
//      *operator-loop* update which traverses each loop once and flips
//      the spins along it with probability 1/2,
//   3) flipping the residual "free" spins (sites with no operators) with
//      probability 1/2.
//
// Step (2) is responsible for the *off-diagonal* updates that make this
// algorithm ergodic and very efficient (typical autocorrelation times of
// a few sweeps even at the critical point of the 2D AF Heisenberg model).

#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "qmc/heisenberg.hpp"
#include "qmc/lattice.hpp"
#include "qmc/operator_string.hpp"
#include "qmc/rng.hpp"
#include "qmc/types.hpp"

namespace qmc {

struct SseStats {
    std::uint64_t diag_inserts_proposed{0};
    std::uint64_t diag_inserts_accepted{0};
    std::uint64_t diag_removes_proposed{0};
    std::uint64_t diag_removes_accepted{0};
    std::uint64_t loops_constructed{0};
    std::uint64_t loops_flipped{0};
    std::uint64_t legs_traversed{0};
    Length        max_n_observed{0};
    Length        cutoff_M{0};

    Real diag_insert_acc() const {
        return diag_inserts_proposed
                   ? static_cast<Real>(diag_inserts_accepted) / diag_inserts_proposed
                   : 0.0;
    }
    Real diag_remove_acc() const {
        return diag_removes_proposed
                   ? static_cast<Real>(diag_removes_accepted) / diag_removes_proposed
                   : 0.0;
    }
    Real loop_flip_frac() const {
        return loops_constructed
                   ? static_cast<Real>(loops_flipped) / loops_constructed
                   : 0.0;
    }
};

class SseEngine {
public:
    SseEngine(const Lattice& lattice,
              const HeisenbergModel& model,
              Real beta,
              Pcg32 rng,
              Length initial_M = 32);

    // Run one Monte Carlo step (diagonal + loop + free-spin flip).
    void mc_step();

    // Trigger adaptive growth of M and clear acceptance counters.
    // Call between thermalization sweeps.
    void thermalize_step();

    const Lattice&         lattice()  const { return lat_; }
    const HeisenbergModel& model()    const { return model_; }
    const SpinConfig&      spins()    const { return spins_; }
    const OperatorString&  operators()const { return ops_; }
    const SseStats&        stats()    const { return stats_; }
    Real                   beta()     const { return beta_; }
    Pcg32&                 rng()            { return rng_; }

    // Observables that depend purely on the operator string can be read
    // off cheaply via these helpers.
    Length n_op() const { return ops_.n_op(); }

private:
    // --- algorithmic primitives ---
    void diagonal_update_();
    void loop_update_();
    void flip_free_spins_();

    // --- members ---
    const Lattice&         lat_;
    HeisenbergModel        model_;
    Real                   beta_{1.0};
    Pcg32                  rng_;

    SpinConfig             spins_;        // ±1
    OperatorString         ops_;
    LinkedVertices         lvl_;
    SseStats               stats_;

    Real                   bond_weight_{0.5};   // J/2 in units of energy
};

inline SseEngine::SseEngine(const Lattice& lattice,
                            const HeisenbergModel& model,
                            Real beta,
                            Pcg32 rng,
                            Length initial_M)
    : lat_(lattice),
      model_(model),
      beta_(beta),
      rng_(rng) {
    model_.validate(lat_);
    if (beta_ <= 0.0) throw std::invalid_argument("SseEngine: beta must be > 0");
    if (initial_M <= 0) initial_M = 32;
    bond_weight_ = model_.bond_matrix_element();   // = J / 2

    // Random initial spin configuration with zero total magnetization
    // (helps avoid early stalls at very low T).
    spins_.assign(lat_.n_sites(), 0);
    for (int s = 0; s < lat_.n_sites(); ++s) {
        spins_[s] = (s & 1) ? +1 : -1;
    }
    // Random shuffle to break the initial Néel order (loop update will
    // restore it).
    for (int s = lat_.n_sites() - 1; s > 0; --s) {
        const int t = static_cast<int>(rng_.uniform_int(s + 1));
        std::swap(spins_[s], spins_[t]);
    }
    ops_.resize(initial_M);
}

// -----------------------------------------------------------------------------
// Diagonal update
// -----------------------------------------------------------------------------
inline void SseEngine::diagonal_update_() {
    const int    Nb     = lat_.n_bonds();
    const Length M      = ops_.size();
    const Real   bweight = bond_weight_;
    // Working spin state: starts from spins_ and is propagated forward
    // through off-diagonal operators as we sweep.
    std::vector<std::int8_t> state = spins_;

    // Pre-compute the constant probability factors. Sandvik's MH ratios:
    //   P_insert(b) = beta * Nb * <a|H_b|a> / (M - n)
    //   P_remove    = (M - n + 1) / (beta * Nb * <a|H_b|a>)
    // since the matrix element is the same for every diagonal operator
    // (J/2) on antiparallel spins, we can factor it out.
    const Real prefactor_insert_base = beta_ * Nb * bweight;

    for (Length p = 0; p < M; ++p) {
        OpCode op = ops_[p];
        if (op_is_identity(op)) {
            // Try to insert a diagonal operator on a random bond.
            const Bond b = static_cast<Bond>(rng_.uniform_int(static_cast<std::uint32_t>(Nb)));
            const auto& bd = lat_.bond(b);
            // Antiparallel spins required for a non-zero matrix element.
            if (state[bd.i] != state[bd.j]) {
                ++stats_.diag_inserts_proposed;
                const Length denom = M - ops_.n_op();
                const Real   acc   = prefactor_insert_base / static_cast<Real>(denom);
                if (acc >= 1.0 || rng_.uniform() < acc) {
                    ops_.set(p, pack_op(b, 0));
                    ++stats_.diag_inserts_accepted;
                }
            }
        } else if (op_is_diagonal(op)) {
            // Try to remove the diagonal operator.
            ++stats_.diag_removes_proposed;
            const Length denom = M - ops_.n_op() + 1;
            const Real   acc   = static_cast<Real>(denom) / prefactor_insert_base;
            if (acc >= 1.0 || rng_.uniform() < acc) {
                ops_.set(p, kIdentity);
                ++stats_.diag_removes_accepted;
            }
        } else {
            // Off-diagonal operator: propagate the state forward.
            const auto& bd = lat_.bond(op_bond(op));
            state[bd.i] = static_cast<std::int8_t>(-state[bd.i]);
            state[bd.j] = static_cast<std::int8_t>(-state[bd.j]);
        }
    }
    if (ops_.n_op() > stats_.max_n_observed) stats_.max_n_observed = ops_.n_op();
}

// -----------------------------------------------------------------------------
// Loop update
// -----------------------------------------------------------------------------
//
// Operator-loop update for the AF Heisenberg model: the deterministic
// "switch" rule on the SSE vertex graph. When a loop enters a vertex
// through leg `e` it exits through leg `e XOR 1` -- i.e. the leg on
// the *same* time-side (top or bottom) but on the *other* site of the
// bond. The two legs are then jumped from via the linked vertex list:
//
//     leg `e` --(visit)--> leg `e^1` --(link)--> next leg
//
// Each vertex is visited by at most two loops (one per time-side).
// When a loop is flipped:
//   * the spin on the visited bottom legs becomes -spin,
//   * the spin on the visited top legs becomes -spin,
//   * if exactly one side of a vertex is visited the operator's diagonal
//     / off-diagonal type is toggled,
//   * if both sides are visited the type is unchanged.
//
// This leaves the SSE weight invariant, hence each loop is flipped with
// probability 1/2.

inline void SseEngine::loop_update_() {
    lvl_.build(ops_, lat_, spins_);
    auto&       link = lvl_.link();
    const int   Nv   = lvl_.n_vertices();
    if (Nv == 0) return;
    const int   total_legs = 4 * Nv;

    // visited[leg] = 0 untouched, 1 in a flipped loop, 2 in a kept loop.
    std::vector<std::uint8_t> visited(static_cast<std::size_t>(total_legs), 0);

    for (int start_leg = 0; start_leg < total_legs; ++start_leg) {
        if (visited[start_leg]) continue;
        ++stats_.loops_constructed;
        const std::uint8_t mark = rng_.bernoulli(0.5) ? 1u : 2u;
        if (mark == 1u) ++stats_.loops_flipped;

        std::int32_t leg = start_leg;
        do {
            ++stats_.legs_traversed;
            visited[leg]      = mark;
            const std::int32_t exit_leg = leg ^ 1;          // switch rule
            visited[exit_leg] = mark;
            const std::int32_t next = link[exit_leg];
            if (next < 0) break;                            // free track (shouldn't happen here)
            leg = next;
        } while (leg != start_leg);
    }

    // Update operator types: a vertex's type toggles iff exactly one of
    // its time-sides (bottom = legs {0,1}, top = legs {2,3}) was
    // touched by a *flipped* loop.
    for (int v = 0; v < Nv; ++v) {
        const bool bot_flipped = (visited[4 * v + 0] == 1) ||
                                 (visited[4 * v + 1] == 1);
        const bool top_flipped = (visited[4 * v + 2] == 1) ||
                                 (visited[4 * v + 3] == 1);
        if (bot_flipped == top_flipped) continue;           // either both or neither
        const Length p  = lvl_.vertex_position(v);
        const OpCode op = ops_[p];
        ops_.set(p, pack_op(op_bond(op), op_type(op) ^ 1));
    }

    // Update the persistent spin state: spins_[s] is the value at
    // imaginary-time position 0 (= the bottom of the periodic operator
    // string). It is flipped iff the very first leg encountered when
    // sweeping the operators at site `s` belongs to a flipped loop.
    const auto& first_leg = lvl_.first_leg();
    for (Site s = 0; s < lat_.n_sites(); ++s) {
        const std::int32_t fl = first_leg[s];
        if (fl >= 0 && visited[fl] == 1) {
            spins_[s] = static_cast<std::int8_t>(-spins_[s]);
        }
    }
}

// -----------------------------------------------------------------------------
// Free-spin flip
// -----------------------------------------------------------------------------
inline void SseEngine::flip_free_spins_() {
    if (lvl_.n_vertices() == 0) {
        // No operators at all -> all spins are free.
        for (auto& s : spins_) {
            if (rng_.bernoulli(0.5)) s = static_cast<std::int8_t>(-s);
        }
        return;
    }
    const auto& first_leg = lvl_.first_leg();
    for (Site s = 0; s < lat_.n_sites(); ++s) {
        if (first_leg[s] < 0) {
            if (rng_.bernoulli(0.5)) spins_[s] = static_cast<std::int8_t>(-spins_[s]);
        }
    }
}

// -----------------------------------------------------------------------------
// MC steps
// -----------------------------------------------------------------------------
inline void SseEngine::mc_step() {
    diagonal_update_();
    loop_update_();
    flip_free_spins_();
}

inline void SseEngine::thermalize_step() {
    mc_step();
    if (ops_.grow_if_needed(1.4)) {
        stats_.cutoff_M = ops_.size();
    } else {
        stats_.cutoff_M = ops_.size();
    }
}

} // namespace qmc
