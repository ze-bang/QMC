// SPDX-License-Identifier: MIT
//
// Lattice abstraction for spin systems on regular lattices.
//
// A `Lattice` is a flat list of `Bond`s plus a sublattice assignment.
// We support:
//   * Chain      (1D, 2 sites/unit cell? -- no, 1 site/unit cell, bipartite)
//   * Square     (2D, bipartite)
//   * Honeycomb  (2D, bipartite, 2 sites/unit cell)
//
// Building further geometries (triangular, kagome, ...) is a matter of
// constructing a Lattice with a custom bond list. Note: SSE in the
// directed-loop / operator-loop form on _frustrated_ lattices suffers
// from the well-known sign problem; we don't claim correctness there.

#pragma once

#include <array>
#include <cassert>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include "qmc/types.hpp"

namespace qmc {

struct BondDesc {
    Site i;       // first site
    Site j;       // second site
    int  type{0}; // generic bond label (e.g. NN vs NNN), unused by core SSE
};

class Lattice {
public:
    Lattice() = default;

    Lattice(std::string name,
            int n_sites,
            std::vector<BondDesc> bonds,
            std::vector<std::int8_t> sublattice = {})
        : name_(std::move(name)),
          n_sites_(n_sites),
          bonds_(std::move(bonds)),
          sublattice_(std::move(sublattice)) {
        if (sublattice_.empty()) sublattice_.assign(n_sites_, 0);
        if (static_cast<int>(sublattice_.size()) != n_sites_) {
            throw std::invalid_argument("Lattice: sublattice size mismatch");
        }
        validate_();
    }

    const std::string&             name()      const { return name_; }
    int                            n_sites()   const { return n_sites_; }
    int                            n_bonds()   const { return static_cast<int>(bonds_.size()); }
    const std::vector<BondDesc>&   bonds()     const { return bonds_; }
    const BondDesc&                bond(Bond b)const { return bonds_[b]; }
    const std::vector<std::int8_t>& sublattice() const { return sublattice_; }

    // Coordination number averaged over sites.
    Real avg_coordination() const {
        return n_sites_ ? 2.0 * static_cast<Real>(n_bonds()) / n_sites_ : 0.0;
    }

    bool is_bipartite() const {
        for (auto sub : sublattice_) {
            if (sub != 0 && sub != 1) return false;
        }
        for (const auto& b : bonds_) {
            if (sublattice_[b.i] == sublattice_[b.j]) return false;
        }
        return true;
    }

    // Build the per-site bond adjacency list. Each entry is the bond
    // index incident on that site. Built lazily.
    const std::vector<std::vector<Bond>>& site_bonds() const {
        if (site_bonds_.empty()) build_site_bonds_();
        return site_bonds_;
    }

private:
    std::string                       name_;
    int                               n_sites_{0};
    std::vector<BondDesc>             bonds_;
    std::vector<std::int8_t>          sublattice_;
    mutable std::vector<std::vector<Bond>> site_bonds_;

    void validate_() const {
        for (std::size_t b = 0; b < bonds_.size(); ++b) {
            const auto& bd = bonds_[b];
            if (bd.i < 0 || bd.j < 0 || bd.i >= n_sites_ || bd.j >= n_sites_) {
                throw std::invalid_argument(
                    "Lattice: bond " + std::to_string(b) +
                    " refers to out-of-range site");
            }
            if (bd.i == bd.j) {
                throw std::invalid_argument(
                    "Lattice: self-loop bond " + std::to_string(b));
            }
        }
    }

    void build_site_bonds_() const {
        site_bonds_.assign(n_sites_, {});
        for (Bond b = 0; b < n_bonds(); ++b) {
            site_bonds_[bonds_[b].i].push_back(b);
            site_bonds_[bonds_[b].j].push_back(b);
        }
    }
};

// -----------------------------------------------------------------------------
// Lattice factories
// -----------------------------------------------------------------------------

// 1D chain with periodic boundary conditions (length L, even => bipartite).
inline Lattice make_chain(int L) {
    if (L < 2) throw std::invalid_argument("make_chain: need L >= 2");
    std::vector<BondDesc>     bonds;
    bonds.reserve(L);
    for (int i = 0; i < L; ++i) bonds.push_back({i, (i + 1) % L});
    std::vector<std::int8_t> sub(L);
    for (int i = 0; i < L; ++i) sub[i] = static_cast<std::int8_t>(i & 1);
    return Lattice("chain[L=" + std::to_string(L) + "]", L, std::move(bonds),
                   std::move(sub));
}

// 2D square lattice with periodic BCs (Lx * Ly sites, bipartite when Lx,Ly even).
inline Lattice make_square(int Lx, int Ly) {
    if (Lx < 2 || Ly < 2)
        throw std::invalid_argument("make_square: need Lx,Ly >= 2");
    auto idx = [Lx, Ly](int x, int y) { return ((y % Ly) * Lx + (x % Lx)); };
    const int N = Lx * Ly;
    std::vector<BondDesc> bonds;
    bonds.reserve(2 * N);
    for (int y = 0; y < Ly; ++y) {
        for (int x = 0; x < Lx; ++x) {
            bonds.push_back({idx(x, y), idx(x + 1, y)});
            bonds.push_back({idx(x, y), idx(x, y + 1)});
        }
    }
    std::vector<std::int8_t> sub(N);
    for (int y = 0; y < Ly; ++y)
        for (int x = 0; x < Lx; ++x)
            sub[idx(x, y)] = static_cast<std::int8_t>((x + y) & 1);
    return Lattice("square[" + std::to_string(Lx) + "x" + std::to_string(Ly) + "]",
                   N, std::move(bonds), std::move(sub));
}

// 2D honeycomb lattice with PBC. 2 sites per unit cell.
//   Unit cell sites: A = (x,y,0), B = (x,y,1)
//   A connects to: B(x,y), B(x-1,y), B(x,y-1)
inline Lattice make_honeycomb(int Lx, int Ly) {
    if (Lx < 2 || Ly < 2)
        throw std::invalid_argument("make_honeycomb: need Lx,Ly >= 2");
    auto idx = [Lx, Ly](int x, int y, int s) {
        return 2 * ((y % Ly) * Lx + (x % Lx)) + s;
    };
    auto wrap = [](int v, int L) { return (v % L + L) % L; };
    const int N = 2 * Lx * Ly;
    std::vector<BondDesc> bonds;
    bonds.reserve(3 * Lx * Ly);
    for (int y = 0; y < Ly; ++y) {
        for (int x = 0; x < Lx; ++x) {
            const int A = idx(x, y, 0);
            bonds.push_back({A, idx(x, y, 1)});
            bonds.push_back({A, idx(wrap(x - 1, Lx), y, 1)});
            bonds.push_back({A, idx(x, wrap(y - 1, Ly), 1)});
        }
    }
    std::vector<std::int8_t> sub(N);
    for (int i = 0; i < N; ++i) sub[i] = static_cast<std::int8_t>(i & 1);
    return Lattice("honeycomb[" + std::to_string(Lx) + "x" + std::to_string(Ly) + "]",
                   N, std::move(bonds), std::move(sub));
}

// Dispatch from a string name (used by config file driver).
inline Lattice make_lattice(const std::string& kind, int Lx, int Ly) {
    if (kind == "chain")     return make_chain(Lx);
    if (kind == "square")    return make_square(Lx, Ly);
    if (kind == "honeycomb") return make_honeycomb(Lx, Ly);
    throw std::invalid_argument("Unknown lattice kind: " + kind);
}

} // namespace qmc
