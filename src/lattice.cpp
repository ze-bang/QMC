/**
 * @file lattice.cpp
 * @brief Implementation of lattice structures
 */

#include "sse/lattice.hpp"
#include <fmt/format.h>
#include <stdexcept>
#include <algorithm>
#include <iostream>

namespace sse {

Lattice::Lattice(SiteIdx n_sites, const std::vector<Bond>& bonds)
    : n_sites_(n_sites), bonds_(bonds) {
    buildSiteBonds();
}

void Lattice::buildSiteBonds() {
    site_bonds_.resize(n_sites_);
    for (auto& sb : site_bonds_) sb.clear();
    
    n_bond_types_ = 0;
    for (BondIdx b = 0; b < numBonds(); ++b) {
        const auto& bond = bonds_[b];
        if (bond.i >= n_sites_ || bond.j >= n_sites_) {
            throw std::runtime_error(fmt::format(
                "Invalid bond ({}, {}): sites must be < {}", 
                bond.i, bond.j, n_sites_));
        }
        site_bonds_[bond.i].push_back(b);
        site_bonds_[bond.j].push_back(b);
        n_bond_types_ = std::max(n_bond_types_, bond.type + 1);
    }
    
    max_coordination_ = 0;
    for (SiteIdx s = 0; s < n_sites_; ++s) {
        max_coordination_ = std::max(max_coordination_, 
                                     static_cast<int>(site_bonds_[s].size()));
    }
}

std::vector<SiteIdx> Lattice::neighbors(SiteIdx site) const {
    std::vector<SiteIdx> result;
    for (BondIdx b : site_bonds_[site]) {
        const auto& bond = bonds_[b];
        result.push_back(bond.i == site ? bond.j : bond.i);
    }
    return result;
}

std::vector<int> Lattice::siteCoords(SiteIdx site) const {
    if (dimensions_.empty()) return {};
    
    std::vector<int> coords(dimensions_.size());
    SiteIdx s = site;
    for (int d = static_cast<int>(dimensions_.size()) - 1; d >= 0; --d) {
        coords[d] = s % dimensions_[d];
        s /= dimensions_[d];
    }
    return coords;
}

SiteIdx Lattice::coordsToSite(const std::vector<int>& coords) const {
    if (coords.size() != dimensions_.size()) {
        throw std::runtime_error("Coordinate dimension mismatch");
    }
    
    SiteIdx site = 0;
    for (size_t d = 0; d < dimensions_.size(); ++d) {
        site = site * dimensions_[d] + coords[d];
    }
    return site;
}

Lattice Lattice::chain(SiteIdx L, bool periodic) {
    std::vector<Bond> bonds;
    bonds.reserve(L);
    
    SiteIdx n_bonds = periodic ? L : L - 1;
    for (SiteIdx i = 0; i < n_bonds; ++i) {
        bonds.emplace_back(i, (i + 1) % L);
    }
    
    Lattice lat(L, bonds);
    lat.type_ = LatticeType::Chain;
    lat.dimensions_ = {L};
    return lat;
}

Lattice Lattice::square(SiteIdx Lx, SiteIdx Ly, bool periodic) {
    SiteIdx N = Lx * Ly;
    std::vector<Bond> bonds;
    bonds.reserve(2 * N);
    
    auto idx = [Lx](SiteIdx x, SiteIdx y) { return y * Lx + x; };
    
    for (SiteIdx y = 0; y < Ly; ++y) {
        for (SiteIdx x = 0; x < Lx; ++x) {
            // Horizontal bond
            if (periodic || x + 1 < Lx) {
                bonds.emplace_back(idx(x, y), idx((x + 1) % Lx, y));
            }
            // Vertical bond
            if (periodic || y + 1 < Ly) {
                bonds.emplace_back(idx(x, y), idx(x, (y + 1) % Ly));
            }
        }
    }
    
    Lattice lat(N, bonds);
    lat.type_ = LatticeType::Square;
    lat.dimensions_ = {Lx, Ly};
    return lat;
}

Lattice Lattice::triangular(SiteIdx Lx, SiteIdx Ly, bool periodic) {
    SiteIdx N = Lx * Ly;
    std::vector<Bond> bonds;
    bonds.reserve(3 * N);
    
    auto idx = [Lx](SiteIdx x, SiteIdx y) { return y * Lx + x; };
    
    for (SiteIdx y = 0; y < Ly; ++y) {
        for (SiteIdx x = 0; x < Lx; ++x) {
            // Horizontal
            if (periodic || x + 1 < Lx) {
                bonds.emplace_back(idx(x, y), idx((x + 1) % Lx, y));
            }
            // Vertical
            if (periodic || y + 1 < Ly) {
                bonds.emplace_back(idx(x, y), idx(x, (y + 1) % Ly));
            }
            // Diagonal
            if ((periodic || (x + 1 < Lx && y + 1 < Ly))) {
                bonds.emplace_back(idx(x, y), idx((x + 1) % Lx, (y + 1) % Ly));
            }
        }
    }
    
    Lattice lat(N, bonds);
    lat.type_ = LatticeType::Triangular;
    lat.dimensions_ = {Lx, Ly};
    return lat;
}

Lattice Lattice::honeycomb(SiteIdx Lx, SiteIdx Ly, bool periodic) {
    // 2 sites per unit cell
    SiteIdx N = 2 * Lx * Ly;
    std::vector<Bond> bonds;
    bonds.reserve(3 * Lx * Ly);
    
    // Site indexing: (x, y, sublattice) -> site index
    auto idx = [Lx, Ly](SiteIdx x, SiteIdx y, int sub) {
        return 2 * ((y % Ly) * Lx + (x % Lx)) + sub;
    };
    
    for (SiteIdx y = 0; y < Ly; ++y) {
        for (SiteIdx x = 0; x < Lx; ++x) {
            // A-B bond within unit cell
            bonds.emplace_back(idx(x, y, 0), idx(x, y, 1));
            
            // A(x,y) to B(x-1, y)
            if (periodic || x > 0) {
                bonds.emplace_back(idx(x, y, 0), idx((x - 1 + Lx) % Lx, y, 1));
            }
            
            // A(x,y) to B(x, y-1)
            if (periodic || y > 0) {
                bonds.emplace_back(idx(x, y, 0), idx(x, (y - 1 + Ly) % Ly, 1));
            }
        }
    }
    
    Lattice lat(N, bonds);
    lat.type_ = LatticeType::Honeycomb;
    lat.dimensions_ = {Lx, Ly};
    return lat;
}

Lattice Lattice::kagome(SiteIdx Lx, SiteIdx Ly, bool periodic) {
    // 3 sites per unit cell
    SiteIdx N = 3 * Lx * Ly;
    std::vector<Bond> bonds;
    bonds.reserve(6 * Lx * Ly);
    
    auto idx = [Lx, Ly](SiteIdx x, SiteIdx y, int sub) {
        return 3 * ((y % Ly) * Lx + (x % Lx)) + sub;
    };
    
    for (SiteIdx y = 0; y < Ly; ++y) {
        for (SiteIdx x = 0; x < Lx; ++x) {
            // Bonds within unit cell (triangle)
            bonds.emplace_back(idx(x, y, 0), idx(x, y, 1));
            bonds.emplace_back(idx(x, y, 1), idx(x, y, 2));
            bonds.emplace_back(idx(x, y, 2), idx(x, y, 0));
            
            // Bonds to neighboring cells
            if (periodic || x + 1 < Lx) {
                bonds.emplace_back(idx(x, y, 1), idx((x + 1) % Lx, y, 0));
            }
            if (periodic || y + 1 < Ly) {
                bonds.emplace_back(idx(x, y, 2), idx(x, (y + 1) % Ly, 0));
            }
            if (periodic || (x + 1 < Lx && y + 1 < Ly)) {
                bonds.emplace_back(idx(x, y, 2), idx((x + 1) % Lx, (y + 1) % Ly, 1));
            }
        }
    }
    
    Lattice lat(N, bonds);
    lat.type_ = LatticeType::Kagome;
    lat.dimensions_ = {Lx, Ly};
    return lat;
}

Lattice Lattice::cubic(SiteIdx Lx, SiteIdx Ly, SiteIdx Lz, bool periodic) {
    SiteIdx N = Lx * Ly * Lz;
    std::vector<Bond> bonds;
    bonds.reserve(3 * N);
    
    auto idx = [Lx, Ly](SiteIdx x, SiteIdx y, SiteIdx z) {
        return z * Lx * Ly + y * Lx + x;
    };
    
    for (SiteIdx z = 0; z < Lz; ++z) {
        for (SiteIdx y = 0; y < Ly; ++y) {
            for (SiteIdx x = 0; x < Lx; ++x) {
                // x-direction
                if (periodic || x + 1 < Lx) {
                    bonds.emplace_back(idx(x, y, z), idx((x + 1) % Lx, y, z));
                }
                // y-direction
                if (periodic || y + 1 < Ly) {
                    bonds.emplace_back(idx(x, y, z), idx(x, (y + 1) % Ly, z));
                }
                // z-direction
                if (periodic || z + 1 < Lz) {
                    bonds.emplace_back(idx(x, y, z), idx(x, y, (z + 1) % Lz));
                }
            }
        }
    }
    
    Lattice lat(N, bonds);
    lat.type_ = LatticeType::Cubic;
    lat.dimensions_ = {Lx, Ly, Lz};
    return lat;
}

Lattice Lattice::custom(SiteIdx n_sites, const std::vector<Bond>& bonds) {
    Lattice lat(n_sites, bonds);
    lat.type_ = LatticeType::Custom;
    return lat;
}

void Lattice::print() const {
    const char* type_names[] = {
        "Chain", "Square", "Triangular", "Honeycomb", "Kagome", "Cubic", "Custom"
    };
    
    std::cout << fmt::format("Lattice: {} with {} sites and {} bonds\n",
                             type_names[static_cast<int>(type_)],
                             n_sites_, bonds_.size());
    
    if (!dimensions_.empty()) {
        std::cout << "  Dimensions: [";
        for (size_t i = 0; i < dimensions_.size(); ++i) {
            std::cout << dimensions_[i];
            if (i + 1 < dimensions_.size()) std::cout << " x ";
        }
        std::cout << "]\n";
    }
    
    std::cout << fmt::format("  Max coordination: {}\n", max_coordination_);
    std::cout << fmt::format("  Bond types: {}\n", n_bond_types_);
}

} // namespace sse
