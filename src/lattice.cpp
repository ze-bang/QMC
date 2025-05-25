#include "../include/lattice.hpp"
#include <cmath>
#include <stdexcept>
#include <algorithm>

namespace qmc {

using Index = Lattice::Index;
using Position = Lattice::Position;

// SquareLattice implementation
SquareLattice::SquareLattice(int L, int W, bool periodic) 
    : L_(L), W_(W), periodic_(periodic) {
    
    if (L <= 0 || W <= 0) {
        throw std::invalid_argument("Lattice dimensions must be positive");
    }
    
    // Create sites
    sites_.resize(L * W);
    for (int y = 0; y < W; ++y) {
        for (int x = 0; x < L; ++x) {
            int idx = siteIndex(x, y);
            sites_[idx].index = idx;
            sites_[idx].position = Position(x, y, 0);
            
            // Add neighbors
            std::vector<Index> neighbors;
            
            // Right neighbor
            if (x < L - 1 || periodic_) {
                neighbors.push_back(siteIndex((x + 1) % L, y));
            }
            
            // Left neighbor
            if (x > 0 || periodic_) {
                neighbors.push_back(siteIndex((x - 1 + L) % L, y));
            }
            
            // Up neighbor
            if (y < W - 1 || periodic_) {
                neighbors.push_back(siteIndex(x, (y + 1) % W));
            }
            
            // Down neighbor
            if (y > 0 || periodic_) {
                neighbors.push_back(siteIndex(x, (y - 1 + W) % W));
            }
            
            sites_[idx].neighbors = neighbors;
        }
    }
}

Index SquareLattice::siteIndex(int x, int y) const {
    return x + y * L_;
}

std::pair<int, int> SquareLattice::coordinates(Index index) const {
    return {index % L_, index / L_};
}

std::vector<std::pair<Index, Index>> SquareLattice::getBonds() const {
    std::vector<std::pair<Index, Index>> bonds;
    bonds.reserve(2 * L_ * W_);  // Approximately 2 bonds per site (right and up)
    
    for (int y = 0; y < W_; ++y) {
        for (int x = 0; x < L_; ++x) {
            Index current = siteIndex(x, y);
            
            // Add bond to the right
            if (x < L_ - 1 || periodic_) {
                bonds.emplace_back(current, siteIndex((x + 1) % L_, y));
            }
            
            // Add bond upwards
            if (y < W_ - 1 || periodic_) {
                bonds.emplace_back(current, siteIndex(x, (y + 1) % W_));
            }
        }
    }
    
    return bonds;
}

// TriangularLattice implementation
TriangularLattice::TriangularLattice(int L, int W, bool periodic) 
    : L_(L), W_(W), periodic_(periodic) {
    
    if (L <= 0 || W <= 0) {
        throw std::invalid_argument("Lattice dimensions must be positive");
    }
    
    // Create sites
    sites_.resize(L * W);
    for (int y = 0; y < W; ++y) {
        for (int x = 0; x < L; ++x) {
            int idx = siteIndex(x, y);
            sites_[idx].index = idx;
            
            // Triangular lattice has a brick-wall representation
            // with even rows shifted by 0.5
            double xpos = x + (y % 2) * 0.5;
            sites_[idx].position = Position(xpos, y, 0);
            
            // Add neighbors
            std::vector<Index> neighbors;
            
            // Right neighbor
            if (x < L - 1 || periodic_) {
                neighbors.push_back(siteIndex((x + 1) % L, y));
            }
            
            // Left neighbor
            if (x > 0 || periodic_) {
                neighbors.push_back(siteIndex((x - 1 + L) % L, y));
            }
            
            // Up-right neighbor
            if (y < W - 1 || periodic_) {
                int new_x = (y % 2 == 0) ? x : (x + 1) % L;
                if (new_x < L || periodic_) {
                    neighbors.push_back(siteIndex(new_x % L, (y + 1) % W));
                }
            }
            
            // Up-left neighbor
            if (y < W - 1 || periodic_) {
                int new_x = (y % 2 == 0) ? (x - 1 + L) % L : x;
                if (new_x >= 0 || periodic_) {
                    neighbors.push_back(siteIndex((new_x + L) % L, (y + 1) % W));
                }
            }
            
            // Down-right neighbor
            if (y > 0 || periodic_) {
                int new_x = (y % 2 == 0) ? x : (x + 1) % L;
                if (new_x < L || periodic_) {
                    neighbors.push_back(siteIndex(new_x % L, (y - 1 + W) % W));
                }
            }
            
            // Down-left neighbor
            if (y > 0 || periodic_) {
                int new_x = (y % 2 == 0) ? (x - 1 + L) % L : x;
                if (new_x >= 0 || periodic_) {
                    neighbors.push_back(siteIndex((new_x + L) % L, (y - 1 + W) % W));
                }
            }
            
            sites_[idx].neighbors = neighbors;
        }
    }
}

Index TriangularLattice::siteIndex(int x, int y) const {
    return x + y * L_;
}

std::pair<int, int> TriangularLattice::coordinates(Index index) const {
    return {index % L_, index / L_};
}

std::vector<std::pair<Index, Index>> TriangularLattice::getBonds() const {
    std::vector<std::pair<Index, Index>> bonds;
    bonds.reserve(3 * L_ * W_);  // Approximately 3 bonds per site
    
    for (int y = 0; y < W_; ++y) {
        for (int x = 0; x < L_; ++x) {
            Index current = siteIndex(x, y);
            
            // Add bond to the right
            if (x < L_ - 1 || periodic_) {
                bonds.emplace_back(current, siteIndex((x + 1) % L_, y));
            }
            
            // Add diagonal bonds
            if (y < W_ - 1 || periodic_) {
                // Up-right
                int new_x = (y % 2 == 0) ? x : (x + 1) % L_;
                if (new_x < L_ || periodic_) {
                    bonds.emplace_back(current, siteIndex(new_x % L_, (y + 1) % W_));
                }
                
                // Up-left
                new_x = (y % 2 == 0) ? (x - 1 + L_) % L_ : x;
                if ((new_x >= 0 && new_x < L_) || periodic_) {
                    bonds.emplace_back(current, siteIndex((new_x + L_) % L_, (y + 1) % W_));
                }
            }
        }
    }
    
    return bonds;
}

// CustomLattice implementation
CustomLattice::CustomLattice(
    const std::vector<Position>& positions,
    const std::vector<std::vector<Index>>& neighbors) {
    
    if (positions.size() != neighbors.size()) {
        throw std::invalid_argument("Number of positions must match number of neighbor lists");
    }
    
    // Create sites
    sites_.resize(positions.size());
    for (size_t i = 0; i < positions.size(); ++i) {
        sites_[i].index = static_cast<Index>(i);
        sites_[i].position = positions[i];
        sites_[i].neighbors = neighbors[i];
    }
    
    // Create bonds list
    bonds_.clear();
    for (size_t i = 0; i < sites_.size(); ++i) {
        for (const auto& j : sites_[i].neighbors) {
            if (i < static_cast<size_t>(j)) { // Avoid duplicates
                bonds_.emplace_back(i, j);
            }
        }
    }
}

std::vector<std::pair<Index, Index>> CustomLattice::getBonds() const {
    return bonds_;
}

} // namespace qmc
