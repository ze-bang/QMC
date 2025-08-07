#include "lattice.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <numeric>

namespace SSE {

// Base Lattice implementation
void Lattice::add_site(int index, const std::vector<double>& position) {
    sites_.emplace_back(index, position);
}

void Lattice::add_bond(int site1, int site2, double coupling, const std::string& type) {
    bonds_.emplace_back(site1, site2, coupling, type);
}

void Lattice::build_neighbor_list() {
    neighbors_.clear();
    neighbors_.resize(sites_.size());
    
    for (const auto& bond : bonds_) {
        neighbors_[bond.site1].push_back(bond.site2);
        neighbors_[bond.site2].push_back(bond.site1);
    }
    
    // Remove duplicates and sort
    for (auto& neighbor_list : neighbors_) {
        std::sort(neighbor_list.begin(), neighbor_list.end());
        neighbor_list.erase(std::unique(neighbor_list.begin(), neighbor_list.end()), 
                           neighbor_list.end());
    }
}

double Lattice::get_coordination_number() const {
    if (neighbors_.empty()) return 0.0;
    
    double total = 0.0;
    for (const auto& neighbor_list : neighbors_) {
        total += neighbor_list.size();
    }
    return total / sites_.size();
}

void Lattice::print_info() const {
    std::cout << "Lattice Information:" << std::endl;
    std::cout << "  Name: " << get_name() << std::endl;
    std::cout << "  Dimension: " << dimension_ << std::endl;
    std::cout << "  Number of sites: " << sites_.size() << std::endl;
    std::cout << "  Number of bonds: " << bonds_.size() << std::endl;
    std::cout << "  Coordination number: " << get_coordination_number() << std::endl;
}

// Square Lattice implementation
SquareLattice::SquareLattice(int Lx, int Ly, bool periodic) 
    : Lattice(2), Lx_(Lx), Ly_(Ly), periodic_(periodic) {
    generate_lattice();
}

void SquareLattice::generate_lattice() {
    sites_.clear();
    bonds_.clear();
    
    // Generate sites
    for (int y = 0; y < Ly_; ++y) {
        for (int x = 0; x < Lx_; ++x) {
            int index = y * Lx_ + x;
            sites_.emplace_back(index, std::vector<double>{static_cast<double>(x), static_cast<double>(y)});
        }
    }
    
    // Generate bonds
    for (int y = 0; y < Ly_; ++y) {
        for (int x = 0; x < Lx_; ++x) {
            int site = get_site_index(x, y);
            
            // Horizontal bonds
            if (x < Lx_ - 1) {
                int neighbor = get_site_index(x + 1, y);
                bonds_.emplace_back(site, neighbor, 1.0, "horizontal");
            } else if (periodic_) {
                int neighbor = get_site_index(0, y);
                bonds_.emplace_back(site, neighbor, 1.0, "horizontal");
            }
            
            // Vertical bonds
            if (y < Ly_ - 1) {
                int neighbor = get_site_index(x, y + 1);
                bonds_.emplace_back(site, neighbor, 1.0, "vertical");
            } else if (periodic_) {
                int neighbor = get_site_index(x, 0);
                bonds_.emplace_back(site, neighbor, 1.0, "vertical");
            }
        }
    }
    
    build_neighbor_list();
}

int SquareLattice::get_site_index(int x, int y) const {
    return y * Lx_ + x;
}

std::pair<int, int> SquareLattice::get_coordinates(int site) const {
    return {site % Lx_, site / Lx_};
}

// Triangular Lattice implementation
TriangularLattice::TriangularLattice(int Lx, int Ly, bool periodic) 
    : Lattice(2), Lx_(Lx), Ly_(Ly), periodic_(periodic) {
    generate_lattice();
}

void TriangularLattice::generate_lattice() {
    sites_.clear();
    bonds_.clear();
    
    // Generate sites
    for (int y = 0; y < Ly_; ++y) {
        for (int x = 0; x < Lx_; ++x) {
            int index = y * Lx_ + x;
            double pos_x = x + 0.5 * (y % 2);
            double pos_y = y * std::sqrt(3.0) / 2.0;
            sites_.emplace_back(index, std::vector<double>{pos_x, pos_y});
        }
    }
    
    // Generate bonds
    for (int y = 0; y < Ly_; ++y) {
        for (int x = 0; x < Lx_; ++x) {
            int site = get_site_index(x, y);
            
            // Horizontal bonds
            if (x < Lx_ - 1) {
                int neighbor = get_site_index(x + 1, y);
                bonds_.emplace_back(site, neighbor, 1.0, "horizontal");
            } else if (periodic_) {
                int neighbor = get_site_index(0, y);
                bonds_.emplace_back(site, neighbor, 1.0, "horizontal");
            }
            
            // Diagonal bonds (up-right and up-left)
            if (y < Ly_ - 1) {
                if (y % 2 == 0) {  // Even rows
                    // Up-right
                    int neighbor = get_site_index(x, y + 1);
                    bonds_.emplace_back(site, neighbor, 1.0, "diagonal1");
                    
                    // Up-left
                    if (x > 0) {
                        neighbor = get_site_index(x - 1, y + 1);
                        bonds_.emplace_back(site, neighbor, 1.0, "diagonal2");
                    } else if (periodic_) {
                        neighbor = get_site_index(Lx_ - 1, y + 1);
                        bonds_.emplace_back(site, neighbor, 1.0, "diagonal2");
                    }
                } else {  // Odd rows
                    // Up-left
                    int neighbor = get_site_index(x, y + 1);
                    bonds_.emplace_back(site, neighbor, 1.0, "diagonal1");
                    
                    // Up-right
                    if (x < Lx_ - 1) {
                        neighbor = get_site_index(x + 1, y + 1);
                        bonds_.emplace_back(site, neighbor, 1.0, "diagonal2");
                    } else if (periodic_) {
                        neighbor = get_site_index(0, y + 1);
                        bonds_.emplace_back(site, neighbor, 1.0, "diagonal2");
                    }
                }
            } else if (periodic_) {
                // Handle periodic boundary in y-direction
                if (y % 2 == 0) {
                    int neighbor = get_site_index(x, 0);
                    bonds_.emplace_back(site, neighbor, 1.0, "diagonal1");
                    
                    if (x > 0) {
                        neighbor = get_site_index(x - 1, 0);
                        bonds_.emplace_back(site, neighbor, 1.0, "diagonal2");
                    } else {
                        neighbor = get_site_index(Lx_ - 1, 0);
                        bonds_.emplace_back(site, neighbor, 1.0, "diagonal2");
                    }
                } else {
                    int neighbor = get_site_index(x, 0);
                    bonds_.emplace_back(site, neighbor, 1.0, "diagonal1");
                    
                    if (x < Lx_ - 1) {
                        neighbor = get_site_index(x + 1, 0);
                        bonds_.emplace_back(site, neighbor, 1.0, "diagonal2");
                    } else {
                        neighbor = get_site_index(0, 0);
                        bonds_.emplace_back(site, neighbor, 1.0, "diagonal2");
                    }
                }
            }
        }
    }
    
    build_neighbor_list();
}

int TriangularLattice::get_site_index(int x, int y) const {
    return y * Lx_ + x;
}

std::pair<int, int> TriangularLattice::get_coordinates(int site) const {
    return {site % Lx_, site / Lx_};
}

// 1D Chain implementation
Chain::Chain(int L, bool periodic) : Lattice(1), L_(L), periodic_(periodic) {
    generate_lattice();
}

void Chain::generate_lattice() {
    sites_.clear();
    bonds_.clear();
    
    // Generate sites
    for (int i = 0; i < L_; ++i) {
        sites_.emplace_back(i, std::vector<double>{static_cast<double>(i)});
    }
    
    // Generate bonds
    for (int i = 0; i < L_ - 1; ++i) {
        bonds_.emplace_back(i, i + 1, 1.0, "nearest_neighbor");
    }
    
    if (periodic_ && L_ > 2) {
        bonds_.emplace_back(L_ - 1, 0, 1.0, "nearest_neighbor");
    }
    
    build_neighbor_list();
}

// Honeycomb Lattice implementation
HoneycombLattice::HoneycombLattice(int Lx, int Ly, bool periodic) 
    : Lattice(2), Lx_(Lx), Ly_(Ly), periodic_(periodic) {
    generate_lattice();
}

void HoneycombLattice::generate_lattice() {
    sites_.clear();
    bonds_.clear();
    
    // Generate sites (2 sites per unit cell)
    for (int y = 0; y < Ly_; ++y) {
        for (int x = 0; x < Lx_; ++x) {
            // A sublattice
            int index_A = get_site_index(x, y, 0);
            double pos_x_A = x * 1.5;
            double pos_y_A = y * std::sqrt(3.0) + (x % 2) * std::sqrt(3.0) / 2.0;
            sites_.emplace_back(index_A, std::vector<double>{pos_x_A, pos_y_A});
            
            // B sublattice
            int index_B = get_site_index(x, y, 1);
            double pos_x_B = x * 1.5 + 0.5;
            double pos_y_B = y * std::sqrt(3.0) + (x % 2) * std::sqrt(3.0) / 2.0;
            sites_.emplace_back(index_B, std::vector<double>{pos_x_B, pos_y_B});
        }
    }
    
    // Generate bonds within unit cells and between neighboring cells
    for (int y = 0; y < Ly_; ++y) {
        for (int x = 0; x < Lx_; ++x) {
            int site_A = get_site_index(x, y, 0);
            int site_B = get_site_index(x, y, 1);
            
            // Bond within unit cell
            bonds_.emplace_back(site_A, site_B, 1.0, "intracell");
            
            // Bonds to neighboring cells
            if (x % 2 == 0) {
                // Even x: connect to upper neighbors
                if (y < Ly_ - 1) {
                    int neighbor_A = get_site_index(x, y + 1, 0);
                    bonds_.emplace_back(site_B, neighbor_A, 1.0, "intercell");
                } else if (periodic_) {
                    int neighbor_A = get_site_index(x, 0, 0);
                    bonds_.emplace_back(site_B, neighbor_A, 1.0, "intercell");
                }
                
                if (x < Lx_ - 1) {
                    int neighbor_A = get_site_index(x + 1, y, 0);
                    bonds_.emplace_back(site_B, neighbor_A, 1.0, "intercell");
                } else if (periodic_) {
                    int neighbor_A = get_site_index(0, y, 0);
                    bonds_.emplace_back(site_B, neighbor_A, 1.0, "intercell");
                }
            } else {
                // Odd x: connect to lower neighbors
                if (y > 0) {
                    int neighbor_A = get_site_index(x, y - 1, 0);
                    bonds_.emplace_back(site_B, neighbor_A, 1.0, "intercell");
                } else if (periodic_) {
                    int neighbor_A = get_site_index(x, Ly_ - 1, 0);
                    bonds_.emplace_back(site_B, neighbor_A, 1.0, "intercell");
                }
                
                if (x < Lx_ - 1) {
                    int neighbor_A = get_site_index(x + 1, y, 0);
                    bonds_.emplace_back(site_B, neighbor_A, 1.0, "intercell");
                } else if (periodic_) {
                    int neighbor_A = get_site_index(0, y, 0);
                    bonds_.emplace_back(site_B, neighbor_A, 1.0, "intercell");
                }
            }
        }
    }
    
    build_neighbor_list();
}

int HoneycombLattice::get_site_index(int x, int y, int sublattice) const {
    return 2 * (y * Lx_ + x) + sublattice;
}

} // namespace SSE
