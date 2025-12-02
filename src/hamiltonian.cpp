/**
 * @file hamiltonian.cpp
 * @brief Implementation of Hamiltonian class
 */

#include "sse/hamiltonian.hpp"
#include <fmt/format.h>
#include <iostream>
#include <algorithm>
#include <cmath>

namespace sse {

Hamiltonian::Hamiltonian(const BondMatrix& bond_matrix) {
    bond_matrices_.push_back(bond_matrix);
    computeEnergyOffsets();
}

Hamiltonian::Hamiltonian(const std::vector<BondMatrix>& bond_matrices)
    : bond_matrices_(bond_matrices) {
    computeEnergyOffsets();
}

void Hamiltonian::computeEnergyOffsets() {
    energy_offsets_.resize(bond_matrices_.size());
    
    for (size_t t = 0; t < bond_matrices_.size(); ++t) {
        // Find minimum diagonal element
        // We need all weights to be non-negative: w = -H + offset >= 0
        // So offset >= max(H_diagonal)
        Real max_diag = bond_matrices_[t](0, 0);
        for (int i = 1; i < 4; ++i) {
            max_diag = std::max(max_diag, bond_matrices_[t](i, i));
        }
        // Add small positive value to ensure strict positivity
        energy_offsets_[t] = max_diag + 0.25;
    }
}

Hamiltonian Hamiltonian::heisenberg(Real J) {
    // H = J * (Sx·Sx + Sy·Sy + Sz·Sz)
    // = J * (1/2 * (S+·S- + S-·S+) + Sz·Sz)
    // = J/2 * (S+·S- + S-·S+) + J * Sz·Sz
    //
    // In basis |00⟩, |01⟩, |10⟩, |11⟩:
    // Sz·Sz gives 1/4 on |00⟩ and |11⟩, -1/4 on |01⟩ and |10⟩
    // S+·S- flips |01⟩ ↔ |10⟩
    
    BondMatrix H = BondMatrix::Zero();
    
    // Diagonal: Sz_i * Sz_j
    // |00⟩: (+1/2)(+1/2) = 1/4  [both down, but convention: 0=up here? Let's use 0=down]
    // Actually, let's be consistent: |0⟩ = |↓⟩, |1⟩ = |↑⟩
    // Sz|0⟩ = -1/2|0⟩, Sz|1⟩ = +1/2|1⟩
    // 
    // |00⟩ -> index 0: Sz_i*Sz_j = (-1/2)(-1/2) = 1/4
    // |01⟩ -> index 1: Sz_i*Sz_j = (-1/2)(+1/2) = -1/4
    // |10⟩ -> index 2: Sz_i*Sz_j = (+1/2)(-1/2) = -1/4  
    // |11⟩ -> index 3: Sz_i*Sz_j = (+1/2)(+1/2) = 1/4
    
    H(0, 0) = J * 0.25;   // |↓↓⟩
    H(1, 1) = -J * 0.25;  // |↓↑⟩
    H(2, 2) = -J * 0.25;  // |↑↓⟩
    H(3, 3) = J * 0.25;   // |↑↑⟩
    
    // Off-diagonal: (S+_i S-_j + S-_i S+_j)/2
    // S+|0⟩ = |1⟩, S+|1⟩ = 0
    // S-|0⟩ = 0, S-|1⟩ = |0⟩
    // 
    // S+_i S-_j |01⟩ = S+_i|0⟩ ⊗ S-_j|1⟩ = |1⟩ ⊗ |0⟩ = |10⟩
    // S-_i S+_j |10⟩ = S-_i|1⟩ ⊗ S+_j|0⟩ = |0⟩ ⊗ |1⟩ = |01⟩
    
    H(1, 2) = J * 0.5;  // |↓↑⟩ ↔ |↑↓⟩
    H(2, 1) = J * 0.5;
    
    return Hamiltonian(H);
}

Hamiltonian Hamiltonian::xxz(Real Jxy, Real Jz) {
    BondMatrix H = BondMatrix::Zero();
    
    // Diagonal: Jz * Sz·Sz
    H(0, 0) = Jz * 0.25;
    H(1, 1) = -Jz * 0.25;
    H(2, 2) = -Jz * 0.25;
    H(3, 3) = Jz * 0.25;
    
    // Off-diagonal: Jxy/2 * (S+S- + S-S+)
    H(1, 2) = Jxy * 0.5;
    H(2, 1) = Jxy * 0.5;
    
    return Hamiltonian(H);
}

Hamiltonian Hamiltonian::xy(Real J) {
    return xxz(J, 0.0);
}

Hamiltonian Hamiltonian::ising(Real J) {
    return xxz(0.0, J);
}

Hamiltonian Hamiltonian::xxzWithField(Real Jxy, Real Jz, Real h) {
    BondMatrix H = BondMatrix::Zero();
    
    // Heisenberg part
    H(0, 0) = Jz * 0.25;
    H(1, 1) = -Jz * 0.25;
    H(2, 2) = -Jz * 0.25;
    H(3, 3) = Jz * 0.25;
    
    H(1, 2) = Jxy * 0.5;
    H(2, 1) = Jxy * 0.5;
    
    // Field: -h*(Sz_i + Sz_j) (distributed evenly to bonds)
    // For bond, add -h/z * (Sz_i + Sz_j) where z is coordination
    // Here we add full field, user should scale if needed
    // |00⟩: -h*(-1/2 - 1/2) = h
    // |01⟩: -h*(-1/2 + 1/2) = 0
    // |10⟩: -h*(+1/2 - 1/2) = 0
    // |11⟩: -h*(+1/2 + 1/2) = -h
    H(0, 0) += h;
    H(3, 3) -= h;
    
    return Hamiltonian(H);
}

Hamiltonian Hamiltonian::xyz(Real Jx, Real Jy, Real Jz) {
    BondMatrix H = BondMatrix::Zero();
    
    // Sx·Sx = 1/4 * (S+ + S-)(S+ + S-) = 1/4 * (S+S+ + S+S- + S-S+ + S-S-)
    // Sy·Sy = -1/4 * (S+ - S-)(S+ - S-) = -1/4 * (S+S+ - S+S- - S-S+ + S-S-)
    //
    // Sx·Sx + Sy·Sy = 1/2 * (S+S- + S-S+)  [XY part]
    // Sx·Sx - Sy·Sy = 1/2 * (S+S+ + S-S-)  [creates pairs]
    //
    // For Sx·Sx:
    // |00⟩↔|11⟩: 1/4
    // |01⟩↔|10⟩: 1/4
    //
    // For real SSE, we need real Hamiltonian. Standard XYZ:
    // H = Jx*Sx·Sx + Jy*Sy·Sy + Jz*Sz·Sz
    //   = (Jx+Jy)/4 * (S+S- + S-S+) + (Jx-Jy)/4 * (S+S+ + S-S-) + Jz*Sz·Sz
    
    // This creates sign problem unless Jx == Jy
    // We proceed with the general case
    
    Real Jpm = (Jx + Jy) / 4.0;  // S+S- + S-S+ coefficient
    Real Jpp = (Jx - Jy) / 4.0;  // S+S+ + S-S- coefficient
    
    // Diagonal: Jz * Sz·Sz
    H(0, 0) = Jz * 0.25;
    H(1, 1) = -Jz * 0.25;
    H(2, 2) = -Jz * 0.25;
    H(3, 3) = Jz * 0.25;
    
    // Off-diagonal from S+S- + S-S+
    H(1, 2) = Jpm;
    H(2, 1) = Jpm;
    
    // Off-diagonal from S+S+ + S-S- (pair creation/annihilation)
    H(0, 3) = Jpp;
    H(3, 0) = Jpp;
    
    return Hamiltonian(H);
}

Hamiltonian Hamiltonian::custom(const BondMatrix& matrix) {
    return Hamiltonian(matrix);
}

Hamiltonian Hamiltonian::heisenbergDM(Real J, Real Dx, Real Dy, Real Dz) {
    // H = J*S·S + D·(S×S)
    // DM term: D·(Si × Sj) = Dx(Sy_i Sz_j - Sz_i Sy_j) + cyclic
    // This is imaginary/antisymmetric, causes sign problem
    // We implement real part only
    
    BondMatrix H = heisenberg(J).getBondMatrix(0);
    
    // For real DM interaction with D along z:
    // Dz*(Sx_i Sy_j - Sy_i Sx_j) = i*Dz/2 * (S+_i S-_j - S-_i S+_j)
    // This is purely imaginary, so in real representation it's zero
    // 
    // DM interaction generally creates sign problem.
    // For now, return just Heisenberg as placeholder.
    (void)Dx; (void)Dy; (void)Dz;
    
    return Hamiltonian(H);
}

const BondMatrix& Hamiltonian::getBondMatrix(int bond_type) const {
    if (bond_type < 0 || bond_type >= static_cast<int>(bond_matrices_.size())) {
        throw std::runtime_error(fmt::format(
            "Invalid bond type {}, have {} types", bond_type, bond_matrices_.size()));
    }
    return bond_matrices_[bond_type];
}

Real Hamiltonian::matrixElement(VertexState in_state, VertexState out_state, int bond_type) const {
    // in_state encodes (s_i_in, s_j_in) in lower 2 bits
    // out_state encodes (s_i_out, s_j_out) in lower 2 bits
    // But vertex state has all 4 bits
    
    int in_idx = in_state & 0x3;   // (s_i_in, s_j_in)
    int out_idx = (out_state >> 2) & 0x3;  // (s_i_out, s_j_out)
    
    return bond_matrices_[bond_type](out_idx, in_idx);
}

bool Hamiltonian::isAllowedTransition(VertexState in_state, VertexState out_state, int bond_type) const {
    return std::abs(matrixElement(in_state, out_state, bond_type)) > TOLERANCE;
}

Real Hamiltonian::getEnergyOffset(int bond_type) const {
    if (bond_type < 0 || bond_type >= static_cast<int>(energy_offsets_.size())) {
        return 0.0;
    }
    return energy_offsets_[bond_type];
}

void Hamiltonian::print() const {
    std::cout << fmt::format("Hamiltonian with {} bond type(s)\n", bond_matrices_.size());
    
    for (size_t t = 0; t < bond_matrices_.size(); ++t) {
        std::cout << fmt::format("\nBond type {}:\n", t);
        for (int i = 0; i < 4; ++i) {
            std::cout << "  [";
            for (int j = 0; j < 4; ++j) {
                std::cout << fmt::format("{:8.4f}", bond_matrices_[t](i, j));
                if (j < 3) std::cout << ", ";
            }
            std::cout << "]\n";
        }
        std::cout << fmt::format("  Energy offset: {:.4f}\n", energy_offsets_[t]);
    }
}

} // namespace sse
