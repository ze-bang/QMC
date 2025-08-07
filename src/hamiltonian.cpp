#include "hamiltonian.h"
#include <iostream>
#include <cmath>
#include <stdexcept>

namespace SSE {

// Base Hamiltonian implementation
Hamiltonian::Hamiltonian(std::shared_ptr<Lattice> lattice) 
    : lattice_(lattice), total_diagonal_weight_(0.0) {
    bonds_ = lattice_->get_bonds();
}

void Hamiltonian::add_term(const std::vector<int>& sites, 
                          const std::vector<std::string>& operators, 
                          double coefficient) {
    if (sites.size() != operators.size()) {
        throw std::invalid_argument("Number of sites must match number of operators");
    }
    terms_.emplace_back(sites, operators, coefficient);
}

double Hamiltonian::diagonal_matrix_element(int bond, const std::vector<int>& spins) const {
    const auto& bond_obj = bonds_[bond];
    int site1 = bond_obj.site1;
    int site2 = bond_obj.site2;
    
    // For most models, find the corresponding term
    for (const auto& term : terms_) {
        if (term.sites.size() == 2 && 
            ((term.sites[0] == site1 && term.sites[1] == site2) ||
             (term.sites[0] == site2 && term.sites[1] == site1))) {
            if (term.is_diagonal) {
                return evaluate_term_diagonal(term, spins);
            }
        }
    }
    return 0.0;
}

double Hamiltonian::offdiagonal_matrix_element(int bond, const std::vector<int>& /* spins */) const {
    const auto& bond_obj = bonds_[bond];
    // For most models, the off-diagonal matrix element is the coupling strength
    return std::abs(bond_obj.coupling);
}

bool Hamiltonian::can_apply_offdiagonal(int bond, const std::vector<int>& spins) const {
    const auto& bond_obj = bonds_[bond];
    int site1 = bond_obj.site1;
    int site2 = bond_obj.site2;
    
    // For spin-1/2 systems, can flip if spins are opposite
    return spins[site1] != spins[site2];
}

void Hamiltonian::apply_offdiagonal(int bond, std::vector<int>& spins) const {
    if (!can_apply_offdiagonal(bond, spins)) return;
    
    const auto& bond_obj = bonds_[bond];
    int site1 = bond_obj.site1;
    int site2 = bond_obj.site2;
    
    // Flip both spins
    spins[site1] = 1 - spins[site1];
    spins[site2] = 1 - spins[site2];
}

double Hamiltonian::get_diagonal_weight(int bond, const std::vector<int>& spins) const {
    return std::abs(diagonal_matrix_element(bond, spins));
}

double Hamiltonian::get_offdiagonal_probability(int /* bond */) const {
    // Simple probability based on coupling strength
    return 0.5;  // For many models, this is a reasonable default
}

void Hamiltonian::calculate_diagonal_weights() {
    total_diagonal_weight_ = 0.0;
    
    // Calculate maximum possible diagonal weight
    for (const auto& bond : bonds_) {
        total_diagonal_weight_ += std::abs(bond.coupling);
    }
}

double Hamiltonian::evaluate_term_diagonal(const HamiltonianTerm& term, const std::vector<int>& spins) const {
    if (term.sites.size() == 2) {
        int s1 = spins[term.sites[0]];
        int s2 = spins[term.sites[1]];
        
        // Convert to spin values: 0 -> -1/2, 1 -> +1/2
        double sz1 = (s1 == 0) ? -0.5 : 0.5;
        double sz2 = (s2 == 0) ? -0.5 : 0.5;
        
        // Evaluate different operator types
        double result = 0.0;
        for (size_t i = 0; i < term.operators.size(); ++i) {
            if (term.operators[i] == "Sz") {
                result += (i == 0 ? sz1 : sz2);
            } else if (term.operators[i] == "SzSz") {
                result = sz1 * sz2;
                break;
            }
        }
        
        return term.coefficient * result;
    }
    
    return 0.0;
}

bool Hamiltonian::can_apply_term_offdiagonal(const HamiltonianTerm& term, const std::vector<int>& spins) const {
    if (term.is_diagonal) return false;
    
    // For S+ and S- operators
    for (size_t i = 0; i < term.operators.size(); ++i) {
        if (term.operators[i] == "S+" && spins[term.sites[i]] == 1) return false;  // Can't raise up spin
        if (term.operators[i] == "S-" && spins[term.sites[i]] == 0) return false;  // Can't lower down spin
    }
    
    return true;
}

void Hamiltonian::apply_term_offdiagonal(const HamiltonianTerm& term, std::vector<int>& spins) const {
    if (!can_apply_term_offdiagonal(term, spins)) return;
    
    for (size_t i = 0; i < term.operators.size(); ++i) {
        if (term.operators[i] == "S+") {
            spins[term.sites[i]] = 1;  // Flip up
        } else if (term.operators[i] == "S-") {
            spins[term.sites[i]] = 0;  // Flip down
        } else if (term.operators[i] == "Sx") {
            spins[term.sites[i]] = 1 - spins[term.sites[i]];  // Flip
        }
    }
}

void Hamiltonian::print_info() const {
    std::cout << "Hamiltonian Information:" << std::endl;
    std::cout << "  Name: " << get_name() << std::endl;
    std::cout << "  Number of terms: " << terms_.size() << std::endl;
    std::cout << "  Number of bonds: " << bonds_.size() << std::endl;
    std::cout << "  Total diagonal weight: " << total_diagonal_weight_ << std::endl;
}

// Heisenberg Model implementation
HeisenbergModel::HeisenbergModel(std::shared_ptr<Lattice> lattice, double J) 
    : Hamiltonian(lattice), J_(J) {
    construct_hamiltonian();
}

void HeisenbergModel::construct_hamiltonian() {
    terms_.clear();
    
    // H = J Σ_⟨i,j⟩ (S^x_i S^x_j + S^y_i S^y_j + S^z_i S^z_j)
    //   = J Σ_⟨i,j⟩ (1/2(S^+_i S^-_j + S^-_i S^+_j) + S^z_i S^z_j)
    
    for (const auto& bond : bonds_) {
        // Diagonal term: J S^z_i S^z_j
        add_term({bond.site1, bond.site2}, {"Sz", "Sz"}, J_);
        
        // Off-diagonal terms: J/2 (S^+_i S^-_j + S^-_i S^+_j)
        add_term({bond.site1, bond.site2}, {"S+", "S-"}, J_ / 2.0);
        add_term({bond.site1, bond.site2}, {"S-", "S+"}, J_ / 2.0);
    }
    
    calculate_diagonal_weights();
}

// Ising Model implementation
IsingModel::IsingModel(std::shared_ptr<Lattice> lattice, double J) 
    : Hamiltonian(lattice), J_(J) {
    construct_hamiltonian();
}

void IsingModel::construct_hamiltonian() {
    terms_.clear();
    
    // H = J Σ_⟨i,j⟩ S^z_i S^z_j
    for (const auto& bond : bonds_) {
        add_term({bond.site1, bond.site2}, {"Sz", "Sz"}, J_);
    }
    
    calculate_diagonal_weights();
}

// XY Model implementation
XYModel::XYModel(std::shared_ptr<Lattice> lattice, double J) 
    : Hamiltonian(lattice), J_(J) {
    construct_hamiltonian();
}

void XYModel::construct_hamiltonian() {
    terms_.clear();
    
    // H = J Σ_⟨i,j⟩ (S^x_i S^x_j + S^y_i S^y_j)
    //   = J Σ_⟨i,j⟩ 1/2(S^+_i S^-_j + S^-_i S^+_j)
    
    for (const auto& bond : bonds_) {
        add_term({bond.site1, bond.site2}, {"S+", "S-"}, J_ / 2.0);
        add_term({bond.site1, bond.site2}, {"S-", "S+"}, J_ / 2.0);
    }
    
    calculate_diagonal_weights();
}

// XXZ Model implementation
XXZModel::XXZModel(std::shared_ptr<Lattice> lattice, double J, double Delta) 
    : Hamiltonian(lattice), J_(J), Delta_(Delta) {
    construct_hamiltonian();
}

void XXZModel::construct_hamiltonian() {
    terms_.clear();
    
    // H = J Σ_⟨i,j⟩ (S^x_i S^x_j + S^y_i S^y_j + Δ S^z_i S^z_j)
    
    for (const auto& bond : bonds_) {
        // Diagonal term: J Δ S^z_i S^z_j
        add_term({bond.site1, bond.site2}, {"Sz", "Sz"}, J_ * Delta_);
        
        // Off-diagonal terms: J/2 (S^+_i S^-_j + S^-_i S^+_j)
        add_term({bond.site1, bond.site2}, {"S+", "S-"}, J_ / 2.0);
        add_term({bond.site1, bond.site2}, {"S-", "S+"}, J_ / 2.0);
    }
    
    calculate_diagonal_weights();
}

// Heisenberg with Field implementation
HeisenbergField::HeisenbergField(std::shared_ptr<Lattice> lattice, double J, double h) 
    : Hamiltonian(lattice), J_(J), h_(h) {
    construct_hamiltonian();
}

void HeisenbergField::construct_hamiltonian() {
    terms_.clear();
    
    // H = J Σ_⟨i,j⟩ S_i·S_j - h Σ_i S^z_i
    
    // Exchange terms
    for (const auto& bond : bonds_) {
        add_term({bond.site1, bond.site2}, {"Sz", "Sz"}, J_);
        add_term({bond.site1, bond.site2}, {"S+", "S-"}, J_ / 2.0);
        add_term({bond.site1, bond.site2}, {"S-", "S+"}, J_ / 2.0);
    }
    
    // Magnetic field terms
    for (int i = 0; i < lattice_->size(); ++i) {
        add_term({i}, {"Sz"}, -h_);
    }
    
    calculate_diagonal_weights();
}

} // namespace SSE
