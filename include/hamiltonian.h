#pragma once

#include <vector>
#include <string>
#include <memory>
#include "lattice.h"

namespace SSE {

struct HamiltonianTerm {
    std::vector<int> sites;
    std::vector<std::string> operators;  // "Sx", "Sy", "Sz", "S+", "S-", "I"
    double coefficient;
    bool is_diagonal;
    
    HamiltonianTerm(std::vector<int> s, std::vector<std::string> ops, double coeff)
        : sites(std::move(s)), operators(std::move(ops)), coefficient(coeff) {
        // Determine if term is diagonal
        is_diagonal = true;
        for (const auto& op : operators) {
            if (op == "S+" || op == "S-" || op == "Sx") {
                is_diagonal = false;
                break;
            }
        }
    }
};

class Hamiltonian {
protected:
    std::shared_ptr<Lattice> lattice_;
    std::vector<HamiltonianTerm> terms_;
    std::vector<Bond> bonds_;  // Copy of relevant bonds from lattice
    double total_diagonal_weight_;
    
public:
    explicit Hamiltonian(std::shared_ptr<Lattice> lattice);
    virtual ~Hamiltonian() = default;
    
    // Pure virtual methods
    virtual void construct_hamiltonian() = 0;
    virtual std::string get_name() const = 0;
    
    // Common interface
    void add_term(const std::vector<int>& sites, 
                  const std::vector<std::string>& operators, 
                  double coefficient);
    
    // Matrix element calculations
    double diagonal_matrix_element(int bond, const std::vector<int>& spins) const;
    double offdiagonal_matrix_element(int bond, const std::vector<int>& spins) const;
    bool can_apply_offdiagonal(int bond, const std::vector<int>& spins) const;
    void apply_offdiagonal(int bond, std::vector<int>& spins) const;
    
    // SSE-specific methods
    double get_diagonal_weight(int bond, const std::vector<int>& spins) const;
    double get_offdiagonal_probability(int bond) const;
    
    // Getters
    int num_bonds() const { return static_cast<int>(bonds_.size()); }
    int num_terms() const { return static_cast<int>(terms_.size()); }
    const Bond& get_bond(int i) const { return bonds_[i]; }
    const HamiltonianTerm& get_term(int i) const { return terms_[i]; }
    const std::vector<Bond>& get_bonds() const { return bonds_; }
    double get_total_diagonal_weight() const { return total_diagonal_weight_; }
    
    // Utility
    void print_info() const;
    
protected:
    void calculate_diagonal_weights();
    double evaluate_term_diagonal(const HamiltonianTerm& term, const std::vector<int>& spins) const;
    bool can_apply_term_offdiagonal(const HamiltonianTerm& term, const std::vector<int>& spins) const;
    void apply_term_offdiagonal(const HamiltonianTerm& term, std::vector<int>& spins) const;
};

// Heisenberg model: H = J Σ_⟨i,j⟩ S_i·S_j
class HeisenbergModel : public Hamiltonian {
private:
    double J_;  // Exchange coupling
    
public:
    HeisenbergModel(std::shared_ptr<Lattice> lattice, double J);
    
    void construct_hamiltonian() override;
    std::string get_name() const override { return "Heisenberg Model"; }
    
    double get_J() const { return J_; }
};

// Ising model: H = J Σ_⟨i,j⟩ S^z_i S^z_j
class IsingModel : public Hamiltonian {
private:
    double J_;  // Exchange coupling
    
public:
    IsingModel(std::shared_ptr<Lattice> lattice, double J);
    
    void construct_hamiltonian() override;
    std::string get_name() const override { return "Ising Model"; }
    
    double get_J() const { return J_; }
};

// XY model: H = J Σ_⟨i,j⟩ (S^x_i S^x_j + S^y_i S^y_j)
class XYModel : public Hamiltonian {
private:
    double J_;  // Exchange coupling
    
public:
    XYModel(std::shared_ptr<Lattice> lattice, double J);
    
    void construct_hamiltonian() override;
    std::string get_name() const override { return "XY Model"; }
    
    double get_J() const { return J_; }
};

// XXZ model: H = J Σ_⟨i,j⟩ (S^x_i S^x_j + S^y_i S^y_j + Δ S^z_i S^z_j)
class XXZModel : public Hamiltonian {
private:
    double J_;     // Exchange coupling
    double Delta_; // Anisotropy parameter
    
public:
    XXZModel(std::shared_ptr<Lattice> lattice, double J, double Delta);
    
    void construct_hamiltonian() override;
    std::string get_name() const override { return "XXZ Model"; }
    
    double get_J() const { return J_; }
    double get_Delta() const { return Delta_; }
};

// Heisenberg model with magnetic field: H = J Σ_⟨i,j⟩ S_i·S_j - h Σ_i S^z_i
class HeisenbergField : public Hamiltonian {
private:
    double J_;  // Exchange coupling
    double h_;  // Magnetic field
    
public:
    HeisenbergField(std::shared_ptr<Lattice> lattice, double J, double h);
    
    void construct_hamiltonian() override;
    std::string get_name() const override { return "Heisenberg Model with Field"; }
    
    double get_J() const { return J_; }
    double get_h() const { return h_; }
};

} // namespace SSE
