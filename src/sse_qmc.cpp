#include "sse_qmc.h"
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <climits>

namespace SSE {

bool Operator::is_diagonal() const {
    if (type == 0) return true;  // Identity is diagonal
    // This will be determined by the specific Hamiltonian term
    return type % 2 == 1;  // Simple convention: odd types are diagonal
}

SSE_QMC::SSE_QMC(std::shared_ptr<Lattice> lattice, 
                 std::unique_ptr<Hamiltonian> hamiltonian,
                 double beta, 
                 int max_operators,
                 unsigned int seed)
    : lattice_(std::move(lattice))
    , hamiltonian_(std::move(hamiltonian))
    , beta_(beta)
    , L_(max_operators)
    , n_(0)
    , rng_(seed)
    , uniform_dist_(0.0, 1.0)
    , int_dist_(0, INT_MAX)
    , total_steps_(0)
    , thermalization_steps_(0)
    , measurement_steps_(0)
    , measurement_interval_(1)
{
    // Initialize spin configuration
    spins_.resize(lattice_->size());
    initialize_spins();
    
    // Initialize operator string
    operators_.resize(L_);
    
    // Initialize vertex lists
    vertex_list_.resize(L_, std::vector<int>(4, -1));
    first_.resize(lattice_->size(), -1);
    last_.resize(lattice_->size(), -1);
    
    // Initialize observables
    observables_ = std::make_unique<Observables>(lattice_.get(), hamiltonian_.get());
    
    std::cout << "SSE QMC initialized with:" << std::endl;
    std::cout << "  Lattice: " << lattice_->get_name() << " (" << lattice_->size() << " sites)" << std::endl;
    std::cout << "  Hamiltonian: " << hamiltonian_->get_name() << std::endl;
    std::cout << "  Beta: " << beta_ << std::endl;
    std::cout << "  Max operators: " << L_ << std::endl;
}

void SSE_QMC::initialize_spins() {
    // Initialize with random spin configuration
    for (int i = 0; i < lattice_->size(); ++i) {
        spins_[i] = (random_real() < 0.5) ? 0 : 1;  // 0 = down, 1 = up
    }
}

void SSE_QMC::run_simulation(long long therm_steps, long long meas_steps, int meas_interval) {
    std::cout << "\nStarting SSE QMC simulation..." << std::endl;
    
    thermalization_steps_ = therm_steps;
    measurement_steps_ = meas_steps;
    measurement_interval_ = meas_interval;
    
    // Thermalization
    std::cout << "Thermalizing for " << therm_steps << " steps..." << std::endl;
    thermalize(therm_steps);
    
    // Measurement
    std::cout << "Measuring for " << meas_steps << " steps..." << std::endl;
    observables_->reset_all();
    measure(meas_steps, meas_interval);
    
    // Calculate final statistics
    observables_->calculate_all_statistics();
    
    std::cout << "Simulation completed!" << std::endl;
    print_statistics();
}

void SSE_QMC::thermalize(long long steps) {
    for (long long step = 0; step < steps; ++step) {
        diagonal_update();
        construct_vertex_list();
        loop_update();
        
        // Adjust cutoff periodically
        if (step % 1000 == 0) {
            adjust_cutoff();
        }
        
        if (step % 10000 == 0 && step > 0) {
            std::cout << "  Thermalization step " << step << "/" << steps 
                     << ", n = " << n_ << std::endl;
        }
    }
}

void SSE_QMC::measure(long long steps, int interval) {
    for (long long step = 0; step < steps; ++step) {
        diagonal_update();
        construct_vertex_list();
        loop_update();
        
        // Measure observables
        if (step % interval == 0) {
            observables_->measure_all(spins_, n_, beta_);
        }
        
        // Adjust cutoff periodically
        if (step % 1000 == 0) {
            adjust_cutoff();
        }
        
        if (step % 10000 == 0 && step > 0) {
            std::cout << "  Measurement step " << step << "/" << steps 
                     << ", n = " << n_ << std::endl;
        }
    }
}

void SSE_QMC::diagonal_update() {
    const int num_bonds = hamiltonian_->num_bonds();
    
    for (int p = 0; p < L_; ++p) {
        if (operators_[p].is_identity()) {
            // Try to insert operator
            int bond = random_bond();
            double prob = beta_ * hamiltonian_->get_diagonal_weight(bond, spins_) * num_bonds / (L_ - n_);
            
            if (random_real() < prob) {
                operators_[p] = Operator(1, bond);  // Insert diagonal operator
                n_++;
            }
        } else if (operators_[p].is_diagonal()) {
            // Try to remove operator
            int bond = operators_[p].bond;
            double prob = (L_ - n_ + 1) / (beta_ * hamiltonian_->get_diagonal_weight(bond, spins_) * num_bonds);
            
            if (random_real() < prob) {
                operators_[p] = Operator();  // Remove operator (set to identity)
                n_--;
            }
        } else {
            // Off-diagonal operator: propagate spins
            propagate_state(p);
        }
    }
}

void SSE_QMC::offdiagonal_update() {
    for (int p = 0; p < L_; ++p) {
        if (!operators_[p].is_identity() && operators_[p].is_diagonal()) {
            int bond = operators_[p].bond;
            
            // Try to convert to off-diagonal
            if (hamiltonian_->can_apply_offdiagonal(bond, spins_)) {
                double prob = hamiltonian_->get_offdiagonal_probability(bond);
                if (random_real() < prob) {
                    operators_[p].type = 2;  // Convert to off-diagonal
                }
            }
        }
        
        if (!operators_[p].is_identity() && !operators_[p].is_diagonal()) {
            propagate_state(p);
        }
    }
}

void SSE_QMC::construct_vertex_list() {
    // Reset vertex list
    std::fill(first_.begin(), first_.end(), -1);
    std::fill(last_.begin(), last_.end(), -1);
    
    for (int p = 0; p < L_; ++p) {
        for (int leg = 0; leg < 4; ++leg) {
            vertex_list_[p][leg] = -1;
        }
    }
    
    // Build linked vertex list
    for (int p = 0; p < L_; ++p) {
        if (!operators_[p].is_identity()) {
            int bond = operators_[p].bond;
            const auto& bond_obj = hamiltonian_->get_bond(bond);
            int site1 = bond_obj.site1;
            int site2 = bond_obj.site2;
            
            // Link vertices
            int v0 = 4 * p;      // vertex at site1, leg 0
            int v1 = 4 * p + 1;  // vertex at site1, leg 1
            int v2 = 4 * p + 2;  // vertex at site2, leg 0
            int v3 = 4 * p + 3;  // vertex at site2, leg 1
            
            if (last_[site1] != -1) {
                vertex_list_[last_[site1] / 4][last_[site1] % 4] = v0;
                vertex_list_[p][0] = last_[site1];
            } else {
                first_[site1] = v0;
            }
            
            if (last_[site2] != -1) {
                vertex_list_[last_[site2] / 4][last_[site2] % 4] = v2;
                vertex_list_[p][2] = last_[site2];
            } else {
                first_[site2] = v2;
            }
            
            last_[site1] = v1;
            last_[site2] = v3;
        }
    }
    
    // Close loops
    for (int site = 0; site < lattice_->size(); ++site) {
        if (first_[site] != -1) {
            vertex_list_[last_[site] / 4][last_[site] % 4] = first_[site];
            vertex_list_[first_[site] / 4][first_[site] % 4] = last_[site];
        }
    }
}

void SSE_QMC::loop_update() {
    std::vector<bool> visited(L_ * 4, false);
    
    for (int p = 0; p < L_; ++p) {
        if (!operators_[p].is_identity() && !operators_[p].is_diagonal()) {
            for (int leg = 0; leg < 4; leg += 2) {  // Start from legs 0 and 2
                int start_vertex = 4 * p + leg;
                if (!visited[start_vertex]) {
                    flip_spins_in_loop(start_vertex);
                    
                    // Mark all vertices in this loop as visited
                    int current = start_vertex;
                    do {
                        visited[current] = true;
                        int next_leg = current % 4;
                        next_leg = (next_leg == 0 || next_leg == 2) ? next_leg + 1 : next_leg - 1;
                        current = 4 * (current / 4) + next_leg;
                        current = vertex_list_[current / 4][current % 4];
                    } while (current != start_vertex);
                }
            }
        }
    }
}

void SSE_QMC::propagate_state(int p) {
    if (!operators_[p].is_identity()) {
        int bond = operators_[p].bond;
        if (hamiltonian_->can_apply_offdiagonal(bond, spins_)) {
            hamiltonian_->apply_offdiagonal(bond, spins_);
        }
    }
}

void SSE_QMC::flip_spins_in_loop(int start_vertex) {
    if (random_real() < 0.5) return;  // Randomly decide whether to flip
    
    int current = start_vertex;
    do {
        int p = current / 4;
        int leg = current % 4;
        int site = (leg < 2) ? hamiltonian_->get_bond(operators_[p].bond).site1 
                             : hamiltonian_->get_bond(operators_[p].bond).site2;
        
        // Flip spin
        spins_[site] = 1 - spins_[site];
        
        // Move to next vertex in loop
        int next_leg = (leg % 2 == 0) ? leg + 1 : leg - 1;
        current = 4 * p + next_leg;
        current = vertex_list_[current / 4][current % 4];
        
    } while (current != start_vertex);
}

void SSE_QMC::adjust_cutoff() {
    if (n_ > 0.7 * L_) {
        // Increase cutoff
        int new_L = std::max(L_ + 100, static_cast<int>(1.2 * L_));
        resize_operator_string();
        L_ = new_L;
    }
}

void SSE_QMC::resize_operator_string() {
    int old_L = operators_.size();
    operators_.resize(L_);
    vertex_list_.resize(L_, std::vector<int>(4, -1));
    
    // Initialize new operators as identity
    for (int i = old_L; i < L_; ++i) {
        operators_[i] = Operator();
    }
}

double SSE_QMC::get_energy() const {
    return observables_->get_observable("energy").value;
}

double SSE_QMC::get_magnetization() const {
    return observables_->get_observable("magnetization").value;
}

void SSE_QMC::print_statistics() const {
    std::cout << "\n=== SSE QMC Results ===" << std::endl;
    observables_->print_results();
    
    std::cout << "\nSimulation Statistics:" << std::endl;
    std::cout << "  Thermalization steps: " << thermalization_steps_ << std::endl;
    std::cout << "  Measurement steps: " << measurement_steps_ << std::endl;
    std::cout << "  Average operator number: " << n_ << std::endl;
    std::cout << "  Cutoff length: " << L_ << std::endl;
    std::cout << "  Fill ratio: " << static_cast<double>(n_) / L_ << std::endl;
}

} // namespace SSE
