#pragma once

#include <vector>
#include <random>
#include <memory>
#include "lattice.h"
#include "hamiltonian.h"
#include "observables.h"

namespace SSE {

struct Operator {
    int type;  // 0: identity, >0: Hamiltonian term index
    int bond;  // bond index for non-identity operators
    
    Operator(int t = 0, int b = -1) : type(t), bond(b) {}
    
    bool is_identity() const { return type == 0; }
    bool is_diagonal() const;
    bool is_offdiagonal() const { return !is_diagonal(); }
};

class SSE_QMC {
private:
    // System parameters
    std::shared_ptr<Lattice> lattice_;
    std::unique_ptr<Hamiltonian> hamiltonian_;
    std::unique_ptr<Observables> observables_;
    
    // QMC parameters
    double beta_;           // inverse temperature
    int L_;                // maximum operator string length
    int n_;                // current number of operators
    
    // State vectors
    std::vector<int> spins_;              // current spin configuration
    std::vector<Operator> operators_;     // operator string
    std::vector<std::vector<int>> vertex_list_;  // linked vertex list
    std::vector<int> first_;              // first vertex on each site
    std::vector<int> last_;               // last vertex on each site
    
    // Random number generation
    std::mt19937 rng_;
    std::uniform_real_distribution<double> uniform_dist_;
    std::uniform_int_distribution<int> int_dist_;
    
    // Statistics
    long long total_steps_;
    long long thermalization_steps_;
    long long measurement_steps_;
    int measurement_interval_;
    
public:
    SSE_QMC(std::shared_ptr<Lattice> lattice, 
            std::unique_ptr<Hamiltonian> hamiltonian,
            double beta, 
            int max_operators = 1000,
            unsigned int seed = std::random_device{}());
    
    ~SSE_QMC() = default;
    
    // Main simulation methods
    void run_simulation(long long therm_steps, long long meas_steps, int meas_interval = 1);
    void thermalize(long long steps);
    void measure(long long steps, int interval = 1);
    
    // Update methods
    void diagonal_update();
    void offdiagonal_update();
    void construct_vertex_list();
    void loop_update();
    
    // Utility methods
    void initialize_spins();
    void adjust_cutoff();
    double get_energy() const;
    double get_magnetization() const;
    void print_statistics() const;
    
    // Getters
    double get_beta() const { return beta_; }
    int get_n_operators() const { return n_; }
    const std::vector<int>& get_spins() const { return spins_; }
    const Lattice& get_lattice() const { return *lattice_; }
    const Hamiltonian& get_hamiltonian() const { return *hamiltonian_; }
    
private:
    // Helper methods
    void resize_operator_string();
    int random_site() { return int_dist_(rng_) % lattice_->size(); }
    int random_bond() { return int_dist_(rng_) % hamiltonian_->num_bonds(); }
    double random_real() { return uniform_dist_(rng_); }
    
    // Vertex operations
    void propagate_state(int p);
    void flip_spins_in_loop(int start_vertex);
};

} // namespace SSE
