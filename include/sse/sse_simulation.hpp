#pragma once

/**
 * @file sse_simulation.hpp
 * @brief Main SSE simulation class with updates and measurements
 * 
 * Implements the complete SSE QMC algorithm:
 * 1. Diagonal update: insert/remove diagonal operators
 * 2. Loop update: modify operator string using directed loops
 * 3. Measurements: compute physical observables
 */

#include <vector>
#include <memory>
#include <string>
#include <functional>
#include <chrono>
#include "sse/types.hpp"
#include "sse/lattice.hpp"
#include "sse/hamiltonian.hpp"
#include "sse/vertex.hpp"
#include "sse/sse_config.hpp"
#include "sse/measurements.hpp"
#include "sse/random.hpp"

namespace sse {

/**
 * @brief Main SSE simulation class
 */
class SSESimulation {
public:
    /**
     * @brief Construct simulation with lattice and Hamiltonian
     */
    SSESimulation(const Lattice& lattice, 
                  const Hamiltonian& hamiltonian,
                  const SimulationParams& params);
    
    /**
     * @brief Initialize simulation
     */
    void initialize();
    
    /**
     * @brief Run thermalization
     */
    void thermalize();
    
    /**
     * @brief Run production sweeps with measurements
     */
    void run();
    
    /**
     * @brief Perform single Monte Carlo sweep
     */
    void sweep();
    
    /**
     * @brief Diagonal update: insert/remove diagonal operators
     */
    void diagonalUpdate();
    
    /**
     * @brief Loop update: directed loop algorithm
     */
    void loopUpdate();
    
    /**
     * @brief Simple cluster update (Swendsen-Wang style)
     */
    void clusterUpdate();
    
    /**
     * @brief Adjust operator string length
     */
    void adjustCutoff();
    
    // Getters
    const SSEConfig& getConfig() const { return config_; }
    const Measurements& getMeasurements() const { return *measurements_; }
    const SimulationParams& getParams() const { return params_; }
    
    int64_t numSweeps() const { return sweep_count_; }
    Real acceptanceRate() const { 
        return n_diagonal_attempts_ > 0 ? 
               static_cast<Real>(n_diagonal_accepts_) / n_diagonal_attempts_ : 0.0;
    }
    
    /**
     * @brief Get timing statistics
     */
    struct TimingStats {
        double diagonal_update_time;
        double loop_update_time;
        double measurement_time;
        double total_time;
    };
    TimingStats getTimingStats() const { return timing_; }
    
    /**
     * @brief Save checkpoint
     */
    void saveCheckpoint(const std::string& filename) const;
    
    /**
     * @brief Load checkpoint
     */
    void loadCheckpoint(const std::string& filename);
    
    /**
     * @brief Print simulation status
     */
    void printStatus() const;
    
private:
    // System definition
    Lattice lattice_;
    Hamiltonian hamiltonian_;
    VertexDataCollection vertex_data_;
    SimulationParams params_;
    
    // Configuration
    SSEConfig config_;
    
    // Measurements
    std::unique_ptr<Measurements> measurements_;
    
    // Random number generator
    Random rng_;
    
    // Statistics
    int64_t sweep_count_ = 0;
    int64_t n_diagonal_attempts_ = 0;
    int64_t n_diagonal_accepts_ = 0;
    int64_t n_loop_updates_ = 0;
    int64_t total_loop_length_ = 0;
    
    // Timing
    TimingStats timing_ = {0, 0, 0, 0};
    
    // Helper functions
    void propagateSpins();
    int loopTraverse(int start_vertex);
    
    // Precomputed quantities
    Real prob_insert_factor_;
    Real prob_remove_factor_;
};

/**
 * @brief Parallel tempering extension
 */
class ParallelTempering {
public:
    ParallelTempering(const Lattice& lattice,
                      const Hamiltonian& hamiltonian,
                      const std::vector<Real>& temperatures,
                      const SimulationParams& base_params);
    
    void initialize();
    void run();
    void attemptSwaps();
    
    const SSESimulation& getSimulation(int replica) const { return *simulations_[replica]; }
    
private:
    std::vector<std::unique_ptr<SSESimulation>> simulations_;
    std::vector<Real> temperatures_;
    int n_replicas_;
    int64_t n_swap_attempts_ = 0;
    int64_t n_swap_accepts_ = 0;
};

} // namespace sse
