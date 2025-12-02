/**
 * @file sse_simulation.cpp
 * @brief Implementation of SSE simulation with diagonal and loop updates
 */

#include "sse/sse_simulation.hpp"
#include <fmt/format.h>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <nlohmann/json.hpp>

namespace sse {

SSESimulation::SSESimulation(const Lattice& lattice,
                             const Hamiltonian& hamiltonian,
                             const SimulationParams& params)
    : lattice_(lattice),
      hamiltonian_(hamiltonian),
      vertex_data_(hamiltonian),
      params_(params),
      config_(lattice, vertex_data_),
      rng_(params.seed) {
    
    // Precompute update probabilities
    // P(insert) ∝ β * N_bonds * w_max
    // P(remove) ∝ 1/β * 1/N_bonds * 1/w_max
    Real beta = params_.beta;
    BondIdx n_bonds = lattice_.numBonds();
    
    prob_insert_factor_ = beta * n_bonds;
    prob_remove_factor_ = 1.0 / (beta * n_bonds);
}

void SSESimulation::initialize() {
    config_.initialize(rng_);
    
    // Initial operator string length
    int init_M = params_.max_expansion_order > 0 
                 ? params_.max_expansion_order 
                 : std::max(10, static_cast<int>(lattice_.numSites() * params_.beta));
    config_.resize(init_M);
    
    measurements_ = std::make_unique<Measurements>(config_, params_.beta, params_.n_bins);
    
    // Reset statistics
    sweep_count_ = 0;
    n_diagonal_attempts_ = 0;
    n_diagonal_accepts_ = 0;
    n_loop_updates_ = 0;
    total_loop_length_ = 0;
}

void SSESimulation::thermalize() {
    auto start = std::chrono::high_resolution_clock::now();
    
    for (int i = 0; i < params_.n_therm; ++i) {
        sweep();
        
        // Adjust cutoff during thermalization
        if ((i + 1) % 100 == 0) {
            adjustCutoff();
        }
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    timing_.total_time += std::chrono::duration<double>(end - start).count();
    
    // Reset statistics after thermalization
    n_diagonal_attempts_ = 0;
    n_diagonal_accepts_ = 0;
}

void SSESimulation::run() {
    auto start = std::chrono::high_resolution_clock::now();
    
    int sweeps_per_bin = params_.n_sweeps / params_.n_bins;
    
    for (int bin = 0; bin < params_.n_bins; ++bin) {
        measurements_->startBin();
        
        for (int s = 0; s < sweeps_per_bin; ++s) {
            sweep();
            
            if ((s + 1) % params_.measure_every == 0) {
                measurements_->measure(config_);
            }
            
            ++sweep_count_;
        }
        
        measurements_->endBin();
        
        // Periodic checkpoint
        if (params_.checkpoint_every > 0 && (bin + 1) % params_.checkpoint_every == 0) {
            saveCheckpoint(fmt::format("checkpoint_{}.json", bin + 1));
        }
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    timing_.total_time += std::chrono::duration<double>(end - start).count();
}

void SSESimulation::sweep() {
    // Diagonal update
    auto t1 = std::chrono::high_resolution_clock::now();
    diagonalUpdate();
    auto t2 = std::chrono::high_resolution_clock::now();
    timing_.diagonal_update_time += std::chrono::duration<double>(t2 - t1).count();
    
    // Build vertex list for loop update
    config_.buildVertexList();
    
    // Loop update
    t1 = std::chrono::high_resolution_clock::now();
    if (params_.use_loop_update) {
        loopUpdate();
    }
    t2 = std::chrono::high_resolution_clock::now();
    timing_.loop_update_time += std::chrono::duration<double>(t2 - t1).count();
}

void SSESimulation::diagonalUpdate() {
    int M = config_.operatorStringLength();
    int64_t n = config_.numOperators();
    BondIdx n_bonds = lattice_.numBonds();
    
    // Propagating spin state
    std::vector<StateIdx> spin = config_.getSpins();
    
    for (int p = 0; p < M; ++p) {
        OperatorCode op = config_.getOperator(p);
        
        if (op.isIdentity()) {
            // Try to insert diagonal operator
            ++n_diagonal_attempts_;
            
            // Choose random bond
            BondIdx bond_idx = rng_.uniformInt(n_bonds);
            const Bond& bond = lattice_.getBond(bond_idx);
            int bond_type = bond.type;
            
            // Get current spins at bond sites
            StateIdx si = spin[bond.i];
            StateIdx sj = spin[bond.j];
            
            // Get vertex data for this bond type
            const VertexData& vd = vertex_data_.get(bond_type);
            VertexState vs = vd.getDiagonalVertex(si, sj);
            
            if (!vd.isAllowed(vs)) continue;
            
            Real weight = vd.getWeight(vs);
            
            // Acceptance probability for insertion
            // P(insert) = β * N_b * w / (M - n)
            Real P_insert = prob_insert_factor_ * weight / (M - n);
            
            if (rng_.uniform01() < P_insert) {
                config_.setOperator(p, OperatorCode::diagonal(bond_idx, vs));
                config_.incrementOperators();
                ++n;
                ++n_diagonal_accepts_;
            }
        } else if (op.isDiagonal()) {
            // Try to remove diagonal operator
            ++n_diagonal_attempts_;
            
            BondIdx bond_idx = op.bond();
            const Bond& bond = lattice_.getBond(bond_idx);
            int bond_type = bond.type;
            
            VertexState vs = op.vertexState();
            const VertexData& vd = vertex_data_.get(bond_type);
            Real weight = vd.getWeight(vs);
            
            // Acceptance probability for removal
            // P(remove) = (M - n + 1) / (β * N_b * w)
            Real P_remove = (M - n + 1) * prob_remove_factor_ / weight;
            
            if (rng_.uniform01() < P_remove) {
                config_.setOperator(p, OperatorCode::identity());
                config_.decrementOperators();
                --n;
                ++n_diagonal_accepts_;
            }
        } else {
            // Off-diagonal operator: propagate spin state
            BondIdx bond_idx = op.bond();
            const Bond& bond = lattice_.getBond(bond_idx);
            VertexState vs = op.vertexState();
            
            // Update spins according to vertex output legs
            spin[bond.i] = getSpinAtLeg(vs, 2);  // Output leg for site i
            spin[bond.j] = getSpinAtLeg(vs, 3);  // Output leg for site j
        }
    }
    
    // Update final spin configuration
    config_.getSpins() = spin;
}

void SSESimulation::loopUpdate() {
    int M = config_.operatorStringLength();
    if (config_.numOperators() == 0) {
        // No operators: flip spins randomly
        for (SiteIdx s = 0; s < lattice_.numSites(); ++s) {
            if (rng_.uniform01() < 0.5) {
                config_.flipSpin(s);
            }
        }
        return;
    }
    
    // Perform multiple loop updates
    int n_loops = std::max(1, static_cast<int>(config_.numOperators() / 2));
    
    for (int loop = 0; loop < n_loops; ++loop) {
        // Choose random starting vertex
        int v0;
        int attempts = 0;
        do {
            int p = rng_.uniformInt(M);
            if (config_.getOperator(p).isIdentity()) continue;
            int leg = rng_.uniformInt(NUM_LEGS);
            v0 = config_.vertexIndex(p, leg);
            ++attempts;
        } while (config_.isUnlinked(v0) && attempts < 100);
        
        if (attempts >= 100) continue;
        
        // Traverse loop
        int loop_length = loopTraverse(v0);
        total_loop_length_ += loop_length;
        ++n_loop_updates_;
    }
    
    // Update spin configuration from operators
    propagateSpins();
}

int SSESimulation::loopTraverse(int start_vertex) {
    int v = start_vertex;
    int loop_length = 0;
    const int MAX_LOOP = 1000000;
    
    // Flip direction (whether we're moving up or down in imaginary time)
    // For directed loop: we traverse through the operator string
    
    do {
        int p = config_.vertexToOperator(v);
        int leg = config_.vertexToLeg(v);
        
        OperatorCode& op = config_.getOperators()[p];
        BondIdx bond_idx = op.bond();
        const Bond& bond = lattice_.getBond(bond_idx);
        int bond_type = bond.type;
        VertexState vs = op.vertexState();
        
        const VertexData& vd = vertex_data_.get(bond_type);
        
        // Sample exit transition
        auto [exit_leg, new_vs] = vd.sampleTransition(vs, leg, rng_.uniform01());
        
        // Update operator
        if (new_vs != vs) {
            op.setVertexState(new_vs);
        }
        
        // Move to linked vertex
        int v_exit = config_.vertexIndex(p, exit_leg);
        v = config_.getLink(v_exit);
        
        if (v < 0) {
            // Unlinked (shouldn't happen in valid config)
            break;
        }
        
        ++loop_length;
    } while (v != start_vertex && loop_length < MAX_LOOP);
    
    return loop_length;
}

void SSESimulation::propagateSpins() {
    std::vector<StateIdx>& spin = config_.getSpins();
    
    // Propagate through operator string to get consistent spin configuration
    int M = config_.operatorStringLength();
    
    for (int p = 0; p < M; ++p) {
        OperatorCode op = config_.getOperator(p);
        if (op.isIdentity()) continue;
        
        BondIdx bond_idx = op.bond();
        const Bond& bond = lattice_.getBond(bond_idx);
        VertexState vs = op.vertexState();
        
        // Set input spins from current configuration
        StateIdx si_in = spin[bond.i];
        StateIdx sj_in = spin[bond.j];
        
        // Check consistency with vertex input legs
        StateIdx vs_si_in = getSpinAtLeg(vs, 0);
        StateIdx vs_sj_in = getSpinAtLeg(vs, 1);
        
        if (si_in != vs_si_in || sj_in != vs_sj_in) {
            // Inconsistency detected - update vertex state
            VertexState new_vs = makeVertexState(si_in, sj_in,
                                                  getSpinAtLeg(vs, 2),
                                                  getSpinAtLeg(vs, 3));
            config_.getOperators()[p].setVertexState(new_vs);
            vs = new_vs;
        }
        
        // Update spins to output state
        spin[bond.i] = getSpinAtLeg(vs, 2);
        spin[bond.j] = getSpinAtLeg(vs, 3);
    }
}

void SSESimulation::adjustCutoff() {
    int64_t n = config_.numOperators();
    int M = config_.operatorStringLength();
    
    // Increase cutoff if we're close to saturation
    if (n > M * 0.9) {
        int new_M = static_cast<int>(M * 1.5) + 10;
        config_.resize(new_M);
    }
    
    // Could also decrease, but usually not necessary
}

void SSESimulation::clusterUpdate() {
    // Swendsen-Wang style cluster update
    // Not used in main SSE, but provided as alternative
    
    // Mark all sites as unvisited
    std::vector<bool> visited(lattice_.numSites(), false);
    std::vector<SiteIdx> cluster;
    
    for (SiteIdx start = 0; start < lattice_.numSites(); ++start) {
        if (visited[start]) continue;
        
        // Build cluster starting from this site
        cluster.clear();
        cluster.push_back(start);
        visited[start] = true;
        
        StateIdx start_spin = config_.getSpin(start);
        
        for (size_t i = 0; i < cluster.size(); ++i) {
            SiteIdx site = cluster[i];
            
            for (SiteIdx neighbor : lattice_.neighbors(site)) {
                if (visited[neighbor]) continue;
                
                if (config_.getSpin(neighbor) == start_spin) {
                    // Add to cluster with some probability
                    // (depends on bond strength)
                    if (rng_.uniform01() < 0.5) {  // Simplified
                        cluster.push_back(neighbor);
                        visited[neighbor] = true;
                    }
                }
            }
        }
        
        // Flip cluster with probability 1/2
        if (rng_.uniform01() < 0.5) {
            for (SiteIdx site : cluster) {
                config_.flipSpin(site);
            }
        }
    }
}

void SSESimulation::saveCheckpoint(const std::string& filename) const {
    nlohmann::json j;
    
    j["sweep_count"] = sweep_count_;
    j["n_operators"] = config_.numOperators();
    j["operator_string_length"] = config_.operatorStringLength();
    
    // Save spin configuration
    j["spins"] = config_.getSpins();
    
    // Save operator string (compact representation)
    std::vector<uint64_t> ops;
    for (const auto& op : config_.getOperators()) {
        ops.push_back(op.code());
    }
    j["operators"] = ops;
    
    // Save RNG state would require more work with PCG
    
    std::ofstream file(filename);
    file << j.dump(2);
}

void SSESimulation::loadCheckpoint(const std::string& filename) {
    std::ifstream file(filename);
    nlohmann::json j;
    file >> j;
    
    sweep_count_ = j["sweep_count"];
    
    // Restore spin configuration
    std::vector<StateIdx> spins = j["spins"];
    for (SiteIdx s = 0; s < spins.size(); ++s) {
        config_.setSpin(s, spins[s]);
    }
    
    // Restore operator string
    std::vector<uint64_t> ops = j["operators"];
    config_.resize(static_cast<int>(ops.size()));
    for (size_t p = 0; p < ops.size(); ++p) {
        // Reconstruct OperatorCode from raw value
        // This is simplified; real implementation needs proper deserialization
        if (ops[p] == 0) {
            config_.setOperator(p, OperatorCode::identity());
        } else {
            // Would need to decode the operator code properly
        }
    }
    
    config_.setNumOperators(j["n_operators"]);
}

void SSESimulation::printStatus() const {
    std::cout << fmt::format("SSE Simulation Status:\n");
    std::cout << fmt::format("  Sweeps: {}\n", sweep_count_);
    std::cout << fmt::format("  Operators: {} / {}\n", 
                             config_.numOperators(), config_.operatorStringLength());
    std::cout << fmt::format("  Diagonal acceptance rate: {:.2f}%\n", 
                             acceptanceRate() * 100);
    
    if (n_loop_updates_ > 0) {
        std::cout << fmt::format("  Average loop length: {:.2f}\n",
                                 static_cast<Real>(total_loop_length_) / n_loop_updates_);
    }
    
    std::cout << fmt::format("\nTiming:\n");
    std::cout << fmt::format("  Diagonal update: {:.2f}s\n", timing_.diagonal_update_time);
    std::cout << fmt::format("  Loop update: {:.2f}s\n", timing_.loop_update_time);
    std::cout << fmt::format("  Total: {:.2f}s\n", timing_.total_time);
}

// ParallelTempering implementation

ParallelTempering::ParallelTempering(const Lattice& lattice,
                                     const Hamiltonian& hamiltonian,
                                     const std::vector<Real>& temperatures,
                                     const SimulationParams& base_params)
    : temperatures_(temperatures), n_replicas_(static_cast<int>(temperatures.size())) {
    
    simulations_.reserve(n_replicas_);
    
    for (int r = 0; r < n_replicas_; ++r) {
        SimulationParams params = base_params;
        params.beta = 1.0 / temperatures[r];
        params.seed = base_params.seed + r * 12345;
        
        simulations_.push_back(std::make_unique<SSESimulation>(lattice, hamiltonian, params));
    }
}

void ParallelTempering::initialize() {
    for (auto& sim : simulations_) {
        sim->initialize();
    }
}

void ParallelTempering::run() {
    // Run sweeps with periodic swap attempts
    // This is a simplified implementation; real HPC would use MPI
    
    for (int i = 0; i < simulations_[0]->getParams().n_sweeps; ++i) {
        // Run one sweep on each replica
        for (auto& sim : simulations_) {
            sim->sweep();
        }
        
        // Attempt swaps between adjacent replicas
        if ((i + 1) % 10 == 0) {
            attemptSwaps();
        }
    }
}

void ParallelTempering::attemptSwaps() {
    for (int r = 0; r < n_replicas_ - 1; ++r) {
        ++n_swap_attempts_;
        
        Real beta1 = simulations_[r]->getParams().beta;
        Real beta2 = simulations_[r + 1]->getParams().beta;
        
        int64_t n1 = simulations_[r]->getConfig().numOperators();
        int64_t n2 = simulations_[r + 1]->getConfig().numOperators();
        
        // Metropolis criterion for swap
        Real delta = (beta2 - beta1) * (n2 - n1);
        
        if (delta <= 0 || std::exp(-delta) > 
            static_cast<Real>(rand()) / RAND_MAX) {
            // Swap configurations
            std::swap(simulations_[r], simulations_[r + 1]);
            ++n_swap_accepts_;
        }
    }
}

} // namespace sse
