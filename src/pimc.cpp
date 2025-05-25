#include "../include/qmc.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <iostream>
#include <random>
#include <chrono>
#include <functional>

namespace qmc {

PathIntegralMonteCarlo::PathIntegralMonteCarlo(std::shared_ptr<SpinModel> model, 
                                               const PIMCParameters& params)
    : QuantumMonteCarlo(model, params), configuration_() {
    if (!model) {
        throw std::invalid_argument("Model cannot be null");
    }
    
    // Initialize the configuration with random spins
    initialize();
}

// Extension of the PathIntegralMonteCarlo class with a more complete implementation
void PathIntegralMonteCarlo::initialize() {
    // Initialize with a random configuration
    const auto& lattice = model_->getLattice();
    size_t n_sites = lattice->size();
    int n_slices = static_cast<const PIMCParameters&>(params_).slices;
    
    // Create 2D configuration: [site][slice]
    configuration_.resize(n_sites);
    for (size_t i = 0; i < n_sites; ++i) {
        configuration_[i].resize(n_slices);
    }
    
    // Fill with random spins
    std::uniform_int_distribution<int> dist(0, 1);
    for (size_t i = 0; i < n_sites; ++i) {
        for (int j = 0; j < n_slices; ++j) {
            // Convert 0/1 to +1/-1 (spin up/down)
            configuration_[i][j] = dist(rng_) == 0 ? 1 : -1;
        }
    }
    
    // Make sure spins are consistent across the time boundary (periodic)
    for (size_t i = 0; i < n_sites; ++i) {
        configuration_[i][n_slices - 1] = configuration_[i][0];
    }
}

void PathIntegralMonteCarlo::warmup() {
    for (int step = 0; step < params_.warmupSteps; ++step) {
        localUpdate();
        
        // Every 10 steps, do a global update
        if (step % 10 == 0) {
            globalUpdate();
        }
    }
}

void PathIntegralMonteCarlo::run() {
    // Initialize measurement variables
    energy_ = 0.0;
    energySquared_ = 0.0;
    magnetization_ = 0.0;
    magnetizationSquared_ = 0.0;
    energyError_ = 0.0;
    magnetizationError_ = 0.0;
    
    for (auto& [name, _] : observableValues_) {
        observableValues_[name] = 0.0;
        observableErrors_[name] = 0.0;
    }
    
    // Number of bins for error estimation
    int numBins = params_.measurementSteps / params_.binSize;
    if (numBins == 0) numBins = 1;
    
    std::vector<double> energyBins(numBins, 0.0);
    std::vector<double> magBins(numBins, 0.0);
    std::map<std::string, std::vector<double>> obsBins;
    
    for (const auto& [name, _] : observableValues_) {
        obsBins[name] = std::vector<double>(numBins, 0.0);
    }
    
    int current_bin = 0;
    int measurements_in_bin = 0;
    
    // Measurement loop
    for (int step = 0; step < params_.measurementSteps; ++step) {
        // Perform Monte Carlo updates
        localUpdate();
        
        // Every 10 steps, do a global update
        if (step % 10 == 0) {
            globalUpdate();
        }
        
        // Measure energy and other observables
        calculateMeasurements();
        
        // Extract slice for measurements (we use the first time slice)
        std::vector<double> slice_config(configuration_.size());
        for (size_t i = 0; i < configuration_.size(); ++i) {
            slice_config[i] = configuration_[i][0];
        }
        
        // Bin the measurements
        double step_energy = PathIntegralMonteCarlo::calculateEnergy();
        energyBins[current_bin] += step_energy;
        
        // Measure magnetization
        double step_mag = 0.0;
        for (size_t i = 0; i < configuration_.size(); ++i) {
            step_mag += 0.5 * configuration_[i][0]; // From first time slice
        }
        step_mag /= configuration_.size();
        magBins[current_bin] += step_mag;
        
        // Measure custom observables
        for (const auto& [name, func] : observableFunctions_) {
            double val = func(slice_config);
            obsBins[name][current_bin] += val;
        }
        
        // Update bin counter
        measurements_in_bin++;
        if (measurements_in_bin >= params_.binSize) {
            // Normalize the bin averages
            energyBins[current_bin] /= measurements_in_bin;
            magBins[current_bin] /= measurements_in_bin;
            
            for (auto& [name, bins] : obsBins) {
                bins[current_bin] /= measurements_in_bin;
            }
            
            // Move to next bin
            current_bin++;
            measurements_in_bin = 0;
            
            // If all bins are full, stop measuring
            if (current_bin >= numBins) break;
        }
    }
    
    // Calculate final averages and errors
    energy_ = std::accumulate(energyBins.begin(), energyBins.end(), 0.0) / numBins;
    magnetization_ = std::accumulate(magBins.begin(), magBins.end(), 0.0) / numBins;
    
    // Calculate squared averages for specific heat and susceptibility
    for (int i = 0; i < numBins; ++i) {
        energySquared_ += energyBins[i] * energyBins[i];
        magnetizationSquared_ += magBins[i] * magBins[i];
    }
    energySquared_ /= numBins;
    magnetizationSquared_ /= numBins;
    
    // Calculate errors (standard error of the mean)
    if (numBins > 1) {
        double energy_var = 0.0;
        double mag_var = 0.0;
        
        for (int i = 0; i < numBins; ++i) {
            energy_var += (energyBins[i] - energy_) * (energyBins[i] - energy_);
            mag_var += (magBins[i] - magnetization_) * (magBins[i] - magnetization_);
        }
        
        energy_var /= (numBins - 1);
        mag_var /= (numBins - 1);
        
        energyError_ = std::sqrt(energy_var / numBins);
        magnetizationError_ = std::sqrt(mag_var / numBins);
    }
    
    // Process custom observables
    for (auto& [name, func] : observableFunctions_) {
        auto& bins = obsBins[name];
        observableValues_[name] = std::accumulate(bins.begin(), bins.end(), 0.0) / numBins;
        
        if (numBins > 1) {
            double var = 0.0;
            for (int i = 0; i < numBins; ++i) {
                var += (bins[i] - observableValues_[name]) * (bins[i] - observableValues_[name]);
            }
            var /= (numBins - 1);
            observableErrors_[name] = std::sqrt(var / numBins);
        }
    }
    
    // Calculate derived quantities like specific heat and susceptibility
    calculateErrorsAndDerivedQuantities();
}

void PathIntegralMonteCarlo::calculateMeasurements() {
    // This would include correlation functions and other observables
    // Here we just calculate energy from the configuration
}

double PathIntegralMonteCarlo::calculateEnergy() {
    // Calculate the energy from the PIMC configuration
    // For a quantum spin system, we need to consider both kinetic and potential terms
    
    const auto& lattice = model_->getLattice();
    const int n_slices = static_cast<const PIMCParameters&>(params_).slices;
    const double beta = 1.0 / params_.temperature;
    const double dtau = beta / n_slices;
    
    // Potential energy (from the model's classical energy function)
    double pot_energy = 0.0;
    
    for (int slice = 0; slice < n_slices; ++slice) {
        // Extract configuration at this time slice
        std::vector<double> slice_config(configuration_.size());
        for (size_t i = 0; i < configuration_.size(); ++i) {
            slice_config[i] = configuration_[i][slice];
        }
        
        // Get energy from the model
        pot_energy += model_->getEnergy(slice_config);
    }
    pot_energy /= n_slices;
    
    // Kinetic energy (from spin changes between slices)
    // For quantum spin systems, this comes from non-commuting terms
    double kin_energy = 0.0;
    
    for (size_t site = 0; site < configuration_.size(); ++site) {
        for (int slice = 0; slice < n_slices; ++slice) {
            int next_slice = (slice + 1) % n_slices;
            if (configuration_[site][slice] != configuration_[site][next_slice]) {
                // There's a spin flip, which contributes to kinetic energy
                kin_energy += 0.5 / dtau;  // Simplified approximation
            }
        }
    }
    kin_energy /= n_slices;
    
    return pot_energy + kin_energy;
}

void PathIntegralMonteCarlo::localUpdate() {
    // Local update for PIMC - single spin flip
    const auto& lattice = model_->getLattice();
    size_t n_sites = lattice->size();
    int n_slices = static_cast<const PIMCParameters&>(params_).slices;
    double beta = 1.0 / params_.temperature;
    double dtau = beta / n_slices;
    
    // Try to flip each spin
    for (size_t site = 0; site < n_sites; ++site) {
        for (int slice = 0; slice < n_slices; ++slice) {
            // Calculate energy change for flipping this spin
            int prev_slice = (slice > 0) ? slice - 1 : n_slices - 1;
            int next_slice = (slice < n_slices - 1) ? slice + 1 : 0;
            
            // Check for kink creation/destruction in imaginary time
            double time_energy_change = 0.0;
            if (configuration_[site][prev_slice] == configuration_[site][slice]) {
                // Creating a kink backward
                time_energy_change += 0.5 / dtau;
            } else {
                // Destroying a kink backward
                time_energy_change -= 0.5 / dtau;
            }
            
            if (configuration_[site][slice] == configuration_[site][next_slice]) {
                // Creating a kink forward
                time_energy_change += 0.5 / dtau;
            } else {
                // Destroying a kink forward
                time_energy_change -= 0.5 / dtau;
            }
            
            // Calculate spatial energy change
            // Extract configuration at this time slice
            std::vector<double> old_config(configuration_.size());
            for (size_t i = 0; i < configuration_.size(); ++i) {
                old_config[i] = configuration_[i][slice];
            }
            
            double old_energy = model_->getLocalEnergy(old_config, site);
            
            // Flip spin and recalculate
            std::vector<double> new_config = old_config;
            new_config[site] *= -1;
            
            double new_energy = model_->getLocalEnergy(new_config, site);
            
            double spatial_energy_change = new_energy - old_energy;
            
            // Total energy change
            double delta_e = time_energy_change + spatial_energy_change;
            
            // Accept or reject based on Metropolis criterion
            bool accept = false;
            if (delta_e <= 0) {
                accept = true;
            } else {
                std::uniform_real_distribution<double> dist(0.0, 1.0);
                if (dist(rng_) < std::exp(-delta_e * dtau)) {
                    accept = true;
                }
            }
            
            if (accept) {
                configuration_[site][slice] *= -1;  // Flip the spin
            }
        }
    }
}

void PathIntegralMonteCarlo::globalUpdate() {
    // Global update for PIMC using cluster algorithm
    const int n_slices = static_cast<const PIMCParameters&>(params_).slices;
    const size_t n_sites = configuration_.size();
    
    // Create a 2D array to track which sites have been visited
    std::vector<std::vector<bool>> visited(n_sites, std::vector<bool>(n_slices, false));
    
    // Pick a random starting point
    std::uniform_int_distribution<size_t> site_dist(0, n_sites - 1);
    std::uniform_int_distribution<int> slice_dist(0, n_slices - 1);
    
    size_t start_site = site_dist(rng_);
    int start_slice = slice_dist(rng_);
    
    // Decide whether to flip the cluster
    std::uniform_real_distribution<double> flip_dist(0.0, 1.0);
    bool do_flip = (flip_dist(rng_) < 0.5);
    
    // Build and potentially flip a cluster
    std::vector<std::pair<size_t, int>> stack;
    stack.emplace_back(start_site, start_slice);
    visited[start_site][start_slice] = true;
    
    int original_spin = configuration_[start_site][start_slice];
    
    while (!stack.empty()) {
        auto [site, slice] = stack.back();
        stack.pop_back();
        
        // Flip the spin if needed
        if (do_flip) {
            configuration_[site][slice] *= -1;
        }
        
        // Add neighbors in spatial dimensions with some probability
        for (const auto& neighbor : model_->getLattice()->getNeighbors(site)) {
            if (!visited[neighbor][slice] && 
                configuration_[neighbor][slice] == original_spin) {
                
                double bond_strength = 1.0;  // This should depend on the model
                double prob = 1.0 - std::exp(-2.0 * bond_strength / params_.temperature);
                
                std::uniform_real_distribution<double> dist(0.0, 1.0);
                if (dist(rng_) < prob) {
                    stack.emplace_back(neighbor, slice);
                    visited[neighbor][slice] = true;
                }
            }
        }
        
        // Add neighbors in time dimension
        int next_slice = (slice + 1) % n_slices;
        if (!visited[site][next_slice] && 
            configuration_[site][next_slice] == original_spin) {
            
            stack.emplace_back(site, next_slice);
            visited[site][next_slice] = true;
        }
        
        int prev_slice = (slice > 0) ? slice - 1 : n_slices - 1;
        if (!visited[site][prev_slice] && 
            configuration_[site][prev_slice] == original_spin) {
            
            stack.emplace_back(site, prev_slice);
            visited[site][prev_slice] = true;
        }
    }
}


} // namespace qmc
