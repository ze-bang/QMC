#include "../include/qmc.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <iostream>
#include <random>
#include <chrono>

namespace qmc {

// StochasticSeriesExpansion implementation
StochasticSeriesExpansion::StochasticSeriesExpansion(
    std::shared_ptr<SpinModel> model, const SSEParameters& params)
    : QuantumMonteCarlo(model, params) {
    
    if (model->getSpinRepresentation() != SpinRepresentation::QUANTUM_S_HALF) {
        throw std::invalid_argument("SSE only works with quantum spin-1/2 models");
    }
    
    // Initialize with empty operator sequence
    operatorSequence_.clear();
    
    // Initialize max order
    currentOrder_ = 0;
}

void StochasticSeriesExpansion::initialize() {
    // Create a random initial spin configuration
    const auto& lattice = model_->getLattice();
    size_t n_sites = lattice->size();
    
    spinConfiguration_.resize(n_sites);
    std::uniform_int_distribution<int> dist(0, 1);
    
    for (size_t i = 0; i < n_sites; ++i) {
        // Convert 0/1 to +1/-1 (spin up/down)
        spinConfiguration_[i] = dist(rng_) == 0 ? 1 : -1;
    }
    
    // Initialize operator sequence with identity operators
    int maxOrder = static_cast<const SSEParameters&>(params_).maxOrder;
    operatorSequence_.resize(maxOrder);
    
    for (int i = 0; i < maxOrder; ++i) {
        operatorSequence_[i] = {0, -1, -1, -1}; // Identity operator
    }
    
    currentOrder_ = 0;
}

void StochasticSeriesExpansion::warmup() {
    for (int step = 0; step < params_.warmupSteps; ++step) {
        diagonalUpdate();
        offDiagonalUpdate();
        
        // Every 100 steps, adjust the maximum expansion order if needed
        if (step % 100 == 0 && static_cast<const SSEParameters&>(params_).adjustMaxOrder) {
            adjustExpansionOrder();
        }
    }
}

void StochasticSeriesExpansion::run() {
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
        diagonalUpdate();
        offDiagonalUpdate();
        
        // Every 100 steps, adjust the maximum expansion order if needed
        if (step % 100 == 0 && static_cast<const SSEParameters&>(params_).adjustMaxOrder) {
            adjustExpansionOrder();
        }
        
        // Measure energy and other observables
        calculateMeasurements();
        
        // Bin the measurements
        double step_energy = -currentOrder_ / (params_.temperature * model_->getLattice()->size());
        energyBins[current_bin] += step_energy;
        
        // Measure magnetization
        double step_mag = 0.0;
        for (size_t i = 0; i < spinConfiguration_.size(); ++i) {
            step_mag += 0.5 * spinConfiguration_[i];
        }
        step_mag /= spinConfiguration_.size();
        magBins[current_bin] += step_mag;
        
        // Measure custom observables
        for (const auto& [name, func] : observableFunctions_) {
            double val = func(spinConfiguration_);
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

void StochasticSeriesExpansion::calculateMeasurements() {
    // Energy can be directly measured from the expansion order
    // E/N = -<n>/(beta*N)
    // Other measurements require sampling of the configurations
}

void StochasticSeriesExpansion::diagonalUpdate() {
    const auto& lattice = model_->getLattice();
    size_t n_sites = lattice->size();
    auto bonds = lattice->getBonds();
    size_t n_bonds = bonds.size();
    
    double beta = 1.0 / params_.temperature;
    int max_order = static_cast<const SSEParameters&>(params_).maxOrder;
    
    // Get the quantum spin model (safe cast since we checked in constructor)
    auto qhModel = std::dynamic_pointer_cast<QuantumHeisenbergModel>(model_);
    
    // Loop through the operator sequence
    for (int p = 0; p < max_order; ++p) {
        auto& op = operatorSequence_[p];
        
        if (op.type == 0) {
            // Current operator is identity, try to insert a diagonal operator
            // with probability proportional to beta * matrix element / (M - n)
            
            // Choose a random bond
            std::uniform_int_distribution<size_t> bond_dist(0, n_bonds - 1);
            size_t bond_idx = bond_dist(rng_);
            auto bond = bonds[bond_idx];
            
            int i = bond.first;
            int j = bond.second;
            
            // Calculate the matrix element for the diagonal operator
            double diag_element = qhModel->getDiagonalMatrixElement(spinConfiguration_, i, j);
            
            // Only non-zero elements can be inserted
            if (std::abs(diag_element) > 1e-10) {
                double prob = beta * diag_element / (max_order - currentOrder_);
                
                std::uniform_real_distribution<double> dist(0.0, 1.0);
                if (dist(rng_) < prob) {
                    // Insert diagonal operator
                    op.type = 1;
                    op.bond = bond_idx;
                    op.site1 = i;
                    op.site2 = j;
                    currentOrder_++;
                }
            }
        } else if (op.type == 1) {
            // Current operator is diagonal, try to remove it
            // with probability proportional to (M - n + 1) / (beta * matrix element)
            
            int i = op.site1;
            int j = op.site2;
            
            // Calculate the matrix element
            double diag_element = qhModel->getDiagonalMatrixElement(spinConfiguration_, i, j);
            
            if (std::abs(diag_element) > 1e-10) {
                double prob = (max_order - currentOrder_ + 1) / (beta * diag_element);
                
                std::uniform_real_distribution<double> dist(0.0, 1.0);
                if (dist(rng_) < prob) {
                    // Remove diagonal operator
                    op.type = 0;
                    op.bond = -1;
                    op.site1 = -1;
                    op.site2 = -1;
                    currentOrder_--;
                }
            }
        }
        // Off-diagonal operators (type 2) are not modified in the diagonal update
    }
}

void StochasticSeriesExpansion::offDiagonalUpdate() {
    // This is the loop update algorithm for the SSE method
    // We identify closed loops in the operator-state configuration space
    // and flip all spins along these loops
    
    // For simplicity, we'll implement a single loop update
    // More advanced implementations would do multiple loops
    
    // Get the quantum spin model
    auto qhModel = std::dynamic_pointer_cast<QuantumHeisenbergModel>(model_);
    
    // First, build linked vertex list for efficient loop construction
    buildOperatorSequence();
    
    // If no operators, nothing to do
    if (currentOrder_ == 0) {
        return;
    }
    
    int max_order = static_cast<const SSEParameters&>(params_).maxOrder;
    
    // Start from a random position
    std::uniform_int_distribution<int> pos_dist(0, max_order - 1);
    int pos = pos_dist(rng_);
    
    // Skip until we find a non-identity operator
    int count = 0;
    while (operatorSequence_[pos].type == 0 && count < max_order) {
        pos = (pos + 1) % max_order;
        count++;
    }
    
    // If all operators are identity, exit
    if (count == max_order) {
        return;
    }
    
    // Start the loop update
    auto& start_op = operatorSequence_[pos];
    int site = start_op.site1; // Start with site1
    int start_site = site;
    bool loopClosed = false;
    
    std::uniform_real_distribution<double> flip_dist(0.0, 1.0);
    bool flip_spins = (flip_dist(rng_) < 0.5);
    
    // Loop until we close the loop
    while (!loopClosed) {
        auto& current_op = operatorSequence_[pos];
        
        if (current_op.type == 1) {
            // Diagonal operator - continue the path
            if (site == current_op.site1) {
                site = current_op.site2;
            } else {
                site = current_op.site1;
            }
        } else if (current_op.type == 2) {
            // Off-diagonal operator - can change the operator type
            // and continue the path
            
            // Determine which leg of the vertex we entered
            bool entered_site1 = (site == current_op.site1);
            
            // Change operator from off-diagonal to diagonal or vice versa
            // with probability depending on the matrix elements
            double offdiag_element = qhModel->getOffDiagonalMatrixElement(
                spinConfiguration_, current_op.site1, current_op.site2);
            
            double diag_element = qhModel->getDiagonalMatrixElement(
                spinConfiguration_, current_op.site1, current_op.site2);
            
            if (flip_spins) {
                // Change the operator type
                current_op.type = (current_op.type == 1) ? 2 : 1;
                
                // Flip the spins at both sites
                spinConfiguration_[current_op.site1] *= -1;
                spinConfiguration_[current_op.site2] *= -1;
            }
            
            // Continue the path through the other site
            site = entered_site1 ? current_op.site2 : current_op.site1;
        }
        
        // Move to the next position (this is simplified, in full implementation
        // we would have linked lists for efficient traversal)
        pos = (pos + 1) % max_order;
        while (operatorSequence_[pos].type == 0) {
            pos = (pos + 1) % max_order;
        }
        
        // Check if we've closed the loop
        if (site == start_site && pos == 0) {
            loopClosed = true;
        }
    }
}

void StochasticSeriesExpansion::adjustExpansionOrder() {
    int maxOrder = static_cast<const SSEParameters&>(params_).maxOrder;
    
    // If the current order is more than 80% of the maximum, double the size
    if (currentOrder_ > 0.8 * maxOrder) {
        int newMaxOrder = 2 * maxOrder;
        operatorSequence_.resize(newMaxOrder, {0, -1, -1, -1});
        
        auto& modifiedParams = const_cast<SSEParameters&>(
            static_cast<const SSEParameters&>(params_));
        modifiedParams.maxOrder = newMaxOrder;
    }
}

void StochasticSeriesExpansion::buildOperatorSequence() {
    // In a full implementation, we would build linked vertex lists
    // to efficiently identify and traverse loops
    // This is a simplified version, just identifying the operators
}

} // namespace qmc
