#include "../include/qmc.hpp"
#include <cmath>
#include <algorithm>
#include <numeric>

namespace qmc {

const ProjectorQuantumMonteCarlo::ProjParameters ProjectorQuantumMonteCarlo::defaultProjParams = {};

ProjectorQuantumMonteCarlo::ProjectorQuantumMonteCarlo(
    std::shared_ptr<SpinModel> model, 
    const ProjParameters& params)
    : QuantumMonteCarlo(model, params) {
    // Initialize random number generator
    std::random_device rd;
    rng_ = std::mt19937(rd());
}

void ProjectorQuantumMonteCarlo::initialize() {
    int nSites = model_->getNumberOfSites();
    const auto& projParams = static_cast<const ProjParameters&>(params_);
    
    // Initialize walkers - start with a reasonable number
    int nWalkers = std::max(100, nSites * 10);
    walkers_.resize(nWalkers);
    
    // Initialize each walker with a random configuration
    std::uniform_int_distribution<> spinDist(0, 1);
    for (auto& walker : walkers_) {
        walker.resize(nSites);
        for (int i = 0; i < nSites; ++i) {
            walker[i] = spinDist(rng_) * 2 - 1; // Convert to ±1
        }
    }
}

void ProjectorQuantumMonteCarlo::warmup() {
    // Projector QMC doesn't need traditional warmup
    // Instead, we evolve in imaginary time to reach ground state
    const auto& projParams = static_cast<const ProjParameters&>(params_);
    double dt = projParams.projectionTime / projParams.timeSlices;
    
    // Project for half the total projection time as warmup
    int warmupSlices = projParams.timeSlices / 2;
    
    for (int slice = 0; slice < warmupSlices; ++slice) {
        // Apply projection operator
        std::vector<std::vector<int>> newWalkers;
        std::vector<double> weights;
        
        for (const auto& walker : walkers_) {
            // Calculate local energy
            double localE = 0.0;
            const auto& bonds = model_->getBonds();
            
            // Diagonal part (Ising terms)
            for (const auto& bond : bonds) {
                if (bond.type == SpinModel::BondType::ISING) {
                    localE += bond.strength * walker[bond.site1] * walker[bond.site2];
                }
            }
            
            // Weight for branching
            double weight = std::exp(-dt * localE);
            weights.push_back(weight);
            
            // Off-diagonal moves (Heisenberg terms)
            std::uniform_real_distribution<> uniform(0.0, 1.0);
            std::vector<int> newWalker = walker;
            
            for (const auto& bond : bonds) {
                if (bond.type == SpinModel::BondType::HEISENBERG) {
                    // Probability to flip pair of spins
                    double pFlip = 0.5 * bond.strength * dt;
                    if (walker[bond.site1] != walker[bond.site2] && 
                        uniform(rng_) < pFlip) {
                        newWalker[bond.site1] *= -1;
                        newWalker[bond.site2] *= -1;
                    }
                }
            }
            
            newWalkers.push_back(newWalker);
        }
        
        // Branching process
        double totalWeight = std::accumulate(weights.begin(), weights.end(), 0.0);
        double avgWeight = totalWeight / weights.size();
        
        walkers_.clear();
        for (size_t i = 0; i < newWalkers.size(); ++i) {
            int nCopies = static_cast<int>(weights[i] / avgWeight + uniform(rng_));
            for (int j = 0; j < nCopies; ++j) {
                walkers_.push_back(newWalkers[i]);
            }
        }
        
        // Population control - keep number of walkers reasonable
        if (walkers_.size() > 10000) {
            // Randomly sample to reduce population
            std::shuffle(walkers_.begin(), walkers_.end(), rng_);
            walkers_.resize(5000);
        } else if (walkers_.size() < 50) {
            // Duplicate walkers if population too small
            size_t originalSize = walkers_.size();
            for (size_t i = 0; i < originalSize && walkers_.size() < 100; ++i) {
                walkers_.push_back(walkers_[i]);
            }
        }
    }
}

void ProjectorQuantumMonteCarlo::run() {
    const auto& projParams = static_cast<const ProjParameters&>(params_);
    double dt = projParams.projectionTime / projParams.timeSlices;
    
    // Measurement accumulators
    std::vector<double> energyMeasurements;
    std::vector<double> magnetizationMeasurements;
    std::vector<std::unordered_map<std::string, double>> customMeasurements;
    
    for (int step = 0; step < params_.measurementSteps; ++step) {
        // Evolve one time step
        std::vector<std::vector<int>> newWalkers;
        std::vector<double> weights;
        
        double stepEnergy = 0.0;
        double stepMagnetization = 0.0;
        std::unordered_map<std::string, double> stepCustomObs;
        
        for (const auto& walker : walkers_) {
            // Calculate local energy and magnetization
            double localE = 0.0;
            const auto& bonds = model_->getBonds();
            
            for (const auto& bond : bonds) {
                if (bond.type == SpinModel::BondType::ISING) {
                    localE += bond.strength * walker[bond.site1] * walker[bond.site2];
                }
            }
            
            stepEnergy += localE;
            stepMagnetization += std::abs(std::accumulate(walker.begin(), walker.end(), 0.0));
            
            // Custom observables
            SpinModel::SpinConfiguration config;
            config.spins = walker;
            for (const auto& [name, func] : observableFunctions_) {
                stepCustomObs[name] += func(config);
            }
            
            // Weight and evolution
            double weight = std::exp(-dt * localE);
            weights.push_back(weight);
            
            // Off-diagonal moves
            std::uniform_real_distribution<> uniform(0.0, 1.0);
            std::vector<int> newWalker = walker;
            
            for (const auto& bond : bonds) {
                if (bond.type == SpinModel::BondType::HEISENBERG) {
                    double pFlip = 0.5 * bond.strength * dt;
                    if (walker[bond.site1] != walker[bond.site2] && 
                        uniform(rng_) < pFlip) {
                        newWalker[bond.site1] *= -1;
                        newWalker[bond.site2] *= -1;
                    }
                }
            }
            
            newWalkers.push_back(newWalker);
        }
        
        // Normalize measurements
        int nWalkers = walkers_.size();
        int nSites = model_->getNumberOfSites();
        energyMeasurements.push_back(stepEnergy / (nWalkers * nSites));
        magnetizationMeasurements.push_back(stepMagnetization / (nWalkers * nSites));
        
        std::unordered_map<std::string, double> normalizedCustomObs;
        for (const auto& [name, value] : stepCustomObs) {
            normalizedCustomObs[name] = value / nWalkers;
        }
        customMeasurements.push_back(normalizedCustomObs);
        
        // Branching
        double totalWeight = std::accumulate(weights.begin(), weights.end(), 0.0);
        double avgWeight = totalWeight / weights.size();
        
        walkers_.clear();
        std::uniform_real_distribution<> uniform(0.0, 1.0);
        for (size_t i = 0; i < newWalkers.size(); ++i) {
            int nCopies = static_cast<int>(weights[i] / avgWeight + uniform(rng_));
            for (int j = 0; j < nCopies; ++j) {
                walkers_.push_back(newWalkers[i]);
            }
        }
        
        // Population control
        if (walkers_.size() > 10000) {
            std::shuffle(walkers_.begin(), walkers_.end(), rng_);
            walkers_.resize(5000);
        } else if (walkers_.size() < 50) {
            size_t originalSize = walkers_.size();
            for (size_t i = 0; i < originalSize && walkers_.size() < 100; ++i) {
                walkers_.push_back(walkers_[i]);
            }
        }
    }
    
    // Calculate averages
    energy_ = std::accumulate(energyMeasurements.begin(), energyMeasurements.end(), 0.0) 
              / energyMeasurements.size();
    magnetization_ = std::accumulate(magnetizationMeasurements.begin(), 
                                   magnetizationMeasurements.end(), 0.0) 
                     / magnetizationMeasurements.size();
    
    // Calculate energy squared for specific heat
    energySquared_ = 0.0;
    for (double e : energyMeasurements) {
        energySquared_ += e * e;
    }
    energySquared_ /= energyMeasurements.size();
    
    // Calculate magnetization squared for susceptibility
    magnetizationSquared_ = 0.0;
    for (double m : magnetizationMeasurements) {
        magnetizationSquared_ += m * m;
    }
    magnetizationSquared_ /= magnetizationMeasurements.size();
    
    // Custom observables
    for (const auto& [name, func] : observableFunctions_) {
        double sum = 0.0;
        for (const auto& measurement : customMeasurements) {
            sum += measurement.at(name);
        }
        observableValues_[name] = sum / customMeasurements.size();
    }
    
    calculateErrorsAndDerivedQuantities();
}

void ProjectorQuantumMonteCarlo::calculateMeasurements() {
    // Measurements are calculated during the run() method
    // This is because projector QMC requires special handling of walkers
}

} // namespace qmc