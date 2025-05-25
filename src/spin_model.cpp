#include "../include/spin_model.hpp"
#include <cmath>
#include <random>
#include <numeric>
#include <algorithm>

namespace qmc {

// Base SpinModel implementation
SpinModel::SpinModel(std::shared_ptr<Lattice> lattice, SpinRepresentation rep)
    : lattice_(lattice), spinRep_(rep) {
    if (!lattice) {
        throw std::invalid_argument("Lattice cannot be null");
    }
}

double SpinModel::getMagnetization(const SpinConfiguration& config) const {
    if (config.size() != lattice_->size()) {
        throw std::invalid_argument("Configuration size does not match lattice size");
    }
    
    double mag = std::accumulate(config.begin(), config.end(), 0.0);
    return mag / lattice_->size();
}

Eigen::Vector3d SpinModel::getMagnetization(const SpinVectorConfiguration& config) const {
    if (config.size() != lattice_->size()) {
        throw std::invalid_argument("Configuration size does not match lattice size");
    }
    
    Eigen::Vector3d mag = std::accumulate(config.begin(), config.end(), 
                                         Eigen::Vector3d::Zero().eval());
    return mag / lattice_->size();
}

int SpinModel::getStateSpaceDimension() const {
    switch (spinRep_) {
        case SpinRepresentation::ISING:
        case SpinRepresentation::QUANTUM_S_HALF:
            return 2;
        case SpinRepresentation::HEISENBERG:
        case SpinRepresentation::XY:
            return 3;  // Continuous, but we return dimension of the space
        default:
            return 1;
    }
}

int SpinModel::getHilbertSpaceDimension() const {
    return static_cast<int>(std::pow(getStateSpaceDimension(), lattice_->size()));
}

// IsingModel implementation
IsingModel::IsingModel(std::shared_ptr<Lattice> lattice, double J, double h)
    : SpinModel(lattice, SpinRepresentation::ISING), J_(J), h_(h) {}

double IsingModel::getEnergy(const SpinConfiguration& config) const {
    if (config.size() != lattice_->size()) {
        throw std::invalid_argument("Configuration size does not match lattice size");
    }
    
    double energy = 0.0;
    
    // Interaction energy
    for (const auto& bond : lattice_->getBonds()) {
        energy -= J_ * config[bond.first] * config[bond.second];
    }
    
    // External field energy
    if (h_ != 0.0) {
        for (size_t i = 0; i < config.size(); ++i) {
            energy -= h_ * config[i];
        }
    }
    
    return energy;
}

double IsingModel::getLocalEnergy(const SpinConfiguration& config, Lattice::Index site) const {
    if (config.size() != lattice_->size() || site >= static_cast<Lattice::Index>(config.size())) {
        throw std::invalid_argument("Invalid configuration or site index");
    }
    
    double energy = 0.0;
    
    // Interaction energy with neighbors
    for (const auto& neighbor : lattice_->getNeighbors(site)) {
        energy -= J_ * config[site] * config[neighbor];
    }
    
    // External field energy
    energy -= h_ * config[site];
    
    return energy;
}

// HeisenbergModel implementation
HeisenbergModel::HeisenbergModel(std::shared_ptr<Lattice> lattice, 
                               double Jx, double Jy, double Jz, double h)
    : SpinModel(lattice, SpinRepresentation::HEISENBERG), 
      Jx_(Jx), Jy_(Jy), Jz_(Jz), h_(h) {}

double HeisenbergModel::getEnergy(const SpinVectorConfiguration& config) const {
    if (config.size() != lattice_->size()) {
        throw std::invalid_argument("Configuration size does not match lattice size");
    }
    
    double energy = 0.0;
    
    // Interaction energy
    for (const auto& bond : lattice_->getBonds()) {
        const auto& S_i = config[bond.first];
        const auto& S_j = config[bond.second];
        
        energy -= Jx_ * S_i.x() * S_j.x();
        energy -= Jy_ * S_i.y() * S_j.y();
        energy -= Jz_ * S_i.z() * S_j.z();
    }
    
    // External field energy
    if (h_ != 0.0) {
        for (const auto& spin : config) {
            energy -= h_ * spin.z();
        }
    }
    
    return energy;
}

double HeisenbergModel::getLocalEnergy(const SpinVectorConfiguration& config, Lattice::Index site) const {
    if (config.size() != lattice_->size() || site >= static_cast<Lattice::Index>(config.size())) {
        throw std::invalid_argument("Invalid configuration or site index");
    }
    
    double energy = 0.0;
    const auto& spin_i = config[site];
    
    // Interaction energy with neighbors
    for (const auto& neighbor : lattice_->getNeighbors(site)) {
        const auto& spin_j = config[neighbor];
        
        energy -= Jx_ * spin_i.x() * spin_j.x();
        energy -= Jy_ * spin_i.y() * spin_j.y();
        energy -= Jz_ * spin_i.z() * spin_j.z();
    }
    
    // External field energy
    energy -= h_ * spin_i.z();
    
    return energy;
}

// XYModel implementation
XYModel::XYModel(std::shared_ptr<Lattice> lattice, double J, double h)
    : SpinModel(lattice, SpinRepresentation::XY), J_(J), h_(h) {}

double XYModel::getEnergy(const SpinVectorConfiguration& config) const {
    if (config.size() != lattice_->size()) {
        throw std::invalid_argument("Configuration size does not match lattice size");
    }
    
    double energy = 0.0;
    
    // Interaction energy
    for (const auto& bond : lattice_->getBonds()) {
        const auto& S_i = config[bond.first];
        const auto& S_j = config[bond.second];
        
        energy -= J_ * (S_i.x() * S_j.x() + S_i.y() * S_j.y());
    }
    
    // External field energy
    if (h_ != 0.0) {
        for (const auto& spin : config) {
            energy -= h_ * spin.x();
        }
    }
    
    return energy;
}

double XYModel::getLocalEnergy(const SpinVectorConfiguration& config, Lattice::Index site) const {
    if (config.size() != lattice_->size() || site >= static_cast<Lattice::Index>(config.size())) {
        throw std::invalid_argument("Invalid configuration or site index");
    }
    
    double energy = 0.0;
    const auto& spin_i = config[site];
    
    // Interaction energy with neighbors
    for (const auto& neighbor : lattice_->getNeighbors(site)) {
        const auto& spin_j = config[neighbor];
        
        energy -= J_ * (spin_i.x() * spin_j.x() + spin_i.y() * spin_j.y());
    }
    
    // External field energy
    energy -= h_ * spin_i.x();
    
    return energy;
}

// QuantumHeisenbergModel implementation
QuantumHeisenbergModel::QuantumHeisenbergModel(std::shared_ptr<Lattice> lattice,
                                             double Jx, double Jy, double Jz, double h)
    : SpinModel(lattice, SpinRepresentation::QUANTUM_S_HALF),
      Jx_(Jx), Jy_(Jy), Jz_(Jz), h_(h) {}

double QuantumHeisenbergModel::getEnergy(const SpinConfiguration& config) const {
    // This is a classical approximation of the quantum energy
    // For QMC, we need matrix elements, not direct energy
    if (config.size() != lattice_->size()) {
        throw std::invalid_argument("Configuration size does not match lattice size");
    }
    
    double energy = 0.0;
    
    // Diagonal (Sz-Sz) interaction energy
    for (const auto& bond : lattice_->getBonds()) {
        energy += 0.25 * Jz_ * config[bond.first] * config[bond.second];
    }
    
    // External field energy
    if (h_ != 0.0) {
        for (size_t i = 0; i < config.size(); ++i) {
            energy -= 0.5 * h_ * config[i];
        }
    }
    
    return energy;
}

double QuantumHeisenbergModel::getLocalEnergy(const SpinConfiguration& config, Lattice::Index site) const {
    // This is a classical approximation 
    if (config.size() != lattice_->size() || site >= static_cast<Lattice::Index>(config.size())) {
        throw std::invalid_argument("Invalid configuration or site index");
    }
    
    double energy = 0.0;
    
    // Diagonal (Sz-Sz) interaction energy with neighbors
    for (const auto& neighbor : lattice_->getNeighbors(site)) {
        energy += 0.25 * Jz_ * config[site] * config[neighbor];
    }
    
    // External field energy
    energy -= 0.5 * h_ * config[site];
    
    return energy;
}

double QuantumHeisenbergModel::getDiagonalMatrixElement(
    const SpinConfiguration& config, Lattice::Index i, Lattice::Index j) const {
    
    // Diagonal term: Jz * Sz_i * Sz_j + h * (Sz_i + Sz_j)/2
    return 0.25 * Jz_ * config[i] * config[j] - 0.5 * h_ * (config[i] + config[j]);
}

double QuantumHeisenbergModel::getOffDiagonalMatrixElement(
    const SpinConfiguration& config, Lattice::Index i, Lattice::Index j) const {
    
    // Off-diagonal term: 0.5 * (Jx * Sx_i * Sx_j + Jy * Sy_i * Sy_j)
    // For spin-1/2: <↑↓|Sx_i Sx_j + Sy_i Sy_j|↓↑> = 0.5 * (Jx + Jy)
    
    // This is only non-zero when the spins are opposite
    if (config[i] != config[j]) {
        return 0.5 * (Jx_ + Jy_);
    }
    
    return 0.0;
}

int QuantumHeisenbergModel::getHilbertSpaceDimension() const {
    return static_cast<int>(std::pow(2.0, static_cast<double>(lattice_->size())));
}

} // namespace qmc
