#pragma once

#include "lattice.hpp"
#include <memory>
#include <vector>
#include <string>
#include <functional>
#include <random>
#include <Eigen/Dense>

namespace qmc {

/**
 * @brief Spin representation enum
 */
enum class SpinRepresentation {
    ISING,      // S = ±1/2 (Ising)
    HEISENBERG,  // S = continuous 3D vector (Heisenberg)
    XY,         // S = continuous 2D vector (XY model)
    QUANTUM_S_HALF  // S = 1/2 quantum spin
};

/**
 * @brief Base class for spin Hamiltonian models
 */
class SpinModel {
public:
    using SpinConfiguration = std::vector<double>;  // For Ising, values are ±1
    using SpinVectorConfiguration = std::vector<Eigen::Vector3d>;  // For Heisenberg/XY

    SpinModel(std::shared_ptr<Lattice> lattice, SpinRepresentation rep);
    virtual ~SpinModel() = default;

    // Virtual methods to be implemented by derived classes
    virtual double getEnergy(const SpinConfiguration& config) const = 0;
    virtual double getLocalEnergy(const SpinConfiguration& config, Lattice::Index site) const = 0;
    
    // For vector spin models (Heisenberg, XY)
    virtual double getEnergy(const SpinVectorConfiguration& config) const {
        throw std::runtime_error("Vector energy not implemented for this model");
    }
    
    virtual double getLocalEnergy(const SpinVectorConfiguration& config, Lattice::Index site) const {
        throw std::runtime_error("Vector local energy not implemented for this model");
    }

    // Virtual methods for calculating observables
    virtual double getMagnetization(const SpinConfiguration& config) const;
    virtual Eigen::Vector3d getMagnetization(const SpinVectorConfiguration& config) const;
    
    // Virtual methods for quantum models
    virtual int getStateSpaceDimension() const;
    virtual int getHilbertSpaceDimension() const;
    
    // Access methods
    std::shared_ptr<Lattice> getLattice() const { return lattice_; }
    SpinRepresentation getSpinRepresentation() const { return spinRep_; }
    virtual std::string getName() const = 0;

protected:
    std::shared_ptr<Lattice> lattice_;
    SpinRepresentation spinRep_;
};

/**
 * @brief Ising model Hamiltonian
 * 
 * Implements the Ising model: H = -J ∑(i,j) S_i S_j - h ∑_i S_i
 * where S_i = ±1
 */
class IsingModel : public SpinModel {
public:
    IsingModel(std::shared_ptr<Lattice> lattice, double J = 1.0, double h = 0.0);
    ~IsingModel() override = default;

    double getEnergy(const SpinConfiguration& config) const override;
    double getLocalEnergy(const SpinConfiguration& config, Lattice::Index site) const override;
    
    std::string getName() const override { return "Ising"; }
    
    // Parameters
    double getJ() const { return J_; }
    double getH() const { return h_; }

private:
    double J_; // Exchange coupling
    double h_; // External field
};

/**
 * @brief Heisenberg model Hamiltonian
 * 
 * Implements the Heisenberg model: H = -J ∑(i,j) S_i · S_j - h ∑_i S^z_i
 * where S_i is a 3D unit vector
 */
class HeisenbergModel : public SpinModel {
public:
    HeisenbergModel(std::shared_ptr<Lattice> lattice, 
                   double Jx = 1.0, double Jy = 1.0, double Jz = 1.0,
                   double h = 0.0);
    ~HeisenbergModel() override = default;

    double getEnergy(const SpinVectorConfiguration& config) const override;
    double getLocalEnergy(const SpinVectorConfiguration& config, Lattice::Index site) const override;
    
    std::string getName() const override { return "Heisenberg"; }
    
    // Parameters
    double getJx() const { return Jx_; }
    double getJy() const { return Jy_; }
    double getJz() const { return Jz_; }
    double getH() const { return h_; }

private:
    double Jx_, Jy_, Jz_; // Exchange couplings
    double h_;            // External field
};

/**
 * @brief XY model Hamiltonian
 * 
 * Implements the XY model: H = -J ∑(i,j) (S^x_i S^x_j + S^y_i S^y_j) - h ∑_i S^x_i
 */
class XYModel : public SpinModel {
public:
    XYModel(std::shared_ptr<Lattice> lattice, double J = 1.0, double h = 0.0);
    ~XYModel() override = default;

    double getEnergy(const SpinVectorConfiguration& config) const override;
    double getLocalEnergy(const SpinVectorConfiguration& config, Lattice::Index site) const override;
    
    std::string getName() const override { return "XY"; }
    
    // Parameters
    double getJ() const { return J_; }
    double getH() const { return h_; }

private:
    double J_; // Exchange coupling
    double h_; // External field
};

/**
 * @brief Quantum Heisenberg model
 * 
 * Implements the quantum Heisenberg model: 
 * H = ∑(i,j) [Jx S^x_i S^x_j + Jy S^y_i S^y_j + Jz S^z_i S^z_j] - h ∑_i S^z_i
 * Where S^a are spin-1/2 operators
 */
class QuantumHeisenbergModel : public SpinModel {
public:
    QuantumHeisenbergModel(std::shared_ptr<Lattice> lattice,
                          double Jx = 1.0, double Jy = 1.0, double Jz = 1.0,
                          double h = 0.0);
    ~QuantumHeisenbergModel() override = default;

    // For quantum models, these methods have different meanings
    double getEnergy(const SpinConfiguration& config) const override;
    double getLocalEnergy(const SpinConfiguration& config, Lattice::Index site) const override;
    
    // For SSE QMC, we need matrix elements
    double getDiagonalMatrixElement(const SpinConfiguration& config, Lattice::Index i, Lattice::Index j) const;
    double getOffDiagonalMatrixElement(const SpinConfiguration& config, Lattice::Index i, Lattice::Index j) const;
    
    // Number of states per site (2 for S=1/2)
    int getStateSpaceDimension() const override { return 2; }
    
    // Total Hilbert space dimension
    int getHilbertSpaceDimension() const override;
    
    std::string getName() const override { return "Quantum Heisenberg"; }
    
    // Parameters
    double getJx() const { return Jx_; }
    double getJy() const { return Jy_; }
    double getJz() const { return Jz_; }
    double getH() const { return h_; }

private:
    double Jx_, Jy_, Jz_; // Exchange couplings
    double h_;            // External field
};

} // namespace qmc
