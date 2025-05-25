#pragma once

#include "spin_model.hpp"
#include <vector>
#include <random>
#include <memory>
#include <functional>
#include <string>
#include <unordered_map>

namespace qmc {

/**
 * @brief Abstract base class for all QMC algorithms
 */
class QuantumMonteCarlo {
public:
    struct Parameters {
        int warmupSteps = 10000;      // Number of warmup steps
        int measurementSteps = 50000; // Number of measurement steps
        int binSize = 100;            // Size of measurement bins for error analysis
        double temperature = 1.0;     // Temperature 
    };
    
    using ObservableFunction = std::function<double(const SpinModel::SpinConfiguration&)>;
    
    QuantumMonteCarlo(std::shared_ptr<SpinModel> model, const Parameters& params);
    virtual ~QuantumMonteCarlo() = default;
    
    // Main methods
    virtual void initialize() = 0;
    virtual void warmup() = 0;
    virtual void run() = 0;
    
    // Register a custom observable to measure
    void registerObservable(const std::string& name, ObservableFunction func);
    
    // Get results
    double getEnergy() const { return energy_; }
    double getEnergyError() const { return energyError_; }
    double getSpecificHeat() const { return specificHeat_; }
    double getMagnetization() const { return magnetization_; }
    double getMagnetizationError() const { return magnetizationError_; }
    double getSusceptibility() const { return susceptibility_; }
    
    // Get a custom observable result
    double getObservable(const std::string& name) const;
    double getObservableError(const std::string& name) const;
    
    // Set parameters
    void setTemperature(double T) { params_.temperature = T; }
    void setWarmupSteps(int steps) { params_.warmupSteps = steps; }
    void setMeasurementSteps(int steps) { params_.measurementSteps = steps; }
    void setBinSize(int size) { params_.binSize = size; }
    
    // Get model
    std::shared_ptr<SpinModel> getModel() const { return model_; }
    
    // Get name of the algorithm
    virtual std::string getName() const = 0;

protected:
    std::shared_ptr<SpinModel> model_;
    Parameters params_;
    std::mt19937 rng_;
    
    // Common measurements
    double energy_ = 0.0;
    double energySquared_ = 0.0;
    double energyError_ = 0.0;
    double magnetization_ = 0.0;
    double magnetizationSquared_ = 0.0;
    double magnetizationError_ = 0.0;
    double specificHeat_ = 0.0;
    double susceptibility_ = 0.0;
    
    // Custom observables
    std::unordered_map<std::string, ObservableFunction> observableFunctions_;
    std::unordered_map<std::string, double> observableValues_;
    std::unordered_map<std::string, double> observableErrors_;
    
    // Helper methods
    void calculateErrorsAndDerivedQuantities();
    virtual void calculateMeasurements() = 0;
};

/**
 * @brief Stochastic Series Expansion (SSE) QMC algorithm
 * 
 * Implementation of the SSE algorithm for quantum spin models.
 */
class StochasticSeriesExpansion : public QuantumMonteCarlo {
public:
    struct SSEParameters : public Parameters {
        int maxOrder = 100;  // Maximum expansion order
        bool adjustMaxOrder = true;  // Auto-adjust expansion order
    };
    
    static const SSEParameters defaultSSEParams;
    StochasticSeriesExpansion(std::shared_ptr<SpinModel> model, 
                             const SSEParameters& params = defaultSSEParams);
    ~StochasticSeriesExpansion() override = default;
    
    void initialize() override;
    void warmup() override;
    void run() override;
    
    std::string getName() const override { return "SSE"; }
    
    // SSE specific methods
    int getExpansionOrder() const { return currentOrder_; }
    void setMaxExpansionOrder(int order) { 
        static_cast<SSEParameters&>(params_).maxOrder = order; 
    }
    
    // Access to parameters
    const SSEParameters& getSSEParameters() const {
        return static_cast<const SSEParameters&>(params_);
    }

protected:
    void calculateMeasurements() override;
    
    // SSE specific moves
    void diagonalUpdate();
    void offDiagonalUpdate();
    
    // Helper methods for SSE
    void adjustExpansionOrder();
    void buildOperatorSequence();
    
    // SSE specific data structures
    struct Operator {
        int type;        // 0: identity, 1: diagonal, 2: off-diagonal
        int bond;        // Index in the bond list
        int site1, site2; // Sites involved
    };
    
    std::vector<Operator> operatorSequence_;
    std::vector<double> spinConfiguration_;
    int currentOrder_ = 0;
};

/**
 * @brief Path Integral Monte Carlo (PIMC) QMC algorithm
 * 
 * Implementation of the PIMC algorithm for quantum spin models.
 */
class PathIntegralMonteCarlo : public QuantumMonteCarlo {
public:
    struct PIMCParameters : public Parameters {
        int slices = 40;  // Number of imaginary time slices
    };
    static const PIMCParameters defaultPIMCParams;
    PathIntegralMonteCarlo(std::shared_ptr<SpinModel> model, 
                          const PIMCParameters& params = defaultPIMCParams);
    ~PathIntegralMonteCarlo() override = default;
    
    void initialize() override;
    void warmup() override;
    void run() override;
    
    std::string getName() const override { return "PIMC"; }

protected:
    void calculateMeasurements() override;
    
    // PIMC specific moves
    void localUpdate();
    void globalUpdate();
    
    // Energy calculation
    double calculateEnergy();
    
    // PIMC specific data structures
    // 2D array: [site][slice]
    std::vector<std::vector<int>> configuration_;
};

/**
 * @brief Zero-temperature Projector QMC algorithm
 * 
 * Implementation of the projector QMC algorithm for ground state properties.
 */
class ProjectorQuantumMonteCarlo : public QuantumMonteCarlo {
public:
    struct ProjParameters : public Parameters {
        double projectionTime = 20.0;  // Total projection time
        int timeSlices = 100;          // Number of time slices
    };
    static const ProjParameters defaultProjParams;
    ProjectorQuantumMonteCarlo(std::shared_ptr<SpinModel> model, 
                              const ProjParameters& params = defaultProjParams);
    ~ProjectorQuantumMonteCarlo() override = default;
    
    void initialize() override;
    void warmup() override;
    void run() override;
    
    std::string getName() const override { return "Projector QMC"; }

protected:
    void calculateMeasurements() override;
    
    // Projector QMC specific data structure
    std::vector<std::vector<int>> walkers_; // Multiple walkers
};

} // namespace qmc
