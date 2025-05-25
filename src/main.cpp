#include "../include/lattice.hpp"
#include "../include/spin_model.hpp"
#include "../include/qmc.hpp"
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <iomanip>
#include <chrono>

void printHeader() {
    std::cout << "=====================================================" << std::endl;
    std::cout << "       Quantum Monte Carlo - SSE Implementation       " << std::endl;
    std::cout << "=====================================================" << std::endl;
}

void printResults(const qmc::QuantumMonteCarlo& qmc) {
    std::cout << "Results:" << std::endl;
    std::cout << "Energy per site: " << qmc.getEnergy() 
              << " ± " << qmc.getEnergyError() << std::endl;
    std::cout << "Magnetization per site: " << qmc.getMagnetization() 
              << " ± " << qmc.getMagnetizationError() << std::endl;
    std::cout << "Specific Heat: " << qmc.getSpecificHeat() << std::endl;
    std::cout << "Susceptibility: " << qmc.getSusceptibility() << std::endl;
}

int main(int argc, char* argv[]) {
    printHeader();
    
    try {
        // Parse command-line arguments (simplified here)
        int L = 8; // Lattice size
        double T = 1.0; // Temperature
        std::string latticeType = "square";
        std::string modelType = "heisenberg";
        int warmupSteps = 10000;
        int measurementSteps = 50000;
        
        std::cout << "Quantum Monte Carlo simulation:" << std::endl;
        std::cout << "  Lattice: " << latticeType << " " << L << "x" << L << std::endl;
        std::cout << "  Model: " << modelType << std::endl;
        std::cout << "  Temperature: " << T << std::endl;
        std::cout << "  Warmup steps: " << warmupSteps << std::endl;
        std::cout << "  Measurement steps: " << measurementSteps << std::endl;
        std::cout << std::endl;
        
        // Create the lattice
        std::shared_ptr<qmc::Lattice> lattice;
        if (latticeType == "square") {
            lattice = std::make_shared<qmc::SquareLattice>(L, L);
        } else if (latticeType == "triangular") {
            lattice = std::make_shared<qmc::TriangularLattice>(L, L);
        } else {
            throw std::invalid_argument("Unknown lattice type: " + latticeType);
        }
        
        std::cout << "Created " << lattice->getName() << " lattice with " 
                  << lattice->size() << " sites" << std::endl;
        
        // Create the model
        std::shared_ptr<qmc::SpinModel> model;
        if (modelType == "ising") {
            model = std::make_shared<qmc::IsingModel>(lattice);
        } else if (modelType == "heisenberg") {
            model = std::make_shared<qmc::QuantumHeisenbergModel>(lattice);
        } else {
            throw std::invalid_argument("Unknown model type: " + modelType);
        }
        
        std::cout << "Created " << model->getName() << " model" << std::endl;
        
        // Create QMC simulation
        qmc::StochasticSeriesExpansion::SSEParameters params;
        params.temperature = T;
        params.warmupSteps = warmupSteps;
        params.measurementSteps = measurementSteps;
        params.maxOrder = 4 * lattice->size(); // Initial maxOrder guess
        
        qmc::StochasticSeriesExpansion simulation(model, params);
        
        std::cout << "Created " << simulation.getName() << " simulation" << std::endl;
        std::cout << std::endl;
        
        // Register custom observables if needed
        simulation.registerObservable("StaggeredMagnetization", [&lattice](const std::vector<double>& config) {
            double m = 0.0;
            for (size_t i = 0; i < config.size(); ++i) {
                auto coords = std::dynamic_pointer_cast<qmc::SquareLattice>(lattice)->coordinates(i);
                double sign = (coords.first + coords.second) % 2 == 0 ? 1.0 : -1.0;
                m += sign * 0.5 * config[i];
            }
            return std::abs(m) / lattice->size();
        });
        
        // Run the simulation
        auto startTime = std::chrono::high_resolution_clock::now();
        
        std::cout << "Initializing simulation..." << std::endl;
        simulation.initialize();
        
        std::cout << "Performing warmup..." << std::endl;
        simulation.warmup();
        
        std::cout << "Running measurements..." << std::endl;
        simulation.run();
        
        auto endTime = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = endTime - startTime;
        
        std::cout << "Simulation completed in " << elapsed.count() << " seconds" << std::endl;
        std::cout << std::endl;
        
        // Print results
        printResults(simulation);
        
        // Print custom observables
        std::cout << "Staggered Magnetization: " << simulation.getObservable("StaggeredMagnetization")
                  << " ± " << simulation.getObservableError("StaggeredMagnetization") << std::endl;
        
        std::cout << std::endl;
        std::cout << "Final expansion order: " << simulation.getExpansionOrder() << std::endl;
        std::cout << "Max expansion order: " << simulation.getSSEParameters().maxOrder << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
