#include "../include/sse_qmc.h"
#include "../include/lattice.h"
#include "../include/hamiltonian.h"
#include <iostream>
#include <memory>
#include <vector>
#include <iomanip>

using namespace SSE;

int main() {
    std::cout << "=== Temperature Dependence Study ===" << std::endl;
    
    // Parameters
    int Lx = 6, Ly = 6;
    double J = 1.0;
    std::vector<double> temperatures = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0};
    
    // Create lattice (shared for all temperatures)
    auto base_lattice = std::make_unique<SquareLattice>(Lx, Ly, true);
    std::cout << "System: " << base_lattice->get_name() 
              << " (" << base_lattice->size() << " sites)" << std::endl;
    std::cout << "Model: Heisenberg with J = " << J << std::endl;
    
    // Results storage
    std::cout << "\nTemperature\tBeta\tEnergy\t\tMagnetization" << std::endl;
    std::cout << std::string(60, '-') << std::endl;
    
    for (double T : temperatures) {
        double beta = 1.0 / T;
        
        // Create fresh lattice and Hamiltonian for each temperature
        auto lattice = std::make_unique<SquareLattice>(Lx, Ly, true);
        auto shared_lattice = std::shared_ptr<Lattice>(lattice.release());
        auto hamiltonian = std::make_unique<HeisenbergModel>(shared_lattice, J);
        
        // Create QMC instance
        auto qmc = std::make_unique<SSE_QMC>(
            shared_lattice,
            std::move(hamiltonian),
            beta,
            1000,
            42
        );
        
        // Run shorter simulation for temperature scan
        qmc->run_simulation(5000, 20000, 10);
        
        // Print results
        std::cout << std::fixed << std::setprecision(3)
                  << T << "\t\t"
                  << beta << "\t"
                  << std::setprecision(6)
                  << qmc->get_energy() << "\t"
                  << qmc->get_magnetization() << std::endl;
    }
    
    std::cout << "\nTemperature scan completed!" << std::endl;
    std::cout << "Note: At low temperatures, expect antiferromagnetic correlations." << std::endl;
    
    return 0;
}
