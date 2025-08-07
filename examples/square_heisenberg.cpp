#include "../include/sse_qmc.h"
#include "../include/lattice.h"
#include "../include/hamiltonian.h"
#include <iostream>
#include <memory>

using namespace SSE;

int main() {
    std::cout << "=== Square Lattice Heisenberg Model Example ===" << std::endl;
    
    // Parameters
    int Lx = 8, Ly = 8;
    double J = 1.0;          // Antiferromagnetic coupling
    double beta = 10.0;      // Low temperature
    bool periodic = true;
    
    // Create square lattice
    auto lattice = std::make_unique<SquareLattice>(Lx, Ly, periodic);
    std::cout << "Created " << lattice->get_name() << " with " << lattice->size() << " sites" << std::endl;
    
    // Create Heisenberg Hamiltonian
    auto shared_lattice = std::shared_ptr<Lattice>(lattice.release());
    auto hamiltonian = std::make_unique<HeisenbergModel>(shared_lattice, J);
    std::cout << "Created " << hamiltonian->get_name() << " with J = " << J << std::endl;
    
    // Create SSE QMC instance
    auto qmc = std::make_unique<SSE_QMC>(
        shared_lattice,
        std::move(hamiltonian),
        beta,
        1000,  // Initial cutoff
        42     // Random seed
    );
    
    // Run simulation
    std::cout << "\nRunning simulation..." << std::endl;
    qmc->run_simulation(
        10000,  // Thermalization steps
        50000,  // Measurement steps
        10      // Measurement interval
    );
    
    std::cout << "\nSimulation completed successfully!" << std::endl;
    
    return 0;
}
