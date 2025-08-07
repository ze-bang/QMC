#include "../include/sse_qmc.h"
#include "../include/lattice.h"
#include "../include/hamiltonian.h"
#include <iostream>
#include <memory>
#include <vector>

using namespace SSE;

void demonstrate_lattice(std::unique_ptr<Lattice> lattice, const std::string& /* model_name */) {
    std::cout << "\n" << std::string(50, '=') << std::endl;
    lattice->print_info();
    
    // Create Heisenberg Hamiltonian for demonstration
    auto shared_lattice = std::shared_ptr<Lattice>(lattice.release());
    auto hamiltonian = std::make_unique<HeisenbergModel>(shared_lattice, 1.0);
    
    // Quick QMC run
    auto qmc = std::make_unique<SSE_QMC>(
        shared_lattice,
        std::move(hamiltonian),
        5.0,   // beta
        500,   // cutoff
        42     // seed
    );
    
    std::cout << "\nRunning short simulation..." << std::endl;
    qmc->run_simulation(2000, 5000, 10);
}

int main() {
    std::cout << "=== Different Lattice Types Demo ===" << std::endl;
    
    // 1. Square Lattice
    demonstrate_lattice(
        std::make_unique<SquareLattice>(6, 6, true),
        "Square Lattice Heisenberg"
    );
    
    // 2. Triangular Lattice
    demonstrate_lattice(
        std::make_unique<TriangularLattice>(4, 4, true),
        "Triangular Lattice Heisenberg"
    );
    
    // 3. 1D Chain
    demonstrate_lattice(
        std::make_unique<Chain>(24, true),
        "1D Chain Heisenberg"
    );
    
    // 4. Honeycomb Lattice
    demonstrate_lattice(
        std::make_unique<HoneycombLattice>(3, 3, true),
        "Honeycomb Lattice Heisenberg"
    );
    
    std::cout << "\n" << std::string(50, '=') << std::endl;
    std::cout << "Lattice comparison completed!" << std::endl;
    std::cout << "Each lattice shows different magnetic properties:" << std::endl;
    std::cout << "- Square: 2D antiferromagnet, Néel order at low T" << std::endl;
    std::cout << "- Triangular: Frustrated, potential spin liquid behavior" << std::endl;
    std::cout << "- Chain: 1D, algebraic correlations, no long-range order" << std::endl;
    std::cout << "- Honeycomb: 2D with sublattice structure" << std::endl;
    
    return 0;
}
