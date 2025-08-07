#include "../include/sse_qmc.h"
#include "../include/lattice.h"
#include "../include/hamiltonian.h"
#include <iostream>
#include <memory>
#include <vector>
#include <iomanip>

using namespace SSE;

void test_model(std::unique_ptr<Hamiltonian> hamiltonian, const std::string& name, double beta) {
    std::cout << "\n" << std::string(40, '-') << std::endl;
    std::cout << "Testing: " << name << std::endl;
    hamiltonian->print_info();
    
    // Create a simple square lattice for testing
    auto lattice = std::make_shared<SquareLattice>(6, 6, true);
    
    auto qmc = std::make_unique<SSE_QMC>(
        lattice,
        std::move(hamiltonian),
        beta,
        500,
        42
    );
    
    std::cout << "Running simulation..." << std::endl;
    qmc->run_simulation(3000, 10000, 10);
}

int main() {
    std::cout << "=== Different Hamiltonian Models Demo ===" << std::endl;
    
    double beta = 5.0;  // Moderate temperature
    
    // Create a shared lattice for all models
    auto create_lattice = []() { return std::make_shared<SquareLattice>(6, 6, true); };
    
    // 1. Heisenberg Model
    test_model(
        std::make_unique<HeisenbergModel>(create_lattice(), 1.0),
        "Heisenberg Model (J=1.0)",
        beta
    );
    
    // 2. Ising Model
    test_model(
        std::make_unique<IsingModel>(create_lattice(), 1.0),
        "Ising Model (J=1.0)",
        beta
    );
    
    // 3. XY Model
    test_model(
        std::make_unique<XYModel>(create_lattice(), 1.0),
        "XY Model (J=1.0)",
        beta
    );
    
    // 4. XXZ Model (Easy-axis)
    test_model(
        std::make_unique<XXZModel>(create_lattice(), 1.0, 2.0),
        "XXZ Model (J=1.0, Δ=2.0, Easy-axis)",
        beta
    );
    
    // 5. XXZ Model (Easy-plane)
    test_model(
        std::make_unique<XXZModel>(create_lattice(), 1.0, 0.5),
        "XXZ Model (J=1.0, Δ=0.5, Easy-plane)",
        beta
    );
    
    // 6. Heisenberg with Field
    test_model(
        std::make_unique<HeisenbergField>(create_lattice(), 1.0, 0.5),
        "Heisenberg with Field (J=1.0, h=0.5)",
        beta
    );
    
    std::cout << "\n" << std::string(50, '=') << std::endl;
    std::cout << "Model comparison completed!" << std::endl;
    std::cout << "Different models show different physics:" << std::endl;
    std::cout << "- Heisenberg: Full SU(2) symmetry, quantum fluctuations" << std::endl;
    std::cout << "- Ising: Classical spins, sharp phase transition" << std::endl;
    std::cout << "- XY: U(1) symmetry, vortex excitations" << std::endl;
    std::cout << "- XXZ: Interpolates between Ising and XY" << std::endl;
    std::cout << "- Field: Breaks spin symmetry, induces magnetization" << std::endl;
    
    return 0;
}
