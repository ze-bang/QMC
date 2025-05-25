# QMC_cpp: Quantum Monte Carlo for Spin Hamiltonians

This library provides a modern C++ implementation of Quantum Monte Carlo (QMC) algorithms for simulating spin Hamiltonians on arbitrary lattices. It focuses on high performance, flexibility, and ease of use.

## Features

- **Arbitrary Lattices**: Support for square, triangular, and custom lattice structures (including Kagome)
- **Multiple Spin Models**:
  - Classical Ising Model
  - Classical XY Model
  - Classical Heisenberg Model
  - Quantum Heisenberg Model (S=1/2)
- **State-of-the-Art QMC Algorithms**:
  - Stochastic Series Expansion (SSE) with operator-loop updates
  - Path Integral Monte Carlo (PIMC) with cluster updates
  - Projector Quantum Monte Carlo for ground state properties
- **Advanced Measurements**:
  - Energy and magnetization with error estimates
  - Specific heat and susceptibility
  - Correlation functions
  - Custom observables
- **Highly Configurable**: Easy to customize parameters and observables
- **Parallelization**: OpenMP support for parallel simulations
- **Modern C++**: Written in C++17 with clean, object-oriented design

## Requirements

- C++17 compiler (GCC 7+, Clang 5+, or MSVC 2019+)
- CMake 3.10 or higher
- Eigen3 library (optional, for enhanced linear algebra operations)
- Boost libraries (optional, for program options and serialization)
- OpenMP (optional, for parallelization)

## Building

You can build the project using the provided build script:

```bash
./build.sh
```

Or manually with CMake:

```bash
mkdir build
cd build
cmake ..
make
```

## Usage Examples

### Basic Example

```cpp
#include "lattice.hpp"
#include "spin_model.hpp"
#include "qmc.hpp"

int main() {
    // Create a 8x8 square lattice with periodic boundary conditions
    auto lattice = std::make_shared<qmc::SquareLattice>(8, 8, true);
    
    // Create a quantum Heisenberg model
    auto model = std::make_shared<qmc::QuantumHeisenbergModel>(lattice);
    
    // Set up SSE parameters
    qmc::StochasticSeriesExpansion::SSEParameters params;
    params.temperature = 0.5;
    params.warmupSteps = 10000;
    params.measurementSteps = 50000;
    
    // Create and run the simulation
    qmc::StochasticSeriesExpansion simulation(model, params);
    simulation.initialize();
    simulation.warmup();
    simulation.run();
    
    // Print results
    std::cout << "Energy per site: " << simulation.getEnergy() << std::endl;
    std::cout << "Magnetization: " << simulation.getMagnetization() << std::endl;
    std::cout << "Specific Heat: " << simulation.getSpecificHeat() << std::endl;
    std::cout << "Susceptibility: " << simulation.getSusceptibility() << std::endl;
    
    return 0;
}
```

## Theoretical Background

### Quantum Monte Carlo

Quantum Monte Carlo methods are powerful numerical techniques for studying quantum many-body systems. Unlike classical Monte Carlo, QMC methods can simulate quantum mechanical effects like tunneling, entanglement, and zero-point motion.

### Stochastic Series Expansion (SSE)

SSE is based on a Taylor series expansion of the partition function:

$Z = \text{Tr}[e^{-\beta H}] = \sum_{\alpha} \sum_{n=0}^{\infty} \frac{\beta^n}{n!} (-1)^n \langle\alpha|H^n|\alpha\rangle$

Where $|\alpha\rangle$ are basis states and $H$ is the Hamiltonian. The algorithm samples terms in this expansion using Metropolis Monte Carlo.

### Path Integral Monte Carlo (PIMC)

PIMC maps a quantum system to a classical system with an extra dimension (imaginary time). It uses the Trotter decomposition:

$e^{-\beta H} \approx (e^{-\beta H/M})^M = (e^{-\beta V/2M} e^{-\beta K/M} e^{-\beta V/2M})^M + O(\beta^3/M^2)$

Where $H = K + V$ is the Hamiltonian split into kinetic and potential terms.

## Example Programs

The repository contains several example programs:

1. **simple_sse** - A basic example of the Stochastic Series Expansion algorithm
2. **heisenberg_afm** - A temperature scan for the antiferromagnetic Heisenberg model
3. **custom_lattice** - Demonstrates simulating a Kagome lattice with geometric frustration
4. **path_integral_mc** - Demonstrates the Path Integral Monte Carlo algorithm

To run an example:

```bash
cd build
./examples/simple_sse
```

## Advanced Usage

### Custom Lattices

The library supports defining custom lattice structures:

```cpp
// Create a Kagome lattice
std::vector<qmc::Lattice::Position> positions = {...};
std::vector<std::vector<qmc::Lattice::Index>> neighbors = {...};
auto lattice = std::make_shared<qmc::CustomLattice>(positions, neighbors);
```

### Custom Observables

You can register custom observables to measure during the simulation:

```cpp
simulation.registerObservable("StaggeredMagnetization", [&lattice](const std::vector<double>& config) {
    double m = 0.0;
    for (size_t i = 0; i < config.size(); ++i) {
        auto coords = std::dynamic_pointer_cast<qmc::SquareLattice>(lattice)->coordinates(i);
        double sign = (coords.first + coords.second) % 2 == 0 ? 1.0 : -1.0;
        m += sign * 0.5 * config[i];
    }
    return std::abs(m) / lattice->size();
});
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

This library implements algorithms described in:
- A.W. Sandvik, "Stochastic series expansion method with operator-loop update", Phys. Rev. B 59, R14157 (1999)
- O.F. Syljuåsen and A.W. Sandvik, "Quantum Monte Carlo with directed loops", Phys. Rev. E 66, 046701 (2002)
- M. Troyer, F. Alet, S. Trebst, "Path Integral Monte Carlo Simulation of Quantum Spin Models", AIP Conference Proceedings 690, 156 (2003)
- N. Kawashima and K. Harada, "Recent Developments of World-Line Monte Carlo Methods", J. Phys. Soc. Jpn. 73, 1379 (2004)