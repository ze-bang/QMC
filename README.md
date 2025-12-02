# SSE Quantum Monte Carlo for Spin-1/2 Systems

A high-performance C++ implementation of the Stochastic Series Expansion (SSE) Quantum Monte Carlo algorithm for spin-1/2 systems with arbitrary lattice geometries and nearest-neighbor interaction Hamiltonians.

## Features

- **Arbitrary Lattices**: Chain, square, triangular, honeycomb, kagome, cubic, and custom lattices
- **General Hamiltonians**: Heisenberg, XXZ, XY, Ising, and custom 4×4 bond matrices
- **Efficient Updates**: Diagonal updates and directed loop algorithm
- **Measurements**: Energy, specific heat, magnetization, susceptibility, correlations
- **High Performance**: 
  - OpenMP parallelization
  - MPI support for parallel tempering and distributed computing
  - PCG random number generator for high-quality statistics
  - Optimized data structures
- **Extensible**: Easy to add new lattices, models, and measurements

## Requirements

- C++20 compiler (GCC 10+, Clang 12+, or MSVC 2019+)
- CMake 3.16+
- Eigen3 (linear algebra)
- OpenMP (parallelization)
- MPI (distributed computing)
- Optional: HDF5 (checkpointing)

## Building

```bash
# Clone the repository
git clone <repository-url>
cd QMC_cpp

# Create build directory
mkdir build && cd build

# Configure with CMake
cmake .. -DCMAKE_BUILD_TYPE=Release

# Build
cmake --build . -j$(nproc)

# Run tests
ctest --output-on-failure

# Install (optional)
sudo cmake --install .
```

### Build Options

- `CMAKE_BUILD_TYPE`: `Release` (optimized), `Debug` (with sanitizers), `RelWithDebInfo`
- Custom compiler: `cmake .. -DCMAKE_CXX_COMPILER=g++-12`

## Usage

### Command Line

```bash
# Simple Heisenberg model on 8x8 square lattice at β=2
./sse_qmc --lattice square --L 8 --beta 2.0 --model heisenberg

# With output file
./sse_qmc --lattice square --L 16 --beta 1.0 --sweeps 100000 --output results.json

# From configuration file
./sse_qmc -c examples/heisenberg_square.json

# Using MPI for parallel runs
mpirun -np 4 ./sse_qmc --lattice square --L 16 --beta 2.0
```

### Configuration File (JSON)

```json
{
    "lattice": {
        "type": "square",
        "Lx": 16,
        "Ly": 16
    },
    "model": {
        "type": "heisenberg",
        "J": 1.0
    },
    "simulation": {
        "beta": 2.0,
        "thermalization": 10000,
        "sweeps": 100000,
        "bins": 100,
        "seed": 42
    },
    "output": "results.json",
    "verbose": true
}
```

### Available Options

**Lattice types:**
- `chain`: 1D chain
- `square`: 2D square lattice
- `triangular`: 2D triangular lattice
- `honeycomb`: 2D honeycomb lattice
- `kagome`: 2D kagome lattice
- `cubic`: 3D cubic lattice

**Model types:**
- `heisenberg`: H = J(S_x·S_x + S_y·S_y + S_z·S_z)
- `xxz`: H = J_xy(S_x·S_x + S_y·S_y) + J_z·S_z·S_z
- `xy`: H = J(S_x·S_x + S_y·S_y)
- `ising`: H = J·S_z·S_z

## Algorithm

The SSE algorithm represents the partition function as:

$$Z = \text{Tr}\left[\sum_{n=0}^{\infty} \frac{\beta^n}{n!} (-H)^n\right]$$

The simulation samples operator strings using:

1. **Diagonal Update**: Insert/remove diagonal operators with Metropolis acceptance
2. **Loop Update**: Modify operator string using directed loop algorithm for efficient sampling

### Key References

1. A. W. Sandvik, "Stochastic series expansion method with operator-loop update", Phys. Rev. B 59, R14157 (1999)
2. O. F. Syljuåsen and A. W. Sandvik, "Quantum Monte Carlo with directed loops", Phys. Rev. E 66, 046701 (2002)
3. A. W. Sandvik, "Computational Studies of Quantum Spin Systems", AIP Conf. Proc. 1297, 135 (2010)

## Library API

```cpp
#include <sse.hpp>

using namespace sse;

// Create lattice
auto lattice = Lattice::square(16, 16);

// Create Hamiltonian
auto H = Hamiltonian::heisenberg(1.0);  // J = 1

// Set up simulation
SimulationParams params;
params.beta = 2.0;
params.n_therm = 10000;
params.n_sweeps = 100000;

SSESimulation sim(lattice, H, params);
sim.initialize();
sim.thermalize();
sim.run();

// Get results
const auto& meas = sim.getMeasurements();
std::cout << "Energy: " << meas.energy() << " ± " << meas.energyError() << "\n";
```

### Custom Hamiltonian

```cpp
// Create custom 4x4 bond Hamiltonian matrix
BondMatrix H_custom;
H_custom << 0.25,  0,    0,    0,
            0,    -0.25, 0.5,  0,
            0,     0.5, -0.25, 0,
            0,     0,    0,    0.25;

auto H = Hamiltonian::custom(H_custom);
```

### Custom Lattice

```cpp
// Define bonds manually
std::vector<Bond> bonds = {
    {0, 1, 0},  // Site 0 to Site 1, type 0
    {1, 2, 0},
    {2, 0, 1},  // Different bond type
    // ...
};

auto lattice = Lattice::custom(n_sites, bonds);
```

## Output Format

Results are saved in JSON format:

```json
{
    "parameters": {
        "lattice": "square",
        "L": 16,
        "beta": 2.0,
        "model": "heisenberg",
        "J": 1.0
    },
    "results": {
        "energy": {"mean": -0.6694, "error": 0.0002},
        "specific_heat": {"mean": 0.127, "error": 0.003},
        "magnetization_sq": {"mean": 0.0156, "error": 0.0001},
        "susceptibility": {"mean": 0.985, "error": 0.015},
        "stag_magnetization": {"mean": 0.0892, "error": 0.0003}
    }
}
```

## Performance

Typical performance on modern hardware (single core):

| System Size | Sweeps/second |
|-------------|---------------|
| 8×8         | ~10,000       |
| 16×16       | ~2,500        |
| 32×32       | ~600          |
| 64×64       | ~150          |

Performance scales approximately as O(N) where N is the number of sites.

## Directory Structure

```
QMC_cpp/
├── CMakeLists.txt          # Build configuration
├── README.md               # This file
├── include/
│   └── sse/
│       ├── types.hpp       # Type definitions
│       ├── random.hpp      # Random number generator
│       ├── lattice.hpp     # Lattice structures
│       ├── hamiltonian.hpp # Hamiltonian definitions
│       ├── vertex.hpp      # Vertex data for SSE
│       ├── sse_config.hpp  # Configuration state
│       ├── measurements.hpp# Physical observables
│       └── sse_simulation.hpp # Main simulation class
├── src/
│   ├── main.cpp            # Command-line interface
│   ├── lattice.cpp
│   ├── hamiltonian.cpp
│   ├── vertex.cpp
│   ├── sse_config.cpp
│   ├── sse_simulation.cpp
│   ├── measurements.cpp
│   └── random.cpp
├── tests/
│   ├── test_main.cpp
│   ├── test_lattice.cpp
│   ├── test_hamiltonian.cpp
│   └── test_simulation.cpp
├── benchmarks/
│   └── benchmark.cpp
└── examples/
    ├── heisenberg_square.json
    └── heisenberg_triangular.json
```

## Extending the Code

### Adding a New Lattice

1. Add static factory method to `Lattice` class in `include/sse/lattice.hpp`
2. Implement in `src/lattice.cpp`

### Adding a New Model

1. Add static factory method to `Hamiltonian` class in `include/sse/hamiltonian.hpp`
2. Implement in `src/hamiltonian.cpp`

### Adding New Measurements

1. Add observable to `Measurements` class in `include/sse/measurements.hpp`
2. Implement measurement in `src/measurements.cpp`

## License

MIT License

## Citation

If you use this code, please cite:

```bibtex
@software{sse_qmc,
    title = {SSE Quantum Monte Carlo for Spin Systems},
    author = {Your Name},
    year = {2024},
    url = {https://github.com/your-repo}
}
```

## Contributing

Contributions are welcome! Please open an issue or submit a pull request.
