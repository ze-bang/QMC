/**
 * @file main.cpp
 * @brief Main executable for SSE QMC simulation
 * 
 * Usage:
 *   sse_qmc [options] config.json
 *   sse_qmc --lattice <type> --L <size> --beta <beta> [options]
 * 
 * Examples:
 *   sse_qmc --lattice square --L 8 --beta 1.0 --model heisenberg
 *   sse_qmc -c config.json
 */

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <nlohmann/json.hpp>
#include <fmt/format.h>
#include <mpi.h>
#include "sse.hpp"

using namespace sse;
using json = nlohmann::json;

void printUsage(const char* prog) {
    std::cout << "SSE Quantum Monte Carlo for Spin-1/2 Systems\n\n";
    std::cout << "Usage:\n";
    std::cout << "  " << prog << " [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  -c, --config <file>     Load configuration from JSON file\n";
    std::cout << "  --lattice <type>        Lattice type: chain, square, triangular, honeycomb, kagome, cubic\n";
    std::cout << "  --L <size>              Linear system size\n";
    std::cout << "  --Lx, --Ly, --Lz <size> System size in each dimension\n";
    std::cout << "  --beta <value>          Inverse temperature\n";
    std::cout << "  --model <type>          Model: heisenberg, xxz, xy, ising\n";
    std::cout << "  --J <value>             Exchange coupling (default: 1.0)\n";
    std::cout << "  --Jz <value>            Z-coupling for XXZ (default: J)\n";
    std::cout << "  --h <value>             Magnetic field (default: 0.0)\n";
    std::cout << "  --therm <sweeps>        Thermalization sweeps (default: 10000)\n";
    std::cout << "  --sweeps <sweeps>       Measurement sweeps (default: 100000)\n";
    std::cout << "  --bins <n>              Number of bins (default: 100)\n";
    std::cout << "  --seed <value>          Random seed (default: 42)\n";
    std::cout << "  --output <file>         Output file for results\n";
    std::cout << "  -v, --verbose           Verbose output\n";
    std::cout << "  -h, --help              Show this help message\n";
}

struct ProgramOptions {
    std::string config_file;
    std::string lattice_type = "square";
    int L = 0;
    int Lx = 0, Ly = 0, Lz = 0;
    Real beta = 1.0;
    std::string model_type = "heisenberg";
    Real J = 1.0;
    Real Jz = 1.0;
    Real h = 0.0;
    int n_therm = 10000;
    int n_sweeps = 100000;
    int n_bins = 100;
    uint64_t seed = 42;
    std::string output_file;
    bool verbose = false;
    bool use_mpi = false;
};

ProgramOptions parseArgs(int argc, char* argv[]) {
    ProgramOptions opts;
    
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        
        if (arg == "-h" || arg == "--help") {
            printUsage(argv[0]);
            std::exit(0);
        } else if (arg == "-c" || arg == "--config") {
            opts.config_file = argv[++i];
        } else if (arg == "--lattice") {
            opts.lattice_type = argv[++i];
        } else if (arg == "--L") {
            opts.L = std::stoi(argv[++i]);
        } else if (arg == "--Lx") {
            opts.Lx = std::stoi(argv[++i]);
        } else if (arg == "--Ly") {
            opts.Ly = std::stoi(argv[++i]);
        } else if (arg == "--Lz") {
            opts.Lz = std::stoi(argv[++i]);
        } else if (arg == "--beta") {
            opts.beta = std::stod(argv[++i]);
        } else if (arg == "--model") {
            opts.model_type = argv[++i];
        } else if (arg == "--J") {
            opts.J = std::stod(argv[++i]);
        } else if (arg == "--Jz") {
            opts.Jz = std::stod(argv[++i]);
        } else if (arg == "--h") {
            opts.h = std::stod(argv[++i]);
        } else if (arg == "--therm") {
            opts.n_therm = std::stoi(argv[++i]);
        } else if (arg == "--sweeps") {
            opts.n_sweeps = std::stoi(argv[++i]);
        } else if (arg == "--bins") {
            opts.n_bins = std::stoi(argv[++i]);
        } else if (arg == "--seed") {
            opts.seed = std::stoull(argv[++i]);
        } else if (arg == "--output") {
            opts.output_file = argv[++i];
        } else if (arg == "-v" || arg == "--verbose") {
            opts.verbose = true;
        } else if (arg == "--mpi") {
            opts.use_mpi = true;
        } else {
            std::cerr << "Unknown option: " << arg << "\n";
            printUsage(argv[0]);
            std::exit(1);
        }
    }
    
    return opts;
}

ProgramOptions loadConfigFromJson(const std::string& filename) {
    ProgramOptions opts;
    
    std::ifstream file(filename);
    if (!file) {
        throw std::runtime_error("Cannot open config file: " + filename);
    }
    
    json j;
    file >> j;
    
    if (j.contains("lattice")) {
        auto& lat = j["lattice"];
        opts.lattice_type = lat.value("type", "square");
        opts.L = lat.value("L", 0);
        opts.Lx = lat.value("Lx", 0);
        opts.Ly = lat.value("Ly", 0);
        opts.Lz = lat.value("Lz", 0);
    }
    
    if (j.contains("model")) {
        auto& model = j["model"];
        opts.model_type = model.value("type", "heisenberg");
        opts.J = model.value("J", 1.0);
        opts.Jz = model.value("Jz", opts.J);
        opts.h = model.value("h", 0.0);
    }
    
    if (j.contains("simulation")) {
        auto& sim = j["simulation"];
        opts.beta = sim.value("beta", 1.0);
        opts.n_therm = sim.value("thermalization", 10000);
        opts.n_sweeps = sim.value("sweeps", 100000);
        opts.n_bins = sim.value("bins", 100);
        opts.seed = sim.value("seed", 42);
    }
    
    opts.output_file = j.value("output", "");
    opts.verbose = j.value("verbose", false);
    
    return opts;
}

Lattice createLattice(const ProgramOptions& opts) {
    int Lx = opts.Lx > 0 ? opts.Lx : opts.L;
    int Ly = opts.Ly > 0 ? opts.Ly : opts.L;
    int Lz = opts.Lz > 0 ? opts.Lz : opts.L;
    
    if (Lx <= 0) {
        throw std::runtime_error("Invalid lattice size");
    }
    
    if (opts.lattice_type == "chain") {
        return Lattice::chain(Lx);
    } else if (opts.lattice_type == "square") {
        return Lattice::square(Lx, Ly > 0 ? Ly : Lx);
    } else if (opts.lattice_type == "triangular") {
        return Lattice::triangular(Lx, Ly > 0 ? Ly : Lx);
    } else if (opts.lattice_type == "honeycomb") {
        return Lattice::honeycomb(Lx, Ly > 0 ? Ly : Lx);
    } else if (opts.lattice_type == "kagome") {
        return Lattice::kagome(Lx, Ly > 0 ? Ly : Lx);
    } else if (opts.lattice_type == "cubic") {
        return Lattice::cubic(Lx, Ly > 0 ? Ly : Lx, Lz > 0 ? Lz : Lx);
    } else {
        throw std::runtime_error("Unknown lattice type: " + opts.lattice_type);
    }
}

Hamiltonian createHamiltonian(const ProgramOptions& opts) {
    if (opts.model_type == "heisenberg") {
        if (opts.h != 0.0) {
            return Hamiltonian::xxzWithField(opts.J, opts.J, opts.h);
        }
        return Hamiltonian::heisenberg(opts.J);
    } else if (opts.model_type == "xxz") {
        if (opts.h != 0.0) {
            return Hamiltonian::xxzWithField(opts.J, opts.Jz, opts.h);
        }
        return Hamiltonian::xxz(opts.J, opts.Jz);
    } else if (opts.model_type == "xy") {
        return Hamiltonian::xy(opts.J);
    } else if (opts.model_type == "ising") {
        return Hamiltonian::ising(opts.J);
    } else {
        throw std::runtime_error("Unknown model type: " + opts.model_type);
    }
}

void saveResults(const std::string& filename, 
                 const ProgramOptions& opts,
                 const Measurements& measurements,
                 const SSESimulation& sim) {
    json j;
    
    // Input parameters
    j["parameters"]["lattice"] = opts.lattice_type;
    j["parameters"]["L"] = opts.L > 0 ? opts.L : opts.Lx;
    j["parameters"]["beta"] = opts.beta;
    j["parameters"]["model"] = opts.model_type;
    j["parameters"]["J"] = opts.J;
    j["parameters"]["Jz"] = opts.Jz;
    j["parameters"]["h"] = opts.h;
    j["parameters"]["sweeps"] = opts.n_sweeps;
    j["parameters"]["thermalization"] = opts.n_therm;
    
    // Results
    auto results = measurements.getResults();
    for (const auto& [name, val] : results) {
        j["results"][name]["mean"] = val.first;
        j["results"][name]["error"] = val.second;
    }
    
    // Statistics
    j["statistics"]["sweeps"] = sim.numSweeps();
    j["statistics"]["acceptance_rate"] = sim.acceptanceRate();
    
    std::ofstream file(filename);
    file << j.dump(2);
}

int main(int argc, char* argv[]) {
    // Initialize MPI
    int mpi_rank = 0, mpi_size = 1;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    
    try {
        // Parse command line arguments
        ProgramOptions opts = parseArgs(argc, argv);
        
        // Load config file if specified
        if (!opts.config_file.empty()) {
            opts = loadConfigFromJson(opts.config_file);
        }
        
        // Validate options
        if (opts.L <= 0 && opts.Lx <= 0) {
            if (mpi_rank == 0) {
                std::cerr << "Error: System size not specified\n";
                printUsage(argv[0]);
            }
            MPI_Finalize();
            return 1;
        }
        
        // Create lattice
        Lattice lattice = createLattice(opts);
        
        // Create Hamiltonian
        Hamiltonian hamiltonian = createHamiltonian(opts);
        
        // Print info on rank 0
        if (mpi_rank == 0 && opts.verbose) {
            std::cout << "=== SSE QMC Simulation ===\n\n";
            lattice.print();
            std::cout << "\n";
            hamiltonian.print();
            std::cout << "\nSimulation parameters:\n";
            std::cout << fmt::format("  β = {:.4f} (T = {:.4f})\n", opts.beta, 1.0/opts.beta);
            std::cout << fmt::format("  Thermalization: {} sweeps\n", opts.n_therm);
            std::cout << fmt::format("  Measurements: {} sweeps\n", opts.n_sweeps);
            std::cout << fmt::format("  MPI processes: {}\n\n", mpi_size);
        }
        
        // Set up simulation parameters
        SimulationParams params;
        params.beta = opts.beta;
        params.n_therm = opts.n_therm;
        params.n_sweeps = opts.n_sweeps / mpi_size;  // Divide work
        params.n_bins = opts.n_bins;
        params.seed = opts.seed + mpi_rank * 12345;  // Different seed per rank
        
        // Create and run simulation
        SSESimulation sim(lattice, hamiltonian, params);
        
        if (mpi_rank == 0 && opts.verbose) {
            std::cout << "Initializing...\n";
        }
        sim.initialize();
        
        if (mpi_rank == 0 && opts.verbose) {
            std::cout << "Thermalizing...\n";
        }
        sim.thermalize();
        
        if (mpi_rank == 0 && opts.verbose) {
            std::cout << "Running production...\n";
        }
        sim.run();
        
        // Collect results from all ranks (simplified - just rank 0 output)
        if (mpi_rank == 0) {
            sim.getMeasurements().print();
            sim.printStatus();
            
            if (!opts.output_file.empty()) {
                saveResults(opts.output_file, opts, sim.getMeasurements(), sim);
                std::cout << "\nResults saved to: " << opts.output_file << "\n";
            }
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        MPI_Finalize();
        return 1;
    }
    
    MPI_Finalize();
    return 0;
}
