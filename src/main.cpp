#include <iostream>
#include <memory>
#include <string>
#include <fstream>
#include <chrono>
#include <iomanip>
#include "sse_qmc.h"
#include "lattice.h"
#include "hamiltonian.h"

using namespace SSE;

void print_usage(const char* program_name) {
    std::cout << "Usage: " << program_name << " [options]\n"
              << "Options:\n"
              << "  --lattice <type>     Lattice type: square, triangular, chain, honeycomb\n"
              << "  --model <type>       Model type: heisenberg, ising, xy, xxz\n"
              << "  --Lx <size>          Lattice size in x direction\n"
              << "  --Ly <size>          Lattice size in y direction (for 2D lattices)\n"
              << "  --periodic           Use periodic boundary conditions\n"
              << "  --J <coupling>       Exchange coupling strength\n"
              << "  --Delta <anisotropy> Anisotropy parameter (for XXZ model)\n"
              << "  --h <field>          Magnetic field strength\n"
              << "  --beta <temp>        Inverse temperature\n"
              << "  --therm <steps>      Thermalization steps\n"
              << "  --meas <steps>       Measurement steps\n"
              << "  --interval <int>     Measurement interval\n"
              << "  --seed <seed>        Random seed\n"
              << "  --output <file>      Output file prefix\n"
              << "  --help               Show this help message\n"
              << std::endl;
}

struct Parameters {
    std::string lattice_type = "square";
    std::string model_type = "heisenberg";
    int Lx = 8;
    int Ly = 8;
    bool periodic = true;
    double J = 1.0;
    double Delta = 1.0;
    double h = 0.0;
    double beta = 10.0;
    long long therm_steps = 10000;
    long long meas_steps = 50000;
    int meas_interval = 1;
    unsigned int seed = 42;
    std::string output_prefix = "sse_qmc";
};

void parse_arguments(int argc, char* argv[], Parameters& params) {
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        
        if (arg == "--help") {
            print_usage(argv[0]);
            exit(0);
        } else if (arg == "--lattice" && i + 1 < argc) {
            params.lattice_type = argv[++i];
        } else if (arg == "--model" && i + 1 < argc) {
            params.model_type = argv[++i];
        } else if (arg == "--Lx" && i + 1 < argc) {
            params.Lx = std::stoi(argv[++i]);
        } else if (arg == "--Ly" && i + 1 < argc) {
            params.Ly = std::stoi(argv[++i]);
        } else if (arg == "--periodic") {
            params.periodic = true;
        } else if (arg == "--open") {
            params.periodic = false;
        } else if (arg == "--J" && i + 1 < argc) {
            params.J = std::stod(argv[++i]);
        } else if (arg == "--Delta" && i + 1 < argc) {
            params.Delta = std::stod(argv[++i]);
        } else if (arg == "--h" && i + 1 < argc) {
            params.h = std::stod(argv[++i]);
        } else if (arg == "--beta" && i + 1 < argc) {
            params.beta = std::stod(argv[++i]);
        } else if (arg == "--therm" && i + 1 < argc) {
            params.therm_steps = std::stoll(argv[++i]);
        } else if (arg == "--meas" && i + 1 < argc) {
            params.meas_steps = std::stoll(argv[++i]);
        } else if (arg == "--interval" && i + 1 < argc) {
            params.meas_interval = std::stoi(argv[++i]);
        } else if (arg == "--seed" && i + 1 < argc) {
            params.seed = static_cast<unsigned int>(std::stoul(argv[++i]));
        } else if (arg == "--output" && i + 1 < argc) {
            params.output_prefix = argv[++i];
        } else {
            std::cerr << "Unknown argument: " << arg << std::endl;
            print_usage(argv[0]);
            exit(1);
        }
    }
}

std::unique_ptr<Lattice> create_lattice(const Parameters& params) {
    if (params.lattice_type == "square") {
        return std::make_unique<SquareLattice>(params.Lx, params.Ly, params.periodic);
    } else if (params.lattice_type == "triangular") {
        return std::make_unique<TriangularLattice>(params.Lx, params.Ly, params.periodic);
    } else if (params.lattice_type == "chain") {
        return std::make_unique<Chain>(params.Lx, params.periodic);
    } else if (params.lattice_type == "honeycomb") {
        return std::make_unique<HoneycombLattice>(params.Lx, params.Ly, params.periodic);
    } else {
        throw std::invalid_argument("Unknown lattice type: " + params.lattice_type);
    }
}

std::unique_ptr<Hamiltonian> create_hamiltonian(std::shared_ptr<Lattice> lattice, const Parameters& params) {
    if (params.model_type == "heisenberg") {
        if (params.h != 0.0) {
            return std::make_unique<HeisenbergField>(lattice, params.J, params.h);
        } else {
            return std::make_unique<HeisenbergModel>(lattice, params.J);
        }
    } else if (params.model_type == "ising") {
        return std::make_unique<IsingModel>(lattice, params.J);
    } else if (params.model_type == "xy") {
        return std::make_unique<XYModel>(lattice, params.J);
    } else if (params.model_type == "xxz") {
        return std::make_unique<XXZModel>(lattice, params.J, params.Delta);
    } else {
        throw std::invalid_argument("Unknown model type: " + params.model_type);
    }
}

void save_parameters(const Parameters& params, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Warning: Could not save parameters to " << filename << std::endl;
        return;
    }
    
    file << "# SSE QMC Parameters\n";
    file << "lattice_type: " << params.lattice_type << "\n";
    file << "model_type: " << params.model_type << "\n";
    file << "Lx: " << params.Lx << "\n";
    file << "Ly: " << params.Ly << "\n";
    file << "periodic: " << (params.periodic ? "true" : "false") << "\n";
    file << "J: " << params.J << "\n";
    file << "Delta: " << params.Delta << "\n";
    file << "h: " << params.h << "\n";
    file << "beta: " << params.beta << "\n";
    file << "therm_steps: " << params.therm_steps << "\n";
    file << "meas_steps: " << params.meas_steps << "\n";
    file << "meas_interval: " << params.meas_interval << "\n";
    file << "seed: " << params.seed << "\n";
    
    file.close();
}

int main(int argc, char* argv[]) {
    std::cout << "=== SSE Quantum Monte Carlo ===" << std::endl;
    std::cout << "Stochastic Series Expansion for Spin-1/2 Systems" << std::endl;
    std::cout << std::string(50, '=') << std::endl;
    
    // Parse command line arguments
    Parameters params;
    parse_arguments(argc, argv, params);
    
    try {
        // Create lattice
        auto lattice = create_lattice(params);
        lattice->print_info();
        
        // Create Hamiltonian
        auto shared_lattice = std::shared_ptr<Lattice>(lattice.release());
        auto hamiltonian = create_hamiltonian(shared_lattice, params);
        hamiltonian->print_info();
        
        // Create SSE QMC instance
        auto qmc = std::make_unique<SSE_QMC>(
            shared_lattice,
            std::move(hamiltonian),
            params.beta,
            1000,  // Initial cutoff
            params.seed
        );
        
        // Save parameters
        save_parameters(params, params.output_prefix + "_params.txt");
        
        // Run simulation
        auto start_time = std::chrono::high_resolution_clock::now();
        
        qmc->run_simulation(params.therm_steps, params.meas_steps, params.meas_interval);
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
        
        std::cout << "\nSimulation completed in " << duration.count() << " seconds" << std::endl;
        
        // Save results
        // Note: The observables are owned by the QMC object, so we'd need to add getter methods
        // For now, results are printed to console
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}

// Example function for running a temperature scan
void run_temperature_scan() {
    std::cout << "\n=== Temperature Scan Example ===" << std::endl;
    
    // Parameters
    int Lx = 8, Ly = 8;
    double J = 1.0;
    std::vector<double> temperatures = {0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0};
    
    // Create lattice and Hamiltonian
    auto lattice = std::make_unique<SquareLattice>(Lx, Ly, true);
    auto shared_lattice = std::shared_ptr<Lattice>(lattice.release());
    auto hamiltonian = std::make_unique<HeisenbergModel>(shared_lattice, J);
    
    std::cout << "Temperature\tEnergy\t\tSpecific_Heat\tSusceptibility" << std::endl;
    std::cout << std::string(70, '-') << std::endl;
    
    for (double T : temperatures) {
        double beta = 1.0 / T;
        
        auto qmc = std::make_unique<SSE_QMC>(
            std::unique_ptr<Lattice>(shared_lattice.get()),
            std::unique_ptr<Hamiltonian>(hamiltonian.get()),
            beta,
            1000,
            42
        );
        
        // Quick simulation for demo
        qmc->run_simulation(5000, 10000, 10);
        
        std::cout << std::fixed << std::setprecision(3)
                  << T << "\t\t"
                  << qmc->get_energy() << "\t\t"
                  << "N/A" << "\t\t"  // Would need to implement getter for specific heat
                  << "N/A" << std::endl;  // Would need to implement getter for susceptibility
    }
    
    // Prevent double delete
    // shared_lattice.release();
    // hamiltonian.release();
}

// Example for different lattice types
void demonstrate_lattices() {
    std::cout << "\n=== Lattice Types Demonstration ===" << std::endl;
    
    std::vector<std::unique_ptr<Lattice>> lattices;
    lattices.push_back(std::make_unique<SquareLattice>(4, 4, true));
    lattices.push_back(std::make_unique<TriangularLattice>(4, 4, true));
    lattices.push_back(std::make_unique<Chain>(16, true));
    lattices.push_back(std::make_unique<HoneycombLattice>(3, 3, true));
    
    for (auto& lattice : lattices) {
        std::cout << "\n" << std::string(40, '-') << std::endl;
        lattice->print_info();
    }
}
