/**
 * @file benchmark.cpp
 * @brief Performance benchmarks for SSE QMC
 */

#include <iostream>
#include <chrono>
#include <iomanip>
#include "sse.hpp"

using namespace sse;
using namespace std::chrono;

void benchmark_lattice_creation() {
    std::cout << "=== Lattice Creation Benchmark ===\n";
    
    for (int L : {16, 32, 64, 128}) {
        auto start = high_resolution_clock::now();
        for (int i = 0; i < 100; ++i) {
            auto lat = Lattice::square(L, L);
        }
        auto end = high_resolution_clock::now();
        double time_us = duration_cast<microseconds>(end - start).count() / 100.0;
        std::cout << "Square " << L << "x" << L << ": " << time_us << " μs\n";
    }
    std::cout << "\n";
}

void benchmark_diagonal_update() {
    std::cout << "=== Diagonal Update Benchmark ===\n";
    
    for (int L : {8, 16, 32}) {
        auto lattice = Lattice::square(L, L);
        auto H = Hamiltonian::heisenberg(1.0);
        
        SimulationParams params;
        params.beta = 1.0;
        params.n_therm = 100;
        params.seed = 42;
        
        SSESimulation sim(lattice, H, params);
        sim.initialize();
        sim.thermalize();
        
        int n_sweeps = 1000;
        auto start = high_resolution_clock::now();
        for (int i = 0; i < n_sweeps; ++i) {
            sim.diagonalUpdate();
        }
        auto end = high_resolution_clock::now();
        
        double time_ms = duration_cast<microseconds>(end - start).count() / 1000.0;
        double time_per_sweep = time_ms / n_sweeps;
        double sweeps_per_sec = n_sweeps / (time_ms / 1000.0);
        
        std::cout << "L=" << L << " (" << lattice.numSites() << " sites): "
                  << std::fixed << std::setprecision(3) 
                  << time_per_sweep << " ms/sweep, "
                  << sweeps_per_sec << " sweeps/s\n";
    }
    std::cout << "\n";
}

void benchmark_loop_update() {
    std::cout << "=== Loop Update Benchmark ===\n";
    
    for (int L : {8, 16, 32}) {
        auto lattice = Lattice::square(L, L);
        auto H = Hamiltonian::heisenberg(1.0);
        
        SimulationParams params;
        params.beta = 2.0;  // Lower T = more operators
        params.n_therm = 200;
        params.seed = 42;
        
        SSESimulation sim(lattice, H, params);
        sim.initialize();
        sim.thermalize();
        
        // Prepare for loop update
        sim.diagonalUpdate();
        sim.getConfig().buildVertexList();
        
        int n_updates = 1000;
        auto start = high_resolution_clock::now();
        for (int i = 0; i < n_updates; ++i) {
            sim.loopUpdate();
        }
        auto end = high_resolution_clock::now();
        
        double time_ms = duration_cast<microseconds>(end - start).count() / 1000.0;
        double time_per_update = time_ms / n_updates;
        
        std::cout << "L=" << L << " (" << lattice.numSites() << " sites): "
                  << std::fixed << std::setprecision(3)
                  << time_per_update << " ms/update\n";
    }
    std::cout << "\n";
}

void benchmark_full_sweep() {
    std::cout << "=== Full Sweep Benchmark ===\n";
    
    for (int L : {8, 16, 32, 64}) {
        auto lattice = Lattice::square(L, L);
        auto H = Hamiltonian::heisenberg(1.0);
        
        SimulationParams params;
        params.beta = 1.0;
        params.n_therm = 100;
        params.n_sweeps = 1000;
        params.n_bins = 10;
        params.seed = 42;
        
        SSESimulation sim(lattice, H, params);
        sim.initialize();
        sim.thermalize();
        
        auto start = high_resolution_clock::now();
        sim.run();
        auto end = high_resolution_clock::now();
        
        double time_s = duration_cast<milliseconds>(end - start).count() / 1000.0;
        double sweeps_per_sec = params.n_sweeps / time_s;
        
        std::cout << "L=" << L << " (" << lattice.numSites() << " sites): "
                  << std::fixed << std::setprecision(1)
                  << sweeps_per_sec << " sweeps/s, "
                  << "total " << time_s << "s\n";
    }
    std::cout << "\n";
}

void benchmark_scaling() {
    std::cout << "=== Scaling with System Size ===\n";
    std::cout << "L\tN\tn_avg\ttime/sweep\n";
    
    for (int L : {4, 8, 12, 16, 24, 32}) {
        auto lattice = Lattice::square(L, L);
        auto H = Hamiltonian::heisenberg(1.0);
        
        SimulationParams params;
        params.beta = 1.0;
        params.n_therm = 100;
        params.n_sweeps = 500;
        params.seed = 42;
        
        SSESimulation sim(lattice, H, params);
        sim.initialize();
        sim.thermalize();
        
        auto start = high_resolution_clock::now();
        sim.run();
        auto end = high_resolution_clock::now();
        
        double time_ms = duration_cast<microseconds>(end - start).count() / 1000.0;
        double time_per_sweep = time_ms / params.n_sweeps;
        double n_avg = sim.getMeasurements().avgOperators();
        
        std::cout << L << "\t" << L*L << "\t" 
                  << std::fixed << std::setprecision(0) << n_avg << "\t"
                  << std::setprecision(3) << time_per_sweep << " ms\n";
    }
}

int main() {
    std::cout << "SSE QMC Performance Benchmarks\n";
    std::cout << "==============================\n\n";
    
    benchmark_lattice_creation();
    benchmark_diagonal_update();
    benchmark_loop_update();
    benchmark_full_sweep();
    benchmark_scaling();
    
    return 0;
}
