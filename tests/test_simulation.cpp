/**
 * @file test_simulation.cpp
 * @brief Tests for SSE simulation
 */

#include <iostream>
#include <cmath>
#include "sse.hpp"

extern int g_test_count;
extern int g_test_passed;
extern int g_test_failed;

#define TEST_ASSERT(condition, message) \
    do { \
        ++g_test_count; \
        if (condition) { \
            ++g_test_passed; \
        } else { \
            ++g_test_failed; \
            std::cerr << "FAIL: " << message << " at " << __FILE__ << ":" << __LINE__ << "\n"; \
        } \
    } while(0)

#define TEST_ASSERT_NEAR(a, b, tol, message) \
    TEST_ASSERT(std::abs((a) - (b)) < (tol), message)

void test_simulation() {
    std::cout << "Testing SSE Simulation...\n";
    
    // Test basic simulation setup
    {
        auto lattice = sse::Lattice::chain(4);
        auto H = sse::Hamiltonian::heisenberg(1.0);
        
        sse::SimulationParams params;
        params.beta = 1.0;
        params.n_therm = 100;
        params.n_sweeps = 100;
        params.n_bins = 10;
        params.seed = 12345;
        
        sse::SSESimulation sim(lattice, H, params);
        sim.initialize();
        
        TEST_ASSERT(sim.getConfig().numOperators() == 0, 
                   "Initial configuration should have 0 operators");
        
        TEST_ASSERT(sim.getConfig().operatorStringLength() >= 4,
                   "Operator string should have minimum length");
    }
    
    // Test diagonal update
    {
        auto lattice = sse::Lattice::chain(4);
        auto H = sse::Hamiltonian::heisenberg(1.0);
        
        sse::SimulationParams params;
        params.beta = 2.0;
        params.n_therm = 100;
        params.n_sweeps = 100;
        params.seed = 12345;
        
        sse::SSESimulation sim(lattice, H, params);
        sim.initialize();
        sim.thermalize();
        
        // After thermalization, should have some operators
        TEST_ASSERT(sim.getConfig().numOperators() > 0,
                   "Should have operators after thermalization");
        
        // Average n should scale with beta * N_bonds
        double expected_n = params.beta * lattice.numBonds() * 0.5;  // Rough estimate
        double actual_n = static_cast<double>(sim.getConfig().numOperators());
        TEST_ASSERT(actual_n > expected_n * 0.1 && actual_n < expected_n * 10,
                   "Number of operators should be reasonable");
    }
    
    // Test short simulation produces measurements
    {
        auto lattice = sse::Lattice::square(4, 4);
        auto H = sse::Hamiltonian::heisenberg(1.0);
        
        sse::SimulationParams params;
        params.beta = 1.0;
        params.n_therm = 100;
        params.n_sweeps = 1000;
        params.n_bins = 10;
        params.seed = 42;
        
        sse::SSESimulation sim(lattice, H, params);
        sim.initialize();
        sim.thermalize();
        sim.run();
        
        const auto& meas = sim.getMeasurements();
        
        // Energy should be negative for antiferromagnet
        TEST_ASSERT(meas.energy() < 0, "Energy should be negative for Heisenberg");
        
        // Magnetization squared should be small for antiferromagnet at finite T
        TEST_ASSERT(meas.magnetizationSquared() >= 0, "m^2 should be non-negative");
        TEST_ASSERT(meas.magnetizationSquared() < 0.5, "m^2 should be small for AFM");
        
        // Error bars should be positive
        TEST_ASSERT(meas.energyError() >= 0, "Error should be non-negative");
    }
    
    // Test different lattices
    {
        sse::SimulationParams params;
        params.beta = 0.5;
        params.n_therm = 50;
        params.n_sweeps = 100;
        params.n_bins = 5;
        
        auto H = sse::Hamiltonian::heisenberg(1.0);
        
        // Square lattice
        {
            auto lat = sse::Lattice::square(4, 4);
            sse::SSESimulation sim(lat, H, params);
            sim.initialize();
            sim.thermalize();
            TEST_ASSERT(sim.getConfig().numOperators() > 0, "Square lattice should work");
        }
        
        // Triangular lattice
        {
            auto lat = sse::Lattice::triangular(4, 4);
            sse::SSESimulation sim(lat, H, params);
            sim.initialize();
            sim.thermalize();
            TEST_ASSERT(sim.getConfig().numOperators() > 0, "Triangular lattice should work");
        }
        
        // Honeycomb lattice
        {
            auto lat = sse::Lattice::honeycomb(4, 4);
            sse::SSESimulation sim(lat, H, params);
            sim.initialize();
            sim.thermalize();
            TEST_ASSERT(sim.getConfig().numOperators() > 0, "Honeycomb lattice should work");
        }
    }
    
    // Test different models
    {
        auto lattice = sse::Lattice::square(4, 4);
        
        sse::SimulationParams params;
        params.beta = 1.0;
        params.n_therm = 50;
        params.n_sweeps = 50;
        
        // XXZ model
        {
            auto H = sse::Hamiltonian::xxz(1.0, 0.5);
            sse::SSESimulation sim(lattice, H, params);
            sim.initialize();
            sim.thermalize();
            TEST_ASSERT(sim.getConfig().numOperators() > 0, "XXZ model should work");
        }
        
        // XY model
        {
            auto H = sse::Hamiltonian::xy(1.0);
            sse::SSESimulation sim(lattice, H, params);
            sim.initialize();
            sim.thermalize();
            TEST_ASSERT(sim.getConfig().numOperators() > 0, "XY model should work");
        }
    }
    
    // Test random number generator
    {
        sse::Random rng(42);
        
        double sum = 0.0;
        int N = 10000;
        for (int i = 0; i < N; ++i) {
            double r = rng.uniform01();
            TEST_ASSERT(r >= 0.0 && r < 1.0, "Random should be in [0,1)");
            sum += r;
        }
        double mean = sum / N;
        TEST_ASSERT_NEAR(mean, 0.5, 0.05, "Random mean should be ~0.5");
    }
    
    std::cout << "Simulation tests completed.\n\n";
}
