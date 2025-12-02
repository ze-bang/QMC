/**
 * @file test_hamiltonian.cpp
 * @brief Tests for Hamiltonian and vertex data
 */

#include <iostream>
#include <cmath>
#include "sse/hamiltonian.hpp"
#include "sse/vertex.hpp"

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

void test_hamiltonian() {
    std::cout << "Testing Hamiltonian...\n";
    
    // Test Heisenberg Hamiltonian
    {
        auto H = sse::Hamiltonian::heisenberg(1.0);
        const auto& mat = H.getBondMatrix(0);
        
        // Check diagonal elements
        TEST_ASSERT_NEAR(mat(0, 0), 0.25, 1e-10, "H(|00>,|00>) should be 1/4");
        TEST_ASSERT_NEAR(mat(1, 1), -0.25, 1e-10, "H(|01>,|01>) should be -1/4");
        TEST_ASSERT_NEAR(mat(2, 2), -0.25, 1e-10, "H(|10>,|10>) should be -1/4");
        TEST_ASSERT_NEAR(mat(3, 3), 0.25, 1e-10, "H(|11>,|11>) should be 1/4");
        
        // Check off-diagonal elements
        TEST_ASSERT_NEAR(mat(1, 2), 0.5, 1e-10, "H(|01>,|10>) should be 1/2");
        TEST_ASSERT_NEAR(mat(2, 1), 0.5, 1e-10, "H(|10>,|01>) should be 1/2");
        
        // Should be Hermitian
        TEST_ASSERT_NEAR(mat(0, 1), mat(1, 0), 1e-10, "H should be Hermitian");
        TEST_ASSERT_NEAR(mat(0, 3), mat(3, 0), 1e-10, "H should be Hermitian");
    }
    
    // Test XXZ Hamiltonian
    {
        auto H = sse::Hamiltonian::xxz(1.0, 0.5);
        const auto& mat = H.getBondMatrix(0);
        
        // Z coupling should be Jz = 0.5
        TEST_ASSERT_NEAR(mat(0, 0), 0.125, 1e-10, "XXZ diagonal should scale with Jz");
        TEST_ASSERT_NEAR(mat(1, 2), 0.5, 1e-10, "XXZ off-diagonal should scale with Jxy");
    }
    
    // Test Ising Hamiltonian
    {
        auto H = sse::Hamiltonian::ising(1.0);
        const auto& mat = H.getBondMatrix(0);
        
        // Should have no off-diagonal elements
        TEST_ASSERT_NEAR(mat(1, 2), 0.0, 1e-10, "Ising should have no off-diagonal");
        TEST_ASSERT_NEAR(mat(2, 1), 0.0, 1e-10, "Ising should have no off-diagonal");
    }
    
    // Test XY Hamiltonian
    {
        auto H = sse::Hamiltonian::xy(1.0);
        const auto& mat = H.getBondMatrix(0);
        
        // Diagonal Sz·Sz should be zero
        TEST_ASSERT_NEAR(mat(0, 0), 0.0, 1e-10, "XY should have zero diagonal");
        TEST_ASSERT_NEAR(mat(3, 3), 0.0, 1e-10, "XY should have zero diagonal");
    }
    
    // Test energy offset
    {
        auto H = sse::Hamiltonian::heisenberg(1.0);
        sse::Real offset = H.getEnergyOffset(0);
        TEST_ASSERT(offset > 0.25, "Energy offset should make all weights positive");
    }
    
    std::cout << "Hamiltonian tests completed.\n\n";
    
    // Test Vertex Data
    std::cout << "Testing VertexData...\n";
    
    {
        auto H = sse::Hamiltonian::heisenberg(1.0);
        sse::VertexData vd(H, 0);
        
        // All diagonal vertices should be allowed
        auto diag_states = vd.getDiagonalStates();
        TEST_ASSERT(diag_states.size() == 4, "Should have 4 diagonal vertex states");
        
        // Weights should be positive
        for (auto vs : vd.getAllowedStates()) {
            TEST_ASSERT(vd.getWeight(vs) > 0, "Vertex weight should be positive");
        }
        
        // Test transition sampling
        sse::VertexState vs = sse::makeVertexState(0, 1, 0, 1);  // |01⟩ → |01⟩
        if (vd.isAllowed(vs)) {
            for (int leg = 0; leg < 4; ++leg) {
                const auto& trans = vd.getTransitions(vs, leg);
                if (!trans.empty()) {
                    TEST_ASSERT_NEAR(trans.back().cumulative_prob, 1.0, 1e-10,
                                    "Transition probabilities should sum to 1");
                }
            }
        }
    }
    
    std::cout << "VertexData tests completed.\n\n";
}
