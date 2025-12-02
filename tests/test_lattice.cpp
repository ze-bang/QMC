/**
 * @file test_lattice.cpp
 * @brief Tests for lattice structures
 */

#include <iostream>
#include <cmath>
#include "sse/lattice.hpp"

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

void test_lattice() {
    std::cout << "Testing Lattice...\n";
    
    // Test 1D chain
    {
        auto lat = sse::Lattice::chain(10);
        TEST_ASSERT(lat.numSites() == 10, "Chain should have 10 sites");
        TEST_ASSERT(lat.numBonds() == 10, "Periodic chain should have 10 bonds");
        TEST_ASSERT(lat.maxCoordination() == 2, "Chain coordination should be 2");
        
        auto neighbors = lat.neighbors(5);
        TEST_ASSERT(neighbors.size() == 2, "Chain site should have 2 neighbors");
    }
    
    // Test open chain
    {
        auto lat = sse::Lattice::chain(10, false);
        TEST_ASSERT(lat.numBonds() == 9, "Open chain should have L-1 bonds");
    }
    
    // Test 2D square lattice
    {
        auto lat = sse::Lattice::square(4, 4);
        TEST_ASSERT(lat.numSites() == 16, "4x4 square should have 16 sites");
        TEST_ASSERT(lat.numBonds() == 32, "4x4 periodic square should have 32 bonds");
        TEST_ASSERT(lat.maxCoordination() == 4, "Square lattice coordination should be 4");
        
        auto neighbors = lat.neighbors(5);
        TEST_ASSERT(neighbors.size() == 4, "Square lattice site should have 4 neighbors");
    }
    
    // Test triangular lattice
    {
        auto lat = sse::Lattice::triangular(4, 4);
        TEST_ASSERT(lat.numSites() == 16, "4x4 triangular should have 16 sites");
        TEST_ASSERT(lat.numBonds() == 48, "4x4 periodic triangular should have 48 bonds");
        TEST_ASSERT(lat.maxCoordination() == 6, "Triangular lattice coordination should be 6");
    }
    
    // Test honeycomb lattice
    {
        auto lat = sse::Lattice::honeycomb(4, 4);
        TEST_ASSERT(lat.numSites() == 32, "4x4 honeycomb should have 32 sites");
        TEST_ASSERT(lat.maxCoordination() == 3, "Honeycomb lattice coordination should be 3");
    }
    
    // Test kagome lattice
    {
        auto lat = sse::Lattice::kagome(4, 4);
        TEST_ASSERT(lat.numSites() == 48, "4x4 kagome should have 48 sites");
        TEST_ASSERT(lat.maxCoordination() == 4, "Kagome lattice coordination should be 4");
    }
    
    // Test cubic lattice
    {
        auto lat = sse::Lattice::cubic(4, 4, 4);
        TEST_ASSERT(lat.numSites() == 64, "4x4x4 cubic should have 64 sites");
        TEST_ASSERT(lat.numBonds() == 192, "4x4x4 periodic cubic should have 192 bonds");
        TEST_ASSERT(lat.maxCoordination() == 6, "Cubic lattice coordination should be 6");
    }
    
    std::cout << "Lattice tests completed.\n\n";
}
