/**
 * @file test_main.cpp
 * @brief Test runner using minimal test framework
 */

#include <iostream>
#include <cstdlib>

// Simple test framework
int g_test_count = 0;
int g_test_passed = 0;
int g_test_failed = 0;

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

// Declare test functions
void test_lattice();
void test_hamiltonian();
void test_simulation();

int main() {
    std::cout << "=== SSE QMC Tests ===\n\n";
    
    test_lattice();
    test_hamiltonian();
    test_simulation();
    
    std::cout << "\n=== Test Summary ===\n";
    std::cout << "Total: " << g_test_count << "\n";
    std::cout << "Passed: " << g_test_passed << "\n";
    std::cout << "Failed: " << g_test_failed << "\n";
    
    return g_test_failed > 0 ? 1 : 0;
}
