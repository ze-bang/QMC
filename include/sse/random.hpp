#pragma once

/**
 * @file random.hpp
 * @brief High-quality random number generation for SSE QMC
 * 
 * Uses PCG (Permuted Congruential Generator) for high quality
 * statistically robust random numbers.
 */

#include <pcg_random.hpp>
#include <random>
#include <array>
#include "sse/types.hpp"

namespace sse {

/**
 * @brief Thread-safe random number generator wrapper
 * 
 * Each thread should have its own instance to avoid synchronization overhead.
 */
class Random {
public:
    using generator_type = pcg64;
    
    /**
     * @brief Construct with seed
     */
    explicit Random(uint64_t seed = 42) 
        : rng_(seed), uniform_(0.0, 1.0), uniform_int_() {}
    
    /**
     * @brief Seed with multiple values for better entropy
     */
    void seed(uint64_t s1, uint64_t s2 = 0) {
        rng_.seed(s1, s2);
    }
    
    /**
     * @brief Generate uniform random double in [0, 1)
     */
    Real uniform01() {
        return uniform_(rng_);
    }
    
    /**
     * @brief Generate uniform random integer in [0, n)
     */
    template<typename T>
    T uniformInt(T n) {
        return std::uniform_int_distribution<T>(0, n - 1)(rng_);
    }
    
    /**
     * @brief Generate random integer in [a, b]
     */
    template<typename T>
    T uniformInt(T a, T b) {
        return std::uniform_int_distribution<T>(a, b)(rng_);
    }
    
    /**
     * @brief Generate random boolean
     */
    bool randomBool() {
        return uniform01() < 0.5;
    }
    
    /**
     * @brief Get raw generator for custom distributions
     */
    generator_type& generator() { return rng_; }
    const generator_type& generator() const { return rng_; }
    
    /**
     * @brief Generate random spin state (0 or 1)
     */
    StateIdx randomSpin() {
        return uniformInt<StateIdx>(2);
    }
    
private:
    generator_type rng_;
    std::uniform_real_distribution<Real> uniform_;
    std::uniform_int_distribution<uint64_t> uniform_int_;
};

/**
 * @brief Create a seeded random generator for parallel simulations
 * 
 * @param base_seed Base seed
 * @param thread_id Thread/rank identifier
 * @return Uniquely seeded Random object
 */
inline Random createParallelRandom(uint64_t base_seed, int thread_id) {
    // Use splitmix64 to generate independent seeds
    auto splitmix = [](uint64_t x) {
        x += 0x9e3779b97f4a7c15;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9;
        x = (x ^ (x >> 27)) * 0x94d049bb133111eb;
        return x ^ (x >> 31);
    };
    
    Random rng;
    rng.seed(splitmix(base_seed + thread_id), 
             splitmix(base_seed * 31 + thread_id * 37));
    return rng;
}

} // namespace sse
