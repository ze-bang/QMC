// SPDX-License-Identifier: MIT
//
// PCG32 (O'Neill 2014, "PCG: A Family of Simple Fast Space-Efficient
// Statistically Good Algorithms for Random Number Generation"). A small,
// header-only, fast pseudo-random number generator with a long period
// (2^64) and excellent statistical properties for Monte Carlo work.
//
// We avoid std::mt19937 because:
//   * it is large (~2.5 KiB state) which is bad for parallel replicas,
//   * it has known statistical issues on some standard QMC tests, and
//   * its seeding from a single 32-bit value is famously poor.
//
// The implementation here is deliberately minimal and only exposes the
// surface that the rest of the library needs.

#pragma once

#include <array>
#include <chrono>
#include <cstdint>
#include <limits>
#include <thread>

#include "qmc/types.hpp"

namespace qmc {

class Pcg32 {
public:
    using result_type = std::uint32_t;

    static constexpr result_type min() noexcept { return 0u; }
    static constexpr result_type max() noexcept {
        return std::numeric_limits<result_type>::max();
    }

    // Default-constructed engines are seeded from a high-resolution clock
    // mixed with the calling thread id. Use an explicit seed for
    // reproducible runs.
    Pcg32() noexcept { seed(default_seed_()); }
    explicit Pcg32(std::uint64_t seed_value, std::uint64_t stream = 0xda3e39cb94b95bdbULL) noexcept {
        seed(seed_value, stream);
    }

    void seed(std::uint64_t seed_value,
              std::uint64_t stream = 0xda3e39cb94b95bdbULL) noexcept {
        state_ = 0u;
        inc_   = (stream << 1u) | 1u;
        next_();
        state_ += seed_value;
        next_();
    }

    // Raw 32-bit uniform draw.
    result_type operator()() noexcept { return next_(); }

    // Uniform real in [0, 1). Uses the 24 highest bits to fit a float
    // mantissa exactly, then promotes; sufficient for QMC weights.
    Real uniform() noexcept {
        // 53-bit double from two 32-bit draws, full mantissa precision.
        const std::uint64_t a = next_() >> 5;       // 27 bits
        const std::uint64_t b = next_() >> 6;       // 26 bits
        return (a * (1ULL << 26) + b) * (1.0 / (1ULL << 53));
    }

    // Uniform integer in [0, n). Lemire (2019) bounded fast method.
    std::uint32_t uniform_int(std::uint32_t n) noexcept {
        std::uint64_t m = static_cast<std::uint64_t>(next_()) * static_cast<std::uint64_t>(n);
        std::uint32_t l = static_cast<std::uint32_t>(m);
        if (l < n) {
            const std::uint32_t t = static_cast<std::uint32_t>(-n) % n;
            while (l < t) {
                m = static_cast<std::uint64_t>(next_()) * static_cast<std::uint64_t>(n);
                l = static_cast<std::uint32_t>(m);
            }
        }
        return static_cast<std::uint32_t>(m >> 32);
    }

    bool bernoulli(Real p) noexcept { return uniform() < p; }

private:
    std::uint64_t state_ = 0;
    std::uint64_t inc_   = 0;

    [[gnu::always_inline]]
    inline std::uint32_t next_() noexcept {
        const std::uint64_t old = state_;
        state_ = old * 6364136223846793005ULL + inc_;
        const std::uint32_t xorshifted =
            static_cast<std::uint32_t>(((old >> 18u) ^ old) >> 27u);
        const std::uint32_t rot = static_cast<std::uint32_t>(old >> 59u);
        return (xorshifted >> rot) | (xorshifted << ((-static_cast<int>(rot)) & 31));
    }

    static std::uint64_t default_seed_() noexcept {
        const auto t = std::chrono::high_resolution_clock::now()
                           .time_since_epoch()
                           .count();
        const std::size_t tid =
            std::hash<std::thread::id>{}(std::this_thread::get_id());
        // SplitMix64 mix step.
        std::uint64_t z = static_cast<std::uint64_t>(t) ^
                          (static_cast<std::uint64_t>(tid) * 0x9E3779B97F4A7C15ULL);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        return z ^ (z >> 31);
    }
};

} // namespace qmc
