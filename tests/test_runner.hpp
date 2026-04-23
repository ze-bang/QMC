// SPDX-License-Identifier: MIT
//
// A trivial test runner. We deliberately avoid pulling in Catch2 /
// doctest / GoogleTest as a dependency: the test surface here is small
// and a 50-line harness is enough.

#pragma once

#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

namespace qmc::test {

struct Case {
    std::string             name;
    std::function<void()>   fn;
};

inline std::vector<Case>& registry() {
    static std::vector<Case> r;
    return r;
}

inline int& fail_count() {
    static int n = 0;
    return n;
}

struct Register {
    Register(const std::string& name, std::function<void()> fn) {
        registry().push_back({name, std::move(fn)});
    }
};

inline int run_all(const std::string& filter = "") {
    int passed = 0;
    int failed = 0;
    for (const auto& c : registry()) {
        if (!filter.empty() && c.name.find(filter) == std::string::npos) continue;
        try {
            const int before = fail_count();
            c.fn();
            if (fail_count() == before) {
                std::printf("[ PASS ] %s\n", c.name.c_str());
                ++passed;
            } else {
                std::printf("[ FAIL ] %s\n", c.name.c_str());
                ++failed;
            }
        } catch (const std::exception& e) {
            std::printf("[THROW ] %s: %s\n", c.name.c_str(), e.what());
            ++failed;
            ++fail_count();
        } catch (...) {
            std::printf("[THROW ] %s: unknown exception\n", c.name.c_str());
            ++failed;
            ++fail_count();
        }
    }
    std::printf("\n%d passed, %d failed.\n", passed, failed);
    return failed == 0 ? 0 : 1;
}

inline void report_failure(const char* file, int line, const std::string& msg) {
    ++fail_count();
    std::printf("    assertion failed at %s:%d: %s\n", file, line, msg.c_str());
}

} // namespace qmc::test

#define QMC_TEST(name) \
    static void test_##name();                                               \
    static ::qmc::test::Register reg_##name{#name, &test_##name};            \
    static void test_##name()

#define QMC_REQUIRE(cond)                                                    \
    do {                                                                     \
        if (!(cond)) {                                                       \
            ::qmc::test::report_failure(__FILE__, __LINE__, #cond);          \
        }                                                                    \
    } while (0)

#define QMC_REQUIRE_NEAR(a, b, tol)                                          \
    do {                                                                     \
        const double _a = (a);                                               \
        const double _b = (b);                                               \
        const double _t = (tol);                                             \
        if (!(std::fabs(_a - _b) <= _t)) {                                   \
            ::qmc::test::report_failure(__FILE__, __LINE__,                  \
                std::string(#a) + " ~= " + #b + ": got " +                   \
                std::to_string(_a) + " expected " + std::to_string(_b) +     \
                " (tol " + std::to_string(_t) + ")");                        \
        }                                                                    \
    } while (0)
