cmake_minimum_required(VERSION 3.20)
project(qmc_sse
    VERSION 0.1.0
    DESCRIPTION "Stochastic Series Expansion QMC for lattice spin systems"
    LANGUAGES CXX
)

# -----------------------------------------------------------------------------
# Build options
# -----------------------------------------------------------------------------
option(QMC_BUILD_TESTS       "Build unit / physics tests"      ON)
option(QMC_BUILD_BENCHMARKS  "Build throughput benchmark"      ON)
option(QMC_ENABLE_OPENMP     "Use OpenMP for replica parallelism"  ON)
option(QMC_ENABLE_LTO        "Enable link-time optimization in Release" ON)
option(QMC_ENABLE_NATIVE     "Tune for the host CPU (-march=native)"    ON)
option(QMC_WERROR            "Treat warnings as errors"        OFF)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
                 Debug Release RelWithDebInfo MinSizeRel)
endif()

# -----------------------------------------------------------------------------
# Header-only library target
# -----------------------------------------------------------------------------
add_library(qmc_sse_lib INTERFACE)
target_include_directories(qmc_sse_lib INTERFACE
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
)
target_compile_features(qmc_sse_lib INTERFACE cxx_std_20)

# Warning flags (gcc / clang / msvc).
if(MSVC)
    target_compile_options(qmc_sse_lib INTERFACE /W4 /permissive-)
    if(QMC_WERROR)
        target_compile_options(qmc_sse_lib INTERFACE /WX)
    endif()
else()
    target_compile_options(qmc_sse_lib INTERFACE
        -Wall -Wextra -Wpedantic -Wshadow -Wnon-virtual-dtor
        -Wcast-align -Wunused -Woverloaded-virtual
        -Wdouble-promotion -Wformat=2
    )
    if(QMC_WERROR)
        target_compile_options(qmc_sse_lib INTERFACE -Werror)
    endif()
    if(QMC_ENABLE_NATIVE)
        target_compile_options(qmc_sse_lib INTERFACE
            $<$<CONFIG:Release,RelWithDebInfo>:-march=native>)
    endif()
endif()

# OpenMP --------------------------------------------------------------------
if(QMC_ENABLE_OPENMP)
    find_package(OpenMP COMPONENTS CXX)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(qmc_sse_lib INTERFACE OpenMP::OpenMP_CXX)
        message(STATUS "qmc_sse: OpenMP enabled (${OpenMP_CXX_VERSION})")
    else()
        message(STATUS "qmc_sse: OpenMP requested but not found, building serial")
    endif()
endif()

# IPO / LTO -----------------------------------------------------------------
if(QMC_ENABLE_LTO)
    include(CheckIPOSupported)
    check_ipo_supported(RESULT _ipo_ok OUTPUT _ipo_err)
    if(_ipo_ok)
        set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
        set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELWITHDEBINFO ON)
    else()
        message(STATUS "qmc_sse: LTO requested but not supported: ${_ipo_err}")
    endif()
endif()

# -----------------------------------------------------------------------------
# CLI driver
# -----------------------------------------------------------------------------
add_executable(qmc_sse apps/qmc_sse.cpp)
target_link_libraries(qmc_sse PRIVATE qmc_sse_lib)

# -----------------------------------------------------------------------------
# Tests
# -----------------------------------------------------------------------------
if(QMC_BUILD_TESTS)
    enable_testing()
    file(GLOB QMC_TEST_SOURCES CONFIGURE_DEPENDS "tests/test_*.cpp")
    add_executable(qmc_tests ${QMC_TEST_SOURCES})
    target_link_libraries(qmc_tests PRIVATE qmc_sse_lib)
    target_include_directories(qmc_tests PRIVATE tests)
    add_test(NAME qmc_tests COMMAND qmc_tests)
    set_tests_properties(qmc_tests PROPERTIES TIMEOUT 600)
endif()

# -----------------------------------------------------------------------------
# Benchmarks
# -----------------------------------------------------------------------------
if(QMC_BUILD_BENCHMARKS)
    add_executable(qmc_bench benchmarks/bench_sse.cpp)
    target_link_libraries(qmc_bench PRIVATE qmc_sse_lib)
endif()

# -----------------------------------------------------------------------------
# Install / export
# -----------------------------------------------------------------------------
include(GNUInstallDirs)
install(TARGETS qmc_sse RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Friendly final message.
message(STATUS "qmc_sse configured: build=${CMAKE_BUILD_TYPE} "
               "tests=${QMC_BUILD_TESTS} bench=${QMC_BUILD_BENCHMARKS}")
