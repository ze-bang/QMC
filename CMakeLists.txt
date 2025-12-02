cmake_minimum_required(VERSION 3.16)
project(SSE_QMC VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Compiler flags
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra -Wpedantic -fsanitize=address,undefined")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG -ffast-math")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g -march=native")

# Find required packages
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)

# Optional: HDF5 for checkpointing
find_package(HDF5 COMPONENTS CXX)

# fmt library for formatting
include(FetchContent)
FetchContent_Declare(
    fmt
    GIT_REPOSITORY https://github.com/fmtlib/fmt.git
    GIT_TAG 10.1.1
)
FetchContent_MakeAvailable(fmt)

# json library for configuration
FetchContent_Declare(
    json
    GIT_REPOSITORY https://github.com/nlohmann/json.git
    GIT_TAG v3.11.3
)
FetchContent_MakeAvailable(json)

# PCG random number generator (high quality RNG)
FetchContent_Declare(
    pcg
    GIT_REPOSITORY https://github.com/imneme/pcg-cpp.git
    GIT_TAG v0.98.1
)
FetchContent_MakeAvailable(pcg)

# Main library
add_library(sse_qmc_lib STATIC
    src/lattice.cpp
    src/hamiltonian.cpp
    src/vertex.cpp
    src/sse_config.cpp
    src/sse_simulation.cpp
    src/measurements.cpp
    src/random.cpp
)

target_include_directories(sse_qmc_lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${pcg_SOURCE_DIR}/include
)

target_link_libraries(sse_qmc_lib PUBLIC
    Eigen3::Eigen
    OpenMP::OpenMP_CXX
    MPI::MPI_CXX
    fmt::fmt
    nlohmann_json::nlohmann_json
)

if(HDF5_FOUND)
    target_compile_definitions(sse_qmc_lib PUBLIC USE_HDF5)
    target_link_libraries(sse_qmc_lib PUBLIC ${HDF5_CXX_LIBRARIES})
    target_include_directories(sse_qmc_lib PUBLIC ${HDF5_INCLUDE_DIRS})
endif()

# Main executable
add_executable(sse_qmc src/main.cpp)
target_link_libraries(sse_qmc PRIVATE sse_qmc_lib)

# Tests
enable_testing()
add_executable(test_sse
    tests/test_lattice.cpp
    tests/test_hamiltonian.cpp
    tests/test_simulation.cpp
    tests/test_main.cpp
)
target_link_libraries(test_sse PRIVATE sse_qmc_lib)
add_test(NAME SSETests COMMAND test_sse)

# Benchmark executable
add_executable(benchmark_sse benchmarks/benchmark.cpp)
target_link_libraries(benchmark_sse PRIVATE sse_qmc_lib)

# Install
install(TARGETS sse_qmc RUNTIME DESTINATION bin)
install(TARGETS sse_qmc_lib ARCHIVE DESTINATION lib)
install(DIRECTORY include/ DESTINATION include)
