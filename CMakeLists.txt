cmake_minimum_required(VERSION 3.12)
project(SSE_QMC VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find required packages
find_package(OpenMP)

# Include directories
include_directories(include)

# Source files (exclude examples)
file(GLOB SOURCES "src/*.cpp")
file(GLOB HEADERS "include/*.h")

# Create the main executable
add_executable(sse_qmc ${SOURCES} ${HEADERS})

# Link libraries
if(OpenMP_CXX_FOUND)
    target_link_libraries(sse_qmc OpenMP::OpenMP_CXX)
endif()

# Compiler flags
target_compile_options(sse_qmc PRIVATE -O3 -march=native -Wall -Wextra)

# Enable debug symbols in debug mode
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(sse_qmc PRIVATE -g -O0)
endif()

# Example executables
file(GLOB EXAMPLE_SOURCES "examples/*.cpp")
foreach(EXAMPLE_SOURCE ${EXAMPLE_SOURCES})
    get_filename_component(EXAMPLE_NAME ${EXAMPLE_SOURCE} NAME_WE)
    add_executable(${EXAMPLE_NAME} ${EXAMPLE_SOURCE} 
                   src/sse_qmc.cpp src/lattice.cpp src/hamiltonian.cpp src/observables.cpp)
    target_include_directories(${EXAMPLE_NAME} PRIVATE include)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(${EXAMPLE_NAME} OpenMP::OpenMP_CXX)
    endif()
    target_compile_options(${EXAMPLE_NAME} PRIVATE -O3 -march=native -Wall -Wextra)
endforeach()
