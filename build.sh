#!/bin/bash

# Build script for SSE QMC

set -e  # Exit on any error

echo "=== Building SSE Quantum Monte Carlo ==="

# Check if CMake is available
if command -v cmake &> /dev/null; then
    echo "Using CMake build system..."
    
    # Create build directory
    mkdir -p build
    cd build
    
    # Configure
    echo "Configuring..."
    cmake ..
    
    # Build
    echo "Building..."
    make -j$(nproc)
    
    echo "Build completed successfully!"
    echo "Executables are in the build/ directory:"
    ls -la sse_qmc *_heisenberg *_scan *_comparison 2>/dev/null || true
    
elif command -v make &> /dev/null; then
    echo "Using Makefile build system..."
    
    # Build using Makefile
    make clean
    make -j$(nproc)
    
    echo "Build completed successfully!"
    echo "Executable: sse_qmc"
    
else
    echo "Error: Neither CMake nor Make found!"
    echo "Please install build tools:"
    echo "  Ubuntu/Debian: sudo apt install build-essential cmake"
    echo "  CentOS/RHEL: sudo yum install gcc-c++ make cmake"
    exit 1
fi

echo ""
echo "=== Quick Test ==="
if [ -f "sse_qmc" ]; then
    echo "Running quick test..."
    ./sse_qmc --help
elif [ -f "build/sse_qmc" ]; then
    echo "Running quick test..."
    cd build
    ./sse_qmc --help
else
    echo "Warning: Executable not found!"
fi

echo ""
echo "=== Usage Examples ==="
echo "Run the main program:"
echo "  ./sse_qmc --lattice square --model heisenberg --Lx 8 --Ly 8 --beta 10"
echo ""
echo "Run example programs:"
echo "  ./square_heisenberg"
echo "  ./temperature_scan"
echo "  ./lattice_comparison"
echo "  ./model_comparison"
echo ""
echo "Build completed successfully!"
