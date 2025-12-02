#pragma once

/**
 * @file types.hpp
 * @brief Common type definitions for SSE QMC simulation
 * 
 * This header defines fundamental types used throughout the SSE simulation.
 */

#include <cstdint>
#include <complex>
#include <vector>
#include <array>
#include <Eigen/Dense>

namespace sse {

// Basic types for efficiency
using StateIdx = uint8_t;       // State index (0 or 1 for spin-1/2)
using BondIdx = uint32_t;       // Bond index
using SiteIdx = uint32_t;       // Site index
using OperIdx = uint64_t;       // Operator code

// Real and complex types
using Real = double;
using Complex = std::complex<double>;

// Matrix types using Eigen
using RealMatrix = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using ComplexMatrix = Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using RealVector = Eigen::Matrix<Real, Eigen::Dynamic, 1>;
using ComplexVector = Eigen::Matrix<Complex, Eigen::Dynamic, 1>;

// Fixed size for spin-1/2 (2x2 matrices)
using SpinMatrix = Eigen::Matrix<Real, 2, 2, Eigen::RowMajor>;
using SpinVector = Eigen::Matrix<Real, 2, 1>;

// 4x4 bond Hamiltonian for two spin-1/2
using BondMatrix = Eigen::Matrix<Real, 4, 4, Eigen::RowMajor>;

/**
 * @brief Spin state enumeration for spin-1/2
 */
enum class Spin : StateIdx {
    Down = 0,   // |↓⟩ or |0⟩
    Up = 1      // |↑⟩ or |1⟩
};

/**
 * @brief Operator type enumeration
 */
enum class OperatorType : uint8_t {
    Identity = 0,       // Identity (null) operator
    Diagonal = 1,       // Diagonal operator
    OffDiagonal = 2     // Off-diagonal operator
};

/**
 * @brief Vertex leg enumeration for 4-leg vertex (2 sites, in/out)
 * 
 * Leg layout for bond (i,j):
 *   leg 0: site i, input (bottom)
 *   leg 1: site j, input (bottom)  
 *   leg 2: site i, output (top)
 *   leg 3: site j, output (top)
 */
constexpr int NUM_LEGS = 4;

/**
 * @brief Vertex state: encodes the 4 spin states at vertex legs
 * 
 * Bits: [s_i_in, s_j_in, s_i_out, s_j_out]
 * For spin-1/2: 16 possible vertex states (2^4)
 */
using VertexState = uint8_t;

/**
 * @brief Compute vertex state from individual spin states
 */
inline VertexState makeVertexState(StateIdx si_in, StateIdx sj_in, 
                                    StateIdx si_out, StateIdx sj_out) {
    return static_cast<VertexState>((si_in) | (sj_in << 1) | 
                                     (si_out << 2) | (sj_out << 3));
}

/**
 * @brief Extract spin at leg from vertex state
 */
inline StateIdx getSpinAtLeg(VertexState vs, int leg) {
    return (vs >> leg) & 1;
}

/**
 * @brief Flip spin at leg in vertex state
 */
inline VertexState flipSpinAtLeg(VertexState vs, int leg) {
    return vs ^ (1 << leg);
}

/**
 * @brief Check if vertex is diagonal (input state == output state)
 */
inline bool isDiagonalVertex(VertexState vs) {
    return (vs & 0x3) == ((vs >> 2) & 0x3);
}

/**
 * @brief Bond structure for storing lattice connectivity
 */
struct Bond {
    SiteIdx i;          // First site
    SiteIdx j;          // Second site
    int type;           // Bond type (for different couplings)
    
    Bond() = default;
    Bond(SiteIdx i_, SiteIdx j_, int type_ = 0) 
        : i(i_), j(j_), type(type_) {}
};

/**
 * @brief Operator code encoding bond index and vertex state
 * 
 * Layout: [bond_index (upper bits) | vertex_state (lower 4 bits) | type (1 bit)]
 * Identity operator: code == 0
 */
class OperatorCode {
public:
    static constexpr OperIdx IDENTITY = 0;
    static constexpr int VERTEX_BITS = 4;
    static constexpr int TYPE_BITS = 1;
    
    OperatorCode() : code_(IDENTITY) {}
    
    static OperatorCode identity() { return OperatorCode(); }
    
    static OperatorCode diagonal(BondIdx bond, VertexState vs) {
        return OperatorCode(1 | (static_cast<OperIdx>(vs) << TYPE_BITS) | 
                           (static_cast<OperIdx>(bond) << (TYPE_BITS + VERTEX_BITS)));
    }
    
    static OperatorCode offDiagonal(BondIdx bond, VertexState vs) {
        return OperatorCode(0 | (static_cast<OperIdx>(vs) << TYPE_BITS) | 
                           (static_cast<OperIdx>(bond) << (TYPE_BITS + VERTEX_BITS)));
    }
    
    bool isIdentity() const { return code_ == IDENTITY; }
    bool isDiagonal() const { return !isIdentity() && (code_ & 1); }
    bool isOffDiagonal() const { return !isIdentity() && !(code_ & 1); }
    
    BondIdx bond() const { 
        return static_cast<BondIdx>(code_ >> (TYPE_BITS + VERTEX_BITS)); 
    }
    
    VertexState vertexState() const { 
        return static_cast<VertexState>((code_ >> TYPE_BITS) & 0xF); 
    }
    
    void setVertexState(VertexState vs) {
        code_ = (code_ & ~(static_cast<OperIdx>(0xF) << TYPE_BITS)) | 
                (static_cast<OperIdx>(vs) << TYPE_BITS);
    }
    
    OperIdx code() const { return code_; }
    
private:
    explicit OperatorCode(OperIdx code) : code_(code) {}
    OperIdx code_;
};

/**
 * @brief Simulation parameters
 */
struct SimulationParams {
    Real beta = 1.0;            // Inverse temperature
    int n_therm = 10000;        // Thermalization sweeps
    int n_sweeps = 100000;      // Measurement sweeps
    int n_bins = 100;           // Number of bins for error analysis
    int measure_every = 1;      // Measure every N sweeps
    int checkpoint_every = 1000;// Checkpoint frequency
    uint64_t seed = 42;         // Random seed
    bool use_loop_update = true;// Use loop (cluster) update
    int max_expansion_order = 0;// Max expansion order (0 = auto)
};

} // namespace sse
