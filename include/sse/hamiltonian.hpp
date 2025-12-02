#pragma once

/**
 * @file hamiltonian.hpp
 * @brief Hamiltonian definitions for spin-1/2 systems
 * 
 * Supports arbitrary nearest-neighbor spin interaction matrices.
 * The Hamiltonian is specified by a 4x4 matrix acting on bond (i,j).
 */

#include <vector>
#include <array>
#include <memory>
#include <functional>
#include <Eigen/Dense>
#include "sse/types.hpp"
#include "sse/lattice.hpp"

namespace sse {

/**
 * @brief Pauli matrices and common spin operators
 */
namespace SpinOp {
    // Pauli matrices (σ/2 convention for spin-1/2)
    inline SpinMatrix Sx() {
        SpinMatrix m;
        m << 0.0, 0.5,
             0.5, 0.0;
        return m;
    }
    
    inline SpinMatrix Sy() {
        // Note: Sy is imaginary, we store real part only for real Hamiltonians
        // For complex case, use ComplexMatrix
        SpinMatrix m;
        m << 0.0, 0.0,
             0.0, 0.0;
        return m;
    }
    
    inline SpinMatrix iSy() {
        // i*Sy for use in real Hamiltonians
        SpinMatrix m;
        m << 0.0, 0.5,
            -0.5, 0.0;
        return m;
    }
    
    inline SpinMatrix Sz() {
        SpinMatrix m;
        m << 0.5, 0.0,
             0.0, -0.5;
        return m;
    }
    
    inline SpinMatrix Sp() {
        // S+ = Sx + i*Sy
        SpinMatrix m;
        m << 0.0, 1.0,
             0.0, 0.0;
        return m;
    }
    
    inline SpinMatrix Sm() {
        // S- = Sx - i*Sy
        SpinMatrix m;
        m << 0.0, 0.0,
             1.0, 0.0;
        return m;
    }
    
    inline SpinMatrix Id() {
        return SpinMatrix::Identity();
    }
}

/**
 * @brief Kronecker product of two matrices
 */
template<typename Derived1, typename Derived2>
auto kronecker(const Eigen::MatrixBase<Derived1>& A, 
               const Eigen::MatrixBase<Derived2>& B) {
    using Scalar = typename Derived1::Scalar;
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> result(
        A.rows() * B.rows(), A.cols() * B.cols());
    
    for (int i = 0; i < A.rows(); ++i) {
        for (int j = 0; j < A.cols(); ++j) {
            result.block(i * B.rows(), j * B.cols(), B.rows(), B.cols()) = 
                A(i, j) * B;
        }
    }
    return result;
}

/**
 * @brief Hamiltonian class for spin-1/2 systems
 * 
 * The bond Hamiltonian is a 4x4 matrix acting on the 4-dimensional
 * Hilbert space of two spin-1/2 sites: {|00⟩, |01⟩, |10⟩, |11⟩}
 * where |0⟩ = |↓⟩ and |1⟩ = |↑⟩.
 */
class Hamiltonian {
public:
    /**
     * @brief Default constructor
     */
    Hamiltonian() = default;
    
    /**
     * @brief Construct Hamiltonian with uniform coupling
     * @param bond_matrix 4x4 bond Hamiltonian matrix
     */
    explicit Hamiltonian(const BondMatrix& bond_matrix);
    
    /**
     * @brief Construct Hamiltonian with different bond types
     * @param bond_matrices Map from bond type to Hamiltonian matrix
     */
    explicit Hamiltonian(const std::vector<BondMatrix>& bond_matrices);
    
    /**
     * @brief Create Heisenberg Hamiltonian: H = J * S_i · S_j
     */
    static Hamiltonian heisenberg(Real J = 1.0);
    
    /**
     * @brief Create XXZ Hamiltonian: H = Jxy*(Sx·Sx + Sy·Sy) + Jz*Sz·Sz
     */
    static Hamiltonian xxz(Real Jxy = 1.0, Real Jz = 1.0);
    
    /**
     * @brief Create XY Hamiltonian: H = J*(Sx·Sx + Sy·Sy)
     */
    static Hamiltonian xy(Real J = 1.0);
    
    /**
     * @brief Create Ising Hamiltonian: H = J*Sz·Sz
     */
    static Hamiltonian ising(Real J = 1.0);
    
    /**
     * @brief Create general XXZ with field: H = Jxy*(Sx·Sx + Sy·Sy) + Jz*Sz·Sz - h*(Sz_i + Sz_j)
     */
    static Hamiltonian xxzWithField(Real Jxy, Real Jz, Real h);
    
    /**
     * @brief Create fully anisotropic XYZ: H = Jx*Sx·Sx + Jy*Sy·Sy + Jz*Sz·Sz
     */
    static Hamiltonian xyz(Real Jx, Real Jy, Real Jz);
    
    /**
     * @brief Create custom Hamiltonian from 4x4 matrix
     */
    static Hamiltonian custom(const BondMatrix& matrix);
    
    /**
     * @brief Create Hamiltonian with DM interaction
     * H = J*S_i·S_j + D·(S_i × S_j)
     */
    static Hamiltonian heisenbergDM(Real J, Real Dx, Real Dy, Real Dz);
    
    /**
     * @brief Get bond Hamiltonian matrix for given bond type
     */
    const BondMatrix& getBondMatrix(int bond_type = 0) const;
    
    /**
     * @brief Get all bond matrices
     */
    const std::vector<BondMatrix>& getBondMatrices() const { return bond_matrices_; }
    
    /**
     * @brief Get number of bond types
     */
    int numBondTypes() const { return static_cast<int>(bond_matrices_.size()); }
    
    /**
     * @brief Get matrix element ⟨out|H|in⟩
     */
    Real matrixElement(VertexState in_state, VertexState out_state, int bond_type = 0) const;
    
    /**
     * @brief Check if transition is allowed (non-zero matrix element)
     */
    bool isAllowedTransition(VertexState in_state, VertexState out_state, int bond_type = 0) const;
    
    /**
     * @brief Get the energy offset needed for positive weights
     */
    Real getEnergyOffset(int bond_type = 0) const;
    
    /**
     * @brief Print Hamiltonian info
     */
    void print() const;
    
private:
    std::vector<BondMatrix> bond_matrices_;
    std::vector<Real> energy_offsets_;
    
    void computeEnergyOffsets();
    static constexpr Real TOLERANCE = 1e-10;
};

} // namespace sse
