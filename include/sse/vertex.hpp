#pragma once

/**
 * @file vertex.hpp
 * @brief Vertex data structures and operations for SSE
 * 
 * A vertex in SSE represents an interaction event at a bond.
 * The vertex has 4 legs (2 sites × 2 time directions).
 * This file handles vertex weights and loop update probabilities.
 */

#include <vector>
#include <array>
#include <tuple>
#include "sse/types.hpp"
#include "sse/hamiltonian.hpp"

namespace sse {

/**
 * @brief Precomputed vertex data for efficient updates
 * 
 * Stores weights and transition probabilities for all vertex configurations.
 * For spin-1/2 with 4 legs, there are 16 possible vertex states.
 */
class VertexData {
public:
    static constexpr int NUM_VERTEX_STATES = 16;  // 2^4 for spin-1/2
    
    /**
     * @brief Construct vertex data from Hamiltonian and bond type
     */
    VertexData(const Hamiltonian& H, int bond_type = 0);
    
    /**
     * @brief Default constructor
     */
    VertexData() = default;
    
    /**
     * @brief Get vertex weight (always non-negative after shift)
     */
    Real getWeight(VertexState vs) const { return weights_[vs]; }
    
    /**
     * @brief Check if vertex state is allowed (non-zero weight)
     */
    bool isAllowed(VertexState vs) const { return allowed_[vs]; }
    
    /**
     * @brief Check if vertex is diagonal
     */
    bool isDiagonal(VertexState vs) const { return isDiagonalVertex(vs); }
    
    /**
     * @brief Get diagonal vertex state for given input spins
     */
    VertexState getDiagonalVertex(StateIdx si, StateIdx sj) const {
        return makeVertexState(si, sj, si, sj);
    }
    
    /**
     * @brief Get all allowed vertex states
     */
    const std::vector<VertexState>& getAllowedStates() const { return allowed_states_; }
    
    /**
     * @brief Get diagonal vertex states
     */
    const std::vector<VertexState>& getDiagonalStates() const { return diagonal_states_; }
    
    /**
     * @brief Get energy offset
     */
    Real getEnergyOffset() const { return energy_offset_; }
    
    /**
     * @brief Get total weight (sum of all diagonal weights for normalization)
     */
    Real getTotalDiagonalWeight() const { return total_diagonal_weight_; }
    
    // Loop update transition probabilities
    
    /**
     * @brief Transition probability for directed loop update
     * 
     * Returns (exit_leg, new_vertex_state, probability) for loop entering at entry_leg.
     * Uses heat-bath algorithm for detailed balance.
     */
    struct Transition {
        int exit_leg;
        VertexState new_state;
        Real cumulative_prob;
    };
    
    /**
     * @brief Get transitions for loop update entering at given leg
     * @param vs Current vertex state
     * @param entry_leg Leg where loop enters (0-3)
     * @return Vector of possible transitions with cumulative probabilities
     */
    const std::vector<Transition>& getTransitions(VertexState vs, int entry_leg) const {
        return transitions_[vs * NUM_LEGS + entry_leg];
    }
    
    /**
     * @brief Sample a transition using directed loop algorithm
     * @param vs Current vertex state
     * @param entry_leg Entry leg
     * @param random Uniform random number in [0,1)
     * @return Tuple of (exit_leg, new_vertex_state)
     */
    std::tuple<int, VertexState> sampleTransition(VertexState vs, int entry_leg, Real random) const;
    
    /**
     * @brief Check if bounce (entering and exiting same leg) is possible
     */
    bool canBounce(VertexState vs, int leg) const;
    
    /**
     * @brief Print vertex data info
     */
    void print() const;
    
private:
    std::array<Real, NUM_VERTEX_STATES> weights_;
    std::array<bool, NUM_VERTEX_STATES> allowed_;
    std::vector<VertexState> allowed_states_;
    std::vector<VertexState> diagonal_states_;
    Real energy_offset_ = 0.0;
    Real total_diagonal_weight_ = 0.0;
    
    // Transition table for loop updates: [vertex_state * 4 + entry_leg]
    std::array<std::vector<Transition>, NUM_VERTEX_STATES * NUM_LEGS> transitions_;
    
    void computeWeights(const Hamiltonian& H, int bond_type);
    void computeTransitions();
    
    static constexpr Real TOLERANCE = 1e-12;
};

/**
 * @brief Collection of vertex data for all bond types
 */
class VertexDataCollection {
public:
    VertexDataCollection() = default;
    
    /**
     * @brief Construct from Hamiltonian
     */
    explicit VertexDataCollection(const Hamiltonian& H);
    
    /**
     * @brief Get vertex data for bond type
     */
    const VertexData& get(int bond_type) const { return data_[bond_type]; }
    
    /**
     * @brief Number of bond types
     */
    int numBondTypes() const { return static_cast<int>(data_.size()); }
    
    /**
     * @brief Get maximum energy offset across all bond types
     */
    Real maxEnergyOffset() const { return max_energy_offset_; }
    
private:
    std::vector<VertexData> data_;
    Real max_energy_offset_ = 0.0;
};

} // namespace sse
