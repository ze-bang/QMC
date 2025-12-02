/**
 * @file vertex.cpp
 * @brief Implementation of vertex data structures
 */

#include "sse/vertex.hpp"
#include <fmt/format.h>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <cmath>

namespace sse {

VertexData::VertexData(const Hamiltonian& H, int bond_type) {
    computeWeights(H, bond_type);
    computeTransitions();
}

void VertexData::computeWeights(const Hamiltonian& H, int bond_type) {
    const BondMatrix& Hmat = H.getBondMatrix(bond_type);
    
    // Compute energy offset to make all weights non-negative
    // Weight = -H_ij + δ_ij * energy_offset
    energy_offset_ = H.getEnergyOffset(bond_type);
    
    allowed_states_.clear();
    diagonal_states_.clear();
    total_diagonal_weight_ = 0.0;
    
    for (int vs = 0; vs < NUM_VERTEX_STATES; ++vs) {
        VertexState vertex_state = static_cast<VertexState>(vs);
        
        // Input state: (s_i_in, s_j_in) -> index in {0,1,2,3}
        int in_idx = (getSpinAtLeg(vertex_state, 0)) | 
                     (getSpinAtLeg(vertex_state, 1) << 1);
        // Output state: (s_i_out, s_j_out)
        int out_idx = (getSpinAtLeg(vertex_state, 2)) | 
                      (getSpinAtLeg(vertex_state, 3) << 1);
        
        // Weight = -⟨out|H|in⟩ for off-diagonal
        // Weight = -⟨in|H|in⟩ + offset for diagonal
        Real w = -Hmat(out_idx, in_idx);
        if (in_idx == out_idx) {
            w += energy_offset_;
        }
        
        if (w > TOLERANCE) {
            weights_[vs] = w;
            allowed_[vs] = true;
            allowed_states_.push_back(vertex_state);
            
            if (isDiagonalVertex(vertex_state)) {
                diagonal_states_.push_back(vertex_state);
                total_diagonal_weight_ += w;
            }
        } else {
            weights_[vs] = 0.0;
            allowed_[vs] = false;
        }
    }
}

void VertexData::computeTransitions() {
    // For each (vertex_state, entry_leg), compute transition probabilities
    // using heat-bath algorithm
    //
    // When a loop enters at leg l, it can:
    // 1. Exit through any of the 4 legs (including bouncing back)
    // 2. The vertex state may change if leg states flip
    //
    // For spin-1/2, the loop can flip the spin at entry leg.
    // Conservation: if entry leg flips, one other leg must also flip
    // to maintain a valid vertex.
    
    for (int vs = 0; vs < NUM_VERTEX_STATES; ++vs) {
        VertexState vertex_state = static_cast<VertexState>(vs);
        
        if (!allowed_[vs]) {
            // No transitions for disallowed vertices
            for (int leg = 0; leg < NUM_LEGS; ++leg) {
                transitions_[vs * NUM_LEGS + leg].clear();
            }
            continue;
        }
        
        for (int entry_leg = 0; entry_leg < NUM_LEGS; ++entry_leg) {
            std::vector<Transition>& trans = transitions_[vs * NUM_LEGS + entry_leg];
            trans.clear();
            
            // Consider all possible exit legs
            Real total_weight = 0.0;
            std::vector<std::pair<int, VertexState>> candidates;
            std::vector<Real> candidate_weights;
            
            for (int exit_leg = 0; exit_leg < NUM_LEGS; ++exit_leg) {
                // New vertex state after flipping entry and exit leg spins
                VertexState new_state = flipSpinAtLeg(vertex_state, entry_leg);
                if (exit_leg != entry_leg) {
                    new_state = flipSpinAtLeg(new_state, exit_leg);
                }
                // For bouncing (entry == exit), flip only once (already done)
                if (exit_leg == entry_leg) {
                    new_state = flipSpinAtLeg(vertex_state, entry_leg);
                }
                
                if (allowed_[new_state]) {
                    Real w = weights_[new_state];
                    candidates.emplace_back(exit_leg, new_state);
                    candidate_weights.push_back(w);
                    total_weight += w;
                }
            }
            
            // Also consider straight-through (no flip at entry)
            // This corresponds to the "directed loop" update where
            // the worm can traverse without flipping
            for (int exit_leg = 0; exit_leg < NUM_LEGS; ++exit_leg) {
                if (exit_leg == entry_leg) continue; // Already considered
                
                // Check if we can exit without flipping entry
                // This requires matching states on entry/exit legs
                if (getSpinAtLeg(vertex_state, entry_leg) == 
                    getSpinAtLeg(vertex_state, exit_leg)) {
                    // Can traverse through
                    candidates.emplace_back(exit_leg, vertex_state);
                    Real w = weights_[vs];
                    candidate_weights.push_back(w);
                    total_weight += w;
                }
            }
            
            // Build cumulative probability distribution
            Real cumsum = 0.0;
            for (size_t i = 0; i < candidates.size(); ++i) {
                cumsum += candidate_weights[i] / total_weight;
                trans.push_back({candidates[i].first, candidates[i].second, cumsum});
            }
            
            // Ensure last probability is exactly 1.0
            if (!trans.empty()) {
                trans.back().cumulative_prob = 1.0;
            }
        }
    }
}

std::tuple<int, VertexState> VertexData::sampleTransition(
    VertexState vs, int entry_leg, Real random) const {
    
    const auto& trans = transitions_[vs * NUM_LEGS + entry_leg];
    
    if (trans.empty()) {
        // No valid transitions; shouldn't happen for valid configurations
        return {entry_leg, vs};  // Bounce back unchanged
    }
    
    // Binary search for transition
    for (const auto& t : trans) {
        if (random < t.cumulative_prob) {
            return {t.exit_leg, t.new_state};
        }
    }
    
    // Fallback to last transition (shouldn't reach here due to cumsum = 1)
    return {trans.back().exit_leg, trans.back().new_state};
}

bool VertexData::canBounce(VertexState vs, int leg) const {
    VertexState new_state = flipSpinAtLeg(vs, leg);
    return allowed_[new_state];
}

void VertexData::print() const {
    std::cout << fmt::format("Vertex Data: {} allowed states, {} diagonal states\n",
                             allowed_states_.size(), diagonal_states_.size());
    std::cout << fmt::format("Energy offset: {:.6f}\n", energy_offset_);
    std::cout << fmt::format("Total diagonal weight: {:.6f}\n", total_diagonal_weight_);
    
    std::cout << "\nAllowed vertices:\n";
    for (VertexState vs : allowed_states_) {
        std::cout << fmt::format("  {:04b}: weight = {:.6f}, diagonal = {}\n",
                                 vs, weights_[vs], isDiagonal(vs));
    }
}

// VertexDataCollection implementation

VertexDataCollection::VertexDataCollection(const Hamiltonian& H) {
    int n_types = H.numBondTypes();
    data_.reserve(n_types);
    
    for (int t = 0; t < n_types; ++t) {
        data_.emplace_back(H, t);
        max_energy_offset_ = std::max(max_energy_offset_, data_.back().getEnergyOffset());
    }
}

} // namespace sse
