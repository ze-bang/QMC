/**
 * @file sse_config.cpp
 * @brief Implementation of SSE configuration management
 */

#include "sse/sse_config.hpp"
#include <algorithm>

namespace sse {

SSEConfig::SSEConfig(const Lattice& lattice, const VertexDataCollection& vertex_data)
    : lattice_(lattice), vertex_data_(vertex_data) {
    
    spins_.resize(lattice.numSites());
    
    // Initial operator string size estimate
    int init_size = std::max(4, static_cast<int>(lattice.numSites()));
    operators_.resize(init_size, OperatorCode::identity());
    
    first_vertex_.resize(lattice.numSites(), -1);
    last_vertex_.resize(lattice.numSites(), -1);
}

void SSEConfig::initialize(Random& rng) {
    // Random initial spin configuration
    for (SiteIdx s = 0; s < lattice_.numSites(); ++s) {
        spins_[s] = rng.randomSpin();
    }
    
    // Clear operator string
    std::fill(operators_.begin(), operators_.end(), OperatorCode::identity());
    n_operators_ = 0;
}

void SSEConfig::ensureCapacity(int min_size) {
    if (static_cast<int>(operators_.size()) < min_size) {
        resize(min_size);
    }
}

void SSEConfig::resize(int new_size) {
    int old_size = static_cast<int>(operators_.size());
    operators_.resize(new_size, OperatorCode::identity());
    
    // Reallocate vertex list
    vertex_links_.resize(new_size * NUM_LEGS, -1);
    
    // Fill new slots with identity
    for (int p = old_size; p < new_size; ++p) {
        operators_[p] = OperatorCode::identity();
    }
}

void SSEConfig::buildVertexList() {
    int M = operatorStringLength();
    
    // Resize and clear vertex list
    vertex_links_.assign(M * NUM_LEGS, -1);
    std::fill(first_vertex_.begin(), first_vertex_.end(), -1);
    std::fill(last_vertex_.begin(), last_vertex_.end(), -1);
    
    // Traverse operator string and link vertices
    for (int p = 0; p < M; ++p) {
        OperatorCode op = operators_[p];
        if (op.isIdentity()) continue;
        
        BondIdx bond_idx = op.bond();
        const Bond& bond = lattice_.getBond(bond_idx);
        
        // Site indices for this bond
        SiteIdx sites[2] = {bond.i, bond.j};
        
        // Link to previous vertex on each site
        for (int leg = 0; leg < 2; ++leg) {  // Bottom legs
            SiteIdx site = sites[leg];
            int v_current = vertexIndex(p, leg);
            
            int v_prev = last_vertex_[site];
            if (v_prev >= 0) {
                // Link to previous vertex's top leg
                int prev_leg = vertexToLeg(v_prev);
                int prev_top = v_prev - prev_leg + leg + 2;  // Corresponding top leg
                
                vertex_links_[v_current] = prev_top;
                vertex_links_[prev_top] = v_current;
            } else {
                // First vertex on this site
                first_vertex_[site] = v_current;
            }
            
            // Update last vertex (now points to current's bottom leg)
            last_vertex_[site] = vertexIndex(p, leg);
        }
    }
    
    // Link first and last vertices (periodic boundary in imaginary time)
    for (SiteIdx site = 0; site < lattice_.numSites(); ++site) {
        int v_first = first_vertex_[site];
        int v_last = last_vertex_[site];
        
        if (v_first >= 0 && v_last >= 0) {
            // Find which legs these correspond to
            int first_p = vertexToOperator(v_first);
            int first_leg = vertexToLeg(v_first);
            int last_p = vertexToOperator(v_last);
            int last_leg = vertexToLeg(v_last);
            
            // Link last's top leg to first's bottom leg
            int first_bottom = vertexIndex(first_p, first_leg);
            int last_top = vertexIndex(last_p, last_leg + 2);
            
            vertex_links_[first_bottom] = last_top;
            vertex_links_[last_top] = first_bottom;
        }
    }
}

} // namespace sse
