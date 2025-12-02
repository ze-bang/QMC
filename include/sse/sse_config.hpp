#pragma once

/**
 * @file sse_config.hpp
 * @brief SSE configuration state (operator string and spin configuration)
 * 
 * The SSE configuration consists of:
 * 1. Operator string: sequence of operators in imaginary time
 * 2. Spin configuration: current spin state at τ=0
 * 3. Linked vertex list: connectivity for loop updates
 */

#include <vector>
#include <array>
#include "sse/types.hpp"
#include "sse/lattice.hpp"
#include "sse/vertex.hpp"

namespace sse {

/**
 * @brief SSE configuration management
 * 
 * Manages the operator string and provides efficient access patterns
 * for diagonal and loop updates.
 */
class SSEConfig {
public:
    /**
     * @brief Construct SSE configuration for given lattice
     */
    SSEConfig(const Lattice& lattice, const VertexDataCollection& vertex_data);
    
    /**
     * @brief Initialize with random spin configuration
     */
    void initialize(Random& rng);
    
    /**
     * @brief Get operator at position p
     */
    OperatorCode getOperator(int p) const { return operators_[p]; }
    
    /**
     * @brief Set operator at position p
     */
    void setOperator(int p, OperatorCode op) { operators_[p] = op; }
    
    /**
     * @brief Get current operator string length (cutoff M)
     */
    int operatorStringLength() const { return static_cast<int>(operators_.size()); }
    
    /**
     * @brief Get number of non-identity operators
     */
    int64_t numOperators() const { return n_operators_; }
    
    /**
     * @brief Set number of non-identity operators
     */
    void setNumOperators(int64_t n) { n_operators_ = n; }
    
    /**
     * @brief Increment number of operators
     */
    void incrementOperators() { ++n_operators_; }
    
    /**
     * @brief Decrement number of operators
     */
    void decrementOperators() { --n_operators_; }
    
    /**
     * @brief Get spin at site
     */
    StateIdx getSpin(SiteIdx site) const { return spins_[site]; }
    
    /**
     * @brief Set spin at site
     */
    void setSpin(SiteIdx site, StateIdx s) { spins_[site] = s; }
    
    /**
     * @brief Flip spin at site
     */
    void flipSpin(SiteIdx site) { spins_[site] ^= 1; }
    
    /**
     * @brief Get spin configuration
     */
    const std::vector<StateIdx>& getSpins() const { return spins_; }
    
    /**
     * @brief Get mutable spin configuration
     */
    std::vector<StateIdx>& getSpins() { return spins_; }
    
    /**
     * @brief Get operator string
     */
    const std::vector<OperatorCode>& getOperators() const { return operators_; }
    
    /**
     * @brief Get mutable operator string
     */
    std::vector<OperatorCode>& getOperators() { return operators_; }
    
    /**
     * @brief Resize operator string if needed
     */
    void ensureCapacity(int min_size);
    
    /**
     * @brief Resize operator string
     */
    void resize(int new_size);
    
    // Linked vertex list operations
    
    /**
     * @brief Build linked vertex list for loop updates
     */
    void buildVertexList();
    
    /**
     * @brief Get linked vertex at position v
     * @return Index of linked vertex, or -1 if not linked
     */
    int getLink(int v) const { return vertex_links_[v]; }
    
    /**
     * @brief Set link between vertices
     */
    void setLink(int v1, int v2) {
        vertex_links_[v1] = v2;
        vertex_links_[v2] = v1;
    }
    
    /**
     * @brief Get first vertex linked to site (at τ=0)
     * @return Vertex index or -1 if none
     */
    int getFirstVertex(SiteIdx site) const { return first_vertex_[site]; }
    
    /**
     * @brief Get last vertex linked to site (at τ=β)
     * @return Vertex index or -1 if none
     */
    int getLastVertex(SiteIdx site) const { return last_vertex_[site]; }
    
    /**
     * @brief Get vertex list (for direct access)
     */
    const std::vector<int>& getVertexLinks() const { return vertex_links_; }
    std::vector<int>& getVertexLinks() { return vertex_links_; }
    
    /**
     * @brief Convert operator position and leg to vertex index
     */
    int vertexIndex(int p, int leg) const {
        return p * NUM_LEGS + leg;
    }
    
    /**
     * @brief Extract operator position from vertex index
     */
    int vertexToOperator(int v) const {
        return v / NUM_LEGS;
    }
    
    /**
     * @brief Extract leg from vertex index
     */
    int vertexToLeg(int v) const {
        return v % NUM_LEGS;
    }
    
    /**
     * @brief Check if vertex is unlinked
     */
    bool isUnlinked(int v) const { return vertex_links_[v] < 0; }
    
    // References
    const Lattice& getLattice() const { return lattice_; }
    const VertexDataCollection& getVertexData() const { return vertex_data_; }
    
private:
    const Lattice& lattice_;
    const VertexDataCollection& vertex_data_;
    
    std::vector<StateIdx> spins_;          // Spin configuration at τ=0
    std::vector<OperatorCode> operators_;  // Operator string
    int64_t n_operators_ = 0;              // Number of non-identity operators
    
    // Linked vertex list for loop updates
    std::vector<int> vertex_links_;        // [p * 4 + leg] -> linked vertex
    std::vector<int> first_vertex_;        // [site] -> first vertex at τ=0
    std::vector<int> last_vertex_;         // [site] -> last vertex at τ=β
};

} // namespace sse
