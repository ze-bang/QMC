#pragma once

/**
 * @file lattice.hpp
 * @brief Lattice structure definitions for arbitrary geometries
 * 
 * Supports various lattice types: chain, square, triangular, honeycomb,
 * kagome, and fully custom lattices defined via bond lists.
 */

#include <vector>
#include <array>
#include <string>
#include <memory>
#include <functional>
#include <unordered_map>
#include "sse/types.hpp"

namespace sse {

/**
 * @brief Enumeration of built-in lattice types
 */
enum class LatticeType {
    Chain,          // 1D chain
    Square,         // 2D square lattice
    Triangular,     // 2D triangular lattice
    Honeycomb,      // 2D honeycomb lattice
    Kagome,         // 2D kagome lattice
    Cubic,          // 3D cubic lattice
    Custom          // User-defined lattice
};

/**
 * @brief Lattice class for arbitrary spin geometries
 * 
 * The lattice stores site positions and bond connectivity.
 * Bonds can have different types for inhomogeneous couplings.
 */
class Lattice {
public:
    /**
     * @brief Default constructor (empty lattice)
     */
    Lattice() = default;
    
    /**
     * @brief Construct lattice from number of sites and bonds
     */
    Lattice(SiteIdx n_sites, const std::vector<Bond>& bonds);
    
    /**
     * @brief Create 1D chain lattice
     * @param L Number of sites
     * @param periodic Use periodic boundary conditions
     */
    static Lattice chain(SiteIdx L, bool periodic = true);
    
    /**
     * @brief Create 2D square lattice
     * @param Lx Size in x direction
     * @param Ly Size in y direction
     * @param periodic Use periodic boundary conditions
     */
    static Lattice square(SiteIdx Lx, SiteIdx Ly, bool periodic = true);
    
    /**
     * @brief Create 2D triangular lattice
     */
    static Lattice triangular(SiteIdx Lx, SiteIdx Ly, bool periodic = true);
    
    /**
     * @brief Create 2D honeycomb lattice
     */
    static Lattice honeycomb(SiteIdx Lx, SiteIdx Ly, bool periodic = true);
    
    /**
     * @brief Create 2D kagome lattice
     */
    static Lattice kagome(SiteIdx Lx, SiteIdx Ly, bool periodic = true);
    
    /**
     * @brief Create 3D cubic lattice
     */
    static Lattice cubic(SiteIdx Lx, SiteIdx Ly, SiteIdx Lz, bool periodic = true);
    
    /**
     * @brief Create custom lattice from bond list
     */
    static Lattice custom(SiteIdx n_sites, const std::vector<Bond>& bonds);
    
    // Accessors
    SiteIdx numSites() const { return n_sites_; }
    BondIdx numBonds() const { return static_cast<BondIdx>(bonds_.size()); }
    
    const Bond& getBond(BondIdx b) const { return bonds_[b]; }
    const std::vector<Bond>& getBonds() const { return bonds_; }
    
    /**
     * @brief Get all bonds connected to a site
     */
    const std::vector<BondIdx>& siteBonds(SiteIdx site) const { 
        return site_bonds_[site]; 
    }
    
    /**
     * @brief Get neighbors of a site
     */
    std::vector<SiteIdx> neighbors(SiteIdx site) const;
    
    /**
     * @brief Get coordination number of a site
     */
    int coordination(SiteIdx site) const { 
        return static_cast<int>(site_bonds_[site].size()); 
    }
    
    /**
     * @brief Get maximum coordination number
     */
    int maxCoordination() const { return max_coordination_; }
    
    /**
     * @brief Get number of bond types
     */
    int numBondTypes() const { return n_bond_types_; }
    
    /**
     * @brief Get lattice dimensions (for regular lattices)
     */
    const std::vector<SiteIdx>& dimensions() const { return dimensions_; }
    
    /**
     * @brief Get site coordinates (for regular lattices)
     */
    std::vector<int> siteCoords(SiteIdx site) const;
    
    /**
     * @brief Convert coordinates to site index
     */
    SiteIdx coordsToSite(const std::vector<int>& coords) const;
    
    /**
     * @brief Get lattice type
     */
    LatticeType type() const { return type_; }
    
    /**
     * @brief Print lattice info
     */
    void print() const;
    
private:
    LatticeType type_ = LatticeType::Custom;
    SiteIdx n_sites_ = 0;
    std::vector<SiteIdx> dimensions_;
    std::vector<Bond> bonds_;
    std::vector<std::vector<BondIdx>> site_bonds_;
    int max_coordination_ = 0;
    int n_bond_types_ = 0;
    
    void buildSiteBonds();
};

} // namespace sse
