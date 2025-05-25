#pragma once

#include <string>
#include <vector>
#include <memory>
#include <map>
#include <unordered_map>
#include <functional>
#include <random>
#include <Eigen/Dense>

namespace qmc {

/**
 * @brief Base class for lattice implementations
 * 
 * Defines the interface for arbitrary lattice structures.
 */
class Lattice {
public:
    using Index = int;
    using Position = Eigen::Vector3d;
    
    struct Site {
        Index index;
        Position position;
        std::vector<Index> neighbors;
    };
    
    Lattice() = default;
    virtual ~Lattice() = default;
    
    // Pure virtual methods that derived classes must implement
    virtual size_t size() const = 0;
    virtual const Site& getSite(Index index) const = 0;
    virtual std::vector<std::pair<Index, Index>> getBonds() const = 0;
    
    // Default implementations
    virtual std::vector<Index> getNeighbors(Index site) const {
        return getSite(site).neighbors;
    }
    
    virtual Position getPosition(Index site) const {
        return getSite(site).position;
    }
    
    virtual std::string getName() const = 0;
    
protected:
    std::vector<Site> sites_;
};

/**
 * @brief Square lattice implementation
 */
class SquareLattice : public Lattice {
public:
    SquareLattice(int L, int W, bool periodic = true);
    ~SquareLattice() override = default;
    
    size_t size() const override { return sites_.size(); }
    const Site& getSite(Index index) const override { return sites_.at(index); }
    std::vector<std::pair<Index, Index>> getBonds() const override;
    std::string getName() const override { return "Square"; }
    
    // Square lattice specific methods
    Index siteIndex(int x, int y) const;
    std::pair<int, int> coordinates(Index index) const;
    int getL() const { return L_; }
    int getW() const { return W_; }
    
private:
    int L_; // Length
    int W_; // Width
    bool periodic_;
};

/**
 * @brief Triangular lattice implementation
 */
class TriangularLattice : public Lattice {
public:
    TriangularLattice(int L, int W, bool periodic = true);
    ~TriangularLattice() override = default;
    
    size_t size() const override { return sites_.size(); }
    const Site& getSite(Index index) const override { return sites_.at(index); }
    std::vector<std::pair<Index, Index>> getBonds() const override;
    std::string getName() const override { return "Triangular"; }
    
    // Triangular lattice specific methods
    Index siteIndex(int x, int y) const;
    std::pair<int, int> coordinates(Index index) const;
    
private:
    int L_; // Length
    int W_; // Width
    bool periodic_;
};

/**
 * @brief Custom lattice implementation
 * 
 * Allows for the creation of arbitrary lattices by specifying
 * sites and their connections manually.
 */
class CustomLattice : public Lattice {
public:
    CustomLattice(const std::vector<Position>& positions, 
                  const std::vector<std::vector<Index>>& neighbors);
    ~CustomLattice() override = default;
    
    size_t size() const override { return sites_.size(); }
    const Site& getSite(Index index) const override { return sites_.at(index); }
    std::vector<std::pair<Index, Index>> getBonds() const override;
    std::string getName() const override { return "Custom"; }
    
private:
    std::vector<std::pair<Index, Index>> bonds_;
};

} // namespace qmc
