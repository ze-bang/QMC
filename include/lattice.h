#pragma once

#include <vector>
#include <string>
#include <unordered_map>

namespace SSE {

struct Site {
    int index;
    std::vector<double> position;
    
    Site(int idx, std::vector<double> pos) : index(idx), position(std::move(pos)) {}
};

struct Bond {
    int site1, site2;
    double coupling;
    std::string type;
    
    Bond(int s1, int s2, double J, std::string t = "exchange") 
        : site1(s1), site2(s2), coupling(J), type(std::move(t)) {}
};

class Lattice {
protected:
    std::vector<Site> sites_;
    std::vector<Bond> bonds_;
    std::vector<std::vector<int>> neighbors_;
    int dimension_;
    
public:
    Lattice(int dim) : dimension_(dim) {}
    virtual ~Lattice() = default;
    
    // Pure virtual methods to be implemented by derived classes
    virtual void generate_lattice() = 0;
    virtual std::string get_name() const = 0;
    
    // Common interface
    void add_site(int index, const std::vector<double>& position);
    void add_bond(int site1, int site2, double coupling, const std::string& type = "exchange");
    void build_neighbor_list();
    
    // Getters
    int size() const { return static_cast<int>(sites_.size()); }
    int num_bonds() const { return static_cast<int>(bonds_.size()); }
    int dimension() const { return dimension_; }
    
    const Site& get_site(int i) const { return sites_[i]; }
    const Bond& get_bond(int i) const { return bonds_[i]; }
    const std::vector<int>& get_neighbors(int site) const { return neighbors_[site]; }
    const std::vector<Site>& get_sites() const { return sites_; }
    const std::vector<Bond>& get_bonds() const { return bonds_; }
    
    // Utility methods
    virtual double get_coordination_number() const;
    void print_info() const;
};

// Square lattice implementation
class SquareLattice : public Lattice {
private:
    int Lx_, Ly_;
    bool periodic_;
    
public:
    SquareLattice(int Lx, int Ly, bool periodic = true);
    
    void generate_lattice() override;
    std::string get_name() const override { return "Square Lattice"; }
    
    // Square lattice specific methods
    int get_site_index(int x, int y) const;
    std::pair<int, int> get_coordinates(int site) const;
    int get_Lx() const { return Lx_; }
    int get_Ly() const { return Ly_; }
    bool is_periodic() const { return periodic_; }
};

// Triangular lattice implementation
class TriangularLattice : public Lattice {
private:
    int Lx_, Ly_;
    bool periodic_;
    
public:
    TriangularLattice(int Lx, int Ly, bool periodic = true);
    
    void generate_lattice() override;
    std::string get_name() const override { return "Triangular Lattice"; }
    
    int get_site_index(int x, int y) const;
    std::pair<int, int> get_coordinates(int site) const;
    int get_Lx() const { return Lx_; }
    int get_Ly() const { return Ly_; }
    bool is_periodic() const { return periodic_; }
};

// 1D Chain implementation
class Chain : public Lattice {
private:
    int L_;
    bool periodic_;
    
public:
    Chain(int L, bool periodic = true);
    
    void generate_lattice() override;
    std::string get_name() const override { return "1D Chain"; }
    
    int get_L() const { return L_; }
    bool is_periodic() const { return periodic_; }
};

// Honeycomb lattice implementation
class HoneycombLattice : public Lattice {
private:
    int Lx_, Ly_;
    bool periodic_;
    
public:
    HoneycombLattice(int Lx, int Ly, bool periodic = true);
    
    void generate_lattice() override;
    std::string get_name() const override { return "Honeycomb Lattice"; }
    
    int get_site_index(int x, int y, int sublattice) const;
    int get_Lx() const { return Lx_; }
    int get_Ly() const { return Ly_; }
    bool is_periodic() const { return periodic_; }
};

} // namespace SSE
