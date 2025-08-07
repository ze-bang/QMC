#pragma once

#include <vector>
#include <string>
#include <unordered_map>
#include <memory>

namespace SSE {

class Lattice;
class Hamiltonian;

struct Observable {
    std::string name;
    double value;
    double error;
    long long count;
    std::vector<double> samples;
    
    Observable(std::string n) : name(std::move(n)), value(0.0), error(0.0), count(0) {}
    
    void add_sample(double sample) {
        samples.push_back(sample);
        count++;
    }
    
    void calculate_statistics();
    void reset() {
        value = 0.0;
        error = 0.0;
        count = 0;
        samples.clear();
    }
};

class Observables {
private:
    std::unordered_map<std::string, Observable> observables_;
    const Lattice* lattice_;
    const Hamiltonian* hamiltonian_;
    
    // Correlation function storage
    std::vector<std::vector<double>> spin_correlations_;
    std::vector<double> correlation_distances_;
    
public:
    Observables(const Lattice* lattice, const Hamiltonian* hamiltonian);
    
    // Basic observables
    void measure_energy(const std::vector<int>& spins, int n_operators, double beta);
    void measure_magnetization(const std::vector<int>& spins);
    void measure_magnetization_squared(const std::vector<int>& spins);
    void measure_susceptibility(const std::vector<int>& spins);
    void measure_specific_heat(const std::vector<int>& spins, int n_operators, double beta);
    
    // Correlation functions
    void measure_spin_correlations(const std::vector<int>& spins);
    void measure_structure_factor(const std::vector<int>& spins);
    
    // Advanced observables
    void measure_staggered_magnetization(const std::vector<int>& spins);
    void measure_uniform_susceptibility(const std::vector<int>& spins);
    void measure_staggered_susceptibility(const std::vector<int>& spins);
    
    // Measurement control
    void measure_all(const std::vector<int>& spins, int n_operators, double beta);
    void calculate_all_statistics();
    void reset_all();
    
    // Data access
    const Observable& get_observable(const std::string& name) const;
    bool has_observable(const std::string& name) const;
    std::vector<std::string> get_observable_names() const;
    
    // Output
    void print_results() const;
    void save_results(const std::string& filename) const;
    void save_correlations(const std::string& filename) const;
    
private:
    Observable& get_or_create_observable(const std::string& name);
    double calculate_distance(int site1, int site2) const;
    std::vector<int> get_momentum_vector(int kx, int ky, int Lx, int Ly) const;
    double calculate_staggered_factor(int site) const;
};

// Autocorrelation analysis
class AutocorrelationAnalysis {
private:
    std::vector<double> data_;
    std::vector<double> autocorr_;
    double tau_int_;
    
public:
    explicit AutocorrelationAnalysis(const std::vector<double>& data);
    
    void calculate_autocorrelation();
    void calculate_integrated_time();
    
    double get_integrated_time() const { return tau_int_; }
    const std::vector<double>& get_autocorrelation() const { return autocorr_; }
    
    void save_autocorrelation(const std::string& filename) const;
};

// Bootstrap error analysis
class BootstrapAnalysis {
private:
    std::vector<double> data_;
    int n_bootstrap_;
    
public:
    BootstrapAnalysis(const std::vector<double>& data, int n_bootstrap = 1000);
    
    std::pair<double, double> calculate_error() const;
    std::vector<double> generate_bootstrap_sample() const;
};

} // namespace SSE
