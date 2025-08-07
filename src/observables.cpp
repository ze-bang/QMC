#include "observables.h"
#include "lattice.h"
#include "hamiltonian.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iomanip>
#include <random>

namespace SSE {

// Observable implementation
void Observable::calculate_statistics() {
    if (samples.empty()) {
        value = 0.0;
        error = 0.0;
        return;
    }
    
    // Calculate mean
    value = std::accumulate(samples.begin(), samples.end(), 0.0) / samples.size();
    
    // Calculate standard error
    if (samples.size() > 1) {
        double variance = 0.0;
        for (double sample : samples) {
            variance += (sample - value) * (sample - value);
        }
        variance /= (samples.size() - 1);
        error = std::sqrt(variance / samples.size());
    } else {
        error = 0.0;
    }
}

// Observables implementation
Observables::Observables(const Lattice* lattice, const Hamiltonian* hamiltonian) 
    : lattice_(lattice), hamiltonian_(hamiltonian) {
    
    // Initialize correlation function storage
    int N = lattice_->size();
    spin_correlations_.resize(N, std::vector<double>(N, 0.0));
    correlation_distances_.resize(N * N);
    
    // Pre-calculate distances
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            correlation_distances_[i * N + j] = calculate_distance(i, j);
        }
    }
}

void Observables::measure_energy(const std::vector<int>& spins, int n_operators, double beta) {
    double energy = 0.0;
    
    // Energy from diagonal terms
    for (int bond = 0; bond < hamiltonian_->num_bonds(); ++bond) {
        energy += hamiltonian_->diagonal_matrix_element(bond, spins);
    }
    
    // Additional contribution from operator number (for SSE)
    energy -= static_cast<double>(n_operators) / beta;
    energy /= lattice_->size();  // Per site
    
    get_or_create_observable("energy").add_sample(energy);
}

void Observables::measure_magnetization(const std::vector<int>& spins) {
    double mag = 0.0;
    for (int spin : spins) {
        mag += (spin == 0) ? -0.5 : 0.5;  // Convert to ±1/2
    }
    mag /= lattice_->size();  // Per site
    
    get_or_create_observable("magnetization").add_sample(mag);
}

void Observables::measure_magnetization_squared(const std::vector<int>& spins) {
    double mag = 0.0;
    for (int spin : spins) {
        mag += (spin == 0) ? -0.5 : 0.5;
    }
    mag /= lattice_->size();
    
    get_or_create_observable("magnetization_squared").add_sample(mag * mag);
}

void Observables::measure_susceptibility(const std::vector<int>& spins) {
    double mag = 0.0;
    for (int spin : spins) {
        mag += (spin == 0) ? -0.5 : 0.5;
    }
    
    double chi = mag * mag / lattice_->size();  // ⟨M²⟩/N
    get_or_create_observable("susceptibility").add_sample(chi);
}

void Observables::measure_specific_heat(const std::vector<int>& spins, int n_operators, double beta) {
    // Specific heat is related to energy fluctuations
    // C = β²(⟨E²⟩ - ⟨E⟩²)
    
    double energy = 0.0;
    for (int bond = 0; bond < hamiltonian_->num_bonds(); ++bond) {
        energy += hamiltonian_->diagonal_matrix_element(bond, spins);
    }
    energy -= static_cast<double>(n_operators) / beta;
    energy /= lattice_->size();
    
    get_or_create_observable("energy_squared").add_sample(energy * energy);
}

void Observables::measure_spin_correlations(const std::vector<int>& spins) {
    int N = lattice_->size();
    
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            double si = (spins[i] == 0) ? -0.5 : 0.5;
            double sj = (spins[j] == 0) ? -0.5 : 0.5;
            spin_correlations_[i][j] += si * sj;
        }
    }
}

void Observables::measure_structure_factor(const std::vector<int>& spins) {
    // For square lattice, calculate S(π,π) and S(0,0)
    const SquareLattice* sq_lattice = dynamic_cast<const SquareLattice*>(lattice_);
    if (!sq_lattice) return;
    
    int Lx = sq_lattice->get_Lx();
    int Ly = sq_lattice->get_Ly();
    
    // S(0,0) - uniform structure factor
    double s_00_real = 0.0;
    for (int y = 0; y < Ly; ++y) {
        for (int x = 0; x < Lx; ++x) {
            int site = y * Lx + x;
            double sz = (spins[site] == 0) ? -0.5 : 0.5;
            s_00_real += sz;
        }
    }
    double s_00 = s_00_real * s_00_real / (Lx * Ly);
    get_or_create_observable("structure_factor_00").add_sample(s_00);
    
    // S(π,π) - staggered structure factor
    double s_pp_real = 0.0;
    for (int y = 0; y < Ly; ++y) {
        for (int x = 0; x < Lx; ++x) {
            int site = y * Lx + x;
            double sz = (spins[site] == 0) ? -0.5 : 0.5;
            double phase = ((x + y) % 2 == 0) ? 1.0 : -1.0;
            s_pp_real += phase * sz;
        }
    }
    double s_pp = s_pp_real * s_pp_real / (Lx * Ly);
    get_or_create_observable("structure_factor_pp").add_sample(s_pp);
}

void Observables::measure_staggered_magnetization(const std::vector<int>& spins) {
    double stag_mag = 0.0;
    
    for (int i = 0; i < lattice_->size(); ++i) {
        double sz = (spins[i] == 0) ? -0.5 : 0.5;
        double stagger_factor = calculate_staggered_factor(i);
        stag_mag += stagger_factor * sz;
    }
    stag_mag /= lattice_->size();
    
    get_or_create_observable("staggered_magnetization").add_sample(stag_mag);
}

void Observables::measure_uniform_susceptibility(const std::vector<int>& spins) {
    // χ_uniform = β Σ_ij ⟨S^z_i S^z_j⟩
    double chi_uniform = 0.0;
    
    for (int i = 0; i < lattice_->size(); ++i) {
        for (int j = 0; j < lattice_->size(); ++j) {
            double si = (spins[i] == 0) ? -0.5 : 0.5;
            double sj = (spins[j] == 0) ? -0.5 : 0.5;
            chi_uniform += si * sj;
        }
    }
    chi_uniform /= lattice_->size();
    
    get_or_create_observable("uniform_susceptibility").add_sample(chi_uniform);
}

void Observables::measure_staggered_susceptibility(const std::vector<int>& spins) {
    // χ_staggered = β Σ_ij (-1)^{i+j} ⟨S^z_i S^z_j⟩
    double chi_stag = 0.0;
    
    for (int i = 0; i < lattice_->size(); ++i) {
        for (int j = 0; j < lattice_->size(); ++j) {
            double si = (spins[i] == 0) ? -0.5 : 0.5;
            double sj = (spins[j] == 0) ? -0.5 : 0.5;
            double stagger_i = calculate_staggered_factor(i);
            double stagger_j = calculate_staggered_factor(j);
            chi_stag += stagger_i * stagger_j * si * sj;
        }
    }
    chi_stag /= lattice_->size();
    
    get_or_create_observable("staggered_susceptibility").add_sample(chi_stag);
}

void Observables::measure_all(const std::vector<int>& spins, int n_operators, double beta) {
    measure_energy(spins, n_operators, beta);
    measure_magnetization(spins);
    measure_magnetization_squared(spins);
    measure_susceptibility(spins);
    measure_specific_heat(spins, n_operators, beta);
    measure_structure_factor(spins);
    measure_staggered_magnetization(spins);
    measure_uniform_susceptibility(spins);
    measure_staggered_susceptibility(spins);
}

void Observables::calculate_all_statistics() {
    for (auto& [name, obs] : observables_) {
        obs.calculate_statistics();
    }
    
    // Calculate specific heat from energy fluctuations
    if (has_observable("energy") && has_observable("energy_squared")) {
        const auto& energy = get_observable("energy");
        const auto& energy_sq = get_observable("energy_squared");
        
        if (energy.samples.size() == energy_sq.samples.size() && !energy.samples.empty()) {
            std::vector<double> specific_heat_samples;
            for (size_t i = 0; i < energy.samples.size(); ++i) {
                double e = energy.samples[i];
                double e2 = energy_sq.samples[i];
                specific_heat_samples.push_back(e2 - e * e);  // β² factor handled elsewhere
            }
            
            auto& spec_heat = get_or_create_observable("specific_heat");
            spec_heat.samples = specific_heat_samples;
            spec_heat.calculate_statistics();
        }
    }
}

void Observables::reset_all() {
    for (auto& [name, obs] : observables_) {
        obs.reset();
    }
    
    // Reset correlation functions
    for (auto& row : spin_correlations_) {
        std::fill(row.begin(), row.end(), 0.0);
    }
}

const Observable& Observables::get_observable(const std::string& name) const {
    auto it = observables_.find(name);
    if (it == observables_.end()) {
        throw std::runtime_error("Observable " + name + " not found");
    }
    return it->second;
}

bool Observables::has_observable(const std::string& name) const {
    return observables_.find(name) != observables_.end();
}

std::vector<std::string> Observables::get_observable_names() const {
    std::vector<std::string> names;
    for (const auto& [name, obs] : observables_) {
        names.push_back(name);
    }
    return names;
}

void Observables::print_results() const {
    std::cout << "\nObservables Results:" << std::endl;
    std::cout << std::string(50, '=') << std::endl;
    
    for (const auto& [name, obs] : observables_) {
        if (obs.count > 0) {
            std::cout << std::setw(25) << std::left << name << ": "
                     << std::setw(12) << std::fixed << std::setprecision(6) << obs.value
                     << " ± " << std::setw(10) << std::setprecision(6) << obs.error
                     << " (" << obs.count << " samples)" << std::endl;
        }
    }
}

void Observables::save_results(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing" << std::endl;
        return;
    }
    
    file << "# Observable Results\n";
    file << "# Name\tValue\tError\tSamples\n";
    
    for (const auto& [name, obs] : observables_) {
        if (obs.count > 0) {
            file << name << "\t" << obs.value << "\t" << obs.error << "\t" << obs.count << "\n";
        }
    }
    
    file.close();
    std::cout << "Results saved to " << filename << std::endl;
}

void Observables::save_correlations(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing" << std::endl;
        return;
    }
    
    file << "# Spin correlation functions\n";
    file << "# i\tj\tdistance\tcorrelation\n";
    
    int N = lattice_->size();
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            double distance = correlation_distances_[i * N + j];
            double correlation = spin_correlations_[i][j];
            file << i << "\t" << j << "\t" << distance << "\t" << correlation << "\n";
        }
    }
    
    file.close();
    std::cout << "Correlations saved to " << filename << std::endl;
}

Observable& Observables::get_or_create_observable(const std::string& name) {
    auto it = observables_.find(name);
    if (it == observables_.end()) {
        it = observables_.emplace(name, Observable(name)).first;
    }
    return it->second;
}

double Observables::calculate_distance(int site1, int site2) const {
    const auto& pos1 = lattice_->get_site(site1).position;
    const auto& pos2 = lattice_->get_site(site2).position;
    
    if (pos1.size() != pos2.size()) return 0.0;
    
    double distance_sq = 0.0;
    for (size_t i = 0; i < pos1.size(); ++i) {
        double diff = pos1[i] - pos2[i];
        distance_sq += diff * diff;
    }
    
    return std::sqrt(distance_sq);
}

double Observables::calculate_staggered_factor(int site) const {
    // For square lattice: (-1)^{x+y}
    const SquareLattice* sq_lattice = dynamic_cast<const SquareLattice*>(lattice_);
    if (sq_lattice) {
        auto [x, y] = sq_lattice->get_coordinates(site);
        return ((x + y) % 2 == 0) ? 1.0 : -1.0;
    }
    
    // For other lattices, use simple alternating pattern
    return (site % 2 == 0) ? 1.0 : -1.0;
}

// AutocorrelationAnalysis implementation
AutocorrelationAnalysis::AutocorrelationAnalysis(const std::vector<double>& data) 
    : data_(data), tau_int_(0.0) {
    calculate_autocorrelation();
    calculate_integrated_time();
}

void AutocorrelationAnalysis::calculate_autocorrelation() {
    size_t N = data_.size();
    autocorr_.resize(N / 2);
    
    double mean = std::accumulate(data_.begin(), data_.end(), 0.0) / N;
    
    for (size_t t = 0; t < N / 2; ++t) {
        double corr = 0.0;
        for (size_t i = 0; i < N - t; ++i) {
            corr += (data_[i] - mean) * (data_[i + t] - mean);
        }
        autocorr_[t] = corr / (N - t);
    }
    
    // Normalize by variance
    if (autocorr_[0] > 0) {
        for (auto& corr : autocorr_) {
            corr /= autocorr_[0];
        }
    }
}

void AutocorrelationAnalysis::calculate_integrated_time() {
    tau_int_ = 0.5;  // Start with 1/2
    
    for (size_t t = 1; t < autocorr_.size(); ++t) {
        if (autocorr_[t] > 0) {
            tau_int_ += autocorr_[t];
        } else {
            break;  // Stop at first negative correlation
        }
        
        // Automatic windowing
        if (t > 6 * tau_int_) break;
    }
}

void AutocorrelationAnalysis::save_autocorrelation(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing" << std::endl;
        return;
    }
    
    file << "# Autocorrelation function\n";
    file << "# t\tC(t)\n";
    
    for (size_t t = 0; t < autocorr_.size(); ++t) {
        file << t << "\t" << autocorr_[t] << "\n";
    }
    
    file.close();
}

// BootstrapAnalysis implementation
BootstrapAnalysis::BootstrapAnalysis(const std::vector<double>& data, int n_bootstrap) 
    : data_(data), n_bootstrap_(n_bootstrap) {}

std::pair<double, double> BootstrapAnalysis::calculate_error() const {
    if (data_.empty()) return {0.0, 0.0};
    
    // Calculate original mean
    double original_mean = std::accumulate(data_.begin(), data_.end(), 0.0) / data_.size();
    
    // Generate bootstrap samples
    std::vector<double> bootstrap_means;
    std::random_device rd;
    std::mt19937 gen(rd());
    
    for (int i = 0; i < n_bootstrap_; ++i) {
        auto bootstrap_sample = generate_bootstrap_sample();
        double bootstrap_mean = std::accumulate(bootstrap_sample.begin(), bootstrap_sample.end(), 0.0) / bootstrap_sample.size();
        bootstrap_means.push_back(bootstrap_mean);
    }
    
    // Calculate bootstrap error
    double bootstrap_var = 0.0;
    double bootstrap_avg = std::accumulate(bootstrap_means.begin(), bootstrap_means.end(), 0.0) / bootstrap_means.size();
    
    for (double mean : bootstrap_means) {
        bootstrap_var += (mean - bootstrap_avg) * (mean - bootstrap_avg);
    }
    bootstrap_var /= (bootstrap_means.size() - 1);
    
    return {original_mean, std::sqrt(bootstrap_var)};
}

std::vector<double> BootstrapAnalysis::generate_bootstrap_sample() const {
    std::vector<double> bootstrap_sample;
    bootstrap_sample.reserve(data_.size());
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, data_.size() - 1);
    
    for (size_t i = 0; i < data_.size(); ++i) {
        bootstrap_sample.push_back(data_[dis(gen)]);
    }
    
    return bootstrap_sample;
}

} // namespace SSE
