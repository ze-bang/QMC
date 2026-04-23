// SPDX-License-Identifier: MIT
//
// Command-line driver for the SSE QMC library.
//
// Usage:
//   qmc_sse [config.conf] [key=value ...]
//
// Either supply a config file (see examples/) or override individual
// keys on the command line. Multiple independent replicas are run in
// parallel via OpenMP and their measurements are merged at the end.

#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "qmc/qmc.hpp"

namespace {

void print_banner(const qmc::Config& cfg, const qmc::Lattice& lat,
                  qmc::Real beta, qmc::Real J, int n_replicas) {
    std::cout << "==============================================================\n"
              << "  qmc_sse  v" << qmc::kVersion << "  (" << qmc::timestamp() << ")\n"
              << "==============================================================\n"
              << "  lattice         : " << lat.name() << "\n"
              << "  N_sites         : " << lat.n_sites() << "\n"
              << "  N_bonds         : " << lat.n_bonds() << "\n"
              << "  bipartite       : " << (lat.is_bipartite() ? "yes" : "no") << "\n"
              << "  model           : antiferromagnetic Heisenberg\n"
              << "  J               : " << J << "\n"
              << "  beta            : " << beta << "\n"
              << "  thermalization  : " << cfg.get_or<int>("thermalization", 5000) << " sweeps\n"
              << "  measurements    : " << cfg.get_or<int>("measurements", 20000) << " sweeps\n"
              << "  measure_every   : " << cfg.get_or<int>("measure_every", 1) << "\n"
              << "  replicas        : " << n_replicas << "\n"
#ifdef _OPENMP
              << "  threads (OpenMP): " << omp_get_max_threads() << "\n"
#else
              << "  threads (OpenMP): 1 (built without OpenMP)\n"
#endif
              << "==============================================================\n";
}

qmc::Config parse_args(int argc, char** argv) {
    qmc::Config cfg;
    int i = 1;
    if (argc > 1 && std::string(argv[1]).find('=') == std::string::npos) {
        cfg = qmc::Config::from_file(argv[1]);
        ++i;
    }
    for (; i < argc; ++i) {
        std::string s = argv[i];
        const auto eq = s.find('=');
        if (eq == std::string::npos) {
            throw std::runtime_error(
                "qmc_sse: extra argument '" + s + "' is not key=value");
        }
        cfg.set(s.substr(0, eq), s.substr(eq + 1));
    }
    return cfg;
}

} // namespace

int main(int argc, char** argv) try {
    qmc::Config cfg = parse_args(argc, argv);

    const std::string lattice_kind = cfg.get_or<std::string>("lattice", std::string("chain"));
    const int Lx                   = cfg.get_or<int>("Lx", 16);
    const int Ly                   = cfg.get_or<int>("Ly", Lx);
    const qmc::Real J              = cfg.get_or<qmc::Real>("J", 1.0);
    const qmc::Real beta           = cfg.get_or<qmc::Real>("beta", 4.0);
    const int therm                = cfg.get_or<int>("thermalization", 5000);
    const int meas                 = cfg.get_or<int>("measurements", 20000);
    const int every                = cfg.get_or<int>("measure_every", 1);
    const int n_replicas           = cfg.get_or<int>("replicas", 1);
    const std::uint64_t base_seed  = cfg.get_or<std::uint64_t>("seed", 0xc0ffee5eed5eed5eULL);
    const std::string out_csv      = cfg.get_or<std::string>("output_csv", std::string{});
    const bool verbose             = cfg.get_or<bool>("verbose", true);

    auto lattice = qmc::make_lattice(lattice_kind, Lx, Ly);
    qmc::HeisenbergModel model{J};
    model.validate(lattice);

    print_banner(cfg, lattice, beta, J, n_replicas);

    // One measurement object per replica, merged at the end via
    // weighted averaging (we use a separate, single-replica final
    // accumulator that re-bins all the samples for the report).
    std::vector<std::unique_ptr<qmc::Measurements>> per_replica;
    std::vector<qmc::SseStats>                      per_replica_stats(n_replicas);
    per_replica.reserve(n_replicas);
    for (int r = 0; r < n_replicas; ++r) per_replica.emplace_back(nullptr);

    std::vector<std::vector<qmc::Real>> energy_series(n_replicas);
    std::vector<std::vector<qmc::Real>> ms_sq_series(n_replicas);

    const auto t0 = std::chrono::steady_clock::now();
    std::atomic<int> replicas_done{0};

    #pragma omp parallel for schedule(dynamic, 1)
    for (int r = 0; r < n_replicas; ++r) {
        qmc::Pcg32 rng(base_seed + 1000003ULL * static_cast<std::uint64_t>(r),
                       0x14057B7EF767814FULL);
        qmc::SseEngine engine(lattice, model, beta, rng, /*initial_M=*/32);

        for (int s = 0; s < therm; ++s) engine.thermalize_step();

        auto m = std::make_unique<qmc::Measurements>(engine);
        for (int s = 0; s < meas; ++s) {
            engine.mc_step();
            if ((s % every) == 0) {
                m->measure();
                energy_series[r].push_back(m->energy().mean());
                ms_sq_series[r].push_back(m->ms_sq().mean());
            }
        }
        per_replica[r]       = std::move(m);
        per_replica_stats[r] = engine.stats();

        const int done = ++replicas_done;
        if (verbose) {
            #pragma omp critical
            std::printf("[replica %3d/%d done]\n", done, n_replicas);
        }
    }

    // Merge: sum means weighted by (independent) samples. Since all
    // replicas were run for the same number of steps we just average
    // their reported means and combine errors in quadrature scaled by
    // 1/sqrt(R).
    auto combine = [&](auto get_obs) {
        qmc::Real mean = 0.0;
        qmc::Real var  = 0.0;
        for (int r = 0; r < n_replicas; ++r) {
            const qmc::Real m = get_obs(*per_replica[r]).mean();
            const qmc::Real e = get_obs(*per_replica[r]).stderr_();
            mean += m;
            var  += e * e;
        }
        mean /= n_replicas;
        const qmc::Real err = std::sqrt(var) / n_replicas;
        return std::pair<qmc::Real, qmc::Real>{mean, err};
    };

    const auto t1   = std::chrono::steady_clock::now();
    const double dt = std::chrono::duration<double>(t1 - t0).count();
    std::cout << "\nrun finished in " << dt << " s ("
              << (n_replicas * meas / dt) << " measurement sweeps / s, summed over replicas)\n\n";

    // Per-replica detailed report.
    for (int r = 0; r < n_replicas; ++r) {
        std::cout << "----- replica " << r << " -----\n";
        per_replica[r]->report(std::cout);
        const auto& st = per_replica_stats[r];
        std::cout << "  [diag insert acc = " << st.diag_insert_acc()
                  << ", remove acc = " << st.diag_remove_acc()
                  << ", loops flipped = " << st.loop_flip_frac()
                  << ", n_max = " << st.max_n_observed
                  << ", M = " << st.cutoff_M << "]\n";
    }

    // Combined summary across replicas (most useful number to report).
    std::cout << "\n===== combined estimates over " << n_replicas << " replica(s) =====\n";
    auto pr = [&](const std::string& name, std::pair<qmc::Real, qmc::Real> v) {
        std::printf("  %-26s  %+.8e  +/-  %.3e\n", name.c_str(), v.first, v.second);
    };
    pr("energy/site",          combine([](const qmc::Measurements& m) -> const qmc::Observable& { return m.energy(); }));
    pr("|m_z|",                combine([](const qmc::Measurements& m) -> const qmc::Observable& { return m.mz_abs(); }));
    pr("m_z^2",                combine([](const qmc::Measurements& m) -> const qmc::Observable& { return m.mz_sq(); }));
    pr("|m_stag|",             combine([](const qmc::Measurements& m) -> const qmc::Observable& { return m.ms_abs(); }));
    pr("m_stag^2",             combine([](const qmc::Measurements& m) -> const qmc::Observable& { return m.ms_sq(); }));
    pr("chi_uniform",          combine([](const qmc::Measurements& m) -> const qmc::Observable& { return m.chi_uniform(); }));
    pr("chi_stag",             combine([](const qmc::Measurements& m) -> const qmc::Observable& { return m.chi_stag(); }));

    // Optional CSV output: per-measurement running mean of energy / m_s^2.
    if (!out_csv.empty()) {
        qmc::CsvWriter w(out_csv, {"replica", "step", "energy_running_mean", "ms_sq_running_mean"});
        for (int r = 0; r < n_replicas; ++r) {
            const auto& es = energy_series[r];
            const auto& ms = ms_sq_series[r];
            for (std::size_t i = 0; i < es.size(); ++i) {
                w.write_row({static_cast<qmc::Real>(r),
                             static_cast<qmc::Real>(i),
                             es[i], ms[i]});
            }
        }
        std::cout << "\nWrote running-mean time series to: " << out_csv << '\n';
    }

    return 0;
} catch (const std::exception& e) {
    std::cerr << "qmc_sse error: " << e.what() << '\n';
    return 1;
}
