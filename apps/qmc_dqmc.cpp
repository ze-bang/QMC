// SPDX-License-Identifier: MIT
//
// Command-line driver for the DQMC engine on the Hubbard model.
//
// Usage:
//   qmc_dqmc [config.conf] [key=value ...]

#include <atomic>
#include <chrono>
#include <cstdio>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "qmc/qmc.hpp"

namespace {

void print_banner(const qmc::Config& cfg, const qmc::Lattice& lat,
                  const qmc::HubbardModel& m, const qmc::DqmcParams& p,
                  int n_replicas) {
    std::cout << "==============================================================\n"
              << "  qmc_dqmc  v" << qmc::kVersion << "  (" << qmc::timestamp() << ")\n"
              << "==============================================================\n"
              << "  lattice         : " << lat.name() << "\n"
              << "  N_sites         : " << lat.n_sites() << "\n"
              << "  bipartite       : " << (lat.is_bipartite() ? "yes" : "no") << "\n"
              << "  model           : Hubbard (Hirsch HS, S^z channel)\n"
              << "  t               : " << m.t  << "\n"
              << "  U               : " << m.U  << "\n"
              << "  mu              : " << m.mu
              << "  (mu = 0 = half-filling on bipartite lattice)\n"
              << "  beta            : " << p.beta << "\n"
              << "  L_tau           : " << p.L_tau
              << "    (dt = "       << p.beta / p.L_tau << ")\n"
              << "  n_stab          : " << p.n_stab << "\n"
              << "  thermalization  : " << cfg.get_or<int>("thermalization", 200) << " sweeps\n"
              << "  measurements    : " << cfg.get_or<int>("measurements", 1000) << " sweeps\n"
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
        if (eq == std::string::npos)
            throw std::runtime_error("qmc_dqmc: extra arg '" + s + "' is not key=value");
        cfg.set(s.substr(0, eq), s.substr(eq + 1));
    }
    return cfg;
}

} // namespace

int main(int argc, char** argv) try {
    qmc::Config cfg = parse_args(argc, argv);

    const std::string lattice_kind = cfg.get_or<std::string>("lattice", std::string("square"));
    const int Lx                   = cfg.get_or<int>("Lx", 4);
    const int Ly                   = cfg.get_or<int>("Ly", Lx);
    const qmc::Real t              = cfg.get_or<qmc::Real>("t",  1.0);
    const qmc::Real U              = cfg.get_or<qmc::Real>("U",  4.0);
    const qmc::Real mu             = cfg.get_or<qmc::Real>("mu", 0.0);

    qmc::DqmcParams params;
    params.beta   = cfg.get_or<qmc::Real>("beta", 4.0);
    params.L_tau  = cfg.get_or<int>("L_tau", 40);
    params.n_stab = cfg.get_or<int>("n_stab", 10);

    const int therm                = cfg.get_or<int>("thermalization", 200);
    const int meas                 = cfg.get_or<int>("measurements", 1000);
    const int every                = cfg.get_or<int>("measure_every", 1);
    const int n_replicas           = cfg.get_or<int>("replicas", 1);
    const std::uint64_t base_seed  = cfg.get_or<std::uint64_t>("seed", 0xfeedfacecafebeefULL);

    auto lattice = qmc::make_lattice(lattice_kind, Lx, Ly);
    qmc::HubbardModel model{t, U, mu};
    model.validate();

    print_banner(cfg, lattice, model, params, n_replicas);

    std::vector<std::unique_ptr<qmc::DqmcMeasurements>> per_replica(n_replicas);
    std::vector<qmc::DqmcStats> per_replica_stats(n_replicas);

    const auto t0 = std::chrono::steady_clock::now();
    std::atomic<int> done{0};

    #pragma omp parallel for schedule(dynamic, 1)
    for (int r = 0; r < n_replicas; ++r) {
        qmc::Pcg32 rng(base_seed + 1000003ULL * static_cast<std::uint64_t>(r),
                       0x14057B7EF767814FULL);
        qmc::DqmcEngine engine(lattice, model, params, rng);
        for (int s = 0; s < therm; ++s) (void)engine.sweep();
        auto m = std::make_unique<qmc::DqmcMeasurements>(engine);
        for (int s = 0; s < meas; ++s) {
            const qmc::Real sgn = engine.sweep();
            if ((s % every) == 0) m->measure(sgn);
        }
        per_replica[r]       = std::move(m);
        per_replica_stats[r] = engine.stats();

        const int d = ++done;
        #pragma omp critical
        std::printf("[replica %3d/%d done]\n", d, n_replicas);
    }

    const auto t1 = std::chrono::steady_clock::now();
    const double dt = std::chrono::duration<double>(t1 - t0).count();
    std::cout << "\nrun finished in " << dt << " s ("
              << (n_replicas * meas / dt) << " sweeps / s, summed over replicas)\n\n";

    for (int r = 0; r < n_replicas; ++r) {
        std::cout << "----- replica " << r << " -----\n";
        per_replica[r]->report(std::cout);
        const auto& st = per_replica_stats[r];
        std::cout << "  [acc = " << st.acceptance()
                  << ", <sign> = " << st.average_sign()
                  << ", max stab drift = " << st.max_recompute_drift << "]\n";
    }

    return 0;
} catch (const std::exception& e) {
    std::cerr << "qmc_dqmc error: " << e.what() << '\n';
    return 1;
}
