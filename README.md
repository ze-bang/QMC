# qmc_sse

A modern, header-only **Stochastic Series Expansion** (SSE) Quantum Monte
Carlo library and command-line driver for spin-1/2 lattice systems.
Implements Sandvik's operator-loop algorithm for the antiferromagnetic
Heisenberg model on bipartite lattices. Designed to be small, fast,
hackable, and dependency-free (no JSON / Boost / FetchContent — just a
modern C++20 compiler and CMake).

```
                  H = J Σ_<ij> S_i · S_j      (J > 0, bipartite lattice)
```

## Features

- **Algorithm** — Sandvik (1999) operator-loop SSE with adaptive operator
  string truncation, free-spin flips, and the Marshall sign
  transformation (sign-problem free on bipartite lattices).
- **Lattices** — built-in 1D chain, 2D square, 2D honeycomb. Adding new
  bipartite lattices = supplying a list of bonds.
- **Observables** — energy, specific heat (jackknife), uniform / staggered
  magnetization (m, m², m⁴), uniform / staggered static
  susceptibilities, Binder cumulant.
- **Statistics** — Flyvbjerg–Petersen logarithmic binning analysis,
  jackknife errors for derived quantities, integrated autocorrelation
  time estimates.
- **Parallelism** — independent replicas in parallel via OpenMP; fast
  PCG32 RNG seeded from a single user-controlled value for full
  reproducibility per replica.
- **Tooling** — CMake build, IPO/LTO, `-march=native`, strict warnings,
  self-contained unit + physics test suite, and a throughput benchmark.

## Build

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

The build has zero external dependencies. OpenMP is auto-detected and
optional (turn off with `-DQMC_ENABLE_OPENMP=OFF`).

## Run a simulation

```bash
# Use a config file …
./build/qmc_sse examples/heisenberg_chain.conf

# …or override individual settings on the command line
./build/qmc_sse examples/heisenberg_square.conf beta=4.0 measurements=200000
```

Config syntax is `key = value` with `#` comments and optional `[section]`
headers (sections are ignored — flat namespace). See
[`examples/`](examples/) for ready-made runs.

### Important keys

| key              | default              | meaning |
|------------------|----------------------|---------|
| `lattice`        | `chain`              | `chain` / `square` / `honeycomb` |
| `Lx`, `Ly`       | `16`, `Lx`           | linear sizes (Ly ignored for chain) |
| `J`              | `1.0`                | exchange (must be > 0) |
| `beta`           | `4.0`                | inverse temperature |
| `thermalization` | `5000`               | warm-up sweeps (with adaptive M) |
| `measurements`   | `20000`              | measurement sweeps |
| `measure_every`  | `1`                  | binning frequency |
| `replicas`       | `1`                  | independent runs (run in parallel) |
| `seed`           | `0xc0ffee5eed5eed5e` | RNG seed (reproducibility) |
| `output_csv`     | (empty)              | optional running-mean time series |

## Example output

```
==============================================================
  qmc_sse  v0.1.0  (2026-04-23T12:34:56)
==============================================================
  lattice         : chain[L=32]
  N_sites         : 32
  N_bonds         : 32
  bipartite       : yes
  model           : antiferromagnetic Heisenberg
  J               : 1
  beta            : 4
  thermalization  : 5000 sweeps
  measurements    : 50000 sweeps
  replicas        : 4
  threads (OpenMP): 8
==============================================================

…

===== combined estimates over 4 replica(s) =====
  energy/site                 -4.42e-01  +/-  3.1e-04
  |m_z|                        2.71e-02  +/-  4.0e-04
  m_stag^2                     5.93e-02  +/-  6.2e-04
  chi_uniform                  1.10e-01  +/-  1.4e-03
  chi_stag                     1.45e+00  +/-  1.6e-02
```

## Project layout

```
include/qmc/        single-header library (everything is inline)
  types.hpp           common typedefs / op-code packing
  rng.hpp             PCG32 RNG + Lemire bounded ints
  config.hpp          dependency-free key=value parser
  logging.hpp         thread-safe logger
  lattice.hpp         lattice abstraction + chain / square / honeycomb
  heisenberg.hpp      AF Heisenberg bond Hamiltonian parameters
  operator_string.hpp SSE operator string + linked-vertex list
  sse_engine.hpp      diagonal + operator-loop + free-spin updates
  observable.hpp      log-binning analysis + jackknife
  measurements.hpp    standard observable set
  output.hpp          CSV writer / timestamps
  qmc.hpp             umbrella header
apps/qmc_sse.cpp    CLI driver (OpenMP replica parallelism)
tests/              minimal home-grown test harness + tests
benchmarks/         throughput benchmark
examples/           ready-to-run config files
docs/               algorithm notes & references
```

## Algorithm

- [`docs/ALGORITHM.md`](docs/ALGORITHM.md) — quick markdown notes on
  the implementation (data structures, update rules, estimators).
- [`docs/sse_qmc.tex`](docs/sse_qmc.tex) — full pedagogical write-up
  in LaTeX (derivation of the SSE expansion, Marshall sign
  transformation, operator-loop construction, observables, error
  analysis). Build with `make -C docs` to produce `docs/sse_qmc.pdf`.
- [`docs/REFERENCES.md`](docs/REFERENCES.md) — the literature the
  implementation actually follows.

## Limitations

- Only the AF Heisenberg model is implemented at present. The data
  structures (operator string, linked vertices, binning analysis) are
  generic; adding e.g. the transverse-field Ising model is mostly a
  matter of writing a new diagonal-update / vertex-table module.
- Frustrated lattices (triangular, kagome, …) are *not* supported by
  the operator-loop update — they suffer from the standard sign
  problem.
- No off-the-shelf parallel tempering / replica exchange yet (each
  replica runs at the same `(beta, J)`).

## License

MIT.  See [`LICENSE`](LICENSE).
