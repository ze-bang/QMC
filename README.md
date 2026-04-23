# qmc_sse / qmc_dqmc

A modern, header-only **Quantum Monte Carlo** library for lattice
many-body systems, with two complementary engines:

| binary       | algorithm                                  | model                              |
|--------------|--------------------------------------------|------------------------------------|
| `qmc_sse`    | Stochastic Series Expansion (Sandvik 1999) | spin-½ AF Heisenberg, bipartite    |
| `qmc_dqmc`   | Determinant QMC (BSS + Hirsch HS)          | single-band Hubbard, bipartite     |

Both are written for clarity, hackability, and reproducibility, share
common infrastructure (lattice, RNG, binning analysis, configuration,
OpenMP replica parallelism), and have **zero external dependencies**
beyond a C++20 compiler and CMake.

## Features

### Common infrastructure
- Self-contained dense linear algebra (GEMM / LU / QR / sym-eig /
  matrix exponential), plain C++20.
- Flyvbjerg–Petersen logarithmic binning analysis with integrated
  autocorrelation time, jackknife errors for derived quantities.
- PCG32 RNG seeded from a single user-controlled value for full
  per-replica reproducibility.
- CMake build, IPO/LTO, `-march=native`, strict warnings, OpenMP
  replica parallelism, and a self-contained test harness.

### `qmc_sse` (SSE for spin-½ Heisenberg)
- Sandvik operator-loop algorithm with adaptive operator-string
  truncation, free-spin flips, Marshall sign transformation
  (sign-problem free on bipartite lattices).
- Lattices: 1D chain, 2D square, 2D honeycomb (PBC).
- Observables: energy, specific heat, uniform/staggered magnetisation
  (m, m², m⁴), uniform/staggered susceptibilities, Binder cumulant.

### `qmc_dqmc` (BSS Determinant QMC for the Hubbard model)
- Hirsch discrete Hubbard–Stratonovich transformation coupling to S^z.
- Single-spin-flip local updates with Sherman–Morrison Green's-function
  updates (O(N²) per accepted move, O(N) per proposal).
- Numerical stabilization by **incremental Householder QR** every
  `n_stab` time slices (Loh–Gubernatis-style scheme); Q-orthogonality
  is preserved throughout, so simulations stay accurate at low
  temperatures.
- Sign-problem-free at half-filling on any bipartite lattice; the
  engine reports the average sign as a diagnostic.
- Observables: density, double occupancy, equal-time spin–spin
  correlations, antiferromagnetic structure factor S(π,π).

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
# SSE — spin-½ Heisenberg
./build/qmc_sse  examples/heisenberg_chain.conf
./build/qmc_sse  examples/heisenberg_square.conf beta=4.0 measurements=200000

# DQMC — half-filled square Hubbard
./build/qmc_dqmc examples/hubbard_square_4x4.conf
./build/qmc_dqmc lattice=square Lx=6 Ly=6 t=1.0 U=4.0 beta=8.0 \
                 L_tau=80 n_stab=10 measurements=2000 replicas=4
```

Config syntax is `key = value` with `#` comments and optional
`[section]` headers (sections are ignored — flat namespace). See
[`examples/`](examples/) for ready-made runs.

### `qmc_sse` keys

| key              | default              | meaning |
|------------------|----------------------|---------|
| `lattice`        | `chain`              | `chain` / `square` / `honeycomb` |
| `Lx`, `Ly`       | `16`, `Lx`           | linear sizes (Ly ignored for chain) |
| `J`              | `1.0`                | exchange (must be > 0) |
| `beta`           | `4.0`                | inverse temperature |
| `thermalization` | `5000`               | warm-up sweeps (with adaptive M) |
| `measurements`   | `20000`              | measurement sweeps |
| `replicas`       | `1`                  | independent runs (run in parallel) |
| `seed`           | `0xc0ffee5eed5eed5e` | RNG seed (reproducibility) |

### `qmc_dqmc` keys

| key              | default              | meaning |
|------------------|----------------------|---------|
| `lattice`        | `square`             | `chain` / `square` / `honeycomb` |
| `Lx`, `Ly`       | `4`, `Lx`            | linear sizes |
| `t`              | `1.0`                | hopping |
| `U`              | `4.0`                | on-site interaction (≥ 0) |
| `mu`             | `0.0`                | chemical potential (`mu=0` = half-filling) |
| `beta`           | `4.0`                | inverse temperature |
| `L_tau`          | `40`                 | Trotter slices (`Δτ = β / Lτ`) |
| `n_stab`         | `10`                 | recompute G every n_stab slices via QR |
| `thermalization` | `200`                | warm-up sweeps |
| `measurements`   | `1000`               | measurement sweeps |
| `replicas`       | `1`                  | independent runs (run in parallel) |
| `seed`           | `0xfeedfacecafebeef` | RNG seed |

## Example DQMC output

```
==============================================================
  qmc_dqmc  v0.1.0  (2026-04-23T...)
==============================================================
  lattice         : square[4x4]
  N_sites         : 16        bipartite : yes
  model           : Hubbard (Hirsch HS, S^z channel)
  t = 1   U = 4   mu = 0  (= half-filling on bipartite lattice)
  beta = 4   L_tau = 40   dt = 0.1   n_stab = 10
==============================================================

----- replica 0 -----
<sign>             1.000000e+00 +/- 0.000e+00
density            1.000000e+00 +/- 7.5e-09
double_occupancy   1.236e-01    +/- 1.6e-03
<m_z^2>            3.66e-03     +/- 6.0e-04
<m_stag^2>         4.10e-02     +/- 5.5e-03
S(pi,pi)           6.55e-01     +/- 8.8e-02
  [acc = 0.64, <sign> = 1.00, max stab drift = 0.00]
```

## Project layout

```
include/qmc/                single-header library (everything inline)
  linalg.hpp                 dense linear algebra (Matrix / GEMM / LU / QR / expm / stable G)
  types.hpp                  common typedefs / SSE op-code packing
  rng.hpp                    PCG32 RNG + Lemire bounded ints
  config.hpp                 dependency-free key=value parser
  logging.hpp                thread-safe logger
  lattice.hpp                lattice abstraction (chain / square / honeycomb)

  ── SSE engine ────────────
  heisenberg.hpp             AF Heisenberg bond Hamiltonian parameters
  operator_string.hpp        SSE operator string + linked-vertex list
  sse_engine.hpp             diagonal + operator-loop + free-spin updates
  measurements.hpp           SSE observables

  ── DQMC engine ───────────
  hubbard.hpp                Hubbard parameters, kinetic K matrix
  dqmc_engine.hpp            B-builder, sweep, Sherman-Morrison, wrap, recompute
  dqmc_measurements.hpp      DQMC equal-time observables (Wick contractions)

  observable.hpp             log-binning analysis + jackknife
  output.hpp                 CSV writer / timestamps
  qmc.hpp                    umbrella header

apps/qmc_sse.cpp           SSE  CLI driver (OpenMP replica parallelism)
apps/qmc_dqmc.cpp          DQMC CLI driver (OpenMP replica parallelism)
tests/                     minimal home-grown test harness + 29 tests
benchmarks/                throughput benchmark
examples/                  ready-to-run config files
docs/                      algorithm notes & full LaTeX write-ups
```

## Algorithm documentation

- [`docs/ALGORITHM.md`](docs/ALGORITHM.md) — SSE algorithm summary.
- [`docs/DQMC.md`](docs/DQMC.md) — DQMC algorithm summary.
- [`docs/sse_qmc.tex`](docs/sse_qmc.tex) — full pedagogical SSE write-up.
- [`docs/dqmc.tex`](docs/dqmc.tex) — full pedagogical DQMC write-up.
- [`docs/REFERENCES.md`](docs/REFERENCES.md) — the literature both engines follow.

Build the LaTeX with `make -C docs` to produce `sse_qmc.pdf` and
`dqmc.pdf`.

## Limitations

- **`qmc_sse`**: only the AF Heisenberg model. Frustrated lattices
  (triangular, kagome, …) suffer from the standard sign problem and
  are not supported by the operator-loop update.
- **`qmc_dqmc`**: only the repulsive single-band Hubbard model with
  the Hirsch S^z-coupling HS field; the engine still runs away from
  half-filling but the average sign decays exponentially with `β U`.
  The dense linear algebra is written for clarity rather than raw
  speed — for production studies you would link against a tuned
  BLAS+LAPACK and add full UDV (singular-value) stabilization.
- No replica exchange / parallel tempering yet (each replica runs at
  the same parameters).

## License

MIT.  See [`LICENSE`](LICENSE).
