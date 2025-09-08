# SSE Quantum Monte Carlo

A stochastic series expansion quantum Monte Carlo implementation for spin-1/2 systems with arbitrary lattice geometry and Hamiltonians.

## Features

- Stochastic Series Expansion (SSE) algorithm
- Support for arbitrary lattice geometries
- Configurable Hamiltonians
- Efficient measurement of observables
- OpenMP parallelization support

## Building

```bash
mkdir build
cd build
cmake ..
make
```

## Usage

Command-line interface (run from build directory):

```bash
./sse_qmc --lattice square --model heisenberg --Lx 8 --Ly 8 --beta 10 --therm 10000 --meas 50000 --interval 10 --seed 42 --output run1
```

Key options:

--lattice square|triangular|chain|honeycomb
--model heisenberg|ising|xy|xxz
--Lx, --Ly lattice sizes (Ly ignored for chain)
--beta inverse temperature β
--therm thermalization sweeps
--meas measurement sweeps
--interval measurement interval
--J coupling (default 1.0)
--Delta anisotropy for XXZ
--h magnetic field for Heisenberg (optional)
--output prefix for saved results

Output files produced:

run1_params.txt (simulation parameters)
run1_observables.txt (mean and error per observable)
run1_correlations.txt (raw spin correlation accumulators)

## Algorithm

This implementation uses the stochastic series expansion formulation of quantum Monte Carlo, which represents the partition function as:

Z = Tr[(-βH)^n / n!] = Σ_n Σ_{S_M} ⟨S_M|(-βH)^n|S_M⟩ / n!

The algorithm performs updates by:
1. Diagonal updates (operator insertion/removal)
2. Off-diagonal updates (operator sequence modification)
3. Loop updates via linked vertex list

Specific heat: the stored observable "specific_heat" contains the variance of energy per site; multiply by β² to obtain C per site.
