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

```bash
./sse_qmc config.json
```

See `examples/` directory for sample configurations.

## Algorithm

This implementation uses the stochastic series expansion formulation of quantum Monte Carlo, which represents the partition function as:

Z = Tr[(-βH)^n / n!] = Σ_n Σ_{S_M} ⟨S_M|(-βH)^n|S_M⟩ / n!

The algorithm performs updates by:
1. Diagonal updates (operator insertion/removal)
2. Off-diagonal updates (operator sequence modification)
3. Linked vertex updates for improved efficiency

## References

- Sandvik, A. W. (2010). Computational studies of quantum spin systems. AIP Conference Proceedings, 1297(1), 135-338.
- Evertz, H. G. (2003). The loop algorithm. Advances in Physics, 52(1), 1-66.
