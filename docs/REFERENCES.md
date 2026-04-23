# References

The implementation follows these foundational papers and reviews.

## SSE engine (`qmc_sse`)

1. A. W. Sandvik and J. Kurkijärvi, *Quantum Monte Carlo simulation
   method for spin systems*, Phys. Rev. B **43**, 5950 (1991).
   Introduces the SSE expansion.

2. A. W. Sandvik, *Stochastic series expansion method with
   operator-loop cluster updates*, Phys. Rev. B **59**, R14157 (1999).
   The operator-loop algorithm used here.

3. O. F. Syljuåsen and A. W. Sandvik, *Quantum Monte Carlo with
   directed loops*, Phys. Rev. E **66**, 046701 (2002). Generalisation
   to arbitrary models via the directed-loop equations (not yet
   implemented in this codebase but the data structures support it).

4. A. W. Sandvik, *Computational Studies of Quantum Spin Systems*,
   AIP Conf. Proc. **1297**, 135 (2010); arXiv:1101.3281. Pedagogical
   review with practical estimators (energy, specific heat, etc.).

## DQMC engine (`qmc_dqmc`)

5. R. Blankenbecler, D. J. Scalapino, R. L. Sugar, *Monte Carlo
   calculations of coupled boson-fermion systems. I*, Phys. Rev. D
   **24**, 2278 (1981). The original BSS Determinant QMC algorithm.

6. J. E. Hirsch, *Discrete Hubbard-Stratonovich transformation for
   fermion lattice models*, Phys. Rev. B **28**, 4059 (1983).
   The discrete S^z-channel HS field used here.

7. E. Y. Loh Jr. and J. E. Gubernatis, *Stable numerical simulations
   of models of interacting electrons in condensed-matter physics*,
   in *Electronic Phase Transitions*, ed. W. Hanke and Yu. V. Kopaev,
   Modern Problems in Condensed Matter Sciences **32**, 177 (1992).
   QR-based stabilisation of the propagator product.

8. F. F. Assaad and H. G. Evertz, *World-line and Determinantal Quantum
   Monte Carlo Methods for Spins, Phonons and Electrons*, Lect. Notes
   Phys. **739**, 277 (2008). Comprehensive modern review.

## Common infrastructure

9. H. Flyvbjerg and H. G. Petersen, *Error estimates on averages of
   correlated data*, J. Chem. Phys. **91**, 461 (1989). The binning
   analysis used to produce error bars.

10. M. E. O'Neill, *PCG: A family of simple fast space-efficient
    statistically good algorithms for random number generation*, Tech.
    Rep. HMC-CS-2014-0905 (Harvey Mudd College, 2014). The RNG.

11. D. Lemire, *Fast random integer generation in an interval*, ACM
    Trans. Model. Comput. Simul. **29**, 1 (2019). The bounded-int
    draw.
