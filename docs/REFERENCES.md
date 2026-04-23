# References

The implementation follows these foundational papers and reviews:

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

5. H. Flyvbjerg and H. G. Petersen, *Error estimates on averages of
   correlated data*, J. Chem. Phys. **91**, 461 (1989). The binning
   analysis used to produce error bars.

6. M. E. O'Neill, *PCG: A family of simple fast space-efficient
   statistically good algorithms for random number generation*, Tech.
   Rep. HMC-CS-2014-0905 (Harvey Mudd College, 2014). The RNG.

7. D. Lemire, *Fast random integer generation in an interval*, ACM
   Trans. Model. Comput. Simul. **29**, 1 (2019). The bounded-int
   draw.
