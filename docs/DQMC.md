# Determinant QMC algorithm — implementation notes

This file is the reference for the `qmc::DqmcEngine` implementation and
companion to the SSE notes in [`ALGORITHM.md`](ALGORITHM.md). For the
detailed pedagogical derivation see [`dqmc.tex`](dqmc.tex).

## Model and conventions

We simulate the single-band Hubbard model on a generic lattice (chain,
square, honeycomb) with periodic boundary conditions:

```
H = -t Σ<ij>,σ (c†_iσ c_jσ + h.c.)
    - μ Σ_iσ n_iσ
    + U Σ_i (n_i↑ - 1/2)(n_i↓ - 1/2)
```

The symmetric `(n - 1/2)` form makes `μ = 0` correspond to half-filling on
bipartite lattices by particle–hole symmetry. The non-interacting part
is encoded in a single-particle matrix `K` with `K_ii = -μ` and
`K_ij = -t` for nearest neighbours.

## Algorithm

1. **Trotter decomposition** at step `Δτ = β / Lτ`:

   `e^{-βH} ≈ ∏_l e^{-Δτ K} e^{-Δτ V_l}`.

2. **Hirsch discrete HS transformation** of the on-site quartic term:

   `e^{-Δτ U (n↑-1/2)(n↓-1/2)} = (e^{-ΔτU/4}/2) Σ_{s=±1} e^{α s (n↑ - n↓)}`

   with `cosh α = e^{Δτ U / 2}`. The Hubbard–Stratonovich field
   `s_{l,i} ∈ {-1, +1}` lives on every site of every time slice.

3. **B-matrices** (built once per slice from the current aux field):

   `B_{l,σ}  = exp(-Δτ K) · diag(exp(σ α s_{l,i}))`
   `B_{l,σ}^{-1} = diag(exp(-σ α s_{l,i})) · exp(+Δτ K)`

4. **Equal-time Green's function** in the convention used here:

   `G_σ(l) = [I + A_σ(l)]^{-1}`, with
   `A_σ(l) = B_{l-1,σ} B_{l-2,σ} ⋯ B_{0,σ} B_{Lτ-1,σ} ⋯ B_{l,σ}`.

   This is the *pre-`B_l`* Green's function — the convention in which
   the local-flip determinant ratio is the simple
   `r = 1 + Δ(1 − G_{ii})` (see step 5 below).

5. **Single-spin-flip update.** A flip `s_{l,i} → −s_{l,i}` yields
   `Δ_↑ = exp(-2α s_{l,i}) − 1` and `Δ_↓ = exp(+2α s_{l,i}) − 1`.
   Determinant ratios:

   ```
   R_↑ = 1 + Δ_↑ (1 - G_↑(i,i))
   R_↓ = 1 + Δ_↓ (1 - G_↓(i,i))
   ```

   Accept with probability `min(1, |R_↑ R_↓|)`. On accept, update the
   Green's function by Sherman–Morrison:

   `G' = G − (Δ / R) (I − G)_{:,i} G_{i,:}`

   in components, `G'_{jk} = G_{jk} − (Δ/R)(δ_{ji} − G_{ji}) G_{ik}`.

6. **Wrapping** to the next slice (after all updates at slice `l`):

   `G_σ ← B_{l,σ} G_σ B_{l,σ}^{-1}`

   With our `B_l = expK · diag(d)`, the correct ordering is
   `expK · (diag(d) · G · diag(1/d)) · expmK` — `diag(1/d)` does **not**
   commute with `expmK`.

7. **Stabilization.** Every `n_stab` slices, re-build `G_σ` from
   scratch using **incremental QR**:

   - Maintain `(Q, R)` with the partial product equal to `Q · R`.
   - Every `n_stab` left-multiplies of `B`, refactorize
     `(B_chunk · Q) = Q_new · R_local` and update `R := R_local · R`.
   - Final `G = (Q^T + R)^{-1} Q^T`.

   This keeps `Q` orthogonal throughout, taming the catastrophic
   conditioning of `exp(-Δτ K)^{Lτ}` for low temperatures. It is the
   standard Loh–Gubernatis-style scheme.

## Sign problem

On any **bipartite** lattice at `μ = 0` (half-filling), particle–hole
symmetry guarantees `det(M_↑) = det(M_↓)` for every aux-field
configuration, so the weight `det(M_↑) det(M_↓) ≥ 0` and the average
sign is exactly `+1`. The engine reports `<sign>` as a diagnostic; on
the test problems (square `Lx × Lx`, honeycomb, chain at `μ = 0`)
this is `1.0 ± 0` to machine precision.

Away from half-filling there is in general a sign problem that grows
exponentially with `β U`. The simulation still runs and reports the
sign; reweight observables by `<O sign> / <sign>` if needed.

## Observables

Computed each measurement from the current `(G_↑, G_↓)`:

| Observable | Formula |
| ---------- | ------- |
| `<n_iσ>`   | `1 - G_σ(i,i)` |
| `<n>`      | `(1/Ns) Σ_i (<n_i↑> + <n_i↓>)` |
| double occ | `(1/Ns) Σ_i (1 - G_↑(i,i)) (1 - G_↓(i,i))` |
| `<S^z_i S^z_j>` | Wick: `¼[(1-G_↑(i,i))(1-G_↑(j,j)) + (↑↔↓) − 2(1-G_↑(i,i))(1-G_↓(j,j)) + (δ_{ij}-G_↑(j,i))G_↑(i,j) + (↑↔↓)]` |
| `<m_z^2>`  | `(1/Ns²) Σ_{ij} <S^z_i S^z_j>` |
| `S(π,π)`   | `(1/Ns) Σ_{ij} (-1)^{a_i+a_j} <S^z_i S^z_j>`, `a_i` = sublattice |

Errors via the same Flyvbjerg–Petersen logarithmic binning as the SSE
code (`qmc::Observable`).

## Verification

The test suite (`tests/test_dqmc.cpp`) exercises:

1. **Free fermion limit (`U = 0`).** With `α = 0` the aux field
   decouples and every configuration gives the closed-form
   `G = (I + e^{-βK})^{-1}`. The engine matches this to ≲10⁻¹³ at
   initialization *and after a full sweep* (so the wrap+stabilize
   cycle is verified to high precision).
2. **Half-filling sum rule.** For `μ = 0` on bipartite lattices,
   `<n> = 1` exactly per configuration (PH symmetry); checked to
   ~10⁻⁹ on the L=4 chain at `U = 4`, `β = 2`.
3. **Average sign.** Exactly `+1` at half-filling.
4. **Double occupancy bounds.** `0 < <n_↑ n_↓> < 1/2` and decreases
   with `U`, as expected.

Sanity runs:

| System            | β | `<n>`   | dbl occ | `S(π,π)` | sign |
| ----------------- | - | ------- | ------- | -------- | ---- |
| 4×4, U=4          | 4 | 1.000   | 0.124   | 0.66     | 1.0  |
| 6×6, U=4          | 8 | 1.000   | 0.120   | 1.58     | 1.0  |
| chain L=8, U=2    | 4 | 1.000   | 0.169   | 0.36     | 1.0  |

`S(π,π)` grows with both system size and β — the expected antiferro-
magnetic ordering tendency of the half-filled square Hubbard model.

## Limitations of this implementation

The code is engineered to be educational, self-contained, and correct
on small to moderate problem sizes. For larger production studies you
would want to:

- Link against a tuned BLAS+LAPACK (the dense ops in `linalg.hpp` are
  written for clarity, not raw speed).
- Replace the QR stabilization with full UDV (singular-value) splitting,
  which is more robust at very low temperatures.
- Add Wolff/global-move type updates and time-displaced measurements
  for spectral function calculations.
- Implement delayed/submatrix updates (Alvarez et al.) to amortize the
  Sherman–Morrison cost.

These are extensions, not corrections — the algorithm is correct as
implemented within its targeted regime (`β U ≲ 30`, lattices up to
~10×10).
