# SSE with operator-loop updates — implementation notes

This file is a short field guide to the algorithm coded up in
`include/qmc/sse_engine.hpp`. It is **not** a tutorial; for that consult
the references in [`REFERENCES.md`](REFERENCES.md).

## 1. Model and basis

We work with the spin-1/2 antiferromagnetic Heisenberg model on a
bipartite lattice,

$$
H = J \sum_{\langle i j\rangle} \mathbf{S}_i \cdot \mathbf{S}_j, \qquad J > 0.
$$

After the **Marshall sign transformation** on one sublattice
($S^\pm \to -S^\pm$ on B sites) the off-diagonal pieces become
non-positive operators, so the SSE expansion has positive weights:

$$
H = -J \sum_b \big( H_{1,b} + H_{2,b} \big) + \tfrac{J}{4} N_b\, \mathbb{1},
$$

with the bond operators

$$
H_{1,b} = \tfrac14 - S^z_i S^z_j, \qquad
H_{2,b} = \tfrac12 (S^+_i S^-_j + S^-_i S^+_j).
$$

Both $H_{1,b}$ and $H_{2,b}$ have matrix element $1/2$ on antiparallel
spin pairs, zero otherwise.

## 2. Operator string

We work in the $S^z$ basis with spins $s_i \in \{-1,+1\}$. The partition
function is expanded as

$$
Z = \sum_{\alpha} \sum_{S_M} \frac{(-J)^n \beta^n (M-n)!}{M!}
    \big\langle \alpha\big| \prod_{p=1}^{M} O_{a_p, b_p} \big|\alpha\big\rangle,
$$

where $S_M = (a_1,b_1; \dots; a_M,b_M)$ is an *operator string* of length
$M$ padded with identities, $n$ counts non-identity operators, and
$O_{0, b} = H_{1,b}$, $O_{1,b} = H_{2,b}$. The sign $(-J)^n$ is positive
on a bipartite lattice — there is **no sign problem**.

Internally we represent the string as `OpCode op = (bond << 1) | type`,
with `kIdentity = -1`. See `include/qmc/operator_string.hpp`.

## 3. Diagonal update (Metropolis sweep)

We sweep $p = 0, 1, \dots, M-1$. For each position:

| current state of `op[p]` | proposal | acceptance |
|--------------------------|----------|------------|
| identity                 | insert diagonal on a random bond $b$ | $\beta\,N_b\,\langle\alpha\|H_{1,b}\|\alpha\rangle / (M - n)$ |
| diagonal                 | remove                                | $(M - n + 1) / (\beta\,N_b\,\langle\alpha\|H_{1,b}\|\alpha\rangle)$ |
| off-diagonal             | propagate the spin state through the bond | (always) |

The matrix element is $1/2$ whenever the two spins on the bond are
antiparallel and zero otherwise; insertion is rejected if the bond
spins are parallel.

## 4. Linked-vertex list

To run the off-diagonal (loop) update we build the **linked-vertex
list** of the operator string. Each non-identity operator is a vertex
with four legs:

```
        leg2  ─────────  leg3            (top, after operator)
                │
              vertex
                │
        leg0  ─────────  leg1            (bottom, before operator)
        (site i)         (site j)
```

For every site we walk forward in imaginary time (operator-string
order) and connect each leg of the current vertex to the previous
unmatched leg on the same site. The first / last legs on each site
are linked across the periodic time boundary.

`link[ℓ]` is the leg connected to leg `ℓ`; it is `-1` for sites with no
operators (free spins). See `LinkedVertices::build` in
`operator_string.hpp`.

## 5. Operator-loop update

Loops live on the leg graph. For the AF Heisenberg model the
deterministic *switch* rule is exact and bounce-free:

> Enter vertex through leg `ℓ`, exit through leg `ℓ XOR 1` — the leg on
> the same time-side (top or bottom) but on the **other** site of the
> bond. Then jump along `link[exit_leg]` to the next vertex.

Each vertex therefore participates in at most two loops (one from its
bottom side and one from its top side). The implementation uses a
single `visited[]` byte array; a value of `1` means "leg belongs to a
loop that we *will* flip", `2` means "loop kept as-is".

Once all loops are constructed we apply two finalization passes:

1. **Operator type update.** Toggle a vertex's diagonal/off-diagonal
   bit iff exactly one of its time-sides was flipped (XOR). Both flipped
   or both kept leaves the type alone.
2. **Spin update.** A site's persistent spin (the value at $p=0$) is
   flipped iff the very first leg encountered for that site
   (`first_leg[s]`) is in a flipped loop.

Each independent loop is flipped with probability $1/2$.

## 6. Free-spin flip

Sites with no operators are free worldlines: their spin can be flipped
independently with probability $1/2$. This step is essential at very
low operator densities (small $\beta$, small $J$) where many sites are
"unvisited" by the operator string.

## 7. Adaptive truncation

During thermalization the cutoff $M$ is grown to roughly $1.4 \times
\max_t n_t$ (Sandvik's recommendation). Once the run reaches the
measurement phase $M$ is frozen — measurements assume a fixed weight
normalization.

## 8. Observables

The estimators implemented in `measurements.hpp` are standard:

- **Energy per site:** $\;E/N = -\langle n\rangle/(\beta N) + J N_b/(4N)$
- **Specific heat (jackknife):** $\;C = \big(\langle n^2\rangle - \langle n\rangle^2 - \langle n\rangle\big)/N$
- **Magnetizations:** $\;|m_z|, m_z^2, |m_s|, m_s^2, m_s^4$ at the
  $p = 0$ slice (sufficient because $S^z_{\text{tot}}$ is conserved by
  the SSE updates).
- **Susceptibilities:** $\;\chi_q = \beta\,\langle (M_q)^2 \rangle / N$
  for the static uniform / staggered components.

Errors are reported via Flyvbjerg–Petersen logarithmic binning; for
`specific_heat()` we use a delete-one jackknife on the raw samples.

## 9. Complexity

A full Monte Carlo step costs $\mathcal{O}(M + N_{\mathrm{sites}})$ for
the diagonal update plus $\mathcal{O}(\langle n\rangle)$ for the loop
update (each leg of the linked-vertex list is touched at most once).
Memory scales as $\mathcal{O}(M + 4 \langle n\rangle)$.
