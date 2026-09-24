# aztec-diagonal-defects

Companion code for

> **C. Won**, *Composite diagonal defects in Aztec-diamond dimer coverings: boundary pinning, Coulomb attraction, and exact asymptotics* (manuscript in preparation).

The model extends domino tilings of the Aztec diamond (the dimer description of Shakti-type artificial spin ice) by a composite defect: a tile joining two diagonally adjacent cells of the same colour. B tiles `(2m,2n)–(2m+1,2n−1)` cover black cells, W tiles `(2m,2n+1)–(2m−1,2n)` cover white cells. They live on complementary period-2 sublattices, never cross, and must occur in equal numbers `k`.

Everything in the paper — exact counts, local statistics, bulk properties, the boundary theorems and their asymptotics, and the thermodynamics — can be reproduced with the scripts in `scripts/`.

## Installation

```bash
pip install -r requirements.txt        # numpy, scipy, sympy, python-flint, matplotlib, networkx, pytest
python -m pytest -q tests              # 23 checks, ~10 s
```

Python ≥ 3.9. `python-flint` supplies fast exact modular and integer determinants.

## Package `aztecdefects`

| module | content |
|---|---|
| `lattice` | cells, colours, diagonal partners, placements A/B of the Aztec diamond, boundary staircases |
| `modular` | Pfaffian, inverse and interpolation modulo primes; Chinese remaindering |
| `pfaffian` | Pfaffian orientation of the non-bipartite covering graph; exact `Z_n(u)=Σ c_k u^k` via `det K(u) = Z(u)^2`; exact k-resolved single-tile and pair counts; double-precision tile probabilities |
| `kasteleyn` | bipartite Kasteleyn matrix; exact tiling counts; k = 1 pair weights in double precision and exactly by multimodular banded elimination (reaches n = 121) |
| `boundary` | Theorem 1 (`N_n(a,b)`), Theorem 2 (the `2k×2k` determinant with the cyclic sign τ), large-n evaluation for k = 1 and k = 2 |
| `bulk` | periodic model on the torus: spectral gap, bulk B–W pair weight, free energy φ(u), pair density, thermodynamics |
| `shuffling` | generalized urban renewal for a face with a chord; face classes |
| `thermo` | entropy S(k) = ln c_k, free-energy landscape, thresholds T_k, k*(T) |
| `bruteforce` | direct enumeration used for validation |

## Reproducing the paper

Run from the repository root (outputs go to `data/` and `figures/`):

| paper | script | runtime* |
|---|---|---|
| Table I, Eq. (6), Sec. III | `python scripts/01_polynomials.py --nmax 13` (`--nmax 17` for Figs. 13–14) | 5 min (17: ~35 min) |
| Sec. IV, Figs. 3–5 | `python scripts/02_n9_statistics.py` | 2 min |
| Sec. V, Fig. 6 (and bulk data for Figs. 11, 13) | `python scripts/03_bulk.py` | 5 min |
| Sec. VI, Fig. 11 | `python scripts/04_equal_energy.py --nmax 61` | 3 min |
| Sec. IX, Fig. 9 | `python scripts/05_k1_large_n.py --n 21 41 61 81 101 121` | 20 min |
| Sec. IX, Fig. 10 | `python scripts/06_boundary_asymptotics.py --nmax 1601` | 1 min |
| Sec. X, Fig. 12 | `python scripts/07_k2.py` | 2 min |
| Sec. VII, Figs. 13–14 | `python scripts/08_thermodynamics.py` | 5 min |
| Figs. 1, 2, 7, 8 | `python scripts/09_schematics.py` | seconds |
| Secs. VIII–X (proof ingredients) | `python scripts/10_verify_proofs.py --nmax 10` | 2 min |

\*single core. `bash make_all.sh` runs everything with the default settings. The cached results in `data/` (polynomials up to n = 17, bulk data, exact k = 1 weights up to n = 81) let the figure scripts run without recomputation; `data/n9_local.npz` is regenerated on the first run of script 02.

## Key results reproduced

* Bulk B–W pair weight `∝ r^{-2.00}`; spectral gap `min|P| = t^4 = u^2` (Sec. V).
* `Z_{2m}(u) = 4^m Z_{2m−1}(u)` (generalized urban renewal, Sec. VIII).
* **Theorem 1**: `M(A_n − b_i − w_j) = 2^{n(n−1)/2} N_n(i, n−1−j)`, `N_n(a,b) = Σ_k 2^k C(n−1−k, a−k) C(n−1−k, b−k)`.
* **Theorem 2**: `c_k/c_0 = 2^{−2kn} |det[τ ν]|` for boundary defects.
* One pair: `c_1/c_0 → (2/π) 4^n n^{-3}` (extrapolated 0.6366165), variance along the edge `n/8`, B–W correlation `4/n`.
* Two pairs: like charges on opposite edges; `c_2/c_0 ≈ (2/π)^2 16^n n^{-6}`, `c_2 c_0/c_1^2 → 1`.
* Equal energies (u = 1): diagonal fraction `0.16264 + 0.61/n`; no frozen corners.
* Thermodynamics: onset `T* ≈ ε_d/(n ln 2)`; staircase `k*(T)` filling the arctic tangency points in doublets.

## Conventions

* Placement A: cells `(i, j)` with `|i+½| + |j+½| ≤ n`; black ⇔ `i+j` even.
* Weight `u` per B tile (hence per B–W pair). With an energy `ε_d` per diagonal tile, `u = exp(−2ε_d/k_BT)`.
* Boundary staircases are indexed clockwise: NE from the N corner, SE from E, SW from S, NW from W.

## License

MIT (see `LICENSE`).
