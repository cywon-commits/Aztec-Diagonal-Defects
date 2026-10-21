# Shakti vertex-model realization (paper Sec. VIII, Table II, Figs. 10-11)

Run everything with `bash shakti/run_all.sh`.

Requires `ortools` (`pip install ortools`). Run scripts from `shakti/code/`.

## Dictionary
* plaquette <-> chessboard cell: a plaquette with centre (X, Y) (Shakti sample shifted by -12) is cell
  (floor((X+Y)/2), floor((Y-X)/2)).  Vertical plaquettes are black cells, horizontal plaquettes white.
* The parallel (block) partner of a plaquette is exactly `aztecdefects.lattice.partner`.
  B = two vertical plaquettes, W = two horizontal plaquettes.
* Covering from islands (`multiplicity.tiling`, `opentiling.tiling_open`): charged block centre -> O tile;
  unhappy z=3 vertex -> pairs the flat plaquette with the plaquette selected by the perpendicular island.

Note: `code/` scripts by the other analysis (cp.py, cp2.py, flip.py, ...) use block parity
"(a+b) odd = horizontal"; the companion vorticity paper uses the opposite parity.

## Results reproduced by the scripts
| script | result |
|---|---|
| `multiplicity.py 8`, `multiplicity.py 12` | periodic samples, random coverings with k = 0, 1, 2: exactly 2 island configurations per covering (global reversal) -> uniform ensemble, g_O = 1 |
| `aztec_sample.py 3 0 2 3`, `aztec_sample.py 2 1 2 3` | parallel-partnered AD_n as open Shakti sample: 2 c0 ground states, 2 c1 lowest one-pair states at 2 eps_d, exact position distribution |
| `opentiling.py 3 0`, `opentiling.py 2 1` | each of these states maps to a covering; every covering appears exactly twice; the set of coverings equals the exact enumeration |
| `aztec_sample.py 2 0 2 3`, `aztec_sample.py 3 1 2 3` | perpendicular-partnered: 2 (not 16) and 160 (not 128) ground states |
| `chordflip.py` | tile flip on a chord face: 2 islands, barrier 2 J1 |
| `fig_shakti_aztec.py` | Fig. 10 (file fig17) |
| `fig_micro_compose.py` | Fig. 11 (file fig16), from `process_figures/` drawn by `draw.py` |
| `path.py 8 ../results/sol_11_8.pkl`, `path.py 8 ../results/sol_22_8.pkl`, `flip.py` | barriers of pair creation, two-pair creation and the chord-free tile flip (Table II) |
| `cp2.py`, `path.py`, `avail.py`, ... | elementary processes and availability probabilities (see `code/README.md`) |

J1, J2 are passed as integers (e.g. 2 3 for r = 1.5).
