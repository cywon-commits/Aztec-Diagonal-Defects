# Code

Requires Python 3 with numpy, matplotlib and ortools (`pip install ortools`).
Run every script from this `code/` directory: the scripts load each other with relative paths.
Pickle files live in `../results/`.

| file | role |
|---|---|
| `v.py` | Island configuration of the user's 4×4 figure (4 O → 4 L pinwheel); the `R` state is the periodic L ground-state unit cell. |
| `shakti.py` | Periodic Shakti torus (size P) built from that unit cell; vertex/arm tables. Used as a module by the other scripts (only the part before `J1=1.0; J2=` is loaded). |
| `cp.py` | CP-SAT: minimal island flips from any ground state to a state with kH horizontal and kV vertical O's. `python3 cp.py 8 1 1 200` |
| `cp2.py` | CP-SAT with prescribed O positions and charge signs (creation, hops, sign flips), plus barriers. Queries are passed as JSON (see the example below). |
| `path.py` | Minimax barrier along the flipped-island set of a saved (X,Y) pair. `python3 path.py 8 ../results/sol_11_8.pkl` |
| `flip.py` | Smallest move between two defect-free ground states (reference barrier 2J2). |
| `cyc.py`, `avail.py` | Chessboard alternating cycles; availability probabilities in the uniform bulk tiling via the infinite-plane K^{-1} (1/(4π²), 1/π⁴). |
| `cp_open.py` | Open samples: E0, minimal energy with O's at given blocks, minimal flips and barrier. `python3 cp_open.py 0 0 8 8 2 3 '{"edge":[[2,0],[2,1]]}' out.pkl` (J1=2, J2=3, i.e. r=1.5) |
| `sitemap.py` | Cost of an O at every boundary block of an L×L open sample. `python3 sitemap.py 12 2 3` |
| `draw.py` | Figures. `python3 draw.py P file.pkl x0 y0 x1 y1 out.png "after title" "before title"` |

Coordinate convention: block (a,b) = [2a,2a+2]×[2b,2b+2] is a horizontal pair when a+b is odd.
This follows the user's figure; the companion vorticity paper uses the opposite parity.

Example `cp2.py` query (12×12 torus; block coordinates; sign 0 = either, ±1 = charge sign):

```
python3 cp2.py 12 '{"hop_diag":[{"1,2":0,"4,4":0},{"2,3":0,"4,4":0}]}' out.pkl
```
