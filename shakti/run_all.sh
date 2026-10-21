#!/usr/bin/env bash
# Shakti vertex-model calculations of Sec. VIII (needs ortools). Exhaustive enumerations: ~1-2 hours on 2 cores.
set -e
cd "$(dirname "$0")/code"
python multiplicity.py 8                 # every covering <-> 2 island configurations (periodic 8x8)
python multiplicity.py 12                # same on 12x12
python aztec_sample.py 3 0 2 3           # A_3 placement A as Shakti sample, r = 1.5
python aztec_sample.py 3 0 5 6           # r = 1.2
python aztec_sample.py 3 0 5 9           # r = 1.8
python aztec_sample.py 2 1 2 3           # A_2 placement B
python aztec_sample.py 2 0 2 3           # perpendicular-partnered counterexamples
python aztec_sample.py 3 1 2 3
python opentiling.py 3 0                 # island states -> coverings, compared with exact enumeration
python opentiling.py 2 1
python chordflip.py                      # chord-face tile flip: 2 islands, barrier 2 J1
python flip.py                           # chord-free face flip: 4 islands, barrier 2 J2
python avail.py                          # availability probabilities 1/(4 pi^2), 1/pi^4
python fig_micro_compose.py              # Fig. 16 (from the process panels drawn by draw.py)
python fig_shakti_aztec.py               # Fig. 17
# elementary processes of Table II (creation, hop, sign flip; results/res12.pkl, log12.txt):
#   python cp2.py 12 '<queries>' ../results/res12.pkl    (queries in code/README.md)
