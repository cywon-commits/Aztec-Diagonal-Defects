#!/usr/bin/env bash
# Reproduce all data and figures of the paper (single core, roughly 1 hour).
set -e
cd "$(dirname "$0")"
python -m pytest -q tests
python scripts/01_polynomials.py --nmax 17
python scripts/02_n9_statistics.py
python scripts/03_bulk.py
python scripts/04_equal_energy.py --nmax 61
python scripts/05_k1_large_n.py --n 21 41 61 81 101 121
python scripts/06_boundary_asymptotics.py --nmax 1601
python scripts/07_k2.py
python scripts/08_thermodynamics.py
python scripts/09_schematics.py
python scripts/10_verify_proofs.py --nmax 10
python scripts/11_placement_B.py
