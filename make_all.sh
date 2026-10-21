#!/usr/bin/env bash
# Reproduce all data and figures of the paper (single core, roughly 1.5 hours with the cached data in data/).
# The Shakti vertex-model calculations (Sec. VIII) are run separately: bash shakti/run_all.sh
set -e
cd "$(dirname "$0")"
python -m pytest -q tests
python scripts/01_polynomials.py --nmax 17
python scripts/02_n9_statistics.py
python scripts/03_bulk.py
python scripts/04_equal_energy.py --nmax 81
python scripts/05_k1_large_n.py --n 21 41 61 81 101 121
python scripts/06_boundary_asymptotics.py --nmax 1601
python scripts/07_k2.py
python scripts/08_thermodynamics.py
python scripts/09_schematics.py
python scripts/10_verify_proofs.py --nmax 10
python scripts/13_placement_local.py
python scripts/11_placement_B.py
python scripts/12_audit_checks.py
