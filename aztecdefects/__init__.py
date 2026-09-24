"""Exact computations for diagonal (B/W) defect tiles in Aztec-diamond domino tilings.

Companion code for C. Won, "Composite diagonal defects in Aztec-diamond dimer
coverings: boundary pinning, Coulomb attraction, and exact asymptotics".

Modules
-------
lattice     : cells, colours, diagonal partners, Aztec-diamond regions, boundary edges
modular     : primes, modular Pfaffian / inverse / interpolation, CRT
pfaffian    : Pfaffian orientation of the non-bipartite covering graph, exact
              generating polynomials Z_n(u), local statistics, float evaluations
kasteleyn   : bipartite (pure-domino) Kasteleyn matrices, exact k=1 pair weights
              (float and multimodular banded elimination)
boundary    : the exact boundary formula (Theorem 1), the determinant formula
              (Theorem 2) and their large-n evaluation
bulk        : periodic model on the torus: spectral gap, bulk pair weight,
              bulk densities and thermodynamics
shuffling   : generalized urban renewal for faces with a chord, face classes
thermo      : free-energy landscape, thresholds T_k, k*(T)
bruteforce  : direct enumeration used for validation
"""
__version__ = "1.0.0"
