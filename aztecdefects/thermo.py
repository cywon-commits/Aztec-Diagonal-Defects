"""Thermodynamics from the exact counts c_k (paper Sec. VII).

Each diagonal tile costs eps_d more than a domino, so E_k = 2 k eps_d,
S(k) = ln c_k, F(k, T) = 2 k eps_d - T S(k), u = exp(-2 eps_d / T).
Temperatures are in units of eps_d / k_B.
"""
import math
import numpy as np


def entropy(c):
    return np.array([math.log(x) for x in c if x > 0])


def free_energy_landscape(c, T):
    S = entropy(c); k = np.arange(len(S))
    return 2 * k - T * S


def k_star(c, T):
    """Minimiser of F(k, T)."""
    return int(np.argmin(free_energy_landscape(c, T)))


def k_mean(c, T):
    S = entropy(c); k = np.arange(len(S)); lw = S - 2 * k / T; w = np.exp(lw - lw.max())
    return float((w * k).sum() / w.sum())


def thresholds(c):
    """Delta_k = ln(c_k/c_{k-1}) and T_k = 2/Delta_k (k = 1, 2, ...), and concavity of ln c_k."""
    S = entropy(c); d = np.diff(S)
    return d, 2 / d, bool(np.all(np.diff(d) <= 0))


def onset_temperature_asymptotic(n):
    """T* = 2 eps_d / ln(c_1/c_0) with ln(c_1/c_0) ~ 2 n ln2 - 3 ln n + ln(2/pi)."""
    return 2 / (2 * n * math.log(2) - 3 * math.log(n) + math.log(2 / math.pi))
