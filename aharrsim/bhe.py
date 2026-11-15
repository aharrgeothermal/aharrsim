"""Borehole thermal resistances of a single U-tube in a grouted borehole.

Line-source (zeroth-order multipole) formulas of Hellström (1991) / Claesson &
Hellström (2011) for two identical pipes at +-x_c from the borehole axis:

    R_b  = 1/(4 pi k_b) [ ln(r_b/r_p) + ln(r_b/(2 x_c)) + sigma ln(r_b^4 / (r_b^4 - x_c^4)) ] + R_p / 2
    R_a  = 1/(pi k_b)   [ ln(2 x_c/r_p) + sigma ln((r_b^2 + x_c^2) / (r_b^2 - x_c^2)) ] + 2 R_p
    sigma = (k_b - k_s) / (k_b + k_s)

R_b is the effective borehole resistance (both pipes at the fluid temperature to
the borehole wall), R_a the internal resistance between the legs, R_p the pipe
wall plus fluid film.  Converted to the delta circuit used by the solver:
leg -> wall  R_1b = R_2b = 2 R_b,  leg <-> leg  R_12 = 4 R_b R_a / (4 R_b - R_a).
Checked against pygfunction's multipole method in validation/verify.py (V8).
"""
from __future__ import annotations

import math

import numpy as np


def pipe_resistance(r_in, r_out, k_pipe, h_film):
    """Pipe wall + fluid film resistance per metre [K m/W]."""
    return math.log(r_out / r_in) / (2 * math.pi * k_pipe) + 1.0 / (2 * math.pi * r_in * np.asarray(h_film, float))


def utube_delta(r_b, r_p_out, x_c, k_grout, k_soil, R_p):
    """Delta-circuit resistances (R_leg_wall, R_leg_leg) per metre for a single U-tube.
    x_c = half the shank spacing; R_p may be an array (per cell)."""
    R_p = np.asarray(R_p, float)
    r_b = np.asarray(r_b, float)
    k_soil = np.asarray(k_soil, float)
    sigma = (k_grout - k_soil) / (k_grout + k_soil)
    R_b = (np.log(r_b / r_p_out) + np.log(r_b / (2 * x_c)) + sigma * np.log(r_b ** 4 / (r_b ** 4 - x_c ** 4))) \
        / (4 * math.pi * k_grout) + 0.5 * R_p
    R_a = (np.log(2 * x_c / r_p_out) + sigma * np.log((r_b ** 2 + x_c ** 2) / (r_b ** 2 - x_c ** 2))) \
        / (math.pi * k_grout) + 2.0 * R_p
    R_lw = 2.0 * R_b
    denom = 4.0 * R_b - R_a
    R_ll = np.where(denom > 1e-9, 4.0 * R_b * R_a / np.maximum(denom, 1e-9), np.inf)
    return R_lw, R_ll, R_b, R_a
