"""Friction, heat-transfer and two-phase flow correlations.

References
----------
Churchill, S.W. (1977) Friction-factor equation spans all fluid-flow regimes.
    Chem. Eng. 84(24):91-92.
Gnielinski, V. (1976) New equations for heat and mass transfer in turbulent
    pipe and channel flow. Int. Chem. Eng. 16:359-368.
Gnielinski, V. (2009) Heat transfer coefficients for turbulent flow in
    concentric annular ducts. Heat Transfer Eng. 30(6):431-436.
Gnielinski, V. (2013) On heat transfer in tubes. Int. J. Heat Mass Transfer
    63:134-140 (laminar-turbulent transition interpolation).
Chen, J.C. (1966) Correlation for boiling heat transfer to saturated fluids in
    convective flow. I&EC Process Des. Dev. 5(3):322-329 (two-phase factor F).
Zuber, N., Findlay, J.A. (1965) Average volumetric concentration in two-phase
    flow systems. J. Heat Transfer 87:453-468.
Hasan, A.R., Kabir, C.S. (2010) Modeling two-phase fluid and heat flows in
    geothermal wells. J. Pet. Sci. Eng. 71:77-86 (drift-flux parameters).
"""
from __future__ import annotations

import numpy as np

G_ACC = 9.80665


# ----------------------------------------------------------------------------
# Friction
# ----------------------------------------------------------------------------
def friction_churchill(Re, rel_rough):
    """Darcy friction factor, all regimes (Churchill 1977)."""
    Re = np.maximum(np.asarray(Re, float), 1e-6)
    A = (2.457 * np.log(1.0 / ((7.0 / Re) ** 0.9 + 0.27 * rel_rough))) ** 16
    B = (37530.0 / Re) ** 16
    return 8.0 * ((8.0 / Re) ** 12 + (A + B) ** (-1.5)) ** (1.0 / 12.0)


# ----------------------------------------------------------------------------
# Single-phase heat transfer
# ----------------------------------------------------------------------------
def _nu_turb_pipe(Re, Pr):
    f = (0.790 * np.log(Re) - 1.64) ** -2
    return (f / 8.0) * (Re - 1000.0) * Pr / (1.0 + 12.7 * np.sqrt(f / 8.0) * (Pr ** (2.0 / 3.0) - 1.0))


def nusselt_pipe(Re, Pr):
    """Fully developed Nusselt number in a circular pipe.

    Laminar 4.36 (uniform flux); Gnielinski (1976) turbulent; linear
    interpolation 2300 < Re < 1e4 as recommended by Gnielinski (2013).
    """
    Re = np.atleast_1d(np.asarray(Re, float))
    Pr = np.atleast_1d(np.asarray(Pr, float)) * np.ones_like(Re)
    nu_lam = 4.364
    nu = np.full(Re.shape, nu_lam)
    turb = Re >= 1.0e4
    nu[turb] = _nu_turb_pipe(Re[turb], Pr[turb])
    tr = (Re > 2300.0) & ~turb
    if np.any(tr):
        g = (Re[tr] - 2300.0) / (1.0e4 - 2300.0)
        nu[tr] = (1.0 - g) * nu_lam + g * _nu_turb_pipe(1.0e4, Pr[tr])
    return np.maximum(nu, nu_lam)


def _nu_turb_annulus(Re, Pr, a, wall):
    lna = np.log(a)
    Re_star = Re * ((1.0 + a ** 2) * lna + (1.0 - a ** 2)) / ((1.0 - a) ** 2 * lna)
    zeta = (1.8 * np.log10(Re_star) - 1.5) ** -2
    k1 = 1.07 + 900.0 / Re - 0.63 / (1.0 + 10.0 * Pr)
    if wall == "inner":
        F = 0.75 * a ** -0.17
    else:
        F = 0.9 - 0.15 * a ** 0.6
    return (zeta / 8.0) * Re * Pr / (k1 + 12.7 * np.sqrt(zeta / 8.0) * (Pr ** (2.0 / 3.0) - 1.0)) * F


def nusselt_annulus(Re, Pr, a, wall):
    """Nusselt number (based on hydraulic diameter) at one wall of a concentric annulus.

    a = d_inner / d_outer; wall = 'inner' (tubing outer surface) or 'outer'
    (casing inner surface).  Gnielinski (2009) for turbulent flow, VDI
    laminar values (one wall heated, the other adiabatic) and linear
    interpolation for 2300 < Re < 1e4.
    """
    Re = np.atleast_1d(np.asarray(Re, float))
    Pr = np.atleast_1d(np.asarray(Pr, float)) * np.ones_like(Re)
    a = np.atleast_1d(np.asarray(a, float)) * np.ones_like(Re)
    a = np.clip(a, 1e-3, 0.999)
    if wall == "inner":
        nu_lam = 3.66 + 1.2 * a ** -0.8
    else:
        nu_lam = 3.66 + 1.2 * a ** 0.5
    nu = nu_lam.copy()
    turb = Re >= 1.0e4
    if np.any(turb):
        nu[turb] = _nu_turb_annulus(Re[turb], Pr[turb], a[turb], wall)
    tr = (Re > 2300.0) & ~turb
    if np.any(tr):
        g = (Re[tr] - 2300.0) / (1.0e4 - 2300.0)
        nu[tr] = (1.0 - g) * nu_lam[tr] + g * _nu_turb_annulus(1.0e4, Pr[tr], a[tr], wall)
    return np.maximum(nu, nu_lam)


def chen_factor(x, rho_l, rho_g, mu_l, mu_g):
    """Chen (1966) convective two-phase enhancement factor F."""
    x = np.clip(np.asarray(x, float), 1e-9, 1 - 1e-9)
    inv_xtt = (x / (1.0 - x)) ** 0.9 * (rho_l / rho_g) ** 0.5 * (mu_g / mu_l) ** 0.1
    F = np.where(inv_xtt <= 0.1, 1.0, 2.35 * (inv_xtt + 0.213) ** 0.736)
    return np.maximum(F, 1.0)


# ----------------------------------------------------------------------------
# Two-phase hydraulics (vertical / inclined upward flow)
# ----------------------------------------------------------------------------
def _void_drift_flux(G, D, x, rho_l, rho_g, sigma, cos_theta):
    """Zuber-Findlay drift-flux void fraction with flow-pattern dependent
    parameters following Hasan & Kabir (2010): bubbly (C0=1.2, Harmathy rise
    velocity), slug/churn (C0=1.2->1.15, Taylor-bubble rise velocity) and a
    smooth transition to no-slip in annular flow (C0=1, Vd=0)."""
    drho = np.maximum(rho_l - rho_g, 1e-6)
    jg = G * x / rho_g
    jl = G * (1.0 - x) / rho_l
    j = jg + jl
    ct = np.clip(cos_theta, 0.0, 1.0)
    st = np.sqrt(np.maximum(1.0 - ct ** 2, 0.0))
    v_bub = 1.53 * (G_ACC * sigma * drho / rho_l ** 2) ** 0.25 * np.sqrt(ct + 1e-12)
    v_tay = 0.35 * np.sqrt(G_ACC * D * drho / rho_l) * np.sqrt(ct + 1e-12) * (1.0 + st) ** 1.2
    v_ann = 3.1 * (G_ACC * sigma * drho / rho_g ** 2) ** 0.25
    # bubbly estimate
    a_b = jg / (1.2 * j + v_bub)
    # slug / churn estimate
    a_s = jg / (1.15 * j + v_tay)
    w_s = np.clip((a_b - 0.20) / 0.10, 0.0, 1.0)          # bubbly -> slug around alpha 0.25
    alpha = (1.0 - w_s) * a_b + w_s * a_s
    # annular: blend to homogeneous (no slip)
    a_h = jg / np.maximum(j, 1e-12)
    w_a = np.clip((jg / v_ann - 0.8) / 0.4, 0.0, 1.0)
    alpha = (1.0 - w_a) * alpha + w_a * a_h
    return np.clip(alpha, 0.0, 1.0 - 1e-9)


def two_phase_gradient(model, G, D, x, rho_l, rho_g, mu_l, mu_g, sigma, rough, cos_theta):
    """Return (rho_grav, fric_grad, rho_mom, alpha) for two-phase upward flow.

    rho_grav  : in-situ mixture density used in the hydrostatic term [kg/m3]
    fric_grad : frictional pressure gradient magnitude [Pa/m]
    rho_mom   : momentum density, momentum flux = G**2 / rho_mom [kg/m3]
    alpha     : void fraction [-]
    model     : 'hem' (homogeneous equilibrium) or 'drift_flux'
    """
    x = np.clip(np.asarray(x, float), 0.0, 1.0)
    v_ns = x / rho_g + (1.0 - x) / rho_l
    rho_ns = 1.0 / v_ns
    mu_m = 1.0 / (x / mu_g + (1.0 - x) / mu_l)                 # McAdams
    Re = np.abs(G) * D / mu_m
    f = friction_churchill(Re, rough / D)
    fric = f * G * G / (2.0 * D * rho_ns)
    if model == "hem":
        alpha = x * rho_l / (x * rho_l + (1.0 - x) * rho_g)
        return rho_ns, fric, rho_ns, alpha
    alpha = _void_drift_flux(np.abs(G), D, x, rho_l, rho_g, sigma, cos_theta)
    alpha = np.maximum(alpha, 0.0)
    rho_m = alpha * rho_g + (1.0 - alpha) * rho_l
    with np.errstate(divide="ignore", invalid="ignore"):
        inv_rm = np.where(x > 0, x ** 2 / (np.maximum(alpha, 1e-12) * rho_g), 0.0) + \
            (1.0 - x) ** 2 / (np.maximum(1.0 - alpha, 1e-12) * rho_l)
    return rho_m, fric, 1.0 / inv_rm, alpha
