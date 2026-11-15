"""Analytical reference solutions used for verification.

1. reverse_circulation_steady
   Steady counter-flow with fluid DOWN the annulus and UP the tubing, linear
   geotherm, constant fluid properties and constant thermal resistances
   (formation R_f and tubing wall R_t, per unit length).  This is the
   reverse-circulation form of the Holmes & Swift (1970) / Kabir et al. (1996)
   solution; derived here from

       w c dTa/dz =  (Tf - Ta)/R_f + (Tt - Ta)/R_t
       w c dTt/dz =  (Tt - Ta)/R_t          (z = depth, tubing flows upward)
       Ta(0) = T_in ,  Tt(L) = Ta(L) ,  Tf = Ts + g z

   With a = w c R_f, b = w c R_t:
       Tt = A exp(l1 z) + B exp(l2 z) + Ts + g z + b g
       Ta = (1 - b l1) A exp(l1 z) + (1 - b l2) B exp(l2 z) + Ts + g z
       l1,2 = [-1 +/- sqrt(1 + 4 a / b)] / (2 a)

2. jaeger_cylinder_flux
   Heat flow per unit length into a cylinder of radius a held at constant
   temperature in an infinite medium initially at uniform temperature
   (Carslaw & Jaeger 1959, Sec. 13.5):
       q' = (8 k dT / pi) * int_0^inf exp(-tau u^2) / (u [J0(u)^2 + Y0(u)^2]) du,
       tau = alpha t / a^2

3. hasan_kabir_f / kabir_reverse_transient
   Quasi-steady formation resistance with the Hasan & Kabir (1994)
   dimensionless time function, giving the classical semi-analytical
   transient approximation for reverse circulation (Kabir et al. 1996 form).
"""
from __future__ import annotations

import math
import numpy as np
from scipy import special, integrate


def reverse_circulation_steady(z, L, w, c, R_f, R_t, T_in, T_s, g):
    """Return (Ta(z), Tt(z)) for the steady reverse-circulation problem."""
    a = w * c * R_f
    b = w * c * R_t
    disc = math.sqrt(1.0 + 4.0 * a / b)
    l1 = (-1.0 + disc) / (2.0 * a)      # > 0
    l2 = (-1.0 - disc) / (2.0 * a)      # < 0
    # unknowns A' = A exp(l1 L) and B; exponentials written to avoid overflow
    M = np.array([[(1.0 - b * l1) * math.exp(-l1 * L), (1.0 - b * l2)],
                  [l1, l2 * math.exp(l2 * L)]])
    rhs = np.array([T_in - T_s, -g])
    Ap, B = np.linalg.solve(M, rhs)
    z = np.asarray(z, float)
    e1 = np.exp(l1 * (z - L))
    e2 = np.exp(l2 * z)
    Tt = Ap * e1 + B * e2 + T_s + g * z + b * g
    Ta = (1.0 - b * l1) * Ap * e1 + (1.0 - b * l2) * B * e2 + T_s + g * z
    return Ta, Tt


def forward_circulation_steady(z, L, w, c, R_f, R_t, T_in, T_s, g):
    """Return (Ta(z), Tt(z)) for the steady forward-circulation problem: fluid
    injected down the centre tubing, returning up the annulus (the annulus
    exchanges with the formation at T_s + g z through R_f, the two streams
    through R_t).  Same method as reverse_circulation_steady.

        w c dTt/dz = (Ta - Tt) / R_t
       -w c dTa/dz = (Tf - Ta) / R_f + (Tt - Ta) / R_t
        Tt(0) = T_in,  Ta(L) = Tt(L)
    """
    a = w * c * R_f
    b = w * c * R_t
    disc = math.sqrt(1.0 + 4.0 * a / b)
    l1 = (1.0 + disc) / (2.0 * a)       # > 0
    l2 = (1.0 - disc) / (2.0 * a)       # < 0
    # Tt = A e^{l1 z} + B e^{l2 z} + Ts + g z - b g ;  Ta = (1 + b l1) A e^{l1 z} + (1 + b l2) B e^{l2 z} + Ts + g z
    # unknowns A' = A exp(l1 L) and B
    M = np.array([[math.exp(-l1 * L), 1.0],
                  [l1, l2 * math.exp(l2 * L)]])
    rhs = np.array([T_in - T_s + b * g, -g])
    Ap, B = np.linalg.solve(M, rhs)
    z = np.asarray(z, float)
    e1 = np.exp(l1 * (z - L))
    e2 = np.exp(l2 * z)
    Tt = Ap * e1 + B * e2 + T_s + g * z - b * g
    Ta = (1.0 + b * l1) * Ap * e1 + (1.0 + b * l2) * B * e2 + T_s + g * z
    return Ta, Tt


def _jaeger_integrand(v, tau):
    u = math.exp(v)
    return math.exp(-tau * u * u) / (special.j0(u) ** 2 + special.y0(u) ** 2)


def jaeger_dimensionless_flux(tau):
    """q' / (2 pi k dT) for the constant-temperature cylinder (Carslaw & Jaeger)."""
    tau = float(tau)
    # integrate in v = ln u; analytic tail for v -> -inf where
    # J0^2 + Y0^2 ~ (4/pi^2) (v - ln2 + gamma)^2
    v_lo = -60.0
    v_hi = 0.5 * math.log(50.0 / tau) if tau < 50.0 else math.log(math.sqrt(50.0 / tau))
    v_hi = max(v_hi, v_lo + 5.0)
    pts = np.linspace(v_lo, v_hi, 40)
    total = 0.0
    for a_, b_ in zip(pts[:-1], pts[1:]):
        val, _ = integrate.quad(_jaeger_integrand, a_, b_, args=(tau,), limit=200, epsabs=0, epsrel=1e-11)
        total += val
    gam = 0.5772156649015329
    tail = (math.pi ** 2 / 4.0) / abs(v_lo - math.log(2.0) + gam)
    integral = total + tail
    q_over = (8.0 / math.pi) * integral          # q' / (k dT)
    return q_over / (2.0 * math.pi)


def jaeger_cylinder_flux(t, k, alpha, a, dT):
    """Heat flow per unit length [W/m] into a cylinder held dT below the far field."""
    tau = alpha * t / a ** 2
    return 2.0 * math.pi * k * dT * jaeger_dimensionless_flux(tau)


def hasan_kabir_f(tD):
    """Hasan & Kabir (1994) dimensionless temperature function f(tD)."""
    tD = np.asarray(tD, float)
    short = 1.1281 * np.sqrt(tD) * (1.0 - 0.3 * np.sqrt(tD))
    long_ = (0.4063 + 0.5 * np.log(tD)) * (1.0 + 0.6 / tD)
    return np.where(tD <= 1.5, short, long_)


def kabir_reverse_transient(t, z, L, w, c, k_e, alpha_e, r_w, U_w, R_t, T_in, T_s, g):
    """Semi-analytical reverse-circulation temperatures at time t (Kabir et al. 1996 form).

    The formation resistance per unit length is (1/(2 pi r_w U_w) + f(tD)/(2 pi k_e)),
    with U_w the wellbore (film + completion) coefficient referred to r_w.
    """
    tD = alpha_e * t / r_w ** 2
    R_f = 1.0 / (2 * math.pi * r_w * U_w) + float(hasan_kabir_f(tD)) / (2 * math.pi * k_e)
    return reverse_circulation_steady(z, L, w, c, R_f, R_t, T_in, T_s, g)
