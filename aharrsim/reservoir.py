"""Reservoir / heat-source elements for the development types that are not pure wellbore
conduction.  Each element gives a production (or return) temperature history; the well
models in aharrsim.simulator take it as a callable T(t).

ConstantReservoir   convective hydrothermal reservoir: T_res, optional linear decline.
DoubletAquifer      cold-water doublet in a confined aquifer (Gringarten & Sauty 1975):
                    streamlines of the source-sink pair, Lauwerier (1955) heat exchange
                    with the confining rock along every streamline.
ParallelFractures   EGS / HDR doublet through N equidistant parallel fractures
                    (Gringarten, Witherspoon & Ohnishi 1975), inversion of the Laplace-domain
                    solution by the Fourier-series method of de Hoog et al. (1982); the same
                    formulation as GEOPHIRES' MPF model.
SingleWellDipole    one well with an injection and a production screen a distance D apart in a porous or
                    densely fractured block (single-well EGS): streamlines of the point source-sink pair,
                    a sharp front or the parallel-fracture response along every streamline.
CyclicFracture      single-well huff-and-puff through one fracture (inject - soak - produce):
                    plug flow along the fracture, 1-D conduction into the rock faces.
HeatStore           lumped, finite heated zone (man-made heat after thermal EOR), swept by a
                    well-mixed doublet, losing heat to the surrounding rock by 1-D conduction.

Verification: validation/verify.py (V10): Lauwerier and erfc limits, loss-free doublet
breakthrough time pi H D^2 (rho c)_a / (3 Q (rho c)_w), inversion vs mpmath (including sharp fronts,
where Gaver-Stehfest smears the response), cyclic model vs Lauwerier in continuous injection,
energy conservation.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.special import erfc

RHOC_W = 4.18e6        # J/m3/K, water near 20-100 C (overridable)


def cbar_w(T_inj_C, T_rock_C, p_Pa=30.0e6):
    """Mean specific heat of water between the injection and the rock temperature [J/kg/K],
    (h(T_rock) - h(T_inj)) / (T_rock - T_inj) at pressure p.  The heat a kilogram of water takes from the rock is an
    enthalpy difference, so this is the heat capacity that sets the speed of a thermal front and the front-crossing
    time t_s = V (rho c)_r / (mdot cbar_w).  Use rhoc_w = cbar_w * rho with Q = mdot / rho in the elements below."""
    import CoolProp.CoolProp as CP
    if abs(T_rock_C - T_inj_C) < 1e-6:
        return float(CP.PropsSI("C", "T", T_rock_C + 273.15, "P", p_Pa, "Water"))
    h1 = CP.PropsSI("H", "T", T_rock_C + 273.15, "P", p_Pa, "Water")
    h0 = CP.PropsSI("H", "T", T_inj_C + 273.15, "P", p_Pa, "Water")
    return float((h1 - h0) / (T_rock_C - T_inj_C))


# ----------------------------------------------------------------------------
class ConstantReservoir:
    """Produced temperature T_res [K], optionally declining linearly (K per year) from t = t_decline_start."""

    def __init__(self, T_res_C, decline_K_per_yr=0.0, T_floor_C=None):
        self.T0 = T_res_C + 273.15
        self.decline = decline_K_per_yr / (365.25 * 86400.0)
        self.T_floor = (T_floor_C + 273.15) if T_floor_C is not None else -np.inf

    def __call__(self, t):
        return float(max(self.T0 - self.decline * t, self.T_floor))

    T_out = None   # not an in-series element


# ----------------------------------------------------------------------------
def stehfest_coefficients(N=14):
    V = np.zeros(N)
    for i in range(1, N + 1):
        s = 0.0
        for k in range((i + 1) // 2, min(i, N // 2) + 1):
            s += (k ** (N // 2) * math.factorial(2 * k)) / (math.factorial(N // 2 - k) * math.factorial(k) * math.factorial(k - 1)
                                                           * math.factorial(i - k) * math.factorial(2 * k - i))
        V[i - 1] = (-1) ** (i + N // 2) * s
    return V


_STEH = stehfest_coefficients(14)


def stehfest(fbar, t, V=_STEH):
    """Gaver-Stehfest inverse Laplace transform of fbar(s) at times t (array)."""
    t = np.asarray(t, float)
    out = np.zeros_like(t)
    ln2 = math.log(2.0)
    for i, Vi in enumerate(V, start=1):
        out += Vi * fbar(i * ln2 / t)
    return out * ln2 / t


def dehoog(fbar, t, M=24, tol=1e-9):
    """Inverse Laplace transform by the accelerated Fourier-series method of de Hoog, Knight & Stokes (1982), vectorised over
    the times t (every time gets its own contour, half period 2 t). fbar must accept a complex array of shape (len(t), 2 M + 1).
    Unlike Gaver-Stehfest it keeps sharp fronts sharp; against 25-digit mpmath it is good to about 1e-7 for the fracture solution."""
    t = np.atleast_1d(np.asarray(t, float)); nt = len(t)
    T = 2.0 * t
    gam = -math.log(tol) / (2.0 * T)
    s = gam[:, None] + 1j * math.pi * np.arange(2 * M + 1)[None, :] / T[:, None]
    a = np.asarray(fbar(s), complex).copy(); a[:, 0] *= 0.5
    e = np.zeros((nt, 2 * M + 1, M + 1), complex); q = np.zeros((nt, 2 * M, M + 1), complex)
    q[:, :, 1] = a[:, 1:2 * M + 1] / a[:, 0:2 * M]
    for r in range(2, M + 2):                                   # quotient-difference table
        n1 = 2 * (M - r + 1) + 1
        e[:, 0:n1, r - 1] = q[:, 1:n1 + 1, r - 1] - q[:, 0:n1, r - 1] + e[:, 1:n1 + 1, r - 2]
        if r < M + 1:
            n2 = 2 * (M - r) + 2
            q[:, 0:n2, r] = q[:, 1:n2 + 1, r - 1] * e[:, 1:n2 + 1, r - 1] / e[:, 0:n2, r - 1]
    d = np.zeros((nt, 2 * M + 1), complex)
    d[:, 0] = a[:, 0]; d[:, 1:2 * M:2] = -q[:, 0, 1:M + 1]; d[:, 2:2 * M + 1:2] = -e[:, 0, 1:M + 1]
    z = np.exp(1j * math.pi * t / T)
    A = np.zeros((2 * M + 2, nt), complex); B = np.zeros((2 * M + 2, nt), complex)
    A[1] = d[:, 0]; B[0] = 1.0; B[1] = 1.0
    for n in range(3, 2 * M + 2):                               # continued fraction
        A[n - 1] = A[n - 2] + d[:, n - 2] * z * A[n - 3]
        B[n - 1] = B[n - 2] + d[:, n - 2] * z * B[n - 3]
    h = 0.5 * (1.0 + (d[:, 2 * M - 1] - d[:, 2 * M]) * z)       # remainder estimate (the "double acceleration")
    rem = -h * (1.0 - np.sqrt(1.0 + d[:, 2 * M] * z / h ** 2))
    A[2 * M + 1] = A[2 * M] + rem * A[2 * M - 1]
    B[2 * M + 1] = B[2 * M] + rem * B[2 * M - 1]
    return np.exp(gam * t) / T * np.real(A[2 * M + 1] / B[2 * M + 1])


class ParallelFractures:
    """Gringarten-Witherspoon-Ohnishi (1975) heat extraction through N parallel equidistant fractures.

    Q         total volumetric flow [m3/s]          N  number of fractures
    L_flow    flow length within each fracture [m]  B  fracture breadth perpendicular to the flow [m]
    S         fracture spacing [m]                  k_r, rhoc_r rock conductivity / volumetric heat capacity
    T_rock    initial rock temperature [C]          T_inj injection temperature [C] (can be a callable of time)
    Dimensionless: t_D = t (rhoc_w q_f)^2 / (4 k_r rhoc_r (B L)^2),  x_ED = S rhoc_w q_f / (4 B L k_r),  q_f = Q/N,
    T_D(s) = exp(-sqrt(s) tanh(x_ED sqrt(s))) / s.  Produced T = T_rock - T_D (T_rock - T_inj).
    """

    def __init__(self, Q, N, L_flow, B, S, k_r, rhoc_r, T_rock_C, T_inj_C, rhoc_w=RHOC_W, impedance_Pa_per_m3s=None):
        self.q_f = Q / N
        self.N, self.L, self.B, self.S = N, L_flow, B, S
        self.k_r, self.rhoc_r, self.rhoc_w = k_r, rhoc_r, rhoc_w
        self.T_rock = T_rock_C + 273.15
        self.T_inj = T_inj_C
        self.gamma2 = (2.0 * B * L_flow * math.sqrt(k_r * rhoc_r) / (rhoc_w * self.q_f)) ** 2   # t = t_D * gamma2
        self.x_ED = S * rhoc_w * self.q_f / (4.0 * B * L_flow * k_r)
        self.Q = Q
        self.impedance = impedance_Pa_per_m3s
        self._cache_t = None

    def T_D(self, t, method="dehoog"):
        """Dimensionless temperature drop at the fracture outlet, 0 at t = 0 -> 1 at long time.
        method 'dehoog' (default) or 'stehfest' (N = 14; versions up to 2.1 used it: it smears the front when the fractures
        are closer than the conduction length, by 0.12 in T_D for fractures 10 m apart and by up to 0.28 on a sharp front)."""
        t = np.atleast_1d(np.asarray(t, float))
        tD = np.maximum(t / self.gamma2, 1e-12)
        xe = self.x_ED
        if method == "stehfest":
            f = lambda s: np.exp(-np.sqrt(s) * np.tanh(np.minimum(xe * np.sqrt(s), 50.0))) / s
            return np.clip(stehfest(f, tD), 0.0, 1.0)

        def f(s):
            r = np.sqrt(s)
            return np.exp(-r * np.tanh(xe * r)) / s

        with np.errstate(all="ignore"):
            out = dehoog(f, tD)
        # the transform underflows where the response is still zero: before the front, or early on an isolated fracture
        out = np.where(np.isfinite(out), out, 0.0)
        return np.clip(out, 0.0, 1.0)

    @property
    def t_sweep(self):
        """Time for the cold front to cross the block if all the heat of the rock were taken [s]: V (rho c)_r / (Q (rho c)_w)."""
        return self.S * self.L * self.B * self.rhoc_r / (self.q_f * self.rhoc_w)

    def groups(self, t):
        """The two groups that fix the response: x = t / t_sweep and y = (S / 2) / sqrt(alpha t)."""
        t = np.asarray(t, float)
        return t / self.t_sweep, 0.5 * self.S / np.sqrt(self.k_r / self.rhoc_r * t)

    def T_out(self, T_in_K, t):
        """Outlet temperature [K] for the current injection temperature (quasi-steady use of the step response)."""
        TD = float(self.T_D(max(t, 1.0))[0])
        return self.T_rock - TD * (self.T_rock - T_in_K)

    def __call__(self, t):
        Tin = self.T_inj(t) if callable(self.T_inj) else self.T_inj
        return self.T_out(Tin + 273.15, t)

    def pumping_dp(self):
        return self.Q * self.impedance if self.impedance else 0.0


# ----------------------------------------------------------------------------
class DoubletAquifer:
    """Gringarten & Sauty (1975) doublet in a confined aquifer: production temperature after cold-water
    reinjection, with conductive heat exchange with the cap and base rock along every streamline.

    Q [m3/s] circulated, D [m] well spacing, H [m] aquifer thickness, phi porosity,
    rhoc_r [J/m3/K] rock grains and confining rock, k_cap [W/m/K] confining rock conductivity,
    T_aq [C] initial aquifer temperature, T_inj [C] injection temperature.
    Loss-free arrival time of the streamline leaving the injector at angle theta:
        tau(theta) = (pi H D^2 (rhoc)_a / (Q (rhoc)_w)) (theta - sin theta cos theta) / (2 sin^3 theta)
    (tau(0) = pi H D^2 (rhoc)_a / (3 Q (rhoc)_w)); along it (Lauwerier)
        T_D = erfc( beta tau / (2 sqrt(t - tau)) ),  beta = 2 sqrt(k_cap rhoc_cap) / ((rhoc)_a H).
    """

    def __init__(self, Q, D, H, phi, rhoc_r, k_cap, T_aq_C, T_inj_C, rhoc_w=RHOC_W, rhoc_cap=None, n_theta=400,
                 k_aq=None, mu=None, r_w=0.1):
        self.Q, self.D, self.H, self.phi = Q, D, H, phi
        self.rhoc_a = phi * rhoc_w + (1 - phi) * rhoc_r
        self.rhoc_w = rhoc_w
        rhoc_cap = rhoc_r if rhoc_cap is None else rhoc_cap
        self.beta = 2.0 * math.sqrt(k_cap * rhoc_cap) / (self.rhoc_a * H)
        self.T_aq = T_aq_C + 273.15
        self.T_inj = T_inj_C
        th = (np.arange(n_theta) + 0.5) / n_theta * math.pi
        self.tau = (math.pi * H * D * D * self.rhoc_a / (Q * rhoc_w)) * (th - np.sin(th) * np.cos(th)) / (2 * np.sin(th) ** 3)
        self.k_aq, self.mu, self.r_w = k_aq, mu, r_w

    @property
    def t_breakthrough(self):
        return float(self.tau.min())

    def T_D(self, t):
        """Dimensionless production temperature drop (0 -> 1)."""
        t = np.asarray(t, float)
        tau = self.tau[None, :]
        tt = t[:, None] if t.ndim else np.array([[t]])
        arg = np.where(tt > tau, self.beta * tau / (2.0 * np.sqrt(np.maximum(tt - tau, 1e-30))), np.inf)
        TD = erfc(arg)
        TD = np.where(tt > tau, TD, 0.0)
        return TD.mean(axis=1)

    def __call__(self, t):
        TD = float(self.T_D(np.array([t]))[0])
        Tin = (self.T_inj(t) if callable(self.T_inj) else self.T_inj) + 273.15
        return self.T_aq - TD * (self.T_aq - Tin)

    def pumping_dp(self):
        """Producer drawdown + injector overpressure from the Thiem equation (doublet, steady), Pa."""
        if self.k_aq is None or self.mu is None:
            return 0.0
        return 2.0 * self.Q * self.mu * math.log(self.D / self.r_w) / (2.0 * math.pi * self.k_aq * self.H)


# ----------------------------------------------------------------------------
# Point source and sink a distance D apart in an unbounded homogeneous medium.  The streamline that leaves the source at an
# angle a from the line to the sink bounds the fraction F = (1 - cos a) / 2 of the flow; its travel time, in units of
# D^3 (rho c)_b / (Q (rho c)_w), is c(F) / (1 - F)^4 with c falling slowly from 0.746 on the axis to 0.718 far out
# (integrated along the streamlines to 1e-6; validation/segs_sandbox.py repeats the integration).
_DIP_F = np.array([0.0, 0.01, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8,
                   0.85, 1.0])
_DIP_C = np.array([0.746012, 0.745435, 0.744579, 0.743176, 0.741803, 0.740462, 0.737877, 0.735426, 0.733115, 0.730952, 0.728941,
                   0.727089, 0.725402, 0.723887, 0.722548, 0.721390, 0.720414, 0.719622, 0.719010, 0.718570, 0.718287, 0.7180])


def dipole_arrival(F):
    """Arrival time of the streamline bounding the fraction F of the flow of a point source-sink pair, in units of
    D^3 (rho c)_b / (Q (rho c)_w).  First arrival (F = 0): 0.746."""
    F = np.asarray(F, float)
    return np.interp(F, _DIP_F, _DIP_C) / np.maximum(1.0 - F, 1e-12) ** 4


class SingleWellDipole:
    """Single-well circulation: water leaves the well through one screen and returns through another a distance D along
    the well, through a porous or densely fractured block that is large compared with D.

    Q [m3/s] circulated, D [m] between the screen centres, rhoc_b [J/m3/K] volumetric heat capacity of the swept medium
    (rock and the water in it), T_rock [C], T_inj [C].  The flow is that of a point source and sink in an unbounded
    homogeneous medium; heat moves with the water along each streamline:

      S is None   local thermal equilibrium (a porous medium, or fractures closer than the conduction length): a sharp
                  front on every streamline, so T_D(t) is the fraction of the flow whose streamlines have been flushed.
                  First arrival at t_bt = 0.746 D^3 (rho c)_b / (Q (rho c)_w); afterwards 1 - T_D = (c / x)^(1/4),
                  x = t Q (rho c)_w / (D^3 (rho c)_b), c = 0.72-0.75: a long tail, half the rise is still there at x = 11.6.
      S [m]       flowing fractures S apart with conduction from the rock between them (k_r [W/m/K]): the parallel-
                  fracture response of Gringarten et al. (1975) along every streamline, with the streamline's own
                  front-crossing time.

    Not in the element: an outer boundary (a bounded block cools faster once the front has reached the boundary), the
    length of the screens, and the effect of the temperature on the flow (cold water is heavier and more viscous; with the
    injection screen above the production screen the cold water sinks towards the producer).  aharrsim.porous solves
    those cases numerically.
    """

    def __init__(self, Q, D, rhoc_b, T_rock_C, T_inj_C, rhoc_w=RHOC_W, S=None, k_r=None, n_tubes=120, impedance_Pa_per_m3s=None):
        self.Q, self.D, self.rhoc_b, self.rhoc_w, self.S, self.k_r = Q, D, rhoc_b, rhoc_w, S, k_r
        self.T_rock = T_rock_C + 273.15
        self.T_inj = T_inj_C
        self.t_ref = D ** 3 * rhoc_b / (Q * rhoc_w)
        self.t_bt = float(_DIP_C[0]) * self.t_ref
        self.impedance = impedance_Pa_per_m3s
        Fg = np.concatenate([np.linspace(0.0, 0.9, 901), 1.0 - np.geomspace(0.1, 1e-6, 400)[1:]])
        self._xg, self._Fg = dipole_arrival(Fg), Fg
        if S is not None:
            if k_r is None:
                raise ValueError("k_r is needed with a fracture spacing S")
            self.t_d = (0.5 * S) ** 2 * rhoc_b / k_r
            Fk = (np.arange(n_tubes) + 0.5) / n_tubes
            self._Rk = dipole_arrival(Fk) * self.t_ref / self.t_d

    def T_D(self, t):
        """Dimensionless temperature drop of the produced water, 0 at first -> 1 at long time."""
        t = np.atleast_1d(np.asarray(t, float))
        if self.S is None:
            return np.interp(t / self.t_ref, self._xg, self._Fg, left=0.0, right=1.0)
        out = np.zeros(len(t))
        R = self._Rk

        def f(s):
            r = np.sqrt(s)
            return np.exp(-R[:, None] * r * np.tanh(r)) / s

        for i, ti in enumerate(t):
            if ti <= 0.0:
                continue
            with np.errstate(all="ignore"):
                v = dehoog(f, np.full(len(R), ti / self.t_d))
            out[i] = float(np.mean(np.clip(np.where(np.isfinite(v), v, 0.0), 0.0, 1.0)))
        return out

    def groups(self, t):
        """x = t / [D^3 (rho c)_b / (Q (rho c)_w)] (first arrival at 0.746) and, with a fracture spacing, y = (S / 2) / sqrt(alpha t)."""
        t = np.asarray(t, float)
        return t / self.t_ref, (None if self.S is None else 0.5 * self.S / np.sqrt(self.k_r / self.rhoc_b * t))

    def T_out(self, T_in_K, t):
        TD = float(self.T_D(max(t, 1.0))[0])
        return self.T_rock - TD * (self.T_rock - T_in_K)

    def __call__(self, t):
        Tin = self.T_inj(t) if callable(self.T_inj) else self.T_inj
        return self.T_out(Tin + 273.15, t)

    def pumping_dp(self):
        return self.Q * self.impedance if self.impedance else 0.0


# ----------------------------------------------------------------------------
class CyclicFracture:
    """Single-well huff-and-puff through one planar fracture (Sage HeatRoot type).

    Plug flow along the fracture (length L_f, breadth B, aperture w) from the well; every fracture
    cell exchanges heat with the rock faces through a 1-D conduction column (both faces by symmetry).
    Schedule per cycle: inject V_cycle at T_inj over t_inj, soak t_soak, produce V_cycle over t_prod.
    Run with simulate(n_cycles); results per cycle: produced temperature history, mean, energy.
    For benchmarking, simulate_continuous(q, t_end) injects continuously (Lauwerier solution).
    """

    def __init__(self, L_f, B, w, k_r, rhoc_r, T_rock_C, T_inj_C, rhoc_w=RHOC_W, n_x=60, n_z=60, z_max=None):
        self.L, self.B, self.w = L_f, B, w
        self.k_r, self.rhoc_r, self.rhoc_w = k_r, rhoc_r, rhoc_w
        self.T_rock, self.T_inj = T_rock_C + 273.15, T_inj_C + 273.15
        self.nx, self.nz = n_x, n_z
        self.dx = L_f / n_x
        self.cell_vol = self.dx * B * w
        alpha = k_r / rhoc_r
        z_max = z_max or 4.0 * math.sqrt(alpha * 10 * 365.25 * 86400)   # ten-year diffusion depth by default
        self.z = np.geomspace(0.002, z_max, n_z)                          # rock column (half-space on each face)
        self.z = np.concatenate([[0.0], self.z])
        self.reset()

    def reset(self):
        self.Tf = np.full(self.nx, self.T_rock)                 # fluid in the fracture cells (initially at rock T)
        self.Tr = np.full((self.nx, self.z.size), self.T_rock)  # rock columns, node 0 = face
        self.t = 0.0

    def _conduct(self, dt):
        """Implicit 1-D conduction in every column with the fluid cell coupled to the face (symmetry: half the
        fluid capacity per face)."""
        z = self.z
        n = z.size
        dz = np.diff(z)
        # node volumes (per unit face area)
        vol = np.concatenate([[0.5 * dz[0]], 0.5 * (dz[:-1] + dz[1:]), [0.5 * dz[-1]]])
        Cr = self.rhoc_r * vol / dt
        G = self.k_r / dz
        Cf = 0.5 * self.rhoc_w * self.w / dt                     # fluid half-aperture per unit face area
        Gf = 2.0 * self.k_r / dz[0]                              # fluid node to first rock node (half cell)
        # unknowns: [Tf, Tr_0..Tr_n-1] per column; tridiagonal
        N = n + 1
        main = np.zeros(N); lower = np.zeros(N - 1); upper = np.zeros(N - 1)
        main[0] = Cf + Gf
        upper[0] = -Gf
        lower[0] = -Gf
        main[1] = Cr[0] + Gf + G[0]
        upper[1] = -G[0]
        for i in range(1, n - 1):
            lower[i] = -G[i - 1]
            main[i + 1] = Cr[i] + G[i - 1] + G[i]
            upper[i + 1] = -G[i]
        lower[n - 1] = -G[n - 2]
        main[n] = Cr[n - 1] + G[n - 2]
        from scipy.linalg import solve_banded
        ab = np.zeros((3, N))
        ab[0, 1:] = upper
        ab[1, :] = main
        ab[2, :-1] = lower
        rhs = np.empty((N, self.nx))
        rhs[0, :] = Cf * self.Tf
        rhs[1:, :] = (Cr[:, None] * self.Tr.T)
        # far node: Dirichlet at the undisturbed rock temperature through the last conductance
        rhs[n, :] += 0.0
        sol = solve_banded((1, 1), ab, rhs)
        self.Tf = sol[0, :].copy()
        self.Tr = sol[1:, :].T.copy()
        self.t += dt

    def _shift_in(self, T_new):
        """Injection of one cell volume: fluid moves one cell away from the well; the last cell leaves the fracture."""
        self.Tf = np.concatenate([[T_new], self.Tf[:-1]])

    def _shift_out(self):
        """Production of one cell volume: fluid moves one cell towards the well; returns the produced temperature.
        The far end refills with water at the local rock-face temperature (fracture tip in contact with rock)."""
        T_prod = self.Tf[0]
        self.Tf = np.concatenate([self.Tf[1:], [self.Tr[-1, 0]]])
        return T_prod

    def simulate_continuous(self, q, t_end, sample_t):
        """Continuous injection at q m3/s; returns the temperature leaving the far end at sample times."""
        self.reset()
        dt = self.cell_vol / q
        out = []
        T_exit = self.T_rock
        for ts in np.sort(sample_t):
            while self.t + dt <= ts + 1e-9:
                self._conduct(dt)
                T_exit = self.Tf[-1]
                self._shift_in(self.T_inj)
            out.append(T_exit)
        return np.array(out)

    def simulate(self, V_cycle, t_inj, t_soak, t_prod, n_cycles, T_inj_C=None, n_sub=4):
        """Run n_cycles of inject / soak / produce.  Returns a list of dicts per cycle."""
        self.reset()
        if T_inj_C is not None:
            self.T_inj = T_inj_C + 273.15
        n_cells = max(1, int(round(V_cycle / self.cell_vol)))
        n_cells = min(n_cells, self.nx)
        q_inj, q_prod = n_cells * self.cell_vol / t_inj, n_cells * self.cell_vol / t_prod
        cycles = []
        for ic in range(n_cycles):
            # inject
            dt = self.cell_vol / q_inj
            for _ in range(n_cells):
                self._conduct(dt)
                self._shift_in(self.T_inj)
            # soak
            for _ in range(n_sub):
                self._conduct(t_soak / n_sub)
            # produce
            dt = self.cell_vol / q_prod
            Tp, tp = [], []
            for k in range(n_cells):
                self._conduct(dt)
                Tp.append(self._shift_out())
                tp.append((k + 1) * dt)
            Tp = np.array(Tp)
            E = self.rhoc_w * n_cells * self.cell_vol * np.mean(Tp - self.T_inj)
            cycles.append(dict(cycle=ic + 1, T_prod_C=Tp - 273.15, t_prod=np.array(tp), T_mean_C=float(Tp.mean() - 273.15),
                               T_first_C=float(Tp[0] - 273.15), T_last_C=float(Tp[-1] - 273.15), E_J=float(E),
                               V=n_cells * self.cell_vol, P_mean_W=float(E / t_prod)))
        return cycles


# ----------------------------------------------------------------------------
class HeatStore:
    """Lumped finite heat store: volume V [m3] of rock+fluid at T_store, swept by a doublet circulating m [kg/s]
    (well mixed), losing heat to surrounding rock at T_0 across area A_b by 1-D conduction (Duhamel).
    Produced temperature = store temperature.  step(dt, m, T_inj) advances the balance."""

    def __init__(self, V, A_b, rhoc_store, k_r, rhoc_r, T_store_C, T0_C):
        self.V, self.A_b, self.rhoc = V, A_b, rhoc_store
        self.k_r, self.alpha = k_r, k_r / rhoc_r
        self.T = T_store_C + 273.15
        self.T0 = T0_C + 273.15
        self.hist_t, self.hist_T = [0.0], [self.T]
        self.t = 0.0

    def _loss(self, dt):
        """Conductive loss over the next dt from the boundary temperature history (piecewise-constant Duhamel)."""
        t_end = self.t + dt
        q = 0.0
        Ts = np.array(self.hist_T + [self.T])
        ts = np.array(self.hist_t + [self.t])
        dT = np.diff(np.concatenate([[self.T0], Ts]))
        # flux from a step dT applied at time ts_i: k dT / sqrt(pi alpha (t - ts_i))
        tau = np.maximum(t_end - ts, 0.5 * dt)
        q = self.k_r * np.sum(dT / np.sqrt(math.pi * self.alpha * tau))
        return max(q, 0.0) * self.A_b

    def step(self, dt, m, T_inj_K, cp_w=4180.0):
        loss = self._loss(dt)
        extract = m * cp_w * (self.T - T_inj_K)
        self.T -= (extract + loss) * dt / (self.rhoc * self.V)
        self.T = max(self.T, self.T0)
        self.t += dt
        self.hist_t.append(self.t)
        self.hist_T.append(self.T)
        return self.T, extract, loss

    def run(self, years, m, T_inj_C, dt_days=10.0):
        out = dict(t=[], T_C=[], Q_MW=[], loss_MW=[])
        n = int(years * 365.25 / dt_days)
        for _ in range(n):
            T, ex, lo = self.step(dt_days * 86400.0, m, T_inj_C + 273.15)
            out["t"].append(self.t); out["T_C"].append(T - 273.15); out["Q_MW"].append(ex / 1e6); out["loss_MW"].append(lo / 1e6)
        return {k: np.array(v) for k, v in out.items()}
