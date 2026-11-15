"""Darcy flow and heat transport in a porous or densely fractured block round a well (axisymmetric r-z).

For single-well systems that inject through one screen and produce through another on the same well: the laboratory
sandbox of a single-well enhanced geothermal system, or a field well with an injection and a production zone.

    u = -(K / mu(T)) (grad p + rho(T) g e_z)            Darcy flux, Boussinesq buoyancy, viscosity of the local water
    div u = q                                           q: screens (uniform flux along each screen)
    (rho c)_b dT/dt + div((rho c)_w u T) = div(lam grad T) + (rho c)_w q T_q        local thermal equilibrium

Finite volumes on a radial grid that is fine at the well; second-order upwind (MUSCL, van Leer slopes, clipped) with a
two-stage Runge-Kutta step, so that a thermal front stays sharp; the pressure field is solved again as the water cools,
because both the weight and the viscosity of the water change the flow. Outer boundaries closed and adiabatic.

Checks (validation/segs_sandbox.py): energy balance to round-off; a piston front; the point source-sink pair in an
unbounded medium (reservoir.SingleWellDipole) reproduced in a large block.
"""
from __future__ import annotations

import math

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve

G = 9.80665


def water_table(p_Pa=2.0e5, T_min_C=1.0, T_max_C=None, n=200):
    """rho(T) [kg/m3], mu(T) [Pa s] and rho c_p (T) [J/m3/K] of liquid water at pressure p, tabulated for np.interp."""
    import CoolProp.CoolProp as CP
    if T_max_C is None:
        T_max_C = min(CP.PropsSI("T", "P", p_Pa, "Q", 0, "Water") - 273.15 - 1.0, 370.0) if p_Pa < 22.0e6 else 370.0
    T = np.linspace(T_min_C, T_max_C, n)
    rho = np.array([CP.PropsSI("D", "T", t + 273.15, "P", p_Pa, "Water") for t in T])
    mu = np.array([CP.PropsSI("V", "T", t + 273.15, "P", p_Pa, "Water") for t in T])
    cp = np.array([CP.PropsSI("C", "T", t + 273.15, "P", p_Pa, "Water") for t in T])
    return dict(T=T, rho=rho, mu=mu, rhoc=rho * cp)


class AxisymmetricBlock:
    """A cylinder of porous rock, radius R and height H, round a well of radius r_w on its axis.

    K permeability [m2] (K_z vertical, if different), phi porosity, rhoc_s volumetric heat capacity of the solid [J/m3/K], lam conductivity of the
    saturated medium [W/m/K], disp thermal dispersivity [m], T0 initial temperature [C], p_Pa pressure for the water
    properties. buoyancy / viscosity switch the two couplings of temperature back on the flow.
    """

    def __init__(self, R, H, r_w, K, phi, rhoc_s, lam, T0, p_Pa=2.0e5, nr=28, nz=120, r1=None, disp=0.0, buoyancy=True,
                 viscosity=True, T_ref=None, rhoc_w=None, water=None, K_z=None):
        self.R, self.H, self.r_w, self.K, self.phi, self.lam, self.disp, self.T0 = R, H, r_w, K, phi, lam, disp, T0
        self.K_z = K if K_z is None else K_z                                 # vertical permeability, if it differs
        self.buoyancy, self.viscosity = buoyancy, viscosity
        # radial cells widening geometrically from the well: first cell dr1 wide (default the larger of 1.5 r_w and R / 40)
        dr1 = (r1 - r_w) if r1 is not None else max(1.5 * r_w, R / 40.0)
        if nr * dr1 >= R - r_w:
            widths = np.full(nr, (R - r_w) / nr)
        else:
            from scipy.optimize import brentq
            g = brentq(lambda g: dr1 * (g ** nr - 1.0) / (g - 1.0) - (R - r_w), 1.0 + 1e-9, 10.0)
            widths = dr1 * g ** np.arange(nr)
        self.rf = r_w + np.concatenate([[0.0], np.cumsum(widths)])           # nr + 1 faces, nr cells
        self.rf[-1] = R
        self.zf = np.linspace(0.0, H, nz + 1)
        self.nr, self.nz = nr, nz
        self.rc = 0.5 * (self.rf[:-1] + self.rf[1:])
        self.zc = 0.5 * (self.zf[:-1] + self.zf[1:])
        self.dz = np.diff(self.zf)
        ring = math.pi * (self.rf[1:] ** 2 - self.rf[:-1] ** 2)
        self.V = ring[:, None] * self.dz[None, :]                            # cell volumes (nr, nz)
        self.Ar = 2.0 * math.pi * self.rf[1:-1, None] * self.dz[None, :]      # radial faces (nr - 1, nz)
        self.Az = np.repeat(ring[:, None], nz - 1, axis=1)                   # vertical faces (nr, nz - 1)
        self.dzc = np.diff(self.zc)
        self.drc = np.diff(self.rc)
        self.TR = 2.0 * math.pi * K * self.dz[None, :] / np.log(self.rc[1:] / self.rc[:-1])[:, None]
        self.TZ = self.K_z * self.Az / self.dzc[None, :]
        self.w = water or water_table(p_Pa)
        self.T_ref = T0 if T_ref is None else T_ref
        self.rhoc_w = float(np.interp(self.T_ref, self.w["T"], self.w["rhoc"])) if rhoc_w is None else rhoc_w
        self.rhoc_b = (1.0 - phi) * rhoc_s + phi * self.rhoc_w
        self.src = np.zeros((nr, nz))            # volumetric source [m3/s] per cell: + injection, - production
        self.T_src = np.zeros((nr, nz))          # temperature of the injected water
        self.screens = []
        self.T = np.full((nr, nz), float(T0))
        self.t = 0.0
        # sparse pattern for the pressure solve
        idx = np.arange(nr * nz).reshape(nz, nr).T                           # idx[i, j] = i + nr j
        self._ra, self._rb = idx[:-1, :].ravel(), idx[1:, :].ravel()
        self._za, self._zb = idx[:, :-1].ravel(), idx[:, 1:].ravel()
        self._kpin = int(idx[-1, -1])

    # ------------------------------------------------------------------ set-up
    def add_screen(self, z0, z1, Q, T_in=None):
        """A screen between heights z0 and z1 [m above the base] carrying Q [m3/s]: Q > 0 injects water at T_in, Q < 0 produces."""
        j = np.where((self.zc > z0) & (self.zc < z1))[0]
        if len(j) == 0:
            raise ValueError("screen falls between cells")
        self.src[0, j] += Q * self.dz[j] / self.dz[j].sum()
        if Q > 0:
            self.T_src[0, j] = T_in
        self.screens.append(dict(z0=z0, z1=z1, Q=Q, T_in=T_in, j=j))
        return self

    # ------------------------------------------------------------------ flow
    def _props(self, T):
        w = self.w
        mu = np.interp(T, w["T"], w["mu"]) if self.viscosity else np.full_like(T, float(np.interp(self.T_ref, w["T"], w["mu"])))
        rho = np.interp(T, w["T"], w["rho"]) if self.buoyancy else np.full_like(T, float(np.interp(self.T_ref, w["T"], w["rho"])))
        return mu, rho

    def solve_flow(self):
        mu, rho = self._props(self.T)
        cr = (self.TR / (0.5 * (mu[:-1, :] + mu[1:, :]))).ravel()
        cz = (self.TZ / (0.5 * (mu[:, :-1] + mu[:, 1:]))).ravel()
        gam = (0.5 * (rho[:, :-1] + rho[:, 1:]) * G * self.dzc[None, :]).ravel()
        n = self.nr * self.nz
        rows = np.concatenate([self._ra, self._rb, self._ra, self._rb, self._za, self._zb, self._za, self._zb])
        cols = np.concatenate([self._ra, self._rb, self._rb, self._ra, self._za, self._zb, self._zb, self._za])
        vals = np.concatenate([cr, cr, -cr, -cr, cz, cz, -cz, -cz])
        A = sp.coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsc()
        b = self.src.T.ravel().copy()                                        # index i + nr j
        np.add.at(b, self._za, cz * gam)
        np.add.at(b, self._zb, -cz * gam)
        A = A + sp.coo_matrix(([cr.max()], ([self._kpin], [self._kpin])), shape=(n, n)).tocsc()
        p = spsolve(A, b)
        P = p.reshape(self.nz, self.nr).T
        self.p = P
        self.FR = (cr.reshape(self.nr - 1, self.nz)) * (P[:-1, :] - P[1:, :])                       # + outward
        self.FZ = (cz.reshape(self.nr, self.nz - 1)) * (P[:, :-1] - P[:, 1:] - gam.reshape(self.nr, self.nz - 1))   # + upward
        out = np.zeros((self.nr, self.nz))
        out[:-1, :] += np.maximum(self.FR, 0.0); out[1:, :] += np.maximum(-self.FR, 0.0)
        out[:, :-1] += np.maximum(self.FZ, 0.0); out[:, 1:] += np.maximum(-self.FZ, 0.0)
        out += np.maximum(-self.src, 0.0)
        # conduction and dispersion on the faces
        lr = self.lam + self.disp * self.rhoc_w * np.abs(self.FR) / self.Ar
        lz = self.lam + self.disp * self.rhoc_w * np.abs(self.FZ) / self.Az
        self.CR = lr * self.Ar / self.drc[:, None]
        self.CZ = lz * self.Az / self.dzc[None, :]
        cond = np.zeros((self.nr, self.nz))
        cond[:-1, :] += self.CR; cond[1:, :] += self.CR; cond[:, :-1] += self.CZ; cond[:, 1:] += self.CZ
        self.dt_max = float(np.min(self.rhoc_b * self.V / (self.rhoc_w * out + 2.0 * cond + 1e-300)))
        return self

    # ------------------------------------------------------------------ transport
    @staticmethod
    def _slope(a, b):
        """van Leer (harmonic) limited slope from the two one-sided gradients."""
        s = np.zeros_like(a)
        m = a * b > 0.0
        s[m] = 2.0 * a[m] * b[m] / (a[m] + b[m])
        return s

    def _rate(self, T):
        nr, nz = self.nr, self.nz
        gr = (T[1:, :] - T[:-1, :]) / self.drc[:, None]
        gz = (T[:, 1:] - T[:, :-1]) / self.dzc[None, :]
        sr = np.zeros_like(T); sz = np.zeros_like(T)
        sr[1:-1, :] = self._slope(gr[:-1, :], gr[1:, :])
        sz[:, 1:-1] = self._slope(gz[:, :-1], gz[:, 1:])
        # face temperatures from the upwind cell, clipped between the two neighbours
        lo = np.minimum(T[:-1, :], T[1:, :]); hi = np.maximum(T[:-1, :], T[1:, :])
        Tl = T[:-1, :] + sr[:-1, :] * (self.rf[1:-1] - self.rc[:-1])[:, None]
        Tr = T[1:, :] + sr[1:, :] * (self.rf[1:-1] - self.rc[1:])[:, None]
        TfR = np.clip(np.where(self.FR >= 0.0, Tl, Tr), lo, hi)
        lo = np.minimum(T[:, :-1], T[:, 1:]); hi = np.maximum(T[:, :-1], T[:, 1:])
        Tb = T[:, :-1] + sz[:, :-1] * (self.zf[1:-1] - self.zc[:-1])[None, :]
        Tt = T[:, 1:] + sz[:, 1:] * (self.zf[1:-1] - self.zc[1:])[None, :]
        TfZ = np.clip(np.where(self.FZ >= 0.0, Tb, Tt), lo, hi)
        HR = self.rhoc_w * self.FR * TfR + self.CR * (T[:-1, :] - T[1:, :])
        HZ = self.rhoc_w * self.FZ * TfZ + self.CZ * (T[:, :-1] - T[:, 1:])
        dE = np.zeros((nr, nz))
        dE[:-1, :] -= HR; dE[1:, :] += HR
        dE[:, :-1] -= HZ; dE[:, 1:] += HZ
        dE += self.rhoc_w * (np.maximum(self.src, 0.0) * self.T_src + np.minimum(self.src, 0.0) * T)
        return dE / (self.rhoc_b * self.V)

    def produced(self, T=None):
        """Flow-weighted temperature of the water entering the production screens [C]."""
        T = self.T if T is None else T
        q = np.maximum(-self.src, 0.0)
        return float((q * T).sum() / q.sum())

    def line(self, r=None, z=None):
        """Temperatures along a vertical line at radius r (returns zc, T) or along a radius at height z (returns rc, T)."""
        if r is not None:
            i = int(np.clip(np.searchsorted(self.rc, r) - 1, 0, self.nr - 2))
            f = (r - self.rc[i]) / (self.rc[i + 1] - self.rc[i])
            return self.zc, (1 - f) * self.T[i, :] + f * self.T[i + 1, :]
        j = int(np.clip(np.searchsorted(self.zc, z) - 1, 0, self.nz - 2))
        f = (z - self.zc[j]) / (self.zc[j + 1] - self.zc[j])
        return self.rc, (1 - f) * self.T[:, j] + f * self.T[:, j + 1]

    def run(self, t_end, dt_out=30.0, flow_update=5.0, cfl=0.45, snapshots=(), lines=None, progress=None):
        """March to t_end [s]. Returns the history of the produced temperature and of the energy balance, the temperature
        fields at the snapshot times and the temperatures along the monitoring lines (dict name -> dict(r=..) or dict(z=..))."""
        steady = not (self.buoyancy or self.viscosity)
        self.solve_flow()
        Qp = float(np.maximum(-self.src, 0.0).sum())
        T_in = float((np.maximum(self.src, 0.0) * self.T_src).sum() / max(np.maximum(self.src, 0.0).sum(), 1e-300))
        hist = dict(t=[self.t], T_prod=[self.produced()], E_out=[0.0], E_block=[0.0])
        snaps, mon = {}, {k: dict(t=[], T=[], x=None, **v) for k, v in (lines or {}).items()}
        snapshots = sorted(snapshots)
        E_out, t_flow, t_next = 0.0, self.t, self.t + dt_out
        E0 = float((self.rhoc_b * self.V * self.T).sum())

        def record():
            hist["t"].append(self.t); hist["T_prod"].append(self.produced()); hist["E_out"].append(E_out)
            hist["E_block"].append(float((self.rhoc_b * self.V * self.T).sum()) - E0)
            for k, m in mon.items():
                x, Tl = self.line(r=m.get("r"), z=m.get("z"))
                m["x"] = x.tolist(); m["t"].append(self.t); m["T"].append(Tl.tolist())

        for k, m in mon.items():
            x, Tl = self.line(r=m.get("r"), z=m.get("z")); m["x"] = x.tolist(); m["t"].append(self.t); m["T"].append(Tl.tolist())
        while self.t < t_end - 1e-9:
            if not steady and self.t - t_flow >= flow_update:
                self.solve_flow(); t_flow = self.t
            stops = [t_end, t_next] + [s for s in snapshots if s > self.t + 1e-9]
            dt = min(cfl * self.dt_max, min(stops) - self.t)
            T0 = self.T
            T1 = T0 + dt * self._rate(T0)
            self.T = 0.5 * T0 + 0.5 * (T1 + dt * self._rate(T1))
            E_out += dt * self.rhoc_w * Qp * (0.5 * (self.produced(T0) + self.produced(T1)) - T_in)
            self.t += dt
            for s in snapshots:
                if abs(self.t - s) < 1e-6 and s not in snaps:
                    snaps[s] = self.T.copy()
            if self.t >= t_next - 1e-9:
                record(); t_next += dt_out
                if progress:
                    progress(self)
        if hist["t"][-1] < self.t - 1e-9:
            record()
        h = {k: np.array(v) for k, v in hist.items()}
        return dict(t=h["t"], T_prod=h["T_prod"], E_out=h["E_out"], E_block=h["E_block"], T_in=T_in, Q=Qp, snapshots=snaps, lines=mon)
