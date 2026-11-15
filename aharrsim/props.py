"""Water / steam properties.

Primary model: IAPWS-95 (thermodynamic) with the IAPWS 2008 viscosity and
IAPWS 2011 thermal-conductivity formulations, all through CoolProp.  For speed
the single-phase states are evaluated with CoolProp's bicubic (p, h) tables
built from IAPWS-95 ("BICUBIC&HEOS"); states the tables cannot handle fall back
to the full Helmholtz equation of state ("HEOS").

Two-phase (flashing) states are handled explicitly with a saturation table
computed once from IAPWS-95, so that quality, phase densities and phase
transport properties are always consistent.

All units SI: p [Pa], h [J/kg], T [K].
"""
from __future__ import annotations

import numpy as np
import CoolProp
from CoolProp import AbstractState

T_CRIT = 647.096        # K
P_CRIT = 22.064e6       # Pa
P_TRIPLE = 611.655      # Pa
T_TRIPLE = 273.16       # K


def surface_tension_iapws(T):
    """IAPWS (2014) surface tension of ordinary water [N/m]."""
    tau = np.clip(1.0 - np.asarray(T, float) / T_CRIT, 0.0, 1.0)
    return 0.2358 * tau ** 1.256 * (1.0 - 0.625 * tau)


class NodeProps:
    """Container for vectorised fluid properties at a set of nodes."""

    __slots__ = ("p", "h", "T", "rho", "cp", "mu", "k", "x", "two_phase",
                 "dTdh", "rho_l", "rho_g", "mu_l", "mu_g", "k_l", "cp_l",
                 "h_l", "h_g", "sigma", "Tsat")

    def __init__(self, n):
        for name in self.__slots__:
            setattr(self, name, np.full(n, np.nan))
        self.two_phase = np.zeros(n, bool)


class WaterIAPWS:
    """IAPWS-95 water/steam (CoolProp)."""

    name = "IAPWS-95 water/steam (CoolProp, bicubic p-h tables)"
    is_constant = False

    def __init__(self, backend: str = "BICUBIC&HEOS"):
        self.backend = backend
        self.AS = AbstractState(backend, "Water")
        self.ASx = AbstractState("HEOS", "Water")      # exact fallback
        self._build_saturation_table()

    # ------------------------------------------------------------------
    def _build_saturation_table(self, n=3000):
        AS = self.ASx
        # log-spaced up to 0.9 p_c, then clustered towards the critical point
        p_lo = np.exp(np.linspace(np.log(P_TRIPLE * 1.0001), np.log(P_CRIT * 0.9), n))
        p_hi = P_CRIT * (1.0 - np.geomspace(0.1, 2e-5, 600))
        p = np.unique(np.concatenate([p_lo, p_hi]))
        n = p.size
        cols = {k: np.empty(n) for k in
                ("T", "h_l", "h_g", "rho_l", "rho_g", "mu_l", "mu_g",
                 "k_l", "k_g", "cp_l", "cp_g", "s_l", "s_g")}
        for i, pi in enumerate(p):
            AS.update(CoolProp.PQ_INPUTS, pi, 0.0)
            cols["T"][i] = AS.T()
            cols["h_l"][i] = AS.hmass()
            cols["rho_l"][i] = AS.rhomass()
            cols["mu_l"][i] = AS.viscosity()
            cols["k_l"][i] = AS.conductivity()
            cols["cp_l"][i] = AS.cpmass()
            cols["s_l"][i] = AS.smass()
            AS.update(CoolProp.PQ_INPUTS, pi, 1.0)
            cols["h_g"][i] = AS.hmass()
            cols["rho_g"][i] = AS.rhomass()
            cols["mu_g"][i] = AS.viscosity()
            cols["k_g"][i] = AS.conductivity()
            cols["cp_g"][i] = AS.cpmass()
            cols["s_g"][i] = AS.smass()
        self._sat_lnp = np.log(p)
        self._sat = cols

    def sat(self, p, name):
        """Saturation property `name` at pressure p (p < p_crit)."""
        lnp = np.log(np.clip(np.asarray(p, float), P_TRIPLE * 1.0001, P_CRIT * (1 - 2e-5)))
        return np.interp(lnp, self._sat_lnp, self._sat[name])

    def Tsat(self, p):
        return self.sat(p, "T")

    # ------------------------------------------------------------------
    def _single(self, p, h):
        """Single-phase state at (p, h): returns T, rho, cp, mu, k."""
        AS = self.AS
        try:
            AS.update(CoolProp.HmassP_INPUTS, h, p)
            T, rho, cp = AS.T(), AS.rhomass(), AS.cpmass()
            # transport properties evaluated exactly (IAPWS 2008 / 2011) at the
            # tabulated (rho, T): the bicubic tables are not accurate enough for
            # viscosity near 0 C (up to 10 % error)
            X = self.ASx
            X.update(CoolProp.DmassT_INPUTS, rho, T)
            out = (T, rho, cp, X.viscosity(), X.conductivity())
            if not np.all(np.isfinite(out)) or cp <= 0.0:
                raise ValueError
            return out
        except Exception:
            AS = self.ASx
            AS.update(CoolProp.HmassP_INPUTS, h, p)
            cp = AS.cpmass()
            if not np.isfinite(cp) or cp <= 0.0:      # numerical trouble right at the critical point
                cp = 1.0e6
            return (AS.T(), AS.rhomass(), cp, AS.viscosity(), AS.conductivity())

    def ph(self, p, h) -> NodeProps:
        """Vectorised properties at nodes given pressure p [Pa] and enthalpy h [J/kg]."""
        p = np.atleast_1d(np.asarray(p, float))
        h = np.atleast_1d(np.asarray(h, float))
        n = p.size
        out = NodeProps(n)
        out.p[:] = p
        out.h[:] = h
        sub = p < P_CRIT * (1 - 2e-5)
        h_l = np.full(n, np.nan)
        h_g = np.full(n, np.nan)
        if np.any(sub):
            h_l[sub] = self.sat(p[sub], "h_l")
            h_g[sub] = self.sat(p[sub], "h_g")
        tp = sub & (h > h_l) & (h < h_g)
        out.two_phase[:] = tp
        # --- two-phase nodes: equilibrium mixture from the saturation table
        if np.any(tp):
            pp = p[tp]
            hl, hg = h_l[tp], h_g[tp]
            x = (h[tp] - hl) / (hg - hl)
            rl, rg = self.sat(pp, "rho_l"), self.sat(pp, "rho_g")
            out.x[tp] = x
            out.T[tp] = self.sat(pp, "T")
            out.Tsat[tp] = out.T[tp]
            out.rho[tp] = 1.0 / (x / rg + (1.0 - x) / rl)
            out.rho_l[tp], out.rho_g[tp] = rl, rg
            out.mu_l[tp], out.mu_g[tp] = self.sat(pp, "mu_l"), self.sat(pp, "mu_g")
            out.k_l[tp], out.cp_l[tp] = self.sat(pp, "k_l"), self.sat(pp, "cp_l")
            out.h_l[tp], out.h_g[tp] = hl, hg
            out.sigma[tp] = surface_tension_iapws(out.T[tp])
            # mixture transport (McAdams viscosity) used only for diagnostics
            out.mu[tp] = 1.0 / (x / out.mu_g[tp] + (1.0 - x) / out.mu_l[tp])
            out.k[tp] = out.k_l[tp]
            out.cp[tp] = out.cp_l[tp]
            out.dTdh[tp] = 0.0
        # --- single-phase nodes
        for i in np.nonzero(~tp)[0]:
            T, rho, cp, mu, k = self._single(p[i], h[i])
            out.T[i], out.rho[i], out.cp[i], out.mu[i], out.k[i] = T, rho, cp, mu, k
            out.dTdh[i] = 1.0 / cp
            if sub[i]:
                out.x[i] = 0.0 if h[i] <= h_l[i] else 1.0
            else:
                out.x[i] = 0.0 if h[i] < 2.1e6 else 1.0   # supercritical: label only
        return out

    def point(self, p, h):
        """Fast single-node state for pressure marching.

        Returns (two_phase, rho, mu, x, rho_l, rho_g, mu_l, mu_g, sigma, T).
        """
        if p < P_CRIT * (1 - 2e-5):
            lnp = np.log(max(p, P_TRIPLE * 1.0001))
            hl = np.interp(lnp, self._sat_lnp, self._sat["h_l"])
            hg = np.interp(lnp, self._sat_lnp, self._sat["h_g"])
            if hl < h < hg:
                x = (h - hl) / (hg - hl)
                S = self._sat
                rl = np.interp(lnp, self._sat_lnp, S["rho_l"])
                rg = np.interp(lnp, self._sat_lnp, S["rho_g"])
                ml = np.interp(lnp, self._sat_lnp, S["mu_l"])
                mg = np.interp(lnp, self._sat_lnp, S["mu_g"])
                T = np.interp(lnp, self._sat_lnp, S["T"])
                rho = 1.0 / (x / rg + (1.0 - x) / rl)
                mu = 1.0 / (x / mg + (1.0 - x) / ml)
                return True, rho, mu, x, rl, rg, ml, mg, float(surface_tension_iapws(T)), T
        T, rho, cp, mu, k = self._single(p, h)
        return False, rho, mu, 0.0, rho, np.nan, mu, np.nan, np.nan, T

    def h_pT(self, p, T):
        p = np.atleast_1d(np.asarray(p, float))
        T = np.atleast_1d(np.asarray(T, float))
        p, T = np.broadcast_arrays(p, T)
        h = np.empty(p.shape)
        for i, (pi, Ti) in enumerate(zip(p.ravel(), T.ravel())):
            try:
                self.AS.update(CoolProp.PT_INPUTS, pi, Ti)
                h.flat[i] = self.AS.hmass()
            except Exception:
                self.ASx.update(CoolProp.PT_INPUTS, pi, Ti)
                h.flat[i] = self.ASx.hmass()
        return h

    def s_ph(self, p, h):
        """Specific entropy [J/kg/K] at (p, h) (exact IAPWS-95)."""
        p = np.atleast_1d(np.asarray(p, float))
        h = np.atleast_1d(np.asarray(h, float))
        s = np.empty(p.size)
        for i in range(p.size):
            self.ASx.update(CoolProp.HmassP_INPUTS, h[i], p[i])
            s[i] = self.ASx.smass()
        return s

    def exact_ph(self, p, h):
        """(T, rho) straight from IAPWS-95 (used to check the tables)."""
        self.ASx.update(CoolProp.HmassP_INPUTS, h, p)
        return self.ASx.T(), self.ASx.rhomass()


class ConstantFluid:
    """Incompressible fluid with constant properties (verification only).

    h = cp * (T - T_ref), no phase change.
    """

    is_constant = True

    def __init__(self, rho=1000.0, cp=4180.0, mu=1.0e-3, k=0.6, T_ref=273.15):
        self.rho, self.cp, self.mu, self.k, self.T_ref = rho, cp, mu, k, T_ref
        self.name = f"constant-property fluid (rho={rho}, cp={cp}, mu={mu}, k={k})"

    def Tsat(self, p):
        return np.full(np.shape(p), np.inf)

    def ph(self, p, h) -> NodeProps:
        p = np.atleast_1d(np.asarray(p, float))
        h = np.atleast_1d(np.asarray(h, float))
        n = p.size
        out = NodeProps(n)
        out.p[:], out.h[:] = p, h
        out.T[:] = self.T_ref + h / self.cp
        out.rho[:] = self.rho
        out.cp[:] = self.cp
        out.mu[:] = self.mu
        out.k[:] = self.k
        out.x[:] = 0.0
        out.dTdh[:] = 1.0 / self.cp
        return out

    def point(self, p, h):
        T = self.T_ref + h / self.cp
        return False, self.rho, self.mu, 0.0, self.rho, np.nan, self.mu, np.nan, np.nan, T

    def h_pT(self, p, T):
        p, T = np.broadcast_arrays(np.atleast_1d(p), np.atleast_1d(np.asarray(T, float)))
        return self.cp * (T - self.T_ref)

    def s_ph(self, p, h):
        T = self.T_ref + np.atleast_1d(h) / self.cp
        return self.cp * np.log(T / self.T_ref)
