"""Coupled wellbore / rock heat-transfer solver for a coaxial closed loop.

Configuration: cold fluid is injected down the annulus (between the centre
tubing and the casing / open hole) and returns up the centre tubing, or, with
flow_direction = 'tubing_down', down the centre tubing and back up the annulus.
There is no mass exchange with the formation; heat reaches the loop by
conduction only.  Internally the two streams are the "feed" (flowing down) and
the "riser" (flowing up); the annulus is whichever of them touches the casing.

Discretisation
--------------
* Along-hole (MD) axis: N cells; the first Nw are the well, the rest extend
  the rock below / beyond TD.
* Fluid: enthalpy at cell faces ("nodes") for both annulus and tubing, a
  weighted box scheme in space (second order, theta=0.5, exponential
  weighting where a cell's transfer number exceeds 1) and variable-step BDF2
  in time (implicit Euler on the first step and after large step-size jumps).  Steady-flow energy equation including gravity work:
      m dh = dQ + m g dTVD    (annulus, downward)
      m dh = dQ - m g dTVD - m d(v^2/2)  (tubing, upward)
  Temperature follows from T(p, h) (IAPWS-95), so friction heating,
  compression/expansion and flashing are captured.
* Solids (casing steel, cement, rock): axisymmetric r-z finite volumes on a
  masked structured grid, same time scheme.  The undisturbed geotherm is made an
  exact discrete steady state by an equilibrium source, so the rock does not
  drift when the well is idle.
* Momentum: steady 1-D pressure marching each iteration (gravity, Churchill
  friction, acceleration; homogeneous or drift-flux two-phase).
* Coupling: all fluid enthalpies and solid temperatures are solved together
  in one sparse linear system, with T(h) linearised (Newton for the
  enthalpy-temperature relation) and iterated with the pressure profile to
  convergence at every time step.
"""
from __future__ import annotations

import math
import time as _time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from .case import Case, DAY, YEAR
from .bhe import pipe_resistance, utube_delta
from .correlations import (G_ACC, friction_churchill, nusselt_pipe, nusselt_annulus,
                           chen_factor, two_phase_gradient)


class Results:
    """Time series and profiles produced by a run.

    Energy balance of every record:  Q_wellhead = Q_rock + Q_res - dE_fluid - Q_ke - Q_pe,
    where Q_wellhead = mdot (h_out - h_in), Q_res is the heat added by a reservoir element
    at TD (coaxial single-well EGS) and Q_pe the potential energy change of the stream
    (zero for a closed loop, + mdot g TVD for a production well, - for an injector)."""

    def __init__(self):
        self.ts = {k: [] for k in (
            "t", "mdot", "T_in", "T_out", "p_in", "p_out", "h_in", "h_out", "x_out",
            "Q_wellhead", "Q_rock", "Q_riser_loss", "dE_fluid", "iters", "flash_md", "choke_dp", "Q_ke",
            "choke_warning", "Q_res", "pump_dp", "pump_kW", "inflow_margin", "Q_pe")}
        self.profiles = []          # list of dicts at output times
        self.meta = {}

    def append(self, **kw):
        for k, v in kw.items():
            self.ts[k].append(v)

    def finalise(self):
        for k in self.ts:
            self.ts[k] = np.asarray(self.ts[k], float)
        return self

    # convenience ---------------------------------------------------------
    @property
    def t_days(self):
        return self.ts["t"] / DAY

    def at(self, t_seconds, key):
        """Linear interpolation of a time series at time t [s]."""
        return float(np.interp(t_seconds, self.ts["t"], self.ts[key]))

    def mean_over(self, t0, t1, key):
        t = self.ts["t"]
        y = self.ts[key]
        # step-wise (implicit Euler): value at t_n applies over (t_{n-1}, t_n]
        tt = np.concatenate([[0.0], t])
        dt = np.diff(tt)
        lo = np.clip(tt[:-1], t0, t1)
        hi = np.clip(tt[1:], t0, t1)
        w = np.maximum(hi - lo, 0.0)
        return float(np.sum(w * y) / max(np.sum(w), 1e-30))

    def to_csv(self, path):
        keys = list(self.ts.keys())
        arr = np.column_stack([self.ts[k] for k in keys])
        np.savetxt(path, arr, delimiter=",", header=",".join(keys), comments="")


class Simulator:
    def __init__(self, case: Case, verbose: bool = False):
        self.case = case
        self.fluid = case.fluid
        self.verbose = verbose
        self.reversed = case.flow_direction == "tubing_down"   # feed = tubing, riser = annulus
        self.mode = case.mode
        self.single = self.mode in ("production", "injection")   # one channel, touching the wall
        self.up = self.mode == "production"                       # single channel flows up
        self.utube = self.mode == "utube"
        self._build_axial()
        self._build_wellbore()
        self.grid = case.rock_model == "grid"
        if self.grid:
            self._build_radial()
            self._build_solid()
        self._init_state()

    # ======================================================================
    # Grid construction
    # ======================================================================
    def _build_axial(self):
        c = self.case
        L = c.md_total
        br = {0.0, L}
        for q in list(c.casing) + list(c.tubing):
            for v in (q.top, q.bottom):
                if 0.0 < v < L:
                    br.add(float(v))
        for v in c.survey.md:
            if 0.0 < v < L:
                br.add(float(v))
        for lay in c.rock_layers:
            for z in (lay.top, lay.bottom):
                m = c.survey.md_at_tvd(z)
                if m is not None and 0.0 < m < L:
                    br.add(float(m))
        br = sorted(br)
        merged = [br[0]]
        for v in br[1:]:
            if v - merged[-1] > 0.5:
                merged.append(v)
        merged[-1] = L
        faces = [0.0]
        self.dz_local_info = []
        for a, b in zip(merged[:-1], merged[1:]):
            dz = self._local_dz(0.5 * (a + b))
            n = max(1, int(math.ceil((b - a) / dz - 1e-9)))
            faces.extend(np.linspace(a, b, n + 1)[1:].tolist())
        s_f = np.array(faces)
        self.Nw = s_f.size - 1
        if c.rock_model == "grid":
            alpha = max(L_.k / (L_.rho * L_.cp) for L_ in c.rock_layers)
            ext = c.ext_length if c.ext_length is not None else max(200.0, 4.0 * math.sqrt(alpha * c.t_end))
            d, tot, extf = c.dz, 0.0, []
            while tot < ext - 1e-9:
                tot += d
                extf.append(L + tot)
                d *= 1.25
            s_f = np.concatenate([s_f, extf])
        self.s_f = s_f
        self.N = s_f.size - 1
        self.ds = np.diff(s_f)
        self.s_c = 0.5 * (s_f[:-1] + s_f[1:])
        self.tvd_f = c.survey.tvd_at(s_f)
        self.tvd_c = c.survey.tvd_at(self.s_c)
        self.dtvd = np.diff(self.tvd_f)
        self.cos = self.dtvd / self.ds
        self.T0_f = np.asarray(c.geotherm(self.tvd_f), float)
        self.T0_c = np.asarray(c.geotherm(self.tvd_c), float)

    def _local_dz(self, s_mid):
        """Axial cell size at MD s_mid.

        Where the centre tubing is poorly insulated, annulus and tubing streams
        exchange heat over a short counter-flow length b = m c R_t; the box
        scheme needs dz << b there, so the grid is refined locally to
        dz <= coupling_frac * b (b evaluated at the lowest non-zero flow rate
        of the schedule with nominal water properties at 80 C).
        """
        c = self.case
        if self.single or self.utube:
            return c.dz
        ts = next((q for q in c.tubing if q.top <= s_mid < q.bottom), None)
        cs = next((q for q in c.casing if q.top <= s_mid < q.bottom), None)
        flows = [p.mdot for p in c.schedule if p.mdot > 0]
        if ts is None or cs is None or not flows or c.coupling_frac <= 0:
            return c.dz
        mu, k, cp = 3.55e-4, 0.67, 4197.0
        Pr = cp * mu / k
        R = ts.wall_resistance()
        bmin = np.inf
        for m in flows:
            if c.h_wall_override is not None:
                h_t = h_a = float(c.h_wall_override)
            else:
                D = 2 * ts.r_in
                Re = 4 * m / (math.pi * D * mu)
                h_t = float(nusselt_pipe(Re, Pr)[0]) * k / D
                Dh = 2 * (cs.r_inner - ts.r_out)
                A = math.pi * (cs.r_inner ** 2 - ts.r_out ** 2)
                Re_a = m * Dh / (A * mu)
                h_a = float(nusselt_annulus(Re_a, Pr, ts.r_out / cs.r_inner, "inner")[0]) * k / Dh
            Rt = R + 1 / (h_t * 2 * math.pi * ts.r_in) + 1 / (h_a * 2 * math.pi * ts.r_out)
            bmin = min(bmin, m * cp * Rt)
        dz = min(c.dz, max(c.dz_min, c.coupling_frac * bmin))
        return dz

    def _build_wellbore(self):
        c = self.case
        Nw = self.Nw
        sc = self.s_c[:Nw]
        self.cas_idx = np.empty(Nw, int)
        self.tub_idx = np.empty(Nw, int)
        for j, s in enumerate(sc):
            self.cas_idx[j] = next(i for i, q in enumerate(c.casing) if q.top <= s < q.bottom)
            if c.mode == "coaxial":
                self.tub_idx[j] = next(i for i, q in enumerate(c.tubing) if q.top <= s < q.bottom)
        self.r_ci = np.array([c.casing[i].r_inner for i in self.cas_idx])
        self.r_hole = np.array([c.casing[i].r_hole for i in self.cas_idx])
        if self.single:
            # one flow channel of radius r_ci; everything outside it is in the solid grid
            self.r_ti = np.zeros(Nw)
            self.r_to = np.zeros(Nw)
            self.Rtw = np.full(Nw, np.inf)
            self.A_a = math.pi * self.r_ci ** 2
            self.Dh_a = 2.0 * self.r_ci
            self.A_t = self.A_a.copy()
            self.D_t = self.Dh_a.copy()
            self.a_ratio = np.zeros(Nw)
        elif self.utube:
            u = c.utube
            self.r_ti = np.full(Nw, u["r_pipe_in"])
            self.r_to = np.full(Nw, u["r_pipe_out"])
            self.Rtw = np.full(Nw, np.inf)             # leg-to-leg exchange comes from the delta circuit
            self.A_a = math.pi * self.r_ti ** 2
            self.Dh_a = 2.0 * self.r_ti
            self.A_t = self.A_a.copy()
            self.D_t = self.Dh_a.copy()
            self.a_ratio = np.zeros(Nw)
            self.x_c = 0.5 * u["shank"]
        else:
            self.r_ti = np.array([c.tubing[i].r_in for i in self.tub_idx])
            self.r_to = np.array([c.tubing[i].r_out for i in self.tub_idx])
            self.Rtw = np.array([c.tubing[i].wall_resistance() for i in self.tub_idx])
            self.A_a = math.pi * (self.r_ci ** 2 - self.r_to ** 2)
            self.Dh_a = 2.0 * (self.r_ci - self.r_to)
            self.A_t = math.pi * self.r_ti ** 2
            self.D_t = 2.0 * self.r_ti
            self.a_ratio = self.r_to / self.r_ci

    def _build_radial(self):
        c = self.case
        ifaces = sorted({round(float(r), 6) for r in self.r_ci} |
                        {round(float(q.r_inner), 6) for q in c.casing} |
                        {round(float(r_out), 6) for q in c.casing for r_out, _ in q.layers})
        r0 = min(self.r_ci)
        ifaces = [r for r in ifaces if r >= r0 - 1e-9]
        g = c.radial_growth
        faces = [ifaces[0]]
        for a, b in zip(ifaces[:-1], ifaces[1:]):
            n = 1 if (b - a) < 0.012 else max(1, int(math.ceil(math.log(b / a) / math.log(g))))
            faces.extend(np.geomspace(a, b, n + 1)[1:].tolist())
        alpha = max(L_.k / (L_.rho * L_.cp) for L_ in c.rock_layers)
        r_far = c.r_far if c.r_far is not None else min(2000.0, max(30.0, 4.0 * math.sqrt(alpha * c.t_end)))
        rb = faces[-1]
        n = max(4, int(math.ceil(math.log(r_far / rb) / math.log(g))))
        faces.extend(np.geomspace(rb, r_far, n + 1)[1:].tolist())
        self.r_f = np.array(faces)
        self.M = self.r_f.size - 1
        self.r_c = np.sqrt(self.r_f[:-1] * self.r_f[1:])
        self.A_ring = math.pi * (self.r_f[1:] ** 2 - self.r_f[:-1] ** 2)

    def _build_solid(self):
        c = self.case
        N, M, Nw = self.N, self.M, self.Nw
        active = np.zeros((N, M), bool)
        kmap = np.zeros((N, M))
        rhoc = np.zeros((N, M))
        self.i0 = np.zeros(Nw, int)
        for j in range(N):
            rock = c.rock_at(self.tvd_c[j])
            if j < Nw:
                cs = c.casing[self.cas_idx[j]]
                self.i0[j] = int(np.argmin(np.abs(self.r_f - cs.r_inner)))
                for i in range(M):
                    rc = self.r_c[i]
                    if rc < cs.r_inner:
                        continue
                    active[j, i] = True
                    mat = None
                    r_prev = cs.r_inner
                    for r_out, m_ in cs.layers:
                        if r_prev <= rc < r_out:
                            mat = m_
                            break
                        r_prev = r_out
                    if mat is None:
                        kmap[j, i], rhoc[j, i] = rock.k, rock.rho * rock.cp
                    else:
                        kmap[j, i], rhoc[j, i] = mat.k, mat.rho * mat.cp
            else:
                active[j, :] = True
                kmap[j, :], rhoc[j, :] = rock.k, rock.rho * rock.cp
        sid = -np.ones((N, M), int)
        sid[active] = np.arange(active.sum())
        self.active, self.kmap, self.rhoc, self.sid = active, kmap, rhoc, sid
        self.Ns = int(active.sum())
        ds = self.ds
        rows, cols, vals = [], [], []

        def add_pair(p, q, G):
            rows.extend([p, q, p, q]); cols.extend([p, q, q, p]); vals.extend([G, G, -G, -G])

        diag = np.zeros(self.Ns)
        # radial pairs
        ln1 = np.log(self.r_f[1:-1] / self.r_c[:-1])
        ln2 = np.log(self.r_c[1:] / self.r_f[1:-1])
        both = active[:, :-1] & active[:, 1:]
        with np.errstate(divide="ignore", invalid="ignore"):
            Rr = (ln1[None, :] / kmap[:, :-1] + ln2[None, :] / kmap[:, 1:]) / (2 * math.pi * ds[:, None])
        jj, ii = np.nonzero(both)
        G = 1.0 / Rr[jj, ii]
        p, q = sid[jj, ii], sid[jj, ii + 1]
        rows_r = np.concatenate([p, q, p, q]); cols_r = np.concatenate([p, q, q, p])
        vals_r = np.concatenate([G, G, -G, -G])
        # axial pairs
        bothz = active[:-1, :] & active[1:, :]
        with np.errstate(divide="ignore", invalid="ignore"):
            Rz = (0.5 * ds[:-1, None] / kmap[:-1, :] + 0.5 * ds[1:, None] / kmap[1:, :]) / self.A_ring[None, :]
        jj, ii = np.nonzero(bothz)
        G = 1.0 / Rz[jj, ii]
        p, q = sid[jj, ii], sid[jj + 1, ii]
        rows_z = np.concatenate([p, q, p, q]); cols_z = np.concatenate([p, q, q, p])
        vals_z = np.concatenate([G, G, -G, -G])
        # boundaries (Dirichlet at undisturbed temperature; the outer radius may be adiabatic = unit cell)
        if c.far_field == "dirichlet":
            jj = np.nonzero(active[:, M - 1])[0]
            Gout = 2 * math.pi * ds[jj] * kmap[jj, M - 1] / math.log(self.r_f[M] / self.r_c[M - 1])
            np.add.at(diag, sid[jj, M - 1], Gout)
        ii = np.nonzero(active[0, :])[0]
        np.add.at(diag, sid[0, ii], kmap[0, ii] * self.A_ring[ii] / (0.5 * ds[0]))
        ii = np.nonzero(active[N - 1, :])[0]
        np.add.at(diag, sid[N - 1, ii], kmap[N - 1, ii] * self.A_ring[ii] / (0.5 * ds[N - 1]))
        Ks = sp.coo_matrix((np.concatenate([vals_r, vals_z, diag]),
                            (np.concatenate([rows_r, rows_z, np.arange(self.Ns)]),
                             np.concatenate([cols_r, cols_z, np.arange(self.Ns)]))),
                           shape=(self.Ns, self.Ns)).tocsr()
        self.Ks = Ks
        jgrid = np.repeat(np.arange(N)[:, None], M, axis=1)
        self.T0_s = self.T0_c[jgrid[active]]
        self.C_s = (rhoc * self.A_ring[None, :] * ds[:, None])[active]
        self.rhs_eq = Ks @ self.T0_s          # = b_bc + S  (equilibrium source)
        Ks_coo = Ks.tocoo()
        self._Ks_rows, self._Ks_cols, self._Ks_vals = Ks_coo.row, Ks_coo.col, Ks_coo.data
        # coupling: annulus <-> first solid cell
        self.sid_wall = np.array([sid[j, self.i0[j]] for j in range(Nw)])
        kw = np.array([kmap[j, self.i0[j]] for j in range(Nw)])
        self.Rhalf = np.log(self.r_c[self.i0] / self.r_ci) / (2 * math.pi * kw)
        # completion resistance (casing + cement, r_ci -> r_hole) and the rock
        # properties at the borehole wall (used only for the scheme weights)
        self.R_comp = np.zeros(Nw)
        self.k_wall_rock = np.zeros(Nw)
        self.alpha_wall = np.zeros(Nw)
        for j in range(Nw):
            cs = c.casing[self.cas_idx[j]]
            r = cs.r_inner
            for r_out, mat in cs.layers:
                self.R_comp[j] += math.log(r_out / r) / (2 * math.pi * mat.k)
                r = r_out
            rk = c.rock_at(self.tvd_c[j])
            self.k_wall_rock[j] = rk.k
            self.alpha_wall[j] = rk.k / (rk.rho * rk.cp)
        # bottom-of-well coupling (rock below TD under the fluid column)
        self.bot_sid, self.bot_R = np.array([], int), np.array([])
        if N > Nw:
            jb = Nw
            ib = np.nonzero(self.r_f[1:] <= self.r_ci[Nw - 1] + 1e-9)[0]
            if ib.size:
                self.bot_sid = sid[jb, ib]
                self.bot_A = self.A_ring[ib]
                self.bot_R = 0.5 * ds[jb] / (kmap[jb, ib] * self.A_ring[ib])

    # ======================================================================
    # State
    # ======================================================================
    def _init_state(self):
        c = self.case
        Nw = self.Nw
        f = self.fluid
        T0n = self.T0_f[:Nw + 1]
        # static pressure with the undisturbed temperature
        p_top = c.pressure_bc.get("outlet", c.pressure_bc.get("inlet"))
        if p_top is None and "bottom" in c.pressure_bc:       # rough static estimate from the bottom pressure
            p_top = max(1.0e5, c.pressure_bc["bottom"] - 900.0 * G_ACC * (self.tvd_f[Nw] - self.tvd_f[0]))
        p = np.empty(Nw + 1)
        p[0] = p_top
        for j in range(Nw):
            rho = f.point(p[j], f.h_pT(p[j], T0n[j])[0])[1]
            p1 = p[j] + rho * G_ACC * self.dtvd[j]
            rho1 = f.point(p1, f.h_pT(p1, T0n[j + 1])[0])[1]
            p[j + 1] = p[j] + 0.5 * (rho + rho1) * G_ACC * self.dtvd[j]
        self.p_a = p.copy()
        self.p_t = p.copy()
        self.h_a = f.h_pT(p, T0n)
        self.h_t = self.h_a.copy()
        self.T_s = self.T0_s.copy() if self.grid else None
        self.dKE = np.zeros(Nw)
        self.dKE_used = np.zeros(Nw)
        self._hist = None
        self.t = 0.0
        if self.mode == "production":
            # warm start: the column holds the feed fluid and the pressures of the flowing well (the rock is
            # still undisturbed), so the run does not begin by filling a steam-filled hole from the bottom
            m0 = next((q.mdot for q in c.schedule if q.mdot > 0), 0.0)
            if m0 > 0:
                self._T_in_current = c.schedule[0].T_in
                self._m_current = m0
                h_feed = self._feed_enthalpy(p[Nw], 0.0)
                self.h_a = np.full(Nw + 1, h_feed)
                self.h_t = self.h_a.copy()
                try:
                    p_flow, _, _ = self._march_single(self.h_a, m0)
                    self.p_a = p_flow.copy()
                    self.p_t = p_flow.copy()
                except RuntimeError:
                    pass
        self.res = Results()
        self.res.meta = dict(name=c.name, Nw=Nw, N=self.N, M=getattr(self, "M", 0),
                             Ns=getattr(self, "Ns", 0), fluid=getattr(f, "name", ""))

    # ======================================================================
    # Heat-transfer coefficients
    # ======================================================================
    def _node_htc(self, nodes, m, side):
        """Film coefficient(s) per cell for 'annulus' (outer wall, inner wall) or 'tubing'."""
        Nw = self.Nw
        npr = nodes
        if side == "annulus":
            G = m / self.A_a
            D = self.Dh_a
        else:
            G = m / self.A_t
            D = self.D_t

        def node_vals(sl):
            mu, k, cp = npr.mu[sl], npr.k[sl], npr.cp[sl]
            tp = npr.two_phase[sl]
            return mu, k, cp, tp

        hs = []
        for sl in (slice(0, Nw), slice(1, Nw + 1)):
            mu, k, cp, tp = node_vals(sl)
            mu = np.where(tp, npr.mu_l[sl], mu)
            k = np.where(tp, npr.k_l[sl], k)
            cp = np.where(tp, npr.cp_l[sl], cp)
            x = np.where(tp, npr.x[sl], 0.0)
            Re = np.abs(G) * D * np.where(tp, 1.0 - x, 1.0) / mu
            Pr = cp * mu / k
            if side == "annulus" and (self.single or self.utube):
                h_o = nusselt_pipe(Re, Pr) * k / D
                h_i = h_o.copy()
                if np.any(tp):
                    F = np.ones(Nw)
                    F[tp] = chen_factor(x[tp], npr.rho_l[sl][tp], npr.rho_g[sl][tp],
                                        npr.mu_l[sl][tp], npr.mu_g[sl][tp])
                    h_o, h_i = h_o * F, h_i * F
                hs.append((h_o, h_i))
            elif side == "annulus":
                h_o = nusselt_annulus(Re, Pr, self.a_ratio, "outer") * k / D
                h_i = nusselt_annulus(Re, Pr, self.a_ratio, "inner") * k / D
                if np.any(tp):
                    F = np.ones(Nw)
                    F[tp] = chen_factor(x[tp], npr.rho_l[sl][tp], npr.rho_g[sl][tp],
                                        npr.mu_l[sl][tp], npr.mu_g[sl][tp])
                    h_o, h_i = h_o * F, h_i * F
                hs.append((h_o, h_i))
            else:
                h = nusselt_pipe(Re, Pr) * k / D
                if np.any(tp):
                    F = np.ones(Nw)
                    F[tp] = chen_factor(x[tp], npr.rho_l[sl][tp], npr.rho_g[sl][tp],
                                        npr.mu_l[sl][tp], npr.mu_g[sl][tp])
                    h = h * F
                hs.append(h)
        if side == "annulus":
            return 0.5 * (hs[0][0] + hs[1][0]), 0.5 * (hs[0][1] + hs[1][1])
        return 0.5 * (hs[0] + hs[1])

    # ======================================================================
    # Pressure marching
    # ======================================================================
    def _grad_terms(self, p, h, G, D, cos, side):
        """(rho_grav, fric_grad, rho_mom, two_phase, x) at a node."""
        if not (1.0e3 <= p <= 400.0e6):
            where = "annulus" if side == "a" else "tubing"
            raise RuntimeError(
                f"pressure in the {where} left the physical range ({p / 1e6:.3g} MPa). The frictional pressure "
                f"drop is too large for this flow rate and pipe size, or the inlet pressure is too low to "
                f"push the flow round the loop (the riser would choke). Reduce the flow rate, enlarge the "
                f"tubing/annulus, or change the pressure boundary condition.")
        tp, rho, mu, x, rl, rg, ml, mg, sig, T = self.fluid.point(p, h)
        if tp:
            model = self.case.two_phase_model
            rho_g_, fric, rho_m, alpha = two_phase_gradient(model, G, D, x, rl, rg, ml, mg, sig,
                                                            self.case.roughness, cos)
            return float(rho_g_), float(fric), float(rho_m), True, x
        Re = abs(G) * D / mu
        f = friction_churchill(Re, self.case.roughness / D)
        return rho, float(f * G * abs(G) / (2.0 * D * rho)), rho, False, x

    def _march(self, h_a, h_t, m):
        """Pressure profiles for the current enthalpies.  The riser (up-flowing
        stream) is marched down from the outlet pressure, then the feed stream
        up to the inlet; if the hot riser drives the loop harder than the pump
        needs, the inlet is held at its minimum pressure and the march is redone
        from the inlet (the excess goes into a wellhead choke)."""
        c = self.case
        Nw = self.Nw
        if not c.include_pressure:
            return self.p_a.copy(), self.p_t.copy(), np.zeros(Nw), None
        rev = self.reversed
        # stream descriptors: enthalpies, mass flux, diameter, old pressures (guesses), direction
        ann = dict(h=h_a.copy(), G=m / self.A_a, D=self.Dh_a, side="a", p_old=self.p_a, up=rev)
        tub = dict(h=h_t.copy(), G=m / self.A_t, D=self.D_t, side="t", p_old=self.p_t, up=not rev)
        F, R = (tub, ann) if rev else (ann, tub)          # feed (down), riser (up)
        R["h"][Nw] = F["h"][Nw]                           # the riser starts with the feed's bottom enthalpy
        pF = np.empty(Nw + 1)
        pR = np.empty(Nw + 1)
        rho_mom = np.empty(Nw + 1)                        # riser momentum density at the nodes
        conv = {"ok": True, "choke": False}

        def cell(st, j, p_known, known_is_top, p_floor=None, riser=False):
            """Unknown end pressure of cell j of stream st.  p_bot - p_top =
            gravity +/- friction (+ for up-flow) + momentum.  With p_floor set,
            returns None as soon as a trial pressure falls below it, so no
            property is evaluated at an unphysical pressure."""
            h_top, h_bot = st["h"][j], st["h"][j + 1]
            G, D, side = st["G"][j], st["D"][j], st["side"]
            sgn = 1.0 if st["up"] else -1.0
            st_k = self._grad_terms(p_known, h_top if known_is_top else h_bot, G, D, self.cos[j], side)
            p_prev = st["p_old"][j + 1] if known_is_top else st["p_old"][j]
            if np.isfinite(p_prev) and p_prev > 1.0e3:
                p_u = p_prev
            else:
                dp0 = st_k[0] * G_ACC * self.dtvd[j] + sgn * st_k[1] * self.ds[j]
                p_u = p_known + dp0 if known_is_top else p_known - dp0
            top = bot = None
            for _ in range(15):
                if p_floor is not None and p_u < p_floor:
                    return None, None, None
                st_u = self._grad_terms(max(p_u, 1000.0), h_bot if known_is_top else h_top, G, D, self.cos[j], side)
                top, bot = (st_k, st_u) if known_is_top else (st_u, st_k)
                dp = 0.5 * (top[0] + bot[0]) * G_ACC * self.dtvd[j] + sgn * 0.5 * (top[1] + bot[1]) * self.ds[j] \
                    + G * G * (1.0 / top[2] - 1.0 / bot[2])
                p_new = p_known + dp if known_is_top else p_known - dp
                if abs(p_new - p_u) < 1.0:
                    p_u = p_new
                    break
                p_u = 0.5 * (p_u + p_new) if (top[3] or bot[3]) else p_new
            else:
                conv["ok"] = False
            if p_floor is not None and p_u < p_floor:
                return None, None, None
            if not (1.0e3 <= p_u <= 400.0e6):
                self._grad_terms(p_u, h_top, G, D, self.cos[j], side)   # raises a clear error
            # near-critical (choking) two-phase flow in the riser: a large fraction of
            # the local pressure is lost over one cell, or the mixture velocity is very high
            if riser:
                p_lo_end = min(p_u, p_known)
                if (top[3] or bot[3]) and (abs(p_u - p_known) > 0.25 * p_lo_end or G / top[2] > 120.0):
                    conv["choke"] = True
            return p_u, top, bot

        turn = lambda rho: c.turn_loss_K * 0.5 * (m / self.A_t[-1]) ** 2 / rho
        state = {"flash": None}

        def riser_cell(j, p_known, known_is_top):
            p, top, bot = cell(R, j, p_known, known_is_top, riser=True)
            rho_mom[j] = top[2]
            rho_mom[j + 1] = bot[2]
            if top[3] and not bot[3]:
                state["flash"] = self.s_f[j + 1]
            return p

        def from_outlet(p_out, p_floor):
            """March down the riser and up the feed; returns False if the feed
            pressure falls below p_floor (loop must be throttled)."""
            pR[0] = p_out
            state["flash"] = None
            for j in range(Nw):
                pR[j + 1] = riser_cell(j, pR[j], True)
            rho_bot = self.fluid.point(pR[Nw], F["h"][Nw])[1]
            pF[Nw] = pR[Nw] + turn(rho_bot)
            for j in range(Nw - 1, -1, -1):
                pj, _, _ = cell(F, j, pF[j + 1], False, p_floor=p_floor)
                if pj is None:
                    return False
                pF[j] = pj
            return True

        def from_inlet(p_in):
            pF[0] = p_in
            state["flash"] = None
            for j in range(Nw):
                pF[j + 1], _, _ = cell(F, j, pF[j], True)
            rho_bot = self.fluid.point(pF[Nw], F["h"][Nw])[1]
            pR[Nw] = pF[Nw] - turn(rho_bot)
            for j in range(Nw - 1, -1, -1):
                pR[j] = riser_cell(j, pR[j + 1], False)

        self.choke_dp = 0.0
        if "outlet" in c.pressure_bc:
            p_spec = c.pressure_bc["outlet"]
            # Thermosiphon stronger than needed to drive the specified flow: the
            # inlet would need (near-)vacuum.  Physically the flow must then be
            # throttled; the inlet is held at its minimum pressure and a wellhead
            # choke takes the excess (isenthalpic, so the heat delivered is unchanged).
            p_min = c.min_inlet_p
            if not self.fluid.is_constant:
                p_min = max(p_min, 1.2 * self._psat(self._T_in_current))
            ok = from_outlet(p_spec, p_min)
            if not ok or np.min(pF) < p_min:
                conv["ok"], conv["choke"] = True, False     # the trial march is discarded
                from_inlet(p_min)
                self.choke_dp = max(0.0, pR[0] - p_spec)
        else:
            from_inlet(c.pressure_bc["inlet"])
        flash_md = state["flash"]
        GR = R["G"]
        v = np.concatenate([GR / rho_mom[:Nw], [GR[-1] / rho_mom[Nw]]])
        dKE = 0.5 * (v[:-1] ** 2 - v[1:] ** 2)               # kinetic-energy gain up each riser cell
        self._march_ok = conv["ok"]
        self.choke_warning = conv["choke"]
        p_a, p_t = (pR, pF) if rev else (pF, pR)
        return p_a, p_t, dKE, flash_md

    # ======================================================================
    # Single-channel wells (production: flow up from a feed zone at TD;
    # injection: flow down into the formation at TD).  The channel is stored in
    # the 'annulus' arrays (h_a, p_a) because it is the stream touching the wall.
    # ======================================================================
    def _feed_enthalpy(self, p, t):
        """Enthalpy of the fluid entering a production well at TD.  A feed given by temperature is liquid at the
        reservoir condition: it is evaluated at the reservoir pressure less the drawdown when p_res (and PI) are
        given, otherwise at the larger of the local pressure and the saturation pressure; it then flashes
        isenthalpically in the well.  A two-phase or steam feed is given by enthalpy."""
        inf = self.case.inflow
        if "h" in inf:
            hv = inf["h"]
            return float(hv(t)) if callable(hv) else float(hv)
        Tv = inf["T"]
        T = float(Tv(t)) if callable(Tv) else float(Tv)
        if "p_res" in inf:
            p_eval = inf["p_res"] - (self._m_current / inf["PI"] if "PI" in inf else 0.0)
        else:
            p_eval = p
        if not self.fluid.is_constant:
            p_eval = max(p_eval, 1.02 * self._psat(T))
        return float(self.fluid.h_pT(p_eval, T)[0])

    def _march_single(self, h, m):
        """Pressure profile of the single channel.  Production: march down from the
        wellhead pressure; with a downhole pump the lower part is marched up from the
        flowing bottom-hole pressure p_res - m/PI and the pump supplies the difference.
        Injection: march down from the wellhead injection pressure.  Returns
        (p, dKE, flash_md) and sets pump_dp, pump_kW, inflow_margin, choke_dp."""
        c = self.case
        Nw = self.Nw
        self.pump_dp = 0.0
        self.pump_kW = 0.0
        self.inflow_margin = np.nan
        self.choke_dp = 0.0
        self.choke_warning = False
        self._march_ok = True
        if not c.include_pressure:
            return self.p_a.copy(), np.zeros(Nw), None
        st = dict(h=h.copy(), G=m / self.A_a, D=self.Dh_a, side="a", p_old=self.p_a, up=self.up)
        pp = np.empty(Nw + 1)
        rho_mom = np.ones(Nw + 1)
        conv = {"ok": True, "choke": False}
        state = {"flash": None}

        def cell(j, p_known, known_is_top, p_floor=None):
            h_top, h_bot = st["h"][j], st["h"][j + 1]
            G, D = st["G"][j], st["D"][j]
            sgn = 1.0 if st["up"] else -1.0
            st_k = self._grad_terms(p_known, h_top if known_is_top else h_bot, G, D, self.cos[j], "a")
            p_prev = st["p_old"][j + 1] if known_is_top else st["p_old"][j]
            if np.isfinite(p_prev) and p_prev > 1.0e3:
                p_u = p_prev
            else:
                dp0 = st_k[0] * G_ACC * self.dtvd[j] + sgn * st_k[1] * self.ds[j]
                p_u = p_known + dp0 if known_is_top else p_known - dp0
            top = bot = None
            for _ in range(15):
                if p_floor is not None and p_u < p_floor:
                    return None
                st_u = self._grad_terms(max(p_u, 1000.0), h_bot if known_is_top else h_top, G, D, self.cos[j], "a")
                top, bot = (st_k, st_u) if known_is_top else (st_u, st_k)
                dp = 0.5 * (top[0] + bot[0]) * G_ACC * self.dtvd[j] + sgn * 0.5 * (top[1] + bot[1]) * self.ds[j] \
                    + G * G * (1.0 / top[2] - 1.0 / bot[2])
                p_new = p_known + dp if known_is_top else p_known - dp
                if abs(p_new - p_u) < 1.0:
                    p_u = p_new
                    break
                p_u = 0.5 * (p_u + p_new) if (top[3] or bot[3]) else p_new
            else:
                conv["ok"] = False
            if p_floor is not None and p_u < p_floor:
                return None
            if not (1.0e3 <= p_u <= 400.0e6):
                self._grad_terms(p_u, h_top, G, D, self.cos[j], "a")
            if st["up"]:
                rho_mom[j], rho_mom[j + 1] = top[2], bot[2]
                if top[3] and not bot[3]:
                    state["flash"] = self.s_f[j + 1]
                p_lo_end = min(p_u, p_known)
                if (top[3] or bot[3]) and (abs(p_u - p_known) > 0.25 * p_lo_end or G / top[2] > 120.0):
                    conv["choke"] = True
            return p_u

        inf = c.inflow or {}
        if self.up and "bottom" in c.pressure_bc:
            pp[Nw] = c.pressure_bc["bottom"]
            for j in range(Nw - 1, -1, -1):
                pj = cell(j, pp[j + 1], False, p_floor=1.0e5)
                if pj is None:                      # the column cannot reach the surface at this bottom pressure
                    pp[:j + 1] = 1.0e5
                    conv["ok"] = False
                    break
                pp[j] = pj
            G = st["G"]
            v = np.concatenate([G / rho_mom[:Nw], [G[-1] / rho_mom[Nw]]])
            dKE = 0.5 * (v[:-1] ** 2 - v[1:] ** 2)
        elif self.up:
            p_wh = c.pressure_bc["outlet"]
            pp[0] = p_wh
            if c.pump is not None:
                jp = int(np.clip(np.searchsorted(self.s_f, c.pump["md"]), 1, Nw - 1))   # pump at face jp
                for j in range(jp):
                    pp[j + 1] = cell(j, pp[j], True)
                p_disch = pp[jp]
                p_wf = inf["p_res"] - m / inf["PI"]
                pp[Nw] = p_wf
                self.pump_dry = False
                for j in range(Nw - 1, jp - 1, -1):
                    pj = cell(j, pp[j + 1], False, p_floor=1.0e5)
                    if pj is None:          # fluid level below the pump: the pump would run dry
                        self.pump_dry = True
                        pp[jp:j + 1] = 1.0e5
                        break
                    pp[j] = pj
                p_suct = pp[jp]
                self.pump_dp = p_disch - p_suct
                rho_p = self.fluid.point(max(p_suct, 1e3), h[jp])[1]
                self.pump_kW = m * max(self.pump_dp, 0.0) / rho_p / c.pump["eta"] / 1e3
                self.pump_two_phase = bool(self.fluid.point(max(p_suct, 1e3), h[jp])[0]) or self.pump_dry
                pp[jp] = p_disch          # node jp holds the discharge pressure (cell jp-1 above it)
                self.p_suction = p_suct
                self.inflow_margin = 0.0
            else:
                for j in range(Nw):
                    pp[j + 1] = cell(j, pp[j], True)
                if "p_res" in inf and "PI" in inf:
                    p_wf_avail = inf["p_res"] - m / inf["PI"]
                    self.inflow_margin = p_wf_avail - pp[Nw]     # > 0: a wellhead choke takes the excess
                    self.choke_dp = max(0.0, self.inflow_margin)
            G = st["G"]
            v = np.concatenate([G / rho_mom[:Nw], [G[-1] / rho_mom[Nw]]])
            dKE = 0.5 * (v[:-1] ** 2 - v[1:] ** 2)
        else:
            pp[0] = c.pressure_bc.get("inlet", c.pressure_bc.get("outlet"))
            for j in range(Nw):
                pp[j + 1] = cell(j, pp[j], True)
            if "p_res" in inf and "PI" in inf:
                self.inflow_margin = pp[Nw] - (inf["p_res"] + m / inf["PI"])   # > 0: injects freely
            dKE = np.zeros(Nw)
        self._march_ok = conv["ok"]
        self.choke_warning = conv["choke"]
        return pp, dKE, state["flash"]

    def _step_single(self, dt, period):
        c = self.case
        f = self.fluid
        Nw = self.Nw
        m = period.mdot
        up = self.up
        self._T_in_current = period.T_in
        self._m_current = m
        gw = 1.0 if c.gravity_work else 0.0
        h_old = self.h_a.copy()
        T_s_old = self.T_s.copy() if self.grid else None
        prev = self._hist
        w_ = dt / prev["dt"] if prev is not None else np.inf
        if c.time_scheme == "bdf2" and prev is not None and w_ <= 2.4:
            a0, a1, a2 = (1 + 2 * w_) / (1 + w_), 1 + w_, w_ * w_ / (1 + w_)
            h_nm1, T_s_nm1 = prev["h_a"], prev["T_s"]
        else:
            a0, a1, a2 = 1.0, 1.0, 0.0
            h_nm1, T_s_nm1 = h_old, T_s_old
        h = self.h_a.copy()
        pp = self.p_a.copy()
        T_s = T_s_old.copy() if self.grid else None
        ds = self.ds[:Nw]
        dtvd = self.dtvd[:Nw]
        flash_md = None
        it_used = 0
        converged = False
        j = np.arange(Nw)
        # unknown columns: up -> nodes 0..Nw-1 (node Nw = feed, known); down -> nodes 1..Nw (node 0 known)
        if up:
            col_j, col_j1 = j, np.where(j + 1 < Nw, j + 1, -1)
        else:
            col_j, col_j1 = j - 1, j
        sgn = -1.0 if up else 1.0          # advection: sgn * m * (h[j+1] - h[j])
        t_new = self.t + dt
        omega, d_prev = 0.6, np.zeros(Nw + 1)
        for it in range(c.max_iter):
            it_used = it + 1
            if up:
                h_known = self._feed_enthalpy(pp[Nw], t_new)
                h[Nw] = h_known
            else:
                inf_ = c.inflow or {}
                if "h" in inf_ or "T" in inf_:
                    h_known = self._feed_enthalpy(pp[0], t_new)
                else:
                    h_known = float(f.h_pT(pp[0], period.T_in)[0])
                h[0] = h_known
            nn = f.ph(pp, h)
            Dn = nn.dTdh.copy()
            cn = nn.T - Dn * h
            if c.h_wall_override is not None:
                h_w = np.full(Nw, float(c.h_wall_override))
            else:
                h_w, _ = self._node_htc(nn, m, "annulus")
            if self.grid:
                Gaw = ds / (1.0 / (h_w * 2 * math.pi * self.r_ci) + self.Rhalf)
            else:
                Gaw = ds / (1.0 / (h_w * 2 * math.pi * self.r_ci) + c.R_formation)
            rho_c = 0.5 * (nn.rho[:-1] + nn.rho[1:])
            Cc = rho_c * self.A_a * ds / dt if (c.fluid_storage and c.single_storage) else np.zeros(Nw)
            if it == 0:
                if c.scheme == "upwind" or m <= 0.0:
                    th = np.ones(Nw)
                else:
                    Dc = 0.5 * (Dn[:-1] + Dn[1:])
                    R_film = 1.0 / (h_w * 2 * math.pi * self.r_ci)
                    if self.grid:
                        tD = np.maximum(self.alpha_wall * max(t_new, dt) / self.r_hole ** 2, 1e-6)
                        f_t = np.where(tD <= 1.5, 1.1281 * np.sqrt(tD) * (1 - 0.3 * np.sqrt(tD)),
                                       (0.4063 + 0.5 * np.log(tD)) * (1 + 0.6 / tD))
                        R_path = R_film + self.R_comp + f_t / (2 * math.pi * self.k_wall_rock)
                    else:
                        R_path = R_film + c.R_formation
                    lam = np.maximum((ds / R_path) * Dc / m, 1e-12)
                    th = np.where(lam <= 1.0, 0.5, 1.0 / (1.0 - np.exp(-np.minimum(lam, 700))) - 1.0 / lam)
                w_j, w_j1 = (th, 1 - th) if up else (1 - th, th)
            Hold = w_j * h_old[:-1] + w_j1 * h_old[1:]
            Hnm1 = w_j * h_nm1[:-1] + w_j1 * h_nm1[1:]
            Sc = a1 * Hold - a2 * Hnm1
            Tbar_c = w_j * cn[:-1] + w_j1 * cn[1:]
            rows, cols, vals = [], [], []
            rhs = np.zeros(Nw + (self.Ns if self.grid else 0))
            rows += [j, j]
            cols += [col_j1, col_j]
            vals += [a0 * Cc * w_j1 + sgn * m + Gaw * w_j1 * Dn[1:],
                     a0 * Cc * w_j - sgn * m + Gaw * w_j * Dn[:-1]]
            rhs[j] += Cc * Sc + sgn * m * G_ACC * dtvd * gw - (m * self.dKE if up else 0.0) - Gaw * Tbar_c
            self.dKE_used = self.dKE.copy()
            if self.grid:
                off = Nw
                sw = off + self.sid_wall
                rows.append(j); cols.append(sw); vals.append(-Gaw)
                rows.append(sw); cols.append(sw); vals.append(Gaw)
                rows += [sw, sw]
                cols += [col_j1, col_j]
                vals += [-Gaw * w_j1 * Dn[1:], -Gaw * w_j * Dn[:-1]]
                rhs[sw] += Gaw * Tbar_c
                rows.append(off + self._Ks_rows); cols.append(off + self._Ks_cols); vals.append(self._Ks_vals)
                rows.append(off + np.arange(self.Ns)); cols.append(off + np.arange(self.Ns)); vals.append(a0 * self.C_s / dt)
                rhs[off:] += self.C_s / dt * (a1 * T_s_old - a2 * T_s_nm1) + self.rhs_eq
            else:
                rhs[j] += Gaw * self.T0_c[:Nw]
            rows = np.concatenate([np.atleast_1d(x) for x in rows])
            cols = np.concatenate([np.atleast_1d(x) for x in cols])
            vals = np.concatenate([np.atleast_1d(x) for x in vals])
            known = cols == -1
            np.add.at(rhs, rows[known], -vals[known] * h_known)
            keep = ~known
            rows, cols, vals = rows[keep], cols[keep], vals[keep]
            nunk = rhs.size
            A = sp.csc_matrix((vals, (rows, cols)), shape=(nunk, nunk))
            x = spla.spsolve(A, rhs)
            if not np.all(np.isfinite(x)):
                raise RuntimeError("linear solve failed (non-finite solution)")
            h_new = h.copy()
            if up:
                h_new[:Nw] = x[:Nw]
            else:
                h_new[1:] = x[:Nw]
            dh = float(np.max(np.abs(h_new - h)))
            # Safeguards for flashing / condensing columns.  T(h) is linearised about the current state, so a cold
            # casing against steam can send the solve across the phase change; and the saturation temperature of a
            # two-phase column depends on the pressure march, which depends on the enthalpy, so the plain
            # fixed-point iteration can oscillate.  Limit the step, under-relax when two-phase states are present
            # (halving the factor when the correction reverses), and keep the enthalpy inside the physical range.
            d_now = h_new - h
            big = np.max(np.abs(d_now))
            if big > 3.0e5:
                d_now *= 3.0e5 / big
            if np.any(nn.two_phase) or np.any(nn.x >= 1.0):
                if it > 0 and np.dot(d_now, d_prev) < 0:
                    omega = max(0.25, 0.5 * omega)
                h_new = h + omega * d_now
            else:
                h_new = h + d_now
                omega = 1.0 if it == 0 else omega
            d_prev = d_now
            if not f.is_constant:
                h_lo = f.h_pT(pp, np.full(Nw + 1, 274.0))
                h_new = np.maximum(h_new, h_lo)
            dT = 0.0
            if self.grid:
                T_s_new = x[Nw:]
                dT = float(np.max(np.abs(T_s_new - T_s)))
                T_s = T_s_new
            h = h_new
            p_new, dKE, flash_md = self._march_single(h, m)
            dp = float(np.max(np.abs(p_new - pp)))
            pp = p_new
            self.p_a = pp
            self.dKE = dKE
            converged = dh < c.tol_h and dT < c.tol_T and dp < c.tol_p and it >= 1 and self._march_ok
            if converged:
                break
        self.h_a, self.p_a = h, pp
        self.h_t, self.p_t = h, pp
        if self.grid:
            self.T_s = T_s
        nn = f.ph(pp, h)
        T_lin = Dn * h + cn
        Tbar = w_j * T_lin[:-1] + w_j1 * T_lin[1:]
        Tw = T_s[self.sid_wall] if self.grid else self.T0_c[:Nw]
        q_rock = Gaw * (Tw - Tbar)
        dE = float(np.sum(Cc * (a0 * (w_j * h[:-1] + w_j1 * h[1:]) - Sc)))
        self._hist = dict(h_a=h_old, h_t=h_old, T_s=T_s_old, dt=dt)
        self._last = dict(na=nn, nt=nn, Ta_bar=Tbar, Tt_bar=Tbar, Tw=Tw, q_rock=q_rock / ds, q_at=np.zeros(Nw),
                          h_ao=h_w, h_ai=h_w, h_ti=h_w, Gaw=Gaw, Gat=np.zeros(Nw))
        i_out, i_in = (0, Nw) if up else (Nw, 0)
        return dict(h_in=h_known, Q_rock=float(np.sum(q_rock)), Q_at=0.0, dE=dE,
                    Q_ke=float(m * np.sum(self.dKE_used)), iters=it_used, flash_md=flash_md,
                    T_out=float(nn.T[i_out]), x_out=float(nn.x[i_out]), T_in_actual=float(nn.T[i_in]),
                    converged=converged, Q_res=0.0)

    def _save_state(self):
        return dict(h_a=self.h_a.copy(), h_t=self.h_t.copy(), p_a=self.p_a.copy(), p_t=self.p_t.copy(),
                    T_s=None if self.T_s is None else self.T_s.copy(), hist=self._hist, dKE=self.dKE.copy(),
                    choke_dp=getattr(self, "choke_dp", 0.0))

    def _restore_state(self, st):
        self.h_a, self.h_t, self.p_a, self.p_t = st["h_a"].copy(), st["h_t"].copy(), st["p_a"].copy(), st["p_t"].copy()
        self.T_s = None if st["T_s"] is None else st["T_s"].copy()
        self._hist, self.dKE, self.choke_dp = st["hist"], st["dKE"].copy(), st["choke_dp"]

    def _psat(self, T):
        """Saturation pressure [Pa] at T [K] (0 above the critical temperature)."""
        if T >= 647.0:
            return 0.0
        lnp = np.interp(T, self.fluid._sat["T"], self.fluid._sat_lnp)
        return float(np.exp(lnp))

    # ======================================================================
    # One time step
    # ======================================================================
    def _step(self, dt, period):
        c = self.case
        f = self.fluid
        Nw = self.Nw
        m = period.mdot
        rev = self.reversed
        self._T_in_current = period.T_in
        gw = 1.0 if c.gravity_work else 0.0
        h_a_old, h_t_old = self.h_a.copy(), self.h_t.copy()
        T_s_old = self.T_s.copy() if self.grid else None
        # time-discretisation coefficients: implicit Euler, or variable-step BDF2
        prev = self._hist
        w_ = dt / prev["dt"] if prev is not None else np.inf
        if c.time_scheme == "bdf2" and prev is not None and w_ <= 2.4:   # zero-stable for ratio < 1 + sqrt(2)
            a0, a1, a2 = (1 + 2 * w_) / (1 + w_), 1 + w_, w_ * w_ / (1 + w_)
            h_a_nm1, h_t_nm1, T_s_nm1 = prev["h_a"], prev["h_t"], prev["T_s"]
        else:
            a0, a1, a2 = 1.0, 1.0, 0.0
            h_a_nm1, h_t_nm1, T_s_nm1 = h_a_old, h_t_old, T_s_old
        h_a, h_t = self.h_a.copy(), self.h_t.copy()
        p_a, p_t = self.p_a.copy(), self.p_t.copy()
        T_s = T_s_old.copy() if self.grid else None
        ds = self.ds[:Nw]
        dtvd = self.dtvd[:Nw]
        flash_md = None
        it_used = 0
        converged = False
        self._march_ok = True
        # Unknown layout: columns 0..Nw-1 hold the feed (down-flowing) stream's nodes
        # 1..Nw (node 0 is the known inlet), columns Nw..2Nw-1 the riser's nodes
        # 0..Nw-1; the riser's node Nw is the feed's node Nw (the bottom of the loop).
        j = np.arange(Nw)
        colF_j, colF_j1 = j - 1, j
        colR_j, colR_j1 = Nw + j, np.where(j + 1 < Nw, Nw + j + 1, Nw - 1)
        if rev:        # feed = tubing, riser = annulus
            colA_j, colA_j1, colT_j, colT_j1 = colR_j, colR_j1, colF_j, colF_j1
            rowA, rowT = Nw + j, j
        else:          # feed = annulus, riser = tubing
            colA_j, colA_j1, colT_j, colT_j1 = colF_j, colF_j1, colR_j, colR_j1
            rowA, rowT = j, Nw + j
        sgn_a = -1.0 if rev else 1.0        # +1: annulus flows down (advection +m h_{j+1} - m h_j)
        sgn_t = -sgn_a
        dh_res = 0.0                        # enthalpy added by a reservoir element at TD
        Gaw_t = np.zeros(Nw)                # riser <-> wall conductance (U-tube only)
        for it in range(c.max_iter):
            it_used = it + 1
            if rev:
                h_in = float(f.h_pT(p_t[0], period.T_in)[0])
                h_t[0] = h_in
                h_a[Nw] = h_t[Nw]
                nt = f.ph(p_t, h_t)
                na = f.ph(np.concatenate([p_a[:Nw], [p_t[Nw]]]), h_a)
            else:
                h_in = float(f.h_pT(p_a[0], period.T_in)[0])
                h_a[0] = h_in
                na = f.ph(p_a, h_a)
                if c.bottom_element is not None:
                    # reservoir element in series at TD: the riser starts at the element's outlet state
                    T_res_out = float(c.bottom_element.T_out(float(na.T[Nw]), self.t + dt))
                    dh_res = float(f.h_pT(p_a[Nw], T_res_out)[0]) - h_a[Nw]
                h_t[Nw] = h_a[Nw] + dh_res
                nt = f.ph(np.concatenate([p_t[:Nw], [p_a[Nw]]]), h_t)
            Da = na.dTdh.copy()
            ca = na.T - Da * h_a
            Dt = nt.dTdh.copy()
            ct = nt.T - Dt * h_t
            if rev:
                Da[Nw], ca[Nw] = Dt[Nw], ct[Nw]
            elif c.bottom_element is None:
                Dt[Nw], ct[Nw] = Da[Nw], ca[Nw]
            # film coefficients ------------------------------------------------
            if c.h_wall_override is not None:
                h_ao = h_ai = h_ti = np.full(Nw, float(c.h_wall_override))
            else:
                h_ao, h_ai = self._node_htc(na, m, "annulus")
                h_ti = self._node_htc(nt, m, "tubing")
            if self.utube:
                u = c.utube
                R_p = 0.5 * (pipe_resistance(u["r_pipe_in"], u["r_pipe_out"], u["k_pipe"], h_ao)
                             + pipe_resistance(u["r_pipe_in"], u["r_pipe_out"], u["k_pipe"], h_ti))
                R_lw, R_ll, _, _ = utube_delta(self.r_ci, u["r_pipe_out"], self.x_c, u["k_grout"], self.k_wall_rock, R_p)
                R_far = self.Rhalf if self.grid else c.R_formation
                Gaw = ds / (R_lw + R_far)
                Gaw_t = ds / (R_lw + R_far)
                Gat = ds / R_ll
                R_film_a = R_lw
            else:
                Gat = ds / (1.0 / (h_ai * 2 * math.pi * self.r_to) + self.Rtw + 1.0 / (h_ti * 2 * math.pi * self.r_ti))
                if self.grid:
                    Gaw = ds / (1.0 / (h_ao * 2 * math.pi * self.r_ci) + self.Rhalf)
                else:
                    Gaw = ds / (1.0 / (h_ao * 2 * math.pi * self.r_ci) + c.R_formation)
                R_film_a = 1.0 / (h_ao * 2 * math.pi * self.r_ci)
            # capacities
            rho_a = 0.5 * (na.rho[:-1] + na.rho[1:])
            rho_t = 0.5 * (nt.rho[:-1] + nt.rho[1:])
            if c.fluid_storage:
                Ca = rho_a * self.A_a * ds / dt
                Ct = rho_t * self.A_t * ds / dt
            else:
                Ca = np.zeros(Nw)
                Ct = np.zeros(Nw)
            # space weights -----------------------------------------------------
            # Box scheme (theta = 0.5) unless a cell's transfer number lambda =
            # G/(m cp) exceeds 1, then exponential weighting (exact for a fixed
            # external temperature, monotone).  For the wall the conductance is
            # taken along the path to the rock (film + completion + transient
            # rock resistance), because the thin casing node itself is not a
            # stiff temperature reservoir.  Frozen after the first iteration.
            # The weight applies to each stream's downstream node.
            if it == 0:
                if c.scheme == "upwind" or m <= 0.0:
                    tha = np.ones(Nw)
                    tht = np.ones(Nw)
                else:
                    Dac = 0.5 * (Da[:-1] + Da[1:])
                    Dtc = 0.5 * (Dt[:-1] + Dt[1:])
                    R_film = R_film_a
                    if self.grid:
                        tD = np.maximum(self.alpha_wall * max(self.t + dt, dt) / self.r_hole ** 2, 1e-6)
                        f_t = np.where(tD <= 1.5, 1.1281 * np.sqrt(tD) * (1 - 0.3 * np.sqrt(tD)),
                                       (0.4063 + 0.5 * np.log(tD)) * (1 + 0.6 / tD))
                        R_path = R_film + self.R_comp + f_t / (2 * math.pi * self.k_wall_rock)
                    else:
                        R_path = R_film + c.R_formation
                    la = (ds / R_path + Gat) * Dac / m
                    lt = (Gat + Gaw_t) * Dtc / m

                    def wexp(lam):
                        lam = np.maximum(lam, 1e-12)
                        return np.where(lam <= 1.0, 0.5, 1.0 / (1.0 - np.exp(-np.minimum(lam, 700))) - 1.0 / lam)
                    tha = wexp(la)
                    tht = wexp(lt)
                # weights on (node j, node j+1): downstream node is j+1 for a down-flowing stream
                wa_j, wa_j1 = (tha, 1 - tha) if rev else (1 - tha, tha)
                wt_j, wt_j1 = (1 - tht, tht) if rev else (tht, 1 - tht)
            Hold_a = wa_j * h_a_old[:-1] + wa_j1 * h_a_old[1:]
            Hold_t = wt_j * h_t_old[:-1] + wt_j1 * h_t_old[1:]
            Hnm1_a = wa_j * h_a_nm1[:-1] + wa_j1 * h_a_nm1[1:]
            Hnm1_t = wt_j * h_t_nm1[:-1] + wt_j1 * h_t_nm1[1:]
            Sa = a1 * Hold_a - a2 * Hnm1_a          # history part of the storage term
            St = a1 * Hold_t - a2 * Hnm1_t
            # assemble ------------------------------------------------------------
            rows, cols, vals = [], [], []
            rhs = np.zeros(2 * Nw + (self.Ns if self.grid else 0))
            Gsum = Gaw + Gat
            Ta_bar_c = wa_j * ca[:-1] + wa_j1 * ca[1:]
            Tt_bar_c = wt_j * ct[:-1] + wt_j1 * ct[1:]
            dKE_a = m * self.dKE if rev else 0.0
            dKE_t = 0.0 if rev else m * self.dKE
            # annulus rows: storage + advection = rock + tubing exchange + gravity work
            rows += [rowA, rowA, rowA, rowA]
            cols += [colA_j1, colA_j, colT_j, colT_j1]
            vals += [a0 * Ca * wa_j1 + sgn_a * m + Gsum * wa_j1 * Da[1:],
                     a0 * Ca * wa_j - sgn_a * m + Gsum * wa_j * Da[:-1],
                     -Gat * wt_j * Dt[:-1],
                     -Gat * wt_j1 * Dt[1:]]
            rhs[rowA] += Ca * Sa + sgn_a * m * G_ACC * dtvd * gw - dKE_a - Gsum * Ta_bar_c + Gat * Tt_bar_c
            # tubing rows
            rows += [rowT, rowT, rowT, rowT]
            cols += [colT_j1, colT_j, colA_j1, colA_j]
            vals += [a0 * Ct * wt_j1 + sgn_t * m + Gat * wt_j1 * Dt[1:],
                     a0 * Ct * wt_j - sgn_t * m + Gat * wt_j * Dt[:-1],
                     -Gat * wa_j1 * Da[1:],
                     -Gat * wa_j * Da[:-1]]
            rhs[rowT] += Ct * St + sgn_t * m * G_ACC * dtvd * gw - dKE_t + Gat * Ta_bar_c - Gat * Tt_bar_c
            if self.utube:
                rows += [rowT, rowT]
                cols += [colT_j1, colT_j]
                vals += [Gaw_t * wt_j1 * Dt[1:], Gaw_t * wt_j * Dt[:-1]]
                rhs[rowT] += -Gaw_t * Tt_bar_c
            if dh_res != 0.0:
                # the riser's bottom node is h_a[Nw] + dh_res: move the known part to the RHS
                cT = (a0 * Ct * wt_j1 + sgn_t * m + Gat * wt_j1 * Dt[1:])[Nw - 1]
                cA = (-Gat * wt_j1 * Dt[1:])[Nw - 1]
                rhs[rowT[Nw - 1]] -= cT * dh_res
                rhs[rowA[Nw - 1]] -= cA * dh_res
            self.dKE_used = self.dKE.copy()
            if self.grid:
                off = 2 * Nw
                sw = off + self.sid_wall
                # annulus <- wall
                rows.append(rowA); cols.append(sw); vals.append(-Gaw)
                # wall rows: extra diagonal + coupling to annulus nodes
                rows.append(sw); cols.append(sw); vals.append(Gaw)
                rows += [sw, sw]
                cols += [colA_j1, colA_j]
                vals += [-Gaw * wa_j1 * Da[1:], -Gaw * wa_j * Da[:-1]]
                rhs[sw] += Gaw * Ta_bar_c
                if self.utube:
                    rows.append(rowT); cols.append(sw); vals.append(-Gaw_t)
                    rows.append(sw); cols.append(sw); vals.append(Gaw_t)
                    rows += [sw, sw]
                    cols += [colT_j1, colT_j]
                    vals += [-Gaw_t * wt_j1 * Dt[1:], -Gaw_t * wt_j * Dt[:-1]]
                    rhs[sw] += Gaw_t * Tt_bar_c
                # bottom coupling (rock below TD <-> bottom node, which is the feed's node Nw, row Nw-1)
                if self.bot_sid.size and c.bottom_element is None:
                    Gb = 1.0 / (self.bot_R + 1.0 / (h_ao[-1] * self.bot_A))
                    sb = off + self.bot_sid
                    rows += [np.array([Nw - 1]), sb, sb, np.full(sb.size, Nw - 1)]
                    cols += [np.array([Nw - 1]), sb, np.full(sb.size, Nw - 1), sb]
                    vals += [np.array([Gb.sum() * Da[Nw]]), Gb, -Gb * Da[Nw], -Gb]
                    rhs[Nw - 1] += -Gb.sum() * ca[Nw]
                    rhs[sb] += Gb * ca[Nw]
                else:
                    Gb = np.zeros(0)
                # solid block
                rows.append(off + self._Ks_rows); cols.append(off + self._Ks_cols); vals.append(self._Ks_vals)
                rows.append(off + np.arange(self.Ns)); cols.append(off + np.arange(self.Ns)); vals.append(a0 * self.C_s / dt)
                rhs[off:] += self.C_s / dt * (a1 * T_s_old - a2 * T_s_nm1) + self.rhs_eq
            else:
                rhs[rowA] += Gaw * self.T0_c[:Nw]
                if self.utube:
                    rhs[rowT] += Gaw_t * self.T0_c[:Nw]
            rows = np.concatenate([np.atleast_1d(x) for x in rows])
            cols = np.concatenate([np.atleast_1d(x) for x in cols])
            vals = np.concatenate([np.atleast_1d(x) for x in vals])
            keep = cols >= 0
            # known inlet node (column -1): move to RHS
            known = cols == -1
            np.add.at(rhs, rows[known], -vals[known] * h_in)
            rows, cols, vals = rows[keep], cols[keep], vals[keep]
            nunk = rhs.size
            A = sp.csc_matrix((vals, (rows, cols)), shape=(nunk, nunk))
            x = spla.spsolve(A, rhs)
            if not np.all(np.isfinite(x)):
                raise RuntimeError("linear solve failed (non-finite solution)")
            h_a_new = h_a.copy()
            h_t_new = h_t.copy()
            if rev:
                h_t_new[1:] = x[:Nw]
                h_a_new[:Nw] = x[Nw:2 * Nw]
                h_a_new[Nw] = h_t_new[Nw]
            else:
                h_a_new[1:] = x[:Nw]
                h_t_new[:Nw] = x[Nw:2 * Nw]
                h_t_new[Nw] = h_a_new[Nw] + dh_res
            dh = max(np.max(np.abs(h_a_new - h_a)), np.max(np.abs(h_t_new - h_t)))
            dT = 0.0
            if self.grid:
                T_s_new = x[2 * Nw:]
                dT = float(np.max(np.abs(T_s_new - T_s)))
                T_s = T_s_new
            h_a, h_t = h_a_new, h_t_new
            p_a_new, p_t_new, dKE, flash_md = self._march(h_a, h_t, m)
            dp = max(np.max(np.abs(p_a_new - p_a)), np.max(np.abs(p_t_new - p_t)))
            p_a, p_t = p_a_new, p_t_new
            self.p_a, self.p_t = p_a, p_t       # used as guesses by _march
            self.dKE = dKE
            converged = dh < c.tol_h and dT < c.tol_T and dp < c.tol_p and it >= 1 and self._march_ok
            if converged:
                break
        # commit ------------------------------------------------------------------
        self.h_a, self.h_t, self.p_a, self.p_t = h_a, h_t, p_a, p_t
        if self.grid:
            self.T_s = T_s
        na = f.ph(p_a, h_a)
        nt = f.ph(p_t, h_t)
        # heat flows from the temperatures the energy equations were solved with
        # (last linearisation), so the discrete balance closes exactly
        Ta_lin = Da * h_a + ca
        Tt_lin = Dt * h_t + ct
        Ta_bar = wa_j * Ta_lin[:-1] + wa_j1 * Ta_lin[1:]
        Tt_bar = wt_j * Tt_lin[:-1] + wt_j1 * Tt_lin[1:]
        if self.grid:
            Tw = T_s[self.sid_wall]
            q_rock = Gaw * (Tw - Ta_bar) + Gaw_t * (Tw - Tt_bar)
            Q_bot = float(np.sum(Gb * (T_s[self.bot_sid] - Ta_lin[Nw]))) if (self.bot_sid.size and c.bottom_element is None) else 0.0
        else:
            Tw = self.T0_c[:Nw]
            q_rock = Gaw * (Tw - Ta_bar) + Gaw_t * (Tw - Tt_bar)
            Q_bot = 0.0
        q_at = Gat * (Tt_bar - Ta_bar)                      # tubing -> annulus
        q_riser = -q_at if rev else q_at                    # riser -> feed (heat short-circuited back down)
        dE = float(np.sum(Ca * (a0 * (wa_j * h_a[:-1] + wa_j1 * h_a[1:]) - Sa)) +
                   np.sum(Ct * (a0 * (wt_j * h_t[:-1] + wt_j1 * h_t[1:]) - St)))
        self._hist = dict(h_a=h_a_old, h_t=h_t_old, T_s=T_s_old, dt=dt)
        self._last = dict(na=na, nt=nt, Ta_bar=Ta_bar, Tt_bar=Tt_bar, Tw=Tw, q_rock=q_rock / ds,
                          q_at=q_riser / ds, h_ao=h_ao, h_ai=h_ai, h_ti=h_ti, Gaw=Gaw, Gat=Gat)
        nR, nF = (na, nt) if rev else (nt, na)
        return dict(h_in=h_in, Q_rock=float(np.sum(q_rock)) + Q_bot, Q_at=float(np.sum(q_riser)),
                    dE=dE, Q_ke=float(m * np.sum(self.dKE_used)), iters=it_used, flash_md=flash_md,
                    T_out=float(nR.T[0]), x_out=float(nR.x[0]), T_in_actual=float(nF.T[0]), converged=converged,
                    Q_res=float(m * dh_res))

    # ======================================================================
    def run(self, progress: bool = False) -> Results:
        c = self.case
        res = self.res
        outs = sorted(set(c.output_times))
        t = 0.0
        wall0 = _time.time()
        for period in c.schedule:
            t_end_p = t + period.duration
            dt = c.dt_initial
            while t < t_end_p - 1e-6:
                step = min(dt, t_end_p - t)
                nxt = next((o for o in outs if o > t + 1e-6), None)
                if nxt is not None and nxt < t + step - 1e-6:
                    step = nxt - t
                # avoid a tiny remainder at the end of a period
                if 0 < t_end_p - (t + step) < 0.2 * step:
                    step = t_end_p - t
                saved = self._save_state()
                stepper = self._step_single if self.single else self._step
                info = stepper(step, period)
                cuts = 0
                while not info["converged"] and cuts < 6 and step > 2.0:
                    self._restore_state(saved)
                    step *= 0.5
                    cuts += 1
                    info = stepper(step, period)
                if not info["converged"]:
                    res.meta["unconverged_steps"] = res.meta.get("unconverged_steps", 0) + 1
                t += step
                self.t = t
                m = period.mdot
                gw = 1.0 if c.gravity_work else 0.0
                if self.single:
                    Nw = self.Nw
                    h_out = self.h_a[0] if self.up else self.h_a[Nw]
                    p_in, p_out = (self.p_a[Nw], self.p_a[0]) if self.up else (self.p_a[0], self.p_a[Nw])
                    # potential energy taken from (production) or given to (injection) the stream
                    Q_pe = m * G_ACC * gw * (self.tvd_f[Nw] - self.tvd_f[0]) * (1.0 if self.up else -1.0)
                else:
                    h_out = self.h_a[0] if self.reversed else self.h_t[0]
                    p_in, p_out = (self.p_t[0], self.p_a[0]) if self.reversed else (self.p_a[0], self.p_t[0])
                    Q_pe = 0.0
                res.append(t=t, mdot=m, T_in=info["T_in_actual"], T_out=info["T_out"],
                           p_in=p_in, p_out=p_out, h_in=info["h_in"], h_out=h_out,
                           x_out=info["x_out"], Q_wellhead=m * (h_out - info["h_in"]),
                           Q_rock=info["Q_rock"], Q_riser_loss=info["Q_at"], dE_fluid=info["dE"],
                           iters=info["iters"], flash_md=(np.nan if info["flash_md"] is None else info["flash_md"]),
                           choke_dp=getattr(self, "choke_dp", 0.0), Q_ke=info["Q_ke"],
                           choke_warning=float(getattr(self, "choke_warning", False)),
                           Q_res=info.get("Q_res", 0.0), pump_dp=getattr(self, "pump_dp", 0.0),
                           pump_kW=getattr(self, "pump_kW", 0.0), inflow_margin=getattr(self, "inflow_margin", np.nan),
                           Q_pe=Q_pe)
                if nxt is not None and abs(t - nxt) < 1e-3:
                    res.profiles.append(self.profile())
                if progress:
                    print(f"  t={t / DAY:10.3f} d  dt={step:9.1f} s  it={info['iters']:2d}  "
                          f"Tout={info['T_out'] - 273.15:7.2f} C  Q={m * (h_out - info['h_in']) / 1e3:9.1f} kW")
                dt = min(step * c.dt_growth, c.dt_max) if step >= dt * 0.999 else dt
        res.meta["wall_time_s"] = _time.time() - wall0
        return res.finalise()

    def profile(self):
        L = self._last
        Nw = self.Nw
        prof = dict(t=self.t, md=self.s_f[:Nw + 1].copy(), tvd=self.tvd_f[:Nw + 1].copy(),
                    T_a=L["na"].T.copy(), T_t=L["nt"].T.copy(), p_a=self.p_a.copy(), p_t=self.p_t.copy(),
                    x_t=L["nt"].x.copy(), x_a=L["na"].x.copy(), T0=self.T0_f[:Nw + 1].copy(),
                    md_c=self.s_c[:Nw].copy(), tvd_c=self.tvd_c[:Nw].copy(), T_wall=np.asarray(L["Tw"]).copy(),
                    q_rock=L["q_rock"].copy(), q_riser=L["q_at"].copy())
        if self.grid:
            prof["T_s"] = self.T_s.copy()
        return prof
