"""Verification of the uneven-flow fracture element (aharrsim.reservoir.ZonedFractures).

  A  even flow (Lorenz coefficient 0) returns the parallel-fracture solution
  B  closely spaced fractures: the closed form theta = 1 - Phi((ln x + sigma^2/2) / sigma)
  C  heat balance: the integral of the retained share over the flushed share is 1
  D  two zones, one without flow: the swept-fraction model (breadth phi B)
  E  convergence with the number of zones
  F  an independent finite-volume solution of conduction between fractures that carry unequal flows, for zones of
     decreasing width (the zoned element assumes zones that do not exchange heat)
  G  the general element (UnevenFractures), which carries that conduction: equal flows, wide zones, finite-volume solution

    python zoned_verify.py      -> results/ZF_verify.json
"""
from __future__ import annotations

import sys
import time

import numpy as np
from scipy.linalg import solve_banded

from common import save
from aharrsim.reservoir import ParallelFractures, ZonedFractures, UnevenFractures, zoned_dense_limit, lorenz_sigma, lognormal_zones, lorenz_of

YEAR = 365.25 * 86400.0
BASE = dict(Q=100.0 / 900.0, N=300, L_flow=137.0, B=143.0, S=1524.0 / 150, k_r=2.8, rhoc_r=2650 * 850.0, T_rock_C=200.0, T_inj_C=60.0, rhoc_w=4160.0 * 900.0)


def fv_fractures(flows, S, L, B, k, rhoc_r, rhoc_w, times, nx_per=16, nz=160, nt_per_decade=120):
    """Finite-volume solution for fractures at spacing S carrying the volume flows `flows` [m3/s] (any order), with one-dimensional
    conduction across the rock between them (as in Gringarten et al.), insulated outer faces, cold water entering at z = 0.
    Second-order upwind marching along the flow, implicit in time.  Returns the flow-weighted outlet T_D at `times`."""
    flows = np.asarray(flows, float); nf = len(flows)
    nx = nf * nx_per; dx = S / nx_per; dz = L / nz
    frac = (np.arange(nf) * nx_per + nx_per // 2).astype(int)              # cell that holds each fracture
    cap = np.zeros(nx); cap[frac] = rhoc_w * flows / B / dz                # advective conductance per unit dz [W/m2/K]
    T = np.zeros((nz + 1, nx))                                             # T_D: 0 = rock temperature, inlet layer handled as 1 in the fractures
    t_edges = np.concatenate([[0.0], np.logspace(np.log10(times.min() * 1e-3), np.log10(times.max()), int(nt_per_decade * (3 + np.log10(times.max() / times.min()))))])
    out = []; rec_t = []
    ab = np.zeros((3, nx))
    cond = k / dx
    for it in range(1, len(t_edges)):
        dt = t_edges[it] - t_edges[it - 1]
        st = rhoc_r * dx / dt
        Tn = T.copy()
        for n in range(1, nz + 1):
            up1 = np.where(cap > 0, 1.0, 0.0) if n == 1 else Tn[n - 1]     # inlet: the fractures carry water at T_D = 1
            if n == 1:
                a_d, rhs_adv = cap, cap * up1                              # first layer: first-order upwind
            else:
                up2 = np.where(cap > 0, 1.0, 0.0) if n == 2 else Tn[n - 2]
                a_d, rhs_adv = 1.5 * cap, cap * (2.0 * up1 - 0.5 * up2)    # second-order upwind
            ab[1] = st + 2.0 * cond + a_d
            ab[1, 0] -= cond; ab[1, -1] -= cond                             # insulated outer faces
            ab[0, 1:] = -cond; ab[2, :-1] = -cond
            Tn[n] = solve_banded((1, 1), ab, st * T[n] + rhs_adv)
        T = Tn
        out.append(float((flows * T[nz, frac]).sum() / flows.sum())); rec_t.append(t_edges[it])
    return np.interp(times, rec_t, out)


def main():
    res = {}
    yrs = np.array([0.5, 1, 2, 3, 5, 7, 10, 15, 20, 25, 30], float); t = yrs * YEAR
    pf = ParallelFractures(**BASE)
    # A
    z0 = ZonedFractures(**BASE, lorenz=0.0)
    res["A_even_flow_max_diff"] = float(np.abs(z0.T_D(t) - pf.T_D(t)).max())
    # B  dense fractures (1 m apart), three Lorenz coefficients
    dense = dict(BASE); dense["N"] = 3048; dense["S"] = 0.5
    B_ = {}
    for Lc in (0.2, 0.5, 0.8):
        z = ZonedFractures(**dense, lorenz=Lc, n_zones=400)
        x = np.array([0.01, 0.02, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 4.0])
        num = 1.0 - z.T_D(x * z.t_sweep); ref = zoned_dense_limit(x, Lc)
        B_[str(Lc)] = dict(x=x.tolist(), element=num.tolist(), closed_form=ref.tolist(), max_diff=float(np.abs(num - ref).max()))
    res["B_dense_limit"] = B_
    # C  heat balance at the base spacing
    C_ = {}
    for Lc in (0.0, 0.2, 0.5, 0.8):
        z = ZonedFractures(**BASE, lorenz=Lc, n_zones=200)
        x = np.logspace(-5, 5, 1201)
        th = 1.0 - z.T_D(x * z.t_sweep)
        C_[str(Lc)] = float(np.trapezoid(th * x, np.log(x)))
    res["C_heat_balance"] = C_
    # D  binary: half of the fracture area takes all the flow
    D_ = {}
    for phi in (0.5, 0.25):
        z = ZonedFractures(**BASE, flow=[1.0, 0.0], area=[phi, 1.0 - phi])
        b = dict(BASE); b["B"] = phi * BASE["B"]
        D_[str(phi)] = dict(max_diff=float(np.abs(z.T_D(t) - ParallelFractures(**b).T_D(t)).max()), lorenz=z.lorenz)
    res["D_binary"] = D_
    # E  convergence with the number of zones
    E_ = {}
    for Lc in (0.2, 0.5, 0.8):
        ref = ZonedFractures(**BASE, lorenz=Lc, n_zones=2000).T_D(t)
        E_[str(Lc)] = {str(K): float(np.abs(ZonedFractures(**BASE, lorenz=Lc, n_zones=K).T_D(t) - ref).max()) for K in (10, 20, 40, 60, 120, 240)}
    res["E_zones"] = E_
    # F  finite-volume cross-check
    F_ = {}
    S, L, B_m = BASE["S"], BASE["L_flow"], BASE["B"]; qf = BASE["Q"] / BASE["N"]
    args = (S, L, B_m, BASE["k_r"], BASE["rhoc_r"], BASE["rhoc_w"], t)
    t0 = time.time()
    num = fv_fractures([qf] * 2, *args)
    F_["even_flow_fv_minus_analytical"] = dict(yr=yrs.tolist(), fv=num.tolist(), analytical=pf.T_D(t).tolist(), max_diff=float(np.abs(num - pf.T_D(t)).max()))
    for tag, hi, lo in (("80_20", 1.6, 0.4), ("95_5", 1.9, 0.1)):
        zoned = ZonedFractures(**BASE, flow=[hi, lo], area=[1.0, 1.0])
        ref = zoned.T_D(t)
        rows = []
        for n in (1, 2, 4, 8, 16):                                          # zone width 2 n S (mirror symmetry at the outer faces)
            fl = [hi * qf] * n + [lo * qf] * n
            num = fv_fractures(fl, *args)
            rows.append(dict(zone_width_m=2 * n * S, fv=num.tolist(), max_diff_zoned=float(np.abs(num - ref).max()), max_diff_even=float(np.abs(num - pf.T_D(t)).max())))
            print(tag, 2 * n * S, rows[-1]["max_diff_zoned"], rows[-1]["max_diff_even"], round(time.time() - t0), flush=True)
        F_[tag] = dict(lorenz=zoned.lorenz, zoned=ref.tolist(), even=pf.T_D(t).tolist(), yr=yrs.tolist(), rows=rows)
    res["F_finite_volume"] = F_
    # G  the general element (UnevenFractures: unequal flows with conduction between the fractures) against its limits and the finite-volume solution
    G_ = dict(even_flow_max_diff=float(np.abs(UnevenFractures(**BASE, mult=np.ones(6)).T_D(t) - pf.T_D(t)).max()), cases=[])
    for tag, hi, lo in (("80_20", 1.6, 0.4), ("95_5", 1.9, 0.1)):
        zref = ZonedFractures(**BASE, flow=[hi, lo], area=[1.0, 1.0]).T_D(t)
        for n in (1, 2, 4, 8, 16, 64):
            a = UnevenFractures(**BASE, mult=[hi] * n + [lo] * n).T_D(t)
            row = dict(split=tag, zone_width_m=2 * n * S, T_D=a.tolist(), max_diff_zoned=float(np.abs(a - zref).max()), max_diff_even=float(np.abs(a - pf.T_D(t)).max()))
            fv = next((r_["fv"] for r_ in F_[tag]["rows"] if abs(r_["zone_width_m"] - 2 * n * S) < 1e-6), None)
            if fv is not None:
                row["max_diff_finite_volume"] = float(np.abs(np.array(fv) - a).max())
            G_["cases"].append(row)
            print("G", tag, round(2 * n * S), {k: round(v, 4) for k, v in row.items() if k.startswith("max")}, flush=True)
    # finer finite-volume grid for one case: the difference is the discretisation error of the finite-volume solution
    a = UnevenFractures(**BASE, mult=[1.6] * 4 + [0.4] * 4).T_D(t)
    fine = fv_fractures([1.6 * qf] * 4 + [0.4 * qf] * 4, *args, nx_per=24, nz=320, nt_per_decade=240)
    G_["finite_volume_fine_grid_max_diff"] = float(np.abs(fine - a).max())
    print("G fine grid", G_["finite_volume_fine_grid_max_diff"], flush=True)
    res["G_uneven"] = G_
    save("ZF_verify", res)
    for k_, v in res.items():
        if k_ == "G_uneven":
            continue
        if k_ == "F_finite_volume":
            print("F even-flow check", v["even_flow_fv_minus_analytical"]["max_diff"]); continue
        print(k_, v if not isinstance(v, dict) or k_ != "B_dense_limit" else {a: b["max_diff"] for a, b in v.items()})


if __name__ == "__main__":
    main()
