"""Calculations added after the editorial review of the manuscript (November 2026). Screening-level unless stated.

    python study_review.py gross      -> results/RV1_gross_need.json       flow and swept volume for gross targets, beside the net targets
    python study_review.py share      -> results/RV2_impedance_share.json  near-well share of the impedance, a = 0.25, 0.5, 0.75
    python study_review.py seeds      -> results/RV3_arrangements.json     20 arrangements of the unequal flows in place of 5
    python study_review.py red        -> results/RV4_red_flow_basis.json   Project Red fit for the produced, mean and injected flow
    python study_review.py outside    -> results/RV5_outside_conduction.json  conduction into the block from the rock around it
    python study_review.py boot       -> results/RV6_percentiles.json      bootstrap of the Latin-hypercube percentiles of the closed loops
"""
from __future__ import annotations

import math
import sys
import time

import numpy as np

from common import save, load
from study_open import SYS, design, strip, geometry, YRS, cbar

YEAR = 365.25 * 86400.0

# figures published for the four fractured examples, and the basis on which each was stated (Table C1 of the manuscript)
TARGETS = {
    10: [dict(value=3.0, basis="not stated", note="lower end of the 3-5 MWe of the classification")],
    14: [dict(value=3.5, basis="gross", note="Project Red, 30-day test, peak (Norbeck and Latimer, 2023)")],
    141: [dict(value=10.0, basis="gross", note="Cape Station, 30-day test (Fervo Energy, 2024)"),
          dict(value=16.0, basis="gross", note="Cape Station, design figure per producer (Fervo Energy, 2026a)")],
    19: [dict(value=50.0, basis="not stated", note="classification")],
}


def need(sid, target, basis="net", a_nw=0.5, hold=0.9, spacing=10.0, m_grid=None):
    """Smallest flow whose 30-year mean output (net or gross) reaches `target` when the swept volume is just large enough for the
    output of year 30 to stay within (1 - hold) of year 1; fractures `spacing` apart. The rule of study_open.need_fine, for either basis."""
    s = SYS[sid]
    k1, k30, km = ("net1", "net30", "mean_net") if basis == "net" else ("gross1", "gross30", "mean_gross")
    for m in (np.arange(20.0, 400.1, 2.5) if m_grid is None else m_grid):
        def held(vf):
            Wv, Lv = geometry(sid, vf)
            d = design(s, N=max(int(Wv / spacing), 2), m=float(m), W=Wv, L=Lv, a_nw=a_nw)
            return d, (d[k30] >= hold * d[k1] and d[k1] > 0)
        d_hi, ok = held(256.0)
        if not ok or d_hi[km] < target:
            continue
        lo, hi = 0.25, 256.0
        for _ in range(40):
            mid = math.sqrt(lo * hi)
            if held(mid)[1]:
                hi = mid
            else:
                lo = mid
        d, _ = held(hi)
        if d[km] >= target:
            return dict(strip(d), vf=hi, basis=basis, target=target, a_nw=a_nw)
    return None


def brief(d):
    return d and dict(m=d["m"], V_km3=round(d["V_km3"], 3), vf=round(d["vf"], 2), N=d["N"], W=round(d["W"]), L=round(d["L"], 1), gross1=round(d["gross1"], 2), gross30=round(d["gross30"], 2),
                      mean_gross=round(d["mean_gross"], 2), mean_net=round(d["mean_net"], 2), net1=round(d["net1"], 2), W_circ=round(d["W_circ"], 2), x_30=round(d["x_30"], 2), y_30=round(d["y_30"], 2))


def gross():
    out = dict(rule="smallest flow whose 30-year mean reaches the target, with the volume just large enough for year 30 to stay within 10 % of year 1; fractures 10 m apart",
               systems={})
    for sid, tg in TARGETS.items():
        rows = []
        for t in tg:
            n = need(sid, t["value"], "net"); g = need(sid, t["value"], "gross")
            rows.append(dict(target=t, net_basis=n, gross_basis=g))
            print(sid, t["value"], "\n   net  ", brief(n), "\n   gross", brief(g), flush=True)
        out["systems"][str(sid)] = dict(name=SYS[sid]["name"], example_V_km3=SYS[sid]["W"] * SYS[sid]["L"] * SYS[sid]["B"] * SYS[sid].get("wings", 1) / 1e9, rows=rows)
        save("RV1_gross_need", out)


def share():
    """The near-well share a of the impedance (Eq. 9) sets how far the impedance falls when fractures are added: I -> a I0 for many fractures.
    Net-target designs and the dense version of each example for a = 0.25, 0.5 and 0.75."""
    out = dict(a=[0.25, 0.5, 0.75], systems={})
    for sid, tg in TARGETS.items():
        s = SYS[sid]; o = dict(name=s["name"], target=tg[0]["value"], rows=[])
        for a in out["a"]:
            n = need(sid, tg[0]["value"], "net", a_nw=a)
            dense = design(s, N=int(s["W"] / 10.0), a_nw=a)
            o["rows"].append(dict(a=a, need=n, dense_example=strip(dense)))
            print(sid, "a", a, "need", brief(n), "| dense example: net", round(dense["mean_net"], 2), "W_circ", round(dense["W_circ"], 2), "I", round(dense["I_MPa_per_Ls"], 3), flush=True)
        out["systems"][str(sid)] = o
        save("RV2_impedance_share", out)


def seeds():
    """Unequal flows are placed along the well at random; the scenario results are means over the arrangements. 5 arrangements (as reported)
    against 20, for the designs that meet the net targets with even flow."""
    from study_rock import SCEN, need_designs
    out = dict(scenarios=SCEN, systems={})
    for sid, nd in need_designs().items():
        s = SYS[sid]; kw = dict(N=nd["N"], m=nd["m"], W=nd["W"], L=nd["L"])
        even = design(s, **kw)
        o = dict(design=nd, even=dict(mean_net=even["mean_net"], T30=even["T30"]), scen={})
        for sc, u in SCEN.items():
            t0 = time.time(); per = []
            for sd in range(20):
                d = design(s, uneven=dict(u, seeds=(sd,)), **kw)
                per.append(dict(seed=sd, mean_net=d["mean_net"], T30=d["T30"], net30=d["net30"]))
            mn = np.array([p["mean_net"] for p in per]); T30 = np.array([p["T30"] for p in per])
            # note: the reported figures average the temperature histories over the arrangements before the plant; this averages the outputs
            o["scen"][sc] = dict(per_seed=per, kept_5=float(mn[:5].mean() / even["mean_net"]), kept_20=float(mn.mean() / even["mean_net"]), kept_sd=float(mn.std(ddof=1) / even["mean_net"]),
                                 T30_5=float(T30[:5].mean()), T30_20=float(T30.mean()), T30_sd=float(T30.std(ddof=1)))
            print(sid, sc, {k: round(v, 3) for k, v in o["scen"][sc].items() if k != "per_seed"}, round(time.time() - t0), "s", flush=True)
        out["systems"][str(sid)] = o
        save("RV3_arrangements", out)


def red():
    """Project Red: the record is plotted against material-balance time and about 70 % of the injected water is recovered (Fervo Energy, 2026a),
    so the flow that swept the fractures lies between the produced 36 kg/s and the injected rate. The fit of field_cases.red() for three flows."""
    import field_cases as fc
    rec = load("FE_project_red_digitised")
    t = np.array(rec["t_yr"]); T = np.array(rec["T_prod_C"]); r = fc.RED
    ok = t > 0.42
    for s in r["shut_ins"]:
        ok &= ~((t >= s - 0.005) & (t < s + 0.07))
    plateau = float(np.median(T[(t > 0.45) & (t < 0.98) & ok])); dT = plateau - r["T_inj"]
    TD_obs = (plateau - T) / dT
    recov = 0.70
    flows = dict(produced=r["m"], mean=r["m"] * 0.5 * (1.0 + 1.0 / recov), injected=r["m"] / recov)
    fs = np.round(np.arange(0.15, 1.001, 0.01), 3); phis = np.round(np.arange(0.20, 1.201, 0.01), 3)
    out = dict(recovered_share=recov, plateau_C=plateau, flows_kg_s=flows, cases={})
    for name, m in flows.items():
        grid = np.array([[float(np.sqrt(np.mean((fc.red_model(t[ok], f, p, m=m)[0] - TD_obs[ok]) ** 2)) * dT) for p in phis] for f in fs])
        i, j = np.unravel_index(np.argmin(grid), grid.shape)
        acc = [(float(fs[a]), float(phis[b])) for a in range(len(fs)) for b in range(len(phis)) if grid[a, b] < 1.0]
        prod = np.array([f * p for f, p in acc])
        jf = int(np.argmin(np.abs(phis - 1.0))); i1 = int(np.argmin(grid[:, jf]))
        j2 = int(np.argmin(grid[-1]))
        ideal = float(grid[-1, jf])
        tt = np.linspace(0.05, 12.0, 240)
        TDi = fc.red_model(tt, 1.0, 1.0, m=m)[0]; TDf = fc.red_model(tt, float(fs[i1]), 1.0, m=m)[0]
        t10 = lambda TD: (float(tt[np.where(TD > 0.1)[0][0]]) if np.any(TD > 0.1) else None)
        out["cases"][name] = dict(m=m, best=dict(f=float(fs[i]), phi=float(phis[j]), rms_K=float(grid[i, j])), ideal_rms_K=ideal,
                                  fit_fractures=dict(f=float(fs[i1]), rms_K=float(grid[i1, jf]), t_10pct_yr=t10(TDf)), fit_volume=dict(phi=float(phis[j2]), rms_K=float(grid[-1, j2])),
                                  f_phi=[float(prod.min()), float(prod.max())] if len(acc) else None, n_acceptable=len(acc), ideal_t_10pct_yr=t10(TDi))
        print(name, round(m, 1), "kg/s:", {k: v for k, v in out["cases"][name].items() if k != "m"}, flush=True)
    save("RV4_red_flow_basis", out)


# ------------------------------------------------------------------------------------------------ conduction from the rock around the block
def block_2d(L, B, v_T, alpha, t_end_yr, dx=4.0, dz=4.0, up=160.0, down=80.0, out_z=160.0, dt_yr=0.02, outside=True, times=None):
    """Closely spaced fractures act as a zone in which the thermal front moves at v_T = L / t_s (water and rock at one temperature).
    Energy balance in the plane of the flow direction (xi) and the height (z), symmetric about mid-height:

        dT/dt + v_T dT/dxi [in the zone] = alpha (d2T/dxi2 + d2T/dz2)

    The zone is 0 < xi < L, |z| < B/2; water enters it at the injection temperature. With outside=False the rock around the zone does not
    conduct (the assumption of the screening model). Returns the dimensionless outlet temperature drop at `times` [yr] and the heat balance.
    Temperature is scaled: 0 = rock, 1 = injection."""
    import scipy.sparse as sp
    from scipy.sparse.linalg import splu
    nx_u, nx_z, nx_d = int(round(up / dx)), int(round(L / dx)), int(round(down / dx))
    nz_z, nz_o = int(round(0.5 * B / dz)), int(round(out_z / dz))
    dx = L / nx_z; dz = 0.5 * B / nz_z
    nx, nz = nx_u + nx_z + nx_d, nz_z + nz_o
    idx = lambda i, j: i * nz + j
    inzone = np.zeros((nx, nz), bool); inzone[nx_u:nx_u + nx_z, :nz_z] = True
    kx = np.ones((nx + 1, nz)); kz = np.ones((nx, nz + 1))                      # face conductances switched off where the rock does not conduct
    if not outside:
        for i in range(nx):
            for j in range(nz):
                if not inzone[i, j]:
                    kx[i, j] = kx[i + 1, j] = 0.0; kz[i, j] = kz[i, j + 1] = 0.0
    kx[0, :] = kx[nx, :] = 0.0; kz[:, 0] = 0.0; kz[:, nz] = 0.0               # symmetry at mid-height; far boundaries closed (far enough not to matter)
    rows, cols, vals = [], [], []
    b_in = np.zeros(nx * nz)
    ax, az = alpha / dx ** 2, alpha / dz ** 2
    for i in range(nx):
        for j in range(nz):
            p = idx(i, j); diag = 0.0
            for (ii, jj, c) in ((i - 1, j, ax * kx[i, j]), (i + 1, j, ax * kx[i + 1, j]), (i, j - 1, az * kz[i, j]), (i, j + 1, az * kz[i, j + 1])):
                if c > 0.0:
                    rows.append(p); cols.append(idx(ii, jj)); vals.append(c); diag -= c
            if inzone[i, j]:                                     # advection, central differences (cell Peclet number below 2); inflow at T = 1
                w = v_T / dx
                left_in = i == nx_u; right_out = i == nx_u + nx_z - 1
                # flux in at the left face
                if left_in:
                    b_in[p] += w * 1.0
                else:
                    rows.append(p); cols.append(idx(i - 1, j)); vals.append(0.5 * w); diag += 0.5 * w
                # flux out at the right face
                if right_out:
                    diag -= w
                else:
                    rows.append(p); cols.append(idx(i + 1, j)); vals.append(-0.5 * w); diag -= 0.5 * w
            rows.append(p); cols.append(p); vals.append(diag)
    A = sp.csc_matrix((vals, (rows, cols)), shape=(nx * nz, nx * nz))
    dt = dt_yr * YEAR
    lu = splu((sp.identity(nx * nz, format="csc") - dt * A).tocsc())
    T = np.zeros(nx * nz)
    n = int(round(t_end_yr / dt_yr))
    t_hist, out_hist = [0.0], [0.0]
    last = [idx(nx_u + nx_z - 1, j) for j in range(nz_z)]
    for step in range(n):
        T = lu.solve(T + dt * b_in)
        t_hist.append((step + 1) * dt_yr); out_hist.append(float(np.mean(T[last])))
    t_hist, out_hist = np.array(t_hist), np.array(out_hist)
    cooled_in = float(T.reshape(nx, nz)[inzone].sum() * dx * dz); cooled_all = float(T.sum() * dx * dz)
    res = dict(TD=np.interp(times, t_hist, out_hist).tolist() if times is not None else None, cooled_zone=cooled_in / (L * 0.5 * B), cooled_outside=(cooled_all - cooled_in) / (L * 0.5 * B))
    return res, t_hist, out_hist


def outside():
    """How much the rock around the swept block adds, for the four net-target designs and for the Project Red block as operated.
    Estimate: heat conducted through a face whose temperature fell by dT at time t_c is rho c dT 2 sqrt(alpha (t - t_c) / pi) per unit area;
    with the front moving at constant speed the top and bottom faces add (4/3) x delta / B of the block's heat and the injection side delta / L,
    delta = 2 sqrt(alpha t / pi). Check: the two-dimensional balance of block_2d with and without conducting surroundings."""
    from study_rock import need_designs
    import field_cases as fc
    cases = {}
    for sid, nd in need_designs().items():
        s = SYS[sid]; wings = s.get("wings", 1)
        V = nd["W"] * nd["L"] * s["B"] * wings
        t_s = V * s["rhoc"] / (nd["m"] * cbar(s)) / YEAR
        cases[str(sid)] = dict(label=s["label"] + " net-target design", L=nd["L"], B=s["B"], W=nd["W"], alpha=s["k"] / s["rhoc"], t_s=t_s, t_end=30.0)
    r = fc.RED
    from aharrsim.reservoir import cbar_w
    cw = cbar_w(r["T_inj"], r["T_rock"] - 6.0, r["p_res"])
    cases["red"] = dict(label="Project Red block at 36 kg/s", L=r["L"], B=r["B"], W=r["W"], alpha=r["k"] / r["rhoc"], t_s=r["W"] * r["L"] * r["B"] * r["rhoc"] / (r["m"] * cw) / YEAR, t_end=8.0)
    out = dict(note="temperature scaled: 0 rock, 1 injection; TD = share of the temperature rise lost at the outlet", cases={})
    for key, c in cases.items():
        t0 = time.time()
        L, B, al, t_s = c["L"], c["B"], c["alpha"], c["t_s"]
        v_T = L / (t_s * YEAR)
        times = np.arange(0.5, 2.5 * max(c["t_end"], t_s) + 0.01, 0.5)
        dx = min(4.0, L / 40.0); dt_yr = 0.02 if t_s > 10 else 0.005
        t_run = float(times[-1])
        on, th, oh = block_2d(L, B, v_T, al, t_run, dx=dx, outside=True, times=times, dt_yr=dt_yr)
        off, th0, oh0 = block_2d(L, B, v_T, al, t_run, dx=dx, outside=False, times=times, dt_yr=dt_yr)
        t10 = lambda th_, oh_: float(np.interp(0.1, oh_, th_)) if oh_[-1] > 0.1 else None
        t10_on, t10_off = t10(th, oh), t10(th0, oh0)
        te = c["t_end"]; x = te / t_s; delta = 2.0 * math.sqrt(al * te * YEAR / math.pi)
        est_tb = (4.0 / 3.0) * min(x, 1.0) * delta / B if x <= 1.0 else None
        est_inj = delta / L
        # the block that lasts as long without help from outside is longer by the ratio of the two times
        res = dict(c, x_end=x, y_note="dense limit", delta_m=delta, est_top_bottom=est_tb, est_injection_side=est_inj,
                   TD_end_with=float(np.interp(te, th, oh)), TD_end_without=float(np.interp(te, th0, oh0)), t_10pct_with=t10_on, t_10pct_without=t10_off,
                   life_ratio=(t10_on / t10_off if t10_on and t10_off else None), heat_from_outside_at_end=None, wall_s=round(time.time() - t0, 1))
        # heat that came from outside the zone by t_end, as a share of the zone's heat: rerun to t_end only
        on_e, _, _ = block_2d(L, B, v_T, al, te, dx=dx, outside=True, dt_yr=dt_yr)
        res["heat_from_outside_at_end"] = on_e["cooled_outside"]
        out["cases"][key] = res
        print(key, c["label"], {k: (round(v, 3) if isinstance(v, float) else v) for k, v in res.items() if k not in ("label",)}, flush=True)
        save("RV5_outside_conduction", out)


def block_slab(L, B, S, v_T, alpha, t_end_yr, nx_z=96, dz=6.0, n_eta=7, up=160.0, down=80.0, out_z=160.0, dt_yr=0.02, long=True, outside=True, mark_yr=None):
    """The swept block with its fractures resolved: one half-slab of rock between a fracture (eta = 0) and the mid-plane of the slab
    (eta = S/2), along the flow direction (xi) and over the height (z, symmetric about mid-height).

        rock:      dT/dt = alpha (d2T/dxi2 + d2T/dz2 + d2T/deta2)
        fracture:  the water of one half-fracture is lumped with the first, thin layer of rock; it moves the heat at v_T (S/2) / h_0 along xi

    long=False switches off conduction along the flow path inside the block and outside=False the rock around the block: with both off this
    is the parallel-fracture element of Gringarten et al. (1975), which serves as the check. Temperature is scaled: 0 = rock, 1 = injection. Lengths in m, times in years (v_T in m/yr, alpha in m2/yr).
    Second-order backward differences in time, central differences in space. Returns the time [yr] and the drop at the outlet."""
    import scipy.sparse as sp
    from scipy.sparse.linalg import splu
    if not outside:
        up = down = out_z = 0.0
    dx = L / nx_z
    nx_u, nx_d = int(round(up / dx)), int(round(down / dx))
    nz_z = max(int(round(0.5 * B / dz)), 2); dz = 0.5 * B / nz_z; nz_o = int(round(out_z / dz))
    nx, nz, ne = nx_u + nx_z + nx_d, nz_z + nz_o, n_eta
    h = 1.6 ** np.arange(ne); h = h / h.sum() * 0.5 * S                          # layers across the half-slab, thin at the fracture
    I, J, M = np.meshgrid(np.arange(nx), np.arange(nz), np.arange(ne), indexing="ij")
    P = (I * nz + J) * ne + M
    zone = (I >= nx_u) & (I < nx_u + nx_z) & (J < nz_z)
    rows, cols, vals = [], [], []
    diag = np.zeros(P.size)

    def couple(a, b, ca, cb, on):
        """Conduction between the cells a and b (index arrays): coefficient ca in the balance of a, cb in that of b."""
        a, b, ca, cb = a[on], b[on], np.broadcast_to(ca, on.shape)[on], np.broadcast_to(cb, on.shape)[on]
        rows.extend([a, b]); cols.extend([b, a]); vals.extend([ca, cb])
        np.add.at(diag, a, -ca); np.add.at(diag, b, -cb)

    both = zone[:-1] & zone[1:]                                                 # along the flow path
    couple(P[:-1], P[1:], alpha / dx ** 2, alpha / dx ** 2, np.where(both, long, outside))
    both = zone[:, :-1] & zone[:, 1:]                                           # over the height
    couple(P[:, :-1], P[:, 1:], alpha / dz ** 2, alpha / dz ** 2, np.where(both, True, outside))
    hm, hp = h[:-1], h[1:]                                                      # across the slab
    couple(P[:, :, :-1], P[:, :, 1:], 2.0 * alpha / (hm * (hm + hp)), 2.0 * alpha / (hp * (hm + hp)), np.ones(P[:, :, :-1].shape, bool))
    rows, cols, vals = np.concatenate(rows), np.concatenate(cols), np.concatenate(vals)
    b_in = np.zeros(P.size)
    w = v_T * 0.5 * S / h[0] / dx                                               # the water: central differences, inflow at T = 1
    f = [list(P[nx_u + i, :nz_z, 0]) for i in range(nx_z)]
    ar, ac, av = [], [], []
    for i in range(nx_z):
        here = np.array(f[i])
        if i == 0:
            b_in[here] += w
        else:
            ar.append(here); ac.append(np.array(f[i - 1])); av.append(np.full(here.size, 0.5 * w)); diag[here] += 0.5 * w
        if i == nx_z - 1:
            diag[here] -= w
        else:
            ar.append(here); ac.append(np.array(f[i + 1])); av.append(np.full(here.size, -0.5 * w)); diag[here] -= 0.5 * w
    rows = np.concatenate([rows] + ar + [np.arange(P.size)]); cols = np.concatenate([cols] + ac + [np.arange(P.size)]); vals = np.concatenate([vals] + av + [diag])
    A = sp.csc_matrix((vals, (rows, cols)), shape=(P.size, P.size))
    dt = dt_yr                                                                  # years: v_T in m/yr and alpha in m2/yr
    eye = sp.identity(P.size, format="csc")
    lu1, lu2 = splu((eye - dt * A).tocsc()), splu((1.5 * eye - dt * A).tocsc())
    last = np.array(f[-1])
    vol = (h[None, None, :] * dx * dz) * np.ones(P.shape); ref = L * 0.5 * B * 0.5 * S

    def cooled(T):                                                              # heat given up by the block and by the rock around it, as shares of the block's heat
        c = T.reshape(P.shape) * vol
        return float(c[zone].sum()) / ref, float(c[~zone].sum()) / ref

    T0 = np.zeros(P.size); T1 = lu1.solve(T0 + dt * b_in)
    t_hist, out_hist = [0.0, dt_yr], [0.0, float(T1[last].mean())]
    info = dict(cells=int(P.size))
    for step in range(1, int(round(t_end_yr / dt_yr))):
        T2 = lu2.solve(2.0 * T1 - 0.5 * T0 + dt * b_in)
        T0, T1 = T1, T2
        t_hist.append((step + 1) * dt_yr); out_hist.append(float(T1[last].mean()))
        if mark_yr is not None and abs((step + 1) * dt_yr - mark_yr) < 0.5 * dt_yr:
            info["cooled_zone_mark"], info["cooled_outside_mark"] = cooled(T1)
    info["cooled_zone"], info["cooled_outside"] = cooled(T1)
    info.update(T_min=float(T1.min()), T_max=float(T1.max()))
    return np.array(t_hist), np.array(out_hist), info


def paths():
    """Two conduction paths that the parallel-fracture element leaves out, added in turn to a block with its fractures resolved (block_slab):
    conduction along the flow path, and the rock around the block. The four net-target designs and the Project Red block at 36 kg/s.
    The first column of the result is the element itself (the analytical solution), the second the same case in the resolved balance."""
    from study_rock import need_designs
    from aharrsim.reservoir import ParallelFractures, cbar_w
    import field_cases as fc
    t10 = lambda t, o: float(np.interp(0.1, o, t)) if o[-1] > 0.1 else None
    cases = {}
    for sid, nd in need_designs().items():
        s = SYS[sid]; wings = s.get("wings", 1); cw = cbar(s)
        cases[str(sid)] = dict(label=s["label"] + " net-target design", N=nd["N"] * wings, m=nd["m"], W=nd["W"], L=nd["L"], B=s["B"], k=s["k"], rhoc=s["rhoc"], T_rock=s["T_rock"],
                               T_inj=s["T_inj"], rhoc_w=cw * s["rho_w"], Q=nd["m"] / s["rho_w"], t_end=45.0)
    r = fc.RED; cw = cbar_w(r["T_inj"], r["T_rock"] - 6.0, r["p_res"]); rho_w = 920.0
    cases["red"] = dict(label="Project Red block at 36 kg/s", N=r["N"], m=r["m"], W=r["W"], L=r["L"], B=r["B"], k=r["k"], rhoc=r["rhoc"], T_rock=r["T_rock"], T_inj=r["T_inj"],
                        rhoc_w=cw * rho_w, Q=r["m"] / rho_w, t_end=10.0)
    out = dict(note="time [yr] at which the outlet has lost a tenth of its temperature rise; scaled temperature drop at 30 years (8 years for Project Red)", cases={})
    for key, c in cases.items():
        t0 = time.time()
        S = c["W"] / (c["N"] / (2 if key == "141" else 1)); al = c["k"] / c["rhoc"]
        v_T = c["Q"] * c["rhoc_w"] / (c["N"] * c["B"] * S * c["rhoc"]) * YEAR                     # m/yr: front speed of the block
        t_s = c["L"] / v_T
        tt = np.arange(0.25, c["t_end"] + 0.01, 0.25)
        pf = ParallelFractures(Q=c["Q"], N=c["N"], L_flow=c["L"], B=c["B"], S=S, k_r=c["k"], rhoc_r=c["rhoc"], T_rock_C=c["T_rock"], T_inj_C=c["T_inj"], rhoc_w=c["rhoc_w"])
        TD_a = np.maximum.accumulate(np.asarray(pf.T_D(tt * YEAR), float))
        dt_yr = 0.02 if t_s > 10 else 0.005
        te = 30.0 if key != "red" else 8.0
        res = dict(label=c["label"], S_m=S, L=c["L"], B=c["B"], t_s=t_s, Pe=v_T * c["L"] / (al * YEAR), x_end=te / t_s, element=dict(t10=t10(tt, TD_a), TD_end=float(np.interp(te, tt, TD_a))))
        for name, kw in (("resolved_as_element", dict(long=False, outside=False)), ("with_path", dict(long=True, outside=False)), ("with_path_and_surroundings", dict(long=True, outside=True))):
            t, o, info = block_slab(c["L"], c["B"], S, v_T, al * YEAR, c["t_end"], dt_yr=dt_yr, dz=6.0 if c["B"] < 200 else 8.0, mark_yr=te, **kw)
            res[name] = dict(t10=t10(t, o), TD_end=float(np.interp(te, t, o)), **info)
            # the time at which this case has lost what the element has lost at the end: at a fixed flow, te over that time is roughly the volume that restores the element's result
            res[name]["t_at_element_end"] = float(np.interp(res["element"]["TD_end"], o, t)) if o[-1] > res["element"]["TD_end"] else None
            if name == "with_path_and_surroundings":
                res[name]["heat_from_outside_at_end"] = info["cooled_outside_mark"]
        res["wall_s"] = round(time.time() - t0, 1)
        out["cases"][key] = res
        print(key, c["label"], "S %.1f Pe %.0f t_s %.1f |" % (S, res["Pe"], t_s), {k: (round(v["t10"], 2) if v.get("t10") else None, round(v["TD_end"], 3), round(v.get("t_at_element_end") or 0, 2)) for k, v in res.items() if isinstance(v, dict)},
              "| outside", round(res["with_path_and_surroundings"]["heat_from_outside_at_end"], 3), "|", res["wall_s"], "s", flush=True)
        save("RV8_conduction_paths", out)


def paths_grid():
    """Grid and time-step check of block_slab on the type 14 net-target design: the element case against the analytical solution (cheap, so
    refined twice), and the case with both conduction paths about a coarser grid than the working one (96 steps along the path)."""
    from study_rock import need_designs
    from aharrsim.reservoir import ParallelFractures
    nd = need_designs()[14]; s = SYS[14]; S = nd["W"] / nd["N"]; al = s["k"] / s["rhoc"] * YEAR
    cw = cbar(s); Q = nd["m"] / s["rho_w"]
    v_T = Q * cw * s["rho_w"] / (nd["N"] * s["B"] * S * s["rhoc"]) * YEAR
    t10 = lambda t, o: float(np.interp(0.1, o, t)) if o[-1] > 0.1 else None
    tt = np.arange(0.25, 45.01, 0.25)
    pf = ParallelFractures(Q=Q, N=nd["N"], L_flow=nd["L"], B=s["B"], S=S, k_r=s["k"], rhoc_r=s["rhoc"], T_rock_C=s["T_rock"], T_inj_C=s["T_inj"], rhoc_w=cw * s["rho_w"])
    TDa = np.maximum.accumulate(np.asarray(pf.T_D(tt * YEAR), float))
    out = dict(case="type 14 net-target design", analytical=dict(t10=t10(tt, TDa), TD_30=float(np.interp(30.0, tt, TDa))), element=[], both_paths=[])
    for lab, g in (("48 steps along the path", dict(nx_z=48)), ("96 steps (working grid)", dict(nx_z=96)), ("192 steps", dict(nx_z=192)), ("96 steps, half the time step", dict(nx_z=96, dt_yr=0.01)),
                   ("96 steps, ten layers across the slab", dict(nx_z=96, n_eta=10)), ("96 steps, half the step over the height", dict(nx_z=96, dz=3.0))):
        t, o, info = block_slab(nd["L"], s["B"], S, v_T, al, 45.0, long=False, outside=False, **g)
        row = dict(grid=lab, t10=t10(t, o), TD_30=float(np.interp(30.0, t, o)), max_diff=float(np.max(np.abs(np.interp(tt, t, o) - TDa))), cells=info["cells"])
        out["element"].append(row); print("element", row, flush=True)
    for lab, g in (("48 steps along the path", dict(nx_z=48)), ("96 steps (working grid)", dict(nx_z=96)), ("48 steps, half the time step", dict(nx_z=48, dt_yr=0.01)),
                   ("48 steps, ten layers across the slab", dict(nx_z=48, n_eta=10)), ("48 steps, half the step over the height", dict(nx_z=48, dz=3.0)),
                   ("48 steps, surroundings twice as far", dict(nx_z=48, up=320.0, down=160.0, out_z=320.0))):
        t0 = time.time(); t, o, info = block_slab(nd["L"], s["B"], S, v_T, al, 45.0, long=True, outside=True, mark_yr=30.0, **g)
        row = dict(grid=lab, t10=t10(t, o), TD_30=float(np.interp(30.0, t, o)), outside_30=info["cooled_outside_mark"], cells=info["cells"], wall_s=round(time.time() - t0, 1))
        out["both_paths"].append(row); print("both paths", row, flush=True)
        save("RV8_conduction_grid", out)


def boot(n_boot=4000, seed=1):
    """Percentiles of the closed-loop Latin-hypercube samples (80, 40 and 40 runs, conductivity falling with temperature): an ordinary
    bootstrap (the run outputs resampled with replacement) of P10, P50 and P90 of the output at year 30; the keys *_90ci hold the 5th and 95th
    percentiles of the resampled values. A descriptive diagnostic of the sampling scatter, not a confidence interval with established coverage:
    Latin-hypercube runs are stratified, not independent."""
    rng = np.random.default_rng(seed)
    MC = [r for r in load("ST1_mc") if r["case"] != 8] + load("ST1_mc_8_m30")
    out = dict(n_boot=n_boot, seed=seed, year=30, note="ordinary bootstrap of the run outputs; *_90ci = 5th and 95th percentiles of the resampled percentile; descriptive, "
               "not a confidence interval (Latin-hypercube runs are stratified)", cases={})
    for c, key in ((7, "MWth"), (8, "MWe"), (20, "MWe")):
        v = np.array([r["30"][key] for r in MC if r["case"] == c and "failed" not in r], float)
        bs = v[rng.integers(0, len(v), size=(n_boot, len(v)))]
        pc = np.percentile(bs, [10, 50, 90], axis=1)
        o = dict(output=key, n=int(len(v)), kT=sorted({bool(r["P"].get("kT")) for r in MC if r["case"] == c}))
        for name, q, row in (("P10", 10, pc[0]), ("P50", 50, pc[1]), ("P90", 90, pc[2])):
            o[name] = float(np.percentile(v, q)); o[name + "_90ci"] = [float(np.percentile(row, 5)), float(np.percentile(row, 95))]
        # the same from each half of the sample, as a plainer check
        h1, h2 = v[::2], v[1::2]
        o["halves"] = dict(P10=[float(np.percentile(h1, 10)), float(np.percentile(h2, 10))], P50=[float(np.percentile(h1, 50)), float(np.percentile(h2, 50))],
                           P90=[float(np.percentile(h1, 90)), float(np.percentile(h2, 90))])
        out["cases"][str(c)] = o
        print(c, {k: (np.round(x, 3).tolist() if k[:1] == "P" else x) for k, x in o.items() if k != "halves"}, "halves", {k: np.round(x, 3).tolist() for k, x in o["halves"].items()}, flush=True)
    save("RV6_percentiles", out)


if __name__ == "__main__":
    {"gross": gross, "share": share, "seeds": seeds, "red": red, "outside": outside, "paths": paths, "paths_grid": paths_grid, "boot": boot}[sys.argv[1]]()
