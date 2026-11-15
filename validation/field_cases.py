"""Field evidence: the published operating record of Fervo's Project Red (two years of production, Fervo Energy, 13 April 2026)
and the 2024 Cape Station test, read with the same model and the same two groups as the screening chart.

Project Red as operated: 36 kg/s (average circulation), injection about 76 C, stimulated volume as published
(975 m of lateral x 110 m between the wells x 143 m of fracture height), 88 flowing fractures (75-100 estimated).
The measured flowing temperature is digitised from the published figure (results/FE_project_red_digitised.json, +/- 1 K).

  f    fraction of the fractures that carry the flow (rock volume kept: the active fractures are S / f apart)
  phi  fraction of the published volume that is swept (rock volume lost)

    python field_cases.py        -> results/FE_field_cases.json
"""
from __future__ import annotations

import math

import numpy as np
import CoolProp.CoolProp as CP

from common import save, load
from aharrsim.reservoir import ParallelFractures, cbar_w
from aharrsim import plant

YEAR = 365.25 * 86400.0
RHO = 900.0

RED = dict(m=36.0, T_inj=76.0, W=975.0, L=110.0, B=143.0, N=88, k=2.8, rhoc=2650 * 850.0, T_rock=191.0, p_res=22.0e6,
           shut_ins=(0.0, 0.47, 0.57, 0.75, 1.0, 1.27), days=614)


def red_model(t_yr, f=1.0, phi=1.0, m=None, r=RED):
    """Dimensionless outlet temperature drop of the Project Red block at times t_yr."""
    m = m or r["m"]
    cw = cbar_w(r["T_inj"], r["T_rock"] - 6.0, r["p_res"])
    Na = f * r["N"]
    pf = ParallelFractures(Q=m / RHO, N=Na, L_flow=r["L"], B=phi * r["B"], S=r["W"] / Na, k_r=r["k"], rhoc_r=r["rhoc"], T_rock_C=r["T_rock"], T_inj_C=r["T_inj"],
                           rhoc_w=cw * RHO)
    t = np.asarray(t_yr, float) * YEAR
    return pf.T_D(t), pf


def red():
    rec = load("FE_project_red_digitised")
    t = np.array(rec["t_yr"]); T = np.array(rec["T_prod_C"])
    r = RED
    ok = t > 0.42
    for s in r["shut_ins"]:                                   # the flowing temperature needs a few weeks to recover after each shut-in
        ok &= ~((t >= s - 0.005) & (t < s + 0.07))
    plateau = float(np.median(T[(t > 0.45) & (t < 0.98) & ok]))
    dT = plateau - r["T_inj"]
    TD_obs = (plateau - T) / dT
    V = r["W"] * r["L"] * r["B"]
    cw = cbar_w(r["T_inj"], r["T_rock"] - 6.0, r["p_res"])
    t_s = V * r["rhoc"] / (r["m"] * cw) / YEAR
    al = r["k"] / r["rhoc"]
    t_end = r["days"] / 365.25
    x_end = t_end / t_s; y_end = 0.5 * (r["W"] / r["N"]) / math.sqrt(al * t_end * YEAR)
    E_out = r["m"] * cw * dT * t_end * YEAR                   # heat produced so far (plateau temperature)
    E_place = V * r["rhoc"] * (r["T_rock"] - r["T_inj"])       # heat in the published volume above the injection temperature
    tt = np.linspace(0.05, 12.0, 240)

    def rms(f, phi):
        TD, _ = red_model(t[ok], f, phi)
        return float(np.sqrt(np.mean((TD - TD_obs[ok]) ** 2)) * dT)

    fs = np.round(np.arange(0.15, 1.001, 0.01), 3); phis = np.round(np.arange(0.20, 1.201, 0.01), 3)
    grid = np.array([[rms(f, p) for p in phis] for f in fs])
    i, j = np.unravel_index(np.argmin(grid), grid.shape)
    best = dict(f=float(fs[i]), phi=float(phis[j]), rms_K=float(grid[i, j]))
    jf = int(np.argmin(np.abs(phis - 1.0))); i1 = int(np.argmin(grid[:, jf])); only_f = dict(f=float(fs[i1]), phi=1.0, rms_K=float(grid[i1, jf]))
    i2 = len(fs) - 1; j2 = int(np.argmin(grid[i2])); only_phi = dict(f=1.0, phi=float(phis[j2]), rms_K=float(grid[i2, j2]))
    ideal = dict(f=1.0, phi=1.0, rms_K=rms(1.0, 1.0))
    # the decline at the end of the record, as drawn (figure) and as stated (text: about 2.5 F), turned into an effective fraction
    end = t > t_end - 0.06
    drop_fig = float(plateau - np.mean(T[end])); drop_txt = 2.5 / 1.8
    eff = {}
    for name, drop in (("figure", drop_fig), ("text", drop_txt)):
        target = drop / dT
        from scipy.optimize import brentq
        ff = brentq(lambda f: float(red_model([t_end], f, 1.0)[0][0]) - target, 0.05, 1.0)
        pp = brentq(lambda p: float(red_model([t_end], 1.0, p)[0][0]) - target, 0.05, 1.5)
        eff[name] = dict(drop_K=drop, T_D=target, f_if_all_volume=ff, active_fractures=ff * r["N"], phi_if_all_fractures=pp)
    cases = dict(ideal=ideal, fit_fractures=only_f, fit_volume=only_phi, fit_both=best)
    curves = {}
    for name, c in cases.items():
        TD, pf = red_model(tt, c["f"], c["phi"])
        Tm = plateau - TD * dT
        below = np.where(TD > 0.1)[0]; cold = np.where(Tm < 150.0)[0]
        c.update(T_C=Tm.tolist(), t_10pct_yr=float(tt[below[0]]) if len(below) else None, t_150C_yr=float(tt[cold[0]]) if len(cold) else None,
                 T_at={str(y): float(np.interp(y, tt, Tm)) for y in (2, 3, 4, 5, 6, 8, 10)}, t_sweep_yr=t_s * c["phi"],
                 x_end=t_end / (t_s * c["phi"]), y_end=y_end / c["f"])
        curves[name] = c
    # acceptable (f, phi) pairs: rms within the digitising accuracy (1 K)
    acc = [(float(fs[a]), float(phis[b])) for a in range(len(fs)) for b in range(len(phis)) if grid[a, b] < 1.0]
    prod = np.array([f * p for f, p in acc])
    # power at the operated conditions with the plant model of the benchmark
    hw = CP.PropsSI("H", "T", plateau + 273.15, "P", 2.0e6, "Water")
    gross = {f"{kind}_{int(Ta)}C": plant.orc_practical(r["m"], hw, 2.0e6, Ta, kind=kind)["W_net"] / 1e6 for kind in ("orc", "sorc") for Ta in (5.0, 15.0, 25.0)}
    hw60 = CP.PropsSI("H", "T", 184.0 + 273.15, "P", 2.0e6, "Water")
    out = dict(source="Fervo Energy, 13 April 2026 (figure digitised from the GEOPHIRES-X documentation copy)", plateau_C=plateau, T_inj_C=r["T_inj"], dT_K=dT,
               t_obs=t.tolist(), T_obs=T.tolist(), used=ok.tolist(), TD_obs=TD_obs.tolist(), V_km3=V / 1e9, t_sweep_yr=t_s, cbar_w=cw, x_end=x_end, y_end=y_end, t_end_yr=t_end,
               recovery_to_date=E_out / E_place, reported=dict(recovery=0.20, gross_MW=2.1, net_MW=1.4, flow_kg_s=36.0, T_prod_C=175.0, peak_gross_MW=3.5, peak_flow=60.0,
                                                                decline_text_F=2.5, onset_days=500),
               effective=eff, t_yr=tt.tolist(), cases=curves, acceptable=dict(n=len(acc), f=[min(a[0] for a in acc), max(a[0] for a in acc)] if acc else None,
                                                                                phi=[min(a[1] for a in acc), max(a[1] for a in acc)] if acc else None,
                                                                                f_phi=[float(prod.min()), float(prod.max())] if acc else None),
               rms_grid=dict(f=fs.tolist(), phi=phis.tolist(), rms_K=grid.tolist()),
               gross_model_MW=gross, gross_model_peak_MW=plant.orc_practical(60.0, hw60, 2.0e6, 15.0)["W_net"] / 1e6)
    return out


def cape():
    """Cape Station: the 2024 test as reported against the benchmark's type 14b, and what the per-well figures need."""
    dt = {float(x["id"]): x for x in load("DT_types")}[14.5]
    a = dt["assumptions"]
    rows = {}
    for label, m, T in (("test, peak flow", 107.0, 195.0), ("test, average", 100.0, 195.0), ("design point 3.0", 100.0, 216.0)):
        h = CP.PropsSI("H", "T", T + 273.15, "P", 2.07e6, "Water")
        rows[label] = dict(m=m, T_wh=T, **{f"{kind}_{int(Ta)}C": plant.orc_practical(m, h, 2.07e6, Ta, kind=kind)["W_net"] / 1e6 for kind in ("orc", "sorc") for Ta in (5.0, 15.0)})
    # what a 16 MWe gross producer needs at 221 C rock (216 C at the wellhead): flow, and rock volume for 10 and 30 years at x = 0.74 (y = 0.15)
    need = {}
    for kind in ("orc", "sorc"):
        h = CP.PropsSI("H", "T", 216.0 + 273.15, "P", 2.07e6, "Water")
        w = plant.orc_practical(1.0, h, 2.07e6, 15.0, kind=kind)["W_net"] / 1e6
        m = 16.0 / w
        cw = cbar_w(60.0, 221.0, 30e6)
        V = lambda yrs: m * cw * yrs * YEAR / (2650 * 850.0 * 0.74) / 1e9
        need[kind] = dict(w_MJ_per_kg=w, flow_kg_s=m, V_km3_10yr=V(10.0), V_km3_30yr=V(30.0), section_m2_per_m_lateral_30yr=V(30.0) * 1e9 / 2286.0)
    y1 = dt["T_out"]["1"]; pl = dt["plant"]
    return dict(reported=dict(flow_peak_kg_s=107.0, flow_range=[93.0, 121.0], T_prod_C=195.0, gross_MW=">10 (Fervo); >9.5 (McClure)", inj_WHP_MPa=[13.8, 15.9],
                              lateral_m=1524.0, spacing_m=137.0, phase1_MW=100.0, in_operation_MW=36.0, per_well_design_MW=16.0),
                model=dict(V_km3=a["V_km3"], t_s_yr=a["t_s_yr"], T_wh_yr1=y1["T_C"], T_out={k: v["T_C"] for k, v in dt["T_out"].items()},
                           pump_kW=dt.get("well_pump_kW"), impedance=a["impedance_MPa_per_L_s"]),
                plant=rows, need_16MW=need)


GER = dict(Lv=4235.0, L_lat=3500.0, pairs_built=6, pairs_working=4, pairs_design=12, T_s=10.0, t_now_yr=0.45, Q_now_MWth=8.5, Q_avg_MWth=11.8, eta_orc=0.043,
           Q_pair_now_MWth=2.0, Q_pair_5yr_MWth=1.3, MWe_reported=0.5, nameplate_MWe=8.2, nameplate_MWth=64.0, loops=4, pair_km_design=6.0, T_in=50.0, T_amb=10.0)


def _ger_loop(n_pairs, L_pair, T_lat, m, t_yr, k, T_in=None, rhoc=2.3e6, g=GER):
    import screening as SC
    T_in = g["T_in"] if T_in is None else T_in
    segs = [dict(L=g["Lv"], T0=g["T_s"], T1=T_lat, r_b=0.156, k=k, rhoc=rhoc), dict(L=L_pair, T0=T_lat, T1=T_lat, r_b=0.108, k=k, rhoc=rhoc, n=n_pairs),
            dict(L=g["Lv"], T0=T_lat, T1=g["T_s"], r_b=0.156, k=k, rhoc=rhoc)]
    T_out, ntu, c = SC.chain(segs, m, T_in, t_yr, 25e6)
    Q = m * c * (T_out - T_in) / 1e6
    W = plant.orc_practical(m, CP.PropsSI("H", "T", T_out + 273.15, "P", 2.0e6, "Water"), 2.0e6, g["T_amb"])["W_net"] / 1e6 if T_out > 71 else 0.0
    return dict(T_out=T_out, Q_MWth=Q, MWe=W, NTU=ntu[1], eta=W / Q if Q > 0 else 0.0)


def _ger_mean(n_pairs, L_pair, T_lat, m, k, t_end):
    """Time averages from start-up to t_end [yr] at constant flow (chart): heat, electricity, outlet temperature."""
    tt = np.concatenate([[6.0 / 24.0 / 365.25], np.geomspace(1.0 / 365.25, t_end, 60)])
    r = [_ger_loop(n_pairs, L_pair, T_lat, m, t, k) for t in tt]
    I = lambda key: float(np.sum(0.5 * (np.array([x[key] for x in r])[1:] + np.array([x[key] for x in r])[:-1]) * np.diff(tt)) / (tt[-1] - tt[0]))
    Q, W = I("Q_MWth"), I("MWe")
    return dict(Q_MWth=Q, MWe=W, T_out=I("T_out"), eta=W / Q if Q > 0 else 0.0)


def geretsried(full=True):
    """Eavor-Loop 1 at Geretsried as reported (Toews, Eavor, 21 May 2026; Whitfield, Drilling Contractor, 22 April 2026) against the conduction
    chart: the flow is not published, so it is set to give the reported heat; the outlet temperature, the conversion and the five-year heat follow."""
    from scipy.optimize import brentq
    g = GER
    rows = []
    for T_lat in (150.0, 160.0, 170.0):
        for k in (2.5, 3.0, 3.5):
            for L_pair in (6000.0, 7000.0):
                f = lambda m: _ger_loop(g["pairs_working"], L_pair, T_lat, m, g["t_now_yr"], k)["Q_MWth"] - g["Q_now_MWth"]
                try:
                    m = brentq(f, 5.0, 400.0)
                except ValueError:
                    continue
                now = _ger_loop(g["pairs_working"], L_pair, T_lat, m, g["t_now_yr"], k); y5 = _ger_loop(g["pairs_working"], L_pair, T_lat, m, 5.0, k)
                rows.append(dict(T_lat=T_lat, k=k, L_pair=L_pair, m=m, T_out=now["T_out"], MWe=now["MWe"], eta=now["eta"], NTU=now["NTU"],
                                 Q_pair_5yr=y5["Q_MWth"] / g["pairs_working"], T_out_5yr=y5["T_out"]))
    central = [r for r in rows if r["T_lat"] == 160.0 and r["k"] == 3.0]
    # the design loop (12 pairs of 6 km at 160 C, k 3): most electricity the chart allows at 5 and 30 years, and the heat at the reported operating point
    design = {}
    for yr in (5.0, 30.0):
        best = max((_ger_loop(12, 6000.0, 160.0, m, yr, 3.0) | dict(m=float(m)) for m in np.arange(10.0, 200.1, 2.5)), key=lambda r: r["MWe"])
        design[str(int(yr))] = best
    # the update also gives the average thermal duty since start-up (11.8 MWth, on which the 4.3 % conversion is stated):
    # (a) the mean the model gives at the flow fitted to today's heat; (b) the constant flow whose mean is the reported average
    c6 = [r for r in central if r["L_pair"] == 6000.0][0]
    mean_fit = _ger_mean(g["pairs_working"], 6000.0, 160.0, c6["m"], 3.0, g["t_now_yr"])
    avg_rows = []
    for T_lat in (150.0, 160.0, 170.0):
        for k in (2.5, 3.0, 3.5):
            for L_pair in (6000.0, 7000.0):
                f = lambda m: _ger_mean(g["pairs_working"], L_pair, T_lat, m, k, g["t_now_yr"])["Q_MWth"] - g["Q_avg_MWth"]
                try:
                    m = brentq(f, 5.0, 300.0, xtol=0.05)
                except ValueError:                  # the reported average is beyond what this rock gives four pairs at a usable flow
                    avg_rows.append(dict(T_lat=T_lat, k=k, L_pair=L_pair, m=None))
                    continue
                mean = _ger_mean(g["pairs_working"], L_pair, T_lat, m, k, g["t_now_yr"])
                now = _ger_loop(g["pairs_working"], L_pair, T_lat, m, g["t_now_yr"], k)
                avg_rows.append(dict(T_lat=T_lat, k=k, L_pair=L_pair, m=m, T_out_mean=mean["T_out"], MWe_mean=mean["MWe"], eta_mean=mean["eta"],
                                     Q_now=now["Q_MWth"], T_out_now=now["T_out"]))
    ok = [r for r in avg_rows if r["m"] is not None]
    average = dict(reported_MWth=g["Q_avg_MWth"], reported_MWe=g["Q_avg_MWth"] * g["eta_orc"], at_fitted_flow=dict(m=c6["m"], **mean_fit),
                   rows=avg_rows, central=[r for r in ok if r["T_lat"] == 160.0 and r["k"] == 3.0],
                   ranges=dict(m=[min(r["m"] for r in ok), max(r["m"] for r in ok)], eta=[min(r["eta_mean"] for r in ok), max(r["eta_mean"] for r in ok)],
                               MWe=[min(r["MWe_mean"] for r in ok), max(r["MWe_mean"] for r in ok)], n_unreachable=len(avg_rows) - len(ok)))
    out = dict(reported=dict(g), rows=rows, central=central, average_since_startup=average,
               ranges=dict(m=[min(r["m"] for r in rows), max(r["m"] for r in rows)], T_out=[min(r["T_out"] for r in rows), max(r["T_out"] for r in rows)],
                           MWe=[min(r["MWe"] for r in rows), max(r["MWe"] for r in rows)], eta=[min(r["eta"] for r in rows), max(r["eta"] for r in rows)],
                           Q_pair_5yr=[min(r["Q_pair_5yr"] for r in rows), max(r["Q_pair_5yr"] for r in rows)]),
               design_loop=design, nameplate_per_loop=dict(MWe=g["nameplate_MWe"] / g["loops"], MWth=g["nameplate_MWth"] / g["loops"]),
               heat_at_5yr_per_design_loop_reported=12 * g["Q_pair_5yr_MWth"], MWe_at_reported_conversion=12 * g["Q_pair_5yr_MWth"] * g["eta_orc"])
    if full:
        from aharrsim.systems import uloop, YEAR as YR
        c = [r for r in central if r["L_pair"] == 6000.0][0]
        sims, means = {}, {}
        for lab, ff, rf in (("isolated laterals", "dirichlet", None), ("laterals 50 m apart", "adiabatic", 25.0)):
            ch = uloop(f"Geretsried {lab}", g["Lv"], 6000.0, 0.216, c["m"], g["T_in"], 2.0e6, [dict(k=3.0, rho=2556.0, cp=900.0)], g["T_s"],
                       dict(gradient_C_per_km=(160.0 - g["T_s"]) / (g["Lv"] / 1000.0)), 6.0, n_lat=g["pairs_working"], layers=None, roughness=1e-4, dz=50.0,
                       lat_far_field=ff, lat_r_far=rf)
            ts = ch["producer"][1].ts
            r = {}
            for yr in (g["t_now_yr"], 1.0, 5.0):
                t = yr * YR
                h = float(np.interp(t, ts["t"], ts["h_out"])); T = float(np.interp(t, ts["t"], ts["T_out"])) - 273.15
                Q = c["m"] * (h - CP.PropsSI("H", "T", g["T_in"] + 273.15, "P", 2.0e6, "Water")) / 1e6
                W = plant.orc_practical(c["m"], h, 2.0e6, g["T_amb"])["W_net"] / 1e6
                r[str(yr)] = dict(T_out=T, Q_MWth=Q, Q_pair=Q / g["pairs_working"], MWe=W, eta=W / Q)
            # mean heat from the first day to today in the full model (the first hours flush the water that stood in the hole)
            tt = np.geomspace(86400.0, g["t_now_yr"] * YR, 80)
            Qt = c["m"] * (np.interp(tt, ts["t"], ts["h_out"]) - CP.PropsSI("H", "T", g["T_in"] + 273.15, "P", 2.0e6, "Water")) / 1e6
            means[lab] = float(np.sum(0.5 * (Qt[1:] + Qt[:-1]) * np.diff(tt)) / (tt[-1] - tt[0]))
            sims[lab] = r
        out["full_chain"] = dict(m=c["m"], cases=sims, mean_since_startup_MWth=means)
    return out


def main():
    out = dict(red=red(), cape=cape(), geretsried=geretsried())
    save("FE_field_cases", out)
    r = out["red"]
    print("Project Red: plateau", round(r["plateau_C"], 1), "dT", round(r["dT_K"], 1), "V km3", round(r["V_km3"], 4), "t_s yr", round(r["t_sweep_yr"], 2), "x_end", round(r["x_end"], 3),
          "y_end", round(r["y_end"], 2), "recovery to date", round(r["recovery_to_date"], 3))
    print("  effective:", {k: {a: round(b, 3) for a, b in v.items()} for k, v in r["effective"].items()})
    for k, c in r["cases"].items():
        print(f"  {k:<14} f {c['f']:.2f} phi {c['phi']:.2f} rms {c['rms_K']:.2f} K  10 % drop at {c['t_10pct_yr']} yr, below 150 C at {c['t_150C_yr']} yr; T at", {y: round(v, 1) for y, v in c["T_at"].items()})
    print("  acceptable (rms < 1 K):", r["acceptable"])
    print("  gross model MW", {k: round(v, 2) for k, v in r["gross_model_MW"].items()}, "peak test 60 kg/s", round(r["gross_model_peak_MW"], 2))
    c = out["cape"]
    print("Cape:", c["model"]["V_km3"], c["model"]["t_s_yr"], {k: round(v, 1) for k, v in c["model"]["T_out"].items()})
    for k, v in c["plant"].items():
        print("  ", k, {a: round(b, 2) for a, b in v.items()})
    print("  16 MW needs", {k: {a: round(b, 3) for a, b in v.items()} for k, v in c["need_16MW"].items()})
    g = out["geretsried"]
    a = out["geretsried"]["average_since_startup"]
    print("Geretsried, average since start-up: reported", a["reported_MWth"], "MWth,", round(a["reported_MWe"], 2), "MWe; at the fitted flow",
          {k: round(v, 2) for k, v in a["at_fitted_flow"].items()}, "| flow for the reported average:",
          [{k: round(v, 2) for k, v in r.items()} for r in a["central"]], "ranges", {k: v for k, v in a["ranges"].items()})
    print("Geretsried: flow for 8.5 MWth", [round(v, 1) for v in g["ranges"]["m"]], "T_out", [round(v) for v in g["ranges"]["T_out"]], "MWe", [round(v, 2) for v in g["ranges"]["MWe"]],
          "eta", [round(v, 3) for v in g["ranges"]["eta"]], "Q per pair at 5 yr", [round(v, 2) for v in g["ranges"]["Q_pair_5yr"]])
    print("  central:", [{k: round(v, 2) for k, v in r.items()} for r in g["central"]])
    print("  design loop, best electricity:", {y: {k: round(v, 2) for k, v in r.items()} for y, r in g["design_loop"].items()}, "nameplate per loop", g["nameplate_per_loop"],
          "| 12 pairs x 1.3 MWth =", g["heat_at_5yr_per_design_loop_reported"], "-> MWe at 4.3 %:", round(g["MWe_at_reported_conversion"], 2))
    if "full_chain" in g:
        for lab, r in g["full_chain"]["cases"].items():
            print("  full chain,", lab, {y: {k: round(v, 2) for k, v in d.items()} for y, d in r.items()},
                  "mean since start-up", round(g["full_chain"]["mean_since_startup_MWth"][lab], 2), "MWth")


if __name__ == "__main__":
    main()
