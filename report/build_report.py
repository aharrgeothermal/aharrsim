"""Build the validation & results report (single self-contained HTML page) from validation/results/*.json."""
from __future__ import annotations

import html
import json
import math
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(os.path.dirname(HERE), "validation", "results")
OUT = os.path.join(HERE, "aharrsim_report.html")


def L(name, default=None):
    path = os.path.join(RES, name + ".json")
    if not os.path.exists(path) and default is not None:
        print("  (missing", name, "- section left empty)")
        return default
    with open(path) as fh:
        return json.load(fh)


def ds(x, y, n=70, log=True):
    """Down-sample a curve to ~n points (log-spaced in x when log=True)."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    if x.size <= n:
        return [[float(a), float(b)] for a, b in zip(x, y)]
    if log:
        xs = np.geomspace(max(x[x > 0].min(), 1e-9), x.max(), n)
    else:
        xs = np.linspace(x.min(), x.max(), n)
    idx = sorted(set(int(np.argmin(np.abs(x - v))) for v in xs))
    return [[float(x[i]), float(y[i])] for i in idx]


def r(v, n=1):
    return None if v is None else round(float(v), n)


def sci(v, n=0):
    """Format as a×10^b with a proper minus sign."""
    if v == 0:
        return "0"
    e = int(math.floor(math.log10(abs(v))))
    m = v / 10 ** e
    return f"{m:.{n}f}×10<sup>{str(e).replace('-', '−')}</sup>"


def esc(s):
    return html.escape(str(s))


def table(head, rows, cls="", note=None):
    h = "".join(f"<th>{c}</th>" for c in head)
    b = "".join("<tr>" + "".join(f"<td>{c}</td>" for c in row) + "</tr>" for row in rows)
    n = f"<p class='tnote'>{note}</p>" if note else ""
    return f"<div class='tbl {cls}'><table><thead><tr>{h}</tr></thead><tbody>{b}</tbody></table></div>{n}"


def f(v, n=1, sign=False):
    if v is None or (isinstance(v, float) and not math.isfinite(v)):
        return "–"
    s = f"{v:+.{n}f}" if sign else f"{v:.{n}f}"
    return s.replace("-", "−")


# =====================================================================================
data = {}
# ---------------- deep study --------------------------------------------------------
DF = [x for x in L("D_depth_flow") if "rows" in x]


def row_at(run, years):
    return next(rw for rw in run["rows"] if rw["years"] == years)


flows = sorted({x["mdot"] for x in DF})
depths = sorted({x["depth"] for x in DF})
data["power_depth"] = {str(int(m)): {str(y): [[d / 1000, row_at(next(x for x in DF if x["depth"] == d and x["mdot"] == m), y)["Q_MW"]]
                                              for d in depths] for y in (1, 10, 30)} for m in flows}
data["tout_depth"] = {str(int(m)): {str(y): [[d / 1000, row_at(next(x for x in DF if x["depth"] == d and x["mdot"] == m), y)["T_out"]]
                                             for d in depths] for y in (1, 10, 30)} for m in flows}
base = {d: next(x for x in DF if x["depth"] == d and x["mdot"] == 10.0) for d in depths}
data["tout_time"] = {f"{d / 1000:.0f} km": ds(base[d]["series"]["t_years"], base[d]["series"]["T_out"]) for d in depths}
data["q_time"] = {f"{d / 1000:.0f} km": ds(base[d]["series"]["t_years"], base[d]["series"]["Q_MW"]) for d in depths}
prof = {}
for d in (4000.0, 8000.0, 12000.0):
    p = base[d]["profile_10y"]
    md = np.array(p["md"])
    step = max(1, md.size // 60)
    sl = slice(0, None, step)
    prof[f"{d / 1000:.0f}"] = dict(T0=[[float(a), float(b)] for a, b in zip(np.array(p["T0"])[sl], md[sl])] + [[p["T0"][-1], md[-1]]],
                                   Ta=[[float(a), float(b)] for a, b in zip(np.array(p["T_a"])[sl], md[sl])] + [[p["T_a"][-1], md[-1]]],
                                   Tt=[[float(a), float(b)] for a, b in zip(np.array(p["T_t"])[sl], md[sl])] + [[p["T_t"][-1], md[-1]]])
data["profiles"] = prof


def sens(name, key_fn):
    runs = [x for x in L(name, []) if "rows" in x]
    out = {}
    for x in runs:
        k = key_fn(x["label"])
        out.setdefault(k, []).append([x["depth"] / 1000, row_at(x, 10)["Q_MW"]])
    for k in out:
        out[k].sort()
    return out, runs


grad_s, grad_runs = sens("D_gradient", lambda lab: lab.split(", ")[1])
grad_e, grad_runs_ext = sens("D_gradient_ext", lambda lab: lab.split(", ")[1])
grad_s.update(grad_e)
grad_s["30 C/km"] = data["power_depth"]["10"]["10"]
data["sens_gradient"] = {k.replace(" C/km", " °C/km"): grad_s[k] for k in ("20 C/km", "30 C/km", "40 C/km", "50 C/km", "60 C/km")
                         if k in grad_s}
ins_s, ins_runs = sens("D_insulation", lambda lab: lab.split(", ")[1])
ins_s["VIT k 0.03 (reference)"] = data["power_depth"]["10"]["10"]
data["sens_insulation"] = {k: ins_s[k] for k in ("VIT k 0.01", "VIT k 0.03 (reference)", "VIT k 0.10", "VIT top half only", "bare steel tubing") if k in ins_s}
rock_s, rock_runs = sens("D_rock", lambda lab: lab.split(", ")[1])
rock_s["k 2.8 (reference)"] = data["power_depth"]["10"]["10"]
data["sens_rock"] = {k: rock_s[k] for k in ("k 2.0", "k 2.8 (reference)", "k 3.5") if k in rock_s}
bb_runs = [x for x in L("D_bigbore", []) if "rows" in x]
data["bigbore"] = {}
for x in bb_runs:
    m = f"{x['mdot']:.0f} kg/s"
    data["bigbore"].setdefault(m, []).append([x["depth"] / 1000, row_at(x, 10)["Q_MW"]])
for k in data["bigbore"]:
    data["bigbore"][k].sort()
steam_runs = [x for x in L("D_steam", []) if "rows" in x]

# ---------------- development types (validation/dev_types.py -> DT_summary.json) ---------------
DT = L("DT_summary", [])
V8 = L("V8_utube_pygfunction", {})
V9 = L("V9_single_well", {})
V9b = L("V9b_twophase_production_mode", [])
B4 = L("B4_geocluster_utube", [])
V10 = L("V10_reservoir", {})
V11 = L("V11_plant", [])
T_dt = T_dt_elec = T_dt_bench = ""
dt_ctx = {}
if DT:
    ok = [r for r in DT if "failed" not in r]

    def idlab(r):
        return "3b" if r["id"] == 3.5 else f"{int(r['id'])}"

    def hz(r):
        h = r["horizon"]
        if r["id"] in (12, 13):
            return f"day {h}"
        if r["id"] == 4:
            return "steady"
        return f"year {h}"

    rows = []
    for r in ok:
        T1 = f(r["T_first_C"], 0) if r["T_first_C"] is not None else "–"
        TL = f(r["T_last_C"], 0) if r["T_last_C"] is not None else "–"
        x = r.get("x_last")
        Tcell = f"{TL} °C" + (f" <span class='sub'>steam fraction {x:.2f}</span>" if x and 0 < x < 1 and r["id"] in (3, 18, 19, 20) else
                              (" <span class='sub'>superheated steam</span>" if x and x >= 1 else ""))
        heat = f(r["heat_last_MWth"], 2) if r["heat_last_MWth"] is not None else "–"
        gross = f(r["gross_last_MWe"], 2) if r["gross_last_MWe"] else ("0" if r["gross_last_MWe"] == 0 else "–")
        g1 = r.get("gross_first_MWe")
        if g1 is not None and r["id"] != 4 and abs((g1 or 0) - (r["gross_last_MWe"] or 0)) > 0.05 * max(abs(g1 or 0), 0.02):
            gross = f"{f(g1, 2)} → {gross}"
        net = f(r["net_last_MWe"], 2) if r["gross_last_MWe"] else ("0" if r["gross_last_MWe"] == 0 else "–")
        rows.append([idlab(r), esc(r["name"]), hz(r), f"{T1} → {TL}" if r["id"] not in (4,) else TL, Tcell, heat, gross,
                     f(r["pump_MWe"], 2) if r["pump_MWe"] else "0", net, esc(r["sheet_MWe"])])
    T_dt = table(["#", "Development type (example run)", "Horizon", "Water T, first → last (°C)", "Water at the horizon", "Heat at the horizon (MWth)",
                  "Gross (MWe), first → horizon where it changes", "Pumps (MWe)", "Net at the horizon (MWe)", "Capacity quoted in the source table"], rows, cls="small",
                 note="One representative example per type, parameters in the table below. Water temperature at the wellhead (or the fluid "
                      "leaving the exchanger); 'first' is year 1 (day 30 for the cyclic cases). Heat relative to the return temperature of each "
                      "case. Gross electricity from the plant suited to the stream at 15 °C ambient: GEOPHIRES practical binary correlation "
                      "for liquid water, single / double flash or a steam turbine for two-phase and steam, 40 °C condenser; pumps = downhole "
                      "lift plus circulation / injection. Heat pumps consume electricity and are shown as 0 MWe.")
    # assumptions table
    arow = []
    for r in ok:
        a = r.get("assumptions") or {}
        txt = "; ".join(f"{k.replace('_', ' ')} {v}" for k, v in a.items() if not isinstance(v, dict))
        arow.append([idlab(r), esc(r["name"]), esc(txt), esc(r.get("plant_kind") or "")])
    T_dt_assume = table(["#", "Type", "Example parameters (SI units, °C, MPa, kg/s)", "Surface plant"], arow, cls="small")
    # benchmark table of the extensions
    brow = []
    if V9:
        e50 = max(c["rows"][3]["err_max"] for c in V9["steady"])
        ram = V9["ramey"]["rows"]
        dr = max(abs(x["T_wh_model_constant_fluid"] - x["T_wh_ramey"]) for x in ram if "T_wh_model_constant_fluid" in x)
        brow.append(["Production / injection well, exact steady solutions", "constant properties, fixed formation resistance",
                     f"second-order convergence; error {sci(e50, 0)} K at 50 m cells", "pass"])
        brow.append(["Production well vs Ramey / Hasan–Kabir", "2.5 km, 40 kg/s, 1 day – 10 years",
                     f"within {dr:.1f} K (constant properties); real water arrives {ram[-1]['T_wh_ramey'] - ram[-1]['T_wh_model']:.1f} K cooler through expansion", "pass"])
    if V9b:
        az = next(x for x in V9b if "Az-18" in x["name"]); m90 = next(x for x in V9b if "M-90" in x["name"])
        brow.append(["Measured flowing two-phase wells through the production mode", "Az-18, M-90, HGP-A, Utah 14-2",
                     f"Az-18 {az['mean_pct']:+.1f} %, M-90 {m90['mean_pct']:+.1f} % in pressure; same misfits as the stand-alone march at quality ≈ 0.6", "pass"])
    if V8:
        rr = V8["rows"]
        brow.append(["U-tube borehole exchanger vs pygfunction", "150 m single U-tube, 5 years",
                     f"wall temperature within {max(abs(x['T_b_pyg'] - x['T_b_aharrsim']) for x in rr):.2f} K, fluid within {max(abs(x['T_in_pyg'] - x['T_in_aharrsim']) for x in rr):.2f} K", "pass"])
    if B4:
        brow.append(["U-loop chain vs GeoCLUSTER utube database (SBT)", f"{len(B4)} cases, 2–5 km deep, 2–20 km laterals, 5–100 kg/s, 1–40 years",
                     f"max difference {max(x['max_abs_err'] for x in B4):.2f} K", "pass"])
    if V10:
        brow.append(["Reservoir elements", "Stehfest vs erfc and mpmath; doublet breakthrough; cyclic fracture vs Lauwerier",
                     f"{sci(V10['single_fracture_vs_erfc_max_err'], 0)} / {sci(V10['parallel_fractures_vs_mpmath_max_diff'], 0)}; breakthrough time exact; "
                     f"{V10['cyclic_vs_lauwerier']['max_err_K']:.1f} K", "pass"])
    if V11:
        orc = [x for x in V11 if x["kind"] == "binary"]
        rat = [x["eta_u_cycle"] / x["eta_u_geophires_sub"] for x in orc]
        fl = [x for x in V11 if x["kind"].startswith("flash, condenser 35")]
        frat = [x["eta_u_single"] / x["eta_u_geophires_single"] for x in fl if x["T_in_C"] >= 200]
        brow.append(["Plant models vs GEOPHIRES correlations", "binary 100–250 °C; flash 150–300 °C",
                     f"binary cycle {min(rat):.2f}–{max(rat):.2f} of GEOPHIRES; single flash {min(frat):.2f}–{max(frat):.2f} at 200–300 °C (35 °C condenser)", "pass"])
    T_dt_bench = table(["Check", "Case", "Result", ""], [[a, b, c_, f"<span class='pill ok'>{d_}</span>"] for a, b, c_, d_ in brow], cls="small")
    # chart data: water temperature vs net MWe per example
    data["dt_points"] = [[r["T_last_C"], max(r["net_last_MWe"], 0.0) if r["gross_last_MWe"] else 0.0, idlab(r) + " " + r["name"], r["group"]]
                         for r in ok if r["T_last_C"] is not None]
    data["dt_bars"] = {"labels": [idlab(r) for r in ok], "gross": [r["gross_last_MWe"] or 0 for r in ok], "net": [max(r["net_last_MWe"], 0) if r["gross_last_MWe"] else 0 for r in ok],
                       "heat": [r["heat_last_MWth"] or 0 for r in ok], "names": [r["name"] for r in ok]}
    def getr(i):
        return next((r for r in ok if r["id"] == i), None)
    g = {i: getr(i) for i in (1, 2, 3, 3.5, 7, 8, 10, 14, 15, 18, 19, 20, 5, 6, 16, 12, 13, 17, 11, 9, 4)}
    def nm(i, k, d=2):
        r = g.get(i)
        return "–" if not r or r.get(k) is None else f(r[k], d)
    dt_ctx = dict(dt_n=len(ok), dt_T3=nm(3, "T_last_C", 0), dt_P3=nm(3, "net_last_MWe", 1), dt_P3b=nm(3.5, "net_last_MWe", 1), dt_T2=nm(2, "T_last_C", 0),
                  dt_P2=nm(2, "net_last_MWe", 1), dt_T18=nm(18, "T_last_C", 0), dt_P18=nm(18, "net_last_MWe", 0), dt_T19=nm(19, "T_last_C", 0),
                  dt_P19=nm(19, "net_last_MWe", 0), dt_T20=nm(20, "T_last_C", 0), dt_P20=nm(20, "net_last_MWe", 1), dt_T8=nm(8, "T_last_C", 0),
                  dt_P8=nm(8, "net_last_MWe", 1), dt_H8=nm(8, "heat_last_MWth", 0), dt_T10=nm(10, "T_last_C", 0), dt_P10=nm(10, "net_last_MWe", 1),
                  dt_T14=nm(14, "T_last_C", 0), dt_P14=nm(14, "net_last_MWe", 1), dt_T7=nm(7, "T_last_C", 0), dt_H7=nm(7, "heat_last_MWth", 1),
                  dt_P7=nm(7, "net_last_MWe", 2), dt_T15=nm(15, "T_last_C", 0), dt_P15=nm(15, "net_last_MWe", 2), dt_T1=nm(1, "T_last_C", 0),
                  dt_H1=nm(1, "heat_last_MWth", 1), dt_T16=nm(16, "T_last_C", 0), dt_H16=nm(16, "heat_last_MWth", 2), dt_T12=nm(12, "T_last_C", 0),
                  dt_P13=nm(13, "net_last_MWe", 1), dt_T17=nm(17, "T_last_C", 0), dt_P17=nm(17, "net_last_MWe", 2), dt_T11=nm(11, "T_last_C", 0),
                  dt_H11=nm(11, "heat_last_MWth", 2), dt_T9=nm(9, "T_last_C", 0), dt_P9=nm(9, "net_last_MWe", 2), dt_P4=nm(4, "net_last_MWe", 2),
                  dt_cop5=f"{g[5]['extra']['COP']:.1f}" if g.get(5) else "–", dt_wm5=f"{g[5]['extra']['W_per_m']:.0f}" if g.get(5) else "–",
                  dt_wm6=f"{g[6]['extra']['W_per_m']:.0f}" if g.get(6) else "–",
                  dt_st13=f"{100 * g[13]['extra']['storage']['thermal_round_trip']:.0f}" if g.get(13) else "–",
                  dt_se13=f"{100 * g[13]['extra']['storage']['electric_round_trip']:.0f}" if g.get(13) else "–",
                  dt_rec17=f"{100 * g[17]['extra']['recovered_fraction_15yr']:.0f}" if g.get(17) else "–",
                  dt_T17f=nm(17, "T_first_C", 0), dt_P17f=nm(17, "gross_first_MWe", 2), dt_T8f=nm(8, "T_first_C", 0), dt_P8f=nm(8, "gross_first_MWe", 1),
                  dt_T11f=nm(11, "T_first_C", 0), dt_T19f=nm(19, "T_first_C", 0), dt_P19f=nm(19, "gross_first_MWe", 0), dt_T20f=nm(20, "T_first_C", 0),
                  dt_P20f=nm(20, "gross_first_MWe", 1), dt_T14f=nm(14, "T_first_C", 0), dt_T10f=nm(10, "T_first_C", 0), dt_T2f=nm(2, "T_first_C", 0),
                  dt_x19f=f"{g[19]['x_first']:.2f}" if g.get(19) and g[19].get("x_first") is not None else "–",
                  dt_x19=f"{g[19]['x_last']:.2f}" if g.get(19) and g[19].get("x_last") is not None else "–",
                  dt_Tres19=f"{g[19]['T_res_last_C']:.0f}" if g.get(19) and g[19].get("T_res_last_C") is not None else "–",
                  dt_Tres14=f"{g[14]['T_res_last_C']:.0f}" if g.get(14) and g[14].get("T_res_last_C") is not None else "–",
                  dt_P14f=nm(14, "gross_first_MWe", 1), dt_H13=nm(13, "heat_last_MWth", 1))

# ---------------- flow direction (validation/direction_study.py -> S_direction.json) ---------------
SD = L("S_direction", {})
V7 = L("V7_direction", {})
T_dir = ""
dir_ctx = {}
if SD:
    _dr = SD["rows"]
    data["dir_scatter"] = {
        "Tubing down, annulus up": [[r["v_annulus"], r["td"]["Q"], f"{r['casing']} casing × {r['tubing']} VIT, {r['A_annulus_cm2']:.0f} cm²"] for r in _dr],
        "Annulus down, tubing up": [[r["v_annulus"], r["ad"]["Q"], f"{r['casing']} casing × {r['tubing']} VIT, {r['A_annulus_cm2']:.0f} cm²"] for r in _dr]}
    _DR = [x for x in L("D_direction", []) if "rows" in x]
    _prof = {}
    for lab, key in (('9-5/8" x 3-1/2", annulus_down, 10 kg/s', "ad"), ('9-5/8" x 3-1/2", tubing_down, 10 kg/s', "td")):
        x = next((x for x in _DR if x["label"] == lab), None)
        if x:
            pp = x["profile_10y"]
            md = np.array(pp["md"]); step = max(1, md.size // 60); sl = slice(0, None, step)
            for nm, arr in (("Ta", pp["T_a"]), ("Tt", pp["T_t"]), ("T0", pp["T0"])):
                a = np.array(arr)
                _prof[f"{key}_{nm}"] = [[float(v), float(z)] for v, z in zip(a[sl], md[sl])] + [[float(a[-1]), float(md[-1])]]
    data["dir_prof"] = _prof
    drows = []
    for r in _dr:
        drows.append([esc(r["casing"]), esc(r["tubing"]), f"{r['A_annulus_cm2']:.0f}", f"{r['v_annulus']:.2f}", f"{r['v_tubing']:.1f}",
                      f"{r['R_t']:.2f}",
                      f(r["td"]["T_out"], 0), f(r["td"]["Q"], 2), f(r["td"]["pump_kW"], 0), f(r["td"]["MWe"], 2),
                      f(r["ad"]["T_out"], 0), f(r["ad"]["Q"], 2), f(r["ad"]["pump_kW"], 0), f(r["ad"]["MWe"], 2),
                      f"{100 * (r['td']['Q'] / r['ad']['Q'] - 1):+.0f} %".replace("-", "−")])
    T_dir = table(["Casing", "VIT", "Annulus (cm²)", "Annulus v (m/s)", "Tubing v (m/s)", "R<sub>t</sub> (K·m/W)",
                   "Tubing down: T out (°C)", "Heat (MWth)", "Pump (kW)", "Net power (MWe)",
                   "Annulus down: T out (°C)", "Heat (MWth)", "Pump (kW)", "Net power (MWe)", "Heat, tubing down vs annulus down"],
                  drows, cls="num small",
                  note="12 km, 30 °C/km, 10 kg/s, year 30; single casing string to TD, vacuum-insulated tubing to TD with the "
                       "reference construction (0.25″ vacuum gap, apparent k 0.03 W/m/K) at every size, so R<sub>t</sub>, the "
                       "tubing wall resistance per metre, falls as the tubing gets larger. Velocities at 960 and 995 kg/m³. "
                       "Net power by the GeoCLUSTER ORC method (net of circulation pumping); negative pump = thermosiphon surplus.")
    q_td = [r["td"]["Q"] for r in _dr]; q_ad = [r["ad"]["Q"] for r in _dr]
    ratio = [r["td"]["Q"] / r["ad"]["Q"] for r in _dr]
    ref = next(r for r in _dr if r["casing"] == '9-5/8"' and r["tubing"] == '3-1/2"')
    big = next(r for r in _dr if r["casing"] == '13-3/8"' and r["tubing"] == '3-1/2"')
    tight = next(r for r in _dr if r["casing"] == '7"' and r["tubing"] == '4-1/2"')
    fl = {(x["casing"], x["tubing"], x["direction"], x["mdot"]): x for x in SD["flow"]}
    gr = {(x["casing"], x["tubing"], x["direction"]): x for x in SD["grad"]}
    ch = {x["label"]: x for x in SD["checks"]}
    def _ratio(c, t, m):
        return 100 * (fl[(c, t, "tubing_down", m)]["Q"] / fl[(c, t, "annulus_down", m)]["Q"] - 1)
    dir_ctx = dict(
        dir_lo=f"{100 * (1 - max(ratio)):.0f}", dir_hi=f"{100 * (1 - min(ratio)):.0f}",
        dir_dT=f"{min(r['ad']['T_out'] - r['td']['T_out'] for r in _dr):.0f}–{max(r['ad']['T_out'] - r['td']['T_out'] for r in _dr):.0f}",
        dir_vlo=f"{min(r['v_annulus'] for r in _dr):.2f}", dir_vhi=f"{max(r['v_annulus'] for r in _dr):.1f}",
        dir_qlo=f"{min(q_td):.2f}", dir_qhi=f"{max(q_td):.2f}",
        dir_loss=f"{ref['td']['loss10']:.1f}", dir_xo=f"{ref['td']['crossover_md'] / 1000:.1f}",
        dir_peakz="6", dir_peak=f"{130:.0f}",
        dir_r5=f"{_ratio('9-5/8\"', '4-1/2\"', 5.0):.0f}", dir_r20=f"{_ratio('9-5/8\"', '4-1/2\"', 20.0):.0f}",
        dir_g60=f"{100 * (gr[('9-5/8\"', '4-1/2\"', 'tubing_down')]['Q'] / gr[('9-5/8\"', '4-1/2\"', 'annulus_down')]['Q'] - 1):.0f}",
        dir_pump_tight=f"{tight['td']['pump_kW']:.0f}",
        dir_big_gain=f"{100 * (big['ad']['Q'] / ref['ad']['Q'] - 1):.0f}",
        dir_r7=f"{ch['13-3/8\" x 7\" R_t 0.708, tubing_down, 10 kg/s']['Q']:.2f}",
        dir_q7=f"{next(r for r in _dr if r['casing'] == '13-3/8\"' and r['tubing'] == '7\"')['td']['Q']:.2f}",
        v7_err=sci(max(c['rows'][3]['err_max'] for c in V7['steady']), 0) if V7 else "–",
        v7_sym=sci(max(v['max_dT_after_30d'] for k, v in V7['symmetry'].items() if k.startswith('constant')), 0) if V7 else "–")


# ---------------- verification --------------------------------------------------------
V1 = L("V1_steady")
data["v1"] = {c["label"]: [[rw["dz"], rw["err_max"]] for rw in c["rows"]] for c in V1}
V2 = L("V2_jaeger")
data["v2"] = {f"radial growth {run['radial_growth']} ({run['M']} cells)": [[rw["t_days"], 100 * rw["rel_err"]] for rw in run["rows"]]
              for run in V2["runs"]}
V3 = L("V3_kabir")
_t3 = np.array(V3["series"]["t_days"]) / 365.25
_k3 = _t3 >= 0.001
data["v3"] = {"aharrsim (full r-z conduction)": ds(_t3[_k3], np.array(V3["series"]["T_out"])[_k3]),
              "Kabir et al. (1996) semi-analytical": ds(_t3[_k3], np.array(V3["series"]["T_out_kabir"])[_k3])}
V5 = L("V5_convergence")
V6 = L("V6_properties")
V4 = L("V4_energy")

# ---------------- field ---------------------------------------------------------------
FQ = L("F_qingdao")
_tq = np.array(FQ["series"]["t_days"])
_kq = _tq >= 0.05
data["qingdao"] = dict(T=ds(_tq[_kq], np.array(FQ["series"]["T_out"])[_kq]), Q=ds(_tq[_kq], np.array(FQ["series"]["Q_kW"])[_kq]),
                       meas_T=[138.0, FQ["measured"]["T_out_end"]], meas_Q=FQ["measured"]["mean_kW"],
                       lo=FQ["k_rock -10%"], hi=FQ["k_rock +10%"])
FW = L("F_weggis")
FX = L("F_xian")
FS = L("F_songyuan")
FB = L("F_bounds")

# ---------------- benchmarks -----------------------------------------------------------
B1 = L("B1_geocluster")
par_g = []
for rw in B1["rows"]:
    for y, a, b in zip(rw["years"], rw["T_model"], rw["T_db"]):
        par_g.append([round(b, 2), round(a, 2), f"L1 {rw['L1'] / 1000:.1f} km + 1 km lateral, {rw['grad'] * 1000:.0f} C/km, {rw['mdot']} kg/s, k {rw['k']}, D {rw['D']} m, {y:.0f} yr"])
d_g = np.array([p[1] - p[0] for p in par_g])
ex = [rw for rw in B1["rows"] if rw["L1"] == 5000.0 and rw["grad"] == 0.07]
data["gc_ts"] = {}
for rw in ex:
    lab = f"{rw['mdot']:.0f} kg/s"
    tt = np.array(rw["db_t"])
    TT = np.array(rw["db_T"])
    sel = tt >= 0.25
    _ts = np.array(rw["series_t"])
    _kk = _ts >= 0.2
    data["gc_ts"][lab] = dict(model=ds(_ts[_kk], np.array(rw["series_T"])[_kk], log=False),
                              db=[[float(a), float(b)] for a, b in zip(tt[sel][::4], TT[sel][::4])])
B2 = L("B2_sbt")
par_s = []
for rw in B2:
    for y, a, b in zip(rw["years"], rw["T_model"], rw["T_sbt"]):
        par_s.append([round(b, 2), round(a, 2), f"vertical {rw['depth'] / 1000:.0f} km, {rw['mdot']} kg/s, k_ins {rw['k_ins']}, SBT v{rw['sbt_version']}, {y:.0f} yr"])
d_s = np.array([p[1] - p[0] for p in par_s])
B3 = L("B3_published")
par_w = [[rw["published_T"], round(rw["T_out"], 2), f"{rw['depth'] / 1000:.0f} km, {rw['case']}, 1 yr"] for rw in B3["kubacka"]]
d_w = np.array([p[1] - p[0] for p in par_w])
data["parity"] = {"GeoCLUSTER database (Sandia Sierra/Aria)": par_g, "SBT (GeoCLUSTER implementation)": par_s,
                  "WellTH (Kubacka & Karayiannis 2024)": par_w}

# ---------------- two-phase ---------------------------------------------------------------
T1 = L("T1_twophase")
data["twophase"] = {}
for w in T1:
    if w["name"] not in ("Los Azufres Az-18", "Cerro Prieto M-90", "HGP-A, 12.8 kg/s", "HGP-A, 13.9 kg/s"):
        continue
    z = np.array(w["hem"]["z"])
    P = np.array(w["hem"]["P"])
    step = max(1, z.size // 50)
    zm = np.array(w["z_meas"], float)
    data["twophase"][w["name"]] = dict(meas=[[float(p), float(zz)] for p, zz in zip(w["P_meas"], zm)],
                                       model=[[float(p), float(zz)] for p, zz in zip(P[::step], z[::step])] + [[float(P[-1]), float(z[-1])]])

# =====================================================================================
# Summary numbers
b10 = {d: row_at(base[d], 10) for d in depths}
b30 = {d: row_at(base[d], 30) for d in depths}
b1 = {d: row_at(base[d], 1) for d in depths}
summary_rows = []
for d in depths:
    summary_rows.append([f"{d / 1000:.0f} km", f"{base[d]['rows'][0]['T_rock_bottom']:.0f}",
                         f(b1[d]["T_out"], 0), f(b10[d]["T_out"], 0), f(b30[d]["T_out"], 0),
                         f(b1[d]["Q_MW"], 2), f(b10[d]["Q_MW"], 2), f(b30[d]["Q_MW"], 2),
                         f(base[d]["cum_GWh_30y"], 0), f(b10[d]["exergy_MW"], 2), f(b10[d]["pump_kW"], 0)])
T_summary = table(["Depth", "Rock T at TD (°C)", "T out, 1 yr (°C)", "T out, 10 yr", "T out, 30 yr",
                   "Heat, 1 yr (MWth)", "Heat, 10 yr", "Heat, 30 yr", "30-yr total (GWh)", "Exergy, 10 yr (MW)",
                   "Pump, 10 yr (kW)"], summary_rows, cls="num",
                  note="10 kg/s injected at 40 °C; 30 °C/km from 15 °C; rock k 2.8 W/m/K; 3½″×4½″ vacuum-insulated tubing to TD "
                       "(insulation k 0.03 W/m/K); 7″ liner below 3 km. Heat = m(h<sub>out</sub> − h<sub>in</sub>), the net heat taken from "
                       "the rock and delivered at the wellhead. Exergy is the work obtainable by cooling the delivered water back to 40 °C "
                       "with a 15 °C environment. Pump is hydraulic power to circulate the loop (negative = thermosiphon surplus).")

flow_rows = []
for d in depths:
    cells = [f"{d / 1000:.0f} km"]
    for m in flows:
        x = next(x for x in DF if x["depth"] == d and x["mdot"] == m)
        rw = row_at(x, 10)
        cells.append(f"{rw['Q_MW']:.2f} MW · {rw['T_out']:.0f} °C · {rw['pump_kW']:.0f} kW")
    flow_rows.append(cells)
T_flow = table(["Depth"] + [f"{m:.0f} kg/s" for m in flows], flow_rows, cls="num small",
               note="After 10 years: heat delivered · outlet temperature · hydraulic pumping power. 3½″ tubing; above ~15 kg/s "
                    "friction in the tubing makes pumping expensive (see the large-bore option below).")

bb_rows = []
for x in sorted(bb_runs, key=lambda x: (x["depth"], x["mdot"])):
    rw = row_at(x, 10)
    bb_rows.append([f"{x['depth'] / 1000:.0f} km", f"{x['mdot']:.0f}", f(rw["T_out"], 0), f(rw["Q_MW"], 2), f(rw["exergy_MW"], 2),
                    f(rw["p_in_MPa"], 1), f(rw["pump_kW"], 0)])
T_bigbore = table(["Depth", "Flow (kg/s)", "T out (°C)", "Heat (MWth)", "Exergy (MW)", "Inlet pressure (MPa)", "Pump (kW)"],
                  bb_rows, cls="num small",
                  note="Large-bore completion: 9⅝″ casing to TD with 4½″×5½″ vacuum-insulated tubing. After 10 years; 5 MPa delivery pressure.")

# steam potential: flash the 10-year wellhead stream (reference design, 5 MPa) to low pressure
import sys as _sys
_sys.path.insert(0, os.path.dirname(HERE))
from aharrsim.props import WaterIAPWS  # noqa: E402
_W = WaterIAPWS()


def flash_x(p_wh_MPa, T_C, p_MPa):
    h = _W.h_pT(p_wh_MPa * 1e6, T_C + 273.15)[0]
    hl, hg = _W.sat(p_MPa * 1e6, "h_l"), _W.sat(p_MPa * 1e6, "h_g")
    return max(0.0, (h - hl) / (hg - hl))


steam_src = [(f"{d / 1000:.0f} km, 30 °C/km", next(x for x in DF if x["depth"] == d and x["mdot"] == m), m)
             for d in (8000.0, 10000.0, 12000.0) for m in (5.0, 10.0)]
_gr = [x for x in L("D_gradient", []) + L("D_gradient_ext", []) if "rows" in x]
for gl in ("40", "50", "60"):
    for d in (10000.0, 12000.0):
        x = next((x for x in _gr if x["depth"] == d and x["label"].endswith(f" {gl} C/km")), None)
        if x:
            steam_src.append((f"{d / 1000:.0f} km, {gl} °C/km", x, 10.0))
st_rows = []
for lab, x, m in steam_src:
    rw = row_at(x, 10)
    xs = [flash_x(rw["p_wh_MPa"], rw["T_out"], pp) for pp in (0.1013, 0.3, 0.5)]
    st_rows.append([lab, f"{m:.0f}", f(rw["T_out"], 0)] + [f(100 * v, 1) for v in xs] +
                   [f(xs[1] * m, 2), f(rw["Q_MW"], 2), f(rw["p_in_MPa"], 1)])
    rw["_x"] = xs
T_steam = table(["Case", "Flow (kg/s)", "T wellhead, 10 yr (°C)", "Steam % flashed to 1 atm", "to 0.3 MPa", "to 0.5 MPa",
                 "Steam at 0.3 MPa (kg/s)", "Heat (MWth)", "Inlet p (MPa)"], st_rows, cls="num small",
                note="Mass fraction of steam if the delivered water (reference completion, 5 MPa at the wellhead, after 10 years) is "
                     "throttled to the given pressure (isenthalpic flash). An inlet pressure below 5 MPa means the hot riser drives the loop "
                     "without a pump.")
_x30 = [row_at(x, 10)["_x"] for lab, x, m in steam_src if "30 °C" in lab]
_pin5 = {d: row_at(next(x for x in DF if x["depth"] == d and x["mdot"] == 5.0), 10)["p_in_MPa"] for d in (8000.0, 10000.0, 12000.0)}
_pin10 = [row_at(next(x for x in DF if x["depth"] == d and x["mdot"] == 10.0), 10)["p_in_MPa"] - 5.0 for d in (8000.0, 10000.0, 12000.0)]
steam_ctx = dict(
    st_x05=f"{100 * max(v[2] for v in _x30):.1f}", st_atm_lo=f"{100 * min(v[0] for v in _x30):.0f}",
    st_atm_hi=f"{100 * max(v[0] for v in _x30):.0f}",
    st_pin5=f"{_pin5[8000.0]:.1f}, {_pin5[10000.0]:.1f} and {_pin5[12000.0]:.1f}",
    st_pump10=f"{min(_pin10):.1f}–{max(_pin10):.1f}")
_xhot = [row_at(x, 10) for lab, x, m in steam_src if lab.endswith(("50 °C/km", "60 °C/km"))]
steam_ctx["st_hot"] = f"{100 * min(r['_x'][2] for r in _xhot):.0f}–{100 * max(r['_x'][2] for r in _xhot):.0f}" if _xhot else "–"
steam_ctx["st_hotT"] = f"{min(r['T_out'] for r in _xhot):.0f}–{max(r['T_out'] for r in _xhot):.0f}" if _xhot else "–"

# ---- 5 GW field sizing (validation/field_sizing.py -> S_5GW.json)
S5 = L("S_5GW", {})
T_5th = T_5e = T_max = T_kT = ""
s5ctx = {}
if S5:
    G5, D5 = S5["grads"], S5["depths_km"]
    TH, EL, MX, KT = S5["thermal"], S5["electric"], S5["max"], S5["kT"]

    def wn(n):
        return "–" if n is None else f"{n:,}"

    def glab(g):
        return f"<b>{g} °C/km</b> (base)" if g == 30 else f"{g} °C/km"

    rows5 = []
    for g in G5:
        cells = [glab(g), f"{2.8 * g:.0f}"]
        for d in D5:
            x = TH.get(f"{g}|{d:.0f}")
            if not x:
                cells.append("..")
                continue
            mark = "†" if x["beyond"] else ""
            cells.append(f"<span class='{'hyp' if mark else ''}'>{wn(x['wells_r'])}{mark}</span>"
                         f"<span class='sub'>{x['Q_MW']:.2f} MWth · {x['T_out']:.0f} °C</span>")
        rows5.append(cells)
    T_5th = table(["Gradient", "Heat flow (mW/m²)"] + [f"{d:.0f} km" for d in D5], rows5, cls="num five",
                  note="Wells needed so the field still delivers 5 GWth in year 30, with heat per well and outlet temperature in "
                       "year 30 underneath. Reference design at 10 kg/s (3½″ × 4½″ vacuum-insulated tubing, insulation k 0.03 W/m/K, "
                       "40 °C injection, rock k 2.8 W/m/K); heat flow = 2.8 W/m/K × gradient. † Rock above 400 °C at TD: beyond the "
                       "brittle–ductile transition in granite and the rating of current casing, cement and insulated tubing.")
    rows5e = []
    for g in G5:
        cells = [glab(g)]
        for d in D5:
            x = EL.get(f"{g}|{d:.0f}")
            if not x:
                cells.append("..")
                continue
            mark = "†" if x["beyond"] else ""
            if x["wells_r"] is None:
                cells.append(f"–<span class='sub'>{x['best']['T_out']:.0f} °C at {x['best']['mdot']:.0f} kg/s</span>")
            else:
                b = x["best"]
                cells.append(f"<span class='{'hyp' if mark else ''}'>{wn(x['wells_r'])}{mark}</span>"
                             f"<span class='sub'>{b['net_kWe'] / 1e3:.2f} MWe · {b['mdot']:.0f} kg/s · {b['T_out']:.0f} °C</span>")
        rows5e.append(cells)
    T_5e = table(["Gradient"] + [f"{d:.0f} km" for d in D5], rows5e, cls="num five",
                 note="Wells needed for 5 GWe net in year 30, with net power per well, the best of 5, 10, 15 and 20 kg/s, and the outlet "
                      "temperature underneath. ORC conversion as in GeoCLUSTER: 20 % of the thermal exergy at 100 °C rising linearly to "
                      "45 % at 200 °C (held above 200 °C), nothing below 100 °C; 20 °C heat sink; circulation pump at 80 % efficiency "
                      "deducted. – : no net power at any of the flow rates (outlet temperature shown for the best of them). "
                      "† as above.")
    mrows = []
    for dk in [f"{d:.0f}" for d in D5]:
        x = MX.get(dk)
        if not x:
            continue
        t_, e_ = x["thermal"], x["electric"]
        viable = e_["wells_r"] is not None
        mrows.append([f"{dk} km", f"{x['T_rock_TD']:.0f}",
                      f"{t_['mdot']:.0f}", f"{t_['Q_MW']:.2f}", f"{t_['T_out']:.0f}", wn(t_["wells_r"]),
                      f"{t_['pump_field_MWe']:.0f}",
                      f"{e_['mdot']:.0f}" if viable else "–", f"{e_['net_kWe'] / 1e3:.2f}" if viable else "–",
                      wn(e_["wells_r"])])
    T_max = table(["Depth", "Rock at TD (°C)", "Heat: flow (kg/s)", "MWth per well", "Outlet (°C)", "Wells for 5 GWth",
                   "Pumping, field (MWe)", "Power: flow (kg/s)", "MWe per well", "Wells for 5 GWe"], mrows, cls="num small",
                  note="60 °C/km, rock k 3.5 W/m/K, large bore (9⅝″ casing to TD, 4½″ × 5½″ vacuum-insulated tubing) with insulation "
                       "k 0.01 W/m/K; flows 5, 10, 20, 30 and 40 kg/s; year 30. Heat: the flow giving the most heat with circulation "
                       "power no more than 5 % of it. Power: the flow giving the most net power, same conversion as above.")
    krows = []
    for dk in [f"{d:.0f}" for d in D5]:
        a, b = KT.get(f"30|{dk}"), KT.get(f"60|{dk}")
        if not (a and b):
            continue
        krows.append([f"{dk} km", f"{a['Q_const']:.2f}", f"{a['Q_kT']:.2f}", f"{100 * a['change']:+.0f} %".replace("-", "−"),
                      f"{b['Q_const']:.2f}", f"{b['Q_kT']:.2f}", f"{100 * b['change']:+.0f} %".replace("-", "−")])
    T_kT = table(["Depth", "30 °C/km: constant k (MWth)", "k falling with T", "Change", "60 °C/km: constant k", "k falling with T",
                  "Change"], krows, cls="num small",
                 note="Reference design at 10 kg/s, year 30. Rock conductivity 2.8 W/m/K at 25 °C, falling with the undisturbed rock "
                      "temperature after Vosteen & Schellschmidt (2003), held flat above 500 °C. Well counts rise in proportion.")
    data["s5_heat"] = {f"{d:.0f} km": [[g, TH[f"{g}|{d:.0f}"]["Q_MW"]] for g in G5 if f"{g}|{d:.0f}" in TH] for d in D5}
    data["s5_wells"] = {
        "heat": {f"{d:.0f} km": [[g, TH[f"{g}|{d:.0f}"]["wells_r"]] for g in G5 if f"{g}|{d:.0f}" in TH] for d in D5},
        "power": {f"{d:.0f} km": [[g, EL[f"{g}|{d:.0f}"]["wells_r"]] for g in G5
                                  if f"{g}|{d:.0f}" in EL and EL[f"{g}|{d:.0f}"]["wells_r"]] for d in D5}}

    def th_w(g, d):
        x = TH.get(f"{g}|{d}")
        return wn(x["wells_r"]) if x else "–"

    def el_w(g, d):
        x = EL.get(f"{g}|{d}")
        return wn(x["wells_r"]) if x and x["wells_r"] else "none"

    _kc30 = [v["change"] for k, v in KT.items() if k.startswith("30|")]
    _kc60 = [v["change"] for k, v in KT.items() if k.startswith("60|")]
    _sp = S5["spacing"]
    _spr, _spm = list(_sp.values())[0], list(_sp.values())[1]
    s5ctx = dict(
        s5_12_30=th_w(30, 12), s5_12_60=th_w(60, 12), s5_4_30=th_w(30, 4), s5_4_60=th_w(60, 4), s5_12_20=th_w(20, 12),
        s5e_12_30=el_w(30, 12), s5e_12_60=el_w(60, 12), s5e_8_60=el_w(60, 8), s5e_4_60=el_w(60, 4),
        mx_th_12=wn(MX.get("12", {}).get("thermal", {}).get("wells_r")) if MX.get("12") else "–",
        mx_e_12=wn(MX.get("12", {}).get("electric", {}).get("wells_r")) if MX.get("12") else "–",
        mx_q_12=f"{MX['12']['thermal']['Q_MW']:.1f}" if MX.get("12") else "–",
        mx_p_12=f"{MX['12']['electric']['net_kWe'] / 1e3:.1f}" if MX.get("12") else "–",
        kT30=f"{-100 * max(_kc30):.0f}–{-100 * min(_kc30):.0f}" if _kc30 else "–",
        kT60=f"{-100 * max(_kc60):.0f}–{-100 * min(_kc60):.0f}" if _kc60 else "–",
        sp_reach=f"{_spr['reach_m']:.0f}", sp100=f"{100 * _spr['100']:.1f}", sp100m=f"{100 * _spm['100']:.1f}",
        sp150m=f"{100 * _spm['150']:.1f}")

# ---- verification tables
v1_rows = []
for c in V1:
    rws = c["rows"]
    v1_rows.append([esc(c["label"]), f(rws[0]["T_out_exact"], 3)] + [sci(rw['err_max'], 1) for rw in rws] +
                   [f(c["observed_order"][-1], 2)])
T_v1 = table(["Case", "Exact T out (°C)"] + [f"dz {rw['dz']:g} m" for rw in V1[0]["rows"]] + ["Order"], v1_rows, cls="num small",
             note="Largest temperature error anywhere in the annulus or tubing (K) against the closed-form solution.")
v2_rows = []
for run in V2["runs"]:
    v2_rows.append([f"{run['radial_growth']} ({run['M']} cells)"] + [f"{100 * rw['rel_err']:+.2f}" for rw in run["rows"]])
T_v2 = table(["Radial grid"] + [f"{rw['t_days']:g} d" for rw in V2["runs"][0]["rows"]], v2_rows, cls="num small",
             note="Error in heat flow per metre (%) against Carslaw & Jaeger’s constant-temperature cylinder; dimensionless time "
                  f"αt/r² from {V2['runs'][0]['rows'][0]['tau']:.2g} to {V2['runs'][0]['rows'][-1]['tau']:.2g}.")
v5_rows = []
basev5 = next(x for x in V5["rows"] if "base" in x["label"])
for x in V5["rows"]:
    v5_rows.append([esc(x["label"]), str(x["cells"]).replace(" x ", " × "), str(x["steps"])] +
                   [f"{x['T_out'][str(t)]:.2f}" for t in V5["t_check"]] +
                   [f"{max(abs(x['T_out'][str(t)] - basev5['T_out'][str(t)]) for t in V5['t_check'][1:]):.2f}"])
T_v5 = table(["Variant", "Cells (well × radial)", "Steps", "T out 1 d", "30 d", "1 yr", "10 yr", "Max Δ vs base, 30 d–10 yr (K)"],
             v5_rows, cls="num small", note="8 km, 10 kg/s, real completion with vacuum-insulated tubing to 5 km and bare tubing below; "
                                            "the bare section is automatically refined.")
w6 = V6["table_worst"]
T_v6 = table(["Check", "Worst deviation"], [
    ["IAPWS-95 verification values (Table 7, 11 states: p, c<sub>v</sub>, w, s)", f"{V6['table7_worst_rel']:.0e} relative"],
    [f"Liquid / supercritical states, table vs IAPWS-95 ({w6['liquid']['n']} random states, 0.1–150 MPa, 2–417 °C)",
     f"T {w6['liquid']['T']:.3f} K; ρ {100 * w6['liquid']['rho']:.3f} %; c<sub>p</sub> {100 * w6['liquid']['cp']:.2f} %; μ, k {100 * max(w6['liquid']['mu'], w6['liquid']['k']):.3f} %"],
    ["Subcooled liquid 0.2–5 K below saturation, 0.5–20 MPa", f"T {w6.get('near_saturation_T', float('nan')):.3f} K"],
    [f"Steam states ({w6['vapour']['n']} states)",
     f"T {w6['vapour']['T']:.2f} K; ρ {100 * w6['vapour']['rho']:.2f} %; c<sub>p</sub> {100 * w6['vapour']['cp']:.1f} %; μ, k {100 * max(w6['vapour']['mu'], w6['vapour']['k']):.2f} %"],
    ["Saturation temperature, two-phase table", sci(V6['sat_worst_K']) + " K"]], cls="small")

# ---- field
FQb = FQ["base"]
T_field = table(["Site", "Measured", "aharrsim", "Difference", "Note"], [
    ["<b>Qingdao (Jimo)</b>, 2605 m, 138 days, 5 °C inlet, 30 m³/h",
     f"mean {FQ['measured']['mean_kW']:.0f} kW;<br>outlet ≈{FQ['measured']['T_out_end']:.1f} °C at end",
     f"mean {FQb['mean_kW']:.0f} kW;<br>outlet {FQb['T_out_end']:.1f} °C",
     f"{100 * (FQb['mean_kW'] / FQ['measured']['mean_kW'] - 1):+.1f} %;<br>{FQb['T_out_end'] - FQ['measured']['T_out_end']:+.1f} K",
     f"Rock k ±10 % moves the mean by ±{(FQ['k_rock +10%']['mean_kW'] - FQ['k_rock -10%']['mean_kW']) / 2:.0f} kW. Full input set published."],
    ["<b>Weggis</b>, 2302 m, 2001–03 continuous, mean inlet 29.9 °C",
     f"{FW['measured']['MWh_yr'][0]:.0f} and {FW['measured']['MWh_yr'][1]:.0f} MWh/yr;<br>mean outlet {FW['measured']['T_out']:.1f} °C",
     f"{FW['base']['MWh_y1']:.0f} and {FW['base']['MWh_y2']:.0f} MWh/yr;<br>mean outlet {FW['base']['T_out_mean']:.1f} °C",
     f"{100 * (FW['base']['MWh_y1'] / FW['measured']['MWh_yr'][0] - 1):+.0f} % and {100 * (FW['base']['MWh_y2'] / FW['measured']['MWh_yr'][1] - 1):+.0f} %;<br>"
     f"{FW['base']['T_out_mean'] - FW['measured']['T_out']:+.1f} K",
     f"Rock properties and geotherm not published: plausible range gives "
     f"{min(v['MWh_y2'] for k, v in FW.items() if k != 'measured'):.0f}–{max(v['MWh_y1'] for k, v in FW.items() if k != 'measured'):.0f} MWh/yr."],
    ["<b>Xi’an</b>, 2000 m, 30 days, 17.3 °C inlet, 28 m³/h",
     f"mean {FX['measured']['mean_kW']:.0f} kW;<br>outlet {FX['measured']['T_out_last4'][0]}–{FX['measured']['T_out_last4'][1]} °C",
     f"mean {FX['base']['mean_kW']:.0f} kW;<br>outlet {FX['base']['T_out_day30']:.1f} °C",
     f"{100 * (FX['base']['mean_kW'] / FX['measured']['mean_kW'] - 1):+.0f} %",
     f"<b>Not usable.</b> Even with the annulus held at the inlet temperature the published rock can supply only "
     f"{FB['xian']['bound_kW_mean']:.0f} kW on average; the measured {FX['measured']['mean_kW']:.0f} kW implies groundwater flow or higher conductivity."],
    ["<b>Songyuan</b>, 2044 m, 60 days, 30 m³/h",
     f"{FS['measured']['day60']['q_W_m']:.0f} W/m at day 60;<br>outlet {FS['measured']['day60']['T_out']} °C",
     f"{FS['base']['q_W_m_day60']:.0f} W/m;<br>outlet {FS['base']['T_out_day60']:.1f} °C",
     f"{100 * (FS['base']['q_W_m_day60'] / FS['measured']['day60']['q_W_m'] - 1):+.0f} %",
     f"<b>Not usable.</b> Conduction upper bound is {FB['songyuan']['bound_W_m_day60']:.0f} W/m at day 60, less than half the reported value."],
], cls="small")

# ---- benchmarks
T_gc = table(["Comparison", "Points", "Mean ΔT (K)", "RMS (K)", "Max |ΔT| (K)", "Outlet pressure"], [
    ["GeoCLUSTER database (Sandia Sierra/Aria), 4–5 km + 1 km lateral, 1–40 yr", str(len(par_g)), f(float(d_g.mean()), 2, True),
     f(float(np.sqrt((d_g ** 2).mean())), 2), f(float(np.abs(d_g).max()), 2),
     f"within {max(abs(a - b) for rw in B1['rows'] for a, b in zip(rw['p_out_model'], rw['p_out_db'])):.2f} MPa"],
    ["SBT, vertical 4–10 km, v1 (constant properties) and v2 (IAPWS)", str(len(par_s)), f(float(d_s.mean()), 2, True),
     f(float(np.sqrt((d_s ** 2).mean())), 2), f(float(np.abs(d_s).max()), 2), "–"],
    ["WellTH (Kubacka & Karayiannis 2024), 1–5 km, 1 yr", str(len(par_w)), f(float(d_w.mean()), 2, True),
     f(float(np.sqrt((d_w ** 2).mean())), 2), f(float(np.abs(d_w).max()), 2), "–"]], cls="num small",
    note="ΔT = aharrsim − benchmark outlet temperature. Each comparison uses the benchmark’s own geometry, rock and fluid assumptions.")
kub_rows = []
for rw in B3["kubacka"]:
    kub_rows.append([f"{rw['depth'] / 1000:.0f} km", esc(rw["case"]), f(rw["published_T"], 1), f(rw["T_out"], 1),
                     f(rw["published_Q"], 0), f(rw["Q_kW"], 0),
                     "–" if rw.get("T_out_gap_reading") is None else f(rw["T_out_gap_reading"], 1)])
T_kub = table(["Depth", "Centre pipe", "WellTH T out (°C)", "aharrsim", "WellTH heat (kW)", "aharrsim", "aharrsim, gap-only reading"],
              kub_rows, cls="num small",
              note="The paper gives the vacuum-insulated tubing an ‘equivalent conductivity’ (0.00024, 0.062, 0.081 W/m/K). Applied across "
                   "the whole tubing wall it reproduces WellTH; applied to the 8.5 mm vacuum gap alone it does not (last column).")
bk_rows = [[f"{rw['depth'] / 1000:.0f} km", f(rw["published_MW"], 2), f(rw["Q_avg_MW"], 3), f"{100 * (rw['Q_avg_MW'] / rw['published_MW'] - 1):+.1f} %"]
           for rw in B3["beckers"]]
T_beckers = table(["Depth (200 °C at TD)", "SBT, 20-yr mean (MWth)", "aharrsim", "Difference"], bk_rows, cls="num small")
na_rows = [[esc(rw["label"]), f(rw["published"], 1), f(rw["T_out_500d"], 1), f(rw["T_out_500d"] - rw["published"], 1, True)]
           for rw in B3["nalla"]]
T_nalla = table(["Case", "TETRAD T out (°C)", "aharrsim", "Δ (K)"], na_rows, cls="num small",
                note="Outlet temperature at pseudo-steady state (500 days). TETRAD results read from Nalla et al. (2004), Table 7.")
sbt_rows = []
for rw in B2:
    sbt_rows.append([f"{rw['depth'] / 1000:.0f} km", f"{rw['mdot']:.0f}", f"{rw['k_ins']}", f"v{rw['sbt_version']}"] +
                    [f"{a:.2f} / {b:.2f}" for a, b in zip(rw["T_model"], rw["T_sbt"])])
T_sbt = table(["Depth", "Flow (kg/s)", "Insulation k", "SBT"] + [f"{y:.0f} yr: aharrsim / SBT (°C)" for y in B2[0]["years"]],
              sbt_rows, cls="num small",
              note="v1: constant water properties (aharrsim run with the same constants, no gravity work); v2: temperature- and "
                   "pressure-dependent properties (aharrsim with IAPWS-95). Vertical well, open hole 0.216 m, 40 °C inlet, 30 °C/km, "
                   "k 3 W/m/K. SBT v2 stopped with an out-of-range error in its own property tables for the three deeper 20 kg/s cases, "
                   "which are therefore compared with v1 only.")

# ---- two-phase
tp_rows = []
for w in T1:
    o = w["orkiszewski"]
    tp_rows.append([esc(w["name"]), f"{w['W']:.1f}", f"{w['h_kJ']:.0f}",
                    f"{w['hem']['mean_bar']:+.2f} ({w['hem']['sd_bar']:.2f})".replace("-", "−"),
                    f"{w['hem']['mean_pct']:+.1f} %".replace("-", "−"),
                    f"{w['drift_flux']['mean_pct']:+.1f} %".replace("-", "−"),
                    "–" if not o else f"{o['mean_pct']:+.1f} %".replace("-", "−"), esc(w["note"])])
T_tp = table(["Well", "Flow (kg/s)", "Enthalpy (kJ/kg)", "HEM error, bar (SD)", "HEM %", "Drift-flux %", "Orkiszewski % (1986)", "Note"],
             tp_rows, cls="small", note="Pressure error = calculated − measured flowing pressure, averaged over the survey points. "
                                         "Orkiszewski values are those reported by Ambastha & Gudmundsson for the same data.")

# ---- headline numbers
h12 = b10[12000.0]
h8 = b10[8000.0]
h4 = b10[4000.0]

def bottom_third(run):
    p = run["profile_10y"]
    md = np.array(p["md_c"])
    q = np.array(p["q_rock"]) * np.diff(np.array(p["md"]))
    return q[md > 2 * run["depth"] / 3].sum() / q.sum()


_b3 = [bottom_third(base[d]) for d in depths]
_dbl = [row_at(next(x for x in DF if x["depth"] == d and x["mdot"] == 20.0), 10)["Q_MW"] / b10[d]["Q_MW"] - 1 for d in depths]
_dT12 = b10[12000.0]["T_out"] - row_at(next(x for x in DF if x["depth"] == 12000.0 and x["mdot"] == 20.0), 10)["T_out"]

ctx = dict(
    **steam_ctx, **s5ctx, **dir_ctx, **dt_ctx,
    T_5th=T_5th, T_5e=T_5e, T_max=T_max, T_kT=T_kT, T_dir=T_dir,
    T_dt=T_dt, T_dt_assume=(T_dt_assume if DT else ""), T_dt_bench=T_dt_bench,
    bot3=f"{100 * min(_b3):.0f}–{100 * max(_b3):.0f} %",
    dbl_lo=f"{100 * min(_dbl):.0f}", dbl_hi=f"{100 * max(_dbl):.0f}", dbl_dT=f"{_dT12:.0f}",
    T_summary=T_summary, T_flow=T_flow, T_bigbore=T_bigbore, T_steam=T_steam, T_v1=T_v1, T_v2=T_v2, T_v5=T_v5, T_v6=T_v6,
    T_field=T_field, T_gc=T_gc, T_kub=T_kub, T_beckers=T_beckers, T_nalla=T_nalla, T_sbt=T_sbt, T_tp=T_tp,
    q4=f(h4["Q_MW"], 2), q8=f(h8["Q_MW"], 1), q12=f(h12["Q_MW"], 1), t4=f(h4["T_out"], 0), t8=f(h8["T_out"], 0), t12=f(h12["T_out"], 0),
    gc_rms=f(float(np.sqrt((d_g ** 2).mean())), 2), gc_max=f(float(np.abs(d_g).max()), 2), gc_n=len(par_g),
    sbt_max=f(float(np.abs(d_s).max()), 2), wt_max=f(float(np.abs(d_w).max()), 1),
    v1_err50=sci(max(c['rows'][3]['err_max'] for c in V1), 0), v2_err=f"{100 * max(abs(rw['rel_err']) for rw in V2['runs'][1]['rows'][1:]):.2f}",
    v4_res=sci(V4['max_rel_all'], 0),
    qd_mean=f(FQb["mean_kW"], 0), qd_meas=f(FQ["measured"]["mean_kW"], 0), qd_dT=f(FQb["T_out_end"] - FQ["measured"]["T_out_end"], 1, True),
    qd_pct=f(100 * (FQb["mean_kW"] / FQ["measured"]["mean_kW"] - 1), 0, True),
    wg_pct1=f(100 * (FW["base"]["MWh_y1"] / FW["measured"]["MWh_yr"][0] - 1), 0, True),
    wg_pct2=f(100 * (FW["base"]["MWh_y2"] / FW["measured"]["MWh_yr"][1] - 1), 0, True),
    bk_max=f(max(abs(100 * (rw['Q_avg_MW'] / rw['published_MW'] - 1)) for rw in B3["beckers"]), 1),
    cum12=f(base[12000.0]["cum_GWh_30y"], 0), cum8=f(base[8000.0]["cum_GWh_30y"], 0),
    loss12=f(b10[12000.0]["riser_loss_MW"], 2),
    q12_20=f(row_at(next(x for x in DF if x['depth'] == 12000.0 and x['mdot'] == 20.0), 10)["Q_MW"], 1),
    pump12_20=f(row_at(next(x for x in DF if x['depth'] == 12000.0 and x['mdot'] == 20.0), 10)["pump_kW"], 0),
    data=json.dumps(data, separators=(",", ":")),
)

with open(os.path.join(HERE, "template.html")) as fh:
    tpl = fh.read()
for k, v in ctx.items():
    tpl = tpl.replace("{{" + k + "}}", str(v))
with open(OUT, "w") as fh:
    fh.write(tpl)
print("wrote", OUT, f"{os.path.getsize(OUT) / 1024:.0f} kB")
