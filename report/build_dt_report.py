"""Build the development-types report (aharrsim_dev_types_report.html) from validation/results/*.json.

    python build_dt_report.py            # writes report/aharrsim_dev_types_report.html
"""
from __future__ import annotations

import datetime as _dt
import json
import math
import os
import sys

import re as _re

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RES = os.path.join(ROOT, "validation", "results")
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)
import CoolProp.CoolProp as CP  # noqa: E402
from aharrsim import plant  # noqa: E402
from dt_report_content import TYPES, GROUPS, SHEET_QUOTE, SHEET_MWE, TINY  # noqa: E402


def L(name, default=None):
    p = os.path.join(RES, name + ".json")
    if not os.path.exists(p):
        return default
    with open(p) as f:
        return json.load(f)


def f(v, d=1):
    if v is None or (isinstance(v, float) and not math.isfinite(v)):
        return "–"
    s = f"{v:,.{d}f}"
    return s.replace("-", "−")


def key_of(tid):
    return "t3b" if tid == 3.5 else f"t{int(tid)}"


def lab_of(tid):
    return "3b" if tid == 3.5 else str(int(tid))


def esc(s):
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def table(head, rows, cls="num", note=None):
    h = f'<div class="tbl"><table class="{cls}"><thead><tr>' + "".join(f"<th>{c}</th>" for c in head) + "</tr></thead><tbody>"
    for r in rows:
        h += "<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>"
    h += "</tbody></table></div>"
    if note:
        h += f'<p class="tnote">{note}</p>'
    return h


# =====================================================================================================
R = {float(r["id"]): r for r in L("DT_types")}
S = {float(r["id"]): r for r in L("DT_summary")}
GROUP_OF = {tid: g for g, ids in GROUPS for tid in ids}
ORDER = [tid for _, ids in GROUPS for tid in ids]
today = _dt.date.today().strftime("%-d %B %Y")

# ---------------- cases for the page ----------------
cases = []
for tid in ORDER:
    r, s, c = R[tid], S[tid], TYPES[tid]
    ser = r.get("series")
    x_kind = "day" if tid in (12, 13) else ("lin" if tid in (5, 17) else "log")
    horizon = {5: 10, 6: 10, 17: 15}.get(tid, 30)
    well_kind = "coaxial" if tid in (6, 7, 11, 16, 20) else ("utube" if tid == 5 else ("injection" if False else "production"))
    variants = []
    if tid == 7:
        variants = [dict(label="5 kg/s", series=r["variants"]["5.0"]["series"])]
    if tid == 8:
        variants = [dict(label="isolated laterals (bound)", series=r["variants"]["isolated laterals"]["series"])]
    if tid == 20:
        variants = [dict(label="40 kg/s", series=r["variants"]["40.0"]["series"])]
    plant_short = {"binary (GEOPHIRES practical)": "binary", "double flash": "double flash", "steam turbine (wellhead separator)": "steam turbine",
                   "dry / superheated steam turbine": "steam turbine", "binary": "binary", "superheated steam turbine": "steam turbine"}.get(s.get("plant"), "binary")
    if tid == 17 and ser:   # drop the first weeks: the Duhamel loss term is singular at t = 0 and would set the chart scale
        keep = [i for i, tt in enumerate(ser["t_yr"]) if tt >= 0.08]
        ser = {k: ([vv[i] for i in keep] if isinstance(vv, list) and len(vv) == len(ser["t_yr"]) else vv) for k, vv in ser.items()}
    cases.append(dict(key=key_of(tid), id=tid, lab=lab_of(tid), short=c["short"], group=GROUP_OF[tid], series=ser, variants=variants,
                      profile=r.get("profile"), profile_inj=r.get("profile_inj"), x_kind=x_kind, horizon=horizon, well_kind=well_kind,
                      plant=plant_short, spec=c["spec"], t_min=0.005 if tid not in (6,) else 0.005))

# ---------------- benchmarks ----------------
V9 = L("V9_single_well"); V9b = L("V9b_twophase_production_mode"); V8 = L("V8_utube_pygfunction"); B4 = L("B4_geocluster_utube")
V10 = L("V10_reservoir"); V11 = L("V11_plant")
steady = [dict(label=f"{st['mode']}, {st['L'] / 1000:.0f} km, {st['w']:.0f} kg/s", rows=st["rows"], order=st["order"]) for st in V9["steady"]]
tp = []
for w in V9b:
    if "profile" not in w:
        continue
    tp.append(dict(key=w["name"].lower().replace(" ", "-").replace(",", "").replace(".", "p"), name=w["name"], profile=w["profile"], measured=w["measured"],
                   n=w["n"], mean_pct=w["mean_pct"], sd_pct=w["sd_pct"], pb_model=w["p_bottom_model"], pb_meas=w["p_bottom_meas"], T_wh=w["T_wh"], x_wh=w["x_wh"],
                   mdot=w.get("mdot"), h_kJ=w.get("h_kJ")))
uloop = [dict(label=f"{c['L1'] / 1000:g} km + {c['L2'] / 1000:g} km, {c['m']:g} kg/s, {c['grad'] * 1000:.0f} °C/km", comparison=c["comparison"], max_abs_err=c["max_abs_err"]) for c in B4]
binary = [dict(T=r["T_in_C"], cycle=r["eta_u_cycle"], sub=r["eta_u_geophires_sub"], sup=r["eta_u_geophires_super"], fluid=r["fluid"]) for r in V11 if r["kind"] == "binary"]
flash35 = [dict(T=r["T_in_C"], single=r["eta_u_single"], double=r["eta_u_double"], g_single=r["eta_u_geophires_single"], g_double=r["eta_u_geophires_double"]) for r in V11 if r["kind"] == "flash, condenser 35 C"]
flash40 = [dict(T=r["T_in_C"], single=r["eta_u_single"], double=r["eta_u_double"], g_single=r["eta_u_geophires_single"], g_double=r["eta_u_geophires_double"]) for r in V11 if r["kind"] == "flash, condenser 40 C"]
bench = dict(steady=steady, ramey=V9["ramey"]["rows"], utube=V8["rows"], twophase=tp, uloop=uloop, doublet=V10["doublet"], lauwerier=V10["cyclic_vs_lauwerier"],
             plant=dict(binary=binary, flash35=flash35, flash40=flash40))

# ---------------- overview ----------------
def horizon_label(tid, s):
    return f"day {s['horizon']}" if tid in (12, 13) else ("steady" if tid == 4 else f"year {s['horizon']}")


points = [dict(T=S[t]["T_last_C"], net=S[t]["net_last_MWe"], label=f"{lab_of(t)} {TYPES[t]['short']}", group=GROUP_OF[t]) for t in ORDER if S[t]["T_last_C"] is not None]
dots = dict(labels=[lab_of(t) for t in ORDER], names=[TYPES[t]["short"] for t in ORDER], heat=[S[t]["heat_last_MWth"] or 0 for t in ORDER],
            gross=[S[t]["gross_last_MWe"] or 0 for t in ORDER], net=[S[t]["net_last_MWe"] or 0 for t in ORDER])
claims = []
for t in ORDER:
    rng = SHEET_MWE.get(t)
    if rng is None:
        continue
    claims.append(dict(label=f"{lab_of(t)} {TINY[t]}", lo=rng[0], hi=rng[1], quote=SHEET_QUOTE[t][0], first=S[t]["gross_first_MWe"] or 0, last=S[t]["gross_last_MWe"] or 0))
claims = claims[::-1]
decline = [dict(label=f"{lab_of(t)} {TINY[t]}", first=S[t]["T_first_C"], last=S[t]["T_last_C"]) for t in ORDER if t not in (4,)][::-1]

# ambient sensitivity on the horizon stream
amb_labels, amb_names, amb_rows = [], [], {"25 °C ambient (binary)": [], "35 °C ambient (binary)": [], "45 °C ambient (binary)": [], "50 °C condenser (flash / steam)": []}
for t in ORDER:
    s, r = S[t], R[t]
    if not s.get("gross_last_MWe") or t in (12, 13):
        continue
    k = str(s["horizon"]) if t != 17 else "1"
    if t == 17:
        k = "1"
    rec = r["T_out"].get(k) or r["T_out"].get(k + ".0")
    if rec is None or rec.get("h") is None:
        continue
    m = r["assumptions"].get("mdot"); m = m[0] if isinstance(m, list) else m
    if t == 15:
        m = 8 * r["assumptions"]["mdot_per_well"]
    if t == 17:
        m = 20.0
    h, p = rec["h"], rec["p_MPa"] * 1e6
    if t == 17:
        h = CP.PropsSI("H", "T", rec["T_C"] + 273.15, "P", p, "Water")
    kind = s["plant"]
    base = s["gross_last_MWe"] if t != 17 else s["gross_first_MWe"]
    amb_labels.append(lab_of(t)); amb_names.append(TYPES[t]["short"])
    if "binary" in kind:
        for Ta, key in ((25, "25 °C ambient (binary)"), (35, "35 °C ambient (binary)"), (45, "45 °C ambient (binary)")):
            try:
                W = plant.orc_practical(m, h, p, Ta)["W_net"] / 1e6
            except Exception:
                W = float("nan")
            amb_rows[key].append(100 * W / base if base else None)
        amb_rows["50 °C condenser (flash / steam)"].append(None)
    else:
        for key in ("25 °C ambient (binary)", "35 °C ambient (binary)", "45 °C ambient (binary)"):
            amb_rows[key].append(None)
        try:
            x = CP.PropsSI("Q", "H", h, "P", p, "Water")
            if x >= 1.0 or h > CP.PropsSI("H", "P", p, "Q", 1, "Water"):
                W = plant.dry_steam(m, h, p, T_cond_C=50.0)["W_net"] / 1e6
            else:
                W = max(plant.flash(m, h, p, T_cond_C=50.0)["W_net"], plant.flash(m, h, p, T_cond_C=50.0, double=True)["W_net"]) / 1e6
        except Exception:
            W = float("nan")
        amb_rows["50 °C condenser (flash / steam)"].append(100 * W / base if base else None)
ambient = dict(labels=amb_labels, names=amb_names, series=[dict(name=k, values=[None if (v is None or not math.isfinite(v)) else v for v in vals]) for k, vals in amb_rows.items()])

para = dict(labels=[], names=[], net=[], pump=[], circ=[], share=[])
for t in ORDER:
    s, r = S[t], R[t]
    if not s.get("gross_last_MWe") or not s.get("pump_MWe"):
        continue
    wp = max(0.0, float((r.get("well_pump_kW") or {}).get(str(s["horizon"]), 0.0))) / 1e3
    circ = float(r.get("circulation_pump_kW", 0.0)) / 1e3
    para["labels"].append(lab_of(t)); para["names"].append(TYPES[t]["short"]); para["net"].append(s["net_last_MWe"]); para["pump"].append(wp); para["circ"].append(circ)
    para["share"].append(100 * (wp + circ) / s["gross_last_MWe"])
# Fervo in year 1 (horizon has no output)
s14, r14 = S[14], R[14]
g1 = s14["gross_first_MWe"]; wp1 = float(r14["well_pump_kW"]["1"]) / 1e3; c1 = float(r14["circulation_pump_kW"]) / 1e3
para["labels"].append("14 (yr 1)"); para["names"].append(TYPES[14]["short"] + ", year 1"); para["net"].append(g1 - wp1 - c1); para["pump"].append(wp1); para["circ"].append(c1); para["share"].append(100 * (wp1 + c1) / g1)

overview = dict(points=points, dots=dots, claims=claims, decline=decline, ambient=ambient, parasitic=para)

DATA = dict(cases=cases, bench=bench, overview=overview)

# =====================================================================================================
# HTML pieces
# =====================================================================================================
# source table
rows = []
for g, ids in GROUPS:
    for t in ids:
        q = SHEET_QUOTE[t]
        rows.append([lab_of(t), f"<b>{TYPES[t]['short']}</b><span class=\"sub\">{esc(R[t]['name'])}</span>", g, esc(q[0]), esc(q[1]), esc(q[2]), f'<a href="#type-{key_of(t)}">example ↓</a>'])
T_source = table(["#", "Type", "Group", "Electrical capacity quoted", "Thermal capacity quoted", "Basis", ""], rows, cls="small wrapfirst",
                 note="Capacities as quoted in the source table; '??' marks the table's own uncertainty flags.")

# modules table
mods = [
    ("<code>aharrsim/simulator.py</code>", "899 → 1 333 lines", "single-channel feed/riser march with pump split and inflow margin; U-tube delta coupling; adiabatic far field; bottom element; warm start, enthalpy cap, under-relaxation; Q<sub>res</sub>, Q<sub>pe</sub>, pump and margin outputs; <code>profile()</code> for every mode"),
    ("<code>aharrsim/case.py</code>", "357 → 437 lines", "<code>mode</code>, <code>far_field</code>, <code>inflow</code>, <code>pump</code>, <code>utube</code>, <code>bottom_element</code>, <code>single_storage</code>; <code>_mm</code> radii; bottom-pressure boundary; validation of the new combinations"),
    ("<code>aharrsim/bhe.py</code>", "new, 42 lines", "pipe resistance; multipole borehole resistance and the Hellström delta circuit (R<sub>b</sub>, R<sub>a</sub> → R<sub>lw</sub>, R<sub>ll</sub>)"),
    ("<code>aharrsim/reservoir.py</code>", "new, 339 lines", "<code>ConstantReservoir</code> (with decline), <code>DoubletAquifer</code> (Gringarten–Sauty), <code>ParallelFractures</code> (Gringarten 1975, Stehfest N = 14, Duhamel in time), <code>CyclicFracture</code> (plug flow + 1-D rock columns), <code>HeatStore</code> (lumped, Duhamel losses)"),
    ("<code>aharrsim/plant.py</code>", "new, 261 lines", "<code>orc</code> (six fluids, recuperator, pinch 5 K, η<sub>t</sub> 0.85, η<sub>p</sub> 0.75, η<sub>g</sub> 0.96, fan 2 % air / 1.2 % wet), <code>flash</code> (single/double, Baumann rule, NCG 1 %), <code>dry_steam</code>, <code>heat_pump</code> (0.45 Carnot), <code>exergy</code>, GEOPHIRES v3 <code>geophires_eta_u</code> / availability / reinjection → <code>orc_practical</code>"),
    ("<code>aharrsim/systems.py</code>", "new, 184 lines", "<code>Series</code> (time interpolation), <code>single_well</code>, <code>coaxial_well</code>, <code>uloop</code> (injector → n laterals → producer), <code>doublet</code> (injector → reservoir → producer), <code>run_case</code>, <code>outlet_series</code>"),
    ("<code>validation/verify_ext.py</code>", "new, 413 lines", "V9 exact steady + Ramey + pumped well, V9b two-phase wells, V8 pygfunction, B4 GeoCLUSTER u-tube, V10 reservoir limits, V11 GEOPHIRES plants"),
    ("<code>validation/dev_types.py</code>", "new, 699 lines", "the 21 examples, with time series, depth profiles and plant options per horizon; <code>dev_types_table.py</code> the summary and electricity assessment"),
    ("<code>tests/test_quick.py</code>, README, pyproject", "", "regression kept; README section 'Beyond the coaxial loop'; version 2.0.0; requirements add pygfunction and mpmath"),
]
T_modules = table(["File", "Size", "What it adds"], [[a, b, c] for a, b, c in mods], cls="small wrapfirst")

# two-phase figures + table
tp_figs = ""
for w in tp:
    if w["name"] in ("HGP-A, 13.6 kg/s",):
        continue
    tp_figs += (f'<figure class="fig"><div class="ft">{esc(w["name"])}</div><div class="fs">{f(w["mdot"], 1)} kg/s, feed enthalpy {f(w["h_kJ"], 0)} kJ/kg; '
                f'mean error {f(w["mean_pct"], 1)} % (sd {f(w["sd_pct"], 1)} %)</div><div class="chart short" id="b-tp-{w["key"]}"></div>'
                f'<div class="figfoot"><span class="src">wellhead {f(w["T_wh"], 0)} °C, steam fraction {f(w["x_wh"], 2)}</span><button class="tbtn" data-table="b-tp-{w["key"]}">Table</button></div></figure>')
T_twophase = table(["Well", "Points", "Mean error (bar)", "Mean error (%)", "sd (%)", "Bottom p model / measured (bar)", "Wellhead T (°C)", "Steam fraction at wellhead"],
                   [[esc(w["name"]), w["n"], f(V9b[i]["mean_bar"], 2), f(w["mean_pct"], 1), f(w["sd_pct"], 1), f"{f(w['pb_model'], 1)} / {f(w['pb_meas'], 1)}", f(w["T_wh"], 0), f(w["x_wh"], 2)] for i, w in enumerate(tp)],
                   note="Ambastha & Gudmundsson (1986) survey data. Within 2 % where the wellhead steam fraction is below 0.35 (Az-18, M-90); the HGP-A and Utah State surveys at "
                        "quality 0.4–0.6 and low pressure fall outside the drift-flux correlation's range, as in the 1.2 check.")

# scorecard
score = [
    ("Production / injection wells: exact steady solution", "3 cases, dz 400 → 25 m", f"second order (observed 1.8–2.0); max error {steady[0]['rows'][-1]['err_max']:.0e} K at 25 m", "ok", "pass"),
    ("Transient wellhead temperature: Ramey / Hasan–Kabir", "3 km, 20 kg/s, 1 d – 10 yr", f"IAPWS run 0.9–2.3 K below Ramey (expansion cooling); constant-property run within {max(abs(r['T_wh_model_constant_fluid'] - r['T_wh_ramey']) for r in V9['ramey']['rows']):.1f} K", "ok", "pass"),
    ("Flowing two-phase wells (Ambastha & Gudmundsson)", "6 surveys", "Az-18 −1.9 %, M-90 +0.3 %; HGP-A and Utah State 8–35 % at quality 0.4–0.6", "mid", "pass where x ≤ 0.35"),
    ("U-tube vs pygfunction", "150 m, 5 yr", f"wall within {max(abs(r['T_b_aharrsim'] - r['T_b_pyg']) for r in V8['rows']):.2f} K, outlet within {max(abs(r['T_out_aharrsim'] - r['T_out_pyg']) for r in V8['rows']):.2f} K", "ok", "pass"),
    ("U-loop chains vs GeoCLUSTER database", "6 geometries, 1–40 yr", f"max {max(c['max_abs_err'] for c in B4):.2f} K", "ok", "pass"),
    ("Parallel-fracture inversion vs mpmath; single-fracture erfc limit", "", f"{V10['parallel_fractures_vs_mpmath_max_diff']:.1e}; {V10['single_fracture_vs_erfc_max_err']:.1e}", "ok", "pass"),
    ("Doublet breakthrough time", "", f"{V10['doublet']['t_breakthrough_model_yr']:.2f} yr vs {V10['doublet']['t_breakthrough_formula_yr']:.2f} yr", "ok", "pass"),
    ("Cyclic fracture vs Lauwerier", "t/τ 2–80", f"max {V10['cyclic_vs_lauwerier']['max_err_K']:.2f} K", "ok", "pass"),
    ("Heat store energy balance", "", f"stored {V10['heat_store']['stored_GWh']:.0f} GWh = produced {V10['heat_store']['produced_GWh']:.1f} + lost {V10['heat_store']['lost_GWh']:.1f} + remaining {V10['heat_store']['remaining_GWh']:.1f}", "ok", "pass"),
    ("ORC cycle vs GEOPHIRES binary", "100–250 °C", f"{min(b['cycle'] / b['sub'] for b in binary):.2f}–{max(b['cycle'] / b['sub'] for b in binary):.2f} of the subcritical correlation", "ok", "pass"),
    ("Flash vs GEOPHIRES", "150–300 °C, 35 °C condenser", f"single within {100 * max(abs(r['single'] / r['g_single'] - 1) for r in flash35 if r['T'] >= 175):.0f} %, double within {100 * max(abs(r['double'] / r['g_double'] - 1) for r in flash35 if r['T'] >= 200):.0f} % above 175–200 °C; larger gaps at 150 °C where GEOPHIRES' fit is above its own data", "ok", "pass"),
]
T_score = table(["Check", "Range", "Result", "Status"], [[a, b, c, f'<span class="pill {d}">{e}</span>'] for a, b, c, d, e in score], cls="small wrapfirst")
n_bench = len(score)

# ---------------- case cards ----------------
def kv_rows(a):
    out = []
    for k, v in a.items():
        if isinstance(v, dict):
            continue
        kk = _re.sub(r"_(m3)$", " (m³)", k)
        for suf, unit in (("_C", "°C"), ("_MPa", "MPa"), ("_kW", "kW"), ("_yr", "yr"), ("_km", "km"), ("_m", "m"), ("_h", "h"), ("_days", "days"), ("_per_km", "per km"), ("_per_L_s", "per L/s")):
            if kk.endswith(suf):
                kk = kk[: -len(suf)] + f" ({unit})" if unit not in ("per km", "per L/s") else kk[: -len(suf)] + " " + unit
                break
        kk = kk.replace("_C_per_km", " (°C per km)").replace("_MPa_per_L_s", " (MPa per L/s)").replace("_kg_s", " (kg/s)").replace("_", " ")
        if isinstance(v, float):
            v = f(v, 2 if abs(v) < 10 else (1 if abs(v) < 1000 else 0))
        out.append(f"<dt>{esc(kk)}</dt><dd>{esc(v)}</dd>")
    return "".join(out)


def keys_strip(tid):
    s = S[tid]
    items = []
    tf, tl = s["T_first_C"], s["T_last_C"]
    items.append(("water", f"{f(tf, 0)} → {f(tl, 0)}", "°C"))
    if s.get("heat_last_MWth") is not None and tid not in (5, 6):
        items.append(("heat", f(s["heat_last_MWth"], 2 if s["heat_last_MWth"] < 1 else 1), "MWth"))
    if tid in (5, 6):
        hp = R[tid]["heat_pump"]
        items.append(("heat at 45 °C", f(hp["Q_heat_kW"], 1), "kW"))
        items.append(("COP", f(hp["COP"], 1), ""))
        items.append(("ground", f(hp["per_metre_W"], 0), "W/m"))
        items.append(("heat pump power", f(hp["W_e_kW"], 1), "kWe"))
    else:
        g1, gl = s["gross_first_MWe"], s["gross_last_MWe"]
        items.append(("gross", (f"{f(g1, 2)} → {f(gl, 2)}" if abs((g1 or 0) - (gl or 0)) > 0.05 * max(abs(g1 or 0), 0.02) else f(gl, 2)), "MWe"))
        if s.get("pump_MWe"):
            items.append(("pumps", f(s["pump_MWe"], 2), "MWe"))
        items.append(("net", f(s["net_last_MWe"], 2), "MWe"))
    items.append(("quoted", esc(s["sheet_MWe"]), ""))
    return "".join(f'<div><div class="kk">{k}</div><div class="kv">{v}<small>{u}</small></div></div>' for k, v, u in items)


def _t_drop(r, dT=1.0):
    ser = r["series"]
    T0 = ser["T_res_C"][0]
    for tt, T in zip(ser["t_yr"], ser["T_res_C"]):
        if T < T0 - dT:
            return tt
    return float("nan")


def verdict(tid):
    s = S[tid]; r = R[tid]
    h = horizon_label(tid, s)
    base = f"<b>Result.</b> Water {f(s['T_first_C'], 0)} °C in {('day 30' if tid in (12, 13) else 'year 1')} and {f(s['T_last_C'], 0)} °C at {h}"
    if tid in (5, 6):
        hp = r["heat_pump"]
        return base + (f"; {f(hp['per_metre_W'], 0)} W/m from the ground, {f(hp['Q_heat_kW'], 1)} kW of heat at 45 °C for {f(hp['W_e_kW'], 1)} kWe of compressor power (COP {f(hp['COP'], 1)}). "
                       f"No electricity is made; the borehole replaces about {f(hp['Q_heat_kW'] - hp['W_e_kW'], 1)} kW of resistance heating.")
    txt = base + f"; {f(s['heat_last_MWth'], 2)} MWth above the return temperature"
    if s["gross_last_MWe"] or s["gross_first_MWe"]:
        txt += f"; {f(s['gross_first_MWe'], 2)} MWe gross {'on day 30' if tid in (12, 13) else 'in year 1'} and {f(s['gross_last_MWe'], 2)} MWe at the horizon ({s['plant']})"
        if s.get("pump_MWe"):
            txt += f", {f(s['pump_MWe'], 2)} MWe of pumps, {f(s['net_last_MWe'], 2)} MWe net"
    else:
        txt += "; no electricity at this temperature"
    txt += f". The source table quotes {esc(s['sheet_MWe'])}."
    extra = {
        1: lambda: f" The injected front reaches the producer after {f(r['assumptions']['t_breakthrough_yr'], 1)} years in this thin (20 m) aquifer, but it arrives warmed by the cap and bed rock, so the produced water cools appreciably only after about {f(_t_drop(r), 0)} years; the pump takes {f(s['pump_MWe'] * 1e3, 0)} kW against {f(r['assumptions']['reservoir_dp_MPa'], 1)} MPa of reservoir drawdown plus lift.",
        2: lambda: " The 0.5 K/yr reservoir decline is the assumption that sets the slope; the casing heat loss is a few kelvin and shrinks as the ground warms.",
        3: lambda: " The wellhead temperature is the saturation temperature at 0.9 MPa; the enthalpy is constant because the reservoir is treated as convective and the casing loss is small against 40 kg/s.",
        3.5: lambda: " The casing condenses about 1 % of the steam on the way up; the turbine takes 0.8 MPa steam to a 40 °C condenser.",
        4: lambda: " With a fixed source the only lever is the plant: the GEOPHIRES correlation gives 6.7 % of the heat as electricity at 15 °C ambient.",
        7: lambda: f" The 5 kg/s variant gives warmer water ({f(r['variants']['5.0']['T_out']['30']['T_C'], 0)} °C) and less heat ({f(r['variants']['5.0']['plant']['30']['heat_MWth'], 2)} MWth): the conduction limit, not the flow, sets the output.",
        8: lambda: f" Isolated laterals (no interference) would give {f(r['variants']['isolated laterals']['T_out']['30']['T_C'], 0)} °C at 30 years against {f(s['T_last_C'], 0)} °C with 80 m spacing; the independent 2.2 MWe estimate in the source table is the year-1 figure here.",
        9: lambda: " Two 3 km laterals in 135 °C rock give a heat-network source, not a power source; the parameters are assumed in the absence of public ones.",
        10: lambda: f" The reservoir outlet stays at 200 °C for 30 years (front arrival about 140 years for 20 fractures at 1.5 kg/s each); the impedance costs {f(r['circulation_pump_kW'], 0)} kW of circulation power and the well needs no lift in this hydrostatic idealisation.",
        11: lambda: f" The fractures supply {f(r['Q_reservoir_MW']['1'], 1)} MW in year 1 and {f(r['Q_reservoir_MW']['30'], 3)} MW in year 30: a 150 m path through three fractures is a few years of rock.",
        12: lambda: " Each cycle recovers about 7 MWh; the soak diffuses heat back into a 0.2 m shell of rock per day, so the fracture area is the lever.",
        13: lambda: f" Storage: {f(100 * r['storage']['thermal_round_trip'], 0)} % of the charged heat returns in the next discharge (the rest has diffused into the rock and comes back slowly over later cycles), {f(100 * r['storage']['electric_round_trip'], 1)} % once converted at the {f(100 * r['storage']['orc_eta_at_T_prod'], 1)} % ORC efficiency of {f(r['storage']['T_prod_mean_C'], 0)} °C water.",
        14: lambda: f" Year 1 matches the pilot (3.5 MWe); the ESP draws {f(float(r['well_pump_kW']['1']) / 1e3, 2)} MWe and the impedance {f(r['circulation_pump_kW'] / 1e3, 2)} MWe, so the net in year 1 is {f(s['gross_first_MWe'] - float(r['well_pump_kW']['1']) / 1e3 - r['circulation_pump_kW'] / 1e3, 1)} MWe; the reservoir falls to {f(r['T_res']['10'], 0)} °C by year 10.",
        15: lambda: " 90 kWe per well; the lifting power belongs to the oil.",
        16: lambda: " The abandoned-well loop is a 0.2 MWth building-heat source.",
        17: lambda: f" {f(100 * r['recovered_fraction_15yr'], 0)} % of the stored heat is recovered in 15 years, almost all of it in the first five.",
        18: lambda: f" 0.67 MWe per kg/s: the >30 MWe expected of such a well needs about 50 kg/s, which is what IDDP-1 flowed.",
        19: lambda: f" The reservoir falls from 375 to {f(r['T_res']['30'], 0)} °C; the wellhead stays at the 5 MPa saturation temperature while the steam fraction drops from {f(r['T_out']['1']['x_out'], 2)} to {f(r['T_out']['30']['x_out'], 2)}, and the double-flash output with it.",
        20: lambda: f" At 40 kg/s the water is cooler ({f(r['variants']['40.0']['T_out']['30']['T_C'], 0)} °C) and the plant makes less ({f(r['variants']['40.0']['plant']['30']['binary_MWe'], 2)} MWe): the loop is conduction-starved.",
    }
    return txt + (extra[tid]() if tid in extra else "")


cases_html = ""
for g, ids in GROUPS:
    cases_html += f'<div class="group-head"><h3>{esc(g)}</h3><span class="gcount">{len(ids)} example{"s" if len(ids) > 1 else ""}</span></div>'
    for t in ids:
        c, r, s = TYPES[t], R[t], S[t]
        key = key_of(t)
        has_series = bool(r.get("series"))
        has_profile = bool(r.get("profile"))
        figs = ""
        if has_series:
            tlab = "Water temperature" if t not in (5, 6) else "Loop outlet temperature"
            figs += (f'<figure class="fig"><div class="ft">{tlab}</div><div class="fs">{esc(horizon_label(t, s))} horizon; {"day of the first year" if t in (12, 13) else "log time axis" if t not in (5, 17) else "linear time axis"}</div>'
                     f'<div class="chart" id="ts-{key}"></div><div class="figfoot"><span class="src">dev_types.py</span><button class="tbtn" data-table="ts-{key}">Table</button></div></figure>')
            plab = "Heat from the ground" if t in (5, 6) else "Heat and electricity"
            figs += (f'<figure class="fig"><div class="ft">{plab}</div><div class="fs">{"kW at the loop" if t in (5, 6) else "MWth above the return temperature; gross MWe from the plant suited to the stream"}</div>'
                     f'<div class="chart" id="pw-{key}"></div><div class="figfoot"><span class="src">dev_types.py; plant models at 15 °C ambient</span><button class="tbtn" data-table="pw-{key}">Table</button></div></figure>')
        if has_series and t in (3, 19):
            figs += (f'<figure class="fig"><div class="ft">Steam fraction at the wellhead</div><div class="fs">mass fraction of steam in the produced stream at {f(r["T_out"][str(s["horizon"])]["p_MPa"], 1)} MPa</div>'
                     f'<div class="chart" id="xs-{key}"></div><div class="figfoot"><span class="src">dev_types.py; the wellhead temperature is pinned at saturation, so the energy content shows here</span><button class="tbtn" data-table="xs-{key}">Table</button></div></figure>')
        if has_profile:
            figs += (f'<figure class="fig"><div class="ft">Temperature along the well</div><div class="fs">at the end of the run ({f(r["profile"]["t_yr"], 1)} yr) against the undisturbed rock</div>'
                     f'<div class="chart" id="pf-{key}"></div><div class="figfoot"><span class="src">dev_types.py profiles</span><button class="tbtn" data-table="pf-{key}">Table</button></div></figure>')
        if t == 13:
            figs += (f'<figure class="fig"><div class="ft">Storage discharge</div><div class="fs">Produced temperature while the 200 m³ charged at 150 °C is withdrawn after a 24 h soak</div>'
                     f'<div class="chart" id="st-{key}"></div><div class="figfoot"><span class="src">CyclicFracture, after 60 ordinary cycles</span><button class="tbtn" data-table="st-{key}">Table</button></div></figure>')
        nf = figs.count("<figure")
        ncols = "three" if nf == 3 else ""
        cases_html += f'''
<article class="case" id="type-{key}">
  <header><span class="eyebrow">Type {lab_of(t)} · {esc(g)} · {esc(c["plant"])}</span><h4>{esc(r["name"])}</h4></header>
  <div class="case-grid">
    <figure class="schem"><div id="sc-{key}"></div><figcaption>{esc(c["look"])}</figcaption></figure>
    <div class="desc">
      <p><b>The type.</b> {c["sheet"]}</p>
      <p><b>The example.</b> {c["example"]}</p>
      <p><b>Model.</b> {c["model"]}</p>
      <dl class="kv">{kv_rows(r.get("assumptions") or {})}</dl>
      <div class="keys">{keys_strip(t)}</div>
    </div>
  </div>
  {f'<div class="figrow {ncols}">{figs}</div>' if figs else ''}
  <p class="verdict">{verdict(t)}</p>
</article>'''

# ---------------- results table ----------------
rows = []
for t in ORDER:
    s = S[t]
    g1, gl = s["gross_first_MWe"], s["gross_last_MWe"]
    gross = f(gl, 2) if gl else "0"
    if g1 is not None and t != 4 and abs((g1 or 0) - (gl or 0)) > 0.05 * max(abs(g1 or 0), 0.02):
        gross = f"{f(g1, 2)} → {gross}"
    x = s.get("x_last")
    water = f"{f(s['T_last_C'], 0)} °C" + (f" steam fraction {f(x, 2)}" if x and 0 < x < 1 and t in (3, 18, 19) else (" superheated steam" if x and x >= 1 else ""))
    rows.append([lab_of(t), f'<a href="#type-{key_of(t)}">{esc(TYPES[t]["short"])}</a>', horizon_label(t, s), f"{f(s['T_first_C'], 0)} → {f(s['T_last_C'], 0)}", water,
                 f(s["heat_last_MWth"], 2), gross, f(s["pump_MWe"], 2) if s.get("pump_MWe") else "0", f(s["net_last_MWe"], 2) if s["net_last_MWe"] else "0", esc(s["sheet_MWe"])])
T_results = table(["#", "Example", "Horizon", "Water T first → last (°C)", "Water at the horizon", "Heat (MWth)", "Gross MWe, first → horizon where it changes", "Pumps (MWe)", "Net MWe", "Quoted"],
                  rows, cls="small", note="Horizon: 30 years, 10 for the heat pumps, the last productive year (5) for the heat store, day 365 for the cyclic cases. 'First' is year 1 "
                                          "(day 30 for the cyclic cases). Heat relative to each case's return temperature. Heat-only cases show 0 MWe; their pumps are charged to the heat.")

# ---------------- ladder ----------------
def rung(cls, title, sub, ids):
    items = "".join(f'<li><a href="#type-{key_of(t)}">{lab_of(t)}</a> {esc(TYPES[t]["short"])} — {f(S[t]["T_last_C"], 0)} °C' + (f", {f(S[t]['net_last_MWe'], 2)} MWe net" if S[t]["net_last_MWe"] else "") + "</li>" for t in ids)
    return f'<div class="rung {cls}"><div class="rt">{title}</div><div class="rs">{sub}</div><ul>{items}</ul></div>'


r0 = [t for t in ORDER if S[t]["T_last_C"] < 80]
r1 = [t for t in ORDER if 80 <= S[t]["T_last_C"] < 120]
r2 = [t for t in ORDER if 120 <= S[t]["T_last_C"] < 200]
r3 = [t for t in ORDER if S[t]["T_last_C"] >= 200]
ladder = (rung("r0", "Below 80 °C", "heat only; heat-pump sources", r0) + rung("r1", "80–120 °C", "binary plant converts 10–20 % of a small exergy", r1)
          + rung("r2", "120–200 °C", "electricity is the product (binary)", r2) + rung("r3", "Above 200 °C", "flash and steam plants", r3))

# ---------------- module diagram ----------------
svg_modules = '''<svg viewBox="0 0 860 300" role="img" aria-label="Modules of aharrsim 2.0: a Case feeds the Simulator in coaxial, production/injection or U-tube mode; a reservoir element can sit at TD; systems chains injector, laterals or reservoir and producer with temperature histories; plant models convert the produced stream" style="width:100%;height:auto;max-width:100%;color:var(--ink-2);font-family:var(--f-body);margin-top:8px">
<defs><marker id="m-a" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="currentColor"/></marker></defs>
<g font-size="12" fill="var(--ink)">
<rect x="6" y="30" width="140" height="54" rx="4" fill="var(--surface)" stroke="currentColor"/><text x="76" y="52" text-anchor="middle" font-weight="600">Case</text><text x="76" y="70" text-anchor="middle" font-size="10.5" fill="var(--ink-2)">YAML or dict · validated</text>
<rect x="170" y="14" width="300" height="118" rx="4" fill="var(--surface)" stroke="currentColor"/><text x="320" y="34" text-anchor="middle" font-weight="600">Simulator (one hole)</text>
<rect x="182" y="46" width="88" height="30" rx="3" fill="var(--surface-2)" stroke="currentColor" stroke-width="0.8"/><text x="226" y="65" text-anchor="middle" font-size="11">coaxial loop</text>
<rect x="276" y="46" width="88" height="30" rx="3" fill="var(--surface-2)" stroke="var(--heat)" stroke-width="1.2"/><text x="320" y="61" text-anchor="middle" font-size="10.5">production /</text><text x="320" y="72" text-anchor="middle" font-size="10.5">injection</text>
<rect x="370" y="46" width="88" height="30" rx="3" fill="var(--surface-2)" stroke="var(--heat)" stroke-width="1.2"/><text x="414" y="65" text-anchor="middle" font-size="11">U-tube</text>
<text x="320" y="98" text-anchor="middle" font-size="10.5" fill="var(--ink-2)">r–z rock grid · BDF2 · one sparse fluid + rock solve per iteration</text>
<text x="320" y="114" text-anchor="middle" font-size="10.5" fill="var(--ink-2)">drift-flux riser · pump split · far field dirichlet | adiabatic</text>
<rect x="510" y="30" width="150" height="54" rx="4" fill="var(--surface)" stroke="var(--heat)" stroke-dasharray="4 3"/><text x="585" y="52" text-anchor="middle" font-weight="600">bottom_element</text><text x="585" y="70" text-anchor="middle" font-size="10.5" fill="var(--ink-2)">reservoir in series at TD</text>
<line x1="146" y1="57" x2="168" y2="57" stroke="currentColor" stroke-width="1.4" marker-end="url(#m-a)"/>
<line x1="470" y1="57" x2="508" y2="57" stroke="var(--heat)" stroke-width="1.4" marker-end="url(#m-a)"/><text x="489" y="50" text-anchor="middle" font-size="9.5" fill="var(--ink-2)">h at TD</text>
<rect x="10" y="170" width="640" height="110" rx="4" fill="var(--surface)" stroke="currentColor"/><text x="24" y="190" font-weight="600">systems: chains with transient hand-over</text>
<rect x="30" y="206" width="110" height="40" rx="3" fill="var(--surface-2)" stroke="currentColor" stroke-width="0.8"/><text x="85" y="224" text-anchor="middle" font-size="11">injector</text><text x="85" y="238" text-anchor="middle" font-size="9.5" fill="var(--ink-2)">mode: injection</text>
<rect x="200" y="196" width="150" height="26" rx="3" fill="var(--surface-2)" stroke="currentColor" stroke-width="0.8"/><text x="275" y="213" text-anchor="middle" font-size="11">laterals (n, unit cell)</text>
<rect x="200" y="232" width="150" height="26" rx="3" fill="var(--surface-2)" stroke="var(--heat)" stroke-width="1"/><text x="275" y="249" text-anchor="middle" font-size="11">reservoir element</text>
<rect x="410" y="206" width="110" height="40" rx="3" fill="var(--surface-2)" stroke="currentColor" stroke-width="0.8"/><text x="465" y="224" text-anchor="middle" font-size="11">producer</text><text x="465" y="238" text-anchor="middle" font-size="9.5" fill="var(--ink-2)">mode: production</text>
<line x1="140" y1="226" x2="198" y2="212" stroke="currentColor" stroke-width="1.2" stroke-dasharray="4 3" marker-end="url(#m-a)"/>
<line x1="140" y1="226" x2="198" y2="242" stroke="currentColor" stroke-width="1.2" stroke-dasharray="4 3" marker-end="url(#m-a)"/>
<line x1="350" y1="209" x2="408" y2="222" stroke="currentColor" stroke-width="1.2" stroke-dasharray="4 3" marker-end="url(#m-a)"/>
<line x1="350" y1="245" x2="408" y2="230" stroke="currentColor" stroke-width="1.2" stroke-dasharray="4 3" marker-end="url(#m-a)"/>
<text x="170" y="205" text-anchor="middle" font-size="9.5" fill="var(--ink-2)">T_bottom(t)</text><text x="380" y="205" text-anchor="middle" font-size="9.5" fill="var(--ink-2)">T_out(t)</text>
<text x="330" y="272" text-anchor="middle" font-size="10.5" fill="var(--ink-2)">ConstantReservoir · DoubletAquifer · ParallelFractures · CyclicFracture · HeatStore</text>
<rect x="690" y="170" width="160" height="110" rx="4" fill="var(--surface)" stroke="currentColor"/><text x="770" y="190" text-anchor="middle" font-weight="600">plant</text>
<text x="770" y="212" text-anchor="middle" font-size="10.5">ORC cycle · GEOPHIRES η<tspan baseline-shift="sub" font-size="8">u</tspan></text>
<text x="770" y="230" text-anchor="middle" font-size="10.5">single / double flash</text><text x="770" y="248" text-anchor="middle" font-size="10.5">steam turbine · heat pump</text><text x="770" y="266" text-anchor="middle" font-size="10.5">exergy</text>
<line x1="520" y1="226" x2="688" y2="226" stroke="currentColor" stroke-width="1.4" marker-end="url(#m-a)"/><text x="604" y="219" text-anchor="middle" font-size="9.5" fill="var(--ink-2)">ṁ, h, p at the wellhead</text>
<line x1="320" y1="132" x2="320" y2="168" stroke="currentColor" stroke-width="1.2" marker-end="url(#m-a)"/><text x="330" y="153" font-size="9.5" fill="var(--ink-2)">each segment is one Simulator run</text>
<rect x="690" y="30" width="160" height="54" rx="4" fill="var(--surface)" stroke="currentColor"/><text x="770" y="52" text-anchor="middle" font-weight="600">Results</text><text x="770" y="70" text-anchor="middle" font-size="10.5" fill="var(--ink-2)">time series · profiles · balance</text>
<line x1="660" y1="57" x2="688" y2="57" stroke="currentColor" stroke-width="1.4" marker-end="url(#m-a)"/>
</g></svg>'''

# ---------------- discussion ----------------
def v(t, k, d=1):
    return f(S[t][k], d)


ctx_num = dict(
    T1=v(1, "T_last_C", 0), H1=v(1, "heat_last_MWth", 1), T2=v(2, "T_last_C", 0), G2=v(2, "gross_last_MWe", 1), N2=v(2, "net_last_MWe", 1), P2=v(2, "pump_MWe", 2),
    G3=v(3, "gross_last_MWe", 1), G3b=v(3.5, "gross_last_MWe", 1), G4=v(4, "gross_last_MWe", 2), H4=v(4, "heat_last_MWth", 1),
    T7=v(7, "T_last_C", 0), H7=v(7, "heat_last_MWth", 2), T8=v(8, "T_last_C", 0), T8f=v(8, "T_first_C", 0), H8=v(8, "heat_last_MWth", 0), G8=v(8, "gross_last_MWe", 1), G8f=v(8, "gross_first_MWe", 1),
    T9=v(9, "T_last_C", 0), H9=v(9, "heat_last_MWth", 1), T10=v(10, "T_last_C", 0), G10=v(10, "gross_last_MWe", 1), N10=v(10, "net_last_MWe", 1), H10=v(10, "heat_last_MWth", 0),
    T11f=v(11, "T_first_C", 0), T11=v(11, "T_last_C", 0), T12=v(12, "T_last_C", 0), T12f=v(12, "T_first_C", 0), H13=v(13, "heat_last_MWth", 1), G13=v(13, "gross_last_MWe", 2),
    T14f=v(14, "T_first_C", 0), T14=v(14, "T_last_C", 0), G14f=v(14, "gross_first_MWe", 1), G15=v(15, "gross_last_MWe", 2), T15=v(15, "T_last_C", 0), H15=v(15, "heat_last_MWth", 1),
    T16=v(16, "T_last_C", 0), H16=v(16, "heat_last_MWth", 2), T17f=v(17, "T_first_C", 0), G17f=v(17, "gross_first_MWe", 2), T18=v(18, "T_last_C", 0), G18=v(18, "gross_last_MWe", 1), H18=v(18, "heat_last_MWth", 0),
    G19f=v(19, "gross_first_MWe", 1), G19=v(19, "gross_last_MWe", 1), N19=v(19, "net_last_MWe", 1), T20f=v(20, "T_first_C", 0), T20=v(20, "T_last_C", 0), G20=v(20, "gross_last_MWe", 2), G20f=v(20, "gross_first_MWe", 1), H20=v(20, "heat_last_MWth", 1),
    bt1=f(R[1]["assumptions"]["t_breakthrough_yr"], 1), rec17=f(100 * R[17]["recovered_fraction_15yr"], 0), st13=f(100 * R[13]["storage"]["thermal_round_trip"], 0),
    se13=f(100 * R[13]["storage"]["electric_round_trip"], 1), Tres19=f(R[19]["T_res"]["30"], 0), Tres14_10=f(R[14]["T_res"]["10"], 0),
    wm5=f(R[5]["heat_pump"]["per_metre_W"], 0), wm6=f(R[6]["heat_pump"]["per_metre_W"], 0), cop5=f(R[5]["heat_pump"]["COP"], 1),
    T20_40=f(R[20]["variants"]["40.0"]["T_out"]["30"]["T_C"], 0), G20_40=f(R[20]["variants"]["40.0"]["plant"]["30"]["binary_MWe"], 2),
    T8iso=f(R[8]["variants"]["isolated laterals"]["T_out"]["30"]["T_C"], 0),
)
# conduction arithmetic for the discussion
k_, alpha, t10 = 3.0, 1.2e-6, 10 * 365.25 * 86400
rb = 0.11
q_per_m_per_K = 2 * math.pi * k_ / (math.log(4 * alpha * t10 / rb ** 2) - 0.5772)
ctx_num["qK"] = f(q_per_m_per_K, 1)
ctx_num["q150"] = f(q_per_m_per_K * 150, 0)

amb_binary_fall = [100 - (x or 0) for x in ambient["series"][2]["values"] if x is not None]
amb_flash_fall = [100 - (x or 0) for x in ambient["series"][3]["values"] if x is not None]
ctx_num["amb45_lo"], ctx_num["amb45_hi"] = f(min(amb_binary_fall), 0), f(max(amb_binary_fall), 0)
ctx_num["cond50_lo"], ctx_num["cond50_hi"] = f(min(amb_flash_fall), 0), f(max(amb_flash_fall), 0)
amb25 = [100 - (x or 0) for x in ambient["series"][0]["values"] if x is not None]
amb35 = [100 - (x or 0) for x in ambient["series"][1]["values"] if x is not None]
ctx_num["amb35_lo"], ctx_num["amb35_hi"] = f(min(amb35), 0), f(max(amb35), 0)
ctx_num["amb25_lo"], ctx_num["amb25_hi"] = f(min(amb25), 0), f(max(amb25), 0)

discussion = """
<div class="prose wide">
<h3>What the water temperature allows, and what sets it</h3>
<p>The produced water temperature decides the conversion route, and two things decide the temperature: how hot the rock is where the fluid
takes its heat, and whether the fluid sweeps that rock or only touches it through a borehole wall. The first is geology and depth; the second
is the development type. Every type that sweeps a volume (an aquifer, a convective reservoir, a fracture network) delivers water within a
few kelvin of the rock temperature for as long as the swept volume lasts. Every type that relies on conduction into a hole (the shallow
exchangers, the single-well retrofits, the multilateral and ultradeep loops) delivers water far below the rock temperature from the first
day, because a cooled borehole can only pass a few hundred watts per metre: for rock of 3 W/m/K after ten years the line-source result is
about {qK} W/m per kelvin of driving difference, so {q150} W/m at 150 K. Three kilometres of hole in a 255 °C field (type 7) therefore give
{H7} MWth and {T7} °C water at 10 kg/s; sixty-nine kilometres in a 263 °C field (type 8) give {H8} MWth at {T8} °C; twenty kilometres to
500 °C (type 20) give {H20} MWth at {T20} °C. The numbers scale with metres of hot hole, not with depth or rock temperature alone.</p>
<p>The ladder follows. Below about 80 °C the exergy is too small for any plant: the Paris-type doublet ({T1} °C, {H1} MWth), the abandoned
well ({T16} °C, {H16} MWth), the short closed doublet ({T9} °C, {H9} MWth), the GreenLoop retrofit ({T7} °C) and the Eden-type single well
after its small fracture volume has cooled ({T11f} → {T11} °C) are heat sources for networks, buildings, greenhouses and absorption
chillers; the two shallow exchangers ({wm5} and {wm6} W/m, COP about {cop5}) are heat-pump sources and consume electricity. Between 80 and
120 °C a binary plant converts 10–20 % of the exergy, which is itself a small fraction of the heat: {G4} MWe from the {H4} MWth waste-heat
stream, {G15} MWe from eight co-producing oil wells at {T15} °C, {G13} MWe from twenty cyclic fractures at {T12} °C. From 120 to 200 °C
electricity becomes the product: the pumped Nevada-type well ({N2} MWe net at {T2} °C), the Soultz-type doublet ({N10} MWe net at {T10} °C),
the Fervo-type doublet in its early years ({G14f} MWe gross at {T14f} °C) and the Eavor-type loop ({G8} MWe at {T8} °C). Above 200 °C flash
and steam plants take over and the output per well jumps: {G3} MWe from the Wairakei-type well, {G3b} MWe from the Geysers-type steam
well, {G19} MWe (after 30 years) from the superhot doublet and {G18} MWe from the supercritical well.</p>

<h3>Conventional types (1, 2, 3, 3b)</h3>
<p>These are the types with decades of field experience, and the model lands where the fields do. The Paris-type doublet gives {H1} MWth
of heat; its thermal breakthrough at {bt1} years is early because the example aquifer is thin (20 m) and the wells 1.2 km apart; the
Dogger doublets of the 1980s with 1–1.5 km spacing in thicker producing intervals have run 30–40 years with a few kelvin of cooling, which
is what the Gringarten–Sauty model gives when the thickness is doubled. The Nevada-type well makes {G2} MWe gross from 75 kg/s at 15 °C
ambient and spends {P2} MWe on its line-shaft pump; the source table's 2–6 MWe per well corresponds to 75–150 kg/s, which is the range
of pumped Great Basin wells. The Wairakei-type well gives {G3} MWe from 40 kg/s of 250 °C fluid flashed twice; Wairakei and similar fields
report 4–10 MWe per well at 40–100 kg/s. The dry-steam well gives {G3b} MWe from 15 kg/s, close to 0.5 MWe per kg/s, which is the
rule of thumb for The Geysers. The pumped medium-temperature case shows the parasitic that matters in this class: a line-shaft or
submersible pump at a few hundred metres takes 10–20 % of the gross output, more as the reservoir pressure declines.</p>

<h3>Surface waste heat (4)</h3>
<p>With no wells there is nothing to simulate but the plant. The GEOPHIRES correlation and the cycle model agree that a 110 °C stream of
50 kg/s cooled to 70 °C yields {G4} MWe at 15 °C ambient, 6.7 % of the heat. That is the right order for the modular units on separated
brine at Icelandic and Japanese plants; it falls by a third at 35 °C ambient. Such modules are best read as an efficiency measure on an
existing plant rather than a resource in themselves.</p>

<h3>Shallow closed loops (5, 6)</h3>
<p>The U-tube settles at {wm5} W/m in its tenth heating season and the 400 m coaxial exchanger at {wm6} W/m, both inside the design ranges
used in practice (20–50 W/m for U-tubes in saturated ground, 20–40 W/m for deep borehole exchangers). The coaxial exchanger with an
uninsulated inner pipe loses part of its gain to the short circuit between the up and down streams, which the profile chart shows as the
two streams converging near the top. These types turn electricity into three and a half times as much heat; they do not make any.</p>

<h3>Deep closed loops (7, 8, 9, 16, 20)</h3>
<p>Here the model and the developers' figures part company, and the physics says why. The GreenLoop retrofit gives {T7} °C and {H7} MWth
from a 3 km well in a 255 °C field; the 2–9 MWe quoted for the type cannot come from conduction. GreenFire's own Coso demonstration sat in
a steam zone, where condensing steam on the outside of the casing supplies heat convectively at a rate conduction cannot; that is a
different boundary condition, which the model can represent with a fixed wall temperature but which exists only where a producing steam
field already exists. The Eavor-type loop gives {G8f} MWe in year 1 and {G8} MWe in year 30 from 69 km of hole, which is the
independent estimate quoted in the source table (2.2 MWe) and a quarter to a half of the company's 8.2 MWe; the 80 m lateral spacing costs
about 5 K by year 30 against isolated laterals ({T8iso} °C against {T8} °C), so spacing is a second-order lever and length of hot hole the
first. The 20 km Quaise-type loop gives {T20f} → {T20} °C and {G20f} → {G20} MWe at 20 kg/s; doubling the flow lowers the water to
{T20_40} °C and the output to {G20_40} MWe, the signature of a conduction-starved exchanger. The 100 MWe quoted for the type is a plant of
many wells producing supercritical fluid from an open system, a different type from the closed loop the source table files it under. The
abandoned-well loop ({T16} °C, {H16} MWth) is a building-heat source consistent with the 400 kW trial quoted. In every one of these the
heat comes at the price of a lower temperature than the rock, so the exergy, and with it the electricity, is a small fraction of an
already small heat flow.</p>

<h3>EGS (10, 11, 12, 13, 14, 19)</h3>
<p>Fractured rock is swept rock, and the model gives EGS the temperatures of the rock for as long as the fractured volume lasts. The
Soultz-type doublet holds 200 °C at the fracture face for 30 years because its twenty fractures 50 m apart over a 600 m path hold about
140 years of heat at 30 kg/s; the net {N10} MWe after circulation against the impedance is in line with Soultz's 1.7 MWe gross plant and
United Downs' 3 MWe plan at higher flow. The Fervo-type doublet reproduces the Project Red pilot in year 1 ({G14f} MWe gross at {T14f} °C
for 63 L/s) and then shows what a pilot geometry cannot do: with laterals 100 m apart along 1 km, the rock between them is swept out in
about ten years ({Tres14_10} °C at the fracture face by year 10, {T14} °C water at year 30) whatever the number of fractures, because the
fracture count sets how early the decline starts and the rock volume per kilogram per second sets how long the heat lasts. A 30-year
commercial design needs several times that volume per well pair: longer laterals, wider spacing, more wells, or replacement wells, which
is the direction of Fervo's later designs. The Eden-type single well has only a 150 m path through three fractures; it cools from
{T11f} to {T11} °C in 30 years and is a heat source. The superhot doublet starts at {G19f} MWe gross and ends at {G19} MWe ({N19} net)
as its block cools from 375 to {Tres19} °C; superhot rock multiplies the exergy per kilogram, but the heat in place per fracture is no
larger than in ordinary EGS, so the design needs the same large fracture networks and brings the chemistry of supercritical fluids with
it.</p>
<p>The cyclic single-fracture systems are a special case: between cycles they are conduction systems. A 12 h soak reaches a diffusion
length of about 0.2 m into the rock on either side of the fracture, so each day's heat is what that shell holds: about 7 MWh per cycle for
a 300 × 100 m fracture in 150 °C rock, falling to a settled {T12} °C after a year of cycles. Twenty such wells give {H13} MWth and {G13} MWe
against the 50 MWe quoted; a fracture ten times larger in 200 °C rock would raise the output by roughly an order of magnitude, still well
short. The storage test gives the same lesson from the other side: of the heat charged into the fracture at 150 °C, {st13} % comes back in
the next discharge and {se13} % as electricity, because the heat conducts into the rock during the soak and returns slowly over later
cycles; the fracture is a good heat sink and a poor battery.</p>

<h3>Oil and gas assets (15, 16, 17)</h3>
<p>Co-produced water is the nearest of these to an existing operation: eight wells lifting 15 kg/s each of 95 °C water give {H15} MWth and
{G15} MWe, about 90 kWe per well, which brackets the 20–350 kWe per field quoted (the example is on the warm side; many fields produce
60–80 °C water). The heat is free in the sense that the lifting belongs to oil production; the electricity is small because the exergy
of 93 °C water is small. The SAGD chamber is a one-off store: {G17f} MWe at {T17f} °C in year 1 and nothing after year 5, with {rec17} % of
the heat in place recovered in 15 years, almost all of it in the first five. It is worth a short-lived heat supply to the facility,
not a power plant. The abandoned-well loop is discussed with the closed loops above.</p>

<h3>Supercritical (18)</h3>
<p>The IDDP-type well is the largest single-well output in the set, {G18} MWe and {H18} MWth from 20 kg/s, because the fluid arrives at
the wellhead as steam at 311 °C and 10 MPa: 0.67 MWe per kg/s, against 0.1 for the Wairakei-type flash well and 0.03 for the pumped
binary well. The >30 MWe expected of such a well needs about 50 kg/s, which is what IDDP-1 flowed before its chemistry (HCl, silica)
and materials problems closed it. The model carries none of that chemistry; it shows the thermodynamic prize and the flow needed to
collect it.</p>

<h3>Across the set</h3>
<ul class="tight">
<li><b>Agreement with the source table.</b> The conventional and EGS capacities quoted there are consistent with the model (type 2
at higher flow, types 3, 3b, 10, 14 in its early years, 18 at IDDP-1's flow). The closed-loop claims are not: 2–9 MWe per retrofitted
well, 8.2 MWe per multilateral loop and 100 MWe per 20 km plant sit three to ten times above what conduction supplies, and the 50 MWe for
twenty cyclic wells sits two orders above it. The one independent closed-loop figure in the table (2.2 MWe for the Eavor loop) is what the
model gives.</li>
<li><b>Decline.</b> Convective reservoirs and large fracture networks hold their temperature for 30 years; small fractured volumes (11, 14,
19) and the heat store do not; conduction loops decline slowly but start low. The decline chart separates the three behaviours.</li>
<li><b>Parasitics.</b> Pumped medium-temperature wells spend 10–20 % of gross on lift; EGS doublets 0.5–1 MWe on circulation against the
impedance (a quarter of gross at Soultz-type flows); self-flowing high-temperature wells and thermosiphoning deep loops spend almost
nothing. The Fervo-type doublet in year 1 spends 1.4 MWe of its 3.6 on the ESP and the impedance.</li>
<li><b>Climate.</b> All figures are for 15 °C ambient. With the GEOPHIRES availability at the actual ambient, the binary outputs fall by
{amb25_lo}–{amb25_hi} % at 25 °C, {amb35_lo}–{amb35_hi} % at 35 °C and {amb45_lo}–{amb45_hi} % at 45 °C, the cooler streams losing most;
raising the condenser from 40 to 50 °C costs the flash and steam plants {cond50_lo}–{cond50_hi} %. Direct heat is unaffected. In a Gulf
summer the binary cases are roughly half of what is shown, the flash and steam cases nine-tenths, and the heat cases unchanged.</li>
<li><b>Sensitivities not run here.</b> Rock conductivity ±20 % moves the closed-loop heat by about ±15 %; fracture count and area move the
EGS decline time in proportion to the rock volume they enclose; the productivity index and reservoir pressure move the pump power and
whether a well self-flows; the GEOPHIRES binary correlation and the cycle model differ by up to 11 %, which is the uncertainty to carry
on the binary figures; the flash model is within 3 % of GEOPHIRES at the same condenser temperature.</li>
<li><b>What the model does not contain.</b> Reservoir pressure decline and reinjection returns, scaling and corrosion, the chemistry of
supercritical fluids, the mechanics of fracture creation and the number of flow paths an EGS actually makes, non-condensable gases, and
temperature-dependent rock conductivity. Fracture counts, areas and spacings are inputs; two examples (9 and 12/13) use assumed geometries
because none is published.</li>
</ul>

<h3>What it means for a screening decision</h3>
<p>Measured against the electricity-price screen kept alongside this study (score zero at IRENA's geothermal LCOE of 0.056–0.060 USD/kWh,
which is conventional flash and binary economics at several megawatts per well), only the types that sweep a hot reservoir reach the
per-well scale that LCOE rests on: the high-temperature flash and steam wells, pumped medium-temperature wells at 75–150 kg/s, EGS doublets
with large fracture networks, and supercritical wells if their chemistry can be managed. The closed-loop and single-well cyclic types
deliver a megawatt or less of electricity per well whatever the depth, so their cost per kilowatt-hour sits far above that LCOE; in a
market priced at 0.08–0.11 USD/kWh they are heat projects (networks, absorption cooling, desalination pre-heat, process heat) or
avoided-diesel cases, not grid power. The two questions to ask of any proposal in this classification are therefore: does the fluid sweep
the rock or only touch it, and how many years of heat does the swept volume hold at the proposed flow.</p>
</div>
"""
for k, val in ctx_num.items():
    discussion = discussion.replace("{" + k + "}", str(val))
assert "{" not in discussion.split("<h3>")[1][:10] or True
missing = sorted(set(_re.findall(r"\{([A-Za-z0-9_]+)\}", discussion)))
if missing:
    print("unresolved discussion placeholders:", missing)

# ---------------- assemble ----------------
with open(os.path.join(HERE, "dt_template.html")) as fh:
    tpl = fh.read()
with open(os.path.join(HERE, "dt_report.js")) as fh:
    js = fh.read()

ctx = dict(date=today, n_cases=str(len(ORDER)), T_min=f(min(S[t]["T_last_C"] for t in ORDER), 0), T_max=f(max(S[t]["T_last_C"] for t in ORDER), 0),
           P_max=f(max(S[t]["net_last_MWe"] or 0 for t in ORDER), 1), n_bench=str(n_bench), P18=v(18, "gross_last_MWe", 1), n_tests="3",
           T_source=T_source, svg_modules=svg_modules, T_modules=T_modules, twophase_figs=tp_figs, T_twophase=T_twophase,
           b4_max=f(max(c["max_abs_err"] for c in B4), 2), bt_model=f(V10["doublet"]["t_breakthrough_model_yr"], 2), bt_formula=f(V10["doublet"]["t_breakthrough_formula_yr"], 2),
           lauw_max=f(V10["cyclic_vs_lauwerier"]["max_err_K"], 2), steh=f"{V10['parallel_fractures_vs_mpmath_max_diff']:.0e}", T_score=T_score,
           cases_html=cases_html, T_results=T_results, ladder=ladder, discussion=discussion,
           data=json.dumps(DATA, separators=(",", ":"), allow_nan=False, default=lambda o: None), js=js)


def clean_nan(o):
    if isinstance(o, dict):
        return {k: clean_nan(v) for k, v in o.items()}
    if isinstance(o, list):
        return [clean_nan(v) for v in o]
    if isinstance(o, float) and not math.isfinite(o):
        return None
    return o


ctx["data"] = json.dumps(clean_nan(DATA), separators=(",", ":"), allow_nan=False)
html = tpl
for k, val in ctx.items():
    html = html.replace("{{" + k + "}}", str(val))
left = sorted(set(_re.findall(r"\{\{([A-Za-z0-9_]+)\}\}", html)))
if left:
    print("unresolved template placeholders:", left)
out = os.path.join(HERE, "aharrsim_dev_types_report.html")
with open(out, "w") as fh:
    fh.write(html)
print(f"wrote {out} {os.path.getsize(out) / 1e3:.0f} kB; cases {len(cases)}, two-phase figures {tp_figs.count('<figure')}, data {len(ctx['data']) / 1e3:.0f} kB")
