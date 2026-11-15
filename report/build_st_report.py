"""Build the uncertainty / open-system study report (aharrsim_uncertainty_report.html) from validation/results/ST*.json."""
from __future__ import annotations

import datetime as _dt
import json
import math
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RES = os.path.join(ROOT, "validation", "results")
sys.path.insert(0, os.path.join(ROOT, "validation")); sys.path.insert(0, ROOT)
import study_cl as CL  # noqa: E402
import study_open as OP  # noqa: E402


def L(n, default=None):
    p = os.path.join(RES, n + ".json")
    return json.load(open(p)) if os.path.exists(p) else default


def f(v, d=1):
    if v is None or (isinstance(v, float) and not math.isfinite(v)):
        return "–"
    return f"{v:,.{d}f}".replace("-", "−")


def table(head, rows, cls="num", note=None):
    h = f'<div class="tbl"><table class="{cls}"><thead><tr>' + "".join(f"<th>{c}</th>" for c in head) + "</tr></thead><tbody>"
    h += "".join("<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>" for r in rows) + "</tbody></table></div>"
    return h + (f'<p class="tnote">{note}</p>' if note else "")


def fig(cid, title, sub, src="", cls=""):
    return (f'<figure class="fig"><div class="ft">{title}</div><div class="fs">{sub}</div><div class="chart {cls}" id="{cid}"></div>'
            f'<div class="figfoot"><span class="src">{src}</span><button class="tbtn" data-table="{cid}">Table</button></div></figure>')


charts = {}
C = {}      # numbers for the prose
today = _dt.date.today().strftime("%-d %B %Y")

# =============================================================================================== part 1
OAT = L("ST1_oat", []); MC = L("ST1_mc", []); OPT = L("ST1_opt", [])
by = {}
for r in OAT:
    if "failed" not in r:
        by[(r["case"], r["tag"])] = r
SHORT = {7: "3 km coaxial retrofit (type 7)", 8: "Multilateral loop, 69 km of hole (type 8)", 20: "20 km coaxial loop (type 20)"}
p1_html = ""
for c in (7, 8, 20):
    base = by.get((c, "base"))
    if not base:
        continue
    b30 = base["30"]
    C[f"b{c}_T"], C[f"b{c}_H"], C[f"b{c}_E"] = b30["T"], b30["MWth"], b30["MWe"]
    rows_h, rows_T, rows_e = [], [], []
    for var, lo, hi, label, kind in CL.OAT[c]:
        a, b = by.get((c, f"oat|{var}|low")), by.get((c, f"oat|{var}|high"))
        if not a or not b:
            continue
        tag = " (design)" if kind == "des" else ""
        rows_h.append(dict(label=label + tag, low=a["30"]["MWth"], high=b["30"]["MWth"], swing=abs(b["30"]["MWth"] - a["30"]["MWth"]), kind=kind, var=var))
        rows_T.append(dict(label=label + tag, low=a["30"]["T"], high=b["30"]["T"], swing=abs(b["30"]["T"] - a["30"]["T"])))
        rows_e.append(dict(label=label + tag, low=a["30"]["MWe"], high=b["30"]["MWe"], swing=abs(b["30"]["MWe"] - a["30"]["MWe"]), kind=kind, var=var))
    for rr in (rows_h, rows_T, rows_e):
        rr.sort(key=lambda r: -r["swing"])
    C[f"rank{c}"] = [r["var"] for r in rows_h if r["kind"] == "unc"]
    C[f"oat{c}"] = {r["var"]: r for r in rows_h}
    C[f"oate{c}"] = {r["var"]: r for r in rows_e}
    charts[f"p1-th-{c}"] = dict(kind="tornado", rows=rows_h, base=b30["MWth"], x="Heat at year 30 (MWth)", unit="MWth")
    charts[f"p1-tt-{c}"] = dict(kind="tornado", rows=rows_T, base=b30["T"], x="Outlet temperature at year 30 (°C)", unit="°C")
    charts[f"p1-te-{c}"] = dict(kind="tornado", rows=rows_e, base=b30["MWe"], x="Gross electricity at year 30 (MWe)", unit="MWe")
    # physics cases
    cats, hv, ev, prow = ["Base (constant k, conduction)"], [b30["MWth"]], [b30["MWe"]], []
    for var, v, label in CL.PHYS[c]:
        r = by.get((c, f"phys|{var}|{v}"))
        if r:
            cats.append(label); hv.append(r["30"]["MWth"]); ev.append(r["30"]["MWe"])
            C[f"phys{c}_{var}_{v}"] = r["30"]
    charts[f"p1-ph-{c}"] = dict(kind="bar", horizontal=True, cats=cats, y="MWth | MWe at year 30", series=[dict(name="Heat (MWth)", data=hv, c="heat"), dict(name="Gross electricity (MWe)", data=ev, c=0)])
    # flow sweep
    var, vals = CL.SWEEP[c]
    sw = [by.get((c, f"sweep|{var}|{v}")) for v in vals]
    sw = [(v, r) for v, r in zip(vals, sw) if r]
    charts[f"p1-sw-{c}"] = dict(kind="xy", x=dict(name="Flow (kg/s)"), y=dict(name="°C  |  MWth ×10  |  MWe ×100", min=0), xl="kg/s", legcols=1,
                                series=[dict(name="Outlet temperature, year 30 (°C)", data=[[v, r["30"]["T"]] for v, r in sw], c=0, sym=True),
                                        dict(name="Heat (MWth ×10)", data=[[v, 10 * r["30"]["MWth"]] for v, r in sw], c="heat", sym=True),
                                        dict(name="Gross electricity (MWe ×100)", data=[[v, 100 * r["30"]["MWe"]] for v, r in sw], c=2, sym=True)])
    C[f"sw{c}"] = [(v, r["30"]["T"], r["30"]["MWth"], r["30"]["MWe"], r["30"]["pump_kW"]) for v, r in sw]
    # Monte Carlo
    mc = [r for r in MC if r["case"] == c and "failed" not in r]
    mc_tbl = ""
    if mc:
        A = {k: np.array([r["30"][k] for r in mc]) for k in ("T", "MWth", "MWe", "MWe35", "MWe45")}
        A1 = {k: np.array([r["1"][k] for r in mc]) for k in ("T", "MWth", "MWe")}
        P = lambda a: np.percentile(a, [10, 50, 90])
        C[f"mc{c}"] = {k: P(v).tolist() for k, v in A.items()}; C[f"mc{c}_n"] = len(mc)
        xs = np.sort(A["MWth"]); ys = (np.arange(len(xs)) + 0.5) / len(xs) * 100
        charts[f"p1-mc-{c}"] = dict(kind="xy", x=dict(name="Heat at year 30 (MWth)"), y=dict(name="Cumulative probability (%)", min=0, max=100), xl="MWth",
                                    series=[dict(name=f"Monte Carlo, {len(mc)} runs", data=[[float(a), float(b)] for a, b in zip(xs, ys)], c=0, step=True)],
                                    marks=[dict(x=b30["MWth"], name="base")])
        rowsmc = [["Outlet temperature, year 1 (°C)"] + [f(v, 0) for v in P(A1["T"])], ["Outlet temperature, year 30 (°C)"] + [f(v, 0) for v in P(A["T"])],
                  ["Heat, year 1 (MWth)"] + [f(v, 2) for v in P(A1["MWth"])], ["Heat, year 30 (MWth)"] + [f(v, 2) for v in P(A["MWth"])],
                  ["Gross electricity, year 1, 15 °C ambient (MWe)"] + [f(v, 2) for v in P(A1["MWe"])],
                  ["Gross electricity, year 30, 15 °C ambient (MWe)"] + [f(v, 2) for v in P(A["MWe"])],
                  ["Gross electricity, year 30, 35 °C ambient (MWe)"] + [f(v, 2) for v in P(A["MWe35"])],
                  ["Gross electricity, year 30, 45 °C ambient (MWe)"] + [f(v, 2) for v in P(A["MWe45"])]]
        mc_tbl = table(["", "P10", "P50", "P90"], rowsmc, note=f"Latin hypercube, {len(mc)} full simulations; triangular distributions over the tornado ranges of gradient, conductivity, heat capacity"
                       + (" and insulation" if c != 8 else "") + "; conductivity falling with temperature in every run; design held at the base case. P10 is the low case.")
    p1_html += f'''
<h3 id="p1-{c}">{SHORT[c]}</h3>
<div class="prose wide"><p>{{{{p1_text_{c}}}}}</p></div>
<div class="figrow">
{fig(f"p1-th-{c}", "What moves the heat", "Year-30 heat at the low and high value of each variable, others at base; blue below base, red above", "study_cl.py oat")}
{fig(f"p1-te-{c}", "What moves the electricity", "Year-30 gross electricity (best plant, 15 °C ambient)", "study_cl.py oat") if b30["MWe"] > 0.01 else fig(f"p1-tt-{c}", "What moves the outlet temperature", "Year-30 outlet temperature; no case in these ranges makes electricity", "study_cl.py oat")}
</div>
<div class="figrow three">
{fig(f"p1-ph-{c}", "Physics beyond constant-k conduction", "Year-30 heat and electricity", "study_cl.py oat", "short")}
{fig(f"p1-sw-{c}", "Flow rate", "Year 30; heat and electricity scaled to share the axis", "study_cl.py oat", "short")}
{fig(f"p1-mc-{c}", "Monte Carlo", "Year-30 heat: % of runs at or below the value", "study_cl.py mc", "short") if mc else ""}
</div>
{mc_tbl}'''

opt_rows = []
for r in sorted([r for r in OPT if "failed" not in r], key=lambda r: (r["case"], r["P"]["nu"], r["P"]["mdot"])):
    y = r["30"]
    opt_rows.append([SHORT[r["case"]].split(" (")[0], f(r["P"]["nu"], 0), f(r["P"]["mdot"], 0), f(y["T"], 0), f(y["MWth"], 1), f(y["MWe"], 2)])
    C[f"opt{r['case']}_{r['P']['nu']:g}_{r['P']['mdot']:g}"] = y
T_opt = table(["Case", "Nu of the zone", "Flow (kg/s)", "Outlet T, year 30 (°C)", "Heat (MWth)", "Gross electricity (MWe)"], opt_rows,
              note="Convecting or steam zone represented as an effective rock conductivity Nu × k over the lower half of the 3 km well "
                   "or around the laterals of the multilateral loop. Nu 3–10 is free convection in a permeable formation; Nu 30–100 approaches a casing held at reservoir temperature by condensing steam. Both loops circulate by thermosiphon at these flows, so no pump power is deducted.") if opt_rows else ""

# range table
rng_rows = []
for c in (7, 8, 20):
    for var, lo, hi, label, kind in CL.OAT[c]:
        if (c, f"oat|{var}|low") not in by or (c, f"oat|{var}|high") not in by:
            continue
        rng_rows.append([SHORT[c].split(" (")[0], label, "uncertainty" if kind == "unc" else "design lever", f"{CL.BASE[c][var]:,g}"])
T_ranges = table(["Case", "Variable and its low / high value", "Kind", "Base"], rng_rows, cls="small wrapfirst")

# =============================================================================================== part 2
SC = L("ST2_screen"); CF = L("ST2_confirm", []); CY = L("ST2_cyclic", [])
p2_html = ""
if SC:
    yrs = SC["years"]
    for sid_s, o in SC["systems"].items():
        sid = int(sid_s); s = OP.SYS[sid]; b = o["base"]
        C[f"s{sid}"] = o
        q = s["quoted"]
        charts[f"p2-T-{sid}"] = dict(kind="xy", x=dict(name="Time (years)", min=0, max=30), y=dict(name="Reservoir outlet temperature (°C)"), xl="year", legcols=1,
                                     series=[dict(name=f"{d['N']:.0f} fractures, {d['S_m']:.1f} m apart" + (", tortuosity 5" if d["tau"] > 1 else ""), data=[[a, bb] for a, bb in zip(yrs, d["T_res"])],
                                                  c=i, dash=d["tau"] > 1) for i, d in enumerate(o["density_curves"])],
                                     marks=[dict(x=min(b["t_sweep_yr"], 30), name=f"piston sweep {b['t_sweep_yr']:.0f} yr")] if b["t_sweep_yr"] < 30 else None)
        charts[f"p2-N-{sid}"] = dict(kind="xy", x=dict(name="Number of fractures in the block", log=True, min=5, max=5000), y=dict(name="30-year mean net electricity (MWe)", min=0), xl="N",
                                     series=[dict(name=f"tortuosity {d['tau']:.0f}", data=[[r["N"], r["mean_net"]] for r in d["rows"]], c=i, sym=True) for i, d in enumerate(o["density"])],
                                     marks=[dict(y=q[0], name=f"quoted {q[0]:g} MWe")] if q[0] < 2.2 * max(r["mean_net"] for d in o["density"] for r in d["rows"]) else None)
        ser = [dict(name="base fractures, base volume", data=[[r["m"], r["mean_net"]] for r in o["flow"][0]["rows"]], c="ref", dash=True, sym=True)]
        for i, v in enumerate([v for v in o["volume"] if v["vf"] in (1, 2, 4, 8, 16)]):
            ser.append(dict(name=f"10 m spacing, volume ×{v['vf']} ({v['V_km3']:.2f} km³)", data=[[r["m"], r["mean_net"]] for r in v["rows"]], c=i, sym=True))
        charts[f"p2-m-{sid}"] = dict(kind="xy", x=dict(name="Flow (kg/s)", min=0, max=300), y=dict(name="30-year mean net electricity (MWe)", min=0), xl="kg/s", legcols=1, series=ser,
                                     marks=[dict(y=qq, name=f"quoted {qq:g} MWe") for qq in sorted(set(q))])
        ch_rows = [[f(r["f"] * 100, 0) + " %", f(r["phi"] * 100, 0) + " %", f(r["T30"], 0), f(r["t_sweep_yr"], 0), f(r["net1"], 2), f(r["net30"], 2), f(r["mean_net"], 2), f(r["dp_MPa"], 1)] for r in o["channel"]]
        T_ch = table(["Active fractures", "Swept volume", "Reservoir T, year 30 (°C)", "Piston life (yr)", "Net year 1 (MWe)", "Net year 30 (MWe)", "30-yr mean net (MWe)", "Reservoir Δp (MPa)"], ch_rows,
                     note=f"Dense design (10 m spacing, {int(s['W'] / 10)} fractures) at the base flow of {s['m']:.0f} kg/s and base volume.")
        cf = [r for r in CF if r.get("sid") == sid and "failed" not in r]
        order = {c["tag"]: i for i, c in enumerate(OP.CONFIRM[sid])}
        cf.sort(key=lambda r: order.get(r["tag"], 9))
        cf_rows = [[r["tag"], f(r["N"], 0), f(r["m"], 0), f(r["W"] * r["L"] * s["B"] / 1e9, 3), f"{f(r['T_wh'][0], 0)} → {f(r['T_wh'][-1], 0)}", f"{f(r['gross'][0], 2)} → {f(r['gross'][-1], 2)}",
                    f"{f(r['screen_gross'][0], 2)} → {f(r['screen_gross'][-1], 2)}", f(r["W_circ"] + r["pump_kW"][-1] / 1e3, 2), f(max(r["gross"][-1] - r["W_circ"] - max(r["pump_kW"][-1], 0) / 1e3, 0.0), 2)] for r in cf]
        C[f"cf{sid}"] = cf
        T_cf = table(["Design", "Fractures", "Flow (kg/s)", "Volume (km³)", "Wellhead T, year 1 → 30 (°C)", "Gross, full chain, year 1 → 30 (MWe)", "Gross, screening (MWe)", "Pumps (MWe)", "Net, year 30 (MWe)"], cf_rows,
                     cls="small wrapfirst", note="Full chain: injector → Gringarten fractures → producer with the wellbore solver (9⅝ in wells above 65 kg/s), 30 years. The screening model is within 3 % of it.") if cf else ""
        need_rows = []
        for nd in o["need"]:
            d = nd["design"]
            need_rows.append([f"{nd['target']:g} MWe", "not reached on the grid (circulation power caps the net)" if not d else f"{d['m']:.0f} kg/s", "" if not d else f(d["V_km3"], 2),
                              "" if not d else f"{d['W']:.0f} × {d['L']:.0f} × {d['B']:.0f} m", "" if not d else f(d["N"], 0), "" if not d else f(d["dp_MPa"], 1), "" if not d else f(d["W_circ"], 1), "" if not d else f(d["mean_net"], 1)])
        T_need = table(["Quoted capacity", "Flow", "Stimulated volume (km³)", "Block W × L × B", "Fractures at 10 m", "Reservoir Δp (MPa)", "Circulation power (MWe)", "30-yr mean net (MWe)"], need_rows, cls="small wrapfirst",
                       note="Smallest flow × volume on the grid that sustains the quoted figure for 30 years (year-30 net at least 90 % of year 1).")
        p2_html += f'''
<h3 id="p2-{sid}">{s["name"]}</h3>
<div class="prose wide"><p>{{{{p2_text_{sid}}}}}</p></div>
<div class="figrow three">
{fig(f"p2-T-{sid}", "Fracture density and thermal decline", f"Reservoir outlet temperature at the base flow ({s['m']:.0f} kg/s) and base volume ({b['V_km3']:.3f} km³)", "study_open.py screen", "short")}
{fig(f"p2-N-{sid}", "Density and tortuosity", "30-year mean net electricity at base flow and volume", "study_open.py screen", "short")}
{fig(f"p2-m-{sid}", "Flow and stimulated volume", "30-year mean net electricity; fractures 10 m apart", "study_open.py screen", "short")}
</div>
{T_need}
<h4>Channelling and partial sweep</h4>
{T_ch}
<h4>Full-chain confirmation</h4>
{T_cf}'''
if CY:
    ser = []
    for i, Tr in enumerate(sorted(set(r["T_rock"] for r in CY))):
        ser.append(dict(name=f"{Tr:.0f} °C rock", data=[[r["area_m2"] / 1e6, r["MWe_20_wells"]] for r in CY if r["T_rock"] == Tr], c=i, sym=True))
    charts["p2-cyc"] = dict(kind="xy", x=dict(name="Fracture wall area per well (km²)", log=True, min=0.05, max=1.2), y=dict(name="Gross electricity, 20 wells (MWe)", min=0), xl="km²", series=ser)
    ser = []
    for i, Tr in enumerate(sorted(set(r["T_rock"] for r in CY))):
        ser.append(dict(name=f"{Tr:.0f} °C rock", data=[[r["area_m2"] / 1e6, r["T_day180"]] for r in CY if r["T_rock"] == Tr], c=i, sym=True))
    charts["p2-cycT"] = dict(kind="xy", x=dict(name="Fracture wall area per well (km²)", log=True, min=0.05, max=1.2), y=dict(name="Produced water, day 180 (°C)"), xl="km²", series=ser)
    C["cy"] = CY

# =============================================================================================== part 3
P3 = L("ST3_open_uncertainty")
p3_html = ""
SEL = {10: 1, 14: 1, 19: 2}
if P3:
    yrs = P3["years"]
    for sid_s, o in P3["systems"].items():
        sid = int(sid_s); s = OP.SYS[sid]
        C[f"p3_{sid}"] = o
        d = o["designs"][SEL[sid]]
        rows = [dict(label=t["label"], low=t["low"]["mean_net"], high=t["high"]["mean_net"], swing=abs(t["high"]["mean_net"] - t["low"]["mean_net"]), var=t["var"]) for t in d["tornado"]]
        for r in rows:
            if r["var"] == "kT":
                r["high"] = d["base"]["mean_net"]; r["swing"] = abs(r["low"] - r["high"])
        rows.sort(key=lambda r: -r["swing"])
        C[f"p3rank{sid}"] = rows
        charts[f"p3-t-{sid}"] = dict(kind="tornado", rows=rows, base=d["base"]["mean_net"], x="30-year mean net electricity (MWe)", unit="MWe")
        charts[f"p3-h-{sid}"] = dict(kind="xy", x=dict(name="Time (years)", min=0, max=30), y=dict(name="Net electricity (MWe)", min=0), xl="year", legcols=1,
                                     series=[dict(name="favourable: gradient +20 %, k 4.0, heat capacity +15 %", data=[[a, b] for a, b in zip(yrs, d["favourable"]["net"])], c=2),
                                             dict(name="base", data=[[a, b] for a, b in zip(yrs, d["base"]["net"])], c=0),
                                             dict(name="base at 45 °C ambient", data=[[a, b] for a, b in zip(yrs, d["hot_ambient"]["net"])], c=3, dash=True),
                                             dict(name="unfavourable: gradient −20 %, k 2.0 falling with T, heat capacity −15 %", data=[[a, b] for a, b in zip(yrs, d["unfavourable"]["net"])], c="heat")],
                                     marks=[dict(y=qq, name=f"quoted {qq:g} MWe") for qq in sorted(set(s["quoted"]))])
        cats = [dd["label"] for dd in o["designs"]]
        mrows = []
        for dd in o["designs"]:
            m15, m35 = dd["mc"]["15"], dd["mc"]["35"]
            mrows.append([dd["label"], f(dd["base"]["mean_net"], 2), f(dd["unfavourable"]["mean_net"], 2), f(dd["favourable"]["mean_net"], 2), " / ".join(f(v, 2) for v in m15["mean_net"]),
                          " / ".join(f(v, 2) for v in m15["net30"]), " / ".join(f(v, 2) for v in m35["mean_net"]), f(100 * m15["p_zero30"], 0) + " %"])
        T_m = table(["Design", "Base mean net (MWe)", "Unfavourable", "Favourable", "Monte Carlo mean net P10 / P50 / P90, 15 °C", "Net at year 30 P10 / P50 / P90, 15 °C", "Mean net P10 / P50 / P90, 35 °C",
                     "Chance of no net output at year 30"], mrows, cls="small wrapfirst",
                    note="Monte Carlo: 3,000 samples of the screening model per design; triangular gradient ±20 %, conductivity 2.0–4.0 W/m/K falling with temperature, heat capacity ±15 %, "
                         "active fractures 10–100 % (mode 50 %), swept volume 40–100 % (mode 80 %), impedance ×0.5–×2. Quoted capacity: " + " – ".join(f"{v:g}" for v in sorted(set(s["quoted"]))) + " MWe.")
        p3_html += f'''
<h3 id="p3-{sid}">{s["name"]}</h3>
<div class="prose wide"><p>{{{{p3_text_{sid}}}}}</p></div>
<div class="figrow">
{fig(f"p3-t-{sid}", "What moves the net electricity", f"30-year mean net, design: {d['label']}", "study_open.py part3")}
{fig(f"p3-h-{sid}", "Favourable and unfavourable ends", f"Net electricity over 30 years, design: {d['label']}", "study_open.py part3")}
</div>
{T_m}'''

# =============================================================================================== text
import st_report_text as TXT  # noqa: E402
texts = TXT.build(C, f)

tpl = open(os.path.join(HERE, "dt_template.html")).read()
style = re.search(r"<style>.*?</style>", tpl, re.S).group(0)
head = '<meta charset="utf-8">\n' + tpl[: tpl.index("<style>")].replace("Twenty Geothermal Development Types", "Geothermal Uncertainty Study")
body = TXT.PAGE
ctx = dict(date=today, p1_html=p1_html, p2_html=p2_html, p3_html=p3_html, T_opt=T_opt, T_ranges=T_ranges,
           cyc_figs=(f'<div class="figrow">{fig("p2-cycT", "Cyclic single fracture: produced temperature", "Day 180 of daily cycles, 200 m³ per cycle; fracture area raised at constant volume (thinner, wider)", "study_open.py cyclic", "short")}'
                     f'{fig("p2-cyc", "Cyclic single fracture: electricity from 20 wells", "Gross binary output on the daily-average stream", "study_open.py cyclic", "short")}</div>') if CY else "")
ctx.update(texts)
style = style.replace("</style>", ".figrow.three .chart.short{height:400px}</style>")
html = head + style + body
for _ in range(3):
    for k, v in ctx.items():
        html = html.replace("{{" + k + "}}", str(v))


def clean(o):
    if isinstance(o, dict):
        return {k: clean(v) for k, v in o.items() if v is not None or k in ("low", "high")}
    if isinstance(o, list):
        return [clean(v) for v in o]
    if isinstance(o, (np.floating, float)):
        return None if not math.isfinite(o) else round(float(o), 5)
    if isinstance(o, np.integer):
        return int(o)
    return o


html = html.replace("{{data}}", json.dumps(clean(dict(charts=charts)), separators=(",", ":"))).replace("{{js}}", open(os.path.join(HERE, "st_report.js")).read())
left = sorted(set(re.findall(r"\{\{([A-Za-z0-9_]+)\}\}", html)))
if left:
    print("unresolved:", left)
out = os.path.join(HERE, "aharrsim_uncertainty_report.html")
open(out, "w").write(html)
print(f"wrote {out} {os.path.getsize(out) / 1e3:.0f} kB, {len(charts)} charts")
