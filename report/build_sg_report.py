"""Build the laboratory-comparison report (aharrsim_segs_report.html) from validation/results/SG_*.json.

The sandbox experiment of Ma et al. (2026) on a single-well enhanced geothermal system against aharrsim: the type 11 example of the
benchmark, the swept-rock chart, a numerical model of the sandbox (aharrsim.porous) and a dipole element (reservoir.SingleWellDipole).
Page skeleton and styles are those of the earlier reports (dt_template.html); charts by sg_report.js (ECharts).
"""
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


def L(n):
    return json.load(open(os.path.join(RES, n + ".json")))


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


lab, cal, fld, chk, rec = L("SG_lab"), L("SG_calib"), L("SG_field"), L("SG_checks"), L("SG_recommended")
REP = lab["reported"]
LAB = lab["lab"]
iso, unc = cal["isotropic"], cal["uncoupled"]
ani = cal.get("anisotropic")
pub = dict(typical=lab["typical"], flow=lab["flow"], T0=lab["T0"], T_in=lab["T_in"], spacing=lab["spacing"])
charts, C = {}, {}
today = _dt.date.today().strftime("%-d %B %Y")
K_fit, a_fit = cal["K_fit"], cal.get("aniso_fit")


def series(r, tmax=90.0, step=2):
    t, T = np.array(r["t_min"]), np.array(r["T_prod"])
    m = t <= tmax + 1e-9
    return [[float(a), float(b)] for a, b in zip(t[m][::step], T[m][::step])]


# =============================================================================================== 1. the nominal case
typ_iso, typ_pub, typ_unc = iso["typical"], pub["typical"], unc["typical"]
T60_T0series = 60.0 * (1 - REP["T0_series"]["60"]["decay60"] / 100.0)
T60_flowseries = 60.0 * (1 - 23.3 / 100.0)
T70_heat = 21.0 + REP["typical"]["heat_at_70min_kW"] * 1e3 / (4.2e6 * 0.5 / 3600.0)
charts["typ"] = dict(x=dict(name="Time since the start of injection (min)", min=0, max=90), y=dict(name="Produced water (°C)", min=25, max=65), legcols=1,
                     series=[dict(name="Model, permeability fitted", data=series(typ_iso), c=0),
                             dict(name="Model, published properties", data=series(typ_pub), c=1),
                             dict(name="Model, flow not affected by temperature", data=series(typ_unc), c=2, dash=True),
                             dict(name="Reported, nominal case", pts=True, symbol="diamond", c="ink",
                                  data=[[30.0, 60.0], [60.0, 37.0], [70.0, T70_heat]]),
                             dict(name="Reported at 60 min in two other series", pts=True, symbol="triangle", c="ink",
                                  data=[[60.0, T60_T0series], [60.0, T60_flowseries]])])
C.update(on_pub=typ_pub["onset_1K_min"], on_fit=typ_iso["onset_1K_min"], on_unc=typ_unc["onset_1K_min"], T60_pub=typ_pub["T60"], T60_fit=typ_iso["T60"], T60_unc=typ_unc["T60"],
         avg_pub=typ_pub["T_avg60"], avg_fit=typ_iso["T_avg60"], avg_unc=typ_unc["T_avg60"], x_pub=typ_pub["x_onset_cyl"], x_fit=typ_iso["x_onset_cyl"], x_unc=typ_unc["x_onset_cyl"],
         rate_fit=typ_iso["max_decline_C_per_min"], rate_pub=typ_pub["max_decline_C_per_min"], K_fit=K_fit, a_fit=a_fit, ts_cyl=typ_pub["t_s_cyl_min"],
         rhoc_b=typ_pub["rhoc_b"] / 1e6, lam=typ_pub["lam"])

# =============================================================================================== 2. plateau against flow
flows = ["0.2", "0.3", "0.4", "0.5", "0.6"]
x_rep = float(np.mean([v["x_cyl"] for v in lab["plateau_scaling"].values()]))
C["x_rep"] = x_rep
C["V_th"] = float(np.mean([v["thermal_volume_m3"] for v in lab["plateau_scaling"].values()]))
qs = np.linspace(0.2, 0.6, 21)
ts_of = lambda q: C["ts_cyl"] * 0.5 / q
charts["onset"] = dict(x=dict(name="Flow rate (m³/h)", min=0.1, max=0.7), y=dict(name="Time on the plateau (min)", min=0, max=80), legcols=1,
                       series=[dict(name="Model, permeability fitted", data=[[float(q), iso["flow"][q]["onset_1K_min"]] for q in flows], c=0, sym=True),
                               dict(name="Model, published properties", data=[[float(q), pub["flow"][q]["onset_1K_min"]] for q in flows], c=1, sym=True),
                               dict(name="Model, flow not affected by temperature", data=[[float(q), unc["flow"][q]["onset_1K_min"]] for q in flows], c=2, dash=True),
                               dict(name=f"Constant flushed share, {x_rep:.2f}", data=[[float(q), x_rep * ts_of(q)] for q in qs], c="ref", w=1.2),
                               dict(name="Reported plateau", pts=True, symbol="diamond", c="ink", data=[[float(q), m] for q, m in REP["plateau_min"].items()])])

# =============================================================================================== 3. reported against model
def dev(a, b):
    return None if a is None or b is None else a - b


rowsT = []
def row(label, reported, getter, d=1, unit=""):
    vals = [getter(v) for v in (pub, iso, ani)]
    rowsT.append([label, (f(reported, d) if isinstance(reported, (int, float)) else reported) + unit] + [f(v, d) for v in vals])
    return vals


devs, devs_a = [], []
for q, m in REP["plateau_min"].items():
    row(f"Plateau at {q} m³/h (min)", m, lambda v, q=q: v["flow"][q]["onset_1K_min"], 0)
for q in ("0.2", "0.6"):
    r_ = REP["avg_flow_series"][q]
    v = row(f"Mean production temperature, 60 min, {q} m³/h (°C)", r_["T"], lambda v, q=q: v["flow"][q]["T_avg60"])
    devs.append(v[1] - r_["T"]); devs_a.append(v[2] - r_["T"])
for T0 in ("50", "60", "70"):
    r_ = REP["T0_series"][T0]
    v = row(f"Mean production temperature, 60 min, pack at {T0} °C (°C)", r_["T_avg"], lambda v, T0=T0: v["T0"][T0]["T_avg60"], 2)
    if T0 != "60":
        devs.append(v[1] - r_["T_avg"]); devs_a.append(v[2] - r_["T_avg"])
for q in ("0.2", "0.6"):
    r_ = REP["avg_flow_series"][q]
    row(f"Mean heat flow, 60 min, {q} m³/h (kW)", r_["kW"], lambda v, q=q: v["flow"][q]["kW_avg60"])
for T0 in ("50", "60", "70"):
    r_ = REP["T0_series"][T0]
    row(f"Mean heat flow, 60 min, pack at {T0} °C (kW)", r_["kW"], lambda v, T0=T0: v["T0"][T0]["kW_avg60"], 2)
for T0 in ("50", "60", "70"):
    r_ = REP["T0_series"][T0]
    row(f"Temperature at 60 min, pack at {T0} °C (°C)", float(T0) * (1 - r_["decay60"] / 100.0), lambda v, T0=T0: v["T0"][T0]["T60"])
row("Injection 21 → 41 °C: change of the mean production temperature (K)", REP["Tin_series"]["dT_avg_21_to_41"], lambda v: v["T_in"]["41"]["T_avg60"] - v["T_in"]["21"]["T_avg60"])
row("Injection 21 → 41 °C: change of the mean heat flow (kW)", REP["Tin_series"]["d_kW_21_to_41"], lambda v: v["T_in"]["41"]["kW_avg60"] - v["T_in"]["21"]["kW_avg60"], 2)
row("Spacing 90 → 120 cm: change of the mean production temperature (K)", REP["spacing"]["dT_avg"], lambda v: v["spacing"]["1.20"]["T_avg60"] - v["spacing"]["0.90"]["T_avg60"])
row("Spacing 90 → 120 cm: change of the mean heat flow (kW)", REP["spacing"]["d_kW"], lambda v: v["spacing"]["1.20"]["kW_avg60"] - v["spacing"]["0.90"]["kW_avg60"])
row("Temperature at 60 min, 90 cm spacing (°C)", 60.0 * (1 - REP["spacing"]["decay_90"] / 100.0), lambda v: v["spacing"]["0.90"]["T60"])
row("Temperature at 60 min, 120 cm spacing (°C)", 60.0 * (1 - REP["spacing"]["decay_120"] / 100.0), lambda v: v["spacing"]["1.20"]["T60"])
row("Steepest decline, nominal case (°C/min)", REP["typical"]["decline_C_per_min"], lambda v: v["typical"]["max_decline_C_per_min"], 2)
row("Heat flow on the plateau, nominal case (kW)", REP["typical"]["heat_first_30min_kW"], lambda v: 4.2e6 * 0.5 / 3600.0 * (60.0 - 21.0) / 1e3)
row("Heat flow at 70 min, nominal case (kW)", REP["typical"]["heat_at_70min_kW"], lambda v: 4.2e6 * 0.5 / 3600.0 * (v["typical"]["T70"] - 21.0) / 1e3)
C["maxdev"] = max(abs(x) for x in devs); C["maxdev_a"] = max(abs(x) for x in devs_a)
C["devs"] = devs
T_cmp = table(["Quantity", "Reported", "Model: published properties", f"Model: permeability {K_fit * 1e10:.1f} × 10⁻¹⁰ m²", f"Model: horizontal permeability {a_fit:.1f} × vertical"], rowsT,
              cls="num wrapfirst",
              note=("Reported temperatures at 60 minutes are computed from the percentage declines of the paper, which it does not define; they are read as a share of the initial temperature in °C. "
                    "The two fitted columns were each adjusted to one number only, the 55.29 °C of the nominal case (pack at 60 °C); every other row is a prediction. "
                    "Model plateau: time until the produced water has lost 1 K."))

# decay percentages of the flow series (two possible readings)
dec = REP["decay60_percent_flow_series"]
rowsD = []
for i, q in enumerate(("0.3", "0.4", "0.5", "0.6")):
    rowsD.append([f"{q} m³/h", f"{f(dec[i], 1)} %", f"{f(dec[i + 1], 1)} %"] + [f"{f(v['flow'][q]['decay60'], 1)} %" for v in (pub, iso, ani)])
T_decay = table(["Decline after 60 min", "Reported, first four values", "Reported, last four values", "Model: published", "Model: permeability fitted", "Model: anisotropy fitted"], rowsD,
                note="The paper prints five declines (8.3, 17.3, 23.3, 29.7 and 32.2 %) for four flow rates. With the last four the 0.5 m³/h value (29.7 %) agrees with the 29.67 % of the temperature series.")

# =============================================================================================== 4. temperature field (fitted model)
g = iso["grid"]
rc, zc = np.array(g["rc"]), np.array(g["zc"])
R_cm = np.arange(-48.0, 48.1, 2.0); Z_cm = np.arange(2.0, 149.0, 2.0)       # full section through the axis, 2 cm cells


def resample(Tij):
    T = np.array(Tij)                       # (nr, nz)
    out = np.zeros((len(Z_cm), len(R_cm)))
    for jz, z in enumerate(Z_cm / 100.0):
        col = np.array([np.interp(z, zc, T[i, :]) for i in range(len(rc))])
        out[jz, :] = np.interp(np.abs(R_cm) / 100.0, rc, col)
    return np.round(out, 1).tolist()


snap = iso["typical"]["snapshots"]
zi = iso["typical"]["z_inj"]
for tmin in ("30", "60"):
    charts[f"field{tmin}"] = dict(kind="field", r=[int(v) for v in R_cm], z=[int(v) for v in Z_cm], T=resample(snap[tmin]), min=21.0, max=60.0,
                                  marks=[dict(x="40", name="vertical sensor line"), dict(y="96", name="injection"), dict(y="56", name="sensor line across the box"), dict(y="4", name="production")])

# fronts along the sensor lines
fr = iso["fronts"]; frp = pub["typical"]["fronts"]
def span(d, t):
    v = d.get(t)
    return "–" if not v else f"{f(v['top_cm'], 0)}–{f(v['bottom_cm'], 0)} cm"
first = lambda d: next((f"{t} min, depth {span(d, t)}" for t in ("10", "20", "30", "40") if d.get(t)), "–")
rowsF = [["Line across the box, 55 cm above the base: reach of a 20 K drop at 40 min", "over 80 % of the monitored distance; little change at the ends",
          f"{f(fr['L55']['40'], 0)} cm from the axis", f"{f(frp['L55']['40'], 0)} cm from the axis"],
         ["Line across the box: outward speed of the 20 K drop", "about 0.75 cm/min", f"{f(fr['L55_cm_per_min'], 2)} cm/min", f"{f(frp['L55_cm_per_min'], 2)} cm/min"],
         ["Vertical line, 10 cm from the wall: first 5 K of cooling", "none for 20 min; at 30 min, depth 50 cm", first(fr["H40"]), first(frp["H40"])],
         ["Vertical line: depth cooled by 5 K at 30 min", "about 50 cm", span(fr["H40"], "30"), span(frp["H40"], "30")],
         ["Vertical line: downward speed of the front, 30–60 min", "about 1.0 cm/min", f"{f(fr['H40_down_cm_per_min'], 2)} cm/min", f"{f(frp['H40_down_cm_per_min'], 2)} cm/min"],
         ["Vertical line: depth reached at 60 min", "about 80 cm at that speed; stable below 110 cm", f"{f(fr['H40']['60']['bottom_cm'], 0)} cm", f"{f(frp['H40']['60']['bottom_cm'], 0)} cm"],
         ["Vertical line: cooled above the injection screen at 60 min", "stable above 30 cm", f"up to {f(fr['H40']['60']['top_cm'], 0)} cm of depth", f"up to {f(frp['H40']['60']['top_cm'], 0)} cm of depth"]]
T_fronts = table(["Along the sensor lines", "Reported", "Model: permeability fitted", "Model: published properties"], rowsF, cls="small wrapfirst",
                 note="Depths are below the top of the box; the injection screen is at 50–60 cm. The model is read at whole 10-minute marks. The paper places the vertical line 10 cm from the wall "
                      "in its methods and 30 cm from the axis in its results; the first is used here, and at the second the model is further ahead. The line across the box has twelve sensors 8 cm apart, "
                      "taken here as 4 to 44 cm either side of the axis. Outward speeds are between 10 and 40 minutes.")
E_60 = REP["T0_series"]["60"]["kW"] * 3600.0 / 1e3            # MJ taken in the hour, nominal case, from the reported mean heat flow
E_place = typ_pub["rhoc_b"] * typ_pub["V_box"] * 39.0 / 1e6     # MJ in the box above the injection temperature
C.update(E_60=E_60, E_place=E_place, pct_out=100.0 * E_60 / E_place)
C.update(r_40=fr["L55"]["40"], r_speed=fr["L55_cm_per_min"], f_first=span(fr["H40"], "30"), f_down=fr["H40_down_cm_per_min"], f_60=fr["H40"]["60"]["bottom_cm"], fp_60=frp["H40"]["60"]["bottom_cm"], fp_down=frp["H40_down_cm_per_min"])

# =============================================================================================== 5. sensitivities of the nominal case
S = lab["sensitivity"]
order = ["as published", "flow not coupled to temperature", "buoyancy only", "viscosity only", "permeability 5.6e-10 m2 (Kozeny-Carman recomputed)",
         "horizontal permeability twice the vertical", "horizontal permeability four times the vertical", "1500 kg/m3 read as the density of the glass", "glass at 2500 kg/m3",
         "spacing read between the near ends of the screens", "conductivity of the pack 1.0 W/m/K", "no dispersion", "dispersivity 5 mm"]
nice = {"as published": "Published properties, cold water heavier and more viscous", "flow not coupled to temperature": "Flow not affected by temperature",
        "buoyancy only": "Weight of the cold water only", "viscosity only": "Viscosity of the cold water only",
        "permeability 5.6e-10 m2 (Kozeny-Carman recomputed)": "Permeability 5.6 × 10⁻¹⁰ m² (Kozeny–Carman for 1 mm beads)",
        "horizontal permeability twice the vertical": "Horizontal permeability twice the vertical", "horizontal permeability four times the vertical": "Horizontal permeability four times the vertical",
        "1500 kg/m3 read as the density of the glass": "1,500 kg/m³ read as the density of the glass, not of the pack", "glass at 2500 kg/m3": "Glass at 2,500 kg/m³",
        "spacing read between the near ends of the screens": "90 cm read between the near ends of the screens", "conductivity of the pack 1.0 W/m/K": "Conductivity of the pack 1.0 W/m/K",
        "no dispersion": "No thermal dispersion", "dispersivity 5 mm": "Dispersivity 5 mm"}
rowsS = [["<b>Reported</b>", "about 30", f"{f(T60_T0series, 1)} (37–46)", "55.29", f"{f(x_rep, 2)}"]]
rowsS.append([f"<b>Permeability fitted, {K_fit * 1e10:.1f} × 10⁻¹⁰ m²</b>", f(typ_iso["onset_1K_min"], 1), f(typ_iso["T60"], 1), f(typ_iso["T_avg60"], 2), f(typ_iso["x_onset_cyl"], 2)])
for k in order:
    r = S[k]
    rowsS.append([nice[k], f(r["onset_1K_min"], 1), f(r["T60"], 1), f(r["T_avg60"], 2), f(r["x_onset_cyl"], 2)])
T_sens = table(["Nominal case (0.5 m³/h, 60 °C, 21 °C, 90 cm)", "Plateau (min)", "T at 60 min (°C)", "Mean T, 60 min (°C)", "Flushed share at onset"], rowsS, cls="num wrapfirst",
               note="Each row changes one thing from the first model row (published properties). Flushed share: of the heat capacity of the bead pack between the screen centres.")
C.update(on_visc=S["viscosity only"]["onset_1K_min"], on_buoy=S["buoyancy only"]["onset_1K_min"], on_a2=S["horizontal permeability twice the vertical"]["onset_1K_min"],
         on_a4=S["horizontal permeability four times the vertical"]["onset_1K_min"], avg_a2=S["horizontal permeability twice the vertical"]["T_avg60"],
         on_glass=S["glass at 2500 kg/m3"]["onset_1K_min"], on_gap=S["spacing read between the near ends of the screens"]["onset_1K_min"])

# =============================================================================================== 6. what a longer test would show
RECi = rec["isotropic"]
names = ["nominal: 0.5 m3/h, 21 C, 90 cm", "0.5 m3/h, 21 C, 120 cm", "recommended: 0.4 m3/h, 31 C, 120 cm"]
labels = ["Nominal: 0.5 m³/h, 21 °C, 90 cm", "0.5 m³/h, 21 °C, 120 cm", "Recommended: 0.4 m³/h, 31 °C, 120 cm"]
ser = []
for i, (n, lb) in enumerate(zip(names, labels)):
    r = RECi[n]
    t, th = np.array(r["t_min"]), np.array(r["theta"])
    ser.append(dict(name=lb, data=[[float(a), float(b)] for a, b in zip(t[::4], th[::4])], c=i))
charts["rec"] = dict(x=dict(name="Time since the start of injection (min)", min=0, max=180), y=dict(name="Share of the temperature rise kept", min=0, max=1.05, hideMax=True), legcols=1,
                     series=ser, bands=[dict(x0=0, x1=60, name="first hour: where the paper compares")])
rr = RECi["recommended: 0.4 m3/h, 31 C, 120 cm"]; rn = RECi["nominal: 0.5 m3/h, 21 C, 90 cm"]; r12 = RECi["0.5 m3/h, 21 C, 120 cm"]
C.update(rec_onset=rr["onset_1K_min"], rec_T60=rr["T60"], rec_decay=rr["decay60"], rec_E=rr["energy_to_onset_MJ"], nom_E=rn["energy_to_onset_MJ"], s12_onset=r12["onset_1K_min"],
         s12_E=r12["energy_to_onset_MJ"], rec_kW=4.2e6 * 0.4 / 3600.0 * (60.0 - 31.0) / 1e3, nom_kW=4.2e6 * 0.5 / 3600.0 * (60.0 - 21.0) / 1e3,
         rec_theta60=rr["theta60"], nom_theta60=rn["theta60"])
rowsR = []
for n, lb, q, Ti in zip(names, labels, (0.5, 0.5, 0.4), (21.0, 21.0, 31.0)):
    r = RECi[n]
    rowsR.append([lb, f(4.2e6 * q / 3600.0 * (60.0 - Ti) / 1e3, 1), f(r["onset_1K_min"], 0), f(r["energy_to_onset_MJ"], 0), f(100 * r["theta60"], 0) + " %", f(r["x_onset_cyl"], 2)])
T_rec = table(["Setting", "Heat flow on the plateau (kW)", "Plateau (min)", "Heat taken on the plateau (MJ)", "Rise kept at 60 min", "Flushed share at onset"], rowsR, cls="num wrapfirst",
              note="Model with the fitted permeability, pack at 60 °C. Heat with the paper's constant of 4.2 MJ/m³/K.")

# =============================================================================================== 7. scale-up curves
cv = fld["curves"]
xb, thb = np.array(cv["box_x"]), np.array(cv["box"])
m = (xb > 0.2) & (xb < 6)
px = cv["piston_x_box"]
charts["curves"] = dict(x=dict(name="x = t Q (ρc)w / (D³ (ρc)b)", log=True, min=0.2, max=40), y=dict(name="Share of the temperature rise kept", min=0, max=1.05, hideMax=True), legcols=1,
                        series=[dict(name="Large block (dipole element)", data=[[a, b] for a, b in zip(cv["x"], cv["unbounded"])], c=0),
                                dict(name="Closed cylinder of radius 0.56 D (the sandbox)", data=[[float(a), float(b)] for a, b in zip(xb[m][::2], thb[m][::2])], c=1),
                                dict(name="Uniform sweep of that cylinder", data=[[0.2, 1.0], [px, 1.0], [px, 0.0], [40.0, 0.0]], c="ref", w=1.2)])
C.update(xbox99=cv["x_at_0.99"]["box"], xbox90=cv["x_at_0.9"]["box"], xbox50=cv["x_at_0.5"]["box"], xun99=cv["x_at_0.99"]["unbounded"], xun90=cv["x_at_0.9"]["unbounded"],
         xun50=cv["x_at_0.5"]["unbounded"], px=px)

# spacing needed
sp = fld["spacing_needed"]
rowsP = []
for mdot in (5.0, 10.0, 20.0, 50.0, 100.0):
    a = [r for r in sp if r["m"] == mdot and r["years"] == 10.0][0]; b = [r for r in sp if r["m"] == mdot and r["years"] == 30.0][0]
    k = (C["xun99"] / C["xbox99"]) ** (1 / 3)        # a sandbox-like bounded block needs a larger spacing for the same first arrival
    rowsP.append([f(mdot, 0), f(a["D_first_arrival_m"], 0), f(b["D_first_arrival_m"], 0), f(b["D_10pct_m"], 0), f(b["D_first_arrival_m"] * (0.746 / C["xbox99"]) ** (1 / 3), 0),
                  f(b["V_ideal_km3"], 3)])
T_space = table(["Flow (kg/s)", "First arrival at 10 yr: spacing (m)", "First arrival at 30 yr (m)", "10 % of the rise lost at 30 yr (m)",
                 "Block bounded like the sandbox, first arrival at 30 yr (m)", "Ideal sweep volume for 30 yr (km³)"], rowsP,
                note="Rock 2,650 kg/m³ × 850 J/kg/K, water 60 → 183 °C at 47 MPa, densely fractured (local thermal equilibrium). Spacing: between the injection and the production zone on the well. "
                     "Ideal sweep: the volume rule of the benchmark, V = ṁ c̄ τ / (0.8 (ρc)r).")
s10 = [r for r in sp if r["m"] == 10.0 and r["years"] == 30.0][0]
C.update(D30=s10["D_first_arrival_m"], D30_10=s10["D_10pct_m"], D30_box=s10["D_first_arrival_m"] * (0.746 / C["xbox99"]) ** (1 / 3), V30=s10["V_ideal_km3"])

# type 11
t11 = fld["type11"]["cases"]
keys = list(t11.keys())
k_bench, k_dip, k_dipS, k_big = keys[0], keys[1], keys[2], keys[3]
t_ref = t11[k_dip]["t_first_yr"] / 0.746012
tt = np.array(t11[k_dip]["series"]["t_yr"])
tt_f = np.linspace(0.25, 30.0, 120)
box_T = 60.0 + np.interp(tt_f / t_ref, xb, thb) * (183.0 - 60.0)
from aharrsim.reservoir import SingleWellDipole, ParallelFractures, cbar_w  # noqa: E402
cw = cbar_w(60.0, 183.0, 47.0e6); Qv = 10.0 / 950.0; rhoc_r = 2650.0 * 850.0
els = {k_bench: ParallelFractures(Q=Qv, N=3, L_flow=150.0, B=150.0, S=30.0, k_r=3.2, rhoc_r=rhoc_r, T_rock_C=183.0, T_inj_C=60.0, rhoc_w=cw * 950.0),
       k_dip: SingleWellDipole(Q=Qv, D=150.0, rhoc_b=rhoc_r, T_rock_C=183.0, T_inj_C=60.0, rhoc_w=cw * 950.0),
       k_big: SingleWellDipole(Q=Qv, D=fld["type11"]["D_holds_m"], rhoc_b=rhoc_r, T_rock_C=183.0, T_inj_C=60.0, rhoc_w=cw * 950.0)}
YR = 365.25 * 86400.0
curve = lambda el: [[float(t), float(el.T_out(60.0 + 273.15, t * YR) - 273.15)] for t in tt_f]
charts["t11"] = dict(x=dict(name="Years of circulation", min=0, max=30), y=dict(name="Water leaving the rock (°C)", min=40, max=200, interval=40), legcols=1,
                     series=[dict(name="Benchmark: three fractures, 0.002 km³", data=curve(els[k_bench]), c=0),
                             dict(name="Sandbox scaled up: 150 m, closed block", data=[[float(a), float(b)] for a, b in zip(tt_f, box_T)], c=1),
                             dict(name="Zones 150 m apart, large block", data=curve(els[k_dip]), c=2),
                             dict(name=f"Zones {fld['type11']['D_holds_m']:.0f} m apart, large block", data=curve(els[k_big]), c=3)])
rowsE = []
for k, lb in ((k_bench, "Benchmark: three fractures, 150 m path, 0.002 km³"), (k_dip, "Dipole, zones 150 m apart, densely fractured"),
              (k_dipS, "Dipole, 150 m, flowing fractures 30 m apart"), (k_big, f"Dipole, zones {fld['type11']['D_holds_m']:.0f} m apart")):
    r = t11[k]
    rowsE.append([lb, f(r["t_first_yr"], 1)] + [f(r["T_out"][y], 0) for y in ("1", "5", "10", "30")] + [f(r["heat_MWth"]["1"], 2) + " → " + f(r["heat_MWth"]["30"], 2),
                                                                                                         f(r["binary_MWe"]["1"], 2) + " → " + f(r["binary_MWe"]["30"], 2)])
T_t11 = table(["Rock at the toe of the type 11 well", "Front crossing or first arrival (yr)", "Wellhead, yr 1 (°C)", "yr 5", "yr 10", "yr 30", "Heat, yr 1 → 30 (MWth)", "Binary, yr 1 → 30 (MWe)"],
              rowsE, cls="num wrapfirst", note="4.8 km coaxial well, 183 °C at the toe, 10 kg/s, 40 °C return: the well of the benchmark's type 11, full wellbore solver, only the rock element changed.")
C.update(b_T1=t11[k_bench]["T_out"]["1"], b_T30=t11[k_bench]["T_out"]["30"], d_T1=t11[k_dip]["T_out"]["1"], d_T10=t11[k_dip]["T_out"]["10"], d_T30=t11[k_dip]["T_out"]["30"],
         d_first=t11[k_dip]["t_first_yr"], g_T30=t11[k_big]["T_out"]["30"], D_holds=fld["type11"]["D_holds_m"], t_ref=t_ref, box_bt=C["xbox99"] * t_ref, box_half=C["xbox50"] * t_ref,
         d_MWe1=t11[k_dip]["binary_MWe"]["1"], d_MWe30=t11[k_dip]["binary_MWe"]["30"], d_H1=t11[k_dip]["heat_MWth"]["1"], d_H30=t11[k_dip]["heat_MWth"]["30"],
         b_first=t11[k_bench]["t_first_yr"])

# =============================================================================================== 8. injection above or below
sb = lab["sensitivity"]["injection through the lower screen"]
ib = iso["inject_bottom"]
charts["updown_lab"] = dict(x=dict(name="Time since the start of injection (min)", min=0, max=75), y=dict(name="Produced water (°C)", min=25, max=65), legcols=1,
                            series=[dict(name="Injection above, fitted permeability", data=series(typ_iso, 75), c=0),
                                    dict(name="Injection below, fitted permeability", data=series(ib, 75), c=0, dash=True),
                                    dict(name="Injection above, published permeability", data=series(typ_pub, 75), c=1),
                                    dict(name="Injection below, published permeability", data=series(sb, 75), c=1, dash=True)])
C.update(ib_on=ib["onset_1K_min"], ib_avg=ib["T_avg60"], sb_on=sb["onset_1K_min"], sb_avg=sb["T_avg60"])
fb = fld["field_buoyancy"]
def fser(k):
    s_ = fb[k]["series"]
    return [[float(a), float(b)] for a, b in zip(s_["t_yr"], s_["theta"])]
charts["updown_field"] = dict(x=dict(name="Years of circulation", min=0, max=30), y=dict(name="Share of the temperature rise kept", min=0, max=1.05, hideMax=True), legcols=1,
                              series=[dict(name="Flow not affected by temperature", data=fser("no coupling, inject above"), c="ref"),
                                      dict(name="10⁻¹⁴ m², injection above", data=fser("K 1e-14, inject above"), c=0),
                                      dict(name="10⁻¹⁴ m², injection below", data=fser("K 1e-14, inject below"), c=0, dash=True),
                                      dict(name="10⁻¹³ m², injection above", data=fser("K 1e-13, inject above"), c=1),
                                      dict(name="10⁻¹³ m², injection below", data=fser("K 1e-13, inject below"), c=1, dash=True)])
rowsB = []
for k, lb in (("no coupling, inject above", "Flow not affected by temperature"), ("element (unbounded)", "The same, dipole element"), ("K 1e-14, inject above", "10⁻¹⁴ m², injection above"),
              ("K 1e-14, inject below", "10⁻¹⁴ m², injection below"), ("K 1e-13, inject above", "10⁻¹³ m², injection above"), ("K 1e-13, inject below", "10⁻¹³ m², injection below"),
              ("K 1e-12, inject above", "10⁻¹² m², injection above"), ("K 1e-12, inject below", "10⁻¹² m², injection below")):
    r = fb[k]
    rowsB.append([lb, f(r.get("t_99_yr"), 1) if r.get("t_99_yr") else "none in 30 yr"] + [f(100 * r["theta"][y], 0) + " %" for y in ("5", "10", "20", "30")])
T_updown = table(["Vertical permeability and arrangement", "1 % of the rise lost after (yr)", "Rise kept at 5 yr", "10 yr", "20 yr", "30 yr"], rowsB, cls="num wrapfirst",
                 note="Zones 150 m apart (screens 20 m long) in a block of 300 m radius and 600 m height, 183 °C, 10 kg/s of 60 °C water at 47 MPa, densely fractured, closed outer boundaries. "
                      "The cold water is 10 % denser and three times as viscous as the hot.")
C.update(fb_a13=fb["K 1e-13, inject above"]["theta"]["30"], fb_b13=fb["K 1e-13, inject below"]["theta"]["30"], fb_a14=fb["K 1e-14, inject above"]["theta"]["30"],
         fb_b14=fb["K 1e-14, inject below"]["theta"]["30"], fb_nc=fb["no coupling, inject above"]["theta"]["30"], fb_el=fb["element (unbounded)"]["theta"]["30"],
         fb_a12_2=fb["K 1e-12, inject above"]["theta"]["2"], fb_a13_first=fb["K 1e-13, inject above"]["t_99_yr"], fb_b13_first=fb["K 1e-13, inject below"]["t_99_yr"])

# =============================================================================================== 9. checks
p200, p400 = chk["piston"]["200"], chk["piston"]["400"]
dp = chk["dipole"]
gr = chk["grid"]
rowsC = [["Energy balance of every run", "heat leaving with the water against heat lost by the block", "closes to 10⁻¹³ of the heat produced"],
         ["Piston front in a thin shell, no conduction", "front at the sweep time; width of the front from the scheme alone",
          f"front at {f(p200['front_over_sweep_time'], 3)} of the sweep time; 10–90 % width {f(100 * p200['width_10_90_over_sweep_time'], 1)} % of the travel time at 200 cells, "
          f"{f(100 * p400['width_10_90_over_sweep_time'], 1)} % at 400"],
         ["Streamline table of the dipole element", "streamlines integrated again (Runge–Kutta, 10⁻¹⁰)", "within 3 × 10⁻⁵" if abs(chk['streamlines']['max_rel'] - 3e-5) < 5e-6 else f"within {chk['streamlines']['max_rel']:.1e}"],
         ["Block model against the dipole element", "screens 1 m apart in a block of 3 m radius and 6 m height, no conduction",
          f"share of the rise lost within {f(dp['max_diff'], 3)} from x = 0.6 to 4 (at x = 2: {f(dp['block'][5], 3)} against {f(dp['element'][5], 3)})"],
         ["Grid of the sandbox model", "20 × 60, 28 × 120 (used), 40 × 180 and 56 × 240 cells",
          f"produced temperature within {f(gr['28x120']['max_diff_from_finest_K'], 2)} K of the finest grid; plateau {f(gr['28x120']['onset_1K_min'], 1)} against {f(gr['56x240']['onset_1K_min'], 1)} min"],
         ["Buoyancy term", "slender cold column in a large hot block", "sinks at K Δρ g / μ within 3 %"],
         ["Type 11 well with the benchmark's element", "the benchmark result run again through the new script", f"{f(C['b_T1'], 1)} → {f(C['b_T30'], 1)} °C, as in the benchmark (135 → 57 °C)"]]
T_checks = table(["Check", "Against", "Result"], rowsC, cls="small wrapfirst")

# reported numbers and their reading
rowsRep = [["Plateau of the production temperature", "about 50 min at 0.3 m³/h, 25 min at 0.6 m³/h; the first 30 min in the nominal case", "taken as read; the three are the same flushed share"],
           ["Decline after the plateau, nominal case", "about 0.8 °C/min, to 37 °C after 30 min; heat flow 22 kW falling to about 10 kW after 70 min", "10 kW is 38 °C at 0.5 m³/h"],
           ["Decline after 60 min, flow series", "8.3, 17.3, 23.3, 29.7 and 32.2 % for 0.3–0.6 m³/h", "five values for four flows; both alignments are shown below"],
           ["Mean production temperature and heat flow, flow series", "59.0 °C and 8.9 kW at 0.2 m³/h; 52.8 °C and 22.3 kW at 0.6 m³/h",
            "pack at 60 °C, injection 21 °C (fixed by the heat flows); averaged over 60 min, as stated for the temperature series"],
           ["Initial temperature of the pack, 0.5 m³/h, 21 °C", "50 / 60 / 70 °C: means 47.72 / 55.29 / 63.44 °C, 15.59 / 20.00 / 24.76 kW, declines 20 / 29.67 / 33.43 %",
            "decline read as a share of the initial temperature in °C: 40.0 / 42.2 / 46.6 °C at 60 min"],
           ["Injection temperature, 0.5 m³/h, 60 °C", "21 → 41 °C: mean production temperature +2.3 K, mean heat flow −10.33 kW", "the same three decline percentages are printed as for the previous series; not used"],
           ["Screen spacing", "90 → 120 cm: mean production temperature +4.4 K, heat flow +2.6 kW, decline 31.5 % → 8.0 %", "flow and temperatures of this pair are not stated; nominal values assumed"],
           ["Temperature field", "vertical line: first cooling at 50 cm depth after 30 min, moving down at about 1.0 cm/min, stable above 30 cm and below 110 cm; line across the box: "
            "20 K drop over 80 % of the distance at 40 min, spreading at about 0.75 cm/min", "the position of the vertical line differs between methods and results"],
           ["The authors' own model", "within 4 % overall, 8 % at most, early in the test; the experiment declines slightly later", "Darcy flow, local thermal equilibrium"]]
T_rep = table(["Quantity", "Reported", "Reading made here"], rowsRep, cls="small wrapfirst",
              note="From the text of the paper; its table of operating conditions could not be retrieved, so the default values of each series are inferred from the numbers.")

# =============================================================================================== text and assembly
import sg_report_text as TXT  # noqa: E402
ctx = dict(date=today, T_cmp=T_cmp, T_decay=T_decay, T_fronts=T_fronts, T_sens=T_sens, T_rec=T_rec, T_space=T_space, T_t11=T_t11, T_updown=T_updown, T_checks=T_checks, T_rep=T_rep)
ctx["F_typ"] = f'<div class="figrow">{fig("typ", "Nominal case: produced water against time", "0.5 m³/h, pack at 60 °C, injection at 21 °C, screens 90 cm apart", "segs_sandbox.py lab, calib", "tall")}' \
               f'{fig("onset", "Time on the plateau against flow rate", "Pack at 60 °C, injection at 21 °C, screens 90 cm apart; model: until 1 K is lost", "segs_sandbox.py lab, calib", "tall")}</div>'
ctx["F_field"] = f'<div class="figrow">{fig("field30", "The sandbox after 30 minutes", "Section through the axis, fitted model, nominal case; the well is in the middle", "segs_sandbox.py calib", "field")}' \
                 f'{fig("field60", "After 60 minutes", "Same section; the cold water has reached the production screen", "segs_sandbox.py calib", "field")}</div>'
ctx["F_rec"] = f'<div class="figrow one">{fig("rec", "Three settings run for three hours", "Share of the temperature rise kept by the produced water; fitted model, pack at 60 °C", "segs_sandbox.py rec")}</div>'
ctx["F_curves"] = f'<div class="figrow">{fig("curves", "One well, two screens: the decline in dimensionless form", "D: spacing of the screens; Q: flow; (ρc)b: heat capacity of the swept medium", "segs_sandbox.py field", "tall")}' \
                  f'{fig("t11", "The type 11 well with four descriptions of the rock at its toe", "Water leaving the rock, 10 kg/s entering at 60 °C, rock at 183 °C", "segs_sandbox.py field", "tall")}</div>'
ctx["F_updown"] = f'<div class="figrow">{fig("updown_lab", "Sandbox: injection above or below", "Nominal case, model only; the experiment injected through the upper screen", "segs_sandbox.py lab, calib", "tall")}' \
                  f'{fig("updown_field", "Field block: injection above or below", "Zones 150 m apart, 10 kg/s, 183 °C rock; by vertical permeability", "segs_sandbox.py field", "tall")}</div>'
texts = TXT.build(C, f, REP)
ctx.update(texts)

tpl = open(os.path.join(HERE, "dt_template.html")).read()
style = re.search(r"<style>.*?</style>", tpl, re.S).group(0)
head = tpl[: tpl.index("<style>")].replace("Twenty Geothermal Development Types", "Single-Well EGS Sandbox Check")
style = re.sub(r"<style>\s*/\*.*?\*/", "<style>\n/* Layout: one reading column (70ch prose) with full-width figure rows and tables. Same identity as the other aharrsim reports. */", style, count=1, flags=re.S)
style = style.replace("</style>", TXT.EXTRA_CSS + "</style>")
html = head + style + TXT.PAGE
for _ in range(3):
    for k, v in ctx.items():
        html = html.replace("{{" + k + "}}", str(v))


def clean(o):
    if isinstance(o, dict):
        return {k: clean(v) for k, v in o.items() if v is not None}
    if isinstance(o, list):
        return [clean(v) for v in o]
    if isinstance(o, (np.floating, float)):
        return None if not math.isfinite(o) else round(float(o), 4)
    if isinstance(o, np.integer):
        return int(o)
    return o


html = html.replace("{{data}}", json.dumps(clean(dict(charts=charts)), separators=(",", ":"), ensure_ascii=False)).replace("{{js}}", open(os.path.join(HERE, "sg_report.js")).read())
left = sorted(set(re.findall(r"\{\{([A-Za-z0-9_]+)\}\}", html)))
if left:
    print("unresolved:", left)
out = os.path.join(HERE, "aharrsim_segs_report.html")
open(out, "w", encoding="utf8").write(html)
print(f"wrote {out} {os.path.getsize(out) / 1e3:.0f} kB, {len(charts)} charts")
print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in C.items() if not isinstance(v, (list, dict))})
