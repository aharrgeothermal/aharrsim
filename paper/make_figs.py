"""Figures 1-11 of the Geothermics manuscript, drawn from validation/results/*.json.

    python make_figs.py            all figures
    python make_figs.py 5 9 11     selected figures
"""
import json
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
# the package root: the first directory above this file that holds validation/results, else AHARRSIM_ROOT
ROOT = next((d for d in (os.path.dirname(HERE), os.path.dirname(os.path.dirname(HERE)))
             if os.path.isdir(os.path.join(d, "validation", "results"))), os.environ.get("AHARRSIM_ROOT", os.getcwd()))
RES = ROOT + "/validation/results"
OUT = os.path.join(HERE, "figs")
os.makedirs(OUT, exist_ok=True)
sys.path.insert(0, ROOT + "/validation"); sys.path.insert(0, ROOT)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8.5, "axes.titlesize": 9, "axes.labelsize": 8.5, "legend.fontsize": 7.2,
                     "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "axes.linewidth": 0.7, "lines.linewidth": 1.3, "figure.dpi": 110,
                     "axes.grid": True, "grid.color": "#dddddd", "grid.linewidth": 0.5, "legend.frameon": False, "axes.spines.top": False,
                     "axes.spines.right": False})
C = ["#1f77b4", "#d62728", "#2ca02c", "#ff7f0e", "#9467bd", "#8c564b", "#17becf", "#7f7f7f"]
YEAR = 365.25 * 86400.0


def L(n):
    with open(os.path.join(RES, n + ".json")) as f:
        return json.load(f)


lab = lambda t: {3.5: "3b", 14.5: "14b"}.get(t, str(int(t)))
SHORT = {1: "HSA doublet", 2: "Pumped MT well", 3: "HT flash well", 3.5: "Dry-steam well", 4: "Waste-heat ORC", 5: "GSHP U-tube", 6: "400 m coaxial BHE",
         7: "Coaxial retrofit", 8: "Multilateral loop", 9: "Closed-loop doublet", 10: "EGS doublet", 11: "Single-well EGS", 12: "Cyclic fracture",
         13: "20 cyclic wells", 14: "Project Red type", 14.5: "Cape Station type", 15: "Co-produced water", 16: "Abandoned-well loop", 17: "SAGD heat store",
         18: "Supercritical well", 19: "Superhot EGS", 20: "20 km loop"}


def save(fig, name):
    fig.savefig(os.path.join(OUT, name), dpi=300, bbox_inches="tight", pad_inches=0.04); plt.close(fig); print("wrote", name)


# ------------------------------------------------------------------ Fig. 1
def fig1():
    R = {float(r["id"]): r for r in L("DT_types")}
    fig, axs = plt.subplots(3, 2, figsize=(7.0, 7.6)); axs = axs.ravel()

    def panel_log(ax, ids, title, dashed_res=False):
        for i, t in enumerate(ids):
            s = R[t]["series"]
            ax.plot(s["t_yr"], s["T_C"], color=C[i % len(C)], label=f"{lab(t)} {SHORT[t]}")
            if dashed_res and s.get("T_res_C"):
                ax.plot(s["t_yr"], s["T_res_C"], color=C[i % len(C)], ls="--", lw=0.9)
        ax.set_xscale("log"); ax.set_xlim(0.005, 30); ax.set_xlabel("Time (years)"); ax.set_ylabel("Produced water (°C)"); ax.set_title(title, loc="left")
        ax.legend(loc="best")

    panel_log(axs[0], [1, 2, 3, 3.5, 15], "(a) Conventional and co-produced water")
    ax = axs[1]
    s5 = R[5]["series"]; t5 = np.array(s5["t_yr"]); T5 = np.array(s5["T_C"]); on = np.array(s5["on"]) == 1
    ax.plot(np.where(on, t5, np.nan), np.where(on, T5, np.nan), color=C[0], label="5 GSHP U-tube, heating seasons")
    s6 = R[6]["series"]; ax.plot(s6["t_yr"], s6["T_C"], color=C[1], label="6 400 m coaxial BHE")
    ax.set_xlim(0, 10); ax.set_xlabel("Time (years)"); ax.set_ylabel("Loop outlet (°C)"); ax.set_title("(b) Shallow exchangers", loc="left"); ax.legend(loc="best")
    panel_log(axs[2], [7, 8, 9, 16, 20], "(c) Deep closed loops (7 at 10 kg/s, 8 at 100, 20 at 20)")
    panel_log(axs[3], [10, 11, 14, 14.5, 18, 19], "(d) Fractured and supercritical (dashed: reservoir outlet)", dashed_res=True)
    axs[3].legend(loc="lower left", fontsize=6.3)
    ax = axs[4]; s = R[12]["series"]
    ax.plot(s["day"], s["T_wellhead_C"], color=C[0], label="wellhead, daily mean")
    ax.plot(s["day"], s["T_fracture_C"], color=C[1], label="leaving the fracture, daily mean")
    ax.plot(s["day"], s["T_first_C"], color=C[2], ls="--", lw=0.9, label="first of the day")
    ax.plot(s["day"], s["T_last_C"], color=C[3], ls="--", lw=0.9, label="last of the day")
    ax.set_xlim(0, 365); ax.set_xlabel("Day of the first year"); ax.set_ylabel("Produced water (°C)"); ax.set_title("(e) Cyclic single fracture (12)", loc="left"); ax.legend(loc="best")
    ax = axs[5]; s = R[17]["series"]
    ax.plot(s["t_yr"], s["T_C"], color=C[0], label="17 chamber outlet")
    ax.axhline(50, color="#888", ls=":", lw=0.9); ax.text(14.8, 52, "50 °C return", ha="right", fontsize=7, color="#555")
    ax.set_xlim(0, 15); ax.set_xlabel("Time (years)"); ax.set_ylabel("Produced water (°C)"); ax.set_title("(f) SAGD heat store (17)", loc="left"); ax.legend(loc="best")
    fig.tight_layout(h_pad=1.2, w_pad=1.0); save(fig, "fig01_temperature.png")


# ------------------------------------------------------------------ Fig. 2
def fig2():
    S = {float(r["id"]): r for r in L("DT_summary")}
    fig, ax = plt.subplots(figsize=(5.4, 4.0))
    for lo, hi, col, name in ((0, 80, "#f2f2f2", "heat only"), (80, 120, "#e6eef7", "binary, 4–7 % of heat"), (120, 200, "#d4e1f1", "binary power"), (200, 350, "#f7e3dc", "flash and steam")):
        ax.axvspan(lo, hi, color=col, lw=0, zorder=0); ax.text((lo + hi) / 2, 70, name, ha="center", fontsize=7, color="#444")
    off = {1: 4, 7: 4, 9: -4, 16: 4, 11: -4, 12: -4, 13: 4, 17: 4, 5: 4, 6: -4, 14: 4, 14.5: -4, 3: -4, 3.5: 4}
    for t, s in S.items():
        if s["T_last_C"] is None:
            continue
        y = max(s["net_last_MWe"] or 0, 0.012)
        hollow = not (s["net_last_MWe"] and s["net_last_MWe"] > 0)
        if hollow:                                   # heat-only cases share the axis floor: stagger them a little so the labels separate
            y *= {1: 1.0, 7: 1.45, 9: 1.0, 16: 1.0, 11: 1.45, 12: 2.1, 14: 2.1, 14.5: 3.0, 17: 1.0, 5: 1.0, 6: 1.45}.get(t, 1.0)
        ax.scatter(s["T_last_C"], y, s=36, facecolor="white" if hollow else C[0], edgecolor=C[0], lw=1.2, zorder=3)
        o = off.get(t, 4)
        ax.annotate(lab(t), (s["T_last_C"], y), xytext=(o, 0), textcoords="offset points", fontsize=7, va="center", ha="left" if o > 0 else "right")
    ax.set_yscale("log"); ax.set_ylim(0.009, 100); ax.set_xlim(0, 340)
    ax.set_xlabel("Produced water temperature at the horizon (°C)"); ax.set_ylabel("Net electricity at the horizon (MWe)")
    ax.text(335, 0.0105, "hollow: no net electricity at the horizon (drawn near the axis floor)", fontsize=6.5, color="#555", va="bottom", ha="right")
    fig.tight_layout(); save(fig, "fig02_scatter.png")


# ------------------------------------------------------------------ Fig. 3
def fig3():
    S = {float(r["id"]): r for r in L("DT_summary")}
    QUOTED = {1: (0.0, 2.0), 2: (2.0, 7.0), 3: (3.0, 25.0), 3.5: (3.0, 25.0), 4: (0.355, 2.0), 7: (2.0, 9.0), 8: (2.05, 2.2), 10: (3.0, 5.0), 13: (50.0, 50.0),
              14: (3.5, 3.5), 14.5: (10.0, 16.0), 15: (0.02, 0.35), 18: (30.0, 30.0), 19: (50.0, 50.0), 20: (100.0, 100.0)}
    LF = {r["m"]: r for r in L("SC_loop_flow")}
    rows = [t for t in S if t in QUOTED]
    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    for i, t in enumerate(rows):
        lo, hi = QUOTED[t]; lo = max(lo, 0.01)
        if hi / lo < 1.15:                           # a single quoted figure: draw it as a short bar
            lo, hi = lo * 0.93, hi * 1.07
        ax.plot([lo, hi], [i, i], color="#bbb", lw=7, solid_capstyle="butt", zorder=1)
        g1, gl = S[t]["gross_first_MWe"] or 0, S[t]["gross_last_MWe"] or 0
        ax.scatter(max(g1, 0.01), i, s=46, facecolor="white", edgecolor=C[0], lw=1.3, zorder=3)
        ax.scatter(max(gl, 0.01), i, s=30, color=C[1], zorder=4)
        if t == 8 and 30.0 in LF:                    # the same loop at the flow that suits electricity
            ax.scatter(LF[30.0]["MWe"][0], i + 0.28, s=30, marker="^", facecolor="white", edgecolor=C[0], lw=1.1, zorder=3)
            ax.scatter(LF[30.0]["MWe"][2], i + 0.28, s=24, marker="^", color=C[1], zorder=4)
    ax.set_yticks(range(len(rows))); ax.set_yticklabels([f"{lab(t)} {SHORT[t]}" for t in rows], fontsize=7.5)
    ax.set_xscale("log"); ax.set_xlim(0.008, 200); ax.set_xlabel("Gross electricity (MWe), log scale"); ax.invert_yaxis()
    fig.legend(handles=[Line2D([], [], color="#bbb", lw=7, label="quoted"), Line2D([], [], marker="o", mfc="white", mec=C[0], ls="", label="model, year 1"),
                        Line2D([], [], marker="o", color=C[1], ls="", label="model, horizon"),
                        Line2D([], [], marker="^", mfc="white", mec=C[0], ls="", label="type 8 at 30 kg/s (year 1; filled: horizon)")], loc="lower center", ncol=4, fontsize=6.8,
               columnspacing=1.2, handletextpad=0.4)
    fig.tight_layout(rect=(0, 0.05, 1, 1)); save(fig, "fig03_claims.png")


# ------------------------------------------------------------------ Fig. 4
def fig4():
    V9 = L("V9_single_well"); V9b = L("V9b_twophase_production_mode"); B4 = L("B4_geocluster_utube"); V11 = L("V11_plant"); V10 = L("V10_reservoir"); CF = L("ST2_confirm")
    fig, axs = plt.subplots(3, 2, figsize=(7.0, 8.4)); axs = axs.ravel()
    ax = axs[0]; rr = V9["ramey"]["rows"]
    ax.plot([r["t_days"] for r in rr], [r["T_wh_model"] for r in rr], "o-", color=C[0], label="aharrsim, IAPWS water")
    ax.plot([r["t_days"] for r in rr], [r["T_wh_model_constant_fluid"] for r in rr], "s-", color=C[2], label="aharrsim, constant properties")
    ax.plot([r["t_days"] for r in rr], [r["T_wh_ramey"] for r in rr], "D--", color="#555", label="Ramey / Hasan–Kabir")
    ax.set_xscale("log"); ax.set_xlabel("Time (days)"); ax.set_ylabel("Wellhead temperature (°C)"); ax.set_title("(a) 3 km producer, 20 kg/s, 160 °C feed", loc="left"); ax.legend()
    ax = axs[1]
    for w, col in ((V9b[0], C[0]), (V9b[4], C[1])):
        ax.plot(np.array(w["profile"]["p_bar"]) / 10.0, w["profile"]["md"], color=col, label=f"{w['name']}, model")
        ax.scatter(np.array(w["measured"]["p_bar"]) / 10.0, w["measured"]["md"], s=14, facecolor="white", edgecolor=col, lw=1, zorder=3, label=f"{w['name']}, measured")
    ax.invert_yaxis(); ax.set_xlabel("Pressure (MPa)"); ax.set_ylabel("Depth (m)"); ax.set_title("(b) Flowing two-phase wells", loc="left"); ax.legend(fontsize=6.5)
    ax = axs[2]
    for i, c in enumerate(B4):
        comp = np.array(c["comparison"])
        ax.plot(comp[:, 0], comp[:, 1], color=C[i % 7], label=f"{c['L1'] / 1000:g}+{c['L2'] / 1000:g} km, {c['m']:g} kg/s, {c['grad'] * 1000:.0f} °C/km")
        ax.scatter(comp[:, 0], comp[:, 2], s=16, facecolor="white", edgecolor=C[i % 7], lw=1, zorder=3)
    ax.set_xscale("log"); ax.set_xlabel("Time (years)"); ax.set_ylabel("Outlet temperature (°C)"); ax.set_title("(c) U-loop chains against GeoCLUSTER (markers)", loc="left"); ax.set_ylim(22, 165); ax.legend(fontsize=5.6, ncol=2, loc="lower left", columnspacing=0.8)
    ax = axs[3]
    b = [r for r in V11 if r["kind"] == "binary"]; f35 = [r for r in V11 if r["kind"] == "flash, condenser 35 C"]; f40 = [r for r in V11 if r["kind"] == "flash, condenser 40 C"]
    ax.plot([r["T_in_C"] for r in b], [r["eta_u_cycle"] for r in b], "o-", color=C[0], label="ORC cycle model")
    ax.plot([r["T_in_C"] for r in b], [r["eta_u_geophires_sub"] for r in b], "D--", color="#555", label="GEOPHIRES binary, subcritical")
    ax.plot([r["T_in_C"] for r in b], [r["eta_u_geophires_super"] for r in b], "D:", color="#999", label="GEOPHIRES binary, supercritical")
    ax.plot([r["T_in_C"] for r in f35], [r["eta_u_double"] for r in f35], "o-", color=C[1], label="double flash, 35 °C condenser")
    ax.plot([r["T_in_C"] for r in f40], [r["eta_u_double"] for r in f40], "o-", color=C[3], label="double flash, 40 °C condenser")
    ax.plot([r["T_in_C"] for r in f35], [r["eta_u_geophires_double"] for r in f35], "D--", color=C[4], label="GEOPHIRES double flash")
    ax.set_xlabel("Inlet temperature (°C)"); ax.set_ylabel("Utilisation efficiency"); ax.set_title("(d) Plant models against GEOPHIRES", loc="left"); ax.legend(fontsize=6)
    ax = axs[4]; sf = V10["sharp_front"]
    ax.plot(sf["x"], sf["mpmath"], "-", color="#555", lw=2.2, alpha=0.5, label="30-digit reference (mpmath)")
    ax.plot(sf["x"], sf["dehoog"], "o", color=C[0], ms=4, label="de Hoog et al. (used here)")
    ax.plot(sf["x"], sf["stehfest"], "s--", color=C[1], ms=3.5, lw=1, label="Gaver–Stehfest, N = 14")
    ax.set_xlabel("x = t / t$_s$"); ax.set_ylabel("Outlet temperature drop, T$_D$"); ax.set_title("(e) Fracture solution at a sharp front", loc="left"); ax.legend(fontsize=6.2, loc="upper left")
    ax = axs[5]
    for sid, col, nm in ((10, C[0], "10"), (14, C[1], "14"), (141, C[2], "14b"), (19, C[3], "19")):
        xs = [g for c in CF if c.get("sid") == sid and "failed" not in c for g in c["gross"]]
        ys = [g for c in CF if c.get("sid") == sid and "failed" not in c for g in c["screen_gross"]]
        ax.scatter(np.maximum(xs, 0.05), np.maximum(ys, 0.05), s=14, facecolor="white", edgecolor=col, lw=1, label=f"type {nm}", zorder=3)
    ax.plot([0.04, 100], [0.04, 100], color="#555", lw=0.8, ls="--")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(0.04, 100); ax.set_ylim(0.04, 100)
    ax.set_xlabel("Gross electricity, full chain (MWe)"); ax.set_ylabel("Gross electricity, screening form (MWe)"); ax.set_title("(f) Screening form against the full chain", loc="left"); ax.legend(fontsize=6.5)
    fig.tight_layout(h_pad=1.2, w_pad=1.0); save(fig, "fig04_verification.png")


# ------------------------------------------------------------------ Fig. 5
def fig5():
    G = L("GX_compare")["rows"]
    fig, axs = plt.subplots(1, 3, figsize=(7.4, 2.9))
    mk = ["o", "s", "D", "^", "v", "P", "X", "*", "<", ">"]
    ax = axs[0]
    for i, r in enumerate(G):
        ax.scatter(r["gx_T"], r["coax_T"], s=22, marker=mk[i], facecolor="white", edgecolor=C[i % 8], lw=1.1, zorder=3, label=f"{r['id']}")
    ax.plot([40, 380], [40, 380], color="#555", lw=0.8, ls="--")
    ax.text(385, 45, "3 and 19: GEOPHIRES passes\nliquid at reservoir temperature\nto the flash plant; aharrsim gives\nthe two-phase wellhead state", fontsize=5.8, color="#444", ha="right", va="bottom")
    ax.set_xlabel("GEOPHIRES (°C)"); ax.set_ylabel("aharrsim (°C)"); ax.set_title("(a) Produced temperature,\ntwo times per case", loc="left"); ax.legend(ncol=2, fontsize=6, loc="upper left", columnspacing=0.6, handletextpad=0.1)
    ax = axs[1]
    for r, col, nm in ((next(g for g in G if g["id"] == "14"), C[0], "14 Project Red type"), (next(g for g in G if g["id"] == "14b"), C[1], "14b Cape Station type")):
        h = r["hist"]
        ax.plot(h["t"], h["coax_T"], color=col, label=nm)
        ax.scatter(h["t"], h["gx_T"], s=12, facecolor="white", edgecolor=col, lw=0.9, zorder=3)
    ax.set_xlim(0, 20); ax.set_xlabel("Time (years)"); ax.set_ylabel("Produced water (°C)"); ax.set_title("(b) Breakthrough: lines aharrsim,\nmarkers GEOPHIRES", loc="left"); ax.legend(fontsize=6.3, loc="upper right")
    ax = axs[2]
    for i, r in enumerate(G):
        ax.scatter(np.maximum(r["gx_out"], 0.03), np.maximum(r["coax_out"], 0.03), s=22, marker=mk[i], facecolor="white", edgecolor=C[i % 8], lw=1.1, zorder=3)
    ax.plot([0.02, 20], [0.02, 20], color="#555", lw=0.8, ls="--")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(0.02, 20); ax.set_ylim(0.02, 20)
    ax.set_xlabel("GEOPHIRES (MW)"); ax.set_ylabel("aharrsim (MW)"); ax.set_title("(c) Gross electricity\n(heat for case 1)", loc="left")
    fig.tight_layout(w_pad=0.8); save(fig, "fig05_geophires.png")


def tornado(ax, rows, base, band, xlabel, title, quoted=None):
    rows = sorted(rows, key=lambda r: abs(r[2] - r[1]))
    for i, (lb, lo, hi, des) in enumerate(rows):
        for v in (lo, hi):
            ax.barh(i, v - base, left=base, height=0.62, color=C[0] if v < base else C[1], alpha=0.45 if des else 0.95, lw=0)
    ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[0] for r in rows])
    if band is not None:
        ax.axvspan(band[0], band[2], color="#999999", alpha=0.18, lw=0); ax.axvline(band[1], color="#555555", lw=0.8, ls=":")
    ax.axvline(base, color="k", lw=0.8, ls="--")
    if quoted:
        for q in quoted:
            ax.axvline(q, color="k", lw=0.8, ls=(0, (1, 2)))
    ax.set_xlabel(xlabel); ax.set_title(title, loc="left"); ax.grid(axis="y", visible=False)


# ------------------------------------------------------------------ Fig. 6
def fig6():
    import study_cl as CL
    by = {(r["case"], r["tag"]): r for r in L("ST1_oat") if "failed" not in r and r["case"] != 8}
    by.update({(r["case"], r["tag"]): r for r in L("ST1_oat_8_m30") if "failed" not in r})          # type 8 at the flow that suits electricity (30 kg/s)
    MC = [r for r in L("ST1_mc") if r["case"] != 8] + L("ST1_mc_8_m30")
    CL.OAT[8] = [("mdot", 20.0, 50.0, "Flow 20 / 50 kg/s", "des") if o[0] == "mdot" else o for o in CL.OAT[8]]
    fig, axs = plt.subplots(1, 3, figsize=(7.4, 3.0))
    for ax, (c, key, xl, ttl) in zip(axs, ((7, "MWth", "Heat (MWth)", "(a) Type 7"), (8, "MWe", "Gross electricity (MWe)", "(b) Type 8 at 30 kg/s"),
                                           (20, "MWe", "Gross electricity (MWe)", "(c) Type 20"))):
        rows = []
        for var, lo, hi, label, kind in CL.OAT[c]:
            a, b = by.get((c, f"oat|{var}|low")), by.get((c, f"oat|{var}|high"))
            if a and b:
                rows.append((label.replace("Rock c", "C").replace("Rock h", "H").replace("Tubing insulation k", "Insulation").replace(" W/m/K", "").replace("TD temperature", "T at depth"),
                             a["30"][key], b["30"][key], kind == "des"))
        kt = by[(c, "phys|kT|True")]["30"][key]; base = by[(c, "base")]["30"][key]
        rows.append(("k falling with T", kt, base, False))
        mc = np.percentile([r["30"][key] for r in MC if r["case"] == c and "failed" not in r], [10, 50, 90])
        tornado(ax, rows, base, mc, xl, ttl)
        ax.tick_params(axis="y", labelsize=6.3)
    fig.tight_layout(); save(fig, "fig06_closed_loop_uncertainty.png")


# ------------------------------------------------------------------ Fig. 7
def fig7():
    import study_open as OP
    S = L("ST2_screen")["systems"]
    fig, axs = plt.subplots(2, 4, figsize=(7.5, 4.7))
    for j, sid in enumerate(("10", "14", "141", "19")):
        o = S[sid]; s = OP.SYS[int(sid)]; q = sorted(set(s["quoted"])); nm = s["label"]
        ax = axs[0, j]
        for i, d in enumerate(o["density"]):
            ax.semilogx([r["N"] for r in d["rows"]], [r["mean_net"] for r in d["rows"]], "o-", ms=2.3, color=C[i], label=f"τ = {d['tau']:.0f}")
        ax.set_title(f"({'abcd'[j]}) Type {nm}: {s['m']:.0f} kg/s,\n{o['base']['V_km3']:.3g} km³", loc="left", fontsize=7.8); ax.set_xlabel("Number of fractures"); ax.set_ylim(bottom=0)
        if j == 0:
            ax.set_ylabel("30-year mean net (MWe)"); ax.legend(fontsize=6.5)
        ax = axs[1, j]
        r0 = o["flow"][0]["rows"]; ax.plot([r["m"] for r in r0], [r["mean_net"] for r in r0], "k--", lw=1, label="example fractures, volume ×1")
        for i, v in enumerate([v for v in o["volume"] if v["vf"] in (1, 2, 4, 8, 16)]):
            ax.plot([r["m"] for r in v["rows"]], [r["mean_net"] for r in v["rows"]], "o-", ms=2.3, color=C[i], label=f"volume ×{v['vf']}")
        for qq in q:
            ax.axhline(qq, color="k", lw=0.8, ls=(0, (1, 2)))
        ax.set_title(f"({'efgh'[j]}) Type {nm}: fractures\n10 m apart", loc="left", fontsize=7.8); ax.set_xlabel("Flow (kg/s)"); ax.set_ylim(bottom=0)
        if j == 0:
            ax.set_ylabel("30-year mean net (MWe)")
    h, l = axs[1, 0].get_legend_handles_labels(); fig.legend(h, l, loc="lower center", ncol=6, fontsize=6.8)
    fig.tight_layout(rect=(0, 0.05, 1, 1), w_pad=0.5); save(fig, "fig07_fractured_sweeps.png")


# ------------------------------------------------------------------ Fig. 8
def fig8():
    import study_open as OP
    P3 = L("ST3_open_uncertainty")["systems"]
    fig, axs = plt.subplots(2, 2, figsize=(7.2, 5.2)); axs = axs.ravel()
    SH = {"T_rock": "Gradient ±20 %", "k": "Conductivity 2.0 / 4.0", "kT": "k falling with T", "rhoc": "Heat capacity ±15 %", "f": "Active fractures 20 %", "phi": "Swept volume 50 %",
          "I0": "Impedance ×2 / ×0.5", "T_amb": "Ambient 45 °C"}
    for ax, (sid, i, ttl) in zip(axs, (("10", 1, "(a) Type 10: 100 fractures, 80 kg/s"), ("14", 2, "(b) Type 14: 16× volume, 80 kg/s"),
                                       ("141", 2, "(c) Type 14b: 8× volume, 200 kg/s"), ("19", 2, "(d) Type 19: 8× volume, 200 kg/s"))):
        d = P3[sid]["designs"][i]; base = d["base"]["mean_net"]
        rows = []
        for t in d["tornado"]:
            lo, hi = t["low"]["mean_net"], t["high"]["mean_net"]
            if t["var"] in ("kT", "f", "phi", "T_amb"):
                lo, hi = min(lo, hi), base
            rows.append((SH.get(t["var"], t["label"]), lo, hi, False))
        tornado(ax, rows, base, d["mc"]["15"]["mean_net"], "Mean net electricity (MWe)", ttl, quoted=sorted(set(OP.SYS[int(sid)]["quoted"]))[:1])
        ax.tick_params(axis="y", labelsize=6.5)
    fig.tight_layout(); save(fig, "fig08_fractured_uncertainty.png")


# ------------------------------------------------------------------ Fig. 9
def fig9():
    SC = L("SC_screening")["swept"]; FE = L("FE_field_cases")["red"]
    X, Y = np.meshgrid(SC["x"], SC["y"]); TH = np.array(SC["theta"])
    fig, ax = plt.subplots(figsize=(5.9, 5.0))
    levels = [0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1.0001]
    cols = ["#3b6fb6", "#6f9bd1", "#a9c5e6", "#d9e6f2", "#f6e7c8", "#f3cf9a", "#eba56a", "#dc7245"]
    cf = ax.contourf(X, Y, TH, levels=levels, colors=cols, alpha=0.85)
    cs = ax.contour(X, Y, TH, levels=[0.1, 0.5, 0.9, 0.99], colors="#333333", linewidths=0.6)
    ax.clabel(cs, fmt={0.1: "10 %", 0.5: "50 %", 0.9: "90 %", 0.99: "99 %"}, fontsize=6.5, inline_spacing=2, manual=[(0.065, 2.6), (0.105, 4.0), (0.3, 4.6), (1.25, 5.6)])
    cb = fig.colorbar(cf, ax=ax, pad=0.02, ticks=[0, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99]); cb.set_label("Retained share of the temperature rise, θ", fontsize=7.5); cb.ax.tick_params(labelsize=6.5)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(0.05, 20); ax.set_ylim(0.05, 20)
    ax.set_xlabel("x = t / t_s   (share of the block's heat capacity flushed)"); ax.set_ylabel("y = (S/2) / √(αt)   (half spacing over conduction length)")
    ax.text(0.056, 12.5, "fractures cool alone:\nθ = 90 % on x·y = 0.43", fontsize=6.2, color="#222")
    ax.text(0.056, 0.058, "block cools as a whole:\nfront at x = 1", fontsize=6.2, color="#222")
    ax.axvline(1.0, color="#333", lw=0.7, ls=":")
    ex = {p["label"]: p for p in SC["points"] if p["kind"] == "example"}
    need = {}
    for p in SC["points"]:
        if p["kind"] == "need":
            need.setdefault(p["label"], p)
    offs = {"10": (-16, 5), "14": (5, 5), "14b": (5, -10), "19": (6, 4)}
    for lb, p in ex.items():
        ax.scatter(p["x"], p["y"], s=34, color="k", zorder=5)
        ax.annotate(lb, (p["x"], p["y"]), xytext=offs[lb], textcoords="offset points", fontsize=7.5, fontweight="bold", zorder=6)
        if lb in need:
            q = need[lb]
            ax.annotate("", xy=(q["x"], q["y"]), xytext=(p["x"], p["y"]), arrowprops=dict(arrowstyle="->", color="k", lw=0.8, shrinkA=3, shrinkB=3), zorder=4)
            ax.scatter(q["x"], q["y"], s=60, marker="*", facecolor="white", edgecolor="k", lw=0.9, zorder=5)
    # Project Red as operated (36 kg/s): path of the published block through time, its position at day 614, and the two readings of the early decline
    t = np.logspace(math.log10(0.3), math.log10(30), 40); ts = FE["t_sweep_yr"]; te = FE["t_end_yr"]
    xr = t / ts; yr = FE["y_end"] * np.sqrt(te / t)
    ax.plot(xr, yr, color=C[1], lw=1.1, ls="--", zorder=4)
    ax.scatter(FE["x_end"], FE["y_end"], s=38, marker="D", color=C[1], zorder=6)
    cf_, cv = FE["cases"]["fit_fractures"], FE["cases"]["fit_volume"]
    ax.scatter([cf_["x_end"], cv["x_end"]], [cf_["y_end"], cv["y_end"]], s=38, marker="D", facecolor="white", edgecolor=C[1], lw=1.1, zorder=6)
    ax.plot([cf_["x_end"], FE["x_end"], cv["x_end"]], [cf_["y_end"], FE["y_end"], cv["y_end"]], color=C[1], lw=0.7, zorder=5)
    ax.annotate("Project Red, day 614\n(published volume)", (FE["x_end"], FE["y_end"]), xytext=(-74, -30), textcoords="offset points", fontsize=6.3, color=C[1],
                arrowprops=dict(arrowstyle="-", color=C[1], lw=0.5))
    ax.annotate("45 % of the\nfractures active", (cf_["x_end"], cf_["y_end"]), xytext=(-58, 4), textcoords="offset points", fontsize=6.3, color=C[1])
    ax.annotate("half the\nvolume swept", (cv["x_end"], cv["y_end"]), xytext=(-14, 9), textcoords="offset points", fontsize=6.3, color=C[1])
    ax.text(xr[-1] * 1.0, yr[-1] * 1.22, "30 yr", fontsize=6.2, color=C[1], ha="center")
    ax.legend(handles=[Line2D([], [], marker="o", color="k", ls="", ms=5, label="example at 30 years"),
                       Line2D([], [], marker="*", mfc="white", mec="k", ls="", ms=8, label="design that holds the quoted capacity"),
                       Line2D([], [], marker="D", color=C[1], ls="--", ms=4.5, label="Project Red as operated, 36 kg/s")], loc="upper right", fontsize=6.5, frameon=True, framealpha=0.85)
    ax.grid(False)
    fig.tight_layout(); save(fig, "fig09_swept_map.png")


# ------------------------------------------------------------------ Fig. 10
def fig10():
    SC = L("SC_screening"); T = SC["touched"]; FE = L("FE_field_cases")["geretsried"]; LF = L("SC_loop_flow")
    S = {float(r["id"]): r for r in L("DT_summary")}
    fig, axs = plt.subplots(1, 2, figsize=(7.3, 3.3))
    ax = axs[0]
    N = np.linspace(0.01, 3.0, 200)
    ax.plot(N, 1 - np.exp(-N), color="k", lw=1.2, label="uniform rock temperature (laterals)")
    ax.plot(N, 1 - (1 - np.exp(-N)) / N, color="k", lw=1.2, ls="--", label="linear gradient, insulated return (coaxial)")
    mk = {"6": ("v", C[6]), "7": ("s", C[0]), "16": ("D", C[4]), "20": ("^", C[1]), "8": ("o", C[2]), "9": ("P", C[3])}
    seen = set()
    for r in T["check"]["rows"]:
        m, col = mk[r["type"]]
        ax.scatter(r["NTU"], r["eff"], s=20, marker=m, facecolor="white", edgecolor=col, lw=1.0, zorder=3, label=f"type {r['type']}" if r["type"] not in seen else None)
        seen.add(r["type"])
    for r in LF:                                         # type 8 at lower flows (full chain), on the lateral curve via the chart's NTU
        pass
    ax.set_xlim(0, 1.4); ax.set_ylim(0, 0.8); ax.set_xlabel("NTU = U′L / (ṁc)"); ax.set_ylabel("Effectiveness"); ax.set_title("(a) Simulated loops on the exchanger curves", loc="left")
    ax.legend(fontsize=6.2, loc="upper left")
    ax = axs[1]; yc = T["yield_curves"]
    for yv, ls, lbl in (("0.5", ":", "0.5 yr"), ("5", "--", "5 yr"), ("30", "-", "30 yr")):
        ax.plot(yc["T_rock"], yc["years"][yv], color="k", lw=1.1, ls=ls, label=f"chart, best flow, {lbl}")
    g1, gl = S[8.0]["gross_first_MWe"], S[8.0]["gross_last_MWe"]
    ax.scatter([262.5, 262.5], [g1 / 60 * 1e3, gl / 60 * 1e3], s=24, marker="o", facecolor=["white", C[2]], edgecolor=C[2], lw=1.1, zorder=3, label="type 8 at 100 kg/s (year 1, 30)")
    lf = {r["m"]: r for r in LF}
    if 30.0 in lf:
        ax.scatter([262.5, 262.5], [lf[30.0]["MWe"][0] / 60 * 1e3, lf[30.0]["MWe"][2] / 60 * 1e3], s=26, marker="^", facecolor=["white", C[2]], edgecolor=C[2], lw=1.1, zorder=3,
                   label="type 8 at 30 kg/s (year 1, 30)")
    rep = FE["reported"]
    ax.errorbar([160], [rep["MWe_reported"] / 30.0 * 1e3], yerr=[[rep["MWe_reported"] / 30.0 * 1e3 - 0.37 / 36.0 * 1e3], [rep["MWe_reported"] / 24.0 * 1e3 - rep["MWe_reported"] / 30.0 * 1e3]],
                xerr=10, fmt="D", ms=4.5, color=C[1], lw=0.9, capsize=2, zorder=4, label="Geretsried loop 1, reported (0.4 yr)")
    ax.scatter([160], [FE["nameplate_per_loop"]["MWe"] / 72.0 * 1e3], s=60, marker="*", facecolor="white", edgecolor=C[1], lw=1.0, zorder=4, label="Geretsried nameplate per loop")
    ax.set_xlim(100, 400); ax.set_ylim(0, 130); ax.set_xlabel("Rock temperature at the laterals (°C)"); ax.set_ylabel("Gross electricity per km of hot hole (kWe/km)")
    ax.set_title("(b) What a kilometre of sealed hole can give", loc="left"); ax.legend(fontsize=6.0, loc="upper left")
    fig.tight_layout(w_pad=1.0); save(fig, "fig10_touched_chart.png")


# ------------------------------------------------------------------ Fig. 11
def fig11():
    FE = L("FE_field_cases")["red"]
    fig, axs = plt.subplots(1, 2, figsize=(7.3, 3.2), gridspec_kw=dict(width_ratios=[1.45, 1]))
    ax = axs[0]
    t = np.array(FE["t_obs"]); T = np.array(FE["T_obs"]); used = np.array(FE["used"])
    ax.scatter(t[~used], T[~used], s=7, color="#bbbbbb", zorder=2, label="record, start-up and after shut-ins")
    ax.scatter(t[used], T[used], s=9, color="k", zorder=3, label="record used (digitised, ±1 K)")
    tt = np.array(FE["t_yr"])
    for key, col, ls, lbl in (("ideal", "#666666", "--", "published volume, 88 fractures"), ("fit_fractures", C[0], "-", "45 % of the fractures active"),
                              ("fit_volume", C[1], "-", "half the volume swept")):
        ax.plot(tt, FE["cases"][key]["T_C"], color=col, ls=ls, lw=1.2, label=lbl)
    ax.axvline(FE["t_end_yr"], color="#888", lw=0.7, ls=":"); ax.text(FE["t_end_yr"] - 0.06, 128, "day 614", fontsize=6.5, color="#555", rotation=90, ha="right", va="center")
    ax.axhline(FE["T_inj_C"], color="#888", lw=0.7, ls=":"); ax.text(7.9, FE["T_inj_C"] + 1.5, "injection", fontsize=6.5, color="#555", ha="right")
    ax.set_xlim(0, 8); ax.set_ylim(55, 185); ax.set_xlabel("Time on production (years)"); ax.set_ylabel("Flowing production temperature (°C)")
    ax.set_title("(a) Project Red: record and forecast at 36 kg/s", loc="left"); ax.legend(fontsize=6.3, loc="upper right")
    ax = axs[1]; g = FE["rms_grid"]
    Fg, Pg = np.meshgrid(g["phi"], g["f"]); Rm = np.array(g["rms_K"])
    cs = ax.contourf(Fg, Pg, Rm, levels=[0, 0.5, 1.0, 2.0, 4.0, 8.0, 100.0], colors=["#dc7245", "#eba56a", "#f3cf9a", "#f6e7c8", "#d9e6f2", "#a9c5e6"], alpha=0.9)
    cb = fig.colorbar(cs, ax=ax, pad=0.02, ticks=[0, 0.5, 1, 2, 4, 8]); cb.set_label("Misfit to the record (K, rms)", fontsize=7.5); cb.ax.tick_params(labelsize=6.5)
    ph = np.linspace(0.3, 1.2, 100)
    for prod, ls in ((FE["acceptable"]["f_phi"][0], ":"), (FE["acceptable"]["f_phi"][1], ":")):
        ax.plot(ph, prod / ph, color="k", lw=0.8, ls=ls)
    ax.text(0.86, 0.66, "f·φ = 0.43–0.53", fontsize=6.5, rotation=-28)
    for key, col, lbl, off in (("fit_fractures", C[0], "45 % of the\nfractures", (-7, -17)), ("fit_volume", C[1], "half the\nvolume", (-7, -17))):
        c = FE["cases"][key]; ax.scatter(c["phi"], c["f"], s=30, color=col, edgecolor="k", lw=0.6, zorder=4, clip_on=False)
        ax.annotate(lbl, (c["phi"], c["f"]), xytext=off, textcoords="offset points", fontsize=6.2, ha="right")
    ax.scatter(1.0, 1.0, s=30, marker="s", color="#666666", edgecolor="k", lw=0.6, zorder=4, clip_on=False)
    ax.annotate("as published", (1.0, 1.0), xytext=(7, -10), textcoords="offset points", fontsize=6.2, ha="left")
    ax.set_xlim(0.3, 1.2); ax.set_ylim(0.2, 1.0); ax.set_xlabel("Swept share of the published volume, φ"); ax.set_ylabel("Active share of the fractures, f")
    ax.set_title("(b) What the record constrains", loc="left"); ax.grid(False)
    fig.tight_layout(w_pad=1.0); save(fig, "fig11_project_red.png")


FIGS = {1: fig1, 2: fig2, 3: fig3, 4: fig4, 5: fig5, 6: fig6, 7: fig7, 8: fig8, 9: fig9, 10: fig10, 11: fig11}
if __name__ == "__main__":
    which = [int(a) for a in sys.argv[1:]] or list(FIGS)
    for n in which:
        FIGS[n]()
