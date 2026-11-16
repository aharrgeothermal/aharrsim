"""Figures 1-13 of the Geothermics manuscript, drawn from validation/results/*.json.

    python make_figs.py                all figures
    python make_figs.py 5 9 11         selected figures
    AUDIT=1 python make_figs.py        also check the lettering of each figure (see audit); AUDIT_JSON=file keeps the findings

The figures are drawn at the journal's widths (188 mm for the full width, 140 mm for one and a half columns) and saved at exactly that
size, with no lettering below 7 points, so that the printed size of the lettering is the size set here.
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
OUT = os.environ.get("FIG_OUT") or os.path.join(HERE, "figs")
os.makedirs(OUT, exist_ok=True)
sys.path.insert(0, ROOT + "/validation"); sys.path.insert(0, ROOT)
FS = 7.0             # smallest lettering, in points at the size the figure is drawn; Elsevier's rule of thumb for the printed size
INK, MUTED = "#222222", "#555555"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8.5, "axes.titlesize": 9, "axes.labelsize": 8.5, "legend.fontsize": FS,
                     "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "axes.linewidth": 0.7, "lines.linewidth": 1.3, "figure.dpi": 110,
                     "axes.grid": True, "grid.color": "#dddddd", "grid.linewidth": 0.5, "legend.frameon": False, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.axisbelow": True, "legend.handlelength": 1.8, "legend.handletextpad": 0.5,
                     "legend.columnspacing": 1.2, "legend.labelspacing": 0.35})
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


def audit(fig, name):
    """Check the lettering of a figure as it is drawn (set AUDIT=1).  Reports: text on text; text reaching into another panel;
    text or a legend lying on a plotted line, marker or axis (pixel test: the pixels a text paints against the pixels the lines and
    markers paint, with filled areas and grid lines left out); lettering below FS points.  Returns the list of findings."""
    import contextlib
    import matplotlib.text as mtext
    from matplotlib.collections import PolyCollection, QuadMesh
    from matplotlib.contour import ContourSet
    from matplotlib.legend import Legend
    from matplotlib.patches import Polygon, Rectangle

    dpi0 = fig.get_dpi(); fig.set_dpi(200)
    drawn = []
    orig = mtext.Text.draw

    def rec(self, renderer, *a, **k):
        if self.get_visible() and self.get_text().strip():
            drawn.append(self)
        return orig(self, renderer, *a, **k)

    def grab():
        fig.canvas.draw()
        return np.asarray(fig.canvas.buffer_rgba())[..., :3].astype(np.int16).copy()

    @contextlib.contextmanager
    def hidden(arts):
        st = [(a_, a_.get_visible()) for a_ in arts]
        for a_ in arts:
            a_.set_visible(False)
        try:
            yield
        finally:
            for a_, v in st:
                a_.set_visible(v)

    mtext.Text.draw = rec
    try:
        full = grab()
    finally:
        mtext.Text.draw = orig
    r = fig.canvas.get_renderer(); H = full.shape[0]
    texts, seen = [], set()
    for t in drawn:
        if id(t) not in seen:
            seen.add(id(t)); texts.append(t)
    box = {id(t): t.get_window_extent(r) for t in texts}
    texts = [t for t in texts if box[id(t)].width > 0 and box[id(t)].height > 0]
    legends = fig.findobj(Legend)
    fills = []
    for ax in fig.axes:
        fills += [c for c in ax.collections if isinstance(c, (PolyCollection, QuadMesh)) or (isinstance(c, ContourSet) and c.filled)]
        fills += [p for p in ax.patches if isinstance(p, (Rectangle, Polygon))]
        fills += list(ax.images)
    with hidden(texts):
        notext = grab()
    with hidden(texts + legends + fills):
        data = grab()
    tmask = np.abs(full - notext).sum(axis=2) > 24          # pixels painted by lettering (and by the boxes and leaders drawn with it)
    ink = data.min(axis=2) < 205                            # lines, markers, axes; the pale grid does not count
    ov = tmask & ink
    from scipy.ndimage import binary_dilation
    near = binary_dilation(tmask, iterations=3) & ink       # lettering within about one point (3 px at 200 dpi) of a line or marker
    short = lambda t: t.get_text().replace("\n", " / ")[:44]
    rows = lambda bb: (slice(max(int(H - bb.y1), 0), max(int(H - bb.y0) + 1, 0)), slice(max(int(bb.x0), 0), max(int(bb.x1) + 1, 0)))
    out, inline, cache = [], [], {}

    def painted(t):                                         # the pixels one text paints, with its background box if it has one
        if id(t) not in cache:
            with hidden([u for u in texts if u is not t]):
                m = np.abs(grab() - notext).sum(axis=2) > 24
            if t.get_bbox_patch() is not None:
                m[rows(t.get_bbox_patch().get_window_extent(r))] = True
            cache[id(t)] = m
        return cache[id(t)]

    for i, ti in enumerate(texts):
        bi = box[id(ti)]
        for tj in texts[i + 1:]:
            bj = box[id(tj)]
            w = min(bi.x1, bj.x1) - max(bi.x0, bj.x0); h = min(bi.y1, bj.y1) - max(bi.y0, bj.y0)
            if w > 1.0 and h > 1.0 and (not (ti.get_rotation() % 90 or tj.get_rotation() % 90) or int((painted(ti) & painted(tj)).sum()) >= 3):
                out.append(f"text on text: '{short(ti)}' x '{short(tj)}'")      # for a slanted label the boxes are not enough: compare the pixels
        for ax in fig.axes:
            if ti.axes is ax or ti.axes is None or not ax.get_visible():
                continue
            ab = ax.get_window_extent(r)
            w = min(bi.x1, ab.x1) - max(bi.x0, ab.x0); h = min(bi.y1, ab.y1) - max(bi.y0, ab.y0)
            if w > 1.0 and h > 1.0:
                out.append(f"text reaches into another panel: '{short(ti)}'")
        n = int(ov[rows(bi)].sum())
        if n >= 4 and ti.get_gid() == "inline":             # a name set on its own curve, as on a contour: looked at by eye, not a fault
            inline.append(short(ti))
        elif n >= 4:
            kind = "boxed label covers a line or marker" if ti.get_bbox_patch() is not None else "text on a line, marker or axis"
            out.append(f"{kind}: '{short(ti)}' ({n} px)")
        elif ti.get_gid() != "inline":
            g = 4                                           # the same test in a box grown by the margin
            m = int(near[max(int(H - bi.y1) - g, 0):max(int(H - bi.y0) + 1 + g, 0), max(int(bi.x0) - g, 0):max(int(bi.x1) + 1 + g, 0)].sum())
            if m >= 4:
                out.append(f"text within a point of a line, marker or axis: '{short(ti)}' ({m} px)")
        if ti.get_fontsize() < FS - 1e-6:
            out.append(f"small lettering, {ti.get_fontsize():.1f} pt: '{short(ti)}'")
        if bi.x0 < 0 or bi.y0 < 0 or bi.x1 > full.shape[1] or bi.y1 > H:
            out.append(f"text cut off at the edge of the figure: '{short(ti)}'")
    for lg in legends:
        if not lg.get_visible():
            continue
        bb = lg.get_window_extent(r)
        if bb.x0 < 0 or bb.y0 < 0 or bb.x1 > full.shape[1] or bb.y1 > H:
            out.append(f"legend cut off at the edge of the figure: legend starting '{lg.get_texts()[0].get_text()[:30] if lg.get_texts() else ''}'")
        sl = (slice(max(int(H - bb.y1) + 2, 0), max(int(H - bb.y0) - 1, 0)), slice(max(int(bb.x0) + 2, 0), max(int(bb.x1) - 1, 0)))
        n = int(ink[sl].sum())
        if n >= 6:
            first = lg.get_texts()[0].get_text()[:30] if lg.get_texts() else ""
            out.append(f"legend lies on a line, marker or axis: legend starting '{first}' ({n} px)")
    fig.set_dpi(dpi0)
    print(f"AUDIT {name}: {len(texts)} texts, {len(legends)} legends, {len(out)} findings" + (f"; names set on their own curve: {', '.join(inline)}" if inline else ""))
    for l in out:
        print("  " + l)
    AUDIT_LOG.append((name, out, inline))
    return out


AUDIT_LOG = []


def save(fig, name):
    fig.canvas.draw(); fig.set_layout_engine("none")        # settle the layout once, so that the check and the file see the same figure
    if os.environ.get("AUDIT"):
        audit(fig, name)
    fig.savefig(os.path.join(OUT, name), dpi=300)            # the file is exactly the size the figure is drawn at
    w, h = fig.get_size_inches()
    plt.close(fig); print(f"wrote {name}  ({w * 25.4:.0f} x {h * 25.4:.0f} mm)")


def panels_with_keys(nrows, ncols, figsize, key_after, key_height=0.24, **kw):
    """A grid of panels with a strip of blank cells under row key_after, one per column, to hold the legends of that row: a legend there
    cannot touch the data and does not push the columns apart.  Returns the figure, the panels (row by row) and the strip cells."""
    rows = [[f"p{r}{c}" for c in range(ncols)] for r in range(nrows)]
    rows.insert(key_after + 1, [f"k{c}" for c in range(ncols)])
    hr = [1.0] * nrows; hr.insert(key_after + 1, key_height)
    fig, d = plt.subplot_mosaic(rows, figsize=figsize, layout="constrained", height_ratios=hr, **kw)
    fig.get_layout_engine().set(h_pad=0.05, w_pad=0.06, hspace=0.0)
    keys = [d[f"k{c}"] for c in range(ncols)]
    for k in keys:
        k.axis("off")
    return fig, [d[f"p{r}{c}"] for r in range(nrows) for c in range(ncols)], keys


def key_legend(fig, cell, handles, x, **kw):
    """A legend in the blank strip cell, with its left edge at x (fraction of the figure width).  It is left out of the layout, so its
    width cannot squeeze the panels; call it once the panels are complete."""
    fig.canvas.draw()
    kw.setdefault("fontsize", FS)
    lg = fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(x, cell.get_position().y1), bbox_transform=fig.transFigure, borderaxespad=0.0, **kw)
    lg.set_in_layout(False)
    return lg


def along(ax, text, p, q, side=1, gap=6.5, **kw):
    """A label written along the straight line from p to q (data coordinates), centred between them and set off the line by gap points:
    side=1 puts it to the left of the direction p -> q, side=-1 to the right.  The slope is taken from the figure as laid out."""
    ax.figure.canvas.draw()
    (x0, y0), (x1, y1) = ax.transData.transform([p, q])
    ang = math.atan2(y1 - y0, x1 - x0); mid = ax.transData.inverted().transform(((x0 + x1) / 2, (y0 + y1) / 2))
    kw.setdefault("fontsize", FS); kw.setdefault("color", INK)
    return ax.annotate(text, mid, xytext=(-side * gap * math.sin(ang), side * gap * math.cos(ang)), textcoords="offset points", rotation=math.degrees(ang),
                       rotation_mode="anchor", ha="center", va="center", **kw)


# ------------------------------------------------------------------ Fig. 1
def fig1():
    R = {float(r["id"]): r for r in L("DT_types")}
    fig, axs, key = panels_with_keys(3, 2, figsize=(7.4, 8.1), key_after=1, key_height=0.31)

    def panel_log(ax, ids, title, flows=None, dashed_res=False):
        hs = []
        for i, t in enumerate(ids):
            sr = R[t]["series"]
            hs += ax.plot(sr["t_yr"], sr["T_C"], color=C[i % len(C)], label=f"{lab(t)} {SHORT[t]}" + (f", {flows[t]} kg/s" if flows and t in flows else ""))
            if dashed_res and sr.get("T_res_C"):
                ax.plot(sr["t_yr"], sr["T_res_C"], color=C[i % len(C)], ls="--", lw=0.9)
        ax.set_xscale("log"); ax.set_xlim(0.005, 30); ax.set_xlabel("Time (years)"); ax.set_ylabel("Produced water (°C)"); ax.set_title(title, loc="left")
        return hs

    panel_log(axs[0], [1, 2, 3, 3.5, 15], "(a) Conventional and co-produced water")
    axs[0].legend(loc="center left", bbox_to_anchor=(0.0, 0.53))
    ax = axs[1]
    s5 = R[5]["series"]; t5 = np.array(s5["t_yr"]); T5 = np.array(s5["T_C"]); on = np.array(s5["on"]) == 1
    ax.plot(np.where(on, t5, np.nan), np.where(on, T5, np.nan), color=C[0], label="5 GSHP U-tube, heating seasons")
    s6 = R[6]["series"]; ax.plot(s6["t_yr"], s6["T_C"], color=C[1], label="6 400 m coaxial BHE")
    ax.set_xlim(0, 10); ax.set_xlabel("Time (years)"); ax.set_ylabel("Loop outlet (°C)"); ax.set_title("(b) Shallow exchangers", loc="left"); ax.legend(loc="upper right")
    hc = panel_log(axs[2], [7, 8, 9, 16, 20], "(c) Deep closed loops", flows={7: 10, 8: 100, 20: 20})
    hd = panel_log(axs[3], [10, 11, 14, 14.5, 18, 19], "(d) Fractured and supercritical", dashed_res=True)
    hd.append(Line2D([], [], color=MUTED, ls="--", lw=0.9, label="reservoir outlet (dashed)"))
    ax = axs[4]; sr = R[12]["series"]
    ax.plot(sr["day"], sr["T_wellhead_C"], color=C[0], label="wellhead, daily mean")
    ax.plot(sr["day"], sr["T_fracture_C"], color=C[1], label="leaving the fracture, daily mean")
    ax.plot(sr["day"], sr["T_first_C"], color=C[2], ls="--", lw=0.9, label="first of the day")
    ax.plot(sr["day"], sr["T_last_C"], color=C[3], ls="--", lw=0.9, label="last of the day")
    ax.set_xlim(0, 365); ax.set_xlabel("Day of the first year"); ax.set_ylabel("Produced water (°C)"); ax.set_title("(e) Cyclic single fracture (12)", loc="left"); ax.legend(loc="upper right")
    ax = axs[5]; sr = R[17]["series"]
    ax.plot(sr["t_yr"], sr["T_C"], color=C[0], label="17 chamber outlet")
    ax.axhline(50, color="#888", ls=":", lw=0.9); ax.text(14.8, 52.5, "50 °C return", ha="right", fontsize=FS, color=MUTED)
    ax.set_xlim(0, 15); ax.set_xlabel("Time (years)"); ax.set_ylabel("Produced water (°C)"); ax.set_title("(f) SAGD heat store (17)", loc="left"); ax.legend(loc="upper right")
    key_legend(fig, key[0], hc, 0.012, ncol=2, columnspacing=1.0)
    key_legend(fig, key[1], hd, 0.525, ncol=2, columnspacing=1.0)
    save(fig, "fig02_temperature.png")


# ------------------------------------------------------------------ Fig. 2
def fig2():
    S = {float(r["id"]): r for r in L("DT_summary")}
    fig, (ax, st) = plt.subplots(2, 1, figsize=(5.5, 4.7), sharex=True, layout="constrained", gridspec_kw=dict(height_ratios=[5.0, 1.45]))
    fig.get_layout_engine().set(hspace=0.0, h_pad=0.035)
    for lo, hi, col, name in ((0, 80, "#f2f2f2", "heat only"), (80, 120, "#e6eef7", "binary,\n4–7 % of\nthe heat"), (120, 200, "#d4e1f1", "binary"),
                              (200, 340, "#f7e3dc", "flash")):
        for a in (ax, st):
            a.axvspan(lo, hi, color=col, lw=0, zorder=0)
        ax.text((lo + hi) / 2, 62, name, ha="center", va="top", fontsize=FS, color="#444", linespacing=1.15)
    # cases with net electricity at the horizon, on the log axis; the label goes on the side named here (default right)
    left = {3, 15}
    for t, r in S.items():
        if r["T_last_C"] is None or not (r["net_last_MWe"] and r["net_last_MWe"] > 0):
            continue
        steam = (r.get("x_last") or 0.0) > 0.01                    # steam or a two-phase mixture at the wellhead: the bands do not class it
        ax.scatter(r["T_last_C"], r["net_last_MWe"], s=30, color=C[1] if steam else C[0], marker="s" if steam else "o", zorder=3)
        o = -4.5 if t in left else 4.5
        ax.annotate(lab(t), (r["T_last_C"], r["net_last_MWe"]), xytext=(o, 0), textcoords="offset points", fontsize=FS, color=INK, va="center", ha="left" if o > 0 else "right")
    ax.set_yscale("log"); ax.set_ylim(0.003, 100); ax.set_ylabel("Electrical output at the horizon (MWe)")     # net of pumping, except types 12 and 13 (pumping not modelled: gross)
    ax.legend(handles=[Line2D([], [], marker="o", color=C[0], ms=5, ls="", label="liquid at the wellhead"),
                       Line2D([], [], marker="s", color=C[1], ms=5, ls="", label="steam or two-phase at the wellhead")],
              loc="lower right", title="bands: conversion route of a liquid stream", title_fontsize=FS, alignment="left", frameon=True, framealpha=1.0, edgecolor="#cccccc")
    # cases with none, in a strip of their own: each at its temperature, on the row and with the label side given here, so that none touch
    rows = {9: (0, "L"), 1: (0, "R"), 5: (0, "R"), 7: (1, "L"), 14: (1, "R"), 6: (1, "R"), 11: (2, "L"), 17: (3, "L"), 14.5: (4, "R"), 16: (4, "R")}
    for t, r in S.items():
        if r["T_last_C"] is None or (r["net_last_MWe"] and r["net_last_MWe"] > 0):
            continue
        row, side = rows[t]
        st.scatter(r["T_last_C"], row, s=26, facecolor="white", edgecolor=C[0], lw=1.1, zorder=3)
        st.annotate(lab(t), (r["T_last_C"], row), xytext=(4.5 if side == "R" else -4.5, 0), textcoords="offset points", fontsize=FS, color=INK, va="center",
                    ha="left" if side == "R" else "right")
    st.set_ylim(4.6, -0.6); st.set_yticks([]); st.grid(axis="y", visible=False)
    st.text(335, 2.0, "no net electricity at the horizon", fontsize=FS, color=MUTED, va="center", ha="right")
    st.set_xlim(0, 340); st.set_xlabel("Wellhead temperature at the horizon (°C)")
    save(fig, "fig03_scatter.png")


# ------------------------------------------------------------------ Fig. 3
def fig3():
    from matplotlib.legend_handler import HandlerTuple
    from matplotlib.patches import Patch
    S = {float(r["id"]): r for r in L("DT_summary")}
    QUOTED = {1: (0.0, 2.0), 2: (2.0, 6.0), 3: (3.0, 25.0), 3.5: (3.0, 25.0), 4: (0.355, 2.0), 7: (2.0, 9.0), 8: (2.05, 2.2), 10: (3.0, 5.0), 13: (50.0, 50.0),
              14: (3.5, 3.5), 14.5: (10.0, 16.0), 15: (0.02, 0.35), 18: (30.0, 30.0), 19: (50.0, 50.0), 20: (100.0, 100.0)}
    LF = {r["m"]: r for r in L("SC_loop_flow")}
    rows = []                                           # label, quoted range, gross in year 1, gross at the horizon
    for t in S:
        if t not in QUOTED:
            continue
        g1, gl = S[t]["gross_first_MWe"] or 0, S[t]["gross_last_MWe"] or 0
        if t == 8 and 30.0 in LF:                       # the same loop at the flow that suits electricity, on a row of its own
            rows.append((f"8 {SHORT[8]}, 100 kg/s", QUOTED[t], g1, gl, t))
            rows.append((f"8 {SHORT[8]}, 30 kg/s", QUOTED[t], LF[30.0]["MWe"][0], LF[30.0]["MWe"][2], t))
        else:
            rows.append((f"{lab(t)} {SHORT[t]}", QUOTED[t], g1, gl, t))
    OTHER = (7, 20)                                     # reported for a different system: steam-zone convection (7), an open multiwell plant (20)
    Z, QG = 0.01, "#b3b3b3"                             # where zero is drawn on the log axis; colour of the quoted figures
    fig, ax = plt.subplots(figsize=(5.5, 5.0), layout="constrained")
    ax.axvspan(0.008, 0.0126, color="#ececec", lw=0, zorder=0)
    ax.text(Z, 2.5, "no electricity", rotation=90, fontsize=FS, color=MUTED, ha="center", va="center")
    for i, (lb, (lo, hi), g1, gl, t) in enumerate(rows):
        if t in OTHER:                                  # hatched outline: a range as it is, one figure as a narrow box
            a, b = (lo, hi) if hi / max(lo, 1e-9) >= 1.15 else (math.sqrt(lo * hi) / 1.16, math.sqrt(lo * hi) * 1.16)
            ax.barh(i, b - a, left=a, height=0.44, facecolor="white", edgecolor=QG, hatch="//////", lw=0.9, zorder=1)
        elif hi / max(lo, 1e-9) < 1.15:                 # one reported figure: a tick tall enough to show behind a marker
            ax.plot([math.sqrt(lo * hi)], [i], marker="|", ms=13, mew=2.6, color=QG, zorder=1)
        else:
            ax.plot([max(lo, Z), hi], [i, i], color=QG, lw=7, solid_capstyle="butt", zorder=1)
        ax.scatter(max(g1, Z), i, s=46, facecolor="none", edgecolor=C[0], lw=1.3, zorder=3)
        ax.scatter(max(gl, Z), i, s=24, color=C[1], zorder=4)
    ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[0] for r in rows], fontsize=7.5)
    ax.set_xscale("log"); ax.set_xlim(0.008, 200); ax.set_ylim(len(rows) - 0.4, -0.6)
    ax.set_xticks([0.01, 0.1, 1, 10, 100]); ax.set_xticklabels(["0", "0.1", "1", "10", "100"]); ax.set_xlabel("Gross electricity (MWe), log scale")
    fb = (math.log10(0.0126) - math.log10(0.008)) / (math.log10(200) - math.log10(0.008))      # the break between zero and the log scale
    for dx in (-0.004, 0.004):
        ax.plot([fb + dx - 0.004, fb + dx + 0.004], [-0.014, 0.014], transform=ax.transAxes, color="k", lw=0.7, clip_on=False, zorder=6)
    fig.legend(handles=[(Patch(facecolor=QG, lw=0), Line2D([], [], marker="|", ms=10, mew=2.4, color=QG, ls="")),
                        Patch(facecolor="white", edgecolor=QG, hatch="//////", lw=0.9),
                        Line2D([], [], marker="o", mfc="none", mec=C[0], mew=1.3, ms=6.5, ls=""), Line2D([], [], marker="o", color=C[1], ms=5, ls="")],
               labels=["reported (bar: range; tick: one figure)", "reported for a different system (types 7 and 20)", "model, year 1", "model, horizon"],
               handler_map={tuple: HandlerTuple(ndivide=None, pad=0.6)}, loc="outside lower center", ncol=2, fontsize=FS, columnspacing=1.6, handlelength=3.0)
    save(fig, "fig04_claims.png")


# ------------------------------------------------------------------ Fig. 4
def fig4():
    V9 = L("V9_single_well"); V9b = L("V9b_twophase_production_mode"); B4 = L("B4_geocluster_utube"); V11 = L("V11_plant"); V10 = L("V10_reservoir"); CF = L("ST2_confirm")
    fig, axs, key = panels_with_keys(3, 2, figsize=(7.4, 8.5), key_after=1, key_height=0.31)
    ax = axs[0]; rr = V9["ramey"]["rows"]
    ax.plot([r["t_days"] for r in rr], [r["T_wh_model"] for r in rr], "o-", ms=4.5, color=C[0], label="aharrsim, IAPWS water")
    ax.plot([r["t_days"] for r in rr], [r["T_wh_model_constant_fluid"] for r in rr], "s-", ms=4.5, color=C[2], label="aharrsim, constant properties")
    ax.plot([r["t_days"] for r in rr], [r["T_wh_ramey"] for r in rr], "D--", ms=4.5, color="#555", label="Ramey / Hasan–Kabir")
    ax.set_xscale("log"); ax.set_ylim(148.9, 157.8); ax.set_xlabel("Time (days)"); ax.set_ylabel("Wellhead temperature (°C)")
    ax.set_title("(a) 3 km producer, 20 kg/s, 160 °C feed", loc="left"); ax.legend(loc="lower right")
    ax = axs[1]                                              # each well named beside its curve; the legend only tells line from marker
    for w, col in ((V9b[0], C[0]), (V9b[4], C[1])):
        ax.plot(np.array(w["profile"]["p_bar"]) / 10.0, w["profile"]["md"], color=col)
        ax.scatter(np.array(w["measured"]["p_bar"]) / 10.0, w["measured"]["md"], s=14, facecolor="white", edgecolor=col, lw=1, zorder=3)
    ax.text(4.9, 1250, V9b[0]["name"], fontsize=FS, color=INK, va="center"); ax.text(5.15, 300, V9b[4]["name"], fontsize=FS, color=INK, va="center")
    ax.invert_yaxis(); ax.set_xlabel("Pressure (MPa)"); ax.set_ylabel("Depth (m)"); ax.set_title("(b) Flowing two-phase wells", loc="left")
    ax.legend(handles=[Line2D([], [], color=INK, label="model"), Line2D([], [], marker="o", ms=3.8, mfc="white", mec=INK, ls="", label="measured")], loc="upper right")
    ax = axs[2]
    for i, c in enumerate(B4):
        comp = np.array(c["comparison"])
        ax.plot(comp[:, 0], comp[:, 1], color=C[i % 7], label=f"{c['L1'] / 1000:g} + {c['L2'] / 1000:g} km, {c['m']:g} kg/s, {c['grad'] * 1000:.0f} °C/km")
        ax.scatter(comp[:, 0], comp[:, 2], s=16, facecolor="white", edgecolor=C[i % 7], lw=1, zorder=3)
    ax.set_xscale("log"); ax.set_xlabel("Time (years)"); ax.set_ylabel("Outlet temperature (°C)"); ax.set_title("(c) U-loop chains against GeoCLUSTER (markers)", loc="left")
    hc = ax.get_legend_handles_labels()[0]
    ax = axs[3]
    b = [r for r in V11 if r["kind"] == "binary"]; f35 = [r for r in V11 if r["kind"] == "flash, condenser 35 C"]; f40 = [r for r in V11 if r["kind"] == "flash, condenser 40 C"]
    h1 = ax.plot([r["T_in_C"] for r in b], [r["eta_u_cycle"] for r in b], "o-", ms=4.5, color=C[0], label="aharrsim")
    h1 += ax.plot([r["T_in_C"] for r in b], [r["eta_u_geophires_sub"] for r in b], "D--", ms=4.5, color="#555", label="GEOPHIRES, subcritical")
    h1 += ax.plot([r["T_in_C"] for r in b], [r["eta_u_geophires_super"] for r in b], "D:", ms=4.5, color="#999", label="GEOPHIRES, supercritical")
    h2 = ax.plot([r["T_in_C"] for r in f35], [r["eta_u_double"] for r in f35], "o-", ms=4.5, color=C[1], label="aharrsim, 35 °C")
    h2 += ax.plot([r["T_in_C"] for r in f40], [r["eta_u_double"] for r in f40], "o-", ms=4.5, color=C[3], label="aharrsim, 40 °C")
    h2 += ax.plot([r["T_in_C"] for r in f35], [r["eta_u_geophires_double"] for r in f35], "D--", ms=4.5, color=C[4], label="GEOPHIRES")
    ax.set_xlabel("Inlet temperature (°C)"); ax.set_ylabel("Utilisation efficiency"); ax.set_title("(d) Plant models against GEOPHIRES", loc="left")
    ax = axs[4]; sf = V10["sharp_front"]
    ax.plot(sf["x"], sf["mpmath"], "-", color="#555", lw=2.2, alpha=0.5, label="30-digit reference")
    ax.plot(sf["x"], sf["dehoog"], "o", color=C[0], ms=4, label="de Hoog et al. (used)")
    ax.plot(sf["x"], sf["stehfest"], "s--", color=C[1], ms=3.5, lw=1, label="Gaver–Stehfest, N = 14")
    ax.set_xlabel("x = t / t$_s$"); ax.set_ylabel("Outlet temperature drop, T$_D$"); ax.set_title("(e) Fracture solution at a sharp front", loc="left"); ax.legend(loc="upper left")
    ax = axs[5]
    for sid, col, nm in ((10, C[0], "10"), (14, C[1], "14"), (141, C[2], "14b"), (19, C[3], "19")):
        xs = [g for c in CF if c.get("sid") == sid and "failed" not in c for g in c["gross"]]
        ys = [g for c in CF if c.get("sid") == sid and "failed" not in c for g in c["screen_gross"]]
        ax.scatter(np.maximum(xs, 0.05), np.maximum(ys, 0.05), s=14, facecolor="white", edgecolor=col, lw=1, label=f"type {nm}", zorder=3)
    ax.plot([0.04, 100], [0.04, 100], color="#555", lw=0.8, ls="--")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(0.04, 100); ax.set_ylim(0.04, 100)
    ax.set_xlabel("Gross electricity, full chain (MWe)"); ax.set_ylabel("Gross electricity,\nscreening form (MWe)"); ax.set_title("(f) Screening form against the full chain", loc="left")
    ax.legend(loc="upper left")
    kw = dict(title_fontsize=FS, alignment="left")
    key_legend(fig, key[0], hc, 0.012, ncol=2, title="depth + lateral, flow, gradient", columnspacing=1.0, handlelength=1.5, **kw)
    key_legend(fig, key[1], h1, 0.535, ncol=1, title="binary (ORC)", **kw)
    key_legend(fig, key[1], h2, 0.775, ncol=1, title="double flash (condenser)", **kw)
    save(fig, "fig01_verification.png")


# ------------------------------------------------------------------ Fig. 5
def fig5():
    G = L("GX_compare")["rows"]
    fig, axs = plt.subplots(1, 3, figsize=(7.4, 3.2), layout="constrained")
    mk = ["o", "s", "D", "^", "v", "P", "X", "*", "<", ">"]
    sty = {r["id"]: (mk[i], C[i % 8]) for i, r in enumerate(G)}             # one marker and colour per example, the same in all three panels
    ax = axs[0]
    for r in G:
        m, col = sty[r["id"]]
        ax.scatter(r["gx_T"], r["coax_T"], s=24, marker=m, facecolor="white", edgecolor=col, lw=1.1, zorder=3)
    ax.plot([25, 400], [25, 400], color="#555", lw=0.8, ls="--")
    for i_ in ("3", "19"):                                                   # the two flashing wells, which lie off the line (see the text)
        r = next(g for g in G if g["id"] == i_)
        ax.annotate(i_, (sum(r["gx_T"]) / 2, r["coax_T"][0]), xytext=(0, -9), textcoords="offset points", fontsize=FS, color=INK, ha="center", va="center")
    ax.set_xlim(25, 400); ax.set_ylim(25, 400)
    ax.set_xlabel("GEOPHIRES (°C)"); ax.set_ylabel("aharrsim (°C)"); ax.set_title("(a) Produced temperature,\ntwo times per case", loc="left")
    ax = axs[1]
    for i_, (xl, yl, ha) in (("14", (3.9, 125.0, "right")), ("14b", (9.9, 150.0, "left"))):
        r = next(g for g in G if g["id"] == i_); m, col = sty[i_]; h = r["hist"]
        ax.plot(h["t"], h["coax_T"], color=col)
        ax.scatter(h["t"], h["gx_T"], s=14, marker=m, facecolor="white", edgecolor=col, lw=0.9, zorder=3)
        ax.text(xl, yl, i_, fontsize=FS, color=INK, ha=ha, va="center")
    ax.set_xlim(0, 20); ax.set_xlabel("Time (years)"); ax.set_ylabel("Produced water (°C)"); ax.set_title("(b) Breakthrough: lines aharrsim,\nmarkers GEOPHIRES", loc="left")
    ax = axs[2]
    LO, FLOOR = 0.015, 0.03; mid = math.sqrt(LO * FLOOR)                     # outputs below FLOOR are drawn inside the shaded bands
    ax.axhspan(LO, FLOOR, color="#ececec", lw=0, zorder=0); ax.axvspan(LO, FLOOR, color="#ececec", lw=0, zorder=0)
    for r in G:
        m, col = sty[r["id"]]
        gx, co = np.array(r["gx_out"]), np.array(r["coax_out"])
        ax.scatter(np.where(gx < FLOOR, mid, gx), np.where(co < FLOOR, mid, co), s=24, marker=m, facecolor="white", edgecolor=col, lw=1.1, zorder=3)
    ax.plot([FLOOR, 20], [FLOOR, 20], color="#555", lw=0.8, ls="--")
    ax.text(17, mid, "below 0.03 MW", fontsize=FS, color=MUTED, ha="right", va="center")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(LO, 20); ax.set_ylim(LO, 20)
    ax.set_xlabel("GEOPHIRES (MW)"); ax.set_ylabel("aharrsim (MW)"); ax.set_title("(c) Gross electricity\n(heat for case 1)", loc="left")
    order = sorted(sty, key=lambda k: float(k.replace("b", ".5")))
    hs = [Line2D([], [], ls="", marker="", label="example")] + [Line2D([], [], ls="", marker=sty[k][0], mfc="white", mec=sty[k][1], mew=1.1, ms=5.5, label=k) for k in order]
    fig.legend(handles=hs, loc="outside lower center", ncol=len(hs), fontsize=FS, columnspacing=1.3, handletextpad=0.2, handlelength=1.2)
    save(fig, "fig06_geophires.png")


BAR_FULL, BAR_PALE, BAND = ("#1f77b4", "#d62728"), ("#9ac2dd", "#ed9e9e"), "#e6e6e6"       # below | above base; paler for a design variable


def tornado(ax, rows, base, band, xlabel, title, quoted=None, wrap=False):
    """One bar per input from the output at its low value to the output at its high value, about the base case."""
    import re
    rows = sorted(rows, key=lambda r: abs(r[2] - r[1]))
    if band is not None:
        ax.axvspan(band[0], band[2], color=BAND, lw=0, zorder=0)
    for i, (lb, lo, hi, des) in enumerate(rows):
        for v in (lo, hi):
            ax.barh(i, v - base, left=base, height=0.62, color=(BAR_PALE if des else BAR_FULL)[v >= base], lw=0, zorder=2)
    two = lambda t: re.sub(r"^(.*?)\s+([−+×±\d].*)$", r"\1\n\2", t) if wrap else t          # name on one line, the two values on the next
    ax.set_yticks(range(len(rows))); ax.set_yticklabels([two(r[0]) for r in rows], fontsize=FS, linespacing=1.05)
    if band is not None:
        ax.axvline(band[1], color="#6e6e6e", lw=1.1, zorder=3)
    ax.axvline(base, color="k", lw=0.9, ls="--", zorder=4)
    for q in quoted or ():
        ax.axvline(q, color="k", lw=1.0, ls=(0, (1, 1.6)), zorder=4)
    ax.set_xlabel(xlabel); ax.set_title(title, loc="left"); ax.grid(axis="y", visible=False)


def tornado_key(fig, design=False, quoted=False, sample="P10–P90 and median of the sample", ncol=2):
    """The key of a tornado figure, under the panels."""
    from matplotlib.legend_handler import HandlerTuple
    from matplotlib.patches import Patch
    hs = [((Patch(facecolor=BAR_FULL[0], lw=0), Patch(facecolor=BAR_FULL[1], lw=0)), "output below | above the base case")]
    if design:
        hs.append(((Patch(facecolor=BAR_PALE[0], lw=0), Patch(facecolor=BAR_PALE[1], lw=0)), "the same for a design variable"))
    hs.append((Line2D([], [], color="k", lw=0.9, ls="--"), "base case"))
    if quoted:
        hs.append((Line2D([], [], color="k", lw=1.0, ls=(0, (1, 1.6))), "published capacity (Table 7)"))
    hs.append(((Patch(facecolor=BAND, lw=0), Line2D([], [], color="#6e6e6e", lw=1.1)), sample))
    fig.legend(handles=[h for h, _ in hs], labels=[l for _, l in hs], handler_map={tuple: HandlerTuple(ndivide=None, pad=0.3)}, loc="outside lower center", ncol=ncol,
               fontsize=FS, columnspacing=2.0, handlelength=2.6)


# ------------------------------------------------------------------ Fig. 6
def fig6():
    import study_cl as CL
    by = {(r["case"], r["tag"]): r for r in L("ST1_oat") if "failed" not in r and r["case"] != 8}
    by.update({(r["case"], r["tag"]): r for r in L("ST1_oat_8_m30") if "failed" not in r})          # type 8 at the flow that suits electricity (30 kg/s)
    MC = [r for r in L("ST1_mc") if r["case"] != 8] + L("ST1_mc_8_m30")
    CL.OAT[8] = [("mdot", 20.0, 50.0, "Flow 20 / 50 kg/s", "des") if o[0] == "mdot" else o for o in CL.OAT[8]]
    fig, axs = plt.subplots(1, 3, figsize=(7.4, 3.9), layout="constrained")
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
        tornado(ax, rows, base, mc, xl, ttl, wrap=True)
    tornado_key(fig, design=True)
    save(fig, "fig07_closed_loop_uncertainty.png")


# ------------------------------------------------------------------ Fig. 7
def fig7():
    import study_open as OP
    S = L("ST2_screen")["systems"]
    G = L("RV1_gross_need")["systems"]
    fig, axs = plt.subplots(2, 4, figsize=(7.4, 5.25), layout="constrained")
    for j, sid in enumerate(("10", "14", "141", "19")):
        o = S[sid]; sy = OP.SYS[int(sid)]; q = sorted(set(sy["quoted"])); nm = sy["label"]
        ax = axs[0, j]
        for i, d in enumerate(o["density"]):
            ax.semilogx([r["N"] for r in d["rows"]], [r["mean_net"] for r in d["rows"]], "o-", ms=2.3, color=C[i], label=f"tortuosity {d['tau']:.0f}")
        ax.set_title(f"({'abcd'[j]}) Type {nm}: {sy['m']:.0f} kg/s,\n{o['base']['V_km3']:.3g} km³", loc="left", fontsize=8); ax.set_xlabel("Number of fractures"); ax.set_ylim(bottom=0)
        if j == 0:
            ax.set_ylabel("30-year mean net (MWe)"); ax.legend(loc="lower right")
        ax = axs[1, j]
        r0 = o["flow"][0]["rows"]; ax.plot([r["m"] for r in r0], [r["mean_net"] for r in r0], "k--", lw=1, label="example fractures, volume ×1")
        for i, v in enumerate([v for v in o["volume"] if v["vf"] in (1, 2, 4, 8, 16)]):
            ax.plot([r["m"] for r in v["rows"]], [r["mean_net"] for r in v["rows"]], "o-", ms=2.3, color=C[i], label=f"fractures 10 m apart, volume ×{v['vf']}")
        for qq in q:
            ax.axhline(qq, color="k", lw=1.0, ls=(0, (1, 1.6)))
        for row in G[sid]["rows"]:                          # the designs of Table 7: the figure as a net target and as a gross target
            for key, face in (("net_basis", "k"), ("gross_basis", "white")):
                if row[key]:
                    ax.plot(row[key]["m"], row[key]["mean_net"], "*", ms=9.5, mfc=face, mec="k", mew=0.9, zorder=6, clip_on=False)
        ax.set_title(f"({'efgh'[j]}) Type {nm}: fractures\n10 m apart", loc="left", fontsize=8); ax.set_xlabel("Flow (kg/s)"); ax.set_ylim(bottom=0)
        if j == 0:
            ax.set_ylabel("30-year mean net (MWe)")
    h, l = axs[1, 0].get_legend_handles_labels()
    l = [t.replace("fractures 10 m apart, ", "") for t in l]
    h.append(Line2D([], [], color="k", lw=1.0, ls=(0, (1, 1.6)))); l.append("published capacity")
    h.append(Line2D([], [], ls="", marker="*", ms=9.5, mfc="k", mec="k", mew=0.9)); l.append("net-target design of Table 7")
    h.append(Line2D([], [], ls="", marker="*", ms=9.5, mfc="white", mec="k", mew=0.9)); l.append("gross-target design of Table 7")
    fig.legend(h, l, loc="outside lower center", ncol=5, fontsize=FS, columnspacing=1.3, handlelength=1.6, title="Bottom row", title_fontsize=FS)
    save(fig, "fig08_fractured_sweeps.png")


# ------------------------------------------------------------------ Fig. 8
def fig8():
    import study_open as OP
    P3 = L("ST3_open_uncertainty")["systems"]
    fig, axs = plt.subplots(2, 2, figsize=(7.4, 5.6), layout="constrained"); axs = axs.ravel()
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
    tornado_key(fig, quoted=True, sample="Monte Carlo P10–P90 and median, 15 °C ambient")
    save(fig, "fig09_fractured_uncertainty.png")


# ------------------------------------------------------------------ Fig. 9
def fig9():
    SC = L("SC_screening")["swept"]; FE = L("FE_field_cases")["red"]
    X, Y = np.meshgrid(SC["x"], SC["y"]); TH = np.array(SC["theta"])
    fig, ax = plt.subplots(figsize=(5.5, 5.35), layout="constrained")
    levels = [0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1.0001]
    cols = ["#3b6fb6", "#6f9bd1", "#a9c5e6", "#d9e6f2", "#f6e7c8", "#f3cf9a", "#eba56a", "#dc7245"]
    cf = ax.contourf(X, Y, TH, levels=levels, colors=cols, alpha=0.85)
    cs = ax.contour(X, Y, TH, levels=[0.1, 0.5, 0.9, 0.99], colors="#333333", linewidths=0.6)
    # each label is set on its own line, where the lines run parallel and nothing else is drawn
    ax.clabel(cs, fmt={0.1: "10 %", 0.5: "50 %", 0.9: "90 %", 0.99: "99 %"}, fontsize=FS, inline_spacing=2, manual=[(0.092, 3.0), (0.0896, 4.8), (0.2, 5.25), (0.451, 12.5)])
    cb = fig.colorbar(cf, ax=ax, pad=0.02, ticks=[0, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99]); cb.set_label("Retained share of the temperature rise, θ", fontsize=7.5); cb.ax.tick_params(labelsize=FS)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(0.05, 20); ax.set_ylim(0.05, 20)
    ax.set_xlabel("x = t / t$_s$   (dimensionless thermal throughput)"); ax.set_ylabel("y = (S/2) / √(αt)   (half spacing over conduction length)")
    note = dict(fontsize=FS, color=INK, ha="right", linespacing=1.15, bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="none", alpha=0.85))
    ax.text(17.5, 7.0, "large y: each fracture cools alone;\nθ = 90 % where x·y = 0.43", va="center", **note)        # the two regimes, written where the map is empty
    ax.text(17.5, 0.073, "small y: the block cools as a whole;\nthe front arrives at x = 1", va="center", **note)
    ax.axvline(1.0, color="#333", lw=0.7, ls=":")
    ex = {p["label"]: p for p in SC["points"] if p["kind"] == "example"}
    need = {}
    for p in SC["points"]:
        if p["kind"] == "need":
            need.setdefault(p["label"], p)
    offs = {"10": (-9, -8), "14": (9, 8), "14b": (0, -10), "19": (-6, 10)}                    # label positions that keep clear of the lines
    for lb, p in ex.items():
        ax.scatter(p["x"], p["y"], s=34, color="k", zorder=5)
        ax.annotate(lb, (p["x"], p["y"]), xytext=offs[lb], textcoords="offset points", fontsize=7.5, fontweight="bold", color="k", ha="center", va="center", zorder=6)
        if lb in need:
            q = need[lb]
            ax.annotate("", xy=(q["x"], q["y"]), xytext=(p["x"], p["y"]), arrowprops=dict(arrowstyle="->", color="k", lw=0.8, shrinkA=3, shrinkB=6), zorder=4)
            ax.scatter(q["x"], q["y"], s=60, marker="*", facecolor="white", edgecolor="k", lw=0.9, zorder=5)
    # Project Red as operated (36 kg/s): path of the published block through time, its position at day 614, and the two readings of the early decline
    t = np.logspace(math.log10(0.3), math.log10(30), 40); ts = FE["t_sweep_yr"]; te = FE["t_end_yr"]
    xr = t / ts; yr = FE["y_end"] * np.sqrt(te / t)
    ax.plot(xr, yr, color=C[1], lw=1.1, ls="--", zorder=4)
    cf_, cv = FE["cases"]["fit_fractures"], FE["cases"]["fit_volume"]
    ax.plot([cf_["x_end"], FE["x_end"], cv["x_end"]], [cf_["y_end"], FE["y_end"], cv["y_end"]], color=C[1], lw=0.7, zorder=5)
    ax.scatter(FE["x_end"], FE["y_end"], s=38, marker="D", color=C[1], zorder=6)
    ax.scatter(cf_["x_end"], cf_["y_end"], s=44, marker="^", facecolor="white", edgecolor=C[1], lw=1.1, zorder=6)     # moved up: fewer fractures, wider slabs
    ax.scatter(cv["x_end"], cv["y_end"], s=44, marker=">", facecolor="white", edgecolor=C[1], lw=1.1, zorder=6)       # moved right: less volume, flushed sooner
    ax.annotate("30 yr", (xr[-1], yr[-1]), xytext=(6, 12), textcoords="offset points", fontsize=FS, color=INK, ha="center", va="center")
    hs = [Line2D([], [], marker="o", color="k", ls="", ms=5, label="example at 30 years, with its type number"),
          Line2D([], [], marker="*", mfc="white", mec="k", ls="", ms=8, label="net-target design of Table 7"),
          Line2D([], [], color=C[1], ls="--", lw=1.1, label="Project Red at 36 kg/s, 0.3 to 30 years"),
          Line2D([], [], marker="D", color=C[1], ls="", ms=4.5, label="Project Red at day 614, as published"),
          Line2D([], [], marker="^", mfc="white", mec=C[1], mew=1.1, ls="", ms=5.5, label="the same, 45 % of the fractures active"),
          Line2D([], [], marker=">", mfc="white", mec=C[1], mew=1.1, ls="", ms=5.5, label="the same, half the volume swept")]
    fig.legend(handles=hs, loc="outside lower center", ncol=2, fontsize=FS, columnspacing=1.6)
    ax.grid(False)
    save(fig, "fig11_swept_map.png")


# ------------------------------------------------------------------ Fig. 10
def fig10():
    from matplotlib.legend_handler import HandlerTuple
    SC = L("SC_screening"); T = SC["touched"]; FE = L("FE_field_cases")["geretsried"]; LF = L("SC_loop_flow")
    S = {float(r["id"]): r for r in L("DT_summary")}
    fig, axs, key = panels_with_keys(1, 2, figsize=(7.4, 4.0), key_after=0, key_height=0.2)
    ax = axs[0]
    N = np.linspace(0.01, 3.0, 200)
    f_lat = lambda n: 1 - math.exp(-n); f_cox = lambda n: 1 - (1 - math.exp(-n)) / n
    ax.plot(N, 1 - np.exp(-N), color="k", lw=1.2)
    ax.plot(N, 1 - (1 - np.exp(-N)) / N, color="k", lw=1.2, ls="--")
    mk = {"6": ("v", C[6]), "7": ("s", C[0]), "16": ("D", C[4]), "20": ("^", C[1]), "8": ("o", C[2]), "9": ("P", C[3])}
    seen = {}
    for r in T["check"]["rows"]:
        m, col = mk[r["type"]]
        h = ax.scatter(r["NTU"], r["eff"], s=20, marker=m, facecolor="white", edgecolor=col, lw=1.0, zorder=3)
        seen.setdefault(r["type"], h)
    ax.set_xlim(0, 1.4); ax.set_ylim(0, 0.8); ax.set_xlabel("NTU = U′L / (ṁc)"); ax.set_ylabel("Effectiveness"); ax.set_title("(a) Simulated loops on the exchanger curves", loc="left")
    # each curve is named along itself, above the line, where no marker lies
    along(ax, "uniform rock temperature (laterals)", (0.62, f_lat(0.62)), (1.22, f_lat(1.22)), side=1, gap=11)
    along(ax, "linear gradient, insulated return (coaxial)", (0.34, f_cox(0.34)), (1.36, f_cox(1.36)), side=1, gap=13)
    ax = axs[1]; yc = T["yield_curves"]
    for yv, ls, lbl in (("0.5", ":", "0.5 yr"), ("5", "--", "5 yr"), ("30", "-", "30 yr")):
        ax.plot(yc["T_rock"], yc["years"][yv], color="k", lw=1.1, ls=ls, label=f"chart, best flow, {lbl}")
    ax.legend(loc="upper left")
    g1, gl = S[8.0]["gross_first_MWe"], S[8.0]["gross_last_MWe"]; lf = {r["m"]: r for r in LF}; ho, fi = dict(facecolor="white"), dict(facecolor=C[2])
    c1 = ax.scatter(262.5, g1 / 60 * 1e3, s=24, marker="o", edgecolor=C[2], lw=1.1, zorder=3, **ho); c30 = ax.scatter(262.5, gl / 60 * 1e3, s=24, marker="o", edgecolor=C[2], lw=1.1, zorder=3, **fi)
    t1 = ax.scatter(262.5, lf[30.0]["MWe"][0] / 60 * 1e3, s=26, marker="^", edgecolor=C[2], lw=1.1, zorder=3, **ho)
    t30 = ax.scatter(262.5, lf[30.0]["MWe"][2] / 60 * 1e3, s=26, marker="^", edgecolor=C[2], lw=1.1, zorder=3, **fi)
    rep = FE["reported"]
    eb = ax.errorbar([160], [rep["MWe_reported"] / 30.0 * 1e3], yerr=[[rep["MWe_reported"] / 30.0 * 1e3 - 0.37 / 36.0 * 1e3], [rep["MWe_reported"] / 24.0 * 1e3 - rep["MWe_reported"] / 30.0 * 1e3]],
                     xerr=10, fmt="D", ms=4.5, color=C[1], lw=0.9, capsize=2, zorder=4)
    st = ax.scatter([160], [FE["nameplate_per_loop"]["MWe"] / 72.0 * 1e3], s=60, marker="*", facecolor="white", edgecolor=C[1], lw=1.0, zorder=4)
    ax.set_xlim(100, 400); ax.set_ylim(0, 130); ax.set_xlabel("Rock temperature at the laterals (°C)"); ax.set_ylabel("Gross electricity per km of lateral (kWe/km)")
    ax.set_title("(b) Conduction limit per kilometre of lateral", loc="left")
    order = sorted(seen, key=int)
    key_legend(fig, key[0], [seen[k] for k in order], 0.075, labels=[f"type {k}" for k in order], ncol=3, columnspacing=1.6, handletextpad=0.2)
    key_legend(fig, key[1], [(c1, c30), (t1, t30), st, eb], 0.555, ncol=1, handler_map={tuple: HandlerTuple(ndivide=None, pad=0.5)}, handlelength=2.2,
               labels=["type 8 at 100 kg/s: year 1 (open) and year 30 (filled)", "type 8 at 30 kg/s: year 1 (open) and year 30 (filled)", "Geretsried, nameplate per loop",
                       "Geretsried loop 1, as reported at 0.4 yr"])
    save(fig, "fig12_touched_chart.png")


# ------------------------------------------------------------------ Fig. 11
def fig11():
    FE = L("FE_field_cases")["red"]
    fig, axs, key = panels_with_keys(1, 2, figsize=(7.4, 3.75), key_after=0, key_height=0.11, width_ratios=[1.45, 1])
    ax = axs[0]
    t = np.array(FE["t_obs"]); T = np.array(FE["T_obs"]); used = np.array(FE["used"])
    hs = [ax.scatter(t[~used], T[~used], s=7, color="#bbbbbb", zorder=2, label="record, start-up and after shut-ins"),
          ax.scatter(t[used], T[used], s=9, color="k", zorder=3, label="record used (digitised, ±1 K)")]
    tt = np.array(FE["t_yr"])
    for k_, col, ls, lbl in (("ideal", "#666666", "--", "published volume, 88 fractures"), ("fit_fractures", C[0], "-", "45 % of the fractures active"),
                             ("fit_volume", C[1], "-", "half the volume swept")):
        hs += ax.plot(tt, FE["cases"][k_]["T_C"], color=col, ls=ls, lw=1.2, label=lbl)
    try:
        cu = L("RQ_red")["fits"]["stage"]["curve"]
        hs += ax.plot(cu["t_yr"], cu["T_C"], color=C[2], ls="-.", lw=1.2, label="unequal flow between stages, $L_c$ = 0.35")
    except (KeyError, FileNotFoundError):
        pass
    ax.axvline(FE["t_end_yr"], color="#888", lw=0.7, ls=":"); ax.text(FE["t_end_yr"] - 0.07, 112, "day 614", fontsize=FS, color=MUTED, rotation=90, ha="right", va="center")
    ax.axhline(FE["T_inj_C"], color="#888", lw=0.7, ls=":"); ax.text(7.9, FE["T_inj_C"] - 2.0, "injection temperature", fontsize=FS, color=MUTED, ha="right", va="top")
    ax.set_xlim(0, 8); ax.set_ylim(55, 185); ax.set_xlabel("Time on production (years)"); ax.set_ylabel("Flowing production temperature (°C)")
    ax.set_title("(a) Project Red: record and model at 36 kg/s", loc="left")
    ax = axs[1]; g = FE["rms_grid"]
    Fg, Pg = np.meshgrid(g["phi"], g["f"]); Rm = np.array(g["rms_K"])
    cs = ax.contourf(Fg, Pg, Rm, levels=[0, 0.5, 1.0, 2.0, 4.0, 8.0, 100.0], colors=["#dc7245", "#eba56a", "#f3cf9a", "#f6e7c8", "#d9e6f2", "#a9c5e6"], alpha=0.9)
    cb = fig.colorbar(cs, ax=ax, pad=0.02, ticks=[0, 0.5, 1, 2, 4, 8]); cb.set_label("Misfit to the record (K, rms)", fontsize=7.5); cb.ax.tick_params(labelsize=FS)
    ph = np.linspace(0.3, 1.2, 100)
    for prod in FE["acceptable"]["f_phi"][:2]:
        ax.plot(ph, prod / ph, color="k", lw=0.8, ls=":")
    lo, hi = FE["acceptable"]["f_phi"][:2]
    along(ax, f"f·φ = {lo:.2f}–{hi:.2f}", (0.86, hi / 0.86), (1.14, hi / 1.14), side=1, gap=10)     # the two dotted lines, named along the upper one
    cF, cV = FE["cases"]["fit_fractures"], FE["cases"]["fit_volume"]
    ax.scatter(cF["phi"], cF["f"], s=30, color=C[0], edgecolor="k", lw=0.6, zorder=4, clip_on=False)
    ax.annotate("45 % of the\nfractures", (cF["phi"], cF["f"]), xytext=(-5, -8), textcoords="offset points", fontsize=FS, color=INK, ha="right", va="top", linespacing=1.1)
    ax.scatter(cV["phi"], cV["f"], s=30, color=C[1], edgecolor="k", lw=0.6, zorder=4, clip_on=False)
    ax.text(0.615, 0.955, "half the\nvolume", fontsize=FS, color=INK, ha="left", va="center", linespacing=1.1)
    ax.annotate("", (cV["phi"], cV["f"]), xytext=(0.607, 0.962), textcoords="data", arrowprops=dict(arrowstyle="-", color=INK, lw=0.5, shrinkA=0, shrinkB=4))    # leader
    ax.scatter(1.0, 1.0, s=30, marker="s", color="#666666", edgecolor="k", lw=0.6, zorder=4, clip_on=False)
    ax.annotate("as published", (1.0, 1.0), xytext=(0, -8), textcoords="offset points", fontsize=FS, color=INK, ha="center", va="top")
    ax.set_xlim(0.3, 1.0); ax.set_ylim(0.2, 1.0); ax.set_xlabel("Swept share of the published volume, φ"); ax.set_ylabel("Active share of the fractures, f")
    ax.set_title("(b) Misfit to the record (K)", loc="left"); ax.grid(False)
    key_legend(fig, key[0], hs, 0.03, ncol=3, columnspacing=1.5)
    save(fig, "fig13_project_red.png")


SCEN_COL = dict(even="#666666", good="#2ca02c", moderate="#ff7f0e", poor="#d62728")
SCEN_LAB = dict(even="Uniform flow", good="Low", moderate="Moderate", poor="High")


def fig_lifetime():
    """Lifetime net electricity and heat per kilometre of hole."""
    LT = L("LT_lifetime"); rows = [r for r in LT["rows"] if r["hole_km"]]
    grp = {1: 0, 2: 0, 3: 0, 3.5: 0, 18: 0, 15: 0, 10: 1, 14: 1, 14.5: 1, 19: 1, 11: 1, 12: 3, 13: 3, 5: 2, 6: 2, 7: 2, 8: 2, 9: 2, 16: 2, 20: 2}
    gcol = ["#1f77b4", "#d62728", "#2ca02c", "#7f7f7f"]; gname = ["wells in a convective reservoir", "flow through fractured rock", "closed loops", "cyclic fracture"]
    order = [1, 2, 3, 3.5, 18, 15, 10, 14, 14.5, 19, 11, 12, 13, 8, 20, 7, 9, 16, 6, 5]
    R = {r["id"]: r for r in rows}
    fig, ax = plt.subplots(figsize=(7.4, 5.0), layout="constrained")
    ylab = []
    for i, tid in enumerate(order):
        r = R[tid]; y = len(order) - 1 - i; c = gcol[grp[tid]]
        if r["GWh_th_per_km"]:
            ax.scatter(r["GWh_th_per_km"], y, s=26, marker="s", facecolor="white", edgecolor=c, lw=1.0, zorder=3)
        if r["GWh_e_per_km"] and r["GWh_e_per_km"] > 0.05:
            ax.scatter(r["GWh_e_per_km"], y, s=34, marker="o", color=c, zorder=4)
            ax.text(r["GWh_e_per_km"] * 0.82, y, f"{r['GWh_e_per_km']:.0f}" if r["GWh_e_per_km"] >= 10 else f"{r['GWh_e_per_km']:.1f}", fontsize=FS, ha="right", va="center", color=INK)
        life = "" if abs(r["life_yr"] - 30) < 0.1 else f", {r['life_yr']:.0f} yr"
        ylab.append(f"{lab(tid)}  {SHORT[tid]} ({r['hole_km']:.3g} km{life})")
    rq = LT.get("rock_quality", {})
    for sid, tid in (("14", 14), ("141", 14.5)):                     # horizontal wells: the hole grows with the block
        o = rq.get(sid, {}).get("need")
        if not o:
            continue
        y = len(order) - 1 - order.index(tid) - 0.28
        ax.plot([o["poor"]["GWh_per_km"], o["good"]["GWh_per_km"]], [y, y], color="#d62728", lw=0.9, zorder=2)
        ax.scatter([o["poor"]["GWh_per_km"], o["moderate"]["GWh_per_km"], o["good"]["GWh_per_km"]], [y] * 3, s=10, marker="|", color="#d62728", zorder=3)
        ax.scatter(o["even"]["GWh_per_km"], y, s=26, marker="D", facecolor="white", edgecolor="#d62728", lw=1.0, zorder=4)
        ax.annotate(f"{o['hole_km']:.1f} km", (o["even"]["GWh_per_km"], y), xytext=(5.5, -2.0), textcoords="offset points", fontsize=FS, ha="left", va="center",
                    color="#d62728")                                   # the resized design has more hole than the example of its row; set just below the heat marker of the row
    ax.set_yticks(range(len(order))); ax.set_yticklabels(ylab[::-1], fontsize=FS)
    ax.set_xscale("log"); ax.set_xlim(0.3, 8000); ax.set_ylim(-0.7, len(order) - 0.4)
    ax.set_xlabel("Energy over the simulated life per kilometre of hole (GWh/km)")
    ax.grid(axis="y", visible=False)
    h = [Line2D([], [], marker="o", ls="", color="#333", ms=5.5, label="net electricity"), Line2D([], [], marker="s", ls="", mfc="white", mec="#333", ms=5, label="heat"),
         Line2D([], [], marker="D", ls="", mfc="white", mec="#d62728", ms=5, label="net electricity, net-target design of Table 7 (uniform flow)"),
         Line2D([], [], marker="|", ls="-", color="#d62728", ms=5, lw=0.9, label="the same design with high, moderate and low flow heterogeneity")]
    h += [Line2D([], [], marker="o", ls="", color=c, ms=5.5, label=n) for c, n in zip(gcol, gname)]
    fig.legend(handles=h, fontsize=FS, loc="outside lower center", ncol=2, columnspacing=2.0)
    save(fig, "fig05_lifetime.png")


def fig_rock():
    """Rock quality: flow allocation, outlet temperature under the three scenarios, effect of the scale of the unevenness."""
    from scipy.stats import norm
    from scipy.special import erfinv
    ND = L("RQ_needs")
    fig, axs = plt.subplots(1, 3, figsize=(7.4, 3.25), layout="constrained", gridspec_kw=dict(width_ratios=[1.17, 1.08, 1]))
    ax = axs[0]
    a = np.linspace(0, 1, 400)
    box = dict(boxstyle="round,pad=0.12", fc="white", ec="none")
    ax.plot(a, a, color=SCEN_COL["even"], ls="--", lw=1.0)
    for key, Lc, xl in (("good", 0.2, 0.66), ("moderate", 0.5, 0.36), ("poor", 0.8, 0.22)):      # each name sits on its own curve
        sig = 2.0 * float(erfinv(Lc)); F = norm.cdf(norm.ppf(np.clip(a, 1e-9, 1 - 1e-9)) + sig); F[0], F[-1] = 0.0, 1.0
        ax.plot(a, F, color=SCEN_COL[key], lw=1.4)
        ax.text(xl, float(np.interp(xl, a, F)), SCEN_LAB[key], fontsize=FS, color=INK, ha="center", va="center", bbox=box, gid="inline")
    pts = [(1 / 27, 0.258, "o", "k", "FORGE, Aug. 2024"), (1 / 27, 0.12, "v", "white", "FORGE, Apr. 2024\n(at most)"), (1 / 3, 2 / 3, "s", "k", "shale wells"),
           (0.12, 0.60, "^", "k", "Beowawe tracer")]
    for x, y, m, fc, t in pts:                                                                    # field measurements: named in the empty corner, not on the curves
        ax.scatter(x, y, s=22, marker=m, facecolor=fc, edgecolor="k", lw=0.8, zorder=5, label=t, clip_on=False)
    ax.legend(loc="lower right", borderaxespad=0.3, borderpad=0.3, handlelength=0.8, handletextpad=0.35, labelspacing=0.22, frameon=True, framealpha=1.0, edgecolor="#c4c4c4",
              fancybox=False).get_frame().set_linewidth(0.5)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_xlabel("Share of the fractures, fastest first"); ax.set_ylabel("Share of the flow")
    ax.set_title("(a) Flow heterogeneity", loc="left")
    along(ax, "uniform flow", (0.70, 0.70), (0.9, 0.9), side=-1, gap=10.5)
    ax = axs[1]
    yrs = np.array(ND["years"]); o = ND["systems"]["141"]
    for key in ("even", "good", "moderate", "poor"):
        T = np.array(o[key]["T_res"]); sd = np.array(o[key]["T_sd_K"])
        ax.plot(yrs, T, color=SCEN_COL[key], lw=1.4, ls="--" if key == "even" else "-")
        if key != "even":
            ax.fill_between(yrs, T - sd, T + sd, color=SCEN_COL[key], alpha=0.18, lw=0)
        ax.text(30.5, T[-1], SCEN_LAB[key], fontsize=FS, color=INK, va="center")
    ax.set_xlim(0, 30); ax.set_ylim(100, 205); ax.set_xlabel("Time (years)"); ax.set_ylabel("Reservoir outlet temperature (°C)")
    ax.set_title("(b) Net-target design of type 14b", loc="left")
    ax = axs[2]
    S2 = L("RQ_scale"); T_r, T_i = 200.0, 60.0; y2 = np.array(S2["years"])
    th = lambda T: float(np.trapezoid((np.array(T) - T_i) / (T_r - T_i), y2) / (y2[-1] - y2[0]))
    for key, Lc in (("good", 0.2), ("moderate", 0.5), ("poor", 0.8)):
        rs = [r for r in S2["rows"] if abs(r["lorenz"] - Lc) < 1e-9]
        w = [r["corridor_m"] for r in rs if r["corridor_m"]]; v = [th(r["T_res"]) for r in rs if r["corridor_m"]]
        if not w:
            continue
        ax.plot(w, v, color=SCEN_COL[key], lw=1.4, marker="o", ms=3)
        lim = [r for r in rs if not r["corridor_m"]]
        end = th(lim[0]["T_res"]) if lim else v[-1]
        if lim:
            ax.plot([560, 1000], [end] * 2, color=SCEN_COL[key], lw=1.6, ls=":")
        ax.text(1080, end, SCEN_LAB[key], fontsize=FS, color=INK, va="center")
    ax.axhline(th(S2["even"]), color=SCEN_COL["even"], lw=1.0, ls="--"); ax.text(1080, th(S2["even"]) + 0.012, "Uniform flow", fontsize=FS, color=INK, va="center")
    ax.axvline(33.7, color="#888", lw=0.7, ls=":"); ax.text(37.5, 0.03, "conduction\nlength, 30 yr", fontsize=FS, color=MUTED, ha="left", va="bottom", linespacing=1.1)
    ax.text(960, 0.46, "dotted:\nseparate\nzones", fontsize=FS, color=MUTED, ha="right", va="center", linespacing=1.1)
    ax.set_xscale("log"); ax.set_xlim(8, 1000); ax.set_ylim(0, 1.08)
    ax.set_xlabel("Width over which the flow varies (m)"); ax.set_ylabel("30-year mean of the\nretained temperature rise")
    ax.set_title("(c) Scale of the heterogeneity", loc="left")
    save(fig, "fig10_rock_quality.png")


FIGS = {1: fig4, 2: fig1, 3: fig2, 4: fig3, 5: fig_lifetime, 6: fig5, 7: fig6, 8: fig7, 9: fig8, 10: fig_rock, 11: fig9, 12: fig10, 13: fig11}
if __name__ == "__main__":
    which = [int(a) for a in sys.argv[1:]] or list(FIGS)
    for n in which:
        FIGS[n]()
    if os.environ.get("AUDIT_JSON"):
        with open(os.environ["AUDIT_JSON"], "w") as f:
            json.dump(AUDIT_LOG, f, indent=1)
