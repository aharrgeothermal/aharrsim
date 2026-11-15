"""Screening charts: the two ways of taking heat from rock in dimensionless form, checked against the simulations.

Swept rock (water flows through equally spaced fractures; Gringarten et al., 1975) depends on two groups only:

    x = t / t_s              t_s = V (rho c)_r / (mdot cbar_w): time for the cold front to cross the block if all its heat were taken
    y = (S / 2) / sqrt(a t)  half the fracture spacing against the conduction length (a = k / (rho c)_r)

    retained fraction  theta = (T_out - T_inj) / (T_rock - T_inj) = 1 - T_D(x, y)
    T_D = L^-1{ exp[-R sqrt(s) tanh sqrt(s)] / s }(tau),   tau = 1 / y^2,   R = 1 / (x y^2)

    limits:  y << 1  the block cools as a whole: a front at x = 1 with relative spread sqrt(2 x y^2 / 3)
             y >> 1  every fracture cools alone:  T_D = erfc[1 / (2 x y)]

Touched rock (water in a sealed hole, heat by conduction only): plug-flow exchanger with the line-source conductance

    U'(t) = 4 pi k / [ln(4 a t / r_b^2) - 0.5772]   [W/m/K],     NTU = U' L / (mdot c)
    uniform rock temperature (laterals):      (T_out - T_in) / (T_rock - T_in) = 1 - exp(-NTU)
    linear gradient along the path:           T(L) = T_r(L) - g A + [T(0) - T_r(0) + g A] exp(-L / A),   A = mdot c / U'

    python screening.py          -> results/SC_screening.json
"""
from __future__ import annotations

import math
import sys

import numpy as np
import mpmath as mp
import CoolProp.CoolProp as CP
from scipy.special import erfc
from scipy.optimize import brentq

from common import save, load
import study_open as O

YEAR = 365.25 * 86400.0
GAMMA = 0.5772156649
IN = 0.0254
mp.mp.dps = 20


# ------------------------------------------------------------------------------------------------ swept rock
def TD(x, y):
    """Dimensionless outlet temperature drop of a block of equally spaced fractures (0 = rock temperature, 1 = injection temperature)."""
    if x <= 0:
        return 0.0
    if y < 0.04:                                    # thin slabs: front at x = 1 with a Gaussian spread (the inversion gets stiff here)
        return float(0.5 * erfc((1.0 - x) / (math.sqrt(2.0) * math.sqrt(2.0 * x * y * y / 3.0))))
    if x * y > 40.0 or y > 60.0:                    # isolated fractures
        return float(erfc(1.0 / (2.0 * x * y)))
    R = 1.0 / (x * y * y); tau = 1.0 / (y * y)
    F = lambda s: mp.exp(-R * mp.sqrt(s) * mp.tanh(mp.sqrt(s))) / s
    try:
        v = float(mp.invertlaplace(F, tau, method="dehoog"))
    except Exception:
        v = float("nan")
    if not np.isfinite(v):
        v = float(erfc(1.0 / (2.0 * x * y))) if y > 3 else float(0.5 * erfc((1.0 - x) / math.sqrt(4.0 * x * y * y / 3.0)))
    return min(max(v, 0.0), 1.0)


def TD_fast(x, y):
    """As TD, for arrays, with the double-precision de Hoog inversion of aharrsim.reservoir (the one the reservoir element uses)."""
    from aharrsim.reservoir import dehoog
    x = np.atleast_1d(np.asarray(x, float)); y = np.atleast_1d(np.asarray(y, float))
    R = 1.0 / (x * y * y); tau = 1.0 / (y * y)

    def F(s):
        r = np.sqrt(s)
        return np.exp(-R[:, None] * r * np.tanh(r)) / s

    with np.errstate(all="ignore"):
        v = dehoog(F, tau)
    v = np.where(np.isfinite(v), v, np.where(x > 1.0, 1.0, 0.0))        # underflow: before the front (0) or long after it (1)
    return np.clip(v, 0.0, 1.0)


def x_at(theta, y):
    """Flushed fraction x at which the retained fraction has fallen to theta, for a given y held at the end of life."""
    f = lambda lx: (1.0 - TD(10 ** lx, y)) - theta
    return 10 ** brentq(f, -3.0, 3.0, xtol=1e-4)


def swept_map(n=241):
    xs = np.logspace(-1.5, 1.5, n); ys = np.logspace(-1.5, 1.5, n)
    X, Y = np.meshgrid(xs, ys)
    th = 1.0 - TD_fast(X.ravel(), Y.ravel()).reshape(X.shape)
    return xs, ys, th


def map_check(n=300, seed=3):
    """The fast inversion against 20-digit mpmath on random points of the map."""
    rng = np.random.default_rng(seed)
    x = 10 ** rng.uniform(-1.5, 1.5, n); y = 10 ** rng.uniform(-1.3, 1.5, n)
    ref = np.array([TD(float(a), float(b)) for a, b in zip(x, y)])
    return dict(n=n, max_abs_diff=float(np.max(np.abs(TD_fast(x, y) - ref))))


def swept_points():
    """Where the examples and their scaled designs sit at 30 years, and the smallest designs that hold the quoted capacity."""
    pts = []
    for sid, s in O.SYS.items():
        for label, kw in O.DESIGNS[sid]:
            d = O.eval_design(sid, kw)
            pts.append(dict(sid=sid, label=s["label"], design=label, kind="example" if label.startswith("example") else "design", x=d["x_30"], y=d["y_30"],
                            theta_model=(d["T30"] - s["T_inj"]) / (s["T_rock"] - s["T_inj"]), theta_map=1.0 - TD(d["x_30"], d["y_30"]), V_km3=d["V_km3"], m=d["m"],
                            mean_net=d["mean_net"]))
    scr = load("ST2_screen")
    for sid, o in scr["systems"].items():
        s = O.SYS[int(sid)]
        for nd in o.get("need_fine") or o["need"]:
            d = nd["design"]
            if d:
                pts.append(dict(sid=int(sid), label=s["label"], design=f"holds {nd['target']:g} MWe", kind="need", x=d["x_30"], y=d["y_30"],
                                theta_model=(d["T30"] - s["T_inj"]) / (s["T_rock"] - s["T_inj"]), theta_map=1.0 - TD(d["x_30"], d["y_30"]), V_km3=d["V_km3"], m=d["m"],
                                mean_net=d["mean_net"]))
    return pts


def swept_check():
    """The map against the full injector -> reservoir -> producer chains of results/ST2_confirm.json (reservoir outlet temperature)."""
    rows = []
    for c in load("ST2_confirm"):
        if "failed" in c:
            continue
        s = O.SYS[c["sid"]]
        V = c["W"] * c["L"] * s["B"] * s.get("wings", 1)
        ts = V * s["rhoc"] / (c["m"] * O.cbar(s)) / YEAR
        al = s["k"] / s["rhoc"]
        for yr, T in zip(c["years"], c["T_res"]):
            x = yr / ts; y = 0.5 * (c["W"] / c["N"]) / math.sqrt(al * yr * YEAR)
            rows.append(dict(sid=c["sid"], tag=c["tag"], year=yr, x=x, y=y, theta_chain=(T - s["T_inj"]) / (s["T_rock"] - s["T_inj"]), theta_map=1.0 - TD(x, y),
                             dT_K=(T - s["T_inj"]) - (1.0 - TD(x, y)) * (s["T_rock"] - s["T_inj"])))
    d = np.array([r["dT_K"] for r in rows])
    return dict(rows=rows, n=len(rows), rms_K=float(np.sqrt(np.mean(d ** 2))), max_K=float(np.max(np.abs(d))), mean_K=float(np.mean(d)))


# ------------------------------------------------------------------------------------------------ touched rock
def U_line(k, alpha, r_b, t_s, R_b=0.0):
    """Conductance of the rock around a hole per unit length and per kelvin [W/m/K] (line source), in series with a completion resistance."""
    return 1.0 / ((math.log(4.0 * alpha * t_s / (r_b * r_b)) - GAMMA) / (4.0 * math.pi * k) + R_b)


def cbar(T1_C, T2_C, p):
    if abs(T2_C - T1_C) < 0.5:
        return CP.PropsSI("C", "T", 0.5 * (T1_C + T2_C) + 273.15, "P", p, "Water")
    return (CP.PropsSI("H", "T", T2_C + 273.15, "P", p, "Water") - CP.PropsSI("H", "T", T1_C + 273.15, "P", p, "Water")) / (T2_C - T1_C)


def chain(segs, m, T_in, t_yr, p=20e6):
    """Outlet temperature of a chain of holes. Each segment: L [m], T0, T1 (undisturbed rock temperature at its start and end along the flow),
    r_b [m], k, rhoc, n (parallel holes sharing the flow), R_b (completion resistance, m K/W). Returns (T_out, list of segment NTU)."""
    T_out = T_in + 20.0
    for _ in range(4):
        c = cbar(T_in, T_out, p)
        T = T_in; ntu = []
        for sg in segs:
            U = U_line(sg["k"], sg["k"] / sg["rhoc"], sg["r_b"], t_yr * YEAR, sg.get("R_b", 0.0))
            A = (m / sg.get("n", 1)) * c / U
            g = (sg["T1"] - sg["T0"]) / sg["L"]
            n = sg["L"] / A; ntu.append(n)
            T = sg["T1"] - g * A + (T - sg["T0"] + g * A) * math.exp(-n)
        T_out = T
    return T_out, ntu, c


def loops():
    """The sealed loops of the benchmark (types 6, 7, 8, 9, 16, 20) as exchanger chains. Coaxial wells: the annulus is the exchanger, the
    insulated inner string is taken as lossless (so the chart is an upper bound for them)."""
    R = lambda od, idd, k: math.log(od / idd) / (2.0 * math.pi * k)
    L = {}
    L["6"] = dict(name="400 m coaxial borehole (type 6)", kind="gradient", T_in=5.0, p=0.3e6, cases={"2": 2.0},
                  segs=[dict(L=400.0, T0=10.0, T1=22.0, r_b=4.25 * IN, k=2.6, rhoc=2500 * 900.0)], T_ref=22.0, T_s=10.0, years=(1, 5, 10), key=lambda r, c: r["T_out"])
    L["7"] = dict(name="3 km coaxial retrofit, 255 °C at depth (type 7)", kind="gradient", T_in=40.0, p=15e6, cases={"5.0": 5.0, "10.0": 10.0},
                  segs=[dict(L=1500.0, T0=15.0, T1=135.0, r_b=6.125 * IN, k=2.5, rhoc=2600 * 900.0), dict(L=1500.0, T0=135.0, T1=255.0, r_b=4.25 * IN, k=2.5, rhoc=2600 * 900.0)],
                  T_ref=255.0, T_s=15.0, years=(1, 10, 30), key=lambda r, c: r["variants"][c]["T_out"])
    L["16"] = dict(name="3 km abandoned gas well, 110 °C at depth (type 16)", kind="gradient", T_in=30.0, p=15e6, cases={"5": 5.0},
                   segs=[dict(L=3000.0, T0=10.0, T1=110.0, r_b=4.25 * IN, k=2.5, rhoc=2500 * 900.0)], T_ref=110.0, T_s=10.0, years=(1, 10, 30), key=lambda r, c: r["T_out"])
    L["20"] = dict(name="20 km coaxial, 500 °C at depth (type 20)", kind="gradient", T_in=40.0, p=100e6, cases={"20.0": 20.0, "40.0": 40.0},
                   segs=[dict(L=20000.0, T0=15.0, T1=500.0, r_b=0.17, k=2.8, rhoc=2700 * 950.0, R_b=R(0.34, 0.30, 1.4))], T_ref=500.0, T_s=15.0, years=(1, 10, 30),
                   key=lambda r, c: r["variants"][c]["T_out"])
    vert = lambda Lz, T_lat, k, rhoc, down: dict(L=Lz, T0=15.0 if down else T_lat, T1=T_lat if down else 15.0, r_b=0.108, k=k, rhoc=rhoc)
    L["8"] = dict(name="12 laterals of 5 km at 263 °C (type 8)", kind="uniform", T_in=40.0, p=25e6, cases={"laterals 80 m apart": 100.0, "isolated laterals": 100.0},
                  segs=[vert(4500.0, 262.5, 2.8, 2600 * 900.0, True), dict(L=5000.0, T0=262.5, T1=262.5, r_b=0.108, k=2.8, rhoc=2600 * 900.0, n=12),
                        vert(4500.0, 262.5, 2.8, 2600 * 900.0, False)], T_ref=262.5, T_s=15.0, years=(1, 10, 30), key=lambda r, c: r["variants"][c]["T_out"],
                  toe=lambda r, c: r["variants"][c]["T_toe"])
    L["9"] = dict(name="2 laterals of 3 km at 135 °C (type 9)", kind="uniform", T_in=40.0, p=15e6, cases={"35": 35.0},
                  segs=[vert(3000.0, 135.0, 2.6, 2600 * 900.0, True), dict(L=3000.0, T0=135.0, T1=135.0, r_b=0.108, k=2.6, rhoc=2600 * 900.0, n=2),
                        vert(3000.0, 135.0, 2.6, 2600 * 900.0, False)], T_ref=135.0, T_s=15.0, years=(1, 10, 30), key=lambda r, c: r["T_out"])
    return L


def touched_check():
    dt = {str(r["id"]).replace(".0", ""): r for r in load("DT_types")}
    rows = []
    for tid, lp in loops().items():
        r = dt[tid]
        for cname, m in lp["cases"].items():
            sim = lp["key"](r, cname)
            for yr in lp["years"]:
                T_sim = sim[str(yr)]["T_C"]
                T_ch, ntu, c = chain(lp["segs"], m, lp["T_in"], yr, lp["p"])
                row = dict(type=tid, name=lp["name"], case=cname, kind=lp["kind"], m=m, year=yr, T_sim=T_sim, T_chart=T_ch, dT_K=T_sim - T_ch, c=c)
                if lp["kind"] == "gradient":
                    N = sum(ntu)                                           # one gradient from surface to depth: NTU adds along the hole
                    # part of the rise that is due to the gradient: remove the (T_s - T_in) term so that all wells fall on one curve
                    e_fix = (lp["T_s"] - lp["T_in"]) * (1.0 - math.exp(-N))
                    row.update(NTU=N, eff=(T_sim - lp["T_in"] - e_fix) / (lp["T_ref"] - lp["T_s"]), eff_raw=(T_sim - lp["T_in"]) / (lp["T_ref"] - lp["T_in"]),
                               eff_curve=1.0 - (1.0 - math.exp(-N)) / N)
                else:
                    N = ntu[1]                                             # the laterals; the verticals carry the whole flow and exchange little
                    T_b, _, _ = chain(lp["segs"][:1], m, lp["T_in"], yr, lp["p"])       # temperature entering the laterals (chart value)
                    T_toe = lp["toe"](r, cname)[str(yr)] if "toe" in lp else None
                    T_end = T_toe if T_toe is not None else T_sim
                    row.update(NTU=N, eff=(T_end - T_b) / (lp["T_ref"] - T_b), eff_raw=(T_sim - lp["T_in"]) / (lp["T_ref"] - lp["T_in"]), eff_curve=1.0 - math.exp(-N),
                               T_toe=T_toe, NTU_vertical=ntu[0])
                rows.append(row)
    d = np.array([r["dT_K"] for r in rows]); rise = np.array([r["T_sim"] - loops()[r["type"]]["T_in"] for r in rows])
    return dict(rows=rows, n=len(rows), rms_K=float(np.sqrt(np.mean(d ** 2))), max_K=float(np.max(np.abs(d))), mean_K=float(np.mean(d)),
                rel_rms=float(np.sqrt(np.mean((d / rise) ** 2))), rel_max=float(np.max(np.abs(d / rise))))


# ------------------------------------------------------------------------------------------------ yardsticks
def specific_work(T_C, T_amb=15.0):
    """Gross electricity from 1 kg/s of pressurised water delivered at T_C [MJ/kg]: best of the practical ORC and flash plants."""
    p = max(2.0e6, CP.PropsSI("P", "T", T_C + 273.15, "Q", 0, "Water") + 0.5e6) if T_C < 370 else 30e6
    h = CP.PropsSI("H", "T", T_C + 273.15, "P", p, "Water")
    return O._specific_work(h, p, T_amb)


_WT = {}


def work_of_T(T_C):
    """specific_work on a 2 K table (the flash optimisation behind it is slow)."""
    if not _WT:
        Ts = np.arange(60.0, 392.1, 2.0)
        _WT["T"] = Ts; _WT["w"] = np.array([specific_work(float(T)) if T > 71 else 0.0 for T in Ts])
    return float(np.interp(T_C, _WT["T"], _WT["w"]))


def touched_yield(T_rock, t_yr, T_in=40.0, k=2.8, rhoc=2650 * 850.0, r_b=0.108):
    """Most electricity a kilometre of sealed hole in rock at T_rock can give at time t_yr (flow chosen freely), with the flow, outlet and heat."""
    U = U_line(k, k / rhoc, r_b, t_yr * YEAR)
    best = None
    for N in np.linspace(0.1, 4.0, 79):
        eff = 1.0 - math.exp(-N)
        T_out = T_in + eff * (T_rock - T_in)
        c = cbar(T_in, T_out, 20e6)
        m_km = U * 1000.0 / (c * N)
        pw = m_km * work_of_T(T_out)
        if best is None or pw > best["MWe_per_km"]:
            best = dict(NTU=float(N), T_out=T_out, kg_s_per_km=m_km, MWe_per_km=pw, MWth_per_km=m_km * c * (T_out - T_in) / 1e6)
    return best


def yield_curves():
    Ts = np.arange(100.0, 401.0, 10.0)
    return dict(T_rock=Ts.tolist(), years={str(y): [touched_yield(float(T), y)["MWe_per_km"] * 1e3 for T in Ts] for y in (0.5, 1, 5, 30)},
                kW_th_per_km_30={str(int(T)): touched_yield(float(T), 30)["MWth_per_km"] * 1e3 for T in (150.0, 200.0, 250.0, 300.0)})


def yardsticks(life_yr=30.0, T_inj_swept=60.0, T_in_loop=40.0, k=2.8, rhoc=2650 * 850.0, r_b=0.108):
    """What a megawatt needs from each contact mode, against rock temperature, for a life of life_yr years."""
    t = life_yr * YEAR; al = k / rhoc
    U = U_line(k, al, r_b, t)
    x90 = {y: x_at(0.9, y) for y in (0.2, 0.5, 1.0)}
    rows = []
    for T in (120.0, 150.0, 175.0, 200.0, 250.0, 300.0, 375.0):
        p = 30e6
        cw = O.cbar_w(T_inj_swept, T, p)
        w = work_of_T(T - 5.0)                                             # about 5 K lost on the way up
        row = dict(T_rock=T, w_MJ_per_kg=w, kg_s_per_MWe=1.0 / w if w > 0 else None, cbar_w=cw)
        for y, x in x90.items():                                           # swept rock: flow a cubic kilometre sustains for the whole life at >= 90 % of the temperature rise
            flow = x * 1e9 * rhoc / (cw * t)
            row[f"swept_kg_s_per_km3_y{y}"] = flow
            row[f"swept_km3_per_MWe_y{y}"] = 1.0 / (flow * w) if w > 0 else None
        # touched rock: laterals at uniform rock temperature, flow per kilometre chosen for the most electricity at the end of life
        best = touched_yield(T, life_yr, T_in_loop, k, rhoc, r_b)
        row.update(touched=best, touched_km_per_MWe=1.0 / best["MWe_per_km"] if best["MWe_per_km"] > 0 else None,
                   touched_kg_s_per_km_NTU1=U * 1000.0 / cbar(T_in_loop, T_in_loop + 0.632 * (T - T_in_loop), 20e6))
        rows.append(row)
    return dict(life_yr=life_yr, U_W_per_m_K=U, U_1yr=U_line(k, al, r_b, YEAR), U_10yr=U_line(k, al, r_b, 10 * YEAR), x90=x90, x95={y: x_at(0.95, y) for y in (0.2, 0.5, 1.0)},
                delta_outside_m=2.0 * math.sqrt(al * t / math.pi), rows=rows, T_inj_swept=T_inj_swept, T_in_loop=T_in_loop, k=k, r_b=r_b)


def loop_flow_runs():
    """The multilateral loop of type 8 at lower flows through the full chain: does the flow the chart picks for electricity do better than
    the example's 100 kg/s?  -> results/SC_loop_flow.json"""
    from aharrsim.systems import uloop
    from dev_types import ts_summary, plant_options
    lp = loops()["8"]
    out = []
    for m in (20.0, 30.0, 40.0, 60.0):
        ch = uloop(f"type 8 at {m:g} kg/s", 4500.0, 5000.0, 0.216, m, 40.0, 2.0e6, [dict(k=2.8, rho=2600.0, cp=900.0)], 15.0, dict(gradient_C_per_km=55.0), 30.0, n_lat=12,
                   layers=None, roughness=1e-4, dz=50.0, lat_far_field="adiabatic", lat_r_far=40.0)
        sm = ts_summary(ch["producer"][1])
        row = dict(m=m, years=[], T_sim=[], T_chart=[], MWe=[], MWth=[], p_wh_MPa=[], unconverged=int(ch["producer"][1].meta.get("unconverged_steps", 0)))
        for y in (1, 10, 30):
            po = plant_options(m, sm[y]["h"], sm[y]["p_MPa"] * 1e6, T_return_C=40.0)
            row["years"].append(y); row["T_sim"].append(sm[y]["T_C"]); row["T_chart"].append(chain(lp["segs"], m, 40.0, y, lp["p"])[0])
            row["MWe"].append(po.get("binary_MWe", 0.0)); row["MWth"].append(po.get("heat_MWth")); row["p_wh_MPa"].append(sm[y]["p_MPa"])
        out.append(row); print(row, flush=True)
        save("SC_loop_flow", out)


def main():
    xs, ys, th = swept_map()
    out = dict(swept=dict(x=xs.tolist(), y=ys.tolist(), theta=np.round(th, 5).tolist(), points=swept_points(), check=swept_check(), inversion_check=map_check(),
                          x_at={str(lv): {str(y): x_at(lv, y) for y in (0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0)} for lv in (0.99, 0.95, 0.9, 0.8, 0.5)}),
               touched=dict(check=touched_check(), yield_curves=yield_curves()), yard=yardsticks())
    save("SC_screening", out)
    c = out["swept"]["check"]; print("swept map vs full chains: n", c["n"], "rms K", round(c["rms_K"], 2), "max K", round(c["max_K"], 2), "mean K", round(c["mean_K"], 2))
    for p in out["swept"]["points"]:
        print(f"  {p['label']:>3} {p['design']:<36} x30 {p['x']:6.2f} y30 {p['y']:5.2f}  theta map {p['theta_map']:.3f} model {p['theta_model']:.3f}")
    print("x at theta:", {lv: {y: round(v, 2) for y, v in d.items()} for lv, d in out["swept"]["x_at"].items()})
    c = out["touched"]["check"]; print("touched chart vs simulated loops: n", c["n"], "rms K", round(c["rms_K"], 2), "max K", round(c["max_K"], 2), "mean K", round(c["mean_K"], 2),
                                       "rel rms", round(c["rel_rms"], 3), "rel max", round(c["rel_max"], 3))
    for r in c["rows"]:
        print(f"  type {r['type']:>2} {r['case']:<20} yr {r['year']:>2}  NTU {r['NTU']:.3f}  T sim {r['T_sim']:6.1f} chart {r['T_chart']:6.1f}  eff {r['eff']:.3f} curve {r['eff_curve']:.3f}")
    y = out["yard"]; print("U' 1/10/30 yr", round(y["U_1yr"], 2), round(y["U_10yr"], 2), round(y["U_W_per_m_K"], 2), "x90", y["x90"], "delta", round(y["delta_outside_m"], 1))
    for r in y["rows"]:
        print(f"  T {r['T_rock']:.0f}  w {r['w_MJ_per_kg']:.4f} MJ/kg  swept {r['swept_kg_s_per_km3_y0.5']:.0f} kg/s per km3 -> {r['swept_km3_per_MWe_y0.5'] or float('nan'):.3f} km3/MWe"
              f" | touched NTU {r['touched']['NTU']:.2f} T_out {r['touched']['T_out']:.0f} {r['touched']['kg_s_per_km']:.2f} kg/s/km {1e3 * r['touched']['MWe_per_km']:.1f} kWe/km"
              f" -> {r['touched_km_per_MWe'] or float('nan'):.0f} km/MWe")


def refresh_points():
    """Recompute the map, the points and the checks, keeping the touched-rock part and the yardsticks (they take minutes)."""
    out = load("SC_screening")
    xs, ys, th = swept_map()
    out["swept"].update(x=xs.tolist(), y=ys.tolist(), theta=np.round(th, 5).tolist(), inversion_check=map_check())
    out["swept"]["points"] = swept_points(); out["swept"]["check"] = swept_check()
    print("fast inversion against mpmath:", out["swept"]["inversion_check"])
    save("SC_screening", out)
    for p in out["swept"]["points"]:
        print(f"  {p['label']:>3} {p['design']:<36} x30 {p['x']:6.2f} y30 {p['y']:5.2f}  theta map {p['theta_map']:.3f} model {p['theta_model']:.3f}")


if __name__ == "__main__":
    loop_flow_runs() if "loopflow" in sys.argv else refresh_points() if "points" in sys.argv else main()
