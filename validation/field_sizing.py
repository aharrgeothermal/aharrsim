"""Wells needed for a 5 GW field, from the deep-study results (deep_study.py).

Basis
  thermal    reference design at 10 kg/s (3-1/2" x 4-1/2" VIT, k_ins 0.03, rock k 2.8 W/m/K, 40 C injection);
             enough wells that the field still delivers 5 GWth in year 30.
  electric   best of 5/10/15/20 kg/s for each depth and gradient.  ORC conversion as in GeoCLUSTER
             (clg_tea_module.py): utilisation efficiency 0.20 of the produced exergy at 100 C rising linearly to
             0.45 at 200 C, nil below 100 C (held at 0.45 above 200 C).  Exergy here is the thermal exergy of the
             delivered stream, dead state 20 C at the delivery pressure.  Circulation pump efficiency 0.8;
             net = gross - pumping; year-30 output.
  max case   "high-high": 60 C/km, rock k 3.5 W/m/K, VIT k 0.01 W/m/K, large bore (9-5/8" casing to TD,
             4-1/2" x 5-1/2" VIT), flows 5-40 kg/s.
  k(T)       reference design at 10 kg/s with rock conductivity falling with temperature
             (Vosteen & Schellschmidt 2003, flat above 500 C), 30 and 60 C/km.

Writes results/S_5GW.json and prints the tables.
"""
from __future__ import annotations

import re

import numpy as np
import CoolProp.CoolProp as CP
from scipy.special import exp1

from common import load, save

TARGET_MW = 5000.0
ETA_PUMP = 0.8
T_DEAD = 293.15
T_INJ = 313.15
YEAR = 30
GRADS = [20, 25, 30, 35, 40, 45, 50, 55, 60]
DEPTHS = [4000.0, 6000.0, 8000.0, 10000.0, 12000.0]
REF_FLOWS = [5.0, 10.0, 15.0, 20.0]
T_LIMIT = 400.0          # rock temperature at TD above which the case is flagged (brittle-ductile, completion limits)


def sig(n, s=3):
    """Round a well count to s significant figures."""
    if n is None or not np.isfinite(n):
        return None
    n = float(n)
    if n < 10 ** s:
        return int(round(n))
    e = int(np.floor(np.log10(n))) - s + 1
    return int(round(n / 10 ** e) * 10 ** e)


def eta_u(T_C):
    return float(np.interp(T_C, [100.0, 200.0], [0.20, 0.45], left=0.0, right=0.45))


def grad_of(label, default=30.0):
    m = re.search(r"(\d+(?:\.\d+)?) C/km", label)
    return float(m.group(1)) if m else default


def row(rec, y=YEAR):
    return next(r for r in rec["rows"] if r["years"] == y)


def h_out_of(rec, r):
    if "h_out" in r:
        return r["h_out"]
    h_in = CP.PropsSI("H", "P", r["p_in_MPa"] * 1e6, "T", T_INJ, "Water")
    return h_in + r["Q_MW"] * 1e6 / rec["mdot"]


def thermal_exergy_W(m, p_MPa, h):
    P = p_MPa * 1e6
    s = CP.PropsSI("S", "P", P, "H", h, "Water")
    h0 = CP.PropsSI("H", "P", P, "T", T_DEAD, "Water")
    s0 = CP.PropsSI("S", "P", P, "T", T_DEAD, "Water")
    return m * ((h - h0) - T_DEAD * (s - s0))


def electric(rec, y=YEAR):
    r = row(rec, y)
    ex = thermal_exergy_W(rec["mdot"], r["p_del_MPa"], h_out_of(rec, r))
    gross = eta_u(r["T_out"]) * ex
    pump = max(0.0, r["pump_kW"]) * 1e3 / ETA_PUMP
    return dict(gross_kWe=gross / 1e3, pump_kWe=pump / 1e3, net_kWe=(gross - pump) / 1e3, T_out=r["T_out"],
                exergy_kW=ex / 1e3, eta_u=eta_u(r["T_out"]))


def interference(k, r_w, spacing, years=30.0, rhoc=2650.0 * 900.0, n=15):
    """Extra temperature drop at a well's wall from an infinite square array of identical line sources,
    as a fraction of the well's own (formation) drop; constant flux, Carslaw & Jaeger line source."""
    a4t = 4.0 * k / rhoc * years * 365.25 * 86400.0
    own = exp1(r_w ** 2 / a4t)
    tot = 0.0
    for i in range(-n, n + 1):
        for j in range(-n, n + 1):
            if i or j:
                tot += exp1((i * i + j * j) * spacing ** 2 / a4t)
    return tot / own, np.sqrt(a4t)


def main():
    runs = {}
    sources = [("D_depth_flow", "ref"), ("D_gradient", "ref"), ("D_gradient_ext", "ref"), ("D_gradient_flow", "ref"),
               ("D_max", "max"), ("D_kT", "kT")]
    missing = []
    for name, var in sources:
        try:
            recs = load(name)
        except FileNotFoundError:
            missing.append(name)
            continue
        for rec in recs:
            if "rows" not in rec:
                print("  failed run:", rec.get("label"), rec.get("failed"))
                continue
            g = 60.0 if var == "max" else grad_of(rec["label"])
            runs[(var, g, rec["depth"], rec["mdot"])] = rec
    flags = [(k, v["flags"]) for k, v in runs.items() if v.get("flags") and
             (v["flags"]["unconverged_steps"] or v["flags"]["choke_warning_steps"])]

    out = dict(basis=__doc__, target_MW=TARGET_MW, year=YEAR, grads=GRADS, depths_km=[d / 1000 for d in DEPTHS],
               T_limit=T_LIMIT, missing=missing, flagged=[[list(k), v] for k, v in flags])

    # ---------------- thermal, reference design at 10 kg/s ----------------
    th = {}
    for g in GRADS:
        for d in DEPTHS:
            rec = runs.get(("ref", float(g), d, 10.0))
            if rec is None:
                continue
            r = row(rec)
            n = TARGET_MW / r["Q_MW"]
            th[f"{g}|{d / 1000:.0f}"] = dict(
                grad=g, depth_km=d / 1000, Q_MW=r["Q_MW"], T_out=r["T_out"], Q1_MW=row(rec, 1)["Q_MW"],
                Q10_MW=row(rec, 10)["Q_MW"], wells=n, wells_r=sig(n), drilled_km=n * d / 1000,
                pump_field_MWe=n * max(0.0, r["pump_kW"]) / ETA_PUMP / 1e3, pump_kW=r["pump_kW"],
                T_rock_TD=r["T_rock_bottom"], beyond=r["T_rock_bottom"] > T_LIMIT,
                x_delivery=r["x_delivery"], p_wh=r["p_wh_MPa"], choke=r["choke_MPa"])
    out["thermal"] = th

    # ---------------- electric, best flow ----------------
    el = {}
    for g in GRADS:
        for d in DEPTHS:
            cand = []
            for m in REF_FLOWS:
                rec = runs.get(("ref", float(g), d, m))
                if rec is not None:
                    e = electric(rec)
                    e["mdot"] = m
                    cand.append(e)
            if not cand:
                continue
            best = max(cand, key=lambda e: e["net_kWe"])
            n = TARGET_MW * 1e3 / best["net_kWe"] if best["net_kWe"] > 0 else None
            el[f"{g}|{d / 1000:.0f}"] = dict(grad=g, depth_km=d / 1000, best=best, wells=n, wells_r=sig(n) if n else None,
                                             flows_run=[e["mdot"] for e in cand],
                                             all={f"{e['mdot']:.0f}": e["net_kWe"] for e in cand},
                                             T_rock_TD=15.0 + g * d / 1000, beyond=15.0 + g * d / 1000 > T_LIMIT)
    out["electric"] = el

    # ---------------- max case ----------------
    mx = {}
    for d in DEPTHS:
        rows_ = []
        for m in (5.0, 10.0, 20.0, 30.0, 40.0):
            rec = runs.get(("max", 60.0, d, m))
            if rec is None:
                continue
            r = row(rec)
            e = electric(rec)
            rows_.append(dict(mdot=m, Q_MW=r["Q_MW"], T_out=r["T_out"], pump_kW=r["pump_kW"],
                              pump_kWe=max(0.0, r["pump_kW"]) / ETA_PUMP, net_kWe=e["net_kWe"], gross_kWe=e["gross_kWe"],
                              p_in=r["p_in_MPa"], choke=r["choke_MPa"], x_delivery=r["x_delivery"],
                              flags=rec.get("flags")))
        if not rows_:
            continue
        # thermal pick: most heat per well with circulation power no more than 5 % of the heat
        ok = [x for x in rows_ if x["pump_kWe"] <= 0.05 * x["Q_MW"] * 1e3] or rows_
        bt = max(ok, key=lambda x: x["Q_MW"])
        be = max(rows_, key=lambda x: x["net_kWe"])
        mx[f"{d / 1000:.0f}"] = dict(depth_km=d / 1000, rows=rows_, T_rock_TD=15.0 + 60.0 * d / 1000,
                                     thermal=dict(mdot=bt["mdot"], Q_MW=bt["Q_MW"], wells=TARGET_MW / bt["Q_MW"],
                                                  wells_r=sig(TARGET_MW / bt["Q_MW"]),
                                                  pump_field_MWe=TARGET_MW / bt["Q_MW"] * bt["pump_kWe"] / 1e3,
                                                  T_out=bt["T_out"]),
                                     electric=dict(mdot=be["mdot"], net_kWe=be["net_kWe"], T_out=be["T_out"],
                                                   wells=TARGET_MW * 1e3 / be["net_kWe"] if be["net_kWe"] > 0 else None,
                                                   wells_r=sig(TARGET_MW * 1e3 / be["net_kWe"]) if be["net_kWe"] > 0 else None))
    out["max"] = mx

    # ---------------- k(T) check ----------------
    kt = {}
    for g in (30.0, 60.0):
        for d in DEPTHS:
            a, b = runs.get(("ref", g, d, 10.0)), runs.get(("kT", g, d, 10.0))
            if a is None or b is None:
                continue
            qa, qb = row(a)["Q_MW"], row(b)["Q_MW"]
            kt[f"{g:.0f}|{d / 1000:.0f}"] = dict(grad=g, depth_km=d / 1000, Q_const=qa, Q_kT=qb, change=qb / qa - 1,
                                                 wells_const=sig(TARGET_MW / qa), wells_kT=sig(TARGET_MW / qb),
                                                 T_out_const=row(a)["T_out"], T_out_kT=row(b)["T_out"])
    out["kT"] = kt

    # ---------------- spacing ----------------
    sp = {}
    for lab, k, rw in (("reference (k 2.8, 8-1/2in hole)", 2.8, 0.108), ("max case (k 3.5, 12-1/4in hole)", 3.5, 0.156)):
        sp[lab] = {str(s): interference(k, rw, s)[0] for s in (75, 100, 150, 200)}
        sp[lab]["reach_m"] = interference(k, rw, 100)[1]
    out["spacing"] = sp
    save("S_5GW", out)
    return out


def report(out):
    G, D = out["grads"], out["depths_km"]
    print("\n5 GWth: wells (reference design, 10 kg/s, year 30); * = rock at TD above %.0f C" % out["T_limit"])
    print("grad   " + "".join(f"{d:>11.0f} km" for d in D))
    for g in G:
        cells = []
        for d in D:
            x = out["thermal"].get(f"{g}|{d:.0f}")
            cells.append(f"{x['wells_r']:>12,}{'*' if x['beyond'] else ' '}" if x else f"{'..':>13}")
        print(f"{g:>3} C/km" + "".join(cells))
    print("\nheat per well, year 30 (MWth) / outlet T (C)")
    for g in G:
        cells = []
        for d in D:
            x = out["thermal"].get(f"{g}|{d:.0f}")
            cells.append(f"{x['Q_MW']:>7.2f}/{x['T_out']:>5.0f}" if x else f"{'..':>13}")
        print(f"{g:>3} C/km " + " ".join(cells))
    print("\n5 GWe: wells (best flow), net kWe per well, flow")
    for g in G:
        cells = []
        for d in D:
            x = out["electric"].get(f"{g}|{d:.0f}")
            if not x:
                cells.append(f"{'..':>24}")
            elif x["wells_r"] is None:
                cells.append(f"{'nil':>24}")
            else:
                cells.append(f"{x['wells_r']:>9,}{'*' if x['beyond'] else ' '}{x['best']['net_kWe']:>7.0f}kWe@{x['best']['mdot']:>2.0f}")
        print(f"{g:>3} C/km" + "".join(cells))
    print("\nMax case (60 C/km, k 3.5, VIT 0.01, large bore), year 30")
    for d, x in out["max"].items():
        print(f"  {d:>2} km  rock {x['T_rock_TD']:.0f} C")
        for r in x["rows"]:
            print(f"     {r['mdot']:>4.0f} kg/s  Q {r['Q_MW']:6.2f} MW  T {r['T_out']:5.0f} C  pump {r['pump_kWe']:7.0f} kWe  "
                  f"net {r['net_kWe']:7.0f} kWe  p_in {r['p_in']:5.1f}  choke {r['choke']:5.1f}  x_del {r['x_delivery']:.3f}"
                  f"  flags {r['flags']}")
        t, e = x["thermal"], x["electric"]
        print(f"     -> 5 GWth: {t['wells_r']:,} wells at {t['mdot']:.0f} kg/s ({t['Q_MW']:.2f} MW, pumping {t['pump_field_MWe']:.0f} MWe);"
              f"  5 GWe: {e['wells_r']} wells at {e['mdot']:.0f} kg/s ({e['net_kWe']:.0f} kWe)")
    print("\nk(T) check (10 kg/s, year 30)")
    for k, x in out["kT"].items():
        print(f"  {x['grad']:.0f} C/km {x['depth_km']:.0f} km: {x['Q_const']:.2f} -> {x['Q_kT']:.2f} MW ({100 * x['change']:+.1f} %)"
              f"  wells {x['wells_const']:,} -> {x['wells_kT']:,}")
    print("\nspacing:", out["spacing"])
    if out["flagged"]:
        print("\nruns with convergence / choke flags:", out["flagged"])
    if out["missing"]:
        print("\nmissing:", out["missing"])


if __name__ == "__main__":
    report(main())
