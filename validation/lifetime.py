"""Lifetime net electricity (and heat) of the development-type examples per kilometre of hole.

    python lifetime.py      -> results/LT_lifetime.json

Net electricity = gross output of the plant less the circulation or well pump, never negative, integrated over the simulated life
of each example (results/DT_types.json; no allowance for outages).  Hole = the wells and laterals of the example as modelled.
"""
from __future__ import annotations

import numpy as np

from common import load, save
from dev_types_table import PLANT_PICK, pick_power

# kilometres of hole in each example as modelled, and a note on what is counted
HOLE = {1: (3.6, "producer and injector"), 2: (2.0, "producer only"), 3: (1.5, "producer only"), 3.5: (2.5, "producer only"), 4: (None, "no wells"),
        5: (0.15, ""), 6: (0.40, ""), 7: (3.0, "existing well"), 8: (69.0, "two verticals and twelve laterals"), 9: (12.0, "two wells and two laterals"),
        10: (10.0, "producer and injector"), 11: (4.8, ""), 12: (3.5, ""), 13: (70.0, "twenty wells"), 14: (6.55, "two wells with 975 m laterals"),
        14.5: (12.34, "three wells with 1524 m laterals"), 15: (16.0, "eight existing wells"), 16: (3.0, "existing well"), 17: (None, "existing wells"),
        18: (4.5, "producer only"), 19: (7.0, "producer and injector"), 20: (20.0, "")}
HOURS = 8766.0


def integ(t, y):
    """Integral of y over time from 0 to t[-1], the first value held back to t = 0."""
    t = np.asarray(t, float); y = np.asarray(y, float)
    return float(np.trapezoid(np.concatenate([[y[0]], y]), np.concatenate([[0.0], t])))


def main():
    R = {float(r["id"]): r for r in load("DT_types")}
    S = {float(r["id"]): r for r in load("DT_summary")}
    rows = []
    for tid in sorted(R):
        r, s = R[tid], S[tid]
        km, note = HOLE[tid]
        se = r.get("series") or {}
        E_e = E_h = None; life = None
        if tid == 4:
            life = 30.0; E_e = s["net_last_MWe"] * life; E_h = None
        elif tid == 5:
            t = np.array(se["t_yr"]); on = np.array(se["on"], float); life = float(t[-1])
            E_e = 0.0; E_h = r["heat_pump"]["Q_heat_kW"] / 1e3 * integ(t, on)
        elif tid in (12, 13):
            d = np.array(se["day"], float) / 365.25; life = float(d[-1])
            E_e = integ(d, np.maximum(np.array(se["field_MWe" if tid == 13 else "MWe"], float), 0.0))
            E_h = integ(d, np.maximum(np.array(se["field_MWth" if tid == 13 else "P_daily_avg_MWth"], float), 0.0))
        else:
            t = np.array(se["t_yr"], float); life = float(t[-1])
            g = np.nan_to_num(np.array([v if v is not None else 0.0 for v in se.get("MWe", [0.0] * len(t))], float))
            wp = r.get("well_pump_kW") or {}
            if wp:
                yk = sorted(wp, key=float); pump = np.interp(t, [float(k) for k in yk], [max(float(wp[k]), 0.0) for k in yk]) / 1e3
            else:
                pump = np.zeros_like(t)
            pump = pump + float(r.get("circulation_pump_kW", 0.0)) / 1e3
            net = np.where(g > 0.0, np.maximum(g - pump, 0.0), 0.0)
            if tid == 6:
                net[:] = 0.0
            E_e = integ(t, net)
            h = np.array([v if v is not None else 0.0 for v in se.get("heat_MWth", [0.0] * len(t))], float)
            E_h = integ(t, np.maximum(h, 0.0))
            # check against the summary table
            assert abs(float(np.interp(float(s["horizon"]), t, net)) - s["net_last_MWe"]) < 0.03 + 0.02 * s["net_last_MWe"], (tid, float(np.interp(float(s["horizon"]), t, net)), s["net_last_MWe"])
        row = dict(id=tid, name=r["name"], life_yr=life, hole_km=km, hole_note=note, net_first_MWe=None, net_last_MWe=s["net_last_MWe"],
                   E_e_GWh=E_e * HOURS / 1e3, mean_net_MWe=E_e / life, E_h_GWh=None if E_h is None else E_h * HOURS / 1e3,
                   mean_heat_MWth=None if E_h is None else E_h / life)
        if km:
            row.update(GWh_e_per_km=row["E_e_GWh"] / km, kWe_per_km=row["mean_net_MWe"] * 1e3 / km,
                       GWh_th_per_km=None if E_h is None else row["E_h_GWh"] / km, kWth_per_km=None if E_h is None else row["mean_heat_MWth"] * 1e3 / km)
        rows.append(row)
    out = dict(rows=rows)
    # the fractured examples under the rock-quality scenarios (screening model), and the designs that hold their quoted capacity
    try:
        ex = load("RQ_examples")["systems"]; nd = load("RQ_needs")["systems"]
        from study_open import SYS, geometry
        rq = {}
        for sid, tid in (("10", 10.0), ("14", 14.0), ("141", 14.5), ("19", 19.0)):
            km = HOLE[tid][0]; s = SYS[int(sid)]
            o = dict(example={k: dict(GWh=ex[sid][k]["GWh"], GWh_per_km=ex[sid][k]["GWh"] / km, mean_net=ex[sid][k]["mean_net"]) for k in ("even", "good", "moderate", "poor")}, hole_km=km)
            if sid in nd:
                d = nd[sid]["design"]; n_w = {10: 2, 14: 2, 141: 3, 19: 2}[int(sid)]
                lat = d["W"] if "stretch" in s else 0.0        # horizontal wells: the lateral runs the length of the block; vertical doublets add no hole
                km_n = n_w * (s["depth"] + lat) / 1e3
                o["need"] = dict(design=d, hole_km=km_n, **{k: dict(GWh=nd[sid][k]["GWh"], GWh_per_km=nd[sid][k]["GWh"] / km_n, mean_net=nd[sid][k]["mean_net"]) for k in ("even", "good", "moderate", "poor") if k in nd[sid]})
            rq[sid] = o
        out["rock_quality"] = rq
    except Exception as e:                                           # the scenario runs may not exist yet
        out["rock_quality_error"] = repr(e)
    save("LT_lifetime", out)
    f = lambda v, d=1: "–" if v is None else f"{v:.{d}f}"
    for w in rows:
        print(f"{w['id']:>5} {w['name'][:40]:40s} life {w['life_yr']:5.1f} km {f(w['hole_km'],2):>6} net_last {w['net_last_MWe']:6.2f} mean {w['mean_net_MWe']:6.2f} GWh_e {w['E_e_GWh']:7.0f} "
              f"GWh_e/km {f(w.get('GWh_e_per_km')):>6} kWe/km {f(w.get('kWe_per_km'),0):>5} | GWh_th {f(w['E_h_GWh'],0):>6} GWh_th/km {f(w.get('GWh_th_per_km')):>6}")
    for sid, o in out.get("rock_quality", {}).items():
        print(sid, "example", {k: (round(v["GWh"]), round(v["GWh_per_km"], 1)) for k, v in o["example"].items()}, "| need", o.get("need") and (round(o["need"]["hole_km"], 1), {k: (round(v["GWh"]), round(v["GWh_per_km"], 1)) for k, v in o["need"].items() if k not in ("design", "hole_km")}))
    print(out.get("rock_quality_error", ""))


if __name__ == "__main__":
    main()
