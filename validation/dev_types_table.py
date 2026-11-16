"""Summary table of the development-type runs (results/DT_types.json -> results/DT_summary.json) with the
electricity assessment per type and the capacity figures quoted in the source spreadsheet."""
from __future__ import annotations

import json
import sys

import numpy as np

from common import load, save

SHEET = {   # from geothermal_development_types.xlsx (Khodayar & Bjornsson 2024): group, stated electrical / thermal capacity, basis
    1: ("Conventional", "<2 MWe", "<5 MWth", "per well"),
    2: ("Conventional", "2–6 MWe (<7)", "<15 MWth", "per well"),
    3: ("Conventional", "3–25 MWe (up to 30–50)", "12.5 MWth typical; 70–125 MWth", "per well"),
    3.5: ("Conventional", "3–25 MWe (dry steam)", "–", "per well"),
    4: ("Surface waste heat", "355 kWe per module; 1–2 MWe stacked", "cooled brine to district heating", "per module"),
    5: ("Closed loop (shallow)", "—", "<10 kW per unit", "per unit"),
    6: ("Closed loop (shallow)", "—", "<10 kW to <5 MWth", "per installation"),
    7: ("Closed loop (deep)", "2–9 MWe claimed; 2 MWe retrofit", "heat for buildings / data centres", "per well"),
    8: ("Closed loop (deep)", "2.2 MWe (classification, attributed to Toews and Holmes, 2021); 8.2 MWe (company, four loops)", "64 MWth (company, four loops)",
        "per loop; the classification attributes the 2.2 MWe to Toews and Holmes (2021), a paper in which the figure could not be located; the company figures are the nameplate of the four-loop Geretsried project (Toews, 2026)"),
    9: ("Closed loop (deep)", "undisclosed", "undisclosed", "–"),
    10: ("EGS / HDR", "3–5 MWe", "~12 MWth", "per well"),
    11: ("EGS / HDR", "small", "small", "per well"),
    12: ("EGS / HDR", "not stated", "not stated", "per well"),
    13: ("EGS / HDR", "up to 50 MWe", "—", "per system (18–20 wells)"),
    14: ("EGS / HDR", "3.5 MWe pilot (63 l/s)", "heat (300,000 homes planned)", "per pilot"),
    14.5: ("EGS / HDR", ">10 MWe per production well (2024 test); about 16 MWe gross (current design)", "–", "per production well (Fervo, 2024, 2026; not in the classification)"),
    15: ("O&G assets", "20–350 kWe", "heating / cooling", "per field"),
    16: ("O&G assets", "—", "<1 MWth (~400 kW trial)", "per well"),
    17: ("O&G assets", "to be proven", "to be proven", "–"),
    18: ("Supercritical", ">30 MWe expected", ">75 MWth expected", "per well"),
    19: ("Supercritical", "up to 50 MWe (??)", "up to 125 MWth (??)", "per well"),
    20: ("Ultradeep", "100 MWe", "—", "per plant"),
}
PLANT_PICK = {1: "binary", 2: "binary", 3: "flash2", 3.5: "dry_steam", 4: "binary", 5: None, 6: None, 7: "binary", 8: "binary",
              9: "binary", 10: "binary", 11: "binary", 12: "binary", 13: "binary", 14: "binary", 14.5: "binary", 15: "binary", 16: "binary",
              17: "binary", 18: "dry_steam", 19: "steam", 20: "best"}


def _h70(p):
    import CoolProp.CoolProp as CP
    return CP.PropsSI("H", "T", 70.0 + 273.15, "P", max(p, 1.0e5), "Water")


def pick_power(po, kind):
    """Gross electrical MW from the plant options for the chosen plant kind."""
    if po is None or kind is None:
        return 0.0, "none"
    if kind == "binary":
        return po.get("binary_MWe", 0.0) or 0.0, "binary (GEOPHIRES practical)"
    if kind == "flash2":
        f2, f1 = po.get("flash2_MWe"), po.get("flash1_MWe")
        if f2 is not None and f1 is not None:
            return (f2, "double flash") if f2 >= f1 else (f1, "single flash")
        return po.get("binary_MWe", 0.0) or 0.0, "binary"
    if kind == "dry_steam":
        if po.get("dry_steam_MWe") is not None:
            return po["dry_steam_MWe"], "dry / superheated steam turbine"
        w, _ = pick_power(po, "flash2")
        return w, "steam turbine (wellhead separator)"
    if kind == "steam":
        if po.get("dry_steam_MWe") is not None:
            return po["dry_steam_MWe"], "superheated steam turbine"
        return pick_power(po, "flash2")
    if kind == "best":
        cands = {"binary": po.get("binary_MWe"), "single flash": po.get("flash1_MWe"), "double flash": po.get("flash2_MWe"),
                 "steam turbine": po.get("dry_steam_MWe")}
        cands = {k: v for k, v in cands.items() if v is not None}
        k = max(cands, key=cands.get)
        return cands[k], k
    return 0.0, "none"


def main():
    R = {float(r["id"]): r for r in load("DT_types")}
    rows = []
    for tid in sorted(R):
        r = R[tid]
        grp, mwe_sheet, mwth_sheet, basis = SHEET[tid if tid != 3.5 else 3.5]
        if "failed" in r:
            rows.append(dict(id=tid, name=r.get("name", ""), group=grp, failed=r["failed"]))
            continue
        T = r.get("T_out", {})
        keys = sorted(T, key=float)
        if tid == 17:   # heat store: report the last year in which it still yields heat (afterwards the circulation heats the ground)
            keys = [k for k in keys if (T[k].get("Q_MWth") or 0.0) > 0.0] or keys
        k_first, k_last = keys[0], keys[-1]
        T_first, T_last = T[k_first].get("T_C"), T[k_last].get("T_C")
        po_last = (r.get("plant") or {}).get(k_last)
        po_first = (r.get("plant") or {}).get(k_first)
        kind = PLANT_PICK.get(tid)
        gross_last, plant_name = pick_power(po_last, kind)
        gross_first, _ = pick_power(po_first, kind)
        pump = 0.0
        if r.get("well_pump_kW"):
            pump += max(0.0, float(r["well_pump_kW"].get(k_last, 0.0))) / 1e3
        pump += float(r.get("circulation_pump_kW", 0.0)) / 1e3
        heat_last = (po_last or {}).get("heat_MWth")
        heat_first = (po_first or {}).get("heat_MWth")
        if heat_last is None and po_last and (r.get("assumptions") or {}).get("mdot"):   # flash / steam wells: heat above 70 C reinjection
            m_ = float(r["assumptions"]["mdot"])
            heat_last = m_ * (T[k_last]["h"] - _h70(T[k_last]["p_MPa"] * 1e6)) / 1e6
            heat_first = m_ * (T[k_first]["h"] - _h70(T[k_first]["p_MPa"] * 1e6)) / 1e6
        ex_last = (po_last or {}).get("exergy_MW")
        extra = {}
        if tid == 5:
            hp = r["heat_pump"]; heat_last = hp["Q_heat_kW"] / 1e3; extra = dict(COP=hp["COP"], W_e_MW=hp["W_e_kW"] / 1e3, W_per_m=hp["per_metre_W"])
            gross_last = gross_first = 0.0; pump = 0.0
        if tid == 6:
            hp = r["heat_pump"]; heat_last = hp["Q_heat_kW"] / 1e3; extra = dict(COP=hp["COP"], W_e_MW=hp["W_e_kW"] / 1e3, W_per_m=hp["per_metre_W"])
            gross_last = gross_first = 0.0
        if tid == 12:
            heat_last = T[k_last]["P_daily_avg_MWth"]; heat_first = T[k_first]["P_daily_avg_MWth"]
        if tid == 13:
            cont = r["continuous"][k_last]; heat_last = cont["P_field_MWth"]; gross_last = cont["MWe_field"]
            heat_first = r["continuous"][k_first]["P_field_MWth"]; gross_first = r["continuous"][k_first]["MWe_field"]
            extra = dict(storage=r["storage"])
        if tid == 17:
            heat_last = T[k_last]["Q_MWth"]; heat_first = T[k_first]["Q_MWth"]; extra = dict(recovered_fraction_15yr=r["recovered_fraction_15yr"])
        net_last = (gross_last - pump) if gross_last > 0 else 0.0   # heat-only cases: pumps are charged to the heat supply
        rows.append(dict(id=tid, name=r["name"], group=grp, horizon=k_last, horizon_first=k_first, T_first_C=T_first, T_last_C=T_last,
                         x_last=(po_last or {}).get("x"), x_first=(po_first or {}).get("x"),
                         T_res_last_C=(r.get("T_res") or {}).get(k_last), T_res_first_C=(r.get("T_res") or {}).get(k_first), heat_first_MWth=heat_first, heat_last_MWth=heat_last, gross_first_MWe=gross_first,
                         gross_last_MWe=gross_last, pump_MWe=pump, net_last_MWe=net_last, exergy_last_MW=ex_last, plant=plant_name,
                         plant_kind=r.get("plant_kind"), sheet_MWe=mwe_sheet, sheet_MWth=mwth_sheet, sheet_basis=basis,
                         assumptions=r.get("assumptions"), flags=r.get("flags"), extra=extra, wall_s=r.get("wall_s"),
                         T_series={k: T[k].get("T_C") for k in keys},
                         binary_last=(po_last or {}).get("binary_MWe"), flash1_last=(po_last or {}).get("flash1_MWe"),
                         flash2_last=(po_last or {}).get("flash2_MWe"), steam_last=(po_last or {}).get("dry_steam_MWe"),
                         cycle_last=(po_last or {}).get("binary_cycle_MWe")))
    save("DT_summary", rows)
    print(f"{'id':>4} {'type':52s} {'T first':>8s} {'T last':>7s} {'MWth':>7s} {'gross MWe':>9s} {'pump':>6s} {'net MWe':>8s}  plant / sheet")
    for r in rows:
        if "failed" in r:
            print(f"{r['id']:>4} {r['name'][:52]:52s} FAILED {r['failed'][:60]}")
            continue
        f = lambda v, d=2: "–" if v is None else f"{v:.{d}f}"
        print(f"{r['id']:>4} {r['name'][:52]:52s} {f(r['T_first_C'],1):>8s} {f(r['T_last_C'],1):>7s} {f(r['heat_last_MWth']):>7s} {f(r['gross_last_MWe']):>9s} "
              f"{f(r['pump_MWe']):>6s} {f(r['net_last_MWe']):>8s}  {r['plant']} | sheet: {r['sheet_MWe']}  flags {r['flags']}")
    return rows


def markdown():
    """The README table, from results/DT_summary.json (and results/SC_loop_flow.json for type 8 at the flow that suits electricity)."""
    lab = lambda t: {3.5: "3b", 14.5: "14b"}.get(t, str(int(t)))
    f = lambda v, d=2: "–" if v is None else f"{v:.{d}f}"
    rng = lambda a, b, d=2: f(b, d) if (a is None or f(a, d) == f(b, d)
                                        or abs(a - (b or 0.0)) < 0.05 * max(abs(a), abs(b or 0.0), 1e-9)) else f"{f(a, d)} → {f(b, d)}"
    T = {float(r["id"]): r for r in load("DT_types")}
    print("| # | Type (example) | Horizon | Water T first → last (°C) | MWth | Gross MWe first → last | Pumps MWe first → last | Net MWe | Quoted |")
    print("|---|---|---|---|---|---|---|---|---|")
    for r in load("DT_summary"):
        if "failed" in r:
            continue
        tid = float(r["id"])
        hor = ("day 365" if tid in (12, 13) else "steady" if tid == 4 else f"year {r['horizon']}")
        wp = T[tid].get("well_pump_kW") or {}
        circ = float(T[tid].get("circulation_pump_kW", 0.0)) / 1e3
        p_first = max(0.0, float(wp.get(r["horizon_first"], 0.0))) / 1e3 + circ if wp else r["pump_MWe"]
        print(f"| {lab(tid)} | {r['name']} | {hor} | {f(r['T_first_C'], 0)} → {f(r['T_last_C'], 0)} | {f(r['heat_last_MWth'])} | "
              f"{rng(r['gross_first_MWe'], r['gross_last_MWe'])} | {rng(p_first, r['pump_MWe'])} | {f(r['net_last_MWe'])} | {r['sheet_MWe']} |")
        if tid == 8:
            try:
                q = next(x for x in load("SC_loop_flow") if abs(x["m"] - 30.0) < 1e-6)
                print(f"| 8 | the same loop at 30 kg/s, the flow that gives most electricity at year 30 | year 30 | {q['T_sim'][0]:.0f} → {q['T_sim'][-1]:.0f} | "
                      f"{q['MWth'][-1]:.2f} | {q['MWe'][0]:.2f} → {q['MWe'][-1]:.2f} | 0.00 | {q['MWe'][-1]:.2f} | as above |")
            except (FileNotFoundError, StopIteration):
                pass


if __name__ == "__main__":
    markdown() if "md" in sys.argv[1:] else main()
