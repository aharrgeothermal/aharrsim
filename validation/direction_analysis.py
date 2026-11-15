"""Tables from results/D_direction.json (direction_study.py): tubing-down vs annulus-down by geometry."""
from __future__ import annotations

import numpy as np

from common import load, save
from field_sizing import electric, row

INCH = 0.0254


def gain_loss(rec, y=10):
    """Heat picked up from the rock and heat given back to it along the annulus (10-yr profile)."""
    p = rec["profile_10y"]
    md = np.array(p["md"]); q = np.array(p["q_rock"]) * np.diff(md)      # W per cell
    gain, loss = q[q > 0].sum(), -q[q < 0].sum()
    mdc = np.array(p["md_c"])
    neg = mdc[q < 0]
    return gain / 1e6, loss / 1e6, (float(neg.max()) if neg.size else 0.0)


def main():
    R = [r for r in load("D_direction") if "rows" in r]
    out = dict(rows=[], flow=[], grad=[], checks=[])
    base = [r for r in R if r["mdot"] == 10.0 and r["grad"] == 30.0 and not r["bare"] and r["R_wall_override"] is None]
    key = lambda r: (r["casing"], r["tubing"])
    geoms = []
    for r in base:
        if key(r) not in geoms:
            geoms.append(key(r))
    print("12 km, 30 C/km, 10 kg/s, year 30.  TD = tubing down / annulus up;  AD = annulus down / tubing up (reference)")
    print(f"{'casing':8s} {'tubing':7s} {'A_ann':>6s} {'v_ann':>5s} {'v_tub':>5s} {'R_t':>5s} | {'TD: T_out':>9s} {'Q':>5s} {'pump':>6s} {'MWe':>6s} {'loss':>5s} {'xover':>5s} | "
          f"{'AD: T_out':>9s} {'Q':>5s} {'pump':>6s} {'MWe':>6s} | {'Q ratio':>7s}")
    for g in geoms:
        td = next((r for r in base if key(r) == g and r["direction"] == "tubing_down"), None)
        ad = next((r for r in base if key(r) == g and r["direction"] == "annulus_down"), None)
        if td is None or ad is None:
            continue
        A = td["A_annulus_m2"]
        v_a = 10.0 / (960.0 * A); v_t = 10.0 / (995.0 * td["A_tubing_m2"])
        rt, ra = row(td, 30), row(ad, 30)
        et, ea = electric(td), electric(ad)
        gain, loss, xo = gain_loss(td)
        rec = dict(casing=g[0], tubing=g[1], A_annulus_cm2=A * 1e4, v_annulus=v_a, v_tubing=v_t, R_t=td["R_tubing_K_m_W"],
                   Dh_cm=td["Dh_annulus_m"] * 100, clearance_in=td["clearance_in"],
                   td=dict(T_out=rt["T_out"], Q=rt["Q_MW"], pump_kW=rt["pump_kW"], MWe=et["net_kWe"] / 1e3, T_out10=row(td, 10)["T_out"],
                           Q10=row(td, 10)["Q_MW"], riser_loss=rt["riser_loss_MW"], gain10=gain, loss10=loss, crossover_md=xo,
                           p_in=rt["p_in_MPa"], T_bottom=rt["T_bottom"], flags=td.get("flags")),
                   ad=dict(T_out=ra["T_out"], Q=ra["Q_MW"], pump_kW=ra["pump_kW"], MWe=ea["net_kWe"] / 1e3, T_out10=row(ad, 10)["T_out"],
                           Q10=row(ad, 10)["Q_MW"], riser_loss=ra["riser_loss_MW"], p_in=ra["p_in_MPa"], T_bottom=ra["T_bottom"],
                           flags=ad.get("flags")))
        out["rows"].append(rec)
        print(f"{g[0]:8s} {g[1]:7s} {A * 1e4:6.0f} {v_a:5.2f} {v_t:5.2f} {td['R_tubing_K_m_W']:5.2f} | "
              f"{rt['T_out']:9.1f} {rt['Q_MW']:5.2f} {rt['pump_kW']:6.0f} {et['net_kWe'] / 1e3:6.3f} {loss:5.2f} {xo:5.0f} | "
              f"{ra['T_out']:9.1f} {ra['Q_MW']:5.2f} {ra['pump_kW']:6.0f} {ea['net_kWe'] / 1e3:6.3f} | {rt['Q_MW'] / ra['Q_MW']:7.2f}")
    print("\nFlow and gradient (year 30)")
    for r in sorted([r for r in R if (r["mdot"] != 10.0 or r["grad"] != 30.0) and not r["bare"] and r["R_wall_override"] is None],
                    key=lambda r: (r["casing"], r["tubing"], r["grad"], r["mdot"], r["direction"])):
        rr = row(r, 30); e = electric(r); gain, loss, xo = gain_loss(r)
        out["flow" if r["grad"] == 30.0 else "grad"].append(dict(label=r["label"], casing=r["casing"], tubing=r["tubing"], direction=r["direction"],
                                                                 mdot=r["mdot"], grad=r["grad"], T_out=rr["T_out"], Q=rr["Q_MW"],
                                                                 pump_kW=rr["pump_kW"], MWe=e["net_kWe"] / 1e3, loss10=loss, crossover_md=xo))
        print(f"  {r['label']:48s} T_out {rr['T_out']:6.1f}  Q {rr['Q_MW']:5.2f}  pump {rr['pump_kW']:7.0f} kW  MWe {e['net_kWe'] / 1e3:6.3f}  "
              f"loss {loss:5.2f}  xover {xo:5.0f} m")
    print("\nChecks (year 30)")
    for r in sorted([r for r in R if r["bare"] or r["R_wall_override"] is not None], key=lambda r: r["label"]):
        rr = row(r, 30); e = electric(r); gain, loss, xo = gain_loss(r)
        out["checks"].append(dict(label=r["label"], direction=r["direction"], T_out=rr["T_out"], Q=rr["Q_MW"], pump_kW=rr["pump_kW"],
                                  MWe=e["net_kWe"] / 1e3, loss10=loss, R_t=r["R_tubing_K_m_W"]))
        print(f"  {r['label']:52s} T_out {rr['T_out']:6.1f}  Q {rr['Q_MW']:5.2f}  pump {rr['pump_kW']:7.0f} kW  MWe {e['net_kWe'] / 1e3:6.3f}  loss {loss:5.2f}")
    flagged = [(r["label"], r["flags"]) for r in R if r.get("flags") and (r["flags"]["unconverged_steps"] or r["flags"]["choke_warning_steps"])]
    if flagged:
        print("\nflagged:", flagged)
    out["flagged"] = flagged
    save("S_direction", out)
    return out


if __name__ == "__main__":
    main()
