"""Command-line interface.

    python -m aharrsim run CASE.yaml [--out DIR] [--progress] [--plot]

Writes to DIR (default: alongside the case file, <name>_results/):
    timeseries.csv   time history at the wellhead (SI units, temperatures in C)
    profiles.csv     depth profiles at each requested output time
    summary.json     key numbers (outlet temperature, thermal power, exergy, pumping)
    summary.png      quick-look plots (with --plot)
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

from . import load_case, Simulator, DAY, YEAR


def _summary(sim, res):
    f = sim.fluid
    c = sim.case
    T_ref = c.T_surface
    ts = res.ts
    out = dict(case=c.name, fluid=getattr(f, "name", ""), grid=dict(well_cells=sim.Nw, axial_cells=sim.N,
               radial_cells=getattr(sim, "M", 0)), wall_time_s=res.meta.get("wall_time_s"),
               points=[])
    times = sorted(set([p["t"] for p in res.profiles] + [float(ts["t"][-1])]))
    for t in times:
        i = int(np.argmin(np.abs(ts["t"] - t)))
        p_in, p_wh, h_in, h_out, m = ts["p_in"][i], ts["p_out"][i], ts["h_in"][i], ts["h_out"][i], ts["mdot"][i]
        choke = ts["choke_dp"][i]
        p_del = p_wh - choke                      # delivery pressure downstream of any wellhead choke
        dl = f.ph(np.array([p_del]), np.array([h_out]))
        # thermal exergy: work obtainable by cooling the delivered stream back to the
        # injection temperature at the delivery pressure, environment at T_surface
        h_r = f.h_pT(p_del, ts["T_in"][i])[0]
        try:
            s_r = f.s_ph(p_del, h_r)[0]
            s_out = f.s_ph(p_del, h_out)[0]
            ex = m * ((h_out - h_r) - T_ref * (s_out - s_r))
        except Exception:
            ex = float("nan")
        q_avail = m * (h_out - h_r)
        if ts["T_in"][i] < T_ref:
            ex = float("nan")    # injection colder than ambient (heat-pump source): exergy not meaningful
        rho_in = f.ph(np.array([p_in]), np.array([h_in])).rho[0]
        out["points"].append(dict(
            t_days=t / DAY, T_in_C=ts["T_in"][i] - 273.15, T_out_C=ts["T_out"][i] - 273.15,
            steam_quality_out=ts["x_out"][i], p_in_MPa=p_in / 1e6, p_wellhead_MPa=p_wh / 1e6,
            choke_dp_MPa=choke / 1e6, p_delivery_MPa=p_del / 1e6,
            T_delivery_C=float(dl.T[0] - 273.15), steam_quality_delivery=float(dl.x[0]) if dl.two_phase[0] else
            (0.0 if dl.x[0] <= 0 else 1.0),
            thermal_power_MW=ts["Q_wellhead"][i] / 1e6, heat_available_MW=q_avail / 1e6, exergy_MW=ex / 1e6,
            riser_heat_loss_MW=ts["Q_riser_loss"][i] / 1e6,
            pump_hydraulic_kW=m * (p_in - p_del) / rho_in / 1e3,
            flash_depth_m=None if not np.isfinite(ts["flash_md"][i]) else ts["flash_md"][i]))
    dt = np.diff(np.concatenate([[0.0], ts["t"]]))
    out["cumulative_heat_GWh"] = float(np.sum(dt * ts["Q_wellhead"]) / 3.6e12)
    out["unconverged_steps"] = int(res.meta.get("unconverged_steps", 0))
    out["steps_with_near_critical_riser_flow"] = int(np.sum(ts["choke_warning"] > 0))
    out["max_energy_balance_residual_kW"] = float(np.max(np.abs(ts["Q_wellhead"] - ts["Q_rock"] - ts["Q_res"] + ts["dE_fluid"] + ts["Q_ke"] + ts["Q_pe"])) / 1e3)
    return out


def _write_profiles(res, path):
    rows = []
    for p in res.profiles:
        for k in range(p["md"].size):
            rows.append([p["t"] / DAY, p["md"][k], p["tvd"][k], p["T0"][k] - 273.15, p["T_a"][k] - 273.15,
                         p["T_t"][k] - 273.15, p["p_a"][k] / 1e6, p["p_t"][k] / 1e6, p["x_t"][k]])
    hdr = "t_days,md_m,tvd_m,T_rock_undisturbed_C,T_annulus_C,T_tubing_C,p_annulus_MPa,p_tubing_MPa,steam_quality_tubing"
    np.savetxt(path, np.array(rows), delimiter=",", header=hdr, comments="", fmt="%.6g")


def _plot(sim, res, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
    t = res.ts["t"] / YEAR
    ax[0].plot(t, res.ts["T_out"] - 273.15, color="#2a78d6", lw=2)
    ax[0].set_xscale("log"); ax[0].set_xlabel("time [years]"); ax[0].set_ylabel("outlet temperature [C]")
    ax[1].plot(t, res.ts["Q_wellhead"] / 1e6, color="#2a78d6", lw=2)
    ax[1].set_xscale("log"); ax[1].set_xlabel("time [years]"); ax[1].set_ylabel("thermal power at wellhead [MW]")
    if res.profiles:
        p = res.profiles[-1]
        ax[2].plot(p["T0"] - 273.15, p["tvd"], color="#898781", lw=1.5, label="undisturbed rock")
        ax[2].plot(p["T_a"] - 273.15, p["tvd"], color="#2a78d6", lw=2, label="annulus (down)")
        ax[2].plot(p["T_t"] - 273.15, p["tvd"], color="#eb6834", lw=2, label="tubing (up)")
        ax[2].invert_yaxis(); ax[2].set_xlabel("temperature [C]"); ax[2].set_ylabel("TVD [m]")
        ax[2].set_title(f"profiles at {p['t'] / DAY:.0f} d", fontsize=10); ax[2].legend(fontsize=8)
    for a in ax:
        a.grid(True, color="#e1e0d9", lw=0.8)
    fig.suptitle(sim.case.name)
    fig.tight_layout()
    fig.savefig(path, dpi=130)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="aharrsim", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="run a case file")
    r.add_argument("case")
    r.add_argument("--out", default=None)
    r.add_argument("--progress", action="store_true")
    r.add_argument("--plot", action="store_true")
    a = ap.parse_args(argv)
    case = load_case(a.case)
    out = a.out or os.path.splitext(a.case)[0] + "_results"
    os.makedirs(out, exist_ok=True)
    sim = Simulator(case)
    res = sim.run(progress=a.progress)
    ts = res.ts
    arr = np.column_stack([ts["t"] / DAY, ts["mdot"], ts["T_in"] - 273.15, ts["T_out"] - 273.15, ts["p_in"] / 1e6,
                           ts["p_out"] / 1e6, ts["choke_dp"] / 1e6, ts["x_out"], ts["Q_wellhead"] / 1e3,
                           ts["Q_rock"] / 1e3, ts["Q_riser_loss"] / 1e3, ts["flash_md"]])
    np.savetxt(os.path.join(out, "timeseries.csv"), arr, delimiter=",", fmt="%.6g",
               header="t_days,mdot_kg_s,T_in_C,T_out_C,p_in_MPa,p_wellhead_MPa,choke_dp_MPa,steam_quality_wellhead,"
                      "thermal_power_kW,heat_from_rock_kW,riser_to_annulus_loss_kW,flash_md_m", comments="")
    _write_profiles(res, os.path.join(out, "profiles.csv"))
    summ = _summary(sim, res)
    def _clean(o):
        if isinstance(o, dict):
            return {k: _clean(v) for k, v in o.items()}
        if isinstance(o, list):
            return [_clean(v) for v in o]
        if isinstance(o, (float, np.floating)):
            return None if not np.isfinite(o) else float(o)
        if isinstance(o, np.integer):
            return int(o)
        return o
    with open(os.path.join(out, "summary.json"), "w") as fh:
        json.dump(_clean(summ), fh, indent=2)
    if a.plot:
        _plot(sim, res, os.path.join(out, "summary.png"))
    for p in summ["points"]:
        print(f"t = {p['t_days']:9.1f} d   T_out = {p['T_out_C']:7.2f} C   Q = {p['thermal_power_MW']:7.3f} MWth   "
              f"exergy = {('%6.3f MW' % p['exergy_MW']) if p['exergy_MW'] == p['exergy_MW'] else '   n/a   '}   "
              f"p_in = {p['p_in_MPa']:6.2f} MPa   p_wh = {p['p_wellhead_MPa']:6.2f} MPa   "
              f"steam fraction at {p['p_delivery_MPa']:.2f} MPa = {p['steam_quality_delivery']:.3f}")
    print(f"cumulative heat {summ['cumulative_heat_GWh']:.1f} GWh;  results in {out}")
    if summ["steps_with_near_critical_riser_flow"]:
        print("WARNING: the flashing riser approaches critical (choked) flow in "
              f"{summ['steps_with_near_critical_riser_flow']} steps; raise the wellhead pressure or enlarge the tubing.")
    if summ["unconverged_steps"]:
        print(f"WARNING: {summ['unconverged_steps']} time steps did not fully converge; check the inputs or reduce dt_max.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
