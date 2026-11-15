"""Heat delivered to surface from 4-12 km coaxial closed loops: base cases and sensitivities.

Reference design (all depths)
  geotherm     15 C surface, 30 C/km (135 C at 4 km ... 375 C at 12 km)
  rock         k 2.8 W/m/K, 2650 kg/m3, 900 J/kg/K (constant)
  completion   9-5/8" 53.5 ppf casing to 3000 m in 12-1/4" hole, 7" 32 ppf liner to TD in
               8-1/2" hole, both cemented (k 1.0)
  centre pipe  3-1/2" x 4-1/2" vacuum-insulated tubing to TD, insulation k_app 0.03 W/m/K
  operation    10 kg/s water injected at 40 C down the annulus, 30 years continuous,
               delivery pressure 5 MPa at the wellhead (single phase)

Reported quantities
  Q_rock      m (h_out - h_in): net heat taken from the rock and delivered to surface
  Q_avail     m (h_out - h(p_del, T_in)): heat released by cooling the delivered stream back to
              the injection temperature (= Q_rock + pump work dissipated by friction)
  exergy      m [(h_out - h_r) - T0 (s_out - s_r)], r = (p_del, T_in), T0 = surface temperature
  pump        hydraulic power m (p_in - p_del) / rho_in (negative = thermosiphon surplus)
"""
from __future__ import annotations

import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from common import save, run, C, DAY, YEAR

YRS = [1, 5, 10, 20, 30]


def k_vosteen(k25, T_C, T_cap=500.0):
    """Vosteen & Schellschmidt (2003) crystalline-rock conductivity at T_C [C] from the 25 C value.
    Held flat above T_cap (the top of their measured range)."""
    k0 = 0.53 * k25 + 0.5 * np.sqrt(1.13 * k25 ** 2 - 0.42 * k25)
    T = np.minimum(T_C, T_cap)
    return k0 / (0.99 + T * (0.0030 - 0.0042 / k0))


def design(depth=6000.0, mdot=10.0, Tin=40.0, grad=30.0, k_rock=2.8, k_ins=0.03, vit_to=None,
           p_out=5.0, years=30.0, Ts=15.0, rhoc=(2650.0, 900.0), dz=50.0, tubing="vit", bigbore=False,
           kT=False):
    """kT=True: rock conductivity falls with the undisturbed temperature (Vosteen & Schellschmidt),
    k_rock being the 25 C value; 250 m layers to 2 km below TD."""
    vit_to = depth if vit_to is None else vit_to
    if bigbore:      # 9-5/8" casing to TD, 4-1/2" x 5-1/2" VIT
        vit = dict(id_in=3.958, layers=[dict(od_in=4.5, k=45.0), dict(od_in=5.0, k=k_ins), dict(od_in=5.5, k=45.0)])
        bare = dict(id_in=3.958, layers=[dict(od_in=4.5, k=45.0)])
    else:
        vit = dict(id_in=2.992, layers=[dict(od_in=3.5, k=45.0), dict(od_in=4.0, k=k_ins), dict(od_in=4.5, k=45.0)])
        bare = dict(id_in=2.992, layers=[dict(od_in=3.5, k=45.0)])
    if tubing == "steel":
        vit = dict(id_in=3.958, layers=[dict(od_in=4.5, k=45.0)])
    tub = [dict(top_m=0, bottom_m=vit_to, **vit)]
    if vit_to < depth:
        tub.append(dict(top_m=vit_to, bottom_m=depth, **bare))
    if bigbore:
        casing = [dict(top_m=0, bottom_m=depth, id_in=8.535,
                       layers=[dict(od_in=9.625, material="steel"), dict(od_in=12.25, material="cement")])]
    else:
        shoe = min(3000.0, 0.5 * depth)
        casing = [dict(top_m=0, bottom_m=shoe, id_in=8.535,
                       layers=[dict(od_in=9.625, material="steel"), dict(od_in=12.25, material="cement")]),
                  dict(top_m=shoe, bottom_m=depth, id_in=6.094,
                       layers=[dict(od_in=7.0, material="steel"), dict(od_in=8.5, material="cement")])]
    if kT:
        edges = np.arange(0.0, depth + 2000.0 + 1.0, 250.0)
        rock_layers = [dict(top_m=float(a), bottom_m=float(b), rho=rhoc[0], cp=rhoc[1],
                            k=float(k_vosteen(k_rock, Ts + grad * 0.5 * (a + b) / 1000.0)))
                       for a, b in zip(edges[:-1], edges[1:])]
        rock_layers[-1]["bottom_m"] = 1e6
    else:
        rock_layers = [dict(k=k_rock, rho=rhoc[0], cp=rhoc[1])]
    return dict(
        name=f"{depth / 1000:.0f} km, {mdot} kg/s",
        well=dict(depth_m=depth, casing=casing, tubing=tub),
        materials=dict(steel=dict(k=45.0, rho=7850.0, cp=490.0), cement=dict(k=1.0, rho=1900.0, cp=1000.0)),
        rock=dict(layers=rock_layers, surface_T_C=Ts,
                  geotherm=dict(gradient_C_per_km=grad)),
        operation=dict(pressure_bc=dict(outlet_MPa=p_out),
                       schedule=[dict(duration_days=years * 365.25, mdot_kg_s=mdot, T_in_C=Tin)]),
        numerics=dict(dz_m=dz, dt_initial_s=300.0, dt_max_days=60.0,
                      output_times_days=[y * 365.25 for y in YRS]),
    )


def summarise(label, d):
    t0 = time.time()
    sim, res = run(d)
    f = sim.fluid
    out = dict(label=label, depth=d["well"]["depth_m"], mdot=d["operation"]["schedule"][0]["mdot_kg_s"])
    T0 = d["rock"]["surface_T_C"] + 273.15
    T_inj = d["operation"]["schedule"][0]["T_in_C"] + 273.15
    rows = []
    for y in YRS:
        t = y * 365.25 * DAY
        i = int(np.argmin(np.abs(res.ts["t"] - t)))
        h_in, h_out = res.ts["h_in"][i], res.ts["h_out"][i]
        p_in, p_wh = res.ts["p_in"][i], res.ts["p_out"][i]
        choke = res.ts["choke_dp"][i]
        p_del = p_wh - choke
        m = res.ts["mdot"][i]
        dl = f.ph(np.array([p_del]), np.array([h_out]))
        x_del = float(dl.x[0]) if dl.two_phase[0] else (0.0 if dl.x[0] <= 0 else 1.0)
        h_r = f.h_pT(p_del, T_inj)[0]
        s_r = f.s_ph(p_del, h_r)[0]
        s_out = f.s_ph(p_del, h_out)[0]
        ex = m * ((h_out - h_r) - T0 * (s_out - s_r))
        rho_in = f.ph(np.array([p_in]), np.array([h_in])).rho[0]
        pump_kW = m * (p_in - p_del) / rho_in / 1e3
        pr = next(p for p in res.profiles if abs(p["t"] - t) < 1.0)
        rows.append(dict(years=y, T_out=float(C(res.ts["T_out"][i])), Q_MW=float(res.ts["Q_wellhead"][i] / 1e6),
                         h_in=float(h_in), h_out=float(h_out),
                         Q_avail_MW=float(m * (h_out - h_r) / 1e6), exergy_MW=float(ex / 1e6),
                         p_in_MPa=float(p_in / 1e6), p_wh_MPa=float(p_wh / 1e6), p_del_MPa=float(p_del / 1e6),
                         choke_MPa=float(choke / 1e6), pump_kW=float(pump_kW),
                         riser_loss_MW=float(res.ts["Q_riser_loss"][i] / 1e6),
                         x_wh=float(res.ts["x_out"][i]), x_delivery=x_del,
                         T_bottom=float(C(pr["T_a"][-1])), T_rock_bottom=float(C(pr["T0"][-1]))))
    out["rows"] = rows
    dt = np.diff(np.concatenate([[0.0], res.ts["t"]]))
    out["cum_GWh_30y"] = float(np.sum(dt * res.ts["Q_wellhead"]) / 3.6e12)
    out["mean_MW_30y"] = float(np.sum(dt * res.ts["Q_wellhead"]) / res.ts["t"][-1] / 1e6)
    keep = res.ts["t"] >= 0.5 * DAY
    out["series"] = dict(t_years=res.ts["t"][keep] / YEAR, T_out=C(res.ts["T_out"][keep]),
                         Q_MW=res.ts["Q_wellhead"][keep] / 1e6)
    pr = res.profiles[YRS.index(10)]
    out["profile_10y"] = dict(md=pr["md"], T_a=C(pr["T_a"]), T_t=C(pr["T_t"]), T0=C(pr["T0"]),
                              p_a=pr["p_a"] / 1e6, p_t=pr["p_t"] / 1e6, md_c=pr["md_c"],
                              q_rock=pr["q_rock"], q_riser=pr["q_riser"], x_t=pr["x_t"])
    out["wall_s"] = time.time() - t0
    out["cells"] = [sim.Nw, sim.N, sim.M]
    ts = res.ts
    out["flags"] = dict(unconverged_steps=int(res.meta.get("unconverged_steps", 0)),
                        choke_warning_steps=int(np.sum(ts["choke_warning"] > 0)),
                        max_energy_residual_kW=float(np.max(np.abs(ts["Q_wellhead"] - ts["Q_rock"] - ts["Q_res"]
                                                                   + ts["dE_fluid"] + ts["Q_ke"] + ts["Q_pe"])) / 1e3))
    r10 = rows[YRS.index(10)]
    print(f"{label:34s} 10y: Tout {r10['T_out']:6.1f} C  Q {r10['Q_MW']:6.2f} MW  Qav {r10['Q_avail_MW']:6.2f}  "
          f"Ex {r10['exergy_MW']:5.2f} MW  pump {r10['pump_kW']:7.1f} kW  x_del {r10['x_delivery']:.3f}  "
          f"({out['wall_s']:.0f}s)", flush=True)
    return out


def _job(args):
    label, kw = args
    try:
        return summarise(label, design(**kw))
    except Exception as e:  # noqa
        print(f"{label}: FAILED {type(e).__name__}: {e}", flush=True)
        return dict(label=label, failed=str(e), kw=kw)


def study(name, jobs, workers=2):
    with ProcessPoolExecutor(max_workers=workers) as ex:
        results = list(ex.map(_job, jobs))
    save(name, results)
    return results


def study_batch(named_jobs, workers=2):
    """Run jobs from several studies in one pool; save each study file as results arrive
    (so a long batch can be monitored and survives an interruption)."""
    from concurrent.futures import as_completed
    out = {name: {} for name, _, _ in named_jobs}
    order = {name: [] for name, _, _ in named_jobs}
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {}
        for name, label, kw in named_jobs:
            order[name].append(label)
            futs[ex.submit(_job, (label, kw))] = (name, label)
        for f in as_completed(futs):
            name, label = futs[f]
            out[name][label] = f.result()
            save(name, [out[name][lb] for lb in order[name] if lb in out[name]])
    return out


DEPTHS = [4000.0, 6000.0, 8000.0, 10000.0, 12000.0]

if __name__ == "__main__":
    which = sys.argv[1:] or ["depth_flow", "gradient", "insulation", "rock", "bigbore", "steam", "extended"]
    workers = 2
    if "depth_flow" in which:
        jobs = [(f"{d / 1000:.0f} km, {m:.0f} kg/s", dict(depth=d, mdot=m)) for d in DEPTHS for m in (5.0, 10.0, 15.0, 20.0)]
        study("D_depth_flow", jobs, workers)
    if "gradient" in which:
        jobs = [(f"{d / 1000:.0f} km, {g:.0f} C/km", dict(depth=d, grad=g)) for d in DEPTHS for g in (20.0, 25.0, 35.0, 40.0)]
        study("D_gradient", jobs, workers)
    if "insulation" in which:
        jobs = []
        for d in DEPTHS:
            jobs.append((f"{d / 1000:.0f} km, VIT k 0.01", dict(depth=d, k_ins=0.01)))
            jobs.append((f"{d / 1000:.0f} km, VIT k 0.10", dict(depth=d, k_ins=0.10)))
            jobs.append((f"{d / 1000:.0f} km, VIT top half only", dict(depth=d, vit_to=0.5 * d)))
            jobs.append((f"{d / 1000:.0f} km, bare steel tubing", dict(depth=d, tubing="steel")))
        study("D_insulation", jobs, workers)
    if "rock" in which:
        jobs = [(f"{d / 1000:.0f} km, k {k}", dict(depth=d, k_rock=k)) for d in DEPTHS for k in (2.0, 3.5)]
        study("D_rock", jobs, workers)
    if "bigbore" in which:
        jobs = [(f"{d / 1000:.0f} km big bore, {m:.0f} kg/s", dict(depth=d, mdot=m, bigbore=True))
                for d in DEPTHS for m in (10.0, 20.0, 40.0)]
        study("D_bigbore", jobs, workers)
    if "steam" in which:
        jobs = [(f"{d / 1000:.0f} km, {m:.0f} kg/s, {p} MPa delivery", dict(depth=d, mdot=m, p_out=p))
                for d in (10000.0, 12000.0) for m in (5.0, 10.0) for p in (0.5, 1.0, 2.0)]
        study("D_steam", jobs, workers)
    if "extended" in which:
        # gradients to 60 C/km (reference design, 10 kg/s); max ("high-high") case; k(T) check;
        # flow sweep at every gradient for the power case.  Deepest first for load balance.
        DEEP_FIRST = DEPTHS[::-1]
        jobs = [("D_gradient_ext", f"{d / 1000:.0f} km, {g:.0f} C/km", dict(depth=d, grad=g))
                for d in DEEP_FIRST for g in (45.0, 50.0, 55.0, 60.0)]
        jobs += [("D_max", f"{d / 1000:.0f} km max, {m:.0f} kg/s",
                  dict(depth=d, mdot=m, grad=60.0, k_rock=3.5, k_ins=0.01, bigbore=True))
                 for d in DEEP_FIRST for m in (5.0, 10.0, 20.0, 30.0, 40.0)]
        jobs += [("D_kT", f"{d / 1000:.0f} km, {g:.0f} C/km, k(T)", dict(depth=d, grad=g, kT=True))
                 for d in DEEP_FIRST for g in (30.0, 60.0)]
        jobs += [("D_gradient_flow", f"{d / 1000:.0f} km, {g:.0f} C/km, {m:.0f} kg/s", dict(depth=d, grad=g, mdot=m))
                 for g in (60.0, 55.0, 50.0, 45.0, 40.0, 35.0, 25.0, 20.0) for d in DEEP_FIRST for m in (5.0, 15.0, 20.0)]
        study_batch(jobs, workers)
