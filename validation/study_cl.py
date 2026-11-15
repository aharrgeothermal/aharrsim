"""Study part 1: uncertainty of the conduction-fed closed loops (types 7, 8, 20 of dev_types.py).

One-at-a-time ranges (tornado), design levers, two physics cases beyond pure conduction (temperature-dependent rock
conductivity; a convecting permeable zone represented by an effective conductivity Nu*k), and a Latin-hypercube
Monte Carlo on the subsurface and completion uncertainties.

    python study_cl.py oat      # one-at-a-time and physics cases  -> results/ST1_oat.json
    python study_cl.py mc       # Latin hypercube                   -> results/ST1_mc.json
"""
from __future__ import annotations

import json
import os
import signal
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import CoolProp.CoolProp as CP

from common import save, load, RESULTS
from aharrsim.systems import coaxial_well, uloop, run_case, Series, YEAR
from aharrsim import plant

IN = 0.0254


def k_of_T(k25, T_C):
    """Vosteen & Schellschmidt (2003) crystalline rock: k(T) = k0 / (0.99 + T (0.0030 - 0.0042 / k0)), k0 at 0 C."""
    k0 = 0.53 * k25 + 0.5 * np.sqrt(1.13 * k25 ** 2 - 0.42 * k25)
    return k0 / (0.99 + T_C * (0.0030 - 0.0042 / k0))


def rock_layers(k, rho, cp, depth, grad, T_s=15.0, kT=False, zone=None, nu=1.0, dz=250.0):
    """Layered rock: optional k(T) evaluated at the undisturbed temperature of each layer, optional convecting zone
    (top, bottom) whose conductivity is multiplied by nu."""
    if not kT and (zone is None or nu == 1.0):
        return [dict(k=k, rho=rho, cp=cp)]
    edges = list(np.arange(0.0, depth + dz, dz))
    if zone:
        edges = sorted(set(edges + [zone[0], zone[1]]))
    out = []
    for a, b in zip(edges[:-1], edges[1:]):
        zm = 0.5 * (a + b)
        kk = float(k_of_T(k, T_s + grad * zm / 1000.0)) if kT else k
        if zone and zone[0] <= zm < zone[1]:
            kk *= nu
        out.append(dict(top_m=a, bottom_m=b if b < edges[-1] else 1e6, k=kk, rho=rho, cp=cp))
    return out


def power(m, h, p, T_amb=15.0):
    """Gross MWe: best of practical binary, flash (condenser at ambient + 25 K) and steam turbine."""
    T = CP.PropsSI("T", "H", h, "P", p, "Water") - 273.15
    w = [0.0]
    if T > 60:
        try:
            w.append(plant.orc_practical(m, h, p, T_amb)["W_net"] / 1e6)
        except Exception:
            pass
    if T > 150:
        try:
            w.append(plant.flash(m, h, p, T_cond_C=T_amb + 25.0, double=True)["W_net"] / 1e6)
        except Exception:
            pass
    return max(w)


def metrics(res, m, T_in):
    ts = res.ts
    out = {}
    for y in (1, 10, 30):
        t = y * YEAR
        h = float(np.interp(t, ts["t"], ts["h_out"])); p = float(np.interp(t, ts["t"], ts["p_out"]))
        T = float(np.interp(t, ts["t"], ts["T_out"])) - 273.15
        h_in = CP.PropsSI("H", "T", T_in + 273.15, "P", p, "Water")
        out[y] = dict(T=T, MWth=m * (h - h_in) / 1e6, MWe=power(m, h, p, 15.0), MWe35=power(m, h, p, 35.0), MWe45=power(m, h, p, 45.0),
                      h=h, p=p, pump_kW=float(np.interp(t, ts["t"], ts["pump_kW"])))
    out["unconverged"] = int(res.meta.get("unconverged_steps", 0))
    return out


# ------------------------------------------------------------------------------------------------ cases
BASE = {
    7: dict(grad=80.0, k=2.5, cp=900.0, rho=2600.0, ins_k=0.03, mdot=10.0, T_in=40.0, depth=3000.0, kT=False, nu=1.0),
    8: dict(grad=55.0, k=2.8, cp=900.0, rho=2600.0, mdot=100.0, T_in=40.0, n_lat=12, L2=5000.0, spacing=80.0, d=0.216, kT=False, nu=1.0, L1=4500.0),
    20: dict(grad=24.25, k=2.8, cp=950.0, rho=2700.0, ins_k=0.03, mdot=20.0, T_in=40.0, depth=20000.0, bore=0.30, kT=False),
}
NAME = {7: "GreenLoop-type coaxial retrofit, 3 km", 8: "Eavor-type multilateral loop", 20: "Quaise-type 20 km coaxial loop"}


def run7(P):
    depth = P["depth"]; half = depth / 2
    casing = [dict(top_m=0, bottom_m=half, id_in=8.535, layers=[dict(od_in=9.625, material="steel"), dict(od_in=12.25, material="cement")]),
              dict(top_m=half, bottom_m=depth, id_in=6.094, layers=[dict(od_in=7.0, material="steel"), dict(od_in=8.5, material="cement")])]
    tubing = [dict(top_m=0, bottom_m=depth, id_in=2.992, layers=[dict(od_in=3.5, k=45), dict(od_in=4.0, k=P["ins_k"]), dict(od_in=4.5, k=45)])]
    if P.get("big"):       # 9-5/8 in casing to TD and 4-1/2 x 5-1/2 in vacuum-insulated tubing for the high-flow convective cases
        casing = [dict(top_m=0, bottom_m=depth, id_in=8.535, layers=[dict(od_in=9.625, material="steel"), dict(od_in=12.25, material="cement")])]
        tubing = [dict(top_m=0, bottom_m=depth, id_in=3.958, layers=[dict(od_in=4.5, k=45), dict(od_in=5.0, k=P["ins_k"]), dict(od_in=5.5, k=45)])]
    rock = rock_layers(P["k"], P["rho"], P["cp"], depth, P["grad"], kT=P["kT"], zone=(half, depth) if P["nu"] != 1.0 else None, nu=P["nu"])
    d = coaxial_well("cl7", depth, casing, tubing, P["mdot"], P["T_in"], 2.0e6, rock, 15.0, dict(gradient_C_per_km=P["grad"]), 30.0, output_years=(1, 30))
    sim, res = run_case(d)
    return metrics(res, P["mdot"], P["T_in"])


def run8(P):
    L1, L2 = P["L1"], P["L2"]
    rock = rock_layers(P["k"], P["rho"], P["cp"], L1 + 500.0, P["grad"], kT=P["kT"], zone=(L1 - 250.0, L1 + 500.0) if P["nu"] != 1.0 else None, nu=P["nu"])
    ch = uloop("cl8", L1, L2, P["d"], P["mdot"], P["T_in"], 2.0e6, rock, 15.0, dict(gradient_C_per_km=P["grad"]), 30.0, n_lat=int(P["n_lat"]),
               layers=None, roughness=1e-4, dz=50.0, lat_far_field="adiabatic", lat_r_far=P["spacing"] / 2.0)
    out = metrics(ch["producer"][1], P["mdot"], P["T_in"])
    out["T_toe30"] = ch["T_toe"](30 * YEAR) - 273.15
    out["hole_km"] = (2 * L1 + int(P["n_lat"]) * L2) / 1000.0
    return out


def run20(P):
    depth = P["depth"]; bore = P["bore"]
    casing = [dict(top_m=0, bottom_m=depth, id_m=bore, layers=[dict(od_m=bore + 0.04, k=1.4, rho=2400.0, cp=900.0, name="glass")])]
    tubing = [dict(top_m=0, bottom_m=depth, id_in=4.892, layers=[dict(od_in=5.5, k=45), dict(od_in=6.0, k=P["ins_k"]), dict(od_in=6.5, k=45)])]
    rock = rock_layers(P["k"], P["rho"], P["cp"], depth, P["grad"], kT=P["kT"], dz=1000.0)
    d = coaxial_well("cl20", depth, casing, tubing, P["mdot"], P["T_in"], 10.0e6, rock, 15.0, dict(gradient_C_per_km=P["grad"]), 30.0, dz=100.0,
                     output_years=(1, 30))
    sim, res = run_case(d)
    return metrics(res, P["mdot"], P["T_in"])


RUN = {7: run7, 8: run8, 20: run20}

# one-at-a-time: (variable, low, high, label, kind)   kind: 'unc' subsurface / completion uncertainty, 'des' design lever
OAT = {
    7: [("grad", 60.0, 100.0, "Gradient 60 / 100 °C/km", "unc"), ("k", 2.0, 3.5, "Rock conductivity 2.0 / 3.5 W/m/K", "unc"),
        ("cp", 765.0, 1035.0, "Rock heat capacity −15 / +15 %", "unc"), ("ins_k", 0.10, 0.01, "Tubing insulation k 0.10 / 0.01 W/m/K", "unc"),
        ("mdot", 5.0, 20.0, "Flow 5 / 20 kg/s", "des"), ("T_in", 60.0, 30.0, "Inlet 60 / 30 °C", "des"), ("depth", 2500.0, 4000.0, "Depth 2.5 / 4 km", "des")],
    8: [("grad", 45.0, 65.0, "Gradient 45 / 65 °C/km", "unc"), ("k", 2.0, 4.0, "Rock conductivity 2.0 / 4.0 W/m/K", "unc"),
        ("cp", 765.0, 1035.0, "Rock heat capacity −15 / +15 %", "unc"), ("mdot", 50.0, 125.0, "Flow 50 / 125 kg/s", "des"),
        ("T_in", 60.0, 30.0, "Inlet 60 / 30 °C", "des"), ("n_lat", 8, 16, "Laterals 8 / 16", "des"), ("L2", 3500.0, 7000.0, "Lateral length 3.5 / 7 km", "des"),
        ("spacing", 40.0, 160.0, "Lateral spacing 40 / 160 m", "des"), ("d", 0.156, 0.311, "Lateral diameter 6⅛ / 12¼ in", "des")],
    20: [("grad", 19.25, 29.25, "TD temperature 400 / 600 °C", "unc"), ("k", 2.0, 4.0, "Rock conductivity 2.0 / 4.0 W/m/K", "unc"),
         ("cp", 807.0, 1093.0, "Rock heat capacity −15 / +15 %", "unc"), ("ins_k", 0.10, 0.01, "Tubing insulation k 0.10 / 0.01 W/m/K", "unc"),
         ("mdot", 10.0, 40.0, "Flow 10 / 40 kg/s", "des"), ("T_in", 60.0, 30.0, "Inlet 60 / 30 °C", "des"), ("depth", 15000.0, 25000.0, "Depth 15 / 25 km", "des"),
         ("bore", 0.20, 0.45, "Bore 0.20 / 0.45 m", "des")],
}
PHYS = {
    7: [("kT", True, "k falls with temperature (Vosteen)"), ("ins_k", 0.3, "Insulation lost (k 0.3)"), ("nu", 3.0, "Convecting zone below 1.5 km, Nu 3"),
        ("nu", 10.0, "Convecting zone, Nu 10"), ("nu", 30.0, "Convecting / steam zone, Nu 30"), ("nu", 100.0, "Steam zone, Nu 100")],
    8: [("kT", True, "k falls with temperature (Vosteen)"), ("nu", 3.0, "Laterals in a convecting formation, Nu 3"), ("nu", 10.0, "Convecting formation, Nu 10")],
    20: [("kT", True, "k falls with temperature (Vosteen)"), ("ins_k", 0.3, "Insulation lost (k 0.3)")],
}
SWEEP = {7: ("mdot", [2.5, 5, 7.5, 10, 15, 20, 30]), 8: ("mdot", [30, 50, 75, 100]), 20: ("mdot", [5, 10, 15, 20, 30, 40, 60])}
# Monte Carlo: triangular(low, mode, high); ins_k log-triangular; T_amb uniform (plant only)
MC = {7: dict(grad=(60, 80, 100), k=(2.0, 2.5, 3.5), cp=(765, 900, 1035), ins_k=(0.01, 0.03, 0.10)),
      8: dict(grad=(45, 55, 65), k=(2.0, 2.8, 4.0), cp=(765, 900, 1035)),
      20: dict(grad=(19.25, 24.25, 29.25), k=(2.0, 2.8, 4.0), cp=(807, 950, 1093), ins_k=(0.01, 0.03, 0.10))}
N_MC = {7: 80, 8: 40, 20: 40}


JOB_TIMEOUT_S = 240


class JobTimeout(BaseException):       # BaseException so the solver's own retry handlers cannot swallow it
    pass


def _timeout(signum, frame):
    raise JobTimeout()


def job(cid, tag, P):
    t0 = time.time()
    try:
        signal.signal(signal.SIGALRM, _timeout); signal.setitimer(signal.ITIMER_REAL, JOB_TIMEOUT_S, 5.0)     # some high-friction loop cases crawl; give up on them
        r = RUN[cid](P)
        signal.setitimer(signal.ITIMER_REAL, 0)
        r.update(case=cid, tag=tag, P={k: (v if not isinstance(v, (np.floating, np.integer)) else float(v)) for k, v in P.items()}, wall_s=time.time() - t0)
        return r
    except JobTimeout:
        signal.setitimer(signal.ITIMER_REAL, 0)
        return dict(case=cid, tag=tag, P=P, failed=f"no result in {JOB_TIMEOUT_S} s")
    except Exception as e:
        signal.setitimer(signal.ITIMER_REAL, 0)
        return dict(case=cid, tag=tag, P=P, failed=f"{type(e).__name__}: {e}", tb=traceback.format_exc()[-800:])


def run_jobs(jobs, name, workers=2):
    out = []
    path = os.path.join(RESULTS, name + ".json")
    if os.path.exists(path) and "--fresh" not in sys.argv:      # resume: keep finished cases
        out = [r for r in json.load(open(path)) if "failed" not in r]
        done = {(r["case"], r["tag"]) for r in out}
        jobs = [j for j in jobs if (j[0], j[1]) not in done]
        print(f"resuming {name}: {len(out)} done, {len(jobs)} to run", flush=True)
    with ProcessPoolExecutor(workers) as ex:
        futs = [ex.submit(job, *j) for j in jobs]
        for i, f in enumerate(as_completed(futs)):
            r = f.result(); out.append(r)
            print(f"[{i + 1}/{len(jobs)}] case {r['case']} {r['tag']}: " + (r.get("failed") or f"T30 {(r.get(30) or r['30'])['T']:.1f} C, {(r.get(30) or r['30'])['MWth']:.2f} MWth, {(r.get(30) or r['30'])['MWe']:.3f} MWe ({r['wall_s']:.0f} s)"), flush=True)
            if (i + 1) % 10 == 0:
                save(name, out)
    save(name, out)
    return out


def tri(u, lo, mode, hi):
    fc = (mode - lo) / (hi - lo)
    return np.where(u < fc, lo + np.sqrt(u * (hi - lo) * (mode - lo)), hi - np.sqrt((1 - u) * (hi - lo) * (hi - mode)))


def main():
    what = sys.argv[1] if len(sys.argv) > 1 else "oat"
    cases = [int(a) for a in sys.argv[2:]] or [7, 8, 20]
    suffix = ""
    if what.endswith("30"):          # type 8 at the flow that suits electricity (30 kg/s) instead of the example's 100 kg/s: "oat30 8", "mc30 8"
        what = what[:-2]; cases = [8]; suffix = "_m30"
        BASE[8]["mdot"] = 30.0
        OAT[8] = [("mdot", 20.0, 50.0, "Flow 20 / 50 kg/s", "des") if o[0] == "mdot" else o for o in OAT[8]]
        SWEEP[8] = ("mdot", [])
    if what == "oat":
        jobs = []
        for c in cases:
            jobs.append((c, "base", dict(BASE[c])))
            for var, lo, hi, label, kind in OAT[c]:
                for side, v in (("low", lo), ("high", hi)):
                    P = dict(BASE[c]); P[var] = v
                    jobs.append((c, f"oat|{var}|{side}", P))
            for var, v, label in PHYS[c]:
                P = dict(BASE[c]); P[var] = v
                jobs.append((c, f"phys|{var}|{v}", P))
            var, vals = SWEEP[c]
            for v in vals:
                P = dict(BASE[c]); P[var] = float(v)
                jobs.append((c, f"sweep|{var}|{v}", P))
        jobs.sort(key=lambda j: {7: 0, 20: 1, 8: 2}[j[0]])
        run_jobs(jobs, "ST1_oat" + ("" if len(cases) == 3 else "_" + "_".join(map(str, cases))) + suffix)
    elif what == "opt":     # convecting / steam zone: does more or less flow use it better? (larger tubing and casing, P["big"], was tried and is worse)
        jobs = []
        for nu in (10.0, 30.0, 100.0):
            for m in (5.0, 10.0, 20.0):
                P = dict(BASE[7]); P.update(nu=nu, mdot=m)
                jobs.append((7, f"opt|nu{nu:g}|m{m:g}", P))
        for nu, m in ((3.0, 100.0), (10.0, 100.0), (30.0, 100.0)):
            P = dict(BASE[8]); P.update(nu=nu, mdot=m)
            jobs.append((8, f"opt|nu{nu:g}|m{m:g}", P))
        run_jobs(jobs, "ST1_opt")
    elif what == "mc":
        rng = np.random.default_rng(20261113)
        jobs = []
        for c in cases:
            n = N_MC[c]; vars_ = list(MC[c])
            U = (np.array([rng.permutation(n) for _ in vars_]).T + rng.random((n, len(vars_)))) / n     # Latin hypercube
            for i in range(n):
                P = dict(BASE[c]); P["kT"] = True
                for j, v in enumerate(vars_):
                    lo, mo, hi = MC[c][v]
                    P[v] = float(np.exp(tri(U[i, j], np.log(lo), np.log(mo), np.log(hi)))) if v == "ins_k" else float(tri(U[i, j], lo, mo, hi))
                jobs.append((c, f"mc|{i}", P))
        jobs.sort(key=lambda j: {7: 0, 20: 1, 8: 2}[j[0]])
        run_jobs(jobs, "ST1_mc" + ("" if len(cases) == 3 else "_" + "_".join(map(str, cases))) + suffix)


if __name__ == "__main__":
    main()
