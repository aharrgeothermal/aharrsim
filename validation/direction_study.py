"""Flow direction and annulus size: inject down the tubing and produce up the annulus
('tubing_down') against the reference arrangement ('annulus_down'), over a matrix of
casing bores and vacuum-insulated tubing sizes.

All runs: 12 km vertical, single casing string to TD (cemented to the hole), VIT to TD with a
0.25 in radial vacuum gap at apparent k 0.03 W/m/K and a 0.25 in outer wall (the reference
3-1/2 x 4-1/2 construction at every size), 30 C/km from 15 C, rock k 2.8 W/m/K, 40 C injection,
5 MPa delivery, 30 years.  Base flow 10 kg/s; 5 and 20 kg/s and 60 C/km for three geometries;
bare steel tubing and a fixed tubing resistance as checks.

Results -> results/D_direction.json (one record per run, as deep_study.summarise plus geometry).
"""
from __future__ import annotations

import math
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

from common import save, INCH
from deep_study import summarise, YRS

K_INS = 0.03
CASINGS = {                      # ID, OD, hole (in)
    "7\"": (6.094, 7.0, 8.5),
    "9-5/8\"": (8.535, 9.625, 12.25),
    "10-3/4\"": (9.760, 10.75, 13.5),
    "13-3/8\"": (12.415, 13.375, 17.5),
}
TUBINGS = {                      # inner pipe ID, OD (in); outer pipe OD = inner OD + 1.0 in
    "2-3/8\"": (1.995, 2.375),
    "2-7/8\"": (2.441, 2.875),
    "3-1/2\"": (2.992, 3.5),
    "4-1/2\"": (3.958, 4.5),
    "5-1/2\"": (4.892, 5.5),
    "7\"": (6.276, 7.0),
}
COMBOS = [("7\"", t) for t in ("2-3/8\"", "2-7/8\"", "3-1/2\"", "4-1/2\"")] + \
         [("9-5/8\"", t) for t in ("2-7/8\"", "3-1/2\"", "4-1/2\"", "5-1/2\"", "7\"")] + \
         [("10-3/4\"", t) for t in ("3-1/2\"", "4-1/2\"", "5-1/2\"", "7\"")] + \
         [("13-3/8\"", t) for t in ("3-1/2\"", "4-1/2\"", "5-1/2\"", "7\"")]
SUBSET = [("7\"", "3-1/2\""), ("9-5/8\"", "4-1/2\""), ("13-3/8\"", "3-1/2\"")]


def geometry(casing, tubing, bare=False, R_wall=None):
    cid, cod, hole = CASINGS[casing]
    tid, tod = TUBINGS[tubing]
    if bare:
        tub = dict(id_in=tid, layers=[dict(od_in=tod, k=45.0)])
    else:
        tub = dict(id_in=tid, layers=[dict(od_in=tod, k=45.0), dict(od_in=tod + 0.5, k=K_INS), dict(od_in=tod + 1.0, k=45.0)])
    if R_wall is not None:
        tub["R_wall_K_m_W"] = R_wall
    r_ci, r_to, r_ti = 0.5 * cid * INCH, 0.5 * (tod + (0.0 if bare else 1.0)) * INCH, 0.5 * tid * INCH
    A_a, A_t = math.pi * (r_ci ** 2 - r_to ** 2), math.pi * r_ti ** 2
    R_t = 0.0 if bare else math.log((tod + 0.5) / tod) / (2 * math.pi * K_INS)
    if R_wall is not None:
        R_t = R_wall
    info = dict(casing=casing, tubing=tubing, casing_id_in=cid, hole_in=hole, tubing_id_in=tid,
                tubing_od_in=tod + (0.0 if bare else 1.0), A_annulus_m2=A_a, A_tubing_m2=A_t,
                Dh_annulus_m=2 * (r_ci - r_to), R_tubing_K_m_W=R_t, clearance_in=0.5 * (cid - tod - (0.0 if bare else 1.0)))
    cas = dict(top_m=0, bottom_m=12000.0, id_in=cid, layers=[dict(od_in=cod, material="steel"), dict(od_in=hole, material="cement")])
    return cas, dict(top_m=0, bottom_m=12000.0, **tub), info


def design(casing, tubing, direction, mdot=10.0, grad=30.0, bare=False, R_wall=None, depth=12000.0,
           Tin=40.0, p_out=5.0, k_rock=2.8, years=30.0):
    cas, tub, info = geometry(casing, tubing, bare, R_wall)
    cas["bottom_m"] = tub["bottom_m"] = depth
    d = dict(
        name=f"{casing} x {tubing} VIT, {direction}, {mdot} kg/s",
        well=dict(depth_m=depth, flow_direction=direction, casing=[cas], tubing=[tub]),
        materials=dict(steel=dict(k=45.0, rho=7850.0, cp=490.0), cement=dict(k=1.0, rho=1900.0, cp=1000.0)),
        rock=dict(layers=[dict(k=k_rock, rho=2650.0, cp=900.0)], surface_T_C=15.0, geotherm=dict(gradient_C_per_km=grad)),
        operation=dict(pressure_bc=dict(outlet_MPa=p_out),
                       schedule=[dict(duration_days=years * 365.25, mdot_kg_s=mdot, T_in_C=Tin)]),
        numerics=dict(dz_m=50.0, dt_initial_s=300.0, dt_max_days=60.0, output_times_days=[y * 365.25 for y in YRS]),
    )
    return d, info


def _job(args):
    label, kw = args
    try:
        d, info = design(**kw)
        out = summarise(label, d)
        out.update(info)
        out["direction"] = kw["direction"]
        out["grad"] = kw.get("grad", 30.0)
        out["bare"] = kw.get("bare", False)
        out["R_wall_override"] = kw.get("R_wall")
        return out
    except Exception as e:  # noqa
        print(f"{label}: FAILED {type(e).__name__}: {e}", flush=True)
        return dict(label=label, failed=str(e), kw=kw)


if __name__ == "__main__":
    jobs = []
    for c, t in COMBOS:
        for direction in ("tubing_down", "annulus_down"):
            jobs.append((f"{c} x {t}, {direction}, 10 kg/s", dict(casing=c, tubing=t, direction=direction)))
    for c, t in SUBSET:
        for direction in ("tubing_down", "annulus_down"):
            for m in (5.0, 20.0):
                jobs.append((f"{c} x {t}, {direction}, {m:.0f} kg/s", dict(casing=c, tubing=t, direction=direction, mdot=m)))
            jobs.append((f"{c} x {t}, {direction}, 10 kg/s, 60 C/km", dict(casing=c, tubing=t, direction=direction, grad=60.0)))
    for direction in ("tubing_down", "annulus_down"):
        jobs.append((f"9-5/8\" x 3-1/2\" bare steel, {direction}, 10 kg/s",
                     dict(casing="9-5/8\"", tubing="3-1/2\"", direction=direction, bare=True)))
    jobs.append(("13-3/8\" x 7\" R_t 0.708, tubing_down, 10 kg/s",
                 dict(casing="13-3/8\"", tubing="7\"", direction="tubing_down", R_wall=0.708)))
    jobs.append(("13-3/8\" x 7\" R_t 0.708, annulus_down, 10 kg/s",
                 dict(casing="13-3/8\"", tubing="7\"", direction="annulus_down", R_wall=0.708)))
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    results, order = {}, [j[0] for j in jobs]
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(_job, j): j[0] for j in jobs}
        for f in as_completed(futs):
            results[futs[f]] = f.result()
            save("D_direction", [results[k] for k in order if k in results])
    print(f"done: {len(results)} runs")
