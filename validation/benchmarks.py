"""Code-to-code benchmarks for deep wells (4-12 km).

B1  GeoCLUSTER / CLGWG tabulated database (Sandia Sierra/Aria model, IAPWS
    water via CoolProp, injection at 20 MPa).  L-shaped coaxial wells with
    L1 = 4-5 km vertical + 1 km lateral, open-hole annulus, inner pipe wall
    19.2 mm thick with k = 0.06 W/m/K, annulus/pipe area ratio 1.
B2  SBT (slender-body theory) model of Beckers et al., as implemented in
    GeoCLUSTER (pnnl/GeoCLUSTER, sbt_v27.py), vertical coaxial wells 4-10 km.
B3  Published single-code results with fully specified inputs:
      Nalla et al. (2004/2005) TETRAD, 5593 m
      Kubacka & Karayiannis (2024) WellTH, 5000 m
      Beckers et al. (2022) SBT, 4000 m (Case 6)
"""
from __future__ import annotations

import math
import os
import sys
import time

import numpy as np

from common import save, run, C, DAY, YEAR

# Public GeoCLUSTER code and database (https://github.com/pnnl/GeoCLUSTER):
#   git clone https://github.com/pnnl/GeoCLUSTER validation/external/GeoCLUSTER
# or point GEOCLUSTER_DIR at an existing clone.
_HERE = os.path.dirname(os.path.abspath(__file__))
GEOCLUSTER = os.environ.get("GEOCLUSTER_DIR", os.path.join(_HERE, "external", "GeoCLUSTER"))
DB = os.path.join(GEOCLUSTER, "dash_app", "data", "clgs_results_final.h5")


# ----------------------------------------------------------------------------
# B1 GeoCLUSTER database
# ----------------------------------------------------------------------------
class GeoClusterDB:
    def __init__(self, path=DB):
        import h5py
        f = h5py.File(path, "r")
        base = "/coaxial/H2O/"
        self.inp = {k: f[base + "input/" + k][:] for k in ("mdot", "L2", "L1", "grad", "D", "T_i", "k_rock", "time")}
        self.fixed = {k: f["/coaxial/fixed_params/" + k][()] for k in f["/coaxial/fixed_params"].keys()}
        shape = (26, 20, 9, 5, 3, 3, 3, 161)
        out = {}
        We = f[base + "output/We"][:]
        valid = np.argwhere(np.isfinite(We.flatten()))[:, 0]
        for name in ("Tout", "Pout"):
            U = f[base + f"output/{name}/U"][:]
            s = f[base + f"output/{name}/sigma"][:]
            Vt = f[base + f"output/{name}/Vt"][:]
            Mk = U @ np.diag(s) @ Vt
            full = np.full((shape[-1], int(np.prod(shape[:-1]))), np.nan, dtype=np.float32)
            full[:, valid] = Mk
            out[name] = full.T.reshape(shape)
        self.Tout, self.Pout = out["Tout"], out["Pout"]
        f.close()

    def idx(self, name, value):
        arr = self.inp[name]
        i = int(np.argmin(np.abs(arr - value)))
        if abs(arr[i] - value) > 1e-6 * max(1.0, abs(value)):
            raise ValueError(f"{name}={value} not on the database grid {arr}")
        return i

    def series(self, mdot, L2, L1, grad, D, Ti_K, k):
        ii = (self.idx("mdot", mdot), self.idx("L2", L2), self.idx("L1", L1), self.idx("grad", grad),
              self.idx("D", D), self.idx("T_i", Ti_K), self.idx("k_rock", k))
        return self.inp["time"], self.Tout[ii], self.Pout[ii]


def geocluster_case(mdot, L1, L2, grad, D, Ti_C, k, years=40.0, dz=50.0):
    R = D / 2.0
    th = 0.0192
    # annulus area = pipe area:  R^2 - (ri + th)^2 = ri^2
    a, b, cc = 2.0, 2.0 * th, th * th - R * R
    ri = (-b + math.sqrt(b * b - 4 * a * cc)) / (2 * a)
    return dict(
        name=f"GeoCLUSTER L1={L1} L2={L2} m={mdot}",
        well=dict(depth_m=L1, lateral_m=L2, roughness_m=2.5e-5,
                  casing=[dict(id_m=D, layers=[])],
                  tubing=[dict(id_m=2 * ri, layers=[dict(od_m=2 * (ri + th), k=0.06)])]),
        rock=dict(layers=[dict(k=k, rho=2750.0, cp=790.0)], surface_T_C=25.0,
                  geotherm=dict(gradient_C_per_km=grad * 1000.0)),
        operation=dict(pressure_bc=dict(inlet_MPa=20.0),
                       schedule=[dict(duration_days=years * 365.0, mdot_kg_s=mdot, T_in_C=Ti_C)]),
        numerics=dict(dz_m=dz, dt_initial_s=300.0, dt_max_days=60.0, r_far_m=150.0,
                      output_times_days=[365.0 * y for y in (1, 5, 10, 20, 30, 40)]),
    )


def b1_geocluster():
    db = GeoClusterDB()
    print("fixed params:", {k: float(v) for k, v in db.fixed.items()})
    points = []
    for L1 in (4000.0, 5000.0):
        for grad in (0.03, 0.05, 0.07):
            for mdot in (5.0, 20.2, 39.2):
                points.append(dict(mdot=mdot, L1=L1, L2=1000.0, grad=grad, D=0.2159, Ti=30.0, k=3.0))
    # rock conductivity and borehole size sensitivity at 4.5 km
    for k in (1.5, 4.5):
        points.append(dict(mdot=12.6, L1=4500.0, L2=1000.0, grad=0.05, D=0.2159, Ti=30.0, k=k))
    points.append(dict(mdot=12.6, L1=4500.0, L2=1000.0, grad=0.05, D=0.3302, Ti=45.0, k=3.0))
    rows = []
    for p in points:
        t_db, T_db, P_db = db.series(p["mdot"], p["L2"], p["L1"], p["grad"], p["D"], p["Ti"] + 273.15, p["k"])
        if not np.all(np.isfinite(T_db[1:])):
            print("  database point pruned (NaN):", p)
            continue
        t0 = time.time()
        sim, res = run(geocluster_case(p["mdot"], p["L1"], p["L2"], p["grad"], p["D"], p["Ti"], p["k"]))
        yrs = np.array([1, 5, 10, 20, 30, 40], float)
        Tm = np.array([C(res.at(y * 365.0 * DAY, "T_out")) for y in yrs])
        Td = np.interp(yrs, t_db, T_db) - 273.15
        Pm = np.array([res.at(y * 365.0 * DAY, "p_out") / 1e6 for y in yrs])
        Pd = np.interp(yrs, t_db, P_db) / 1e6
        T_in = p["Ti"]
        dT_rel = (Tm - T_in) / (Td - T_in) - 1.0
        row = dict(p, years=yrs, T_model=Tm, T_db=Td, dT_K=Tm - Td, rise_rel=dT_rel, p_out_model=Pm,
                   p_out_db=Pd, series_t=res.t_days / 365.0, series_T=C(res.ts["T_out"]),
                   db_t=t_db, db_T=T_db - 273.15)
        rows.append(row)
        print(f"  L1={p['L1']:.0f} grad={p['grad']:.2f} m={p['mdot']:5.1f} k={p['k']} D={p['D']}: "
              + " ".join(f"{y:.0f}y {a:6.2f}/{b:6.2f}" for y, a, b in zip(yrs, Tm, Td))
              + f"   pout {Pm[-1]:.2f}/{Pd[-1]:.2f} MPa  ({time.time() - t0:.0f}s)", flush=True)
    save("B1_geocluster", dict(fixed={k: float(v) for k, v in db.fixed.items()}, rows=rows))


# ----------------------------------------------------------------------------
# B2 SBT (GeoCLUSTER implementation), vertical coaxial
# ----------------------------------------------------------------------------
def _import_sbt():
    app = os.path.join(GEOCLUSTER, "dash_app")
    os.environ["deployment_type"] = "local"
    if app not in sys.path:
        sys.path.insert(0, app)
    import paths  # noqa
    for k, v in list(paths.inpath_dict.items()):
        if isinstance(v, str) and v.startswith("data/"):
            paths.inpath_dict[k] = os.path.join(app, v)
    import sbt_v27
    return sbt_v27


def sbt_case(depth, mdot, k_ins, Tin=40.0, Ts=15.0, grad=0.03, k=3.0, rho=2700.0, cp=850.0,
             D=0.2159, ri=0.0381, th=0.0191, version=1, p_in_bar=100.0):
    """Run GeoCLUSTER's SBT (v1: constant water properties; v2: T,p-dependent)."""
    sbt = _import_sbt()
    hp = {}
    if version == 2:
        hp = dict(HYPERPARAM1=p_in_bar, HYPERPARAM2=1e-6, HYPERPARAM3=1, HYPERPARAM4=1e-5, HYPERPARAM5=15)
    out = sbt.run_sbt(sbt_version=version, mesh_fineness=0, accuracy=3, clg_configuration=1,
                      mdot=mdot, Tinj=Tin, fluid=1, DrillingDepth_L1=depth / 1000.0, HorizontalExtent_L2=0.0,
                      Diameter1=D, Diameter2=2 * ri, PipeParam3=th, PipeParam4=k_ins,
                      PipeParam5="Inject in Annulus", Tsurf=Ts, GeoGradient=grad, k_m=k, c_m=cp, rho_m=rho, **hp)
    t_years = np.asarray(out[0], float).ravel()
    T_K = np.asarray(out[1], float).ravel()
    return t_years, T_K - 273.15


def coax_vertical(depth, mdot, k_ins, Tin=40.0, Ts=15.0, grad=0.03, k=3.0, rho=2700.0, cp=850.0,
                  D=0.2159, ri=0.0381, th=0.0191, years=40.0, constant=True, p_in_MPa=10.0):
    d = dict(
        name=f"vertical {depth} m",
        well=dict(depth_m=depth, roughness_m=1e-6, casing=[dict(id_m=D, layers=[])],
                  tubing=[dict(id_m=2 * ri, layers=[dict(od_m=2 * (ri + th), k=k_ins)])]),
        rock=dict(layers=[dict(k=k, rho=rho, cp=cp)], surface_T_C=Ts, geotherm=dict(gradient_C_per_km=grad * 1000)),
        operation=dict(pressure_bc=dict(inlet_MPa=p_in_MPa),
                       schedule=[dict(duration_days=years * 365.0, mdot_kg_s=mdot, T_in_C=Tin)]),
        numerics=dict(dz_m=50.0, dt_initial_s=300.0, dt_max_days=60.0, r_far_m=200.0,
                      output_times_days=[365.0 * y for y in (1, 5, 10, 20, 40)]),
    )
    if constant:
        d["fluid"] = {"constant": {"rho": 1000.0, "cp": 4200.0, "mu": 6e-4, "k": 0.68}}
        d["options"] = dict(gravity_work=False, include_pressure=False)
    return d


def b2_sbt():
    rows = []
    yrs = np.array([1, 5, 10, 20, 40], float)
    cases = [(d, m, kins) for d in (4000.0, 6000.0, 8000.0, 10000.0) for (m, kins) in ((10.0, 0.06), (20.0, 0.06), (10.0, 0.5))]
    for depth, mdot, k_ins in cases:
        for version in (1, 2):
            if version == 2 and k_ins != 0.06:
                continue
            t0 = time.time()
            try:
                t_s, T_s = sbt_case(depth, mdot, k_ins, version=version)
                ok = np.all(np.isfinite(T_s)) and T_s.min() > -50 and T_s.max() < 600
            except Exception as e:  # noqa
                print(f"  SBT v{version} failed depth={depth} m={mdot} k={k_ins}: {type(e).__name__}: {e}")
                continue
            if not ok:
                print(f"  SBT v{version} returned invalid temperatures depth={depth} m={mdot} k={k_ins}")
                continue
            sim, res = run(coax_vertical(depth, mdot, k_ins, constant=(version == 1)))
            Tm = np.array([C(res.at(y * 365.0 * DAY, "T_out")) for y in yrs])
            Tsb = np.interp(yrs, t_s, T_s)
            rows.append(dict(depth=depth, mdot=mdot, k_ins=k_ins, sbt_version=version, years=yrs,
                             T_model=Tm, T_sbt=Tsb, dT=Tm - Tsb, rise_rel=(Tm - 40.0) / (Tsb - 40.0) - 1.0,
                             sbt_t=t_s, sbt_T=T_s, model_t=res.t_days / 365.0, model_T=C(res.ts["T_out"])))
            print(f"  SBT v{version} {depth:6.0f} m m={mdot:4.1f} k_ins={k_ins:4.2f}: " +
                  " ".join(f"{y:.0f}y {a:6.2f}/{b:6.2f}" for y, a, b in zip(yrs, Tm, Tsb)) +
                  f"  ({time.time() - t0:.0f}s)", flush=True)
    save("B2_sbt", rows)


# ----------------------------------------------------------------------------
# B3 published single-code results
# ----------------------------------------------------------------------------
def nalla_case(mdot=6.3, k_ins=0.07, d_hole=0.311, years=3.0, depth=5593.0, tub_od=0.0889, constant=False):
    """Nalla, Shook, Mines & Bloomfield (2004) INEEL/EXT-03-01433 (Geothermics 34:330, 2005), TETRAD.

    5593 m vertical; casing 228.6/244.5 mm cemented (k 0.87) over the top 1524 m,
    open hole (311 mm) below; tubing 76.2/88.9 mm steel (k 44.83) insulated to
    101.6 mm OD (k 0.07); rock k 1.73, rho*c 1810.7 kJ/m3K; surface 26.7 C,
    350 C at TD (gradient 0.1 W/m2 / 1.73); 6.3 kg/s at 26.7 C down the annulus.
    """
    ins_od = tub_od + 2 * 0.00635
    d = dict(
        name="Nalla et al.",
        well=dict(depth_m=depth,
                  casing=[dict(top_m=0, bottom_m=min(1524.0, depth - 1), id_m=0.2286,
                               layers=[dict(od_m=0.244475, material="steel"), dict(od_m=d_hole, material="cement")]),
                          dict(top_m=min(1524.0, depth - 1), bottom_m=depth, id_m=d_hole, layers=[])],
                  tubing=[dict(id_m=tub_od - 2 * 0.00635, layers=[dict(od_m=tub_od, k=44.83), dict(od_m=ins_od, k=max(k_ins, 1e-6))])]),
        materials=dict(steel=dict(k=44.83, rho=7850.0, cp=3836.93e3 / 7850.0),
                       cement=dict(k=0.87, rho=1900.0, cp=1260.8e3 / 1900.0)),
        rock=dict(layers=[dict(k=1.73, rho=2000.0, cp=1810.7e3 / 2000.0)], surface_T_C=26.7,
                  geotherm=dict(gradient_C_per_km=(350.0 - 26.7) / depth * 1000.0)),
        operation=dict(pressure_bc=dict(outlet_MPa=1.0),
                       schedule=[dict(duration_days=years * 365.0, mdot_kg_s=mdot, T_in_C=26.7)]),
        numerics=dict(dz_m=50.0, dt_initial_s=30.0, output_times_days=[0.1, 500.0]),
    )
    if constant:
        d["fluid"] = {"constant": {"rho": 1000.0, "cp": 4186.8, "mu": 3.0e-4, "k": 0.65}}
        d["options"] = dict(gravity_work=False)
    return d


def kubacka_case(depth=5000.0, pipe="hdpe", k_gap=0.06224, years=1.0):
    """Kubacka & Karayiannis (2024) Int. J. Low-Carbon Technol. 19:2068-2085 (WellTH).

    Vertical; 100 L/min water at 10 C down the annulus; surface 15 C, 30 C/km;
    rock 2750 kg/m3, 800 J/kgK, 2.5 W/mK; hole 215.9 mm; casing OD 198.52 mm,
    wall 10.36 mm (k 45), cement k 1.1 to the hole wall; centre pipe either HDPE
    (ID 88.3 / OD 139.7 mm, k 0.54) or VIT: steel 88.3/101.6 mm, gap to
    118.62 mm with the paper's equivalent gap conductivity at 40 C (Table 7),
    steel 118.62/139.7 mm.  1 year.
    """
    if pipe == "hdpe":
        layers = [dict(od_m=0.1397, k=0.54)]
    elif pipe == "vit_gap":     # alternative reading: k_eq applies to the 8.51 mm gap only
        layers = [dict(od_m=0.1016, k=45.0), dict(od_m=0.11862, k=k_gap), dict(od_m=0.1397, k=45.0)]
    else:                       # k_eq read as the apparent conductivity of the whole VIT wall (D1-D4)
        layers = [dict(od_m=0.1397, k=k_gap)]
    return dict(
        name=f"Kubacka {depth / 1000:.0f} km ({pipe})",
        well=dict(depth_m=depth,
                  casing=[dict(id_m=0.17780, layers=[dict(od_m=0.19852, material="steel"),
                                                     dict(od_m=0.2159, material="cement")])],
                  tubing=[dict(id_m=0.0883, layers=layers)]),
        materials=dict(cement=dict(k=1.1, rho=1900.0, cp=1000.0)),
        rock=dict(layers=[dict(k=2.5, rho=2750.0, cp=800.0)], surface_T_C=15.0,
                  geotherm=dict(gradient_C_per_km=30.0)),
        operation=dict(pressure_bc=dict(outlet_MPa=1.0),
                       schedule=[dict(duration_days=years * 365.0, mdot_kg_s=100.0 / 60.0 * 0.9997, T_in_C=10.0)]),
        numerics=dict(dz_m=50.0, dt_initial_s=300.0, output_times_days=[365.0]),
    )


def beckers_case(depth=4000.0, T_bottom=200.0, years=20.0):
    """Beckers et al. (2022) Geothermics 100:102318, vertical coaxial, SBT.

    Surface 20 C, rock k 2.83, rho 2875, c 825; well radius 0.1143 m; centre pipe
    ID radius 0.0635 m, wall 0.0127 m, vacuum insulated k 0.006 W/mK; 20 kg/s,
    inlet 20 C, 20 years; bottom-hole temperature 200 C.
    """
    return dict(
        name=f"Beckers 2022 {depth} m",
        well=dict(depth_m=depth, casing=[dict(id_m=0.2286, layers=[])],
                  tubing=[dict(id_m=0.127, layers=[dict(od_m=0.127 + 2 * 0.0127, k=0.006)])]),
        rock=dict(layers=[dict(k=2.83, rho=2875.0, cp=825.0)], surface_T_C=20.0,
                  geotherm=dict(gradient_C_per_km=(T_bottom - 20.0) / depth * 1000.0)),
        operation=dict(pressure_bc=dict(outlet_MPa=2.0),
                       schedule=[dict(duration_days=years * 365.0, mdot_kg_s=20.0, T_in_C=20.0)]),
        numerics=dict(dz_m=50.0, dt_initial_s=300.0, dt_max_days=60.0, output_times_days=[365.0 * years]),
    )


def b3_published():
    out = {}
    # --- Kubacka & Karayiannis 2024 (WellTH)
    rows = []
    pub = {  # depth: {case: (T_out, Q_kW)}
        1000.0: {"hdpe": (16.7, 46.7), 0.00024: (17.8, 53.9), 0.06224: (17.6, 52.8), 0.08091: (17.6, 52.6)},
        3000.0: {"hdpe": (31.6, 149.7), 0.00024: (54.6, 309.6), 0.06224: (49.6, 275.2), 0.08091: (48.4, 266.3)},
        5000.0: {"hdpe": (37.9, 193.5), 0.00024: (104.0, 655.0), 0.06224: (84.4, 517.4), 0.08091: (79.9, 486.3)},
    }
    labels = {"hdpe": "HDPE pipe", 0.00024: "VIT 0.01 Pa", 0.06224: "VIT 10 Pa", 0.08091: "VIT >=10 kPa"}
    for depth, cases in pub.items():
        for key, (refT, refQ) in cases.items():
            if key == "hdpe":
                sim, res = run(kubacka_case(depth, "hdpe"))
                Tg = Qg = None
            else:
                sim, res = run(kubacka_case(depth, "vit", key))
                simg, resg = run(kubacka_case(depth, "vit_gap", key))
                Tg = float(C(resg.at(365 * DAY, "T_out")))
                Qg = resg.at(365 * DAY, "Q_wellhead") / 1e3
            T = float(C(res.at(365 * DAY, "T_out")))
            Q = res.at(365 * DAY, "Q_wellhead") / 1e3
            rows.append(dict(depth=depth, case=labels[key], T_out=T, Q_kW=Q, published_T=refT, published_Q=refQ,
                             T_out_gap_reading=Tg, Q_kW_gap_reading=Qg))
            print(f"  Kubacka {depth / 1000:.0f} km {labels[key]:13s} model {T:6.1f} C {Q:6.1f} kW   "
                  f"WellTH {refT} C {refQ} kW", flush=True)
    out["kubacka"] = rows
    # --- Beckers et al. 2022 (SBT)
    rows = []
    for depth, ref in ((2000.0, 0.57), (3000.0, 0.85), (4000.0, 1.10)):
        sim, res = run(beckers_case(depth))
        Qavg = res.mean_over(0, 20 * 365 * DAY, "Q_wellhead") / 1e6
        rows.append(dict(depth=depth, Q_avg_MW=Qavg, published_MW=ref))
        print(f"  Beckers {depth:.0f} m: model 20-yr mean {Qavg:.3f} MWth   SBT {ref} MWth", flush=True)
    out["beckers"] = rows
    # --- Nalla et al. 2004 (TETRAD)
    rows = []
    variants = [("base, 6.3 kg/s", {}, 84.0), ("1.26 kg/s", dict(mdot=1.26), 123.0),
                ("31.5 kg/s", dict(mdot=31.5), 40.0), ("660 mm hole", dict(d_hole=0.6604), 95.0),
                ("insulation k 0.007", dict(k_ins=0.007), 86.0), ("insulation k 0", dict(k_ins=1e-5), 86.6),
                ("tubing OD 127 mm", dict(tub_od=0.127), 83.0),
                ("0.2 W/m2, 2796.5 m", dict(depth=2796.5), 59.0), ("0.5 W/m2, 1118.6 m", dict(depth=1118.6), 41.0)]
    for lab, kw, ref in variants:
        sim, res = run(nalla_case(**kw))
        T500 = float(C(res.at(500 * DAY, "T_out")))
        Tpk = float(C(np.max(res.ts["T_out"])))
        tpk = float(res.t_days[int(np.argmax(res.ts["T_out"]))])
        rows.append(dict(label=lab, T_out_500d=T500, published=ref, peak=Tpk, t_peak_days=tpk))
        print(f"  Nalla {lab:20s} model {T500:6.1f} C  TETRAD {ref} C   (peak {Tpk:.0f} C at {tpk:.2f} d)", flush=True)
    out["nalla"] = rows
    save("B3_published", out)


if __name__ == "__main__":
    which = sys.argv[1:] or ["b1", "b3"]
    for w in which:
        dict(b1=b1_geocluster, b2=b2_sbt, b3=b3_published)[w]()
