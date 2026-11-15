"""The benchmark examples that GEOPHIRES-X can represent, run in GEOPHIRES-X and set beside the aharrsim results.

Run the GEOPHIRES part with an interpreter that has GEOPHIRES-X installed (it is not on PyPI:
    git clone https://github.com/NREL/GEOPHIRES-X && pip install ./GEOPHIRES-X ):
    <geophires python> geophires_compare.py run        -> results/GX_runs.json
then, with any interpreter:
    python geophires_compare.py table                   -> results/GX_compare.json (and a printed table)

Each case mirrors the corresponding example of dev_types.py: same depth, gradient, flow, diameters, rock, injection temperature,
ambient (15 C), pump efficiency (0.7) and, for the fractured systems, the same Gringarten geometry.  GEOPHIRES is given the
reservoir temperature history of the convective and aquifer examples through its linear-drawdown model.  Plant: GEOPHIRES
subcritical ORC (the correlation aharrsim.plant.orc_practical reproduces) or double flash.
"""
from __future__ import annotations

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")

COMMON = {"Ambient Temperature": 15, "Circulation Pump Efficiency": 0.7, "Utilization Factor": 0.9, "Plant Lifetime": 30, "Time steps per year": 4,
          "Print Output to Console": 0, "Water Loss Fraction": 0, "Maximum Drawdown": 1, "Maximum Temperature": 500, "Number of Segments": 1,
          "Number of Production Wells": 1, "Number of Injection Wells": 1, "Injection Wellbore Temperature Gain": 0}


def conv(depth_km, grad, T_surf, m, d_in, T_inj, rock, drawdown=0.0, plant=1, end_use=1, **kw):
    """Convective / aquifer reservoir: linear drawdown of (T_res - T_inj) per year, Ramey production wellbore."""
    d = dict(COMMON)
    d.update({"Reservoir Model": 4, "Drawdown Parameter": drawdown, "Reservoir Depth": depth_km, "Gradient 1": grad, "Surface Temperature": T_surf,
              "Production Flow Rate per Well": m, "Production Well Diameter": d_in, "Injection Well Diameter": d_in, "Ramey Production Wellbore Model": 1,
              "Injection Temperature": T_inj, "Reservoir Volume Option": 4, "Reservoir Volume": 1e9, "Reservoir Density": rock[0], "Reservoir Heat Capacity": rock[1],
              "Reservoir Thermal Conductivity": rock[2], "End-Use Option": end_use, "Power Plant Type": plant})
    d.update(kw)
    return d


def mpf(depth_km, grad, m, T_inj, rock, N, height, width, sep, impedance, plant=1, d_in=6.094, **kw):
    """Gringarten parallel fractures (GEOPHIRES 'multiple parallel fractures'), rectangular fractures height x width (flow length)."""
    d = dict(COMMON)
    d.update({"Reservoir Model": 1, "Reservoir Volume Option": 1, "Fracture Shape": 4, "Fracture Height": height, "Fracture Width": width,
              "Number of Fractures": N, "Fracture Separation": sep, "Reservoir Depth": depth_km, "Gradient 1": grad, "Surface Temperature": 15,
              "Production Flow Rate per Well": m, "Production Well Diameter": d_in, "Injection Well Diameter": d_in, "Ramey Production Wellbore Model": 1,
              "Injection Temperature": T_inj, "Reservoir Impedance": impedance, "Reservoir Density": rock[0], "Reservoir Heat Capacity": rock[1],
              "Reservoir Thermal Conductivity": rock[2], "End-Use Option": 1, "Power Plant Type": plant})
    d.update(kw)
    return d


def sbt(L_vert, L_lat, n_lat, spacing, grad, m, T_inj, rock, d_lat=0.216, d_vert_in=8.535, **kw):
    """Closed U-loop with n_lat horizontal laterals between two verticals (GEOPHIRES slender-body-theory model)."""
    d = dict(COMMON)
    d.update({"Reservoir Model": 8, "Is AGS": True, "Well Geometry Configuration": 1, "Has Nonvertical Section": True, "Multilaterals Cased": False,
              "Number of Multilateral Sections": n_lat, "Vertical Section Length": L_vert, "Vertical Wellbore Spacing": L_lat, "Lateral Spacing": spacing,
              "Junction Depth": L_vert, "Lateral Endpoint Depth": L_vert, "Nonvertical Wellbore Diameter": d_lat, "Discretization Length": 250,
              "SBT Accuracy Desired": 5, "Gradient 1": grad, "Surface Temperature": 15, "Reservoir Depth": L_vert / 1000.0, "Reservoir Volume Option": 4,
              "Reservoir Volume": 1e9, "Reservoir Density": rock[0], "Reservoir Heat Capacity": rock[1], "Reservoir Thermal Conductivity": rock[2],
              "Production Flow Rate per Well": m, "Injection Temperature": T_inj, "Production Well Diameter": d_vert_in, "Injection Well Diameter": d_vert_in,
              "Ramey Production Wellbore Model": 0, "Reservoir Impedance": 1e-4, "Heat Transfer Fluid": 1, "End-Use Option": 1, "Power Plant Type": 1,
              "SBT Generate Wireframe Graphics": False})
    d.update(kw)
    return d


def cases(dt=None):
    """GEOPHIRES inputs per benchmark example id (as strings: '14b' is the Cape Station example)."""
    dd1 = 0.0
    if dt and 1.0 in dt:                      # type 1: linear drawdown that ends at aharrsim's year-30 reservoir temperature
        Tr = dt[1.0]["T_res"]; dd1 = (70.0 - Tr["30"]) / (70.0 - 40.0) / 30.0
    c = {
        "1": conv(1.8, (70.0 - 12.0) / 1.8, 12, 68, 8.535, 40, (2400, 900, 2.2), drawdown=dd1, end_use=2,
                  **{"End-Use Efficiency Factor": 1.0, "Productivity Index": 15.0, "Injectivity Index": 1000.0, "Reservoir Hydrostatic Pressure": 17500, "Production Wellhead Pressure": 500,
                     "Plant Outlet Pressure": 500}),
        "2": conv(2.0, 75.0, 15, 75, 8.535, 70, (2600, 900, 2.5), drawdown=0.5 / (165.0 - 70.0),
                  **{"Productivity Index": 4.0, "Injectivity Index": 1000.0, "Reservoir Hydrostatic Pressure": 18000, "Production Wellhead Pressure": 1200,
                     "Plant Outlet Pressure": 1200}),
        "3": conv(1.5, (250.0 - 15.0) / 1.5, 15, 40, 8.535, 70, (2500, 950, 2.2), plant=4,
                  **{"Productivity Index": 6.0, "Injectivity Index": 1000.0, "Reservoir Hydrostatic Pressure": 9000}),
        "15": conv(2.0, 40.0, 15, 15, 2.441, 60, (2500, 900, 2.3),
                   **{"Number of Production Wells": 8, "Number of Injection Wells": 8, "Productivity Index": 1000.0, "Injectivity Index": 1000.0,
                      "Production Wellhead Pressure": 500, "Plant Outlet Pressure": 500}),
        "10": mpf(5.0, (200.0 - 15.0) / 5.0, 30, 60, (2650, 850, 3.0), 20, 400, 600, 50, 0.3),
        "14": mpf(2.3, (190.0 - 15.0) / 2.3, 57.96, 60, (2650, 850, 2.8), 88, 143, 110, 975.0 / 88, 0.15,
                  **{"Well Geometry Configuration": 4, "Multilaterals Cased": True, "Number of Multilateral Sections": 2,
                     "Nonvertical Length per Multilateral Section": 975}),
        "14b": mpf(2.59, (200.0 - 15.0) / 2.59, 100, 60, (2650, 850, 2.8), 300, 143, 137, 1524.0 / 150, 0.0894,
                   # one injector here: with two, GEOPHIRES' impedance mode charges the pumps for twice the circulating flow
                   **{"Well Geometry Configuration": 4, "Multilaterals Cased": True, "Number of Multilateral Sections": 2,
                      "Nonvertical Length per Multilateral Section": 1524}),
        "19": mpf(3.5, (375.0 - 15.0) / 3.5, 40, 60, (2700, 950, 2.5), 20, 400, 300, 40, 0.3, plant=4),
        "8": sbt(4500, 5000, 12, 80, 55.0, 100, 40, (2600, 900, 2.8)),
        "9": sbt(3000, 3000, 2, 120, 40.0, 35, 40, (2600, 900, 2.6)),
    }
    return c


NOT_REPRESENTED = {
    "3b": "dry-steam producer: no steam wellbore or steam-turbine model",
    "4": "surface stream only: the plant correlation is the same one, so there is nothing to compare",
    "5": "ground-source heat pump on a U-tube: no shallow borehole-exchanger model",
    "6": "shallow coaxial exchanger with heat pump: as 5",
    "7": "coaxial loop: the slender-body model has no coaxial option in this version and the tabulated closed-loop database stops at 70 C/km and 0.216 m bore",
    "11": "single well with a fractured path at the toe of a coaxial completion: no such element",
    "12": "cyclic single fracture: no huff-and-puff element",
    "13": "as 12",
    "16": "coaxial loop in 7 in casing: outside the tabulated database (bore below 0.216 m); no coaxial slender-body option",
    "17": "heat store after thermal recovery: no such element",
    "18": "supercritical feed: no supercritical or steam wellbore",
    "20": "20 km coaxial loop: no coaxial slender-body option; database stops at 5 km",
}


def run():
    from geophires_x_client import GeophiresXClient, GeophiresInputParameters
    import geophires_x
    try:
        dt = {float(r["id"]): r for r in json.load(open(os.path.join(RES, "DT_types.json")))}
    except Exception:
        dt = None
    client = GeophiresXClient()
    out = {"geophires_version": getattr(geophires_x, "__version__", "?"), "cases": {}}
    for cid, params in cases(dt).items():
        try:
            r = client.get_geophires_result(GeophiresInputParameters({k: str(v) for k, v in params.items()})).result
            prof = r.get("POWER GENERATION PROFILE") or r.get("HEAT AND/OR ELECTRICITY EXTRACTION AND GENERATION PROFILE")
            keep = {}
            for sec in ("SUMMARY OF RESULTS", "ENGINEERING PARAMETERS", "RESERVOIR PARAMETERS", "RESERVOIR SIMULATION RESULTS", "SURFACE EQUIPMENT SIMULATION RESULTS"):
                keep[sec] = {k: (v.get("value") if isinstance(v, dict) else v) for k, v in (r.get(sec) or {}).items() if v is not None}
            out["cases"][cid] = dict(inputs=params, profile_header=prof[0] if prof else None, profile=prof[1:] if prof else None,
                                     profiles={k: r[k] for k in r if "PROFILE" in k and k != "SAM CASH FLOW PROFILE" and r[k]}, **keep)
            row1, rowN = prof[1], prof[-1]
            print(cid, "ok", prof[0], row1, rowN, flush=True)
        except Exception as e:
            out["cases"][cid] = dict(inputs=params, failed=f"{type(e).__name__}: {str(e)[:400]}")
            print(cid, "FAILED", str(e)[:300], flush=True)
    json.dump(out, open(os.path.join(RES, "GX_runs.json"), "w"))


NAMES = {"1": "Aquifer doublet (heat)", "2": "Pumped well, binary", "3": "Flash well", "15": "Co-produced water, 8 wells", "10": "EGS doublet, 200 C",
         "14": "Project Red type", "14b": "Cape Station type", "19": "Superhot EGS doublet", "8": "Multilateral loop", "9": "Closed-loop doublet"}
DT = {"14b": 14.5}


def table():
    import numpy as np
    G = json.load(open(os.path.join(RES, "GX_runs.json")))
    R = {float(r["id"]): r for r in json.load(open(os.path.join(RES, "DT_types.json")))}
    S = {float(r["id"]): r for r in json.load(open(os.path.join(RES, "DT_summary.json")))}
    rows = []
    for cid in ("1", "2", "3", "15", "10", "14", "14b", "19", "8", "9"):
        g = G["cases"][cid]
        if "failed" in g:
            rows.append(dict(id=cid, name=NAMES[cid], failed=g["failed"])); continue
        tid = DT.get(cid) or float(cid); r = R[tid]; ser = r["series"]; sm = S[tid]
        hdr, P = g["profile_header"], g["profile"]
        iT = [i for i, h in enumerate(hdr) if "TEMPERATURE" in h][0]
        iP = [i for i, h in enumerate(hdr) if "PUMP" in h][0]
        iN = [i for i, h in enumerate(hdr) if "NET" in h][0]
        heat_case = "NET HEAT" in hdr[iN]
        wp = r.get("well_pump_kW") or {}
        ts = (2, 29) if cid in ("8", "9") else (1, 29)      # GEOPHIRES' closed-loop profile gives an early-transient value in its year-1 row
        rsr = g["RESERVOIR SIMULATION RESULTS"]
        row = dict(id=cid, name=NAMES[cid], heat_case=heat_case, t=list(ts), coax_T=[], gx_T=[], coax_out=[], gx_out=[], coax_pump=[], gx_pump=[],
                   gx_dp_MPa={k.replace("Average ", "").replace(" Pressure Drop", "").lower(): v / 1e3 for k, v in rsr.items() if "Pressure Drop" in k and isinstance(v, (int, float))},
                   gx_dT_prod=g["ENGINEERING PARAMETERS"].get("Average production well temperature drop"), gx_T_inj=g["ENGINEERING PARAMETERS"].get("Injection temperature"))
        for t in ts:
            x = P[t]
            row["coax_T"].append(float(np.interp(t, ser["t_yr"], ser["T_C"])))
            row["gx_T"].append(float(x[iT]))
            if heat_case:
                row["coax_out"].append(float(np.interp(t, ser["t_yr"], ser["heat_MWth"]))); row["gx_out"].append(float(x[iN]))
            else:
                if cid == "19":       # steam / flash on the two-phase stream: values of the summary (years 1 and 30)
                    row["coax_out"].append(sm["gross_first_MWe"] if t <= 2 else sm["gross_last_MWe"])
                else:
                    row["coax_out"].append(float(np.interp(t, ser["t_yr"], ser["MWe"])))
                row["gx_out"].append(float(x[iN]) + float(x[iP]))               # gross = net + pumps
            ys = sorted(wp, key=float)
            row["coax_pump"].append(float(np.interp(t, [float(y) for y in ys], [max(wp[y], 0.0) for y in ys])) / 1e3 if ys else 0.0)
            row["gx_pump"].append(float(x[iP]))
        # whole history (GEOPHIRES row i is i years from the start), and for the cases that break through the time at which the produced
        # water has lost half of its rise above the injection temperature
        tg = np.arange(len(P), dtype=float); Tg = np.array([float(x[iT]) for x in P])
        tc = np.array(ser["t_yr"], float); Tc = np.array(ser["T_C"], float)
        k0 = 2 if cid in ("8", "9") else 1
        row["hist"] = dict(t=tg[k0:].tolist(), gx_T=Tg[k0:].tolist(), coax_T=np.interp(tg[k0:], tc, Tc).tolist())
        if cid in ("14", "14b"):
            T_inj = float(row["gx_T_inj"] or 60.0)

            def t_half(t, T):
                mid = 0.5 * (T[np.searchsorted(t, 1.0)] + T_inj)
                i = int(np.argmax(T < mid))
                return float(t[i - 1] + (T[i - 1] - mid) / (T[i - 1] - T[i]) * (t[i] - t[i - 1])) if i > 0 else None

            tf = np.linspace(0.5, 29.0, 571)
            row["t_half"] = dict(coax=t_half(tf, np.interp(tf, tc, Tc)), gx=t_half(tg, Tg))
            d = np.interp(tg[k0:], tc, Tc) - Tg[k0:]
            row["hist_max_diff_K"] = float(np.max(np.abs(d)))
        rows.append(row)
    out = dict(geophires_version=G["geophires_version"], rows=rows, not_represented=NOT_REPRESENTED)
    json.dump(out, open(os.path.join(RES, "GX_compare.json"), "w"), indent=1)
    print(f"GEOPHIRES-X {G['geophires_version']}:  T at year 1 (2 for the closed loops) and 29 (aharrsim | GEOPHIRES), output MW (gross electricity or heat), pumps MW")
    for r in rows:
        if "failed" in r:
            print(f"{r['id']:>4} {r['name']:<28} FAILED {r['failed'][:80]}"); continue
        f = lambda v: "/".join(f"{x:6.2f}" for x in v)
        print(f"{r['id']:>4} {r['name']:<28} T {f(r['coax_T'])} | {f(r['gx_T'])}   out {f(r['coax_out'])} | {f(r['gx_out'])}   pump {f(r['coax_pump'])} | {f(r['gx_pump'])}   gx dT_prod {r['gx_dT_prod']}")
    return out


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "table"
    dict(run=run, table=table)[what]()
