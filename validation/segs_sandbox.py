"""Single-well enhanced geothermal system (SEGS): the laboratory sandbox of Ma et al. (2026) set against aharrsim.

Ma, B., Li, A., Chen, W., Gao, X., Han, W., Lei, Y., Wang, L., 2026. Experimental investigation on heat transfer characteristics
of single-well enhanced geothermal system. Geothermal Energy 14, 9. doi:10.1186/s40517-026-00379-4

The experiment: a closed, insulated steel cylinder 1.5 m high and 1 m across, packed with 0.8-1.2 mm glass beads and filled with
hot water (50-70 C). One well on the axis: cold water (21-41 C) leaves through an upper screen (10 cm long, 16 mm across), sweeps
the bead pack and returns through a lower screen at the base, 90 or 120 cm away; 0.2-0.6 m3/h for one to two hours. It is the
development type 11 of the benchmark (single-well EGS) at 1 : 150 or so, and a case of swept rock in the limit of local thermal
equilibrium.

    python segs_sandbox.py checks    the numerical block model against a piston front and the unbounded source-sink pair
    python segs_sandbox.py lab       the experiment: flow, temperature and spacing series, sensitivities   -> results/SG_lab.json
    python segs_sandbox.py calib     the one property not measured (vertical permeability) fitted to one number -> results/SG_calib.json
    python segs_sandbox.py field     what it means at field scale; the type 11 example with a dipole sweep -> results/SG_field.json
    python segs_sandbox.py           all of them

The numbers reported by the paper are in REPORTED, with the reading made of each where the text is ambiguous.
"""
from __future__ import annotations

import math
import sys
import time

import numpy as np

from common import save, load
from aharrsim.porous import AxisymmetricBlock, water_table
from aharrsim.reservoir import SingleWellDipole, ParallelFractures, dipole_arrival, cbar_w

MIN = 60.0
YEAR = 365.25 * 86400.0
CRHO_PAPER = 4.2e6        # the paper's heat flow is C rho q (T_out - T_in) with C rho = 4.2e6 J/m3/K (its averages reproduce to 0.1 %)

# ------------------------------------------------------------------------------------------------ the experiment as published
LAB = dict(R=0.5, H=1.5, r_w=0.008, screen=0.10, phi=0.35, rho_bead=1500.0, c_bead=800.0, lam_bead=2.5, K=8.92e-10, d_bead=1.0e-3,
           T0=60.0, T_in=21.0, q=0.5, spacing=0.90)

REPORTED = dict(
    source="Ma et al. (2026), Geothermal Energy 14:9; text of Sections 3-4 (the table of operating conditions was not accessible)",
    plateau_min={"0.3": 50.0, "0.5": 30.0, "0.6": 25.0},          # 0.5: typical case and the injection-temperature series
    typical=dict(q=0.5, T0=60.0, T_in=21.0, plateau_min=30.0, decline_C_per_min=0.8, T_after_30min_of_decline=37.0,
                 heat_first_30min_kW=22.0, heat_at_70min_kW=10.0),
    decay60_percent_flow_series=[8.3, 17.3, 23.3, 29.7, 32.2],     # five values printed for 0.3, 0.4, 0.5, 0.6 m3/h
    avg_flow_series={"0.2": dict(T=59.0, kW=8.9), "0.6": dict(T=52.8, kW=22.3)},
    T0_series={"50": dict(T_avg=47.72, kW=15.59, decay60=20.0), "60": dict(T_avg=55.29, kW=20.00, decay60=29.67),
               "70": dict(T_avg=63.44, kW=24.76, decay60=33.43)},          # 0.5 m3/h, 21 C injection, averages over 60 min
    Tin_series=dict(dT_avg_21_to_41=2.3, d_kW_21_to_41=-10.33),            # 0.5 m3/h, 60 C
    spacing=dict(decay_90=31.5, decay_120=8.0, dT_avg=4.4, d_kW=2.6),
    fronts=dict(H1_first_depth_cm=50.0, H1_first_time_min=30.0, down_cm_per_min=1.0, stable_above_depth_cm=30.0, stable_below_depth_cm=110.0,
                cooling_50_80cm_C_per_min=0.42, L1_radial_cm_per_min=0.75, L1_20K_at_40min_fraction=0.8),
    model_of_the_paper=dict(max_deviation_percent=8.0, overall_percent=4.0, note="experiment declines slightly later; rates converge after 40 min"),
    notes=["decay is taken as (T0 - T at 60 min) / T0 with temperatures in C: the three T0-series values then give 40.0, 42.2 and 46.6 C",
           "the flow series lists five decay values for four flows; 0.5 m3/h is 23.3 % there and 29.67 % in the T0 series, and the typical case"
           " reaches 37 C (38 %): the 60-minute temperature of the nominal case is known to within 36-46 C",
           "line H1 is 10 cm from the wall in the methods and 30 cm from the axis in the results; line L1 is 55 cm from the base",
           "the same three decay values are printed for the T0 series and for the injection-temperature series"])

_W = {}


def water(p=2.0e5):
    if p not in _W:
        _W[p] = water_table(p)
    return _W[p]


def lam_pack(lam_s, lam_w=0.64, phi=0.35):
    """Conductivity of the water-filled bead pack: geometric mean of the two phases."""
    return lam_s ** (1.0 - phi) * lam_w ** phi


def sandbox(q=None, T0=None, T_in=None, spacing=None, centre=True, bulk_density=True, rhoc_extra=0.0, K=None, aniso=1.0, buoyancy=True,
            viscosity=True, lam=None, disp=None, inject="top", nr=28, nz=120, t_end=75 * MIN, lines=True, snapshots=(), flow_update=10.0):
    """One run of the sandbox. spacing: between the screen centres (centre=True) or between their near ends.
    bulk_density: the 1500 kg/m3 of the paper is the density of the bead pack (True) or of the glass (False).
    aniso = K_horizontal / K_vertical with the geometric mean kept. Returns the produced temperature and the measures of the paper."""
    L = LAB
    q = L["q"] if q is None else q; T0 = L["T0"] if T0 is None else T0; T_in = L["T_in"] if T_in is None else T_in
    spacing = L["spacing"] if spacing is None else spacing
    K = L["K"] if K is None else K
    rhoc_pack = L["rho_bead"] * L["c_bead"] * (1.0 if bulk_density else (1.0 - L["phi"]))     # J/K per m3 of pack, solid part
    rhoc_s = rhoc_pack / (1.0 - L["phi"]) + rhoc_extra / (1.0 - L["phi"])
    lam = lam_pack(L["lam_bead"]) if lam is None else lam
    disp = L["d_bead"] if disp is None else disp
    b = AxisymmetricBlock(R=L["R"], H=L["H"], r_w=L["r_w"], K=K * math.sqrt(aniso), K_z=K / math.sqrt(aniso), phi=L["phi"], rhoc_s=rhoc_s, lam=lam,
                          T0=T0, nr=nr, nz=nz, disp=disp, buoyancy=buoyancy, viscosity=viscosity, T_ref=0.5 * (T0 + T_in), water=water())
    Q = q / 3600.0
    z_low = (0.0, L["screen"])
    zc = 0.5 * L["screen"] + spacing if centre else 1.5 * L["screen"] + spacing
    z_up = (zc - 0.5 * L["screen"], zc + 0.5 * L["screen"])
    if inject == "top":
        b.add_screen(*z_low, -Q); b.add_screen(*z_up, +Q, T_in=T_in)
    else:
        b.add_screen(*z_up, -Q); b.add_screen(*z_low, +Q, T_in=T_in)
    mon = dict(H30=dict(r=0.30), H40=dict(r=0.40), L55=dict(z=0.55), Linj=dict(z=zc)) if lines else None
    t0 = time.time()
    r = b.run(t_end, dt_out=30.0, flow_update=flow_update, cfl=0.6, lines=mon, snapshots=snapshots)
    t, T = r["t"], r["T_prod"]
    V_cyl = math.pi * L["R"] ** 2 * (zc - 0.5 * L["screen"])                 # the cylinder between the screen centres
    V_box = math.pi * L["R"] ** 2 * L["H"]
    out = dict(q=q, T0=T0, T_in=T_in, spacing=spacing, z_inj=list(z_up if inject == "top" else z_low), z_prod=list(z_low if inject == "top" else z_up),
               rhoc_b=b.rhoc_b, rhoc_w=b.rhoc_w, lam=lam, t_min=(t / MIN).tolist(), T_prod=T.tolist(), wall_s=time.time() - t0,
               balance=float(np.max(np.abs(r["E_out"] + r["E_block"])) / max(np.max(np.abs(r["E_out"])), 1e-30)),
               t_s_cyl_min=V_cyl * b.rhoc_b / (b.rhoc_w * Q) / MIN, t_s_box_min=V_box * b.rhoc_b / (b.rhoc_w * Q) / MIN,
               V_cyl=V_cyl, V_box=V_box)
    out.update(measures(t / MIN, T, T0, T_in, q))
    out["x_onset_cyl"] = out["onset_1K_min"] / out["t_s_cyl_min"] if out["onset_1K_min"] else None
    out["x_onset_box"] = out["onset_1K_min"] / out["t_s_box_min"] if out["onset_1K_min"] else None
    if lines:
        out["fronts"] = fronts(r["lines"], T0, L["H"])
    if snapshots:
        out["snapshots"] = {str(int(s / MIN)): np.round(r["snapshots"][s], 2).tolist() for s in r["snapshots"]}
        out["grid"] = dict(rc=b.rc.tolist(), zc=b.zc.tolist(), rf=b.rf.tolist(), zf=b.zf.tolist())
    return out


def measures(tm, T, T0, T_in, q):
    """The measures the paper reports, from a produced-temperature history (tm in minutes)."""
    def onset(d):
        k = np.where(T <= T0 - d)[0]
        if len(k) == 0:
            return None
        k = k[0]
        return float(tm[k]) if k == 0 else float(np.interp(T0 - d, [T[k], T[k - 1]], [tm[k], tm[k - 1]]))
    at = lambda m: float(np.interp(m, tm, T)) if tm[-1] >= m - 1e-9 else None
    out = dict(onset_05K_min=onset(0.5), onset_1K_min=onset(1.0), onset_2K_min=onset(2.0), T30=at(30), T40=at(40), T50=at(50), T60=at(60), T70=at(70))
    if tm[-1] >= 60 - 1e-9:
        m = tm <= 60 + 1e-9
        avg = float(np.sum(0.5 * (T[m][1:] + T[m][:-1]) * np.diff(tm[m])) / 60.0)
        out.update(T_avg60=avg, kW_avg60=CRHO_PAPER * q / 3600.0 * (avg - T_in) / 1e3, decay60=100.0 * (T0 - out["T60"]) / T0,
                   theta60=(out["T60"] - T_in) / (T0 - T_in), theta_avg60=(avg - T_in) / (T0 - T_in),
                   kW_start=CRHO_PAPER * q / 3600.0 * (T0 - T_in) / 1e3)
        # steepest decline over ten minutes
        g = [(at(a) - at(a + 10)) / 10.0 for a in np.arange(0, 50.1, 1.0)]
        out["max_decline_C_per_min"] = float(max(g))
    return out


def fronts(lines, T0, H):
    """Front of the 5 K drop along the vertical lines (depth below the top, cm) and of the 20 K drop along the radial lines (cm from the axis)."""
    out = {}
    for name in ("H30", "H40"):
        m = lines[name]; z = np.array(m["x"]); rows = {}
        for t, Tl in zip(m["t"], m["T"]):
            tmn = t / MIN
            if abs(tmn - round(tmn)) > 1e-6 or int(round(tmn)) % 10:
                continue
            cold = np.where(np.array(Tl) <= T0 - 5.0)[0]
            rows[str(int(round(tmn)))] = None if len(cold) == 0 else dict(top_cm=100.0 * (H - z[cold[-1]]), bottom_cm=100.0 * (H - z[cold[0]]))
        out[name] = rows
        d = [(int(k), v["bottom_cm"]) for k, v in rows.items() if v and 30 <= int(k) <= 60]
        out[name + "_down_cm_per_min"] = (d[-1][1] - d[0][1]) / (d[-1][0] - d[0][0]) if len(d) >= 2 and d[-1][0] > d[0][0] else None
    for name in ("L55", "Linj"):
        m = lines[name]; r = np.array(m["x"]); rows = {}
        for t, Tl in zip(m["t"], m["T"]):
            tmn = t / MIN
            if abs(tmn - round(tmn)) > 1e-6 or int(round(tmn)) % 10:
                continue
            cold = np.where(np.array(Tl) <= T0 - 20.0)[0]
            rows[str(int(round(tmn)))] = None if len(cold) == 0 else 100.0 * float(r[cold[-1]])
        out[name] = rows
        d = [(int(k), v) for k, v in rows.items() if v is not None and 10 <= int(k) <= 40]
        out[name + "_cm_per_min"] = (d[-1][1] - d[0][1]) / (d[-1][0] - d[0][0]) if len(d) >= 2 and d[-1][0] > d[0][0] else None
    return out


def brief(r):
    f = lambda v, d=1: "–" if v is None else f"{v:.{d}f}"
    return (f"q {r['q']:.1f} T0 {r['T0']:.0f} Tin {r['T_in']:.0f} D {r['spacing']:.2f}: onset(1 K) {f(r['onset_1K_min'])} min, T60 {f(r.get('T60'))} C, "
            f"decay {f(r.get('decay60'))} %, avg {f(r.get('T_avg60'), 2)} C, {f(r.get('kW_avg60'), 2)} kW, x_onset {f(r['x_onset_cyl'], 2)} (cyl) "
            f"{f(r['x_onset_box'], 2)} (box), balance {r['balance']:.0e}, {r['wall_s']:.0f} s")


# ------------------------------------------------------------------------------------------------ checks of the block model
def checks():
    out = {}
    # (1) piston: a thin shell (6 mm), flow from a one-cell screen at the top to one at the base; no conduction, so only the scheme
    #     and the residence time in the two screen cells (an exponential delay of one cell-transit time each) spread the front
    rw = 0.494
    piston = {}
    for nz in (200, 400):
        b = AxisymmetricBlock(R=0.50, H=1.0, r_w=rw, K=1e-9, phi=0.35, rhoc_s=1.2e6 / 0.65, lam=0.0, T0=60.0, nr=3, nz=nz, disp=0.0,
                              buoyancy=False, viscosity=False, T_ref=40.0, water=water())
        Q = 1.0e-5
        b.add_screen(0.0, 1.0 / nz, -Q); b.add_screen(1.0 - 1.0 / nz, 1.0, +Q, T_in=20.0)
        t_s = math.pi * (0.5 ** 2 - rw ** 2) * 1.0 * b.rhoc_b / (b.rhoc_w * Q)
        r = b.run(1.3 * t_s, dt_out=t_s / 500.0)
        th = (r["T_prod"] - 20.0) / 40.0
        t50, t90, t10 = (float(np.interp(v, th[::-1], r["t"][::-1])) for v in (0.5, 0.9, 0.1))
        piston[str(nz)] = dict(front_over_sweep_time=t50 / t_s, width_10_90_over_sweep_time=(t10 - t90) / t_s, width_cells=(t10 - t90) / t_s * nz,
                               balance=float(np.max(np.abs(r["E_out"] + r["E_block"])) / np.max(np.abs(r["E_out"]))))
        print("piston", nz, {k: round(v, 4) if v > 1e-6 else v for k, v in piston[str(nz)].items()})
    out["piston"] = piston
    # (2) the streamline table of reservoir.SingleWellDipole, integrated again
    from scipy.integrate import solve_ivp

    def vel(t, x):
        r_, z = x
        r1 = math.hypot(r_, z - 0.5); r2 = math.hypot(r_, z + 0.5)
        return [(r_ / r1 ** 3 - r_ / r2 ** 3) / (4 * math.pi), ((z - 0.5) / r1 ** 3 - (z + 0.5) / r2 ** 3) / (4 * math.pi)]

    def tau(F, eps=5e-4):
        a = math.acos(1 - 2 * F)
        hit = lambda t, x: math.hypot(x[0], x[1] + 0.5) - eps
        hit.terminal = True; hit.direction = -1
        s = solve_ivp(vel, [0, 1e8], [eps * math.sin(a), 0.5 - eps * math.cos(a)], events=hit, rtol=1e-10, atol=1e-13, method="DOP853")
        return s.t_events[0][0] + 2 * (4 * math.pi / 3) * eps ** 3
    rows = [dict(F=F, integrated=tau(F), table=float(dipole_arrival(F))) for F in (0.02, 0.12, 0.33, 0.52, 0.72)]
    out["streamlines"] = dict(rows=rows, max_rel=max(abs(r_["integrated"] / r_["table"] - 1) for r_ in rows))
    print("streamline table: max relative difference", out["streamlines"]["max_rel"])
    # (3) the source-sink pair in a block much larger than the spacing, against the element (D = 1 m, short screens)
    D = 1.0
    b = AxisymmetricBlock(R=3.0, H=6.0, r_w=0.01, K=1e-9, phi=0.35, rhoc_s=1.2e6 / 0.65, lam=0.0, T0=60.0, nr=60, nz=240, r1=0.035, disp=0.0,
                          buoyancy=False, viscosity=False, T_ref=40.0, water=water())
    Q = 1.0 / 3600.0
    b.add_screen(2.5 - 0.025, 2.5 + 0.025, -Q); b.add_screen(3.5 - 0.025, 3.5 + 0.025, +Q, T_in=20.0)
    el = SingleWellDipole(Q=Q, D=D, rhoc_b=b.rhoc_b, T_rock_C=60.0, T_inj_C=20.0, rhoc_w=b.rhoc_w)
    r = b.run(4.0 * el.t_ref, dt_out=el.t_ref / 100.0)
    x = r["t"] / el.t_ref
    TD_num = (60.0 - r["T_prod"]) / 40.0
    TD_el = el.T_D(r["t"])
    xs = [0.6, 0.746, 0.9, 1.2, 1.5, 2.0, 3.0, 4.0]
    out["dipole"] = dict(x=xs, block=[float(np.interp(v, x, TD_num)) for v in xs], element=[float(np.interp(v, x, TD_el)) for v in xs],
                         max_diff=float(np.max(np.abs(TD_num - TD_el))), R_over_D=3.0, H_over_D=6.0, screen_over_D=0.05,
                         series=dict(x=x[::4].tolist(), block=TD_num[::4].tolist(), element=TD_el[::4].tolist()),
                         balance=float(np.max(np.abs(r["E_out"] + r["E_block"])) / np.max(np.abs(r["E_out"]))))
    print("dipole in a large block: T_D block", np.round(out["dipole"]["block"], 3), "element", np.round(out["dipole"]["element"], 3), "max diff", round(out["dipole"]["max_diff"], 4))
    # (4) grid refinement of the sandbox itself (flow not coupled to temperature, so that only the grid changes)
    ref = {}
    for nr, nz in ((20, 60), (28, 120), (40, 180), (56, 240)):
        r_ = sandbox(buoyancy=False, viscosity=False, nr=nr, nz=nz, lines=False)
        ref[f"{nr}x{nz}"] = dict(onset_1K_min=r_["onset_1K_min"], T40=r_["T40"], T60=r_["T60"], T_avg60=r_["T_avg60"], wall_s=r_["wall_s"], T=r_["T_prod"], t=r_["t_min"])
        print("grid", nr, nz, brief(r_))
    fine = np.array(ref["56x240"]["T"])
    for k, v in ref.items():
        v["max_diff_from_finest_K"] = float(np.max(np.abs(np.array(v["T"]) - fine)))
        del v["T"], v["t"]
    out["grid"] = ref
    save("SG_checks", out)
    return out


# ------------------------------------------------------------------------------------------------ the experiment
def lab():
    out = dict(reported=REPORTED, lab=LAB)
    # typical case, two hours, with the temperature field at the times the paper discusses
    typ = sandbox(t_end=120 * MIN, snapshots=(20 * MIN, 30 * MIN, 40 * MIN, 60 * MIN))
    print("typical       ", brief(typ))
    out["typical"] = typ
    lin = sandbox(buoyancy=False, viscosity=False, t_end=120 * MIN)
    print("typical linear", brief(lin))
    out["typical_linear"] = lin
    # flow series (60 C, 21 C)
    out["flow"] = {}
    for q in (0.2, 0.3, 0.4, 0.5, 0.6):
        r = typ if q == 0.5 else sandbox(q=q, t_end=(150 if q < 0.35 else 75) * MIN)
        out["flow"][f"{q:.1f}"] = {k: v for k, v in r.items() if k not in ("snapshots", "grid")}
        print("flow          ", brief(r))
    # initial temperature of the pack (0.5 m3/h, 21 C)
    out["T0"] = {}
    for T0 in (50.0, 60.0, 70.0):
        r = typ if T0 == 60.0 else sandbox(T0=T0)
        out["T0"][f"{T0:.0f}"] = {k: v for k, v in r.items() if k not in ("snapshots", "grid")}
        print("T0            ", brief(r))
    # injection temperature (0.5 m3/h, 60 C)
    out["T_in"] = {}
    for Ti in (21.0, 31.0, 41.0):
        r = typ if Ti == 21.0 else sandbox(T_in=Ti)
        out["T_in"][f"{Ti:.0f}"] = {k: v for k, v in r.items() if k not in ("snapshots", "grid")}
        print("T_in          ", brief(r))
    # spacing
    out["spacing"] = {}
    for D in (0.90, 1.20):
        r = typ if D == 0.90 else sandbox(spacing=D, t_end=120 * MIN)
        out["spacing"][f"{D:.2f}"] = {k: v for k, v in r.items() if k not in ("snapshots", "grid")}
        print("spacing       ", brief(r))
    # what the model is sensitive to, on the typical case
    sens = {}
    for name, kw in (("as published", {}),
                     ("flow not coupled to temperature", dict(buoyancy=False, viscosity=False)),
                     ("buoyancy only", dict(viscosity=False)),
                     ("viscosity only", dict(buoyancy=False)),
                     ("permeability 5.6e-10 m2 (Kozeny-Carman recomputed)", dict(K=5.6e-10)),
                     ("1500 kg/m3 read as the density of the glass", dict(bulk_density=False)),
                     ("glass at 2500 kg/m3", dict(rhoc_extra=(2500.0 * 0.65 - 1500.0) * 800.0)),
                     ("spacing read between the near ends of the screens", dict(centre=False)),
                     ("conductivity of the pack 1.0 W/m/K", dict(lam=1.0)),
                     ("no dispersion", dict(disp=0.0)),
                     ("dispersivity 5 mm", dict(disp=5e-3)),
                     ("horizontal permeability twice the vertical", dict(aniso=2.0)),
                     ("horizontal permeability four times the vertical", dict(aniso=4.0)),
                     ("injection through the lower screen", dict(inject="bottom"))):
        r = typ if not kw else sandbox(lines=False, **kw)
        sens[name] = {k: r.get(k) for k in ("onset_05K_min", "onset_1K_min", "onset_2K_min", "T30", "T40", "T50", "T60", "T70", "T_avg60", "kW_avg60", "decay60",
                                            "theta60", "theta_avg60", "x_onset_cyl", "x_onset_box", "rhoc_b", "max_decline_C_per_min", "balance")}
        sens[name]["t_min"] = r["t_min"]; sens[name]["T_prod"] = r["T_prod"]
        print(f"{name:52s}", brief(r))
    out["sensitivity"] = sens
    # the plateau as a flushed share: reported plateaus against flow
    rep = REPORTED["plateau_min"]
    out["plateau_scaling"] = {q: dict(reported_min=m, thermal_volume_m3=float(q) / 3600.0 * m * MIN * typ["rhoc_w"] / typ["rhoc_b"],
                                      x_cyl=m / out["flow"][f"{float(q):.1f}"]["t_s_cyl_min"], x_box=m / out["flow"][f"{float(q):.1f}"]["t_s_box_min"],
                                      model_onset_1K_min=out["flow"][f"{float(q):.1f}"]["onset_1K_min"]) for q, m in rep.items()}
    for q, v in out["plateau_scaling"].items():
        print("plateau", q, {k: round(x, 3) for k, x in v.items()})
    save("SG_lab", out)
    return out


# ------------------------------------------------------------------------------------------------ one free parameter
KEYS = ("onset_05K_min", "onset_1K_min", "onset_2K_min", "T30", "T40", "T50", "T60", "T70", "T_avg60", "kW_avg60", "decay60", "theta60", "theta_avg60",
        "x_onset_cyl", "x_onset_box", "rhoc_b", "max_decline_C_per_min", "balance", "t_s_cyl_min")


def pick(r, series=True):
    o = {k: r.get(k) for k in KEYS}
    if series:
        o["t_min"] = r["t_min"]; o["T_prod"] = r["T_prod"]
    return o


def calibrate():
    """With the published properties the model declines early, because the cold water sinks under its own weight faster than the
    experiment shows. The vertical permeability is the one property that controls that and was not measured (the published value is a
    Kozeny-Carman estimate). Two ways of lowering it are fitted to ONE number, the 60-minute mean production temperature of the nominal
    case (55.29 C), and then run through every other case of the paper without further adjustment:
      (a) an isotropic permeability below the published one;   (b) the published mean with the horizontal above the vertical."""
    from scipy.optimize import brentq
    target = REPORTED["T0_series"]["60"]["T_avg"]
    out = dict(target_T_avg60=target, scan={})
    for name, kw in (("K 8.92e-10 m2 (published)", {}), ("K 5.6e-10 m2", dict(K=5.6e-10)), ("K 3e-10 m2", dict(K=3e-10)), ("K 2e-10 m2", dict(K=2e-10)),
                     ("K 1e-10 m2", dict(K=1e-10)), ("K 3e-11 m2", dict(K=3e-11)), ("Kh/Kv 1.5", dict(aniso=1.5)), ("Kh/Kv 2", dict(aniso=2.0)),
                     ("Kh/Kv 3", dict(aniso=3.0)), ("Kh/Kv 4", dict(aniso=4.0))):
        r = sandbox(lines=False, **kw)
        out["scan"][name] = pick(r, series=False)
        print(f"scan {name:28s}", brief(r))
    fK = lambda lg: sandbox(lines=False, K=10.0 ** lg)["T_avg60"] - target
    fA = lambda a: sandbox(lines=False, aniso=a)["T_avg60"] - target
    try:
        lgK = brentq(fK, math.log10(3e-11), math.log10(8.92e-10), xtol=0.02)
        out["K_fit"] = 10.0 ** lgK
    except ValueError:
        out["K_fit"] = None
    try:
        out["aniso_fit"] = brentq(fA, 1.0, 4.0, xtol=0.03)
    except ValueError:
        out["aniso_fit"] = None
    print("fitted: isotropic K", out["K_fit"], "m2; or Kh/Kv", out["aniso_fit"])
    variants = {}
    if out["K_fit"]:
        variants["isotropic"] = dict(K=out["K_fit"])
    if out["aniso_fit"]:
        variants["anisotropic"] = dict(aniso=out["aniso_fit"])
    variants["uncoupled"] = dict(buoyancy=False, viscosity=False)
    for vname, kw in variants.items():
        v = dict(kw=kw, flow={}, T0={}, T_in={}, spacing={})
        typ = sandbox(t_end=120 * MIN, snapshots=(20 * MIN, 30 * MIN, 40 * MIN, 60 * MIN) if vname == "isotropic" else (), **kw)
        v["typical"] = {k: x for k, x in typ.items() if k not in ("grid",)} if vname == "isotropic" else pick(typ)
        if vname == "isotropic":
            v["grid"] = typ["grid"]
        v["fronts"] = typ.get("fronts")
        for q in (0.2, 0.3, 0.4, 0.5, 0.6):
            r = typ if q == 0.5 else sandbox(q=q, t_end=(150 if q < 0.35 else 75) * MIN, lines=False, **kw)
            v["flow"][f"{q:.1f}"] = pick(r)
        for T0 in (50.0, 60.0, 70.0):
            r = typ if T0 == 60.0 else sandbox(T0=T0, lines=False, **kw)
            v["T0"][f"{T0:.0f}"] = pick(r)
        for Ti in (21.0, 31.0, 41.0):
            r = typ if Ti == 21.0 else sandbox(T_in=Ti, lines=False, **kw)
            v["T_in"][f"{Ti:.0f}"] = pick(r)
        for D in (0.90, 1.20):
            r = typ if D == 0.90 else sandbox(spacing=D, t_end=120 * MIN, lines=False, **kw)
            v["spacing"][f"{D:.2f}"] = pick(r)
        v["inject_bottom"] = pick(sandbox(inject="bottom", lines=False, t_end=120 * MIN, **kw))
        out[vname] = v
        print(f"--- {vname} {kw}")
        for grp in ("flow", "T0", "T_in", "spacing"):
            for k, r in v[grp].items():
                f_ = lambda x, d=1: "–" if x is None else f"{x:.{d}f}"
                print(f"   {grp:8s} {k:5s} onset {f_(r['onset_1K_min'])} min  T60 {f_(r['T60'])}  avg {f_(r['T_avg60'], 2)}  kW {f_(r['kW_avg60'], 2)}  decay {f_(r['decay60'])} %  x_onset {f_(r['x_onset_cyl'], 2)}")
        r = v["inject_bottom"]
        print(f"   injection below: onset {r['onset_1K_min']} T60 {r['T60']} avg {r['T_avg60']}")
    save("SG_calib", out)
    return out


def recommended():
    """The setting the paper recommends (0.4 m3/h, 31 C injection, 120 cm spacing) and the nominal one, run for three hours with each
    variant of the model, so that the whole decline is seen and not only the first hour."""
    cal = load("SG_calib")
    out = {}
    for vname in ("isotropic", "anisotropic", "uncoupled", "published"):
        kw = {} if vname == "published" else (cal.get(vname) or {}).get("kw")
        if kw is None:
            continue
        rows = {}
        for name, a in (("nominal: 0.5 m3/h, 21 C, 90 cm", {}), ("0.5 m3/h, 21 C, 120 cm", dict(spacing=1.20)),
                        ("recommended: 0.4 m3/h, 31 C, 120 cm", dict(q=0.4, T_in=31.0, spacing=1.20)),
                        ("0.2 m3/h, 21 C, 90 cm", dict(q=0.2))):
            r = sandbox(lines=False, t_end=180 * MIN, **a, **kw)
            rows[name] = pick(r)
            rows[name]["theta"] = [(T - r["T_in"]) / (r["T0"] - r["T_in"]) for T in r["T_prod"]]
            rows[name]["energy_to_onset_MJ"] = (CRHO_PAPER * r["q"] / 3600.0 * (r["T0"] - r["T_in"]) * r["onset_1K_min"] * MIN / 1e6) if r["onset_1K_min"] else None
            rows[name]["t_end_min"] = 180.0
            print(f"{vname:12s} {name:38s}", brief(r), "| energy to onset", None if rows[name]["energy_to_onset_MJ"] is None else round(rows[name]["energy_to_onset_MJ"], 1), "MJ")
        out[vname] = rows
    save("SG_recommended", out)
    return out


# ------------------------------------------------------------------------------------------------ field scale
def curves(t_end_min=300.0):
    """One curve for each geometry: retained share of the rise against x_D = t Q (rho c)_w / (D^3 (rho c)_b)."""
    L = LAB
    lin = sandbox(buoyancy=False, viscosity=False, lam=0.0, disp=0.0, t_end=t_end_min * MIN, lines=False)        # the box, advection only
    Q = L["q"] / 3600.0
    D = L["spacing"]
    xD = np.array(lin["t_min"]) * MIN * Q * lin["rhoc_w"] / (D ** 3 * lin["rhoc_b"])
    th_box = (np.array(lin["T_prod"]) - L["T_in"]) / (L["T0"] - L["T_in"])
    xg = np.geomspace(0.2, 40.0, 240)
    el = SingleWellDipole(Q=1.0, D=1.0, rhoc_b=1.0, T_rock_C=1.0, T_inj_C=0.0, rhoc_w=1.0)
    c = dict(x=xg.tolist(), unbounded=(1.0 - el.T_D(xg)).tolist(), box_x=xD[::2].tolist(), box=th_box[::2].tolist(),
             piston_x_box=math.pi * L["R"] ** 2 * D / D ** 3, R_over_D=L["R"] / D,
             note="box: the sandbox geometry (R/D = 0.56, closed 0.5 D above the injection screen); piston: uniform sweep of the cylinder between the screens")
    for th in (0.99, 0.9, 0.5):
        c[f"x_at_{th}"] = dict(unbounded=float(np.interp(1 - th, 1 - np.array(c["unbounded"]), xg)), box=float(np.interp(th, th_box[::-1], xD[::-1])))
    print("x_D at theta 0.99 / 0.9 / 0.5:", {k: {a: round(b, 3) for a, b in v.items()} for k, v in c.items() if k.startswith("x_at")}, "| box curve to x_D", round(float(xD[-1]), 2))
    return c


def refresh_curves():
    """Recompute the curves only and put them into the saved field results."""
    out = load("SG_field")
    out["curves"] = curves()
    save("SG_field", out)
    return out


def field():
    out = {}
    out["curves"] = curves()
    # (2) spacing a single well needs: first arrival and 10 % of the rise lost, 30 years
    rhoc_r = 2650.0 * 850.0
    rows = []
    for m in (5.0, 10.0, 20.0, 50.0, 100.0):
        cw = cbar_w(60.0, 183.0, 47.0e6)
        for yrs in (10.0, 30.0):
            D_bt = (m * cw * yrs * YEAR / (0.746 * rhoc_r)) ** (1 / 3)
            D_10 = (m * cw * yrs * YEAR / (float(dipole_arrival(0.1)) * rhoc_r)) ** (1 / 3)
            rows.append(dict(m=m, years=yrs, D_first_arrival_m=D_bt, D_10pct_m=D_10, V_ideal_km3=m * cw * yrs * YEAR / (0.8 * rhoc_r) / 1e9))
    out["spacing_needed"] = rows
    for r in rows:
        print("spacing needed", {k: round(v, 3) for k, v in r.items()})
    # (3) the type 11 example of the benchmark, with the same well, and three ways of sweeping the rock at its toe
    out["type11"] = type11_variants()
    # (4) does the weight of the cold water matter at field scale? the 150 m dipole in a block, both ways up
    out["field_buoyancy"] = field_buoyancy()
    save("SG_field", out)
    return out


def type11_variants():
    """The Eden-type example (4.8 km coaxial well, 183 C at the toe, 10 kg/s): the benchmark's three parallel fractures against a
    dipole between two zones 150 m apart (densely fractured, and with flowing fractures 30 m apart), and against the spacing that holds."""
    import dev_types as DT
    from aharrsim.systems import coaxial_well
    depth, m = 4800.0, 10.0
    T_rock = 15.0 + 35.0 * 4.8
    rhoc_r = 2650.0 * 850.0
    cw = cbar_w(60.0, T_rock, 47.0e6)
    Q = m / 950.0
    casing = [dict(top_m=0, bottom_m=3000, id_in=8.535, layers=[dict(od_in=9.625, material="steel"), dict(od_in=12.25, material="cement")]),
              dict(top_m=3000, bottom_m=depth, id_in=6.094, layers=[dict(od_in=7.0, material="steel"), dict(od_in=8.5, material="cement")])]
    tubing = [dict(top_m=0, bottom_m=depth, id_in=2.992, layers=[dict(od_in=3.5, k=45), dict(od_in=4.0, k=0.03), dict(od_in=4.5, k=45)])]
    D30 = (m * cw * 30.0 * YEAR / (float(dipole_arrival(0.1)) * rhoc_r)) ** (1 / 3)
    cases = {
        "three fractures, 150 m path (benchmark)": ParallelFractures(Q=Q, N=3, L_flow=150.0, B=150.0, S=30.0, k_r=3.2, rhoc_r=rhoc_r, T_rock_C=T_rock,
                                                                    T_inj_C=60.0, rhoc_w=cw * 950.0),
        "dipole 150 m, densely fractured": SingleWellDipole(Q=Q, D=150.0, rhoc_b=rhoc_r, T_rock_C=T_rock, T_inj_C=60.0, rhoc_w=cw * 950.0),
        "dipole 150 m, fractures 30 m apart": SingleWellDipole(Q=Q, D=150.0, rhoc_b=rhoc_r, T_rock_C=T_rock, T_inj_C=60.0, rhoc_w=cw * 950.0, S=30.0, k_r=3.2),
        f"dipole {D30:.0f} m, densely fractured": SingleWellDipole(Q=Q, D=D30, rhoc_b=rhoc_r, T_rock_C=T_rock, T_inj_C=60.0, rhoc_w=cw * 950.0),
    }
    out = dict(D_holds_m=D30, cases={})
    for name, el in cases.items():
        d = coaxial_well("Eden single-well EGS", depth, casing, tubing, m, 40.0, 1.0e6, [dict(k=3.2, rho=2650.0, cp=850.0)], 12.0,
                         dict(gradient_C_per_km=35.0), 30.0, bottom_element=el)
        sim, res = DT.run_case(d)
        s = DT.ts_summary(res)
        Qres = DT.Series(res.ts["t"], res.ts["Q_res"])
        po = {y: DT.plant_options(m, s[y]["h"], s[y]["p_MPa"] * 1e6, T_return_C=40.0) for y in s}
        yrs = sorted(s, key=float)
        row = dict(T_out={str(y): s[y]["T_C"] for y in yrs}, heat_MWth={str(y): (po[y] or {}).get("heat_MWth") for y in yrs},
                   binary_MWe={str(y): (po[y] or {}).get("binary_MWe") for y in yrs}, Q_res_MW={str(y): float(Qres(float(y) * YEAR)) / 1e6 for y in yrs},
                   T_res={str(y): float(el.T_out(60.0 + 273.15, float(y) * YEAR)) - 273.15 for y in yrs},
                   t_first_yr=(getattr(el, "t_bt", None) or getattr(el, "t_sweep", 0.0)) / YEAR, flags=DT.flags(res))
        tt = np.array([0.25, 0.5, 1, 2, 3, 5, 7, 10, 15, 20, 25, 30])
        row["series"] = dict(t_yr=tt.tolist(), T_res=[float(el.T_out(60.0 + 273.15, t * YEAR)) - 273.15 for t in tt])
        out["cases"][name] = row
        print(f"type 11 | {name:42s} T_out", {y: round(v, 1) for y, v in row["T_out"].items()}, "heat", {y: (round(v, 2) if v is not None else None) for y, v in row["heat_MWth"].items()},
              "first arrival / sweep", round(row["t_first_yr"], 2), "yr")
    return out


def field_buoyancy():
    """A 150 m dipole in a densely fractured block (radius 300 m, 600 m high), 10 kg/s, 183 C rock, 60 C water at 47 MPa:
    injection above or below the production zone, for permeabilities of 1e-14 to 1e-12 m2, against the flow that ignores temperature."""
    Wt = water_table(47.0e6, 20.0, 200.0)
    rows = {}
    D, Qm = 150.0, 10.0
    for K in (None, 1e-14, 1e-13, 1e-12):
        for inject in ("top", "bottom"):
            if K is None and inject == "bottom":
                continue
            b = AxisymmetricBlock(R=300.0, H=600.0, r_w=0.11, K=(K or 1e-13), phi=0.01, rhoc_s=2650.0 * 850.0, lam=3.2, T0=183.0, nr=40, nz=120, r1=4.0,
                                  buoyancy=K is not None, viscosity=K is not None, T_ref=121.5, water=Wt,
                                  rhoc_w=cbar_w(60.0, 183.0, 47.0e6) * float(np.interp(60.0, Wt["T"], Wt["rho"])))
            Q = Qm / float(np.interp(60.0, Wt["T"], Wt["rho"]))
            zu, zl = (300.0 + 0.5 * D - 10.0, 300.0 + 0.5 * D + 10.0), (300.0 - 0.5 * D - 10.0, 300.0 - 0.5 * D + 10.0)
            if inject == "top":
                b.add_screen(*zl, -Q); b.add_screen(*zu, +Q, T_in=60.0)
            else:
                b.add_screen(*zu, -Q); b.add_screen(*zl, +Q, T_in=60.0)
            t0 = time.time()
            r = b.run(30.0 * YEAR, dt_out=YEAR / 12.0, flow_update=YEAR / 50.0, cfl=0.6)
            ty = r["t"] / YEAR; T = r["T_prod"]
            th = (T - 60.0) / 123.0
            first = float(ty[np.argmax(th < 0.99)]) if np.any(th < 0.99) else None
            key = ("no coupling" if K is None else f"K {K:.0e}") + (", inject above" if inject == "top" else ", inject below")
            rows[key] = dict(K=K, inject=inject, t_99_yr=first, t_90_yr=(float(np.interp(0.9, th[::-1], ty[::-1])) if th.min() < 0.9 else None),
                             theta={str(y): float(np.interp(y, ty, th)) for y in (1, 2, 5, 10, 20, 30)}, wall_s=time.time() - t0,
                             series=dict(t_yr=ty[::3].tolist(), theta=th[::3].tolist()),
                             balance=float(np.max(np.abs(r["E_out"] + r["E_block"])) / np.max(np.abs(r["E_out"]))))
            print("field block |", f"{key:28s}", "1 % lost at", None if first is None else round(first, 2), "yr; theta", {k: round(v, 3) for k, v in rows[key]["theta"].items()},
                  f"{rows[key]['wall_s']:.0f} s")
    el = SingleWellDipole(Q=Qm / float(np.interp(60.0, Wt["T"], Wt["rho"])), D=D, rhoc_b=0.99 * 2650.0 * 850.0 + 0.01 * 4.3e6, T_rock_C=183.0, T_inj_C=60.0,
                          rhoc_w=cbar_w(60.0, 183.0, 47.0e6) * float(np.interp(60.0, Wt["T"], Wt["rho"])))
    rows["element (unbounded)"] = dict(t_99_yr=el.t_bt / YEAR, theta={str(y): float(1 - el.T_D(y * YEAR)[0]) for y in (1, 2, 5, 10, 20, 30)})
    print("field block | element: first arrival", round(el.t_bt / YEAR, 2), "yr; theta", {k: round(v, 3) for k, v in rows["element (unbounded)"]["theta"].items()})
    return rows


if __name__ == "__main__":
    what = sys.argv[1:] or ["checks", "lab", "calib", "rec", "field"]
    for w in what:
        dict(checks=checks, lab=lab, calib=calibrate, rec=recommended, field=field, curves=refresh_curves)[w]()
