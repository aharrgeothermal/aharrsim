"""Verification against analytical solutions, conservation and convergence.

V1  steady reverse circulation (exact analytical solution)  -> wellbore solver
V2  constant-temperature cylinder (Carslaw & Jaeger)        -> rock conduction solver
V3  Kabir et al. (1996) semi-analytical transient           -> classical cross-check
V4  energy conservation of the coupled solution
V5  grid / time-step convergence on a realistic 8 km case
V6  water/steam properties: IAPWS-95 check values and table accuracy
"""
from __future__ import annotations

import math
import sys
import time

import numpy as np

from common import save, run, C, DAY, YEAR
from aharrsim.analytic import (reverse_circulation_steady, forward_circulation_steady, jaeger_cylinder_flux,
                              jaeger_dimensionless_flux, kabir_reverse_transient)

CONST = {"constant": {"rho": 1000.0, "cp": 4180.0, "mu": 1.0e-3, "k": 0.6}}


# ----------------------------------------------------------------------------
def v1_steady():
    cases = [
        dict(label="6 km, insulated tubing (R_t = 5 K m/W)", L=6000, w=10.0, Rf=0.5, Rt=5.0),
        dict(label="6 km, bare steel tubing (R_t = 0.05 K m/W)", L=6000, w=10.0, Rf=0.5, Rt=0.05),
        dict(label="12 km, low flow (5 kg/s), R_t = 2 K m/W", L=12000, w=5.0, Rf=0.3, Rt=2.0),
    ]
    Ts, Tin, g, cp = 15.0, 40.0, 0.03, 4180.0
    out = []
    for cs in cases:
        rows = []
        for dz in [400.0, 200.0, 100.0, 50.0, 25.0, 12.5]:
            d = dict(
                name="V1", fluid=CONST,
                well=dict(depth_m=cs["L"],
                          casing=[dict(id_m=0.16, layers=[])],
                          tubing=[dict(id_m=0.076, layers=[dict(od_m=0.1143, k=1.0)], R_wall_K_m_W=cs["Rt"])]),
                rock=dict(layers=[dict(k=2.5, rho=2500, cp=1000)], surface_T_C=Ts,
                          geotherm=dict(gradient_C_per_km=g * 1000)),
                operation=dict(pressure_bc=dict(outlet_MPa=1.0),
                               schedule=[dict(duration_s=1.0, mdot_kg_s=cs["w"], T_in_C=Tin)]),
                numerics=dict(dz_m=dz, dt_initial_s=1.0),
                options=dict(gravity_work=False, include_pressure=False, rock_model="resistance",
                             R_formation_K_m_W=cs["Rf"], h_wall_override_W_m2K=1e12, fluid_storage=False),
            )
            sim, res = run(d)
            md = sim.s_f[:sim.Nw + 1]
            Ta_ex, Tt_ex = reverse_circulation_steady(md, cs["L"], cs["w"], cp, cs["Rf"], cs["Rt"], Tin, Ts, g)
            Ta = C(sim.fluid.ph(sim.p_a, sim.h_a).T)
            Tt = C(sim.fluid.ph(sim.p_t, sim.h_t).T)
            rows.append(dict(dz=dz, T_out=float(Tt[0]), T_out_exact=float(Tt_ex[0]),
                             err_out=float(Tt[0] - Tt_ex[0]),
                             err_max=float(max(np.max(np.abs(Ta - Ta_ex)), np.max(np.abs(Tt - Tt_ex))))))
            if dz == 50.0:
                prof = dict(md=md, Ta=Ta, Tt=Tt, Ta_exact=Ta_ex, Tt_exact=Tt_ex)
        errs = np.array([r["err_max"] for r in rows])
        order = np.log2(errs[:-1] / errs[1:])
        cs_out = dict(cs, rows=rows, observed_order=order.tolist(), profile=prof)
        out.append(cs_out)
        print(f"V1 {cs['label']}: Tout exact {rows[0]['T_out_exact']:.4f} C")
        for r in rows:
            print(f"    dz={r['dz']:6.1f}  Tout={r['T_out']:.5f}  err_out={r['err_out']:+.2e}  err_max={r['err_max']:.2e}")
        print("    observed order:", np.round(order, 2))
    save("V1_steady", out)


# ----------------------------------------------------------------------------
def v2_jaeger():
    k, rhoc, a = 2.5, 2.5e6, 0.1
    alpha = k / rhoc
    t_out_days = [0.01, 0.1, 1.0, 10.0, 100.0, 1000.0, 10950.0]
    out = []
    for growth in [1.3, 1.15, 1.08]:
        d = dict(
            name="V2", fluid=CONST,
            well=dict(depth_m=1000.0,
                      casing=[dict(id_m=2 * a, layers=[])],
                      tubing=[dict(id_m=0.076, layers=[dict(od_m=0.1143, k=1.0)], R_wall_K_m_W=1e9)]),
            rock=dict(layers=[dict(k=k, rho=2500.0, cp=rhoc / 2500.0)], surface_T_C=100.0,
                      geotherm=dict(gradient_C_per_km=0.0)),
            operation=dict(pressure_bc=dict(outlet_MPa=1.0),
                           schedule=[dict(duration_days=10950.0, mdot_kg_s=2000.0, T_in_C=20.0)]),
            numerics=dict(dz_m=100.0, r_far_m=400.0, radial_growth=growth, dt_initial_s=10.0, dt_growth=1.15,
                          dt_max_days=20.0, output_times_days=t_out_days),
            options=dict(gravity_work=False, include_pressure=False, h_wall_override_W_m2K=1e12,
                         fluid_storage=False),
        )
        sim, res = run(d)
        j = int(np.argmin(np.abs(sim.s_c[:sim.Nw] - 500.0)))
        rows = []
        for pr in res.profiles:
            Tabar = 0.5 * (pr["T_a"][j] + pr["T_a"][j + 1])
            dT = 373.15 - Tabar
            q = pr["q_rock"][j]
            q_ex = jaeger_cylinder_flux(pr["t"], k, alpha, a, dT)
            rows.append(dict(t_days=pr["t"] / DAY, tau=alpha * pr["t"] / a ** 2, q=q, q_exact=q_ex,
                             rel_err=q / q_ex - 1.0))
        out.append(dict(radial_growth=growth, M=sim.M, rows=rows))
        print(f"V2 radial growth {growth} (M={sim.M} radial cells)")
        for r in rows:
            print(f"    t={r['t_days']:9.2f} d  tau={r['tau']:10.3e}  q={r['q']:8.2f}  exact={r['q_exact']:8.2f}  rel={r['rel_err']:+.4f}")
    # continuous exact curve for plotting
    taus = np.logspace(-1, 6, 60)
    curve = [dict(tau=float(t), F=jaeger_dimensionless_flux(t)) for t in taus]
    save("V2_jaeger", dict(k=k, rhoc=rhoc, a=a, runs=out, curve=curve))


# ----------------------------------------------------------------------------
def v3_kabir():
    L, w, cp, k, rhoc, rw, Rt = 6000.0, 10.0, 4180.0, 2.5, 2.4e6, 0.108, 2.0
    Ts, Tin, g = 15.0, 40.0, 0.03
    times = [1, 10, 30, 100, 365, 1825, 3650, 10950]
    d = dict(
        name="V3", fluid=CONST,
        well=dict(depth_m=L, casing=[dict(id_m=2 * rw, layers=[])],
                  tubing=[dict(id_m=0.076, layers=[dict(od_m=0.1143, k=1.0)], R_wall_K_m_W=Rt)]),
        rock=dict(layers=[dict(k=k, rho=2400.0, cp=rhoc / 2400.0)], surface_T_C=Ts,
                  geotherm=dict(gradient_C_per_km=g * 1000)),
        operation=dict(pressure_bc=dict(outlet_MPa=1.0),
                       schedule=[dict(duration_days=10950, mdot_kg_s=w, T_in_C=Tin)]),
        numerics=dict(dz_m=50.0, dt_initial_s=60.0, output_times_days=times),
        options=dict(gravity_work=False, include_pressure=False, h_wall_override_W_m2K=1e12),
    )
    sim, res = run(d)
    rows = []
    for t in times:
        Tm = C(res.at(t * DAY, "T_out"))
        Ta_k, Tt_k = kabir_reverse_transient(t * DAY, np.array([0.0]), L, w, cp, k, k / rhoc, rw, 1e12, Rt,
                                             Tin, Ts, g)
        rows.append(dict(t_days=t, T_out_model=float(Tm), T_out_kabir=float(Tt_k[0]),
                         diff=float(Tm - Tt_k[0])))
        print(f"V3 t={t:6d} d  model {Tm:7.3f} C   Kabir {Tt_k[0]:7.3f} C   diff {Tm - Tt_k[0]:+.3f} K")
    tt = res.t_days
    ts = dict(t_days=tt, T_out=C(res.ts["T_out"]))
    kab = [float(kabir_reverse_transient(t * DAY, np.array([0.0]), L, w, cp, k, k / rhoc, rw, 1e12, Rt,
                                         Tin, Ts, g)[1][0]) for t in tt]
    ts["T_out_kabir"] = kab
    save("V3_kabir", dict(rows=rows, series=ts))


# ----------------------------------------------------------------------------
def realistic_8km(dz=50.0, growth=1.15, dtg=1.3, dtmax=30.0, dt0=300.0, scheme="box", tmax=3650.0,
                  outs=(1, 30, 365, 3650)):
    return dict(
        name="8 km reference",
        well=dict(depth_m=8000.0,
                  casing=[dict(top_m=0, bottom_m=3000, id_in=8.681, layers=[dict(od_in=9.625, material="steel"),
                                                                            dict(od_in=12.25, material="cement")]),
                          dict(top_m=3000, bottom_m=8000, id_in=6.184, layers=[dict(od_in=7.0, material="steel"),
                                                                               dict(od_in=8.5, material="cement")])],
                  tubing=[dict(top_m=0, bottom_m=5000, id_in=2.992,
                               layers=[dict(od_in=3.5, k=45), dict(od_in=4.0, k=0.03), dict(od_in=4.5, k=45)]),
                          dict(top_m=5000, bottom_m=8000, id_in=2.992, layers=[dict(od_in=3.5, k=45)])]),
        materials=dict(steel=dict(k=45, rho=7850, cp=490), cement=dict(k=1.0, rho=1900, cp=1000)),
        rock=dict(layers=[dict(top_m=0, bottom_m=3000, k=2.4, rho=2550, cp=900),
                          dict(top_m=3000, bottom_m=1e5, k=3.0, rho=2700, cp=850)],
                  surface_T_C=20.0, geotherm=dict(gradient_C_per_km=30.0)),
        operation=dict(pressure_bc=dict(outlet_MPa=2.0),
                       schedule=[dict(duration_days=tmax, mdot_kg_s=10.0, T_in_C=50.0)]),
        numerics=dict(dz_m=dz, radial_growth=growth, dt_growth=dtg, dt_max_days=dtmax, dt_initial_s=dt0,
                      scheme=scheme, output_times_days=list(outs)),
    )


def v4_v5_convergence():
    t_check = [1, 30, 365, 3650]
    variants = [
        ("dz 200 m", dict(dz=200.0)), ("dz 100 m", dict(dz=100.0)), ("dz 50 m (base)", dict()),
        ("dz 25 m", dict(dz=25.0)),
        ("radial growth 1.30", dict(growth=1.30)), ("radial growth 1.08", dict(growth=1.08)),
        ("dt growth 1.6, dt_max 90 d", dict(dtg=1.6, dtmax=90.0)),
        ("dt growth 1.1, dt_max 5 d", dict(dtg=1.1, dtmax=5.0)),
        ("upwind scheme", dict(scheme="upwind")),
    ]
    out = []
    base = None
    bal = []
    for label, kw in variants:
        t0 = time.time()
        sim, res = run(realistic_8km(**kw))
        vals = {t: float(C(res.at(t * DAY, "T_out"))) for t in t_check}
        Q = {t: float(res.at(t * DAY, "Q_wellhead") / 1e3) for t in t_check}
        resid = res.ts["Q_wellhead"] - (res.ts["Q_rock"] - res.ts["dE_fluid"] - res.ts["Q_ke"])
        rel = float(np.max(np.abs(resid)) / np.max(np.abs(res.ts["Q_wellhead"])))
        bal.append(rel)
        row = dict(label=label, T_out=vals, Q_kW=Q, cells=f"{sim.Nw} x {sim.M}", steps=len(res.ts["t"]),
                   wall_s=time.time() - t0, energy_residual_rel=rel)
        if base is None and "base" in label:
            base = row
        out.append(row)
        print(f"V5 {label:28s} " + "  ".join(f"{t}d:{vals[t]:7.3f}" for t in t_check) +
              f"   steps={row['steps']} wall={row['wall_s']:.1f}s  balance={rel:.1e}")
    save("V5_convergence", dict(rows=out, t_check=t_check))
    # time series of the base run for energy-balance plot
    sim, res = run(realistic_8km())
    save("V4_energy", dict(t_days=res.t_days, Q_wellhead=res.ts["Q_wellhead"] / 1e3,
                           Q_rock=res.ts["Q_rock"] / 1e3, dE=res.ts["dE_fluid"] / 1e3,
                           residual=(res.ts["Q_wellhead"] - res.ts["Q_rock"] + res.ts["dE_fluid"] + res.ts["Q_ke"]) / 1e3,
                           max_rel_all=max(bal)))


# ----------------------------------------------------------------------------
def v6_properties():
    import CoolProp
    from CoolProp import AbstractState
    from aharrsim.props import WaterIAPWS
    AS = AbstractState("HEOS", "Water")
    table7 = [(300, 0.9965560e3, 0.992418352e-1, 0.413018112e1, 0.150151914e4, 0.393062643),
              (300, 0.1005308e4, 0.200022515e2, 0.406798347e1, 0.153492501e4, 0.387405401),
              (300, 0.1188202e4, 0.700004704e3, 0.346135580e1, 0.244357992e4, 0.132609616),
              (500, 0.4350000, 0.999679423e-1, 0.150817541e1, 0.548314253e3, 0.794488271e1),
              (500, 0.4532000e1, 0.999938125, 0.166991025e1, 0.535739001e3, 0.682502725e1),
              (500, 0.8380250e3, 0.100003858e2, 0.322106219e1, 0.127128441e4, 0.256690919e1),
              (500, 0.1084564e4, 0.700000405e3, 0.307437693e1, 0.241200877e4, 0.203237509e1),
              (647, 0.3580000e3, 0.220384756e2, 0.618315728e1, 0.252145078e3, 0.432092307e1),
              (900, 0.2410000, 0.100062559, 0.175890657e1, 0.724027147e3, 0.916653194e1),
              (900, 0.5261500e2, 0.200000690e2, 0.193510526e1, 0.698445674e3, 0.659070225e1),
              (900, 0.8707690e3, 0.700000006e3, 0.266422350e1, 0.201933608e4, 0.417223802e1)]
    worst95 = 0.0
    for T, rho, p, cv, w, s in table7:
        AS.update(CoolProp.DmassT_INPUTS, rho, T)
        errs = [AS.p() / 1e6 / p - 1, AS.cvmass() / 1e3 / cv - 1, AS.speed_sound() / w - 1, AS.smass() / 1e3 / s - 1]
        worst95 = max(worst95, max(abs(e) for e in errs))
    # table (bicubic) vs exact IAPWS-95 over the operating envelope
    W = WaterIAPWS()
    rng = np.random.default_rng(1)
    n = 4000
    ps = 10 ** rng.uniform(np.log10(0.1e6), np.log10(150e6), n)
    Ts = rng.uniform(275.0, 690.0, n)
    worst = {r: dict(T=0.0, rho=0.0, cp=0.0, mu=0.0, k=0.0, n=0) for r in ("liquid", "vapour")}
    for p, T in zip(ps, Ts):
        AS.update(CoolProp.PT_INPUTS, p, T)
        ph = AS.phase()
        if ph == CoolProp.iphase_twophase:
            continue
        h = AS.hmass()
        ref = (AS.T(), AS.rhomass(), AS.cpmass(), AS.viscosity(), AS.conductivity())
        if abs(p - 22.064e6) < 1.5e6 and abs(T - 647.1) < 15:
            continue
        npr = W.ph(np.array([p]), np.array([h]))
        if npr.two_phase[0]:
            continue
        reg = "vapour" if ph in (CoolProp.iphase_gas, CoolProp.iphase_supercritical_gas) else "liquid"
        wr = worst[reg]
        wr["n"] += 1
        wr["T"] = max(wr["T"], abs(npr.T[0] - ref[0]))
        wr["rho"] = max(wr["rho"], abs(npr.rho[0] / ref[1] - 1))
        wr["cp"] = max(wr["cp"], abs(npr.cp[0] / ref[2] - 1))
        wr["mu"] = max(wr["mu"], abs(npr.mu[0] / ref[3] - 1))
        wr["k"] = max(wr["k"], abs(npr.k[0] / ref[4] - 1))
    cnt = sum(w_["n"] for w_ in worst.values())
    # compressed liquid close to saturation (0.5-20 MPa, 0.2-5 K subcooled)
    near = 0.0
    for p in np.geomspace(0.5e6, 20e6, 120):
        AS.update(CoolProp.PQ_INPUTS, p, 0.0)
        Tsat = AS.T()
        for dT in (0.2, 1.0, 5.0):
            AS.update(CoolProp.PT_INPUTS, p, Tsat - dT)
            h = AS.hmass()
            npr = W.ph(np.array([p]), np.array([h]))
            near = max(near, abs(npr.T[0] - AS.T()))
    worst["near_saturation_T"] = near
    # saturation table check
    worst_sat = 0.0
    for p in np.geomspace(1e4, 21e6, 300):
        AS.update(CoolProp.PQ_INPUTS, p, 0.0)
        worst_sat = max(worst_sat, abs(W.Tsat(p) - AS.T()))
    print(f"V6 IAPWS-95 Table 7 worst relative deviation {worst95:.1e}")
    for reg, wr in worst.items():
        if reg == "near_saturation_T":
            print(f"V6 subcooled liquid within 5 K of saturation: worst T error {wr:.3f} K")
            continue
        print(f"V6 bicubic tables vs IAPWS-95, {reg} ({wr['n']} states): " +
              ", ".join(f"{k} {v:.1e}" for k, v in wr.items() if k != "n"))
    print(f"V6 saturation-table Tsat worst deviation {worst_sat:.1e} K")
    save("V6_properties", dict(table7_worst_rel=worst95, table_states=cnt, table_worst=worst,
                               sat_worst_K=worst_sat))


# ----------------------------------------------------------------------------
def v7_direction():
    """Reversed flow (injection down the tubing, return up the annulus).

    (a) steady forward-circulation problem against its closed-form solution, grid convergence;
    (b) with a uniform rock temperature the two flow directions are mirror images of each other,
        so with constant fluid properties they must give the same outlet temperature (with
        IAPWS water they differ slightly because the pressure, and so the properties, are not
        mirror-symmetric)."""
    cases = [dict(label="6 km, insulated tubing (R_t = 5 K m/W)", L=6000, w=10.0, Rf=0.5, Rt=5.0),
             dict(label="6 km, bare steel tubing (R_t = 0.05 K m/W)", L=6000, w=10.0, Rf=0.5, Rt=0.05),
             dict(label="12 km, low flow (5 kg/s), R_t = 2 K m/W", L=12000, w=5.0, Rf=0.3, Rt=2.0)]
    Ts, Tin, g, cp = 15.0, 40.0, 0.03, 4180.0
    steady = []
    for cs in cases:
        rows = []
        for dz in [400.0, 200.0, 100.0, 50.0, 25.0, 12.5]:
            d = dict(name="V7", fluid=CONST,
                     well=dict(depth_m=cs["L"], flow_direction="tubing_down",
                               casing=[dict(id_m=0.16, layers=[])],
                               tubing=[dict(id_m=0.076, layers=[dict(od_m=0.1143, k=1.0)], R_wall_K_m_W=cs["Rt"])]),
                     rock=dict(layers=[dict(k=2.5, rho=2500, cp=1000)], surface_T_C=Ts,
                               geotherm=dict(gradient_C_per_km=g * 1000)),
                     operation=dict(pressure_bc=dict(outlet_MPa=1.0),
                                    schedule=[dict(duration_s=1.0, mdot_kg_s=cs["w"], T_in_C=Tin)]),
                     numerics=dict(dz_m=dz, dt_initial_s=1.0),
                     options=dict(gravity_work=False, include_pressure=False, rock_model="resistance",
                                  R_formation_K_m_W=cs["Rf"], h_wall_override_W_m2K=1e12, fluid_storage=False))
            sim, res = run(d)
            md = sim.s_f[:sim.Nw + 1]
            Ta_ex, Tt_ex = forward_circulation_steady(md, cs["L"], cs["w"], cp, cs["Rf"], cs["Rt"], Tin, Ts, g)
            Ta = C(sim.fluid.ph(sim.p_a, sim.h_a).T)
            Tt = C(sim.fluid.ph(sim.p_t, sim.h_t).T)
            rows.append(dict(dz=dz, T_out=float(Ta[0]), T_out_exact=float(Ta_ex[0]), err_out=float(Ta[0] - Ta_ex[0]),
                             err_max=float(max(np.max(np.abs(Ta - Ta_ex)), np.max(np.abs(Tt - Tt_ex))))))
        errs = np.array([r["err_max"] for r in rows])
        order = np.log2(errs[:-1] / errs[1:])
        steady.append(dict(cs, rows=rows, observed_order=order.tolist()))
        print(f"V7 {cs['label']}: Tout exact {rows[0]['T_out_exact']:.4f} C")
        for r in rows:
            print(f"    dz={r['dz']:6.1f}  Tout={r['T_out']:.5f}  err_max={r['err_max']:.2e}")
        print("    observed order:", np.round(order, 2))

    def sym_case(direction, rock_model, fluid):
        d = dict(name="V7 symmetry", fluid=fluid,
                 well=dict(depth_m=4000.0, flow_direction=direction,
                           casing=[dict(id_in=8.535, layers=[dict(od_in=9.625, material="steel"), dict(od_in=12.25, material="cement")])],
                           tubing=[dict(id_in=2.992, layers=[dict(od_in=3.5, k=45.0), dict(od_in=4.0, k=0.03), dict(od_in=4.5, k=45.0)])]),
                 materials=dict(steel=dict(k=45.0, rho=7850.0, cp=490.0), cement=dict(k=1.0, rho=1900.0, cp=1000.0)),
                 rock=dict(layers=[dict(k=2.8, rho=2650.0, cp=900.0)], surface_T_C=150.0, geotherm=dict(gradient_C_per_km=0.0)),
                 operation=dict(pressure_bc=dict(outlet_MPa=5.0), schedule=[dict(duration_days=3653.0, mdot_kg_s=10.0, T_in_C=40.0)]),
                 numerics=dict(dz_m=50.0, dt_initial_s=300.0, dt_max_days=60.0, output_times_days=[365.25, 3653.0]),
                 options=dict(rock_model=rock_model, R_formation_K_m_W=0.3, gravity_work=fluid != CONST,
                              include_pressure=fluid != CONST))
        return d
    sym = {}
    for fluid, label in ((CONST, "constant properties"), ("iapws", "IAPWS-95")):
        for rock_model in ("resistance", "grid"):
            r = {}
            for direction in ("annulus_down", "tubing_down"):
                sim, res = run(sym_case(direction, rock_model, fluid))
                r[direction] = res
            ta, tb = r["annulus_down"].ts, r["tubing_down"].ts
            late = ta["t"] > 30 * DAY
            dT = np.abs(ta["T_out"] - tb["T_out"])[late]
            dQ = np.abs(ta["Q_wellhead"] - tb["Q_wellhead"])[late]
            print(f"V7 symmetry, {label}, {rock_model}: after 30 d max |dT_out| {dT.max():.2e} K, max |dQ| {dQ.max():.3g} W; "
                  f"T_out at end {C(ta['T_out'][-1]):.4f} / {C(tb['T_out'][-1]):.4f} C")
            sym[f"{label} / {rock_model}"] = dict(max_dT_after_30d=float(dT.max()), max_dQ_after_30d=float(dQ.max()),
                                                  T_end=[float(C(ta["T_out"][-1])), float(C(tb["T_out"][-1]))])
    save("V7_direction", dict(steady=steady, symmetry=sym))


if __name__ == "__main__":
    which = sys.argv[1:] or ["v1", "v2", "v3", "v5", "v6", "v7"]
    for w in which:
        dict(v1=v1_steady, v2=v2_jaeger, v3=v3_kabir, v5=v4_v5_convergence, v6=v6_properties, v7=v7_direction)[w]()
