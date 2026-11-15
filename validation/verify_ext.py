"""Verification of the extended well models (v2.0):

V9   single-channel production / injection wells: exact steady solutions (second-order convergence), transient
     production well against Ramey / Hasan-Kabir, pumped well, and the measured two-phase wells of
     twophase.py run through the production mode.
V8   U-tube borehole heat exchanger against pygfunction (finite-line-source g-function + multipole pipe model).
B4   U-loop (injector - lateral - producer chain) against the GeoCLUSTER 'utube' database (SBT model, Sandia / NREL).

    python verify_ext.py            # all (~20 min)
    python verify_ext.py v9 v8 b4   # selected
V10  reservoir elements (inversion vs erfc / mpmath incl. sharp fronts, doublet breakthrough, cyclic vs Lauwerier, heat-store balance).
V11  plant models vs the GEOPHIRES v3 utilisation-efficiency correlations.
Results -> results/V9_single_well.json, V9b_twophase_production_mode.json, V8_utube_pygfunction.json,
           B4_geocluster_utube.json, V10_reservoir.json, V11_plant.json
"""
from __future__ import annotations

import sys


def v9_single_well():
    import sys, math
    import numpy as np
    from common import run, C, save
    from verify import CONST
    from aharrsim.analytic import hasan_kabir_f
    import CoolProp.CoolProp as CP

    out = {}
    # ---------- (a) exact steady solutions with a fixed formation resistance -----------
    Ts, g, cp = 15.0, 0.03, 4180.0
    def case(mode, L, w, Rf, dz, T_feed=None, T_in=40.0):
        d = dict(name="V9", fluid=CONST,
                 well=dict(depth_m=L, mode=mode, casing=[dict(id_m=0.16, layers=[])]),
                 rock=dict(layers=[dict(k=2.5, rho=2500, cp=1000)], surface_T_C=Ts, geotherm=dict(gradient_C_per_km=g * 1000)),
                 operation=dict(pressure_bc=dict(outlet_MPa=1.0) if mode == "production" else dict(inlet_MPa=1.0),
                                schedule=[dict(duration_s=1.0, mdot_kg_s=w, T_in_C=T_in)]),
                 numerics=dict(dz_m=dz, dt_initial_s=1.0),
                 options=dict(gravity_work=False, include_pressure=False, rock_model="resistance",
                              R_formation_K_m_W=Rf, h_wall_override_W_m2K=1e12, fluid_storage=False))
        if mode == "production":
            d["well"]["inflow"] = dict(T_C=T_feed)
        return d
    steady = []
    for mode, L, w, Rf in (("production", 3000.0, 20.0, 0.3), ("production", 2000.0, 5.0, 0.15), ("injection", 3000.0, 20.0, 0.3)):
        rows = []
        for dz in (400.0, 200.0, 100.0, 50.0, 25.0):
            T_feed = Ts + g * L + 20.0            # feed 20 K above the geotherm at TD
            sim, res = run(case(mode, L, w, Rf, dz, T_feed=T_feed))
            A = w * cp * Rf
            md = sim.s_f[:sim.Nw + 1]
            if mode == "production":
                zeta = L - md
                T_ex = (Ts + g * L) - g * zeta + g * A + (T_feed - (Ts + g * L) - g * A) * np.exp(-zeta / A)
            else:
                T_ex = Ts + g * md - g * A + (40.0 - Ts + g * A) * np.exp(-md / A)
            T = C(sim.fluid.ph(sim.p_a, sim.h_a).T)
            rows.append(dict(dz=dz, err_max=float(np.max(np.abs(T - T_ex))), T_out=float(T[0] if mode == "production" else T[-1]),
                             T_exact=float(T_ex[0] if mode == "production" else T_ex[-1])))
        errs = np.array([r["err_max"] for r in rows]); order = np.log2(errs[:-1] / errs[1:])
        print(f"{mode} L {L:.0f} w {w}: exact T_out {rows[0]['T_exact']:.4f}; errors", [f"{r['err_max']:.2e}" for r in rows], "order", np.round(order, 2))
        steady.append(dict(mode=mode, L=L, w=w, Rf=Rf, rows=rows, order=order.tolist()))
    out["steady"] = steady

    # ---------- (b) transient production well vs Ramey / Hasan-Kabir ----------------
    # 2500 m well, 9-5/8" casing in 12-1/4" hole, feed at 160 C, 40 kg/s; compare wellhead T at 1 d .. 10 yr
    L, w, T_feed = 2500.0, 40.0, 160.0
    d = dict(name="V9b", fluid="iapws",
             well=dict(depth_m=L, mode="production", inflow=dict(T_C=T_feed, p_res_MPa=22.0, PI_kg_s_per_MPa=20.0),
                       casing=[dict(id_in=8.535, layers=[dict(od_in=9.625, material="steel"), dict(od_in=12.25, material="cement")])]),
             materials=dict(steel=dict(k=45.0, rho=7850.0, cp=490.0), cement=dict(k=1.0, rho=1900.0, cp=1000.0)),
             rock=dict(layers=[dict(k=2.5, rho=2600.0, cp=900.0)], surface_T_C=15.0, geotherm=dict(gradient_C_per_km=58.0)),
             operation=dict(pressure_bc=dict(outlet_MPa=1.5), schedule=[dict(duration_days=3653.0, mdot_kg_s=w, T_in_C=40.0)]),
             numerics=dict(dz_m=50.0, dt_initial_s=60.0, dt_max_days=30.0, output_times_days=[1, 30, 365, 3653]),
             options=dict(gravity_work=True))
    sim, res = run(d)
    ts = res.ts
    # Ramey: T_wh = T_geo(0) + g A (1 - e^{-L/A}) ... with A(t) = w cp [k + r_co U f(t)] / (2 pi r_co U k), U from casing+cement, f = Hasan-Kabir
    k_r, rho_r, cp_r = 2.5, 2600.0, 900.0
    alpha = k_r / (rho_r * cp_r)
    r_ci, r_co, r_h = 0.5 * 8.535 * 0.0254, 0.5 * 9.625 * 0.0254, 0.5 * 12.25 * 0.0254
    R_comp = math.log(r_co / r_ci) / (2 * math.pi * 45.0) + math.log(r_h / r_co) / (2 * math.pi * 1.0)
    cpw = 4300.0
    ramey = []
    for tdays in (1, 30, 365, 3653):
        t = tdays * 86400.0
        tD = alpha * t / r_h ** 2
        f_t = hasan_kabir_f(tD)
        R_tot = R_comp + f_t / (2 * math.pi * k_r)          # film neglected
        A = w * cpw * R_tot
        gg = 0.058
        T_wh = (15.0) + gg * A + (T_feed - (15.0 + gg * L) - gg * A) * math.exp(-L / A)
        i = int(np.argmin(np.abs(ts["t"] - t)))
        ramey.append(dict(t_days=tdays, T_wh_model=float(C(ts["T_out"][i])), T_wh_ramey=float(T_wh), p_wf_MPa=float(ts["p_in"][i] / 1e6),
                          margin_MPa=float(ts["inflow_margin"][i] / 1e6), x_out=float(ts["x_out"][i])))
        print(f"  t {tdays:5d} d: model T_wh {C(ts['T_out'][i]):7.2f}  Ramey/HK {T_wh:7.2f}  p_wf {ts['p_in'][i]/1e6:.2f} MPa  margin {ts['inflow_margin'][i]/1e6:+.2f} MPa  x_out {ts['x_out'][i]:.3f}")
    eb = np.max(np.abs(ts["Q_wellhead"] - ts["Q_rock"] - ts["Q_res"] + ts["dE_fluid"] + ts["Q_ke"] + ts["Q_pe"]))
    print("  energy residual", eb, "W; unconverged", res.meta.get("unconverged_steps", 0), "steps", ts["t"].size)
    # like-for-like with Ramey (pure heat exchange): constant properties, no pressure effects
    d_cf = dict(d); d_cf["fluid"] = {"constant": {"rho": 917.0, "cp": 4300.0, "mu": 1.8e-4, "k": 0.68}}
    d_cf["options"] = dict(gravity_work=False)
    sim_cf, res_cf = run(d_cf)
    for r in ramey:
        i = int(np.argmin(np.abs(res_cf.ts["t"] - r["t_days"] * 86400.0)))
        r["T_wh_model_constant_fluid"] = float(C(res_cf.ts["T_out"][i]))
        print(f"  t {r['t_days']:5d} d: constant-property model {r['T_wh_model_constant_fluid']:7.2f} vs Ramey/HK {r['T_wh_ramey']:7.2f}"
              f"  (IAPWS with expansion cooling {r['T_wh_model']:7.2f})")
    out["ramey"] = dict(rows=ramey, energy_residual_W=float(eb))

    # ---------- (c) pumped well: 2000 m, 150 C feed, pump at 400 m ------------------
    d2 = dict(name="V9c", fluid="iapws",
              well=dict(depth_m=2000.0, mode="production", inflow=dict(T_C=150.0, p_res_MPa=17.0, PI_kg_s_per_MPa=10.0),
                        pump=dict(md_m=1000.0, efficiency=0.7),
                        casing=[dict(id_in=8.535, layers=[dict(od_in=9.625, material="steel"), dict(od_in=12.25, material="cement")])]),
              materials=dict(steel=dict(k=45.0, rho=7850.0, cp=490.0), cement=dict(k=1.0, rho=1900.0, cp=1000.0)),
              rock=dict(layers=[dict(k=2.5, rho=2600.0, cp=900.0)], surface_T_C=15.0, geotherm=dict(gradient_C_per_km=67.5)),
              operation=dict(pressure_bc=dict(outlet_MPa=1.0), schedule=[dict(duration_days=365.0, mdot_kg_s=60.0, T_in_C=40.0)]),
              numerics=dict(dz_m=50.0, dt_initial_s=60.0, dt_max_days=30.0, output_times_days=[365]))
    sim2, res2 = run(d2)
    ts2 = res2.ts
    p = res2.profiles[-1]
    print(f"  pumped well: T_wh {C(ts2['T_out'][-1]):.1f} C, pump dp {ts2['pump_dp'][-1]/1e6:.2f} MPa, pump {ts2['pump_kW'][-1]:.0f} kW, "
          f"p_wf {ts2['p_in'][-1]/1e6:.2f} MPa (p_res - m/PI = {17 - 60/10:.2f}), suction {sim2.p_suction/1e6:.2f} MPa, two-phase at suction {sim2.pump_two_phase}")
    # hand check: liquid column 0-1000 m at ~150 C (rho ~ 917): dp ~ 917*9.81*400 = 3.6 MPa + friction
    print("  pressure profile (MD, p MPa):", [(round(float(a)), round(float(b) / 1e6, 2)) for a, b in zip(p["md"][::8], p["p_a"][::8])])
    out["pumped"] = dict(T_wh=float(C(ts2["T_out"][-1])), pump_dp_MPa=float(ts2["pump_dp"][-1] / 1e6), pump_kW=float(ts2["pump_kW"][-1]),
                         p_wf_MPa=float(ts2["p_in"][-1] / 1e6), p_suction_MPa=float(sim2.p_suction / 1e6))
    save("V9_single_well", out)


def v9b_twophase():
    import sys, math
    import numpy as np
    from common import run, C, load, save
    from twophase import WELLS, FT
    import CoolProp.CoolProp as CP

    out = []
    for well in WELLS:
        z0 = well["z"][0]
        zmax = well["z_max"]
        # the model well starts at the shallowest survey point z0 (wellhead pressure known there)
        L = zmax - z0
        strings = []
        for a, b, d_ in well["strings"]:
            if b <= z0 or a >= zmax:
                continue
            strings.append(dict(top_m=max(0.0, a - z0), bottom_m=min(b, zmax) - z0, id_m=d_,
                                layers=[dict(od_m=d_ + 0.02, material="steel"), dict(od_m=d_ + 0.08, material="cement")]))
        strings[0]["top_m"] = 0.0
        strings[-1]["bottom_m"] = L
        # feed enthalpy: the measured surveys are adiabatic with h + g z = const; the model does the same
        h_feed = well["h"] + 9.81 * L
        # rock at the flowing temperature (long-flowing well): no heat exchange, as in the adiabatic survey analysis
        if well["T"] is not None:
            geo = dict(table=[[float(zz - z0), float(TT)] for zz, TT in zip(well["z"], well["T"]) if zz >= z0])
            Tsurf = float(well["T"][0])
        else:
            geo = dict(gradient_C_per_km=0.0)
            Tsurf = float(CP.PropsSI("T", "P", (well["P"][0] + 1.01325) * 1e5, "Q", 0, "Water") - 273.15)
        d = dict(name=well["name"], fluid="iapws",
                 well=dict(depth_m=L, mode="production", inflow=dict(h_J_kg=h_feed), casing=strings, roughness_m=0.0006 * FT),
                 materials=dict(steel=dict(k=45.0, rho=7850.0, cp=490.0), cement=dict(k=1.0, rho=1900.0, cp=1000.0)),
                 rock=dict(layers=[dict(k=2.5, rho=2600.0, cp=900.0)], surface_T_C=Tsurf, geotherm=geo),
                 operation=dict(pressure_bc=dict(outlet_MPa=(well["P"][0] + 1.01325) / 10.0),
                                schedule=[dict(duration_days=1.0, mdot_kg_s=well["W"], T_in_C=40.0)]),
                 numerics=dict(dz_m=10.0, dt_initial_s=60.0, dt_max_days=1.0, output_times_days=[1.0], two_phase_model="drift_flux"))
        sim, res = run(d)
        prof = res.profiles[-1]
        md, pa = prof["md"] + z0, prof["p_a"] / 1e5 - 1.01325
        z = np.array(well["z"], float); P = np.array(well["P"], float)
        sel = (z > z0) & (z <= zmax + 1e-6)
        pm = np.interp(z[sel], md, pa)
        err = pm - P[sel]
        pct = 100 * err / P[sel]
        row = dict(name=well["name"], n=int(sel.sum()), mean_bar=float(err.mean()), sd_bar=float(err.std()), mean_pct=float(pct.mean()),
                   sd_pct=float(pct.std()), p_bottom_model=float(pm[-1]), p_bottom_meas=float(P[sel][-1]), T_wh=float(C(res.ts["T_out"][-1])),
                   x_wh=float(res.ts["x_out"][-1]), flags=dict(unconv=int(res.meta.get("unconverged_steps", 0)), choke=int(np.sum(res.ts["choke_warning"] > 0))))
        if well["T"] is not None:
            Tm = np.interp(z[sel], md, prof["T_a"] - 273.15)
            row["T_mean_err"] = float((Tm - np.array(well["T"], float)[sel]).mean())
        step = max(1, len(md) // 60)
        row["profile"] = dict(md=[float(v) for v in md[::step]] + [float(md[-1])], p_bar=[float(v) for v in pa[::step]] + [float(pa[-1])],
                              T_C=[float(v) for v in (prof["T_a"] - 273.15)[::step]] + [float(prof["T_a"][-1] - 273.15)],
                              x=[float(v) for v in prof["x_a"][::step]] + [float(prof["x_a"][-1])])
        row["measured"] = dict(md=[float(v) for v in z[sel]], p_bar=[float(v) for v in P[sel]],
                               T_C=([float(v) for v in np.array(well["T"], float)[sel]] if well["T"] is not None else None))
        row["mdot"] = float(well["W"]); row["h_kJ"] = float(well["h"] / 1e3)
        out.append(row)
        print(f"{well['name']:22s} n {row['n']:2d}  dp mean {row['mean_bar']:+6.2f} bar ({row['mean_pct']:+5.1f} %)  sd {row['sd_pct']:4.1f} %  "
              f"p_bottom {row['p_bottom_model']:6.1f} vs {row['p_bottom_meas']:6.1f}  T_wh {row['T_wh']:6.1f}  x {row['x_wh']:.3f}  flags {row['flags']}")
    prev = load("T1_twophase")
    print("\nearlier standalone march (T1):")
    for w in (prev["wells"] if isinstance(prev, dict) and "wells" in prev else prev):
        try:
            print(f"  {w['name']:22s} {w.get('drift_flux', w).get('mean_pct', w.get('mean_pct')):+5.1f} %")
        except Exception as e:
            print("  ", w.get("name"), {k: w[k] for k in list(w)[:6]})
    save("V9b_twophase_production_mode", out)


def v8_utube():
    import sys, math
    import numpy as np
    import pygfunction as gt
    from common import run, C, save

    H, r_b, r_pi, r_po, shank, k_g, k_s, rhoc, T_g = 150.0, 0.075, 0.0131, 0.016, 0.06, 1.5, 2.5, 2.3e6, 12.0
    k_p, m, T_in, cp = 0.4, 0.3, 2.0, 4190.0
    years = 5.0
    d = dict(name="V8 utube", fluid="iapws",
             well=dict(depth_m=H, mode="utube", casing=[dict(id_m=2 * r_b, layers=[])],
                       utube=dict(pipe_id_m=2 * r_pi, pipe_od_m=2 * r_po, shank_spacing_m=shank, k_pipe=k_p, k_grout=k_g)),
             rock=dict(layers=[dict(k=k_s, rho=2300.0, cp=rhoc / 2300.0)], surface_T_C=T_g, geotherm=dict(gradient_C_per_km=0.0)),
             operation=dict(pressure_bc=dict(outlet_MPa=0.3), schedule=[dict(duration_days=years * 365.25, mdot_kg_s=m, T_in_C=T_in)]),
             numerics=dict(dz_m=5.0, dt_initial_s=60.0, dt_max_days=10.0, output_times_days=[1, 30, 365, years * 365.25]),
             options=dict(r_far_m=200.0, gravity_work=False))
    sim, res = run(d)
    ts = res.ts
    t = ts["t"]; Q_ext = ts["Q_wellhead"]           # W extracted from the ground (= m (h_out - h_in), fluid warms up)
    print(f"aharrsim: N cells {sim.Nw}, steps {t.size}, Q extracted 1 d {Q_ext[np.argmin(abs(t-86400))]:.0f} W, 1 yr {Q_ext[np.argmin(abs(t-3.156e7))]:.0f} W, end {Q_ext[-1]:.0f} W; T_out end {C(ts['T_out'][-1]):.3f} C")

    # ---- pygfunction: g-function of the single borehole, then T_b by temporal superposition of the aharrsim load,
    #      then the inlet temperature the pipe model needs to extract that load from a borehole at T_b
    alpha = k_s / rhoc
    bh = gt.boreholes.Borehole(H=H, D=0.5, r_b=r_b, x=0.0, y=0.0)
    # time grid: piecewise-constant load on the aharrsim time steps (thinned to ~400 steps)
    idx = np.unique(np.round(np.geomspace(1, t.size - 1, 400)).astype(int))
    idx = np.concatenate([[0], idx])
    tq = t[idx]; Qq = Q_ext[idx]
    t_g = np.concatenate([[1e-3], tq[1:]])
    g = gt.gfunction.gFunction(bh, alpha, time=tq[1:], options={"nSegments": 8}, method="similarities", boundary_condition="UBWT").gFunc
    def g_at(tau):
        tau = np.asarray(tau, float)
        return np.where(tau <= 0, 0.0, np.interp(np.log(np.maximum(tau, 1e-9)), np.log(tq[1:]), g))
    # superposition: load steps dQ_i at times tq[i]; T_b(t) = T_g - sum dQ_i g(t - t_i) / (2 pi k H)
    dQ = np.diff(np.concatenate([[0.0], Qq]))
    pipe_R = gt.pipes.conduction_thermal_resistance_circular_pipe(r_pi, r_po, k_p)
    rows = []
    for tt in (86400.0, 30 * 86400.0, 365.25 * 86400.0, years * 365.25 * 86400.0):
        i = int(np.argmin(np.abs(t - tt)))
        Tb = T_g - np.sum(dQ[tq <= t[i]] * g_at(t[i] - tq[tq <= t[i]])) / (2 * math.pi * k_s * H)
        Tf_mean = 0.5 * (C(ts["T_in"][i]) + C(ts["T_out"][i]))
        # fluid properties at the mean temperature
        import CoolProp.CoolProp as CP
        Tm = Tf_mean + 273.15
        mu, rho, kf, cpf = (CP.PropsSI(k_, "T", Tm, "P", 3e5, "Water") for k_ in ("V", "D", "L", "C"))
        h_f = gt.pipes.convective_heat_transfer_coefficient_circular_pipe(m, r_pi, mu, rho, kf, cpf, 1e-6)
        R_fp = pipe_R + 1.0 / (h_f * 2 * math.pi * r_pi)
        ut = gt.pipes.SingleUTube([(-shank / 2, 0.0), (shank / 2, 0.0)], r_pi, r_po, bh, k_s, k_g, R_fp, J=3)
        # sign convention check: pygfunction Q_f is the heat extracted BY the fluid from the borehole (positive -> fluid warms)
        T_in_pg = ut.get_inlet_temperature(Q_ext[i], Tb, m, cpf)
        T_out_pg = ut.get_outlet_temperature(T_in_pg, Tb, m, cpf)
        # aharrsim's own borehole-wall temperature (mean over depth of the wall node)
        prof = next((p for p in res.profiles if abs(p["t"] - t[i]) < 1.0), None)
        Tw_cx = float(np.mean(prof["T_wall"]) - 273.15) if prof else float("nan")
        rows.append(dict(t_days=t[i] / 86400, Q_W=float(Q_ext[i]), T_b_pyg=float(Tb), T_b_aharrsim=Tw_cx, T_in_pyg=float(T_in_pg),
                         T_in_aharrsim=float(C(ts["T_in"][i])), T_out_pyg=float(T_out_pg), T_out_aharrsim=float(C(ts["T_out"][i])),
                         Rb_eff_pyg=float(ut.effective_borehole_thermal_resistance(m, cpf))))
        print(f"t {t[i]/86400:8.1f} d  Q {Q_ext[i]:6.0f} W | T_b pyg {Tb:6.3f} aharrsim {Tw_cx:6.3f} | T_in pyg {T_in_pg:6.3f} (aharrsim {C(ts['T_in'][i]):.3f}) "
              f"| T_out pyg {T_out_pg:6.3f} aharrsim {C(ts['T_out'][i]):6.3f}  (R_b* {rows[-1]['Rb_eff_pyg']:.4f})")
    save("V8_utube_pygfunction", dict(setup=dict(H=H, r_b=r_b, r_pipe_in=r_pi, r_pipe_out=r_po, shank=shank, k_grout=k_g, k_soil=k_s, rhoc=rhoc,
                                                   T_ground=T_g, mdot=m, T_in=T_in, years=years), rows=rows))


def b4_uloop():
    import sys, time, h5py, math
    import numpy as np
    from benchmarks import DB
    from common import save, C
    from aharrsim.systems import uloop

    class UtubeDB:
        def __init__(self, path=DB):
            f = h5py.File(path, "r")
            base = "/utube/H2O/"
            self.inp = {k: f[base + "input/" + k][:] for k in ("mdot", "L2", "L1", "grad", "D", "T_i", "k_rock", "time")}
            shape = (26, 20, 9, 5, 3, 3, 3, 161)
            We = f[base + "output/We"][:]
            valid = np.argwhere(np.isfinite(We.flatten()))[:, 0]
            U = f[base + "output/Tout/U"][:]; s = f[base + "output/Tout/sigma"][:]; Vt = f[base + "output/Tout/Vt"][:]
            Mk = U @ np.diag(s) @ Vt
            full = np.full((shape[-1], int(np.prod(shape[:-1]))), np.nan, dtype=np.float32)
            full[:, valid] = Mk
            self.Tout = full.T.reshape(shape)
            f.close()
        def idx(self, name, value):
            arr = self.inp[name]; i = int(np.argmin(np.abs(arr - value))); assert abs(arr[i] - value) < 1e-3 * max(1, abs(value)), (name, value); return i
        def series(self, mdot, L2, L1, grad, D, Ti_K, k):
            ii = (self.idx("mdot", mdot), self.idx("L2", L2), self.idx("L1", L1), self.idx("grad", grad), self.idx("D", D), self.idx("T_i", Ti_K), self.idx("k_rock", k))
            return self.inp["time"], self.Tout[ii]

    db = UtubeDB()
    cases = [dict(L1=3500.0, L2=5000.0, m=24.0, grad=0.05, D=0.33, Ti=303.15, k=3.0),
             dict(L1=2000.0, L2=2000.0, m=12.6, grad=0.05, D=0.33, Ti=303.15, k=3.0),
             dict(L1=5000.0, L2=10000.0, m=50.6, grad=0.05, D=0.33, Ti=303.15, k=3.0),
             dict(L1=5000.0, L2=20000.0, m=100.0, grad=0.07, D=0.444, Ti=333.15, k=4.5),
             dict(L1=3500.0, L2=5000.0, m=5.0, grad=0.03, D=0.216, Ti=303.15, k=1.5),
             dict(L1=4500.0, L2=12000.0, m=35.4, grad=0.06, D=0.33, Ti=318.15, k=3.0)]
    rows = []
    for cs in cases:
        t0 = time.time()
        ch = uloop("GC utube", cs["L1"], cs["L2"], cs["D"], cs["m"], cs["Ti"] - 273.15, 20e6,
                   [dict(k=cs["k"], rho=2750.0, cp=790.0)], 25.0, dict(gradient_C_per_km=cs["grad"] * 1000), 40.0,
                   roughness=2.5e-5, dz=50.0 if cs["L2"] <= 10000 else 100.0)
        tdb, Tdb = db.series(cs["m"], cs["L2"], cs["L1"], cs["grad"], cs["D"], cs["Ti"], cs["k"])
        ok = np.isfinite(Tdb)
        comp = []
        for y in (1, 2, 5, 10, 20, 30, 40):
            if y <= tdb[ok].max():
                Tm = ch["T_wh"](y * 365.25 * 86400) - 273.15
                Td = float(np.interp(y, tdb[ok], Tdb[ok])) - 273.15
                comp.append((y, Tm, Td))
        err = [a - b for _, a, b in comp]
        flags = [r.meta.get("unconverged_steps", 0) for r in (ch["injector"][1], ch["lateral"][1], ch["producer"][1])]
        print(f"L1 {cs['L1']:.0f} L2 {cs['L2']:.0f} m {cs['m']} grad {cs['grad']} D {cs['D']} k {cs['k']}: " +
              " ".join(f"{y}y {a:.1f}/{b:.1f}" for y, a, b in comp) + f"  | max|err| {max(abs(e) for e in err):.2f} K  unconv {flags}  ({time.time()-t0:.0f}s)")
        rows.append(dict(cs, comparison=comp, max_abs_err=float(max(abs(e) for e in err)), mean_err=float(np.mean(err)), unconverged=flags))
    save("B4_geocluster_utube", rows)



def v10_reservoir():
    """Reservoir elements: Laplace inversion vs the erfc single-fracture limit and mpmath; doublet breakthrough time and
    loss-free limit; cyclic fracture model vs Lauwerier in continuous injection; heat store energy balance."""
    import math
    import numpy as np
    from scipy.special import erfc
    from common import save
    from aharrsim.reservoir import ParallelFractures, DoubletAquifer, CyclicFracture, HeatStore
    YR = 365.25 * 86400
    out = {}
    pf = ParallelFractures(Q=0.05, N=1, L_flow=500.0, B=500.0, S=1e6, k_r=2.8, rhoc_r=2.4e6, T_rock_C=200, T_inj_C=60)
    t = np.array([0.5, 1, 2, 5, 10, 20, 30]) * YR
    TD = pf.T_D(t); ex = erfc(1 / (2 * np.sqrt(t / pf.gamma2)))
    out["single_fracture_vs_erfc_max_err"] = float(np.max(np.abs(TD - ex)))
    pf2 = ParallelFractures(Q=0.05, N=20, L_flow=500.0, B=300.0, S=50.0, k_r=2.8, rhoc_r=2.4e6, T_rock_C=200, T_inj_C=60)
    TD2 = pf2.T_D(t)
    try:
        import mpmath as mp
        mp.mp.dps = 30
        f = lambda s_: mp.exp(-mp.sqrt(s_) * mp.tanh(pf2.x_ED * mp.sqrt(s_))) / s_
        mpv = np.array([float(mp.invertlaplace(f, float(x), method="dehoog")) for x in t / pf2.gamma2])
        out["parallel_fractures_vs_mpmath_max_diff"] = float(np.max(np.abs(TD2 - mpv)))
        # sharp front: fractures 2.5 m apart (far closer than the conduction length), the block flushed once in 35 years.
        # de Hoog (the model) and Gaver-Stehfest N = 14 (used up to version 2.1) against 30-digit mpmath
        pf3 = ParallelFractures(Q=0.063, N=390, L_flow=880.0, B=143.0, S=2.5, k_r=2.8, rhoc_r=2.25e6, T_rock_C=190, T_inj_C=60, rhoc_w=3.88e6)
        t3 = np.array([10, 15, 20, 25, 28, 30, 32, 35, 40, 50]) * YR
        f3 = lambda s_: mp.exp(-mp.sqrt(s_) * mp.tanh(pf3.x_ED * mp.sqrt(s_))) / s_
        mp3 = np.array([float(mp.invertlaplace(f3, float(x), method="dehoog")) for x in t3 / pf3.gamma2])
        x3, y3 = pf3.groups(t3)
        out["sharp_front"] = dict(years=(t3 / YR).tolist(), x=x3.tolist(), y=y3.tolist(), mpmath=mp3.tolist(), dehoog=pf3.T_D(t3).tolist(),
                                  stehfest=pf3.T_D(t3, method="stehfest").tolist(),
                                  dehoog_max_err=float(np.max(np.abs(pf3.T_D(t3) - mp3))),
                                  stehfest_max_err=float(np.max(np.abs(pf3.T_D(t3, method="stehfest") - mp3))))
        # the two limits of the solution: isolated fractures (erfc) and thin slabs (a front at t_sweep with spread sqrt(2 t_s t_d / 3))
        xs = np.array([0.5, 0.8, 0.9, 1.0, 1.1, 1.2, 1.5]); ts_ = pf3.t_sweep
        td_ = (0.5 * pf3.S) ** 2 / (pf3.k_r / pf3.rhoc_r)
        front = 0.5 * erfc((ts_ - xs * ts_) / (np.sqrt(2.0) * np.sqrt(2.0 * xs * ts_ * td_ / 3.0)))
        out["thin_slab_limit"] = dict(x=xs.tolist(), model=pf3.T_D(xs * ts_).tolist(), front=front.tolist(), max_err=float(np.max(np.abs(pf3.T_D(xs * ts_) - front))))
    except ImportError:
        out["parallel_fractures_vs_mpmath_max_diff"] = None
    da = DoubletAquifer(Q=0.05, D=1000.0, H=50.0, phi=0.15, rhoc_r=2.2e6, k_cap=2.0, T_aq_C=70, T_inj_C=40)
    tb = math.pi * 50 * 1000 ** 2 * da.rhoc_a / (3 * 0.05 * da.rhoc_w)
    da0 = DoubletAquifer(Q=0.05, D=1000.0, H=50.0, phi=0.15, rhoc_r=2.2e6, k_cap=1e-12, T_aq_C=70, T_inj_C=40)
    chk = []
    for yrs in (20, 30, 40, 60, 100):
        tt = yrs * YR
        chk.append(dict(years=yrs, TD_no_loss=float(da0.T_D(np.array([tt]))[0]), arrived_fraction=float(np.mean(da0.tau < tt)),
                        TD_with_caprock=float(da.T_D(np.array([tt]))[0]), T_prod_C=da(tt) - 273.15))
    out["doublet"] = dict(t_breakthrough_model_yr=da.t_breakthrough / YR, t_breakthrough_formula_yr=tb / YR, rows=chk)
    cf = CyclicFracture(L_f=200.0, B=100.0, w=0.005, k_r=2.8, rhoc_r=2.4e6, T_rock_C=200, T_inj_C=50, n_x=80, n_z=80, z_max=30.0)
    q = 0.01
    tau = cf.L * cf.B * cf.w / q
    beta = 2 * math.sqrt(2.8 * 2.4e6) / (cf.rhoc_w * cf.w)
    ts = np.array([2, 5, 10, 20, 40, 80]) * tau
    Tm = cf.simulate_continuous(q, ts[-1], ts) - 273.15
    TD_l = np.where(ts > tau, erfc(beta * tau / (2 * np.sqrt(np.maximum(ts - tau, 1e-30)))), 0.0)
    T_l = 200 - TD_l * (200 - 50)
    out["cyclic_vs_lauwerier"] = dict(t_over_tau=[2, 5, 10, 20, 40, 80], model_C=Tm.tolist(), lauwerier_C=T_l.tolist(),
                                      max_err_K=float(np.max(np.abs(Tm - T_l))))
    hs = HeatStore(V=1e6, A_b=6e4, rhoc_store=2.6e6, k_r=1.5, rhoc_r=2.4e6, T_store_C=200, T0_C=20)
    r = hs.run(15, 20.0, 60.0)
    E0 = 1e6 * 2.6e6 * 180
    Eout = np.sum(r["Q_MW"] * 1e6 * 10 * 86400); Eloss = np.sum(r["loss_MW"] * 1e6 * 10 * 86400)
    Eleft = 1e6 * 2.6e6 * (r["T_C"][-1] - 20)
    out["heat_store"] = dict(stored_GWh=E0 / 3.6e12, produced_GWh=Eout / 3.6e12, lost_GWh=Eloss / 3.6e12, remaining_GWh=Eleft / 3.6e12,
                             balance_error=float((E0 - Eout - Eloss - Eleft) / E0))
    for k, v in out.items():
        print("V10", k, v if not isinstance(v, dict) or len(str(v)) < 300 else {kk: vv for kk, vv in v.items() if kk != "rows"})
    save("V10_reservoir", out)


def v11_plant():
    """Plant models against the GEOPHIRES v3 utilisation-efficiency correlations (Aspen-derived)."""
    import numpy as np
    import CoolProp.CoolProp as CP
    from common import save
    from aharrsim import plant
    m, T_amb = 100.0, 15.0
    rows = []
    for T in (100, 120, 150, 175, 200, 225, 250):
        p = max(5e6, 1.3 * CP.PropsSI("P", "T", T + 273.15, "Q", 0, "Water"))
        h = CP.PropsSI("H", "T", T + 273.15, "P", p, "Water")
        Treinj = float(plant.geophires_reinjection_T("orc", T, T_amb))
        r = plant.orc(m, h, p, T_cond_C=30.0, T_b_min_C=Treinj, cooling="air", dT_ev=2.0)
        av = plant.geophires_availability(T_amb, T, T_amb)
        rows.append(dict(T_in_C=T, kind="binary", eta_u_cycle=r["W_net"] / (m * av), eta_u_geophires_sub=float(plant.geophires_eta_u("orc", T, T_amb)),
                         eta_u_geophires_super=float(plant.geophires_eta_u("sorc", T, T_amb)), fluid=r["fluid"], recuperator=r["recuperator"],
                         W_per_kg_kJ=r["W_net"] / m / 1e3))
    for T in (150, 175, 200, 225, 250, 275, 300):
        p_in = CP.PropsSI("P", "T", T + 273.15, "Q", 0, "Water") * 1.001
        h = CP.PropsSI("H", "T", T + 273.15, "P", p_in, "Water")
        av = plant.geophires_availability(T_amb, T, T_amb)
        for Tc in (35.0, 40.0):
            s1 = plant.flash(m, h, p_in, T_cond_C=Tc); s2 = plant.flash(m, h, p_in, T_cond_C=Tc, double=True)
            rows.append(dict(T_in_C=T, kind=f"flash, condenser {Tc:.0f} C", eta_u_single=s1["W_net"] / (m * av), eta_u_double=s2["W_net"] / (m * av),
                             eta_u_geophires_single=float(plant.geophires_eta_u("flash1", T, T_amb)),
                             eta_u_geophires_double=float(plant.geophires_eta_u("flash2", T, T_amb)), p_sep_bar=s1["p_sep"] / 1e5,
                             W_single_per_kg_kJ=s1["W_net"] / m / 1e3, W_double_per_kg_kJ=s2["W_net"] / m / 1e3))
    for r in rows:
        print("V11", {k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items()})
    save("V11_plant", rows)


if __name__ == "__main__":
    which = sys.argv[1:] or ["v9", "v9b", "v8", "b4", "v10", "v11"]
    for w in which:
        dict(v9=v9_single_well, v9b=v9b_twophase, v8=v8_utube, b4=b4_uloop, v10=v10_reservoir, v11=v11_plant)[w]()
