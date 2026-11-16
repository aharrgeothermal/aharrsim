"""One representative example of each of the 20 development types in geothermal_development_types.xlsx
(Khodayar & Bjornsson 2024), run with the extended aharrsim: well models (coaxial loop, single-channel
production / injection wells, U-tube), reservoir elements (aharrsim.reservoir) and surface plants
(aharrsim.plant).  Reports the produced water temperature and the electricity that can be made from it.

Reference conditions for every case: ambient 15 C; binary plant = GEOPHIRES practical utilisation
(air-cooled, brine returned at the GEOPHIRES reinjection temperature), cycle model as a check; flash and
dry-steam plants from the thermodynamic models with a 40 C condenser; well pumps 70 % efficient; EGS
circulation impedance applied as Q x impedance.  Rock k and rho c as stated per case.

    python dev_types.py            # all types (~15 min on 2 cores)      -> results/DT_types.json
    python dev_types.py 1 7 20     # selected ids
"""
from __future__ import annotations

import math
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import CoolProp.CoolProp as CP

from common import save, load
from aharrsim.systems import (single_well, coaxial_well, doublet, uloop, run_case, outlet_series, Series, YEAR, at_years)
from aharrsim.reservoir import ConstantReservoir, DoubletAquifer, ParallelFractures, CyclicFracture, HeatStore, cbar_w
from aharrsim import plant

IN = 0.0254
T_AMB = 15.0
STEEL = dict(k=45.0, rho=7850.0, cp=490.0)
CEMENT = dict(k=1.0, rho=1900.0, cp=1000.0)
YRS = (1, 5, 10, 20, 30)


def cased(id_in, od_in, hole_in):
    return dict(id_in=id_in, layers=[dict(od_in=od_in, material="steel"), dict(od_in=hole_in, material="cement")])


def _T(h, p):
    return CP.PropsSI("T", "H", h, "P", p, "Water") - 273.15


def _hT(T_C, p):
    return CP.PropsSI("H", "T", T_C + 273.15, "P", p, "Water")


def plant_options(m, h, p, T_return_C=None):
    """Electricity from a stream m kg/s at (h, p) under the reference conditions; all values in MW."""
    T = _T(h, p)
    x = CP.PropsSI("Q", "H", h, "P", p, "Water")
    out = dict(T_C=T, x=float(x) if 0 <= x <= 1 else (1.0 if h > CP.PropsSI("H", "P", p, "Q", 1, "Water") else 0.0),
               exergy_MW=plant.exergy(m, h, p, T_AMB) / 1e6)
    pr = plant.orc_practical(m, h, p, T_AMB)
    out.update(binary_MWe=pr["W_net"] / 1e6, binary_T_return=pr["T_b_out"], binary_eta_u=pr["eta_u"])
    try:
        cyc = plant.orc(m, h, p, T_cond_C=T_AMB + 15.0, T_b_min_C=max(40.0, pr["T_b_out"] if T > 150 else 40.0), cooling="air", dT_ev=3.0)
        out["binary_cycle_MWe"] = cyc["W_net"] / 1e6
        out["binary_cycle_fluid"] = cyc["fluid"]
    except Exception:
        out["binary_cycle_MWe"] = float("nan")
    if T > 140.0 or out["x"] > 0:
        try:
            f1 = plant.flash(m, h, p, T_cond_C=40.0)
            f2 = plant.flash(m, h, p, T_cond_C=40.0, double=True)
            out.update(flash1_MWe=f1["W_net"] / 1e6, flash1_p_sep_bar=(f1["p_sep"] or 0) / 1e5, flash1_T_brine=f1["T_brine_out"],
                       flash2_MWe=f2["W_net"] / 1e6, flash2_p_sep_bar=(f2["p_sep"] or 0) / 1e5)
        except Exception:
            pass
    if out["x"] >= 1.0:
        ds = plant.dry_steam(m, h, p, T_cond_C=40.0)
        out["dry_steam_MWe"] = ds["W_net"] / 1e6
    if T_return_C is not None:
        out["heat_MWth"] = m * (h - _hT(T_return_C, p)) / 1e6
    return out


def ts_summary(res, years=YRS):
    ts = res.ts
    T = Series(ts["t"], ts["T_out"]); P = Series(ts["t"], ts["p_out"]); H = Series(ts["t"], ts["h_out"])
    out = {}
    for y in years:
        t = y * YEAR
        if t <= ts["t"][-1] + 1:
            out[y] = dict(T_C=T(t) - 273.15, p_MPa=P(t) / 1e6, h=H(t),
                          pump_kW=float(np.interp(t, ts["t"], ts["pump_kW"])), margin_MPa=float(np.interp(t, ts["t"], ts["inflow_margin"])) / 1e6,
                          x_out=float(np.interp(t, ts["t"], ts["x_out"])))
    return out


def pipe_friction_Pa(m, d, L, T_C, p, rough=4.6e-5):
    """Single-phase friction pressure drop [Pa] of m kg/s through L metres of pipe of diameter d (Churchill)."""
    from aharrsim.correlations import friction_churchill
    rho = CP.PropsSI("D", "T", T_C + 273.15, "P", p, "Water"); mu = CP.PropsSI("V", "T", T_C + 273.15, "P", p, "Water")
    v = m / (rho * math.pi * d * d / 4.0)
    f = float(friction_churchill(rho * v * d / mu, rough / d))
    return f * L / d * 0.5 * rho * v * v


def loop_pump(ch, m, imp_Pa_per_m3s, T_inj_C, years, rho_res=900.0, dp_plant=0.3e6, eta=0.7):
    """Circulation power of an injector - fractures - producer loop from the closed-loop pressure balance on the simulated columns:
        pump dp = I Q + (BHP - WHP)_producer - (BHP - WHP)_injector + dp_plant,   W = (m / rho_surface) max(dp, 0) / eta.
    The cold, dense injection column pushes and the hot, light production column lifts (thermosiphon), so the pumps supply
    only what the fracture impedance, the wellbore friction and the plant take beyond that.  It does not matter for the total
    whether the pump sits in the producer or at the injection wellhead.  Returns {year: dict}."""
    ti, tp = ch["injector"][1].ts, ch["producer"][1].ts
    rho_s = CP.PropsSI("D", "T", T_inj_C + 273.15, "P", 2.0e6, "Water")
    out = {}
    for y in years:
        t = y * YEAR
        col_i = float(np.interp(t, ti["t"], ti["p_out"] - ti["p_in"]))
        col_p = float(np.interp(t, tp["t"], tp["p_in"] - tp["p_out"]))
        dp_res = imp_Pa_per_m3s * m / rho_res
        dp = dp_res + col_p - col_i + dp_plant
        out[y] = dict(kW=max(dp, 0.0) * m / rho_s / eta / 1e3, dp_MPa=dp / 1e6, dp_res_MPa=dp_res / 1e6, col_inj_MPa=col_i / 1e6, col_prod_MPa=col_p / 1e6,
                      thermosiphon_MPa=(col_i - col_p) / 1e6)
    return out


def sealed_loop_pump(res, m, T_in_C, years, dp_plant=0.3e6, eta=0.7):
    """Circulation power of a sealed loop solved as one well (coaxial or U-tube): the balance of loop_pump without a reservoir term,
        pump dp = p_inlet - p_outlet + dp_plant,   W = (m / rho_surface) max(dp, 0) / eta,
    with the inlet pressure that the pressure march needs to deliver the outlet pressure, so that the thermosiphon is credited in full.
    Returns {year: dict}."""
    ts = res.ts
    rho_s = CP.PropsSI("D", "T", T_in_C + 273.15, "P", 2.0e6, "Water")
    out = {}
    for y in years:
        t = y * YEAR
        if t <= ts["t"][-1] + 1:
            loss = float(np.interp(t, ts["t"], ts["p_in"] - ts["p_out"]))
            dp = loss + dp_plant
            out[y] = dict(kW=max(dp, 0.0) * m / rho_s / eta / 1e3, dp_MPa=dp / 1e6, loop_loss_MPa=loss / 1e6)
    return out


def uloop_pump(s, m, T_in_C, p_in=2.0e6, dp_plant=0.3e6, eta=0.7):
    """The same balance for an injector - laterals - producer loop run as a chain from a fixed injection pressure: the loop loses
    p_in less the pressure that arrives at the production wellhead (s: the summary of the producer). Returns {year: dict}."""
    rho_s = CP.PropsSI("D", "T", T_in_C + 273.15, "P", 2.0e6, "Water")
    return {y: dict(dp_MPa=p_in / 1e6 - s[y]["p_MPa"] + dp_plant / 1e6, kW=max(p_in / 1e6 - s[y]["p_MPa"] + dp_plant / 1e6, 0.0) * 1e6 * m / rho_s / eta / 1e3,
                    loop_loss_MPa=p_in / 1e6 - s[y]["p_MPa"]) for y in s}


def flags(*results):
    f = dict(unconverged=0, choke=0, max_residual_kW=0.0, feed_floor=0, steps=0)
    for r in results:
        ts = r.ts
        f["unconverged"] += int(r.meta.get("unconverged_steps", 0))
        f["choke"] += int(np.sum(ts["choke_warning"] > 0))
        f["feed_floor"] += int(np.sum(ts["feed_floor"] > 0))          # steps on which a temperature-specified feed was evaluated at 1.02 p_sat
        f["steps"] += int(len(ts["t"]))
        f["max_residual_kW"] = max(f["max_residual_kW"], float(np.max(np.abs(ts["Q_wellhead"] - ts["Q_rock"] - ts["Q_res"]
                                                                              + ts["dE_fluid"] + ts["Q_ke"] + ts["Q_pe"])) / 1e3))
    return f


def cheap_power(m, h, p, kind):
    """Gross MWe from the fast plant models (no ORC cycle search) for the time-series charts."""
    if kind is None or m <= 0:
        return 0.0
    T = _T(h, p)
    xq = CP.PropsSI("Q", "H", h, "P", p, "Water")
    x = float(xq) if 0 <= xq <= 1 else (1.0 if h > CP.PropsSI("H", "P", p, "Q", 1, "Water") else 0.0)
    w = {}
    if kind in ("binary", "best"):
        try:
            w["binary"] = plant.orc_practical(m, h, p, T_AMB)["W_net"] / 1e6 if T > 60 else 0.0
        except Exception:
            w["binary"] = 0.0
    if kind in ("flash2", "steam", "best") and (T > 140 or x > 0):
        try:
            w["flash1"] = plant.flash(m, h, p, T_cond_C=40.0)["W_net"] / 1e6
            w["flash2"] = plant.flash(m, h, p, T_cond_C=40.0, double=True)["W_net"] / 1e6
        except Exception:
            pass
    if kind in ("steam", "best") and x >= 1.0:
        w["steam"] = plant.dry_steam(m, h, p, T_cond_C=40.0)["W_net"] / 1e6
    return max(w.values()) if w else 0.0


def ts_series(res, m, T_return_C=None, kind="binary", n=70, t_min_days=2.0, m_plant=None):
    """Decimated (log-spaced) time series of the produced stream for charts, with heat and gross electricity."""
    ts = res.ts
    t = ts["t"]
    tg = np.unique(np.concatenate([np.geomspace(max(t_min_days * 86400.0, t[1]), t[-1], n), [t[-1]]]))
    mp = m if m_plant is None else m_plant
    out = dict(t_yr=[], T_C=[], p_MPa=[], x=[], pump_kW=[], Q_wh_MW=[], heat_MWth=[], MWe=[])
    for tt in tg:
        h = float(np.interp(tt, t, ts["h_out"])); pp = float(np.interp(tt, t, ts["p_out"]))
        out["t_yr"].append(tt / YEAR); out["T_C"].append(float(np.interp(tt, t, ts["T_out"])) - 273.15)
        out["p_MPa"].append(pp / 1e6); out["x"].append(float(np.interp(tt, t, ts["x_out"])))
        out["pump_kW"].append(float(np.interp(tt, t, ts["pump_kW"]))); out["Q_wh_MW"].append(float(np.interp(tt, t, ts["Q_wellhead"])) / 1e6)
        out["heat_MWth"].append(mp * (h - _hT(T_return_C, pp)) / 1e6 if T_return_C is not None else None)
        out["MWe"].append(cheap_power(mp, h, pp, kind))
    return out


def last_profile(res, n=160, idx=-1):
    """Depth profile at an output time (default the last), decimated."""
    if not res.profiles:
        return None
    pr = res.profiles[idx]
    N = len(pr["md"])
    idx = np.unique(np.linspace(0, N - 1, min(N, n)).astype(int))
    out = dict(t_yr=float(pr["t"]) / YEAR)
    for k in ("md", "tvd", "T0", "T_a", "T_t", "p_a", "p_t", "x_t", "x_a"):
        if k in pr:
            v = np.asarray(pr[k], float)
            out[k] = (v[idx] - (273.15 if k in ("T0", "T_a", "T_t") else 0.0)).tolist() if len(v) == N else None
    if "q_rock" in pr:
        q = np.asarray(pr["q_rock"], float); Nc = len(q)
        ic = np.unique(np.linspace(0, Nc - 1, min(Nc, n)).astype(int))
        out["md_c"] = np.asarray(pr["md_c"], float)[ic].tolist(); out["q_rock_W_m"] = q[ic].tolist()
    return out


def fn_series(fn, t_end, n=70, t_min_days=2.0, offset_K=273.15):
    """Sample a callable of time (K) on the chart grid."""
    tg = np.unique(np.concatenate([np.geomspace(t_min_days * 86400.0, t_end, n), [t_end]]))
    return dict(t_yr=(tg / YEAR).tolist(), T_C=[float(fn(tt)) - offset_K for tt in tg])


# =============================================================================================
def t01_lt_hsa():
    """Paris Basin Dogger doublet: 1.8 km, 70 C aquifer, 68 kg/s, 1.2 km spacing, 40 C return (district heating)."""
    depth, m, T_aq, T_inj = 1800.0, 68.0, 70.0, 40.0
    rock = [dict(k=2.2, rho=2400.0, cp=900.0)]
    geo = dict(gradient_C_per_km=(T_aq - 12.0) / depth * 1000)
    Q = m / 985.0
    fac = lambda T_bot: DoubletAquifer(Q=Q, D=1200.0, H=20.0, phi=0.15, rhoc_r=2.2e6, k_cap=2.0, T_aq_C=T_aq,
                                        T_inj_C=lambda t: T_bot(t) - 273.15, k_aq=2e-12, mu=4e-4)
    ch = doublet("LT doublet", fac, depth, 8.535 * IN, m, T_inj, rock, 12.0, geo, 30.0,
                 layers=[dict(od_in=9.625, material="steel"), dict(od_in=12.25, material="cement")],
                 p_inj=0.8e6, p_wh=0.5e6, pump=dict(md_m=250.0, efficiency=0.7), p_res=17.5e6, PI=150e-6)
    res_p = ch["producer"][1]
    s = ts_summary(res_p)
    po = {y: plant_options(m, s[y]["h"], s[y]["p_MPa"] * 1e6, T_return_C=T_inj) for y in s}
    dp_res = ch["reservoir"].pumping_dp()
    series = ts_series(res_p, m, T_inj, "binary"); series["T_res_C"] = fn_series(ch["T_res"], 30 * YEAR)["T_C"]
    return dict(id=1, name="Low-temperature / hot sedimentary aquifer doublet (Paris Basin type)", series=series, profile=last_profile(res_p),
                assumptions=dict(depth_m=depth, T_aquifer_C=T_aq, H_m=20, phi=0.15, spacing_m=1200, mdot=m, T_return_C=T_inj,
                                 t_breakthrough_yr=ch["reservoir"].t_breakthrough / YEAR, reservoir_dp_MPa=dp_res / 1e6),
                T_out=s, T_res={y: ch["T_res"](y * YEAR) - 273.15 for y in s}, plant=po, well_pump_kW={y: s[y]["pump_kW"] for y in s},
                flags=flags(ch["injector"][1], res_p), plant_kind="direct use; binary marginal", profile_inj=last_profile(ch["injector"][1]))


def t02_mt_binary():
    """Nevada-type pumped binary: 2.0 km, 165 C convective reservoir (0.5 K/yr decline), 75 kg/s per well, LSP at 450 m."""
    depth, m = 2000.0, 75.0
    rock = [dict(k=2.5, rho=2600.0, cp=900.0)]
    geo = dict(gradient_C_per_km=75.0)
    resv = ConstantReservoir(165.0, decline_K_per_yr=0.5)
    d = single_well("MT pumped well", "production", depth, 8.535 * IN, m, rock, 15.0, geo, 30.0,
                    layers=[dict(od_in=9.625, material="steel"), dict(od_in=12.25, material="cement")],
                    inflow=dict(T_K=resv, p_res_MPa=18.0, PI_kg_s_per_MPa=40.0), p_wh=1.2e6, pump=dict(md_m=450.0, efficiency=0.7))
    sim, res = run_case(d)
    s = ts_summary(res)
    po = {y: plant_options(m, s[y]["h"], s[y]["p_MPa"] * 1e6, T_return_C=70.0) for y in s}
    return dict(id=2, name="Medium-temperature pumped well + binary plant (Nevada type)", series=ts_series(res, m, 70.0, "binary"), profile=last_profile(res),
                assumptions=dict(depth_m=depth, T_res_C=165, decline_K_per_yr=0.5, mdot=m, pump_md_m=450, p_wh_MPa=1.2),
                T_out=s, plant=po, well_pump_kW={y: s[y]["pump_kW"] for y in s}, flags=flags(res), plant_kind="binary (ORC)")


def t03_ht_flash():
    """Wairakei-type liquid-dominated HT: 250 C feed at 1.5 km, self-flowing two-phase, 40 kg/s, 0.9 MPa wellhead; flash plant."""
    depth, m = 1500.0, 40.0
    rock = [dict(k=2.2, rho=2500.0, cp=950.0)]
    geo = dict(gradient_C_per_km=(250.0 - 15.0) / 1.5)
    d = single_well("HT flash well", "production", depth, 8.535 * IN, m, rock, 15.0, geo, 30.0,
                    layers=[dict(od_in=9.625, material="steel"), dict(od_in=12.25, material="cement")],
                    inflow=dict(T_C=250.0, p_res_MPa=9.0, PI_kg_s_per_MPa=60.0), p_wh=0.9e6)
    sim, res = run_case(d)
    s = ts_summary(res)
    po = {y: plant_options(m, s[y]["h"], s[y]["p_MPa"] * 1e6) for y in s}
    return dict(id=3, name="High-temperature self-flowing well, flash plant (Wairakei type)", series=ts_series(res, m, 70.0, "flash2"), profile=last_profile(res),
                assumptions=dict(depth_m=depth, T_res_C=250, mdot=m, p_wh_MPa=0.9, p_res_MPa=9.0, PI=60),
                T_out=s, plant=po, flags=flags(res), plant_kind="single / double flash")


def t03b_dry_steam():
    """Geysers-type dry steam: 240 C saturated steam reservoir at 2.5 km, 3.3 MPa, 15 kg/s per well, 0.8 MPa wellhead."""
    depth, m = 2500.0, 15.0
    p_res = CP.PropsSI("P", "T", 240 + 273.15, "Q", 1, "Water")
    h_res = CP.PropsSI("H", "P", p_res, "Q", 1, "Water")
    rock = [dict(k=2.5, rho=2600.0, cp=900.0)]
    geo = dict(gradient_C_per_km=(240.0 - 15.0) / 2.5)
    d = single_well("dry steam well", "production", depth, 8.535 * IN, m, rock, 15.0, geo, 30.0,
                    layers=[dict(od_in=9.625, material="steel"), dict(od_in=12.25, material="cement")],
                    inflow=dict(h_J_kg=h_res, p_res_MPa=p_res / 1e6, PI_kg_s_per_MPa=25.0), p_wh=0.8e6)
    sim, res = run_case(d)
    s = ts_summary(res)
    po = {y: plant_options(m, s[y]["h"], s[y]["p_MPa"] * 1e6) for y in s}
    return dict(id=3.5, name="High-temperature vapour-dominated (dry steam, Geysers type)", series=ts_series(res, m, 70.0, "steam"), profile=last_profile(res),
                assumptions=dict(depth_m=depth, T_res_C=240, p_res_MPa=p_res / 1e6, mdot=m, p_wh_MPa=0.8, PI=25),
                T_out=s, plant=po, flags=flags(res), plant_kind="dry steam")


def t04_waste_heat():
    """Climeon-type modular ORC on a 110 C brine stream of 50 kg/s (e.g. separated brine), no wells."""
    m, T = 50.0, 110.0
    p = 1.0e6
    h = _hT(T, p)
    po = plant_options(m, h, p, T_return_C=70.0)
    return dict(id=4, name="Modular ORC on surface waste heat (Climeon type)", assumptions=dict(T_stream_C=T, mdot=m, p_MPa=1.0),
                T_out={0: dict(T_C=T, h=h, p_MPa=1.0, pump_kW=0.0)}, plant={0: po}, flags=dict(unconverged=0, choke=0, max_residual_kW=0.0),
                plant_kind="modular ORC")


def t05_gshp():
    """Ground-source heat pump: 150 m single U-tube, 12 C ground, 0.3 kg/s; heating season Oct-Mar at 2 C inlet for 10 years."""
    H = 150.0
    sched = []
    for yr in range(10):
        sched.append(dict(duration_days=182.0, mdot_kg_s=0.3, T_in_C=2.0, label="heating"))
        sched.append(dict(duration_days=183.25, mdot_kg_s=0.0, T_in_C=12.0, label="off"))
    d = dict(name="GSHP U-tube", fluid="iapws",
             well=dict(depth_m=H, mode="utube", casing=[dict(id_m=0.15, layers=[])],
                       utube=dict(pipe_id_m=0.0262, pipe_od_m=0.032, shank_spacing_m=0.06, k_pipe=0.4, k_grout=1.5)),
             rock=dict(layers=[dict(k=2.5, rho=2300.0, cp=1000.0)], surface_T_C=12.0, geotherm=dict(gradient_C_per_km=0.0)),
             operation=dict(pressure_bc=dict(outlet_MPa=0.3), schedule=sched),
             numerics=dict(dz_m=5.0, dt_initial_s=60.0, dt_max_days=5.0, r_far_m=120.0, output_times_days=[100, 3650]),
             options=dict(gravity_work=False))
    sim, res = run_case(d)
    ts = res.ts
    on = ts["mdot"] > 0
    # season means per year
    seasons = []
    for yr in range(10):
        a, b = yr * 365.25 * 86400, (yr * 365.25 + 182) * 86400
        sel = on & (ts["t"] > a) & (ts["t"] <= b)
        Q = ts["Q_wellhead"][sel]; To = ts["T_out"][sel] - 273.15
        seasons.append(dict(year=yr + 1, Q_mean_kW=float(np.mean(Q) / 1e3), Q_end_kW=float(Q[-1] / 1e3), T_out_mean_C=float(np.mean(To)),
                            T_out_end_C=float(To[-1])))
    last = seasons[-1]
    loss = float(np.mean((ts["p_in"] - ts["p_out"])[on]))                           # loop loss while circulating
    circ = dict(loop_loss_MPa=loss / 1e6, W=max(loss, 0.0) * 0.3 / 1000.0 / 0.7, note="borehole loop only, no heat-pump evaporator; not in the heat-pump consumption")
    hp = plant.heat_pump(last["Q_mean_kW"] * 1e3, 0.5 * (2.0 + last["T_out_mean_C"]), 45.0)
    tg = np.arange(1.0, ts["t"][-1] / 86400.0, 4.0) * 86400.0
    series = dict(t_yr=(tg / YEAR).tolist(), T_C=(np.interp(tg, ts["t"], ts["T_out"]) - 273.15).tolist(),
                  Q_kW=(np.interp(tg, ts["t"], ts["Q_wellhead"]) / 1e3).tolist(), on=(np.interp(tg, ts["t"], ts["mdot"]) > 0).astype(int).tolist())
    return dict(id=5, name="Ground-source heat pump, 150 m U-tube borehole", series=series, profile=last_profile(res, idx=0),   # day 100: mid heating season
                assumptions=dict(H_m=H, T_ground_C=12, mdot=0.3, T_in_C=2,
                heating_days=182, pipe="32 mm SDR11 single U, grout k 1.5", sink_C=45),
                seasons=seasons, heat_pump=dict(COP=hp["COP"], W_e_kW=hp["W_e"] / 1e3, Q_heat_kW=hp["Q_heat"] / 1e3, Q_ground_kW=last["Q_mean_kW"],
                                                per_metre_W=last["Q_mean_kW"] * 1e3 / H),
                T_out={10: dict(T_C=last["T_out_mean_C"], h=float("nan"), p_MPa=0.3, pump_kW=0.0)}, plant={}, flags=flags(res),
                loop_circulation=circ, plant_kind="heat pump (consumes electricity)")


def t06_shallow_coaxial():
    """Deep borehole heat exchanger 400 m coaxial (6-5/8 in casing, 90 mm PE inner pipe), 2 kg/s at 5 C inlet, 10 years, heat pump."""
    depth = 400.0
    casing = [dict(top_m=0, bottom_m=depth, id_in=5.921, layers=[dict(od_in=6.625, material="steel"), dict(od_in=8.5, material="cement")])]
    tubing = [dict(top_m=0, bottom_m=depth, id_mm=73.6, layers=[dict(od_mm=90.0, k=0.4)])]
    d = coaxial_well("DBHE 400 m", depth, casing, tubing, 2.0, 5.0, 0.3e6, [dict(k=2.6, rho=2500.0, cp=900.0)], 10.0,
                     dict(gradient_C_per_km=30.0), 10.0, dz=10.0, output_years=(1, 5, 10))
    d["numerics"]["r_far_m"] = 150.0
    sim, res = run_case(d)
    s = ts_summary(res, years=(1, 5, 10))
    Q = {y: 2.0 * (s[y]["h"] - _hT(5.0, 0.3e6)) / 1e3 for y in s}
    hp = plant.heat_pump(Q[10] * 1e3, 0.5 * (5.0 + s[10]["T_C"]), 45.0)
    lc = sealed_loop_pump(res, 2.0, 5.0, list(s), dp_plant=0.0)[10]
    circ = dict(loop_loss_MPa=lc["loop_loss_MPa"], W=lc["kW"] * 1e3, note="borehole loop only, no heat-pump evaporator; not in the heat-pump consumption")
    return dict(id=6, name="Shallow geothermal: 400 m coaxial borehole heat exchanger + heat pump", series=ts_series(res, 2.0, 5.0, None), profile=last_profile(res),
                assumptions=dict(depth_m=depth, mdot=2.0, T_in_C=5.0, surface_T_C=10, gradient_C_per_km=30),
                T_out=s, Q_ground_kW=Q, heat_pump=dict(COP=hp["COP"], W_e_kW=hp["W_e"] / 1e3, Q_heat_kW=hp["Q_heat"] / 1e3,
                                                        per_metre_W=Q[10] * 1e3 / depth),
                loop_circulation=circ, plant={}, flags=flags(res), plant_kind="heat pump (consumes electricity)")


def t07_greenloop():
    """Deep coaxial retrofit in a hot field (Coso type): 3.0 km, 80 C/km (255 C at TD), 3-1/2 x 4-1/2 VIT in 9-5/8 / 7 in casing, 10 kg/s."""
    depth = 3000.0
    casing = [dict(top_m=0, bottom_m=1500, id_in=8.535, layers=[dict(od_in=9.625, material="steel"), dict(od_in=12.25, material="cement")]),
              dict(top_m=1500, bottom_m=depth, id_in=6.094, layers=[dict(od_in=7.0, material="steel"), dict(od_in=8.5, material="cement")])]
    tubing = [dict(top_m=0, bottom_m=depth, id_in=2.992, layers=[dict(od_in=3.5, k=45), dict(od_in=4.0, k=0.03), dict(od_in=4.5, k=45)])]
    out = {}
    for m in (5.0, 10.0):
        d = coaxial_well(f"GreenLoop {m} kg/s", depth, casing, tubing, m, 40.0, 2.0e6, [dict(k=2.5, rho=2600.0, cp=900.0)], 15.0,
                         dict(gradient_C_per_km=80.0), 30.0)
        sim, res = run_case(d)
        s = ts_summary(res)
        po = {y: plant_options(m, s[y]["h"], s[y]["p_MPa"] * 1e6, T_return_C=40.0) for y in s}
        out[m] = dict(T_out=s, plant=po, flags=flags(res), series=ts_series(res, m, 40.0, "binary"), profile=last_profile(res),
                      loop_pump=sealed_loop_pump(res, m, 40.0, list(s)))
    m = 10.0
    return dict(id=7, name="Single-well deep coaxial retrofit in a hot field (GreenLoop type)", series=out[m]["series"], profile=out[m]["profile"],
                assumptions=dict(depth_m=depth, gradient_C_per_km=80, T_TD_C=255, mdot=[5, 10], T_in_C=40, VIT="3-1/2 x 4-1/2, k 0.03"),
                T_out=out[m]["T_out"], plant=out[m]["plant"], flags=out[m]["flags"], variants={str(k): v for k, v in out.items()},
                plant_kind="binary (ORC)", well_pump_kW={y: out[m]["loop_pump"][y]["kW"] for y in out[m]["loop_pump"]})


def t08_eavor():
    """Eavor-Loop 1 type: two 4.5 km verticals, 12 laterals x 5 km (60 km) at 55 C/km (263 C), 8-1/2 in open-hole laterals,
    9-5/8 in cased verticals, 100 kg/s, 40 C injection; laterals 80 m apart (unit-cell far field) and, as a bound, isolated."""
    L1, L2, n = 4500.0, 5000.0, 12
    rock = [dict(k=2.8, rho=2600.0, cp=900.0)]
    geo = dict(gradient_C_per_km=55.0)
    out = {}
    for lab, ff, rf in (("laterals 80 m apart", "adiabatic", 40.0), ("isolated laterals", "dirichlet", None)):
        ch = uloop(f"Eavor {lab}", L1, L2, 0.216, 100.0, 40.0, 2.0e6, rock, 15.0, geo, 30.0, n_lat=n,
                   layers=None, roughness=1e-4, dz=50.0, lat_far_field=ff, lat_r_far=rf)
        res_p = ch["producer"][1]
        s = ts_summary(res_p)
        po = {y: plant_options(100.0, s[y]["h"], s[y]["p_MPa"] * 1e6, T_return_C=40.0) for y in s}
        ts_i = ch["injector"][1].ts
        ser = ts_series(res_p, 100.0, 40.0, "binary"); ser["T_toe_C"] = fn_series(ch["T_toe"], 30 * YEAR)["T_C"]
        # closed loop: the circulation pump supplies what the loop loses between the two wellheads plus the plant (Eq. 6 without a reservoir term)
        lp = uloop_pump(s, 100.0, 40.0)
        out[lab] = dict(T_out=s, plant=po, flags=flags(ch["injector"][1], ch["lateral"][1], res_p),
                        T_toe={y: ch["T_toe"](y * YEAR) - 273.15 for y in s}, p_wh_MPa={y: s[y]["p_MPa"] for y in s}, loop_pump=lp,
                        thermosiphon_note="injection at 2 MPa; production wellhead pressure from the march; pump = (2 MPa - wellhead + 0.3 MPa plant) m / (rho 0.7)", series=ser,
                        profile=last_profile(res_p), profile_lat=last_profile(ch["lateral"][1]), profile_inj=last_profile(ch["injector"][1]))
    main = out["laterals 80 m apart"]
    return dict(id=8, name="Multilateral deep closed loop (Eavor-Loop 1 type)", series=main["series"], profile=main["profile"],
                profile_lat=main["profile_lat"], profile_inj=main["profile_inj"],
                assumptions=dict(L1_m=L1, L2_m=L2, n_laterals=n, total_km=(2 * L1 + n * L2) / 1000, gradient_C_per_km=55, T_rock_C=263, mdot=100,
                                 T_in_C=40, lateral_id_m=0.216, lateral_spacing_m=80),
                T_out=main["T_out"], plant=main["plant"], flags=main["flags"], variants=out, plant_kind="binary (ORC)",
                well_pump_kW={y: main["loop_pump"][y]["kW"] for y in main["loop_pump"]})


def t09_chevron_acl():
    """Doublet + horizontals advanced closed loop (undisclosed): assumed 3.0 km verticals, 2 x 3 km laterals, 40 C/km (135 C), 35 kg/s."""
    L1, L2, n = 3000.0, 3000.0, 2
    ch = uloop("ACL doublet", L1, L2, 0.216, 35.0, 40.0, 2.0e6, [dict(k=2.6, rho=2600.0, cp=900.0)], 15.0,
               dict(gradient_C_per_km=40.0), 30.0, n_lat=n, roughness=1e-4, dz=50.0, lat_far_field="adiabatic", lat_r_far=60.0)
    res_p = ch["producer"][1]
    s = ts_summary(res_p)
    po = {y: plant_options(35.0, s[y]["h"], s[y]["p_MPa"] * 1e6, T_return_C=40.0) for y in s}
    ser = ts_series(res_p, 35.0, 40.0, "binary"); ser["T_toe_C"] = fn_series(ch["T_toe"], 30 * YEAR)["T_C"]
    lp = uloop_pump(s, 35.0, 40.0)                                    # closed loop: what it loses between the wellheads, plus the plant
    return dict(id=9, name="Doublet with horizontal closed-loop connections (Chevron-type, parameters assumed)", series=ser,
                profile=last_profile(res_p), profile_lat=last_profile(ch["lateral"][1]), profile_inj=last_profile(ch["injector"][1]),
                assumptions=dict(L1_m=L1, L2_m=L2, n_laterals=n, gradient_C_per_km=40, T_rock_C=135, mdot=35, T_in_C=40, lateral_spacing_m=120),
                T_out=s, plant=po, flags=flags(ch["injector"][1], ch["lateral"][1], res_p), plant_kind="binary (ORC)",
                p_wh_MPa={y: s[y]["p_MPa"] for y in s}, loop_pump=lp, well_pump_kW={y: lp[y]["kW"] for y in lp})


def t10_egs_doublet():
    """Classic EGS doublet (Soultz / UDDGP type): 5.0 km, 200 C granite, 20 stimulated fractures 50 m apart, 600 m between wells,
    30 kg/s, 60 C injection, fracture impedance 0.3 MPa/(L/s); pumping from the closed-loop pressure balance (loop_pump)."""
    depth, m, T_rock = 5000.0, 30.0, 200.0
    rock = [dict(k=3.0, rho=2650.0, cp=850.0)]
    geo = dict(gradient_C_per_km=(T_rock - 15.0) / 5.0)
    Q = m / 900.0
    cw = cbar_w(60.0, T_rock, 49.0e6)      # mean specific heat of the water between injection and rock temperature
    fac = lambda T_bot: ParallelFractures(Q=Q, N=20, L_flow=600.0, B=400.0, S=50.0, k_r=3.0, rhoc_r=2650 * 850, T_rock_C=T_rock,
                                           T_inj_C=lambda t: T_bot(t) - 273.15, rhoc_w=cw * 900.0)
    layers = [dict(od_in=7.0, material="steel"), dict(od_in=8.5, material="cement")]
    ch = doublet("EGS doublet", fac, depth, 6.094 * IN, m, 60.0, rock, 15.0, geo, 30.0, layers=layers, p_inj=1.0e6, p_wh=2.0e6,
                 p_res=49.0e6, PI=20e-6)   # 2 MPa wellhead keeps the stream liquid; pumping from the closed-loop balance below
    res_p = ch["producer"][1]
    s = ts_summary(res_p)
    imp = 0.3e6 / 1e-3        # Pa per m3/s
    lp = loop_pump(ch, m, imp, 60.0, list(s))
    po = {y: plant_options(m, s[y]["h"], s[y]["p_MPa"] * 1e6, T_return_C=60.0) for y in s}
    ser = ts_series(res_p, m, 60.0, "binary"); ser["T_res_C"] = fn_series(ch["T_res"], 30 * YEAR)["T_C"]
    return dict(id=10, name="Doublet EGS, hot dry rock (Soultz / UDDGP type)", series=ser, profile=last_profile(res_p), profile_inj=last_profile(ch["injector"][1]),
                assumptions=dict(depth_m=depth, T_rock_C=T_rock, N_fractures=20, L_flow_m=600, B_m=400, S_m=50, mdot=m, T_inj_C=60, p_wh_MPa=2,
                                 impedance_MPa_per_L_s=0.3, V_km3=1000 * 600 * 400 / 1e9, t_s_yr=1000 * 600 * 400 * 2650 * 850 / (m * cw) / YEAR, cbar_w=cw,
                                 loop=lp, friction_MPa=dict(injector=pipe_friction_Pa(m, 6.094 * IN, depth, 80.0, 25e6) / 1e6,
                                                            producer=pipe_friction_Pa(m, 6.094 * IN, depth, 190.0, 25e6) / 1e6)),
                T_out=s, T_res={y: ch["T_res"](y * YEAR) - 273.15 for y in s}, plant=po, well_pump_kW={y: lp[y]["kW"] for y in s},
                flags=flags(ch["injector"][1], res_p), plant_kind="binary (ORC)")


def t11_hybrid_eden():
    """Hybrid single-well EGS (Eden Project type): 4.8 km coaxial well, 35 C/km (183 C at TD), injection down the annulus into
    an upper fracture zone, return up 3-1/2 x 4-1/2 VIT from a lower zone; fracture path 150 m through 3 fractures; 10 kg/s."""
    depth = 4800.0
    casing = [dict(top_m=0, bottom_m=3000, id_in=8.535, layers=[dict(od_in=9.625, material="steel"), dict(od_in=12.25, material="cement")]),
              dict(top_m=3000, bottom_m=depth, id_in=6.094, layers=[dict(od_in=7.0, material="steel"), dict(od_in=8.5, material="cement")])]
    tubing = [dict(top_m=0, bottom_m=depth, id_in=2.992, layers=[dict(od_in=3.5, k=45), dict(od_in=4.0, k=0.03), dict(od_in=4.5, k=45)])]
    m = 10.0
    T_rock = 15.0 + 35.0 * 4.8
    pf = ParallelFractures(Q=m / 950.0, N=3, L_flow=150.0, B=150.0, S=30.0, k_r=3.2, rhoc_r=2650 * 850, T_rock_C=T_rock, T_inj_C=60.0,
                           rhoc_w=cbar_w(60.0, T_rock, 47.0e6) * 950.0)
    d = coaxial_well("Eden single-well EGS", depth, casing, tubing, m, 40.0, 1.0e6, [dict(k=3.2, rho=2650.0, cp=850.0)], 12.0,
                     dict(gradient_C_per_km=35.0), 30.0, bottom_element=pf)
    sim, res = run_case(d)
    s = ts_summary(res)
    po = {y: plant_options(m, s[y]["h"], s[y]["p_MPa"] * 1e6, T_return_C=40.0) for y in s}
    Qres = Series(res.ts["t"], res.ts["Q_res"])
    ser = ts_series(res, m, 40.0, "binary"); ser["Q_res_MW"] = [float(Qres(tt * YEAR)) / 1e6 for tt in ser["t_yr"]]
    lp = sealed_loop_pump(res, m, 40.0, list(s))                      # circulation through the coaxial well, between its two wellheads
    return dict(id=11, name="Hybrid single-well EGS with surface heat exchanger (Eden Project type)", series=ser, profile=last_profile(res),
                assumptions=dict(depth_m=depth, gradient_C_per_km=35, T_TD_C=T_rock, mdot=m, T_in_C=40, fractures="3 x 150 m path, 30 m apart"),
                T_out=s, Q_reservoir_MW={y: Qres(y * YEAR) / 1e6 for y in s}, plant=po, flags=flags(res), plant_kind="direct use; binary marginal",
                loop_pump=lp, well_pump_kW={y: lp[y]["kW"] for y in lp},
                loop_pump_note="wellbore only: the pressure drop of the fracture path at the toe is not specified in this example and not included")


def t12_sage_heatroot():
    """Single-well cyclic fracture loop (Sage HeatRoot type): 3.5 km, 150 C rock, one fracture 300 m x 100 m x 8 mm; daily cycle:
    inject 200 m3 at 40 C over 4 h, soak 12 h, produce over 8 h; the 3.5 km well adds its own losses (production leg modelled)."""
    cf = CyclicFracture(L_f=300.0, B=100.0, w=0.008, k_r=3.0, rhoc_r=2.4e6, T_rock_C=150.0, T_inj_C=40.0, n_x=60, n_z=60, z_max=25.0)
    cycles = cf.simulate(V_cycle=200.0, t_inj=4 * 3600, t_soak=12 * 3600, t_prod=8 * 3600, n_cycles=365)
    T_cycle = np.array([c["T_mean_C"] for c in cycles]); E = np.array([c["E_J"] for c in cycles])
    # production leg: 3.5 km, 7 in casing, mean rate during production 200 m3 / 8 h = 6.9 kg/s, feed at the cycle-mean temperature
    m_prod = 200.0 * 950.0 / (8 * 3600)
    Tin = Series(np.arange(len(cycles)) * 86400.0, T_cycle + 273.15)
    d = single_well("HeatRoot production leg", "production", 3500.0, 6.094 * IN, m_prod, [dict(k=3.0, rho=2600.0, cp=920.0)], 15.0,
                    dict(gradient_C_per_km=(150.0 - 15.0) / 3.5), 1.0, layers=[dict(od_in=7.0, material="steel"), dict(od_in=8.5, material="cement")],
                    inflow=dict(T_K=Tin), p_wh=0.5e6, output_years=(1,), dt_max_days=5.0)
    sim, res = run_case(d)
    ts = res.ts
    Twh = Series(ts["t"], ts["T_out"])
    days = (30, 90, 180, 365)
    s = {dd: dict(T_C=Twh(dd * 86400) - 273.15, T_fracture_C=float(T_cycle[dd - 1]), T_first_C=cycles[dd - 1]["T_first_C"], T_last_C=cycles[dd - 1]["T_last_C"],
                  E_MWh_per_cycle=float(E[dd - 1] / 3.6e9), P_mean_prod_MW=cycles[dd - 1]["P_mean_W"] / 1e6,
                  P_daily_avg_MWth=float(E[dd - 1] / 86400 / 1e6)) for dd in days}
    p = 0.5e6
    po = {dd: plant_options(m_prod / 3.0, _hT(s[dd]["T_C"], p), p, T_return_C=40.0) for dd in days}   # daily-average flow = 8/24 of m_prod
    dsel = np.unique(np.concatenate([np.arange(1, 31), np.arange(35, 366, 5), [365]]))
    series = dict(day=dsel.tolist(), T_fracture_C=[float(T_cycle[d - 1]) for d in dsel], T_wellhead_C=[float(Twh(d * 86400) - 273.15) for d in dsel],
                  T_first_C=[float(cycles[d - 1]["T_first_C"]) for d in dsel], T_last_C=[float(cycles[d - 1]["T_last_C"]) for d in dsel],
                  E_MWh=[float(E[d - 1] / 3.6e9) for d in dsel], P_daily_avg_MWth=[float(E[d - 1] / 86400 / 1e6) for d in dsel],
                  MWe=[cheap_power(m_prod / 3.0, _hT(float(Twh(d * 86400) - 273.15), p), p, "binary") for d in dsel])
    return dict(id=12, name="Single-well cyclic fracture loop (Sage HeatRoot type)", series=series, profile=last_profile(res),
                assumptions=dict(depth_m=3500, T_rock_C=150, fracture="300 m x 100 m x 8 mm", V_cycle_m3=200, t_inj_h=4, t_soak_h=12, t_prod_h=8,
                                 T_inj_C=40, m_prod_kg_s=m_prod),
                T_out=s, plant=po, flags=flags(res), plant_kind="binary on daily-average stream (batch)")


def t13_sage_heatcycle():
    """Multi-well cyclic EGS + thermal storage (Sage HeatCycle/Battery+ type): 20 HeatRoot wells phased (continuous output),
    plus storage: charge the fracture with 150 C water (surplus power through electric heaters), recover after 24 h."""
    base = t12_sage_heatroot()
    s = base["T_out"]
    cont = {dd: dict(T_C=s[dd]["T_C"], P_field_MWth=20 * s[dd]["P_daily_avg_MWth"], MWe_field=20 * base["plant"][dd]["binary_MWe"]) for dd in s}
    # storage round trip: fracture already cycled for a year (150 C rock), charge with 150 C water, soak 24 h, produce
    cf = CyclicFracture(L_f=300.0, B=100.0, w=0.008, k_r=3.0, rhoc_r=2.4e6, T_rock_C=150.0, T_inj_C=40.0, n_x=60, n_z=60, z_max=25.0)
    cf.simulate(V_cycle=200.0, t_inj=4 * 3600, t_soak=12 * 3600, t_prod=8 * 3600, n_cycles=60)
    rock_state = (cf.Tf.copy(), cf.Tr.copy(), cf.t)
    # charge: 200 m3 at 150 C over 4 h, soak 24 h, produce 8 h
    cyc = cf.simulate.__func__  # noqa (keep reference)
    cf.Tf, cf.Tr, cf.t = rock_state
    n_cells = int(round(200.0 / cf.cell_vol)); dt = cf.cell_vol / (n_cells * cf.cell_vol / (4 * 3600))
    for _ in range(n_cells):
        cf._conduct(dt); cf._shift_in(150.0 + 273.15)
    for _ in range(6):
        cf._conduct(24 * 3600 / 6)
    dt = cf.cell_vol / (n_cells * cf.cell_vol / (8 * 3600))
    Tp = []
    for _ in range(n_cells):
        cf._conduct(dt); Tp.append(cf._shift_out())
    Tp = np.array(Tp) - 273.15
    Tp_base = np.array([])
    E_in = 4.18e6 * 200.0 * (150.0 - 40.0)
    E_out = 4.18e6 * 200.0 * np.mean(Tp - 40.0)
    # what the same fracture would have given without charging (the 60-cycle baseline)
    base_c = CyclicFracture(L_f=300.0, B=100.0, w=0.008, k_r=3.0, rhoc_r=2.4e6, T_rock_C=150.0, T_inj_C=40.0, n_x=60, n_z=60, z_max=25.0)
    bc = base_c.simulate(V_cycle=200.0, t_inj=4 * 3600, t_soak=12 * 3600, t_prod=8 * 3600, n_cycles=61)
    E_base = bc[-1]["E_J"]
    p = 0.5e6
    eta_e = plant.orc_practical(1.0, _hT(float(np.mean(Tp)), p), p, T_AMB)["W_net"] / (4.18e3 * (float(np.mean(Tp)) - 40.0))
    ser = dict(base["series"]); ser["field_MWth"] = [20 * v for v in ser["P_daily_avg_MWth"]]; ser["field_MWe"] = [20 * v for v in ser["MWe"]]
    ser["storage_discharge_T_C"] = Tp.tolist(); ser["storage_discharge_frac"] = (np.arange(1, len(Tp) + 1) / len(Tp)).tolist()
    return dict(id=13, name="Multi-well cyclic EGS with subsurface thermal storage (Sage HeatCycle / Battery+ type)", series=ser,
                assumptions=dict(n_wells=20, per_well=base["assumptions"], storage_charge_T_C=150, charge_m3=200, soak_h=24),
                continuous=cont, storage=dict(T_prod_mean_C=float(np.mean(Tp)), T_prod_first_C=float(Tp[0]), T_prod_last_C=float(Tp[-1]),
                                               E_in_MWh=E_in / 3.6e9, E_out_MWh=E_out / 3.6e9, E_baseline_MWh=E_base / 3.6e9,
                                               thermal_round_trip=(E_out - E_base) / E_in, electric_round_trip=(E_out - E_base) / E_in * eta_e,
                                               orc_eta_at_T_prod=eta_e),
                T_out={dd: dict(T_C=s[dd]["T_C"]) for dd in s}, plant={dd: dict(binary_MWe=cont[dd]["MWe_field"]) for dd in s},
                flags=base["flags"], plant_kind="binary; storage round trip")


def t14_fervo():
    """Fervo Project Red type: horizontal doublet at 2.3 km TVD in 190 C rock, at the 2023 test point of 63 L/s with 60 C injection.
    Stimulated volume as published by Fervo in April 2026: 975 m of lateral x 110 m between the wells x 143 m of fracture height
    (0.0153 km3); 88 flowing fractures (the middle of the 75-100 estimated; 16 stages x 6 clusters were pumped on the injector);
    7 in casing; the producer is not pumped (as at Project Red): the injection pumps supply the loop; fracture impedance 0.15 MPa per L/s (assumed)."""
    tvd, lat, m = 2300.0, 975.0, 63.0 * 0.92
    N, sp, B, T_rock = 88, 110.0, 143.0, 190.0
    rock = [dict(k=2.8, rho=2650.0, cp=850.0)]
    geo = dict(gradient_C_per_km=(T_rock - 15.0) / 2.3)
    Q = m / 920.0
    cw = cbar_w(60.0, T_rock, 22.0e6)
    fac = lambda T_bot: ParallelFractures(Q=Q, N=N, L_flow=sp, B=B, S=lat / N, k_r=2.8, rhoc_r=2650 * 850, T_rock_C=T_rock,
                                           T_inj_C=lambda t: T_bot(t) - 273.15, rhoc_w=cw * 920.0)
    survey = dict(md_m=[0.0, tvd, tvd + lat], tvd_m=[0.0, tvd, tvd])
    layers = [dict(od_in=7.0, material="steel"), dict(od_in=8.5, material="cement")]
    ch = doublet("Fervo", fac, tvd, 6.094 * IN, m, 60.0, rock, 15.0, geo, 30.0, layers=layers, p_inj=3.0e6, p_wh=2.0e6,
                 p_res=22.0e6, PI=15e-6, prod_survey=survey, inj_survey=survey)   # 2 MPa: liquid at the wellhead; the producer is not pumped
    res_p = ch["producer"][1]
    s = ts_summary(res_p)
    imp = 0.15e6 / 1e-3
    lp = loop_pump(ch, m, imp, 60.0, list(s))
    po = {y: plant_options(m, s[y]["h"], s[y]["p_MPa"] * 1e6, T_return_C=60.0) for y in s}
    ser = ts_series(res_p, m, 60.0, "binary"); ser["T_res_C"] = fn_series(ch["T_res"], 30 * YEAR)["T_C"]
    V = lat * sp * B
    return dict(id=14, name="Horizontal doublet EGS (Fervo Project Red type)", series=ser, profile=last_profile(res_p), profile_inj=last_profile(ch["injector"][1]),
                assumptions=dict(tvd_m=tvd, lateral_m=lat, N_fractures="88 (75-100 estimated)", well_spacing_m=sp, B_m=B, S_m=lat / N, T_rock_C=T_rock,
                                 V_km3=V / 1e9, t_s_yr=V * 2650 * 850 / (m * cw) / YEAR, cbar_w=cw,
                                 mdot=m, T_inj_C=60, p_wh_MPa=2, impedance_MPa_per_L_s=0.15, loop=lp,
                                 friction_MPa=dict(injector=pipe_friction_Pa(m, 6.094 * IN, tvd + lat, 80.0, 12e6) / 1e6,
                                                   producer=pipe_friction_Pa(m, 6.094 * IN, tvd + lat, 185.0, 12e6) / 1e6)),
                T_out=s, T_res={y: ch["T_res"](y * YEAR) - 273.15 for y in s}, plant=po, well_pump_kW={y: lp[y]["kW"] for y in s},
                flags=flags(ch["injector"][1], res_p), plant_kind="binary (ORC)")


def t14b_cape():
    """Fervo Cape Station type (the Phase I pattern flow-tested in 2024): one producer between two injectors, 1,524 m (5,000 ft)
    laterals 137 m (450 ft) apart at 2.59 km TVD in 200 C rock, 100 kg/s to the producer, 60 C injection, 7 in casing, no downhole
    pump: the injection pumps drive the circulation at the reported 2,000-2,300 psi (15 MPa) wellhead pressure against a 2.07 MPa
    (300 psi) production wellhead.  Assumed: 150 flowing fractures (30 ft clusters, about nine in ten flowing), each with one wing
    towards each injector, 143 m high (the height published for Project Red)."""
    tvd, lat, m, T_rock = 2590.0, 1524.0, 100.0, 200.0
    N, sp, B = 150, 137.0, 143.0
    rock = [dict(k=2.8, rho=2650.0, cp=850.0)]
    geo = dict(gradient_C_per_km=(T_rock - 15.0) / (tvd / 1000.0))
    Q = m / 900.0
    cw = cbar_w(60.0, T_rock, 28.0e6)
    # 2N flow paths of length sp (one wing per injector), each draining a slab lat / N thick: rock volume 2 x lat x sp x B
    fac = lambda T_bot: ParallelFractures(Q=Q, N=2 * N, L_flow=sp, B=B, S=lat / N, k_r=2.8, rhoc_r=2650 * 850, T_rock_C=T_rock,
                                           T_inj_C=lambda t: T_bot(t) - 273.15, rhoc_w=cw * 900.0)
    survey = dict(md_m=[0.0, tvd, tvd + lat / 3.0], tvd_m=[0.0, tvd, tvd])    # inflow spread along the lateral: friction of the full flow over a third of it
    layers = [dict(od_in=7.0, material="steel"), dict(od_in=8.5, material="cement")]
    p_inj, p_wh = 15.0e6, 2.07e6
    ch = doublet("Cape", fac, tvd, 6.094 * IN, m, 60.0, rock, 15.0, geo, 30.0, layers=layers, p_inj=p_inj, p_wh=p_wh,
                 p_res=30.0e6, PI=50e-6, prod_survey=survey, inj_survey=survey, mdot_inj=m / 2.0)
    res_p = ch["producer"][1]
    s = ts_summary(res_p)
    ti, tp = ch["injector"][1].ts, res_p.ts
    col = lambda ts_, a, b: float(np.interp(YEAR, ts_["t"], ts_[a] - ts_[b]))
    imp = ((p_inj - p_wh) - (col(tp, "p_in", "p_out") - col(ti, "p_out", "p_in"))) / Q     # fracture impedance that reproduces the reported wellhead pressures
    lp = loop_pump(ch, m, imp, 60.0, list(s))
    po = {y: plant_options(m, s[y]["h"], s[y]["p_MPa"] * 1e6, T_return_C=60.0) for y in s}
    ser = ts_series(res_p, m, 60.0, "binary"); ser["T_res_C"] = fn_series(ch["T_res"], 30 * YEAR)["T_C"]
    V = 2 * lat * sp * B
    return dict(id=14.5, name="Horizontal EGS, one producer between two injectors (Fervo Cape Station type)", series=ser, profile=last_profile(res_p),
                profile_inj=last_profile(ch["injector"][1]),
                assumptions=dict(tvd_m=tvd, lateral_m=lat, N_fractures="150 (assumed)", well_spacing_m=sp, B_m=B, S_m=lat / N, T_rock_C=T_rock,
                                 V_km3=V / 1e9, t_s_yr=V * 2650 * 850 / (m * cw) / YEAR, cbar_w=cw, mdot=m, T_inj_C=60, p_wh_MPa=p_wh / 1e6,
                                 p_inj_MPa=p_inj / 1e6, impedance_MPa_per_L_s=imp / 1e9, loop=lp,
                                 friction_MPa=dict(injector=pipe_friction_Pa(m / 2.0, 6.094 * IN, tvd + lat / 3.0, 80.0, 25e6) / 1e6,
                                                   producer=pipe_friction_Pa(m, 6.094 * IN, tvd + lat / 3.0, 195.0, 15e6) / 1e6)),
                T_out=s, T_res={y: ch["T_res"](y * YEAR) - 273.15 for y in s}, plant=po, well_pump_kW={y: lp[y]["kW"] for y in s},
                flags=flags(ch["injector"][1], res_p), plant_kind="binary (ORC)")


def t15_coproduction():
    """Co-produced water from oil wells: 2.0 km, 95 C reservoir, 15 kg/s liquid per well up 2-7/8 in tubing (static annulus), 8 wells;
    lifting power is part of oil production.  Water at surface to an ORC."""
    depth, m = 2000.0, 15.0
    layers = [dict(od_in=2.875, material="steel"), dict(od_in=6.094, k=1.5, rho=1000.0, cp=4180.0, name="static annulus water"),
              dict(od_in=7.0, material="steel"), dict(od_in=8.5, material="cement")]
    resv = ConstantReservoir(95.0)
    d = single_well("oil well co-production", "production", depth, 2.441 * IN, m, [dict(k=2.3, rho=2500.0, cp=900.0)], 15.0,
                    dict(gradient_C_per_km=40.0), 30.0, layers=layers, inflow=dict(T_K=resv), p_wh=0.5e6)
    sim, res = run_case(d)
    s = ts_summary(res)
    po = {y: plant_options(8 * m, s[y]["h"], s[y]["p_MPa"] * 1e6, T_return_C=60.0) for y in s}
    return dict(id=15, name="Fluid co-production from producing oil / gas wells (8 wells)", series=ts_series(res, m, 60.0, "binary", m_plant=8 * m), profile=last_profile(res),
                assumptions=dict(depth_m=depth, T_res_C=95, mdot_per_well=m, n_wells=8, tubing="2-7/8 in, static water annulus"),
                T_out=s, plant=po, flags=flags(res), plant_kind="binary (ORC) on the field's water")


def t16_abandoned_well_loop():
    """Closed loop in an abandoned 3 km gas well (CeraPhi / N. Yorkshire type): 7 in casing, 110 C at TD, 3-1/2 x 4-1/2 VIT, 5 kg/s at 30 C."""
    depth = 3000.0
    casing = [dict(top_m=0, bottom_m=depth, id_in=6.094, layers=[dict(od_in=7.0, material="steel"), dict(od_in=8.5, material="cement")])]
    tubing = [dict(top_m=0, bottom_m=depth, id_in=2.992, layers=[dict(od_in=3.5, k=45), dict(od_in=4.0, k=0.03), dict(od_in=4.5, k=45)])]
    m = 5.0
    d = coaxial_well("abandoned well DBHE", depth, casing, tubing, m, 30.0, 1.0e6, [dict(k=2.5, rho=2500.0, cp=900.0)], 10.0,
                     dict(gradient_C_per_km=(110.0 - 10.0) / 3.0), 30.0)
    sim, res = run_case(d)
    s = ts_summary(res)
    po = {y: plant_options(m, s[y]["h"], s[y]["p_MPa"] * 1e6, T_return_C=30.0) for y in s}
    lp = sealed_loop_pump(res, m, 30.0, list(s))
    return dict(id=16, name="Closed loop in a single abandoned hydrocarbon well (CeraPhi / wells2watts type)", series=ts_series(res, m, 30.0, "binary"), profile=last_profile(res),
                assumptions=dict(depth_m=depth, T_TD_C=110, mdot=m, T_in_C=30, VIT="3-1/2 x 4-1/2, k 0.03"),
                T_out=s, plant=po, flags=flags(res), plant_kind="direct use (heat network)", loop_pump=lp, well_pump_kW={y: lp[y]["kW"] for y in lp})


def t17_thermal_eor():
    """Recovery of injected steam heat after SAGD: chamber 1000 x 100 x 30 m at 300 m depth, 200 C, swept by a cold-water doublet at
    20 kg/s (HT ESP), 15 years; lumped store with conductive losses; the 300 m wells lose little."""
    hs = HeatStore(V=1000 * 100 * 30.0, A_b=2 * (1000 * 100 + 1000 * 30 + 100 * 30), rhoc_store=2.3e6, k_r=1.7, rhoc_r=2.3e6, T_store_C=200.0, T0_C=12.0)
    r = hs.run(15.0, 20.0, 50.0)
    yrs = (1, 2, 5, 10, 15)
    s = {y: dict(T_C=float(np.interp(y * YEAR, r["t"], r["T_C"])), Q_MWth=float(np.interp(y * YEAR, r["t"], r["Q_MW"])),
                 loss_MW=float(np.interp(y * YEAR, r["t"], r["loss_MW"]))) for y in yrs}
    E0 = hs.V * hs.rhoc * (200.0 - 12.0)
    E_prod = np.sum(r["Q_MW"] * 1e6 * 10 * 86400)
    p = 1.0e6
    po = {y: plant_options(20.0, _hT(max(s[y]["T_C"], 51.0), p), p, T_return_C=50.0) for y in yrs}
    tg = np.unique(np.concatenate([np.geomspace(10 * 86400.0, r["t"][-1], 70), [r["t"][-1]]]))
    series = dict(t_yr=(tg / YEAR).tolist(), T_C=np.interp(tg, r["t"], r["T_C"]).tolist(), heat_MWth=np.interp(tg, r["t"], r["Q_MW"]).tolist(),
                  loss_MW=np.interp(tg, r["t"], r["loss_MW"]).tolist(),
                  MWe=[cheap_power(20.0, _hT(max(float(np.interp(tt, r["t"], r["T_C"])), 51.0), p), p, "binary") if np.interp(tt, r["t"], r["Q_MW"]) > 0 else 0.0 for tt in tg])
    return dict(id=17, name="Recovery of man-made heat after thermal EOR (SAGD chamber, HT ESP circulation)", series=series,
                assumptions=dict(chamber_m="1000 x 100 x 30", depth_m=300, T_store_C=200, mdot=20, T_inj_C=50, years=15),
                T_out=s, plant=po, recovered_fraction_15yr=float(E_prod / E0), flags=dict(unconverged=0, choke=0, max_residual_kW=0.0),
                plant_kind="binary (ORC), declining")


def t18_supercritical_field():
    """Supercritical under an existing field (IDDP-2 type): 4.5 km, 450 C / 30 MPa feed, 20 kg/s, 100 C/km, wellhead 10 MPa."""
    depth, m = 4500.0, 20.0
    p_res = 30.0e6
    h_res = CP.PropsSI("H", "T", 450 + 273.15, "P", p_res, "Water")
    casing = [dict(top_m=0, bottom_m=3000, id_in=8.535, layers=[dict(od_in=9.625, material="steel"), dict(od_in=12.25, material="cement")]),
              dict(top_m=3000, bottom_m=depth, id_in=6.094, layers=[dict(od_in=7.0, material="steel"), dict(od_in=8.5, material="cement")])]
    d = dict(name="supercritical well", fluid="iapws",
             well=dict(mode="production", depth_m=depth, casing=casing, inflow=dict(h_J_kg=h_res, p_res_MPa=30.0, PI_kg_s_per_MPa=10.0)),
             materials=dict(steel=STEEL, cement=CEMENT),
             rock=dict(layers=[dict(k=2.2, rho=2700.0, cp=950.0)], surface_T_C=15.0, geotherm=dict(gradient_C_per_km=100.0)),
             operation=dict(pressure_bc=dict(outlet_MPa=10.0), schedule=[dict(duration_days=30 * 365.25, mdot_kg_s=m, T_in_C=40.0)]),
             numerics=dict(dz_m=50.0, dt_initial_s=300.0, dt_max_days=60.0, output_times_days=[365.25 * y for y in YRS]))
    sim, res = run_case(d)
    s = ts_summary(res)
    po = {y: plant_options(m, s[y]["h"], s[y]["p_MPa"] * 1e6) for y in s}
    return dict(id=18, name="Supercritical fluid produced from below a volcanic field (IDDP-2 type)", series=ts_series(res, m, 70.0, "steam"), profile=last_profile(res),
                assumptions=dict(depth_m=depth, T_res_C=450, p_res_MPa=30, mdot=m, p_wh_MPa=10, gradient_C_per_km=100),
                T_out=s, plant=po, flags=flags(res), plant_kind="superheated-steam condensing turbine")


def t19_superhot_egs():
    """Superhot EGS (Newberry / AltaRock type): 3.5 km, 375 C rock, doublet through 20 fractures 40 m apart, 40 kg/s, 60 C injection,
    production up 7 in casing to a 5 MPa wellhead (flashes to superheated steam)."""
    depth, m, T_rock = 3500.0, 40.0, 375.0
    rock = [dict(k=2.5, rho=2700.0, cp=950.0)]
    geo = dict(gradient_C_per_km=(T_rock - 15.0) / 3.5)
    Q = m / 900.0
    cw = cbar_w(60.0, T_rock, 34.0e6)
    fac = lambda T_bot: ParallelFractures(Q=Q, N=20, L_flow=300.0, B=400.0, S=40.0, k_r=2.5, rhoc_r=2700 * 950, T_rock_C=T_rock,
                                           T_inj_C=lambda t: T_bot(t) - 273.15, rhoc_w=cw * 900.0)
    layers = [dict(od_in=7.0, material="steel"), dict(od_in=8.5, material="cement")]
    ch = doublet("superhot EGS", fac, depth, 6.094 * IN, m, 60.0, rock, 15.0, geo, 30.0, layers=layers, p_inj=5.0e6, p_wh=5.0e6,
                 p_res=34.0e6, PI=10e-6)
    res_p = ch["producer"][1]
    s = ts_summary(res_p)
    imp = 0.3e6 / 1e-3
    lp = loop_pump(ch, m, imp, 60.0, list(s))
    po = {y: plant_options(m, s[y]["h"], s[y]["p_MPa"] * 1e6) for y in s}
    ser = ts_series(res_p, m, 70.0, "steam"); ser["T_res_C"] = fn_series(ch["T_res"], 30 * YEAR)["T_C"]
    return dict(id=19, name="Superhot-rock EGS doublet (Newberry / AltaRock type)", series=ser, profile=last_profile(res_p), profile_inj=last_profile(ch["injector"][1]),
                assumptions=dict(depth_m=depth, T_rock_C=T_rock, N_fractures=20, L_flow_m=300, B_m=400, S_m=40, mdot=m, T_inj_C=60, p_wh_MPa=5,
                                 impedance_MPa_per_L_s=0.3, V_km3=800 * 300 * 400 / 1e9, t_s_yr=800 * 300 * 400 * 2700 * 950 / (m * cw) / YEAR, cbar_w=cw, loop=lp,
                                 friction_MPa=dict(injector=pipe_friction_Pa(m, 6.094 * IN, depth, 80.0, 20e6) / 1e6)),
                T_out=s, T_res={y: ch["T_res"](y * YEAR) - 273.15 for y in s}, plant=po, well_pump_kW={y: lp[y]["kW"] for y in s},
                flags=flags(ch["injector"][1], res_p), plant_kind="steam turbine on the flashed / superheated stream")


def t20_quaise():
    """Millimetre-wave drilled closed loop (Quaise type): 20 km coaxial, 500 C at TD (24.25 C/km), 0.30 m vitrified bore
    (2 cm glass, k 1.4), 5-1/2 x 7 in VIT, 20 kg/s at 40 C, 10 MPa delivery."""
    depth = 20000.0
    casing = [dict(top_m=0, bottom_m=depth, id_m=0.30, layers=[dict(od_m=0.34, k=1.4, rho=2400.0, cp=900.0, name="glass")])]
    tubing = [dict(top_m=0, bottom_m=depth, id_in=4.892, layers=[dict(od_in=5.5, k=45), dict(od_in=6.0, k=0.03), dict(od_in=6.5, k=45)])]
    out = {}
    for m in (20.0, 40.0):
        d = coaxial_well(f"Quaise {m} kg/s", depth, casing, tubing, m, 40.0, 10.0e6, [dict(k=2.8, rho=2700.0, cp=950.0)], 15.0,
                         dict(gradient_C_per_km=(500.0 - 15.0) / 20.0), 30.0, dz=100.0)
        sim, res = run_case(d)
        s = ts_summary(res)
        po = {y: plant_options(m, s[y]["h"], s[y]["p_MPa"] * 1e6, T_return_C=40.0) for y in s}
        out[m] = dict(T_out=s, plant=po, flags=flags(res), series=ts_series(res, m, 40.0, "best"), profile=last_profile(res),
                      loop_pump=sealed_loop_pump(res, m, 40.0, list(s)))
    m = 20.0
    return dict(id=20, name="Ultradeep millimetre-wave drilled closed loop (Quaise type)", series=out[m]["series"], profile=out[m]["profile"],
                assumptions=dict(depth_m=depth, T_TD_C=500, gradient_C_per_km=24.25, bore_m=0.30, mdot=[20, 40], T_in_C=40, p_out_MPa=10,
                                 VIT="5-1/2 x 7 in, k 0.03"),
                T_out=out[m]["T_out"], plant=out[m]["plant"], flags=out[m]["flags"], variants={str(k): v for k, v in out.items()},
                plant_kind="flash / supercritical ORC", well_pump_kW={y: out[m]["loop_pump"][y]["kW"] for y in out[m]["loop_pump"]})


UNAIDED = (3, 3.5, 18)        # single production wells that flow without a pump: checked for deliverability (inflow_short_MPa)
TYPES = {1: t01_lt_hsa, 2: t02_mt_binary, 3: t03_ht_flash, 3.5: t03b_dry_steam, 4: t04_waste_heat, 5: t05_gshp, 6: t06_shallow_coaxial,
         7: t07_greenloop, 8: t08_eavor, 9: t09_chevron_acl, 10: t10_egs_doublet, 11: t11_hybrid_eden, 12: t12_sage_heatroot,
         13: t13_sage_heatcycle, 14: t14_fervo, 14.5: t14b_cape, 15: t15_coproduction, 16: t16_abandoned_well_loop, 17: t17_thermal_eor,
         18: t18_supercritical_field, 19: t19_superhot_egs, 20: t20_quaise}


def _run(tid):
    t0 = time.time()
    try:
        out = TYPES[tid]()
        out["wall_s"] = time.time() - t0
        if tid in UNAIDED and isinstance(out.get("flags"), dict):
            mg = [v.get("margin_MPa") for v in out.get("T_out", {}).values() if isinstance(v, dict) and v.get("margin_MPa") is not None]
            out["flags"]["inflow_short_MPa"] = max(0.0, -min(mg)) if mg else None
        T = out.get("T_out", {})
        key = max(T) if T else None
        print(f"type {tid}: done in {out['wall_s']:.0f} s; T_out at {key}: {T[key]['T_C'] if key is not None else 'n/a'}", flush=True)
        return out
    except Exception as e:
        print(f"type {tid}: FAILED {type(e).__name__}: {e}", flush=True)
        traceback.print_exc()
        return dict(id=tid, failed=f"{type(e).__name__}: {e}", wall_s=time.time() - t0)


def _clean(o):
    if isinstance(o, dict):
        return {str(k): _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    if isinstance(o, np.ndarray):
        return _clean(o.tolist())
    if isinstance(o, (np.floating, float)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (str, int, bool)) or o is None:
        return o
    return str(o)


if __name__ == "__main__":
    ids = [float(x) if "." in x else int(x) for x in sys.argv[1:]] or list(TYPES)
    workers = 2
    results = {}
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(_run, tid): tid for tid in ids}
        for f in as_completed(futs):
            results[futs[f]] = _clean(f.result())
            try:
                prev = {r["id"]: r for r in load("DT_types")}
            except Exception:
                prev = {}
            prev.update({results[k]["id"]: results[k] for k in results})
            save("DT_types", [prev[k] for k in sorted(prev, key=lambda x: float(x))])
    print("done")
