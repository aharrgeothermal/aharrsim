"""Chains of well segments and development-type templates.

A development is a chain of one-way coupled components: each segment is a aharrsim
run (coaxial loop, single production or injection channel, U-tube), a reservoir
element (aharrsim.reservoir) or a surface plant (aharrsim.plant).  Segments are run in
order; the outlet temperature history of one feeds the inlet of the next through a
callable T(t), so transients propagate without iteration (the flow is one-way).
"""
from __future__ import annotations

import copy

import numpy as np

from .case import case_from_dict
from .simulator import Simulator

YEAR = 365.25 * 86400.0


class Series:
    """Piecewise-linear interpolant of a time series, clamped at the ends."""

    def __init__(self, t, y):
        self.t = np.asarray(t, float)
        self.y = np.asarray(y, float)

    def __call__(self, t):
        return float(np.interp(t, self.t, self.y))

    def shifted(self, dT):
        return Series(self.t, self.y + dT)

    def mapped(self, fn):
        return Series(self.t, np.array([fn(v) for v in self.y]))


def run_case(d: dict, progress=False):
    """Run a case dictionary; returns (sim, res)."""
    sim = Simulator(case_from_dict(copy.deepcopy(d)))
    res = sim.run(progress=progress)
    return sim, res


def outlet_series(res):
    """Outlet temperature [K] and pressure [Pa] histories of a run as Series."""
    ts = res.ts
    return Series(ts["t"], ts["T_out"]), Series(ts["t"], ts["p_out"])


def vertical_well(mode, depth, d_in, layers=None, roughness=4.6e-5, top=0.0):
    """Casing list for a vertical single-channel well of inner diameter d_in [m]."""
    return [dict(top_m=top, bottom_m=depth, id_m=d_in, layers=layers or [])]


def lateral_case(name, L_lat, tvd, d_in, mdot, T_in, p_in, rock, surface_T_C, geotherm, years, layers=None,
                 roughness=4.6e-5, dz=50.0, far_field="dirichlet", r_far=None, output_years=(1, 5, 10, 20, 30, 40)):
    """A horizontal segment at depth tvd carrying mdot from md 0 to L_lat (one-way, no gravity)."""
    d = dict(name=name, fluid="iapws",
             well=dict(mode="injection", survey=dict(md_m=[0.0, L_lat], tvd_m=[tvd, tvd]), roughness_m=roughness,
                       casing=[dict(top_m=0.0, bottom_m=L_lat, id_m=d_in, layers=layers or [])],
                       inflow=dict(T_K=T_in)),
             rock=dict(layers=rock, surface_T_C=surface_T_C, geotherm=geotherm),
             operation=dict(pressure_bc=dict(inlet_MPa=p_in / 1e6),
                            schedule=[dict(duration_days=years * 365.25, mdot_kg_s=mdot, T_in_C=20.0)]),
             numerics=dict(dz_m=dz, dt_initial_s=300.0, dt_max_days=60.0,
                           output_times_days=[365.25 * y for y in output_years if y <= years]),
             options=dict(far_field=far_field))
    if r_far is not None:
        d["numerics"]["r_far_m"] = r_far
    return d


def uloop(name, L1, L2, d_in, mdot, T_in_C, p_inj, rock, surface_T_C, geotherm, years, n_lat=1, layers=None,
          roughness=4.6e-5, dz=50.0, lat_far_field="dirichlet", lat_r_far=None, progress=False):
    """Injector (down) -> n_lat identical laterals (each mdot/n_lat) -> producer (up).  Returns dict of runs and
    the production wellhead temperature series."""
    out_years = [y for y in (1, 5, 10, 20, 30, 40) if y <= years]
    inj = dict(name=f"{name}: injector", fluid="iapws",
               well=dict(mode="injection", depth_m=L1, roughness_m=roughness,
                         casing=[dict(top_m=0.0, bottom_m=L1, id_m=d_in, layers=layers or [])]),
               rock=dict(layers=rock, surface_T_C=surface_T_C, geotherm=geotherm),
               operation=dict(pressure_bc=dict(inlet_MPa=p_inj / 1e6),
                              schedule=[dict(duration_days=years * 365.25, mdot_kg_s=mdot, T_in_C=T_in_C)]),
               numerics=dict(dz_m=dz, dt_initial_s=300.0, dt_max_days=60.0, output_times_days=[365.25 * y for y in out_years]))
    sim_i, res_i = run_case(inj, progress)
    T_bot, p_bot = outlet_series(res_i)
    lat = lateral_case(f"{name}: lateral", L2, L1, d_in, mdot / n_lat, T_bot, p_bot.y[-1], rock, surface_T_C, geotherm, years,
                       layers=layers, roughness=roughness, dz=dz, far_field=lat_far_field, r_far=lat_r_far, output_years=out_years)
    sim_l, res_l = run_case(lat, progress)
    T_toe, p_toe = outlet_series(res_l)
    prod = dict(name=f"{name}: producer", fluid="iapws",
                well=dict(mode="production", depth_m=L1, roughness_m=roughness, inflow=dict(T_K=T_toe),
                          casing=[dict(top_m=0.0, bottom_m=L1, id_m=d_in, layers=layers or [])]),
                rock=dict(layers=rock, surface_T_C=surface_T_C, geotherm=geotherm),
                operation=dict(pressure_bc=dict(bottom_MPa=p_toe.y[-1] / 1e6),
                               schedule=[dict(duration_days=years * 365.25, mdot_kg_s=mdot, T_in_C=T_in_C)]),
                numerics=dict(dz_m=dz, dt_initial_s=300.0, dt_max_days=60.0, output_times_days=[365.25 * y for y in out_years]))
    sim_p, res_p = run_case(prod, progress)
    T_wh, p_wh = outlet_series(res_p)
    return dict(injector=(sim_i, res_i), lateral=(sim_l, res_l), producer=(sim_p, res_p), T_wh=T_wh, p_wh=p_wh,
                T_bottom_injector=T_bot, T_toe=T_toe)


# ----------------------------------------------------------------------------
# Generic single-channel wells and chains used by the development-type templates
# ----------------------------------------------------------------------------
def single_well(name, mode, depth, d_in, mdot, rock, surface_T_C, geotherm, years, layers=None, T_in_C=40.0,
                inflow=None, p_wh=None, p_inj=None, p_bottom=None, pump=None, survey=None, dz=50.0, roughness=4.6e-5,
                far_field="dirichlet", r_far=None, schedule=None, output_years=(1, 5, 10, 20, 30), fluid="iapws",
                two_phase_model="drift_flux", dt_max_days=60.0):
    """Case dictionary for a production ('production') or injection ('injection') well.
    inflow: dict for production wells, e.g. {'T_K': callable or K} or {'h_J_kg': ...}, optional p_res_MPa, PI_kg_s_per_MPa.
    pump: {'md_m': depth, 'efficiency': 0.7} (production only).  schedule: list of periods overriding the single period."""
    well = dict(mode=mode, roughness_m=roughness,
                casing=[dict(top_m=0.0, bottom_m=(survey["md_m"][-1] if survey else depth), id_m=d_in, layers=layers or [])])
    if survey:
        well["survey"] = survey
    else:
        well["depth_m"] = depth
    if inflow:
        well["inflow"] = inflow
    if pump:
        well["pump"] = pump
    if mode == "production":
        pbc = dict(bottom_MPa=p_bottom / 1e6) if p_bottom else dict(outlet_MPa=(p_wh or 1.0e6) / 1e6)
    else:
        pbc = dict(inlet_MPa=(p_inj or 0.5e6) / 1e6)
    sched = schedule or [dict(duration_days=years * 365.25, mdot_kg_s=mdot, T_in_C=T_in_C)]
    d = dict(name=name, fluid=fluid, well=well,
             materials=dict(steel=dict(k=45.0, rho=7850.0, cp=490.0), cement=dict(k=1.0, rho=1900.0, cp=1000.0)),
             rock=dict(layers=rock, surface_T_C=surface_T_C, geotherm=geotherm),
             operation=dict(pressure_bc=pbc, schedule=sched),
             numerics=dict(dz_m=dz, dt_initial_s=300.0, dt_max_days=dt_max_days, two_phase_model=two_phase_model,
                           output_times_days=[365.25 * y for y in output_years if y <= years]),
             options=dict(far_field=far_field))
    if r_far is not None:
        d["numerics"]["r_far_m"] = r_far
    return d


def coaxial_well(name, depth, casing, tubing, mdot, T_in_C, p_out, rock, surface_T_C, geotherm, years,
                 bottom_element=None, far_field="dirichlet", r_far=None, dz=50.0, output_years=(1, 5, 10, 20, 30),
                 schedule=None, materials=None, flow_direction="annulus_down"):
    d = dict(name=name, fluid="iapws",
             well=dict(depth_m=depth, casing=casing, tubing=tubing, flow_direction=flow_direction),
             materials=materials or dict(steel=dict(k=45.0, rho=7850.0, cp=490.0), cement=dict(k=1.0, rho=1900.0, cp=1000.0)),
             rock=dict(layers=rock, surface_T_C=surface_T_C, geotherm=geotherm),
             operation=dict(pressure_bc=dict(outlet_MPa=p_out / 1e6),
                            schedule=schedule or [dict(duration_days=years * 365.25, mdot_kg_s=mdot, T_in_C=T_in_C)]),
             numerics=dict(dz_m=dz, dt_initial_s=300.0, dt_max_days=60.0, two_phase_model="drift_flux",
                           output_times_days=[365.25 * y for y in output_years if y <= years]),
             options=dict(far_field=far_field))
    if bottom_element is not None:
        d["well"]["bottom_element"] = bottom_element
    if r_far is not None:
        d["numerics"]["r_far_m"] = r_far
    return d


def doublet(name, reservoir_factory, depth, d_in, mdot, T_inj_C, rock, surface_T_C, geotherm, years, layers=None,
            p_inj=0.5e6, p_wh=1.0e6, pump=None, p_res=None, PI=None, dz=50.0, prod_survey=None, inj_survey=None,
            prod_depth=None, progress=False, mdot_inj=None):
    """Injector -> reservoir element -> producer.  reservoir_factory(T_inj_series) returns a callable T(t) [K].
    mdot_inj: flow in one injector when it differs from the producer's (a producer fed by two injectors).
    Returns dict(injector, producer, reservoir, T_wh Series, T_res Series)."""
    inj = single_well(f"{name}: injector", "injection", depth, d_in, mdot_inj or mdot, rock, surface_T_C, geotherm, years, layers=layers,
                      T_in_C=T_inj_C, p_inj=p_inj, dz=dz, survey=inj_survey)
    sim_i, res_i = run_case(inj, progress)
    T_bot, _ = outlet_series(res_i)
    resv = reservoir_factory(T_bot)
    inflow = dict(T_K=resv)
    if p_res is not None and PI is not None:
        inflow.update(p_res_MPa=p_res / 1e6, PI_kg_s_per_MPa=PI * 1e6)
    prod = single_well(f"{name}: producer", "production", prod_depth or depth, d_in, mdot, rock, surface_T_C, geotherm, years,
                       layers=layers, inflow=inflow, p_wh=p_wh, pump=pump, dz=dz, survey=prod_survey)
    sim_p, res_p = run_case(prod, progress)
    T_wh, p_wh_s = outlet_series(res_p)
    T_res = Series(res_p.ts["t"], [resv(t) for t in res_p.ts["t"]])
    return dict(injector=(sim_i, res_i), producer=(sim_p, res_p), reservoir=resv, T_wh=T_wh, p_wh=p_wh_s, T_res=T_res,
                T_bottom_injector=T_bot)


def at_years(series, years):
    return [series(y * YEAR) for y in years]
