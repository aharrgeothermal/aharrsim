"""Validation against measured field data from deep coaxial borehole heat exchangers.

Every input is taken from the cited (secondary) sources; anything the sources do
not give is an explicit assumption, listed in ASSUMPTIONS and varied in the
sensitivity runs.
"""
from __future__ import annotations

import copy
import sys

import numpy as np

from common import save, run, C, DAY, YEAR

STEEL = dict(k=45.0, rho=7850.0, cp=490.0)


def qingdao(k_rock=3.49, k_cem=0.73, k_pipe=0.21, rhoc_cem=(1900.0, 1000.0)):
    """Jimo (Qingdao) 2605 m DBHE, 138-day constant-inlet test.

    Source: Bu, Ran & Zhang (2019) Renewable Energy 143:1902-1909, as tabulated by
    Chen & Tomac (2023, eScholarship qt4s97c8n2) Sec. 2.6 / Table 1 and Kolo et
    al. (2024) Geothermal Energy, Table 4.
    """
    return dict(
        name="Qingdao (Jimo) 2605 m",
        well=dict(depth_m=2605.0,
                  casing=[dict(id_m=0.1778 - 2 * 0.00691,
                               layers=[dict(od_m=0.1778, material="steel"),
                                       dict(od_m=0.1778 + 2 * 0.01905, material="cement")])],
                  tubing=[dict(id_m=0.110 - 2 * 0.010, layers=[dict(od_m=0.110, k=k_pipe)])]),
        materials=dict(steel=STEEL, cement=dict(k=k_cem, rho=rhoc_cem[0], cp=rhoc_cem[1])),
        rock=dict(layers=[dict(k=k_rock, rho=2800.0, cp=920.0)], surface_T_C=15.0,
                  geotherm=dict(gradient_C_per_km=27.8)),
        operation=dict(pressure_bc=dict(outlet_MPa=0.3),
                       schedule=[dict(duration_days=138.0, mdot_kg_s=30.0 / 3.6, T_in_C=5.0)]),
        numerics=dict(dz_m=25.0, dt_initial_s=60.0, output_times_days=[1, 30, 138]),
    )


def songyuan(k_pipe=0.40, rhoc=(2500.0, 920.0), k_cem=1.0, r_pipe_od=0.110):
    """Songyuan (Jilin) 2044 m DBHE, 60-day heat-pump operation.

    Source: Huang et al. (2020) Energy 210:118643, as tabulated by Chen & Tomac
    (2023) Sec. 2.7 and Kolo et al. (2024) Table 4.  The inlet temperature fell
    from 36.7 C (day 1) to 29.1 C (day 60); the record between is graphical, so
    a linear decline is assumed.
    """
    sched = []
    for d in range(60):
        Tin = 36.7 + (29.1 - 36.7) * d / 59.0
        sched.append(dict(duration_days=1.0, mdot_kg_s=30.0 / 3.6 * 0.99, T_in_C=Tin))
    return dict(
        name="Songyuan 2044 m",
        well=dict(depth_m=2044.0,
                  casing=[dict(id_m=2 * 0.0809,
                               layers=[dict(od_m=0.1778, material="steel"),
                                       dict(od_m=0.2159, material="cement")])],
                  tubing=[dict(id_m=r_pipe_od - 2 * 0.0123, layers=[dict(od_m=r_pipe_od, k=k_pipe)])]),
        materials=dict(steel=STEEL, cement=dict(k=k_cem, rho=1900.0, cp=1000.0)),
        rock=dict(layers=[dict(k=2.5, rho=rhoc[0], cp=rhoc[1])], surface_T_C=107.3 - 50.7 * 2.044,
                  geotherm=dict(gradient_C_per_km=50.7)),
        operation=dict(pressure_bc=dict(outlet_MPa=0.3), schedule=sched),
        numerics=dict(dz_m=25.0, dt_initial_s=60.0, dt_max_days=1.0, output_times_days=[1, 60]),
    )


def xian(k_pipe=0.45, k_cem=1.0, d_hole=0.2159, cas_id=0.170, pipe=(0.110, 0.010)):
    """Xi'an 2000 m demonstration DBHE, one month continuous operation.

    Source: Wang Z. et al. (2017) ECM 153:603-615 and Cai et al. (2019) ATE
    149:975-986, as tabulated by Zhang J. et al. (2022) Processes 10:121/10:888
    and Kolo et al. (2024) Table 4.  Borehole, grout and pipe dimensions are not
    given in the accessible sources: assumed (casing ID 0.170 m from the 0.085 m
    model radius in Zhang et al.).
    """
    od, wall = pipe
    return dict(
        name="Xi'an 2000 m",
        well=dict(depth_m=2000.0,
                  casing=[dict(id_m=cas_id, layers=[dict(od_m=0.1778, material="steel"),
                                                     dict(od_m=d_hole, material="cement")])],
                  tubing=[dict(id_m=od - 2 * wall, layers=[dict(od_m=od, k=k_pipe)])]),
        materials=dict(steel=dict(k=40.0, rho=7850.0, cp=498.0), cement=dict(k=k_cem, rho=1900.0, cp=1000.0)),
        rock=dict(layers=[dict(top_m=0, bottom_m=500, k=1.59, rho=1760, cp=1433),
                          dict(top_m=500, bottom_m=720, k=1.65, rho=1860, cp=1025),
                          dict(top_m=720, bottom_m=1540, k=1.76, rho=2070, cp=878),
                          dict(top_m=1540, bottom_m=1e5, k=1.88, rho=2270, cp=848)],
                  surface_T_C=15.6, geotherm=dict(gradient_C_per_km=30.0)),
        operation=dict(pressure_bc=dict(outlet_MPa=0.3),
                       schedule=[dict(duration_days=30.0, mdot_kg_s=28.0 / 3.6 * 0.999, T_in_C=17.3)]),
        numerics=dict(dz_m=25.0, dt_initial_s=60.0, output_times_days=[26, 30]),
    )


def weggis(k_rock=2.6, rhoc=(2450.0, 950.0), bht=75.5, T_s=10.0, prior_load=0.5):
    """Weggis (Switzerland) 2302 m DBHE, continuous operation 2001-2003.

    Sources: Eugster & Fueglister, BFE project 38'330 final report (2003) Tables
    1, 5, 10 (measured energies and mean temperatures); geometry from Geothermal
    Energy (2024) s40517-024-00315-4 Table 1 (citing Eugster & Fueglister 2001 and
    Kohl et al. 2002).  Rock properties and the undisturbed temperature profile
    are not in the accessible sources: assumed Molasse values, varied below.
    Operation 1994-2001 was intermittent (about 220-245 MWh/yr); represented
    here as continuous circulation at a reduced flow (prior_load x the 2001-03
    flow) before two years at the 2001-03 mean conditions.
    """
    m = 1.58          # kg/s, from 420 MWh/yr at a mean 7.3 K rise
    return dict(
        name="Weggis 2302 m",
        well=dict(depth_m=2302.0,
                  casing=[dict(id_m=0.177, layers=[dict(od_m=0.197, material="steel"),
                                                    dict(od_m=0.279, material="cement")])],
                  tubing=[dict(top_m=0, bottom_m=1780, id_m=0.073,
                               layers=[dict(od_m=0.105, k=0.09)]),
                          dict(top_m=1780, bottom_m=2302, id_m=0.073,
                               layers=[dict(od_m=0.105, k=50.0)])]),
        materials=dict(steel=dict(k=50.0, rho=7850.0, cp=490.0), cement=dict(k=2.0, rho=1900.0, cp=1000.0)),
        rock=dict(layers=[dict(k=k_rock, rho=rhoc[0], cp=rhoc[1])], surface_T_C=T_s,
                  geotherm=dict(gradient_C_per_km=(bht - T_s) / 2.302)),
        operation=dict(pressure_bc=dict(outlet_MPa=0.3),
                       schedule=[dict(duration_days=7 * 365.25, mdot_kg_s=m * prior_load, T_in_C=33.3, label="1994-2001"),
                                 dict(duration_days=2 * 365.25, mdot_kg_s=m, T_in_C=29.9, label="2001-2003")]),
        numerics=dict(dz_m=25.0, dt_initial_s=300.0, dt_max_days=30.0, output_times_days=[]),
    )


MEASURED = {
    "qingdao": dict(mean_kW=448.49, T_out_end=17.5,
                    note="138-day mean extraction 448.49 kW (172 W/m); outlet 'about 17.5 C' at test end"),
    "songyuan": dict(day1=dict(T_in=36.7, T_out=58.2), day60=dict(T_in=29.1, T_out=43.0, q_W_m=238.0),
                     note="inlet/outlet 36.7/58.2 C day 1 and 29.1/43.0 C day 60; 238 W/m at day 60"),
    "xian": dict(T_out_last4=(25.4, 26.8), mean_kW=286.4,
                 note="mean inlet 17.3 C at 28 m3/h; outlet 25.4-26.8 C over the last 4 days; mean 286.4 kW"),
    "weggis": dict(MWh_yr=(420.5, 419.8), T_in=29.9, T_out=37.2,
                   note="net heat from the DBHE 420.5 MWh (2001/02) and 419.8 MWh (2002/03); mean inlet 29.9 C, outlet 37.2 C"),
}


def run_qingdao():
    out = {}
    variants = {"base": {}, "k_rock -10%": dict(k_rock=3.49 * 0.9), "k_rock +10%": dict(k_rock=3.49 * 1.1),
                "grout k 1.0": dict(k_cem=1.0), "pipe k 0.30": dict(k_pipe=0.30)}
    for lab, kw in variants.items():
        sim, res = run(qingdao(**kw))
        mean_kW = res.mean_over(0, 138 * DAY, "Q_wellhead") / 1e3
        mean_kW_d1 = res.mean_over(1 * DAY, 138 * DAY, "Q_wellhead") / 1e3
        Tend = float(C(res.at(138 * DAY, "T_out")))
        out[lab] = dict(mean_kW=mean_kW, mean_kW_after_day1=mean_kW_d1, T_out_end=Tend,
                        Q_end_kW=res.at(138 * DAY, "Q_wellhead") / 1e3)
        if lab == "base":
            out["series"] = dict(t_days=res.t_days, T_out=C(res.ts["T_out"]), Q_kW=res.ts["Q_wellhead"] / 1e3)
            pr = res.profiles[-1]
            out["profile"] = dict(md=pr["md"], T_a=C(pr["T_a"]), T_t=C(pr["T_t"]), T0=C(pr["T0"]))
        print(f"Qingdao {lab:12s} mean {mean_kW:7.1f} kW (after day 1 {mean_kW_d1:7.1f})  Tout(138 d) {Tend:6.2f} C")
    out["measured"] = MEASURED["qingdao"]
    save("F_qingdao", out)


def run_songyuan():
    out = {}
    variants = {"base": {}, "pipe k 0.30": dict(k_pipe=0.30), "pipe k 0.50": dict(k_pipe=0.50),
                "rho c 2.0e6": dict(rhoc=(2300.0, 870.0)), "rho c 2.6e6": dict(rhoc=(2600.0, 1000.0)),
                "pipe OD 90 mm": dict(r_pipe_od=0.090)}
    for lab, kw in variants.items():
        sim, res = run(songyuan(**kw))
        T1 = float(C(res.at(1 * DAY, "T_out")))
        T60 = float(C(res.at(60 * DAY, "T_out")))
        q60 = res.at(60 * DAY, "Q_wellhead") / 2044.0
        out[lab] = dict(T_out_day1=T1, T_out_day60=T60, q_W_m_day60=q60)
        if lab == "base":
            out["series"] = dict(t_days=res.t_days, T_out=C(res.ts["T_out"]), T_in=C(res.ts["T_in"]),
                                 Q_kW=res.ts["Q_wellhead"] / 1e3)
        print(f"Songyuan {lab:14s} Tout day1 {T1:6.2f} C  day60 {T60:6.2f} C  q60 {q60:6.1f} W/m")
    out["measured"] = MEASURED["songyuan"]
    save("F_songyuan", out)


def run_xian():
    out = {}
    variants = {"base": {}, "grout k 1.5": dict(k_cem=1.5), "pipe k 0.30": dict(k_pipe=0.30),
                "hole 0.25 m": dict(d_hole=0.25), "pipe 90x8 mm": dict(pipe=(0.090, 0.008))}
    for lab, kw in variants.items():
        sim, res = run(xian(**kw))
        T26 = float(C(res.at(26 * DAY, "T_out")))
        T30 = float(C(res.at(30 * DAY, "T_out")))
        mean_kW = res.mean_over(0, 30 * DAY, "Q_wellhead") / 1e3
        out[lab] = dict(T_out_day26=T26, T_out_day30=T30, mean_kW=mean_kW)
        if lab == "base":
            out["series"] = dict(t_days=res.t_days, T_out=C(res.ts["T_out"]), Q_kW=res.ts["Q_wellhead"] / 1e3)
        print(f"Xi'an {lab:14s} Tout day26-30 {T26:6.2f}-{T30:6.2f} C  mean {mean_kW:6.1f} kW")
    out["measured"] = MEASURED["xian"]
    save("F_xian", out)


def run_weggis():
    out = {}
    variants = {"base": {}, "k 2.2": dict(k_rock=2.2), "k 3.0": dict(k_rock=3.0),
                "BHT 73 C": dict(bht=73.0), "BHT 78 C": dict(bht=78.0),
                "prior load 0.3": dict(prior_load=0.3), "prior load 0.7": dict(prior_load=0.7)}
    t0 = 7 * 365.25 * DAY
    for lab, kw in variants.items():
        sim, res = run(weggis(**kw))
        y1 = res.mean_over(t0, t0 + 365.25 * DAY, "Q_wellhead") * 8766.0 / 1e6     # MWh/yr
        y2 = res.mean_over(t0 + 365.25 * DAY, t0 + 2 * 365.25 * DAY, "Q_wellhead") * 8766.0 / 1e6
        Tout = res.mean_over(t0, t0 + 2 * 365.25 * DAY, "T_out") - 273.15
        out[lab] = dict(MWh_y1=y1, MWh_y2=y2, T_out_mean=Tout)
        print(f"Weggis {lab:16s} {y1:6.1f} / {y2:6.1f} MWh/yr   mean Tout {Tout:6.2f} C")
    out["measured"] = MEASURED["weggis"]
    save("F_weggis", out)


def run_bounds():
    """Conduction upper bound: annulus held at the (lowest) inlet temperature over the
    whole well (very high flow, riser perfectly insulated).  No conduction-only
    design with the published rock properties can extract more heat than this."""
    out = {}
    # Xi'an: 30 days at 17.3 C
    d = xian()
    d["operation"]["schedule"] = [dict(duration_days=30.0, mdot_kg_s=400.0, T_in_C=17.3)]
    d["well"]["tubing"][0]["layers"] = [dict(od_m=0.110, k=1e-6)]
    d["options"] = dict(include_pressure=False)
    sim, res = run(d)
    out["xian"] = dict(bound_kW_day30=res.at(30 * DAY, "Q_rock") / 1e3,
                       bound_kW_mean=res.mean_over(0, 30 * DAY, "Q_rock") / 1e3, measured_mean_kW=286.4)
    print("Xi'an bound:", out["xian"])
    # Songyuan: 60 days at the day-60 inlet 29.1 C (lower than the actual history -> higher bound)
    d = songyuan()
    d["operation"]["schedule"] = [dict(duration_days=60.0, mdot_kg_s=400.0, T_in_C=29.1)]
    d["well"]["tubing"][0]["layers"] = [dict(od_m=0.110, k=1e-6)]
    d["numerics"]["dt_max_days"] = 30.0
    d["options"] = dict(include_pressure=False)
    sim, res = run(d)
    out["songyuan"] = dict(bound_W_m_day60=res.at(60 * DAY, "Q_rock") / 2044.0,
                           bound_W_m_day1=res.at(1 * DAY, "Q_rock") / 2044.0, measured_W_m_day60=238.0)
    print("Songyuan bound:", out["songyuan"])
    save("F_bounds", out)


if __name__ == "__main__":
    which = sys.argv[1:] or ["qingdao", "songyuan", "xian", "weggis"]
    for w in which:
        dict(qingdao=run_qingdao, songyuan=run_songyuan, xian=run_xian, weggis=run_weggis, bounds=run_bounds)[w]()
