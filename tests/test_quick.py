"""Quick regression tests (about 1 minute):  python -m pytest -q tests/"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from aharrsim import case_from_dict, Simulator, DAY  # noqa: E402
from aharrsim.analytic import reverse_circulation_steady, jaeger_cylinder_flux  # noqa: E402

CONST = {"constant": {"rho": 1000.0, "cp": 4180.0, "mu": 1.0e-3, "k": 0.6}}


def _steady(dz):
    d = dict(fluid=CONST,
             well=dict(depth_m=6000, casing=[dict(id_m=0.16, layers=[])],
                       tubing=[dict(id_m=0.076, layers=[dict(od_m=0.1143, k=1.0)], R_wall_K_m_W=0.2)]),
             rock=dict(layers=[dict(k=2.5, rho=2500, cp=1000)], surface_T_C=15, geotherm=dict(gradient_C_per_km=30)),
             operation=dict(pressure_bc=dict(outlet_MPa=1.0), schedule=[dict(duration_s=1.0, mdot_kg_s=10.0, T_in_C=40)]),
             numerics=dict(dz_m=dz, dt_initial_s=1.0),
             options=dict(gravity_work=False, include_pressure=False, rock_model="resistance", R_formation_K_m_W=0.5,
                          h_wall_override_W_m2K=1e12, fluid_storage=False))
    sim = Simulator(case_from_dict(d))
    sim.run()
    md = sim.s_f[:sim.Nw + 1]
    Ta, Tt = reverse_circulation_steady(md, 6000, 10.0, 4180.0, 0.5, 0.2, 40.0, 15.0, 0.03)
    Tt_m = sim.fluid.ph(sim.p_t, sim.h_t).T - 273.15
    return float(np.max(np.abs(Tt_m - Tt)))


def test_steady_analytical_second_order():
    e1, e2 = _steady(200.0), _steady(100.0)
    assert e2 < 2e-3
    assert 3.5 < e1 / e2 < 4.5


def test_jaeger_cylinder():
    d = dict(fluid=CONST,
             well=dict(depth_m=600, casing=[dict(id_m=0.2, layers=[])],
                       tubing=[dict(id_m=0.076, layers=[dict(od_m=0.1143, k=1.0)], R_wall_K_m_W=1e9)]),
             rock=dict(layers=[dict(k=2.5, rho=2500, cp=1000)], surface_T_C=100, geotherm=dict(gradient_C_per_km=0)),
             operation=dict(pressure_bc=dict(outlet_MPa=1.0), schedule=[dict(duration_days=100, mdot_kg_s=2000, T_in_C=20)]),
             numerics=dict(dz_m=100, r_far_m=60, dt_initial_s=10, dt_growth=1.15, output_times_days=[1, 100]),
             options=dict(gravity_work=False, include_pressure=False, h_wall_override_W_m2K=1e12, fluid_storage=False))
    sim = Simulator(case_from_dict(d))
    res = sim.run()
    j = int(np.argmin(np.abs(sim.s_c[:sim.Nw] - 300.0)))
    for pr in res.profiles:
        Ta = 0.5 * (pr["T_a"][j] + pr["T_a"][j + 1])
        q_ex = jaeger_cylinder_flux(pr["t"], 2.5, 1e-6, 0.1, 373.15 - Ta)
        assert abs(pr["q_rock"][j] / q_ex - 1) < 0.005


def test_energy_balance_iapws():
    d = dict(well=dict(depth_m=4000,
                       casing=[dict(id_in=6.094, layers=[dict(od_in=7.0, material="steel"), dict(od_in=8.5, material="cement")])],
                       tubing=[dict(id_in=2.992, layers=[dict(od_in=3.5, k=45), dict(od_in=4.0, k=0.03), dict(od_in=4.5, k=45)])]),
             rock=dict(layers=[dict(k=2.8, rho=2650, cp=900)], surface_T_C=15, geotherm=dict(gradient_C_per_km=30)),
             operation=dict(pressure_bc=dict(outlet_MPa=2.0), schedule=[dict(duration_days=60, mdot_kg_s=10, T_in_C=40)]),
             numerics=dict(dz_m=100, output_times_days=[60]))
    sim = Simulator(case_from_dict(d))
    res = sim.run()
    ts = res.ts
    resid = np.abs(ts["Q_wellhead"] - (ts["Q_rock"] - ts["dE_fluid"]))
    assert resid.max() / np.abs(ts["Q_wellhead"]).max() < 1e-4
    assert 45.0 < ts["T_out"][-1] - 273.15 < 60.0


def test_parallel_fractures_inversion():
    """Fracture element (Gringarten et al., 1975) against its closed-form limits.  The sharp-front check fails with the
    Gaver-Stehfest inversion of versions up to 2.1, which is kept as method='stehfest'."""
    from scipy.special import erfc
    from aharrsim.reservoir import ParallelFractures
    yr = 365.25 * 86400.0
    # isolated fractures (wide spacing): T_D = erfc[1 / (2 sqrt(t_D))]
    far = ParallelFractures(Q=0.03, N=5, L_flow=600.0, B=400.0, S=2000.0, k_r=2.8, rhoc_r=2.4e6, T_rock_C=200, T_inj_C=60)
    t = np.array([1.0, 3.0, 10.0, 30.0]) * yr
    assert np.max(np.abs(far.T_D(t) - erfc(0.5 / np.sqrt(t / far.gamma2)))) < 1e-6
    # thin slabs (fractures 2.5 m apart, far inside the conduction length): a sharp front at the sweep time
    pf = ParallelFractures(Q=0.063, N=390, L_flow=880.0, B=143.0, S=2.5, k_r=2.8, rhoc_r=2.25e6, T_rock_C=190, T_inj_C=60, rhoc_w=3.88e6)
    ts = pf.t_sweep
    x, y = pf.groups(ts)
    assert abs(float(x) - 1.0) < 1e-12 and float(y) < 0.05
    assert float(pf.T_D(0.85 * ts)[0]) < 1e-3                       # still at rock temperature
    assert float(pf.T_D(0.85 * ts, method="stehfest")[0]) > 0.05    # the old inversion had already started the decline
    assert abs(float(pf.T_D(ts)[0]) - 0.5) < 0.02                   # half-way at the front
    assert float(pf.T_D(1.15 * ts)[0]) > 0.999
    # all the heat of the block comes out: the integral of 1 - T_D over time is the sweep time
    tt = np.linspace(1e-3, 1.6, 1601) * ts
    td = pf.T_D(tt)
    kept = 1.0 - td
    assert abs(float(np.sum(0.5 * (kept[1:] + kept[:-1]) * np.diff(tt))) / ts - 1.0) < 2e-3
    assert np.diff(td).min() > -1e-5                                # monotonic, no ringing


def test_single_well_dipole_element():
    """Two zones on one well in a large block: first arrival at 0.746 D^3 (rho c)_b / (Q (rho c)_w), half the rise lost at 11.6, the
    long tail (c / x)^(1/4), and the fractured form tending to the sharp-front form when the fractures are close."""
    from aharrsim.reservoir import SingleWellDipole, dipole_arrival
    d = SingleWellDipole(Q=1.0, D=1.0, rhoc_b=1.0, T_rock_C=100.0, T_inj_C=0.0, rhoc_w=1.0)
    assert abs(d.t_bt - 0.746012) < 1e-6
    assert float(d.T_D(0.74)[0]) == 0.0 and float(d.T_D(0.80)[0]) > 0.01
    assert abs(float(d.T_D(float(dipole_arrival(0.5)))[0]) - 0.5) < 1e-6 and abs(float(dipole_arrival(0.5)) - 11.582) < 0.01
    x = np.array([2.0, 5.0, 20.0, 200.0])
    assert np.max(np.abs((1.0 - d.T_D(x)) - (np.interp(d.T_D(x), [0, 1], [0.746, 0.718]) / x) ** 0.25)) < 2e-3
    close = SingleWellDipole(Q=1.0, D=1.0, rhoc_b=1.0, T_rock_C=100.0, T_inj_C=0.0, rhoc_w=1.0, S=0.02, k_r=1.0)
    assert np.max(np.abs(close.T_D(x) - d.T_D(x))) < 5e-3
    wide = SingleWellDipole(Q=1.0, D=1.0, rhoc_b=1.0, T_rock_C=100.0, T_inj_C=0.0, rhoc_w=1.0, S=2.0, k_r=1.0)
    assert float(wide.T_D(0.5)[0]) > 0.03                      # fractures far apart cool before the sweep time
    assert abs(d.T_out(273.15, 11.582) - 323.15) < 0.05         # half-way between 0 and 100 C


def test_porous_block():
    """Block model: energy balance to round-off, no flow in a block at rest, and the unbounded source-sink pair of the element
    reproduced on a coarse grid."""
    from aharrsim.porous import AxisymmetricBlock, water_table
    from aharrsim.reservoir import SingleWellDipole
    W = water_table(2.0e5)
    b = AxisymmetricBlock(R=3.0, H=6.0, r_w=0.01, K=1e-9, phi=0.35, rhoc_s=1.2e6 / 0.65, lam=0.0, T0=60.0, nr=30, nz=120, r1=0.06,
                          buoyancy=False, viscosity=False, T_ref=40.0, water=W)
    Q = 1.0 / 3600.0
    b.add_screen(2.45, 2.55, -Q); b.add_screen(3.45, 3.55, +Q, T_in=20.0)
    el = SingleWellDipole(Q=Q, D=1.0, rhoc_b=b.rhoc_b, T_rock_C=60.0, T_inj_C=20.0, rhoc_w=b.rhoc_w)
    r = b.run(2.0 * el.t_ref, dt_out=el.t_ref / 20.0)
    assert np.max(np.abs(r["E_out"] + r["E_block"])) < 1e-9 * np.max(np.abs(r["E_out"]))
    TD = (60.0 - r["T_prod"]) / 40.0
    assert np.max(np.abs(TD - el.T_D(r["t"]))) < 0.03
    # a block at one temperature is at rest, with or without buoyancy
    c = AxisymmetricBlock(R=0.5, H=1.5, r_w=0.008, K=1e-9, phi=0.35, rhoc_s=1.2e6 / 0.65, lam=1.5, T0=60.0, nr=12, nz=30, water=W)
    c.solve_flow()
    assert np.max(np.abs(c.FZ)) < 1e-15 and np.max(np.abs(c.FR)) < 1e-15
