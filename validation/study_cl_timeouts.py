"""The two multilateral-loop cases of the first one-at-a-time sweep (type 8 at 100 kg/s) that gave no result within the time limit of
that batch: laterals and verticals of 6 1/8 in at 100 kg/s, and 125 kg/s through 8 1/2 in. Two questions:

    python study_cl_timeouts.py limit [seconds]   -> results/RV7_type8_timeouts.json
        run again with the inlet pressure of the study (2 MPa) and a longer limit: are they only slow?
    python study_cl_timeouts.py pressure          -> results/RV7_type8_inlet_pressure.json
        run with 10 and 20 MPa at the injection wellhead and record the pressure that arrives at the production wellhead: the loop of the
        study is one bore throughout, so the two verticals carry the whole flow, and where their friction exceeds the buoyancy of the
        loop plus the inlet pressure the water cannot reach the surface as liquid at the stated flow.
    python study_cl_timeouts.py loads             -> adds to results/RV7_type8_inlet_pressure.json the circulation power that follows from
        the stored losses (loss + 0.3 MPa of plant, pump efficiency 0.7); no solver run
    python study_cl_timeouts.py flows             -> results/RV7_type8_flows.json
        the same loop balance at 20, 30 and 50 kg/s with 2 MPa at the inlet: the flows of the sensitivity table (Table 6)
"""
from __future__ import annotations

import signal
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed

import CoolProp.CoolProp as CP

from common import save
import study_cl as cl
from aharrsim.systems import uloop, YEAR

CASES = {"oat|d|low": dict(d=0.156), "oat|mdot|high": dict(mdot=125.0), "base (100 kg/s)": {},
         "flow 20 kg/s": dict(mdot=20.0), "flow 30 kg/s": dict(mdot=30.0), "flow 50 kg/s": dict(mdot=50.0)}
DP_PLANT, ETA_P = 0.3, 0.7          # MPa and pump efficiency of the loop balance (Eq. 6 of the manuscript)


def circ_MWe(loss_MPa, mdot, T_in_C):
    """Circulation power [MWe] for a loop that loses loss_MPa between its two wellheads: (loss + plant) m / (rho eta), zero if the loop has pressure to spare."""
    rho_s = CP.PropsSI("D", "T", T_in_C + 273.15, "P", 2.0e6, "Water")
    return max(loss_MPa + DP_PLANT, 0.0) * mdot / rho_s / ETA_P


class Stop(BaseException):
    pass


def _stop(signum, frame):
    raise Stop()


def one(args):
    tag, p_inj, limit = args
    P = dict(cl.BASE[8]); P.update(CASES[tag])
    t0 = time.time()
    try:
        signal.signal(signal.SIGALRM, _stop); signal.setitimer(signal.ITIMER_REAL, limit)
        rock = cl.rock_layers(P["k"], P["rho"], P["cp"], P["L1"] + 500.0, P["grad"], kT=P["kT"])
        ch = uloop("cl8", P["L1"], P["L2"], P["d"], P["mdot"], P["T_in"], p_inj, rock, 15.0, dict(gradient_C_per_km=P["grad"]), 30.0, n_lat=int(P["n_lat"]),
                   layers=None, roughness=1e-4, dz=50.0, lat_far_field="adiabatic", lat_r_far=P["spacing"] / 2.0)
        signal.setitimer(signal.ITIMER_REAL, 0)
        r = cl.metrics(ch["producer"][1], P["mdot"], P["T_in"])
        for y in (1, 30):
            T = float(ch["T_wh"](y * YEAR)); p = float(ch["p_wh"](y * YEAR))
            r[f"p_wh_{y}_MPa"] = p / 1e6
            r[f"p_sat_wh_{y}_MPa"] = CP.PropsSI("P", "T", T, "Q", 0, "Water") / 1e6
        r.update(tag=tag, P=P, p_inj_MPa=p_inj / 1e6, loop_loss_1_MPa=p_inj / 1e6 - r["p_wh_1_MPa"], loop_loss_30_MPa=p_inj / 1e6 - r["p_wh_30_MPa"], wall_s=time.time() - t0)
        for y in (1, 30):
            r[f"circ_MWe_{y}"] = circ_MWe(r[f"loop_loss_{y}_MPa"], P["mdot"], P["T_in"])
        return r
    except Stop:
        return dict(tag=tag, P=P, p_inj_MPa=p_inj / 1e6, failed=f"no result in {limit:.0f} s", wall_s=time.time() - t0)
    except Exception as e:                                           # noqa
        signal.setitimer(signal.ITIMER_REAL, 0)
        return dict(tag=tag, P=P, p_inj_MPa=p_inj / 1e6, failed=f"{type(e).__name__}: {e}", tb=traceback.format_exc()[-1500:], wall_s=time.time() - t0)


def run(jobs, name):
    out = []
    with ProcessPoolExecutor(2) as ex:
        futs = {ex.submit(one, j): j for j in jobs}
        for f in as_completed(futs):
            r = f.result(); out.append(r)
            brief = r.get("failed") or {k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items() if isinstance(k, str) and k.startswith(("p_wh", "p_sat", "loop_loss"))} | \
                {k: round(v, 3) for k, v in r.get(30, {}).items() if isinstance(v, float)}
            print(r["tag"], f"inlet {r['p_inj_MPa']:.0f} MPa ->", brief, "| wall %.0f s" % r["wall_s"], flush=True)
            save(name, sorted(out, key=lambda q: (q["tag"], q["p_inj_MPa"])))


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "limit"
    if cmd == "limit":
        limit = float(sys.argv[2]) if len(sys.argv) > 2 else 2400.0
        run([(t, 2.0e6, limit) for t in ("oat|d|low", "oat|mdot|high")], "RV7_type8_timeouts")
    elif cmd == "pressure":
        run([("base (100 kg/s)", 2.0e6, 900.0)] + [(t, p, 900.0) for p in (10.0e6, 20.0e6) for t in ("oat|d|low", "oat|mdot|high", "base (100 kg/s)")], "RV7_type8_inlet_pressure")
    elif cmd == "flows":                             # the loop balance at the flows of the sensitivity table, 2 MPa at the inlet
        run([(t, 2.0e6, 900.0) for t in ("flow 20 kg/s", "flow 30 kg/s", "flow 50 kg/s")], "RV7_type8_flows")
    elif cmd == "loads":                             # circulation power from the stored losses of the pressure runs (no solver run)
        from common import load
        rows = load("RV7_type8_inlet_pressure")
        for r in rows:
            for y in (1, 30):
                if f"loop_loss_{y}_MPa" in r:
                    r[f"circ_MWe_{y}"] = circ_MWe(r[f"loop_loss_{y}_MPa"], r["P"]["mdot"], r["P"]["T_in"])
            print(r["tag"], f"inlet {r['p_inj_MPa']:.0f} MPa:", r.get("failed") or {k: round(r[k], 3) for k in r if k.startswith(("loop_loss", "circ_MWe"))})
        save("RV7_type8_inlet_pressure", rows)
