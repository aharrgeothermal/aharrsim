"""Rock quality: the fractured systems (types 10, 14, 14b, 19) with uneven flow between the fractures.

Three scenarios, set by the Lorenz coefficient of the flow per fracture and by the scale over which the flow varies
(aharrsim.reservoir.UnevenFractures: conduction across the rock between fractures that carry unequal flows):

  good      L_c = 0.2, from fracture to fracture     well-connected network, even sweep
  moderate  L_c = 0.5, corridors 50 m wide           several corridors take most of the flow, partial sweep
  poor      L_c = 0.8, single fractures              a few dominant fractures, early breakthrough, bypassed rock

    python study_rock.py examples   -> results/RQ_examples.json    the four examples under each scenario
    python study_rock.py needs      -> results/RQ_needs.json       the designs that hold their quoted capacity with even flow
    python study_rock.py volume     -> results/RQ_volume.json      volume that holds 90 % of the temperature rise for 30 years
    python study_rock.py scale      -> results/RQ_scale.json       the same flows spread over corridors of increasing width
    python study_rock.py red        -> results/RQ_red.json         the Project Red record read as uneven flow
    python study_rock.py red_curve  -> adds fine curves of the fitted reading to results/RQ_red.json
    python study_rock.py chain      -> results/RQ_chain.json       full injector-reservoir-producer chain for one case
"""
from __future__ import annotations

import math
import sys
import time

import numpy as np

from common import save, load
from aharrsim.reservoir import ParallelFractures, UnevenFractures, ZonedFractures, zoned_dense_limit, flow_pattern, lorenz_of, lorenz_curve, cbar_w
from study_open import SYS, design, strip, geometry, YRS, cbar

YEAR = 365.25 * 86400.0
SCEN = dict(good=dict(lorenz=0.2, corridor_m=0.0), moderate=dict(lorenz=0.5, corridor_m=50.0), poor=dict(lorenz=0.8, corridor_m=0.0))
SEEDS = (0, 1, 2, 3, 4)
KEEP = ("T_res", "net", "gross", "T_sd_K", "net1", "net30", "mean_net", "gross1", "T1", "T30", "GWh", "W_circ", "t_sweep_yr", "x_30", "y_30", "V_km3", "N", "m", "S_m", "t_10pct_yr")


def run(s, scen=None, **kw):
    u = None if scen is None else dict(SCEN[scen], seeds=SEEDS)
    d = design(s, uneven=u, **kw)
    return {k: d[k] for k in KEEP}


def examples():
    out = dict(years=YRS.tolist(), scenarios=SCEN, seeds=list(SEEDS), systems={})
    for sid, s in SYS.items():
        t0 = time.time()
        o = dict(name=s["name"], even=run(s))
        for sc in SCEN:
            o[sc] = run(s, sc)
        out["systems"][str(sid)] = o
        print(sid, {k: (round(o[k]["T30"], 1), round(o[k]["mean_net"], 2), round(o[k]["GWh"])) for k in ("even", *SCEN)}, round(time.time() - t0), "s", flush=True)
        save("RQ_examples", out)


def need_designs():
    scr = load("ST2_screen")
    res = {}
    for sid, s in SYS.items():
        for n in scr["systems"][str(sid)].get("need_fine", []):
            if n["design"]:
                d = n["design"]
                res[sid] = dict(target=n["target"], N=int(d["N"]), m=float(d["m"]), W=float(d["W"]), L=float(d["L"]), vf=float(d["vf"]))
                break
    return res


def needs():
    out = dict(years=YRS.tolist(), scenarios=SCEN, seeds=list(SEEDS), systems={})
    for sid, nd in need_designs().items():
        s = SYS[sid]; t0 = time.time()
        kw = dict(N=nd["N"], m=nd["m"], W=nd["W"], L=nd["L"])
        o = dict(name=s["name"], design=nd, even=run(s, **kw))
        for sc in SCEN:
            o[sc] = run(s, sc, **kw)
        out["systems"][str(sid)] = o
        print(sid, nd["target"], {k: (round(o[k]["T30"], 1), round(o[k]["mean_net"], 2), round(o[k]["net30"], 2)) for k in ("even", *SCEN)}, round(time.time() - t0), "s", flush=True)
        save("RQ_needs", out)


def theta30(s, vf, m, scen, sid, seeds=SEEDS, spacing=10.0):
    """Retained share of the temperature rise after 30 years for vf times the example's volume, fractures `spacing` apart."""
    W, L = geometry(sid, vf)
    N = max(int(W / spacing), 2); Q = m / s["rho_w"]; wings = s.get("wings", 1); cw = cbar(s)
    args = dict(Q=Q, N=N * wings, L_flow=L, B=s["B"], S=W / N, k_r=s["k"], rhoc_r=s["rhoc"], T_rock_C=s["T_rock"], T_inj_C=s["T_inj"], rhoc_w=cw * s["rho_w"])
    t = np.array([30.0 * YEAR])
    if scen is None:
        return 1.0 - float(ParallelFractures(**args).T_D(t)[0])
    sc = SCEN[scen]; cor = max(int(round(sc["corridor_m"] / (W / N))), 1)
    return 1.0 - float(np.mean([UnevenFractures(**args, lorenz=sc["lorenz"], corridor=cor, seed=sd, n_pattern=min(N, 96)).T_D(t)[0] for sd in seeds]))


def volume():
    """Volume at which the outlet keeps 90 % of its temperature rise for 30 years, at the flow of the design that holds the quoted capacity."""
    out = dict(scenarios=SCEN, systems={})
    for sid, nd in need_designs().items():
        s = SYS[sid]; o = dict(design=nd)
        for sc in (None, *SCEN):
            lo, hi = nd["vf"] * 0.3, nd["vf"] * 400.0
            if theta30(s, hi, nd["m"], sc, sid) < 0.9:
                o[sc or "even"] = None; continue
            for _ in range(22):
                mid = math.sqrt(lo * hi)
                if theta30(s, mid, nd["m"], sc, sid) >= 0.9:
                    hi = mid
                else:
                    lo = mid
            W, L = geometry(sid, hi)
            V = W * L * s["B"] * s.get("wings", 1) / 1e9
            o[sc or "even"] = dict(vf=hi, V_km3=V, W=W, L=L, x_30=30.0 / (V * 1e9 * s["rhoc"] / (nd["m"] * cbar(s)) / YEAR))
            print(sid, sc or "even", round(V, 3), "km3  x30", round(o[sc or "even"]["x_30"], 3), flush=True)
        ev = o["even"]["V_km3"]
        for sc in SCEN:
            if o[sc]:
                o[sc]["multiple"] = o[sc]["V_km3"] / ev
        out["systems"][str(sid)] = o
        save("RQ_volume", out)
    # closely spaced fractures in wide zones: closed form
    from scipy.stats import norm
    from aharrsim.reservoir import lorenz_sigma
    out["dense_segregated"] = {str(L): dict(x_90=float(math.exp(norm.ppf(0.1) * lorenz_sigma(L) - 0.5 * lorenz_sigma(L) ** 2)), top10=float(lorenz_curve(0.1, L))) for L in (0.2, 0.5, 0.8)}
    save("RQ_volume", out)


def scale():
    """The Cape Station-type design that holds 10 MWe with even flow: the same six flow multipliers (a lognormal allocation of the given
    Lorenz coefficient) spread over corridors of increasing width, so that only the scale of the unevenness changes."""
    nd = need_designs()[141]; s = SYS[141]
    Q = nd["m"] / s["rho_w"]; cw = cbar(s); S = nd["W"] / nd["N"]
    args = dict(Q=Q, N=nd["N"] * 2, L_flow=nd["L"], B=s["B"], S=S, k_r=s["k"], rhoc_r=s["rhoc"], T_rock_C=s["T_rock"], T_inj_C=s["T_inj"], rhoc_w=cw * s["rho_w"])
    t = YRS * YEAR; dT = s["T_rock"] - s["T_inj"]
    out = dict(years=YRS.tolist(), design=nd, spacing_m=S, n_corridors=6, even=(s["T_rock"] - ParallelFractures(**args).T_D(t) * dT).tolist(), rows=[])
    for Lc in (0.2, 0.5, 0.8):
        mult = flow_pattern(6, Lc, 1, 0)
        z = ZonedFractures(**args, flow=mult, area=np.ones(6))
        out["rows"].append(dict(lorenz=Lc, corridor_m=None, mult=sorted(mult.tolist()), T_res=(s["T_rock"] - z.T_D(t) * dT).tolist()))
        for k in (1, 2, 4, 8, 16, 32):
            TD = np.mean([UnevenFractures(**args, mult=flow_pattern(6 * k, Lc, k, sd)).T_D(t) for sd in (0, 1, 2)], axis=0)
            out["rows"].append(dict(lorenz=Lc, corridor_m=k * S, T_res=(s["T_rock"] - TD * dT).tolist()))
            print(Lc, round(k * S), np.round(out["rows"][-1]["T_res"][-1], 1), flush=True)
            save("RQ_scale", out)


def red_curve():
    """Fine curves of the uneven-flow reading of the Project Red record (Lorenz coefficient 0.35) for the figure."""
    from field_cases import RED, RHO
    fe = load("FE_field_cases")["red"]; out = load("RQ_red")
    r = RED; cw = cbar_w(r["T_inj"], r["T_rock"] - 6.0, r["p_res"]); dT = fe["dT_K"]; plateau = fe["plateau_C"]
    tt = np.linspace(0.05, 8.0, 54)
    for tag, cor in (("fracture", 1), ("stage", 6)):
        TD = np.mean([UnevenFractures(Q=r["m"] / RHO, N=r["N"], L_flow=r["L"], B=r["B"], S=r["W"] / r["N"], k_r=r["k"], rhoc_r=r["rhoc"], T_rock_C=r["T_rock"],
                                      T_inj_C=r["T_inj"], rhoc_w=cw * RHO, lorenz=0.35, corridor=cor, seed=sd).T_D(tt * YEAR) for sd in (0, 1, 2, 3, 4)], axis=0)
        out["fits"][tag]["curve"] = dict(lorenz=0.35, t_yr=tt.tolist(), T_C=(plateau - TD * dT).tolist())
        save("RQ_red", out)
    print("red curves done")


def red():
    """The Project Red record (614 days) read as uneven flow: the Lorenz coefficient that fits it, at the scale of a stage and of a fracture."""
    from field_cases import RED, red_model, RHO
    fe = load("FE_field_cases")["project_red"] if "project_red" in load("FE_field_cases") else load("FE_field_cases")["red"]
    t = np.array(fe["t_obs"]); ok = np.array(fe["used"], bool); TD_obs = np.array(fe["TD_obs"]); dT = fe["dT_K"]; plateau = fe["plateau_C"]
    r = RED; cw = cbar_w(r["T_inj"], r["T_rock"] - 6.0, r["p_res"])
    tg = np.concatenate([np.linspace(0.1, t[ok].max(), 28), [2, 3, 4, 5, 6, 8, 10, 15, 20, 30]])

    def model(Lc, cor, seeds=(0, 1, 2, 3, 4)):
        TD = np.mean([UnevenFractures(Q=r["m"] / RHO, N=r["N"], L_flow=r["L"], B=r["B"], S=r["W"] / r["N"], k_r=r["k"], rhoc_r=r["rhoc"], T_rock_C=r["T_rock"],
                                      T_inj_C=r["T_inj"], rhoc_w=cw * RHO, lorenz=Lc, corridor=cor, seed=sd).T_D(tg * YEAR) for sd in seeds], axis=0)
        return TD

    out = dict(t_yr=tg.tolist(), fits={})
    for tag, cor in (("fracture", 1), ("stage", 6)):
        rows = []
        for Lc in (0.0, 0.1, 0.2, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.7, 0.8):
            TD = model(Lc, cor) if Lc > 0 else red_model(tg, 1.0, 1.0)[0]
            rms = float(np.sqrt(np.mean((np.interp(t[ok], tg[:28], TD[:28]) - TD_obs[ok]) ** 2)) * dT)
            rows.append(dict(lorenz=Lc, rms_K=rms, T_at={str(int(y)): float(plateau - np.interp(y, tg, TD) * dT) for y in (2, 3, 5, 10, 20, 30)}, T_C=(plateau - TD * dT).tolist()))
            print(tag, Lc, round(rms, 2), {k: round(v) for k, v in rows[-1]["T_at"].items()}, flush=True)
        out["fits"][tag] = dict(corridor_fractures=cor, corridor_m=cor * r["W"] / r["N"], rows=rows)
        save("RQ_red", out)
    # the two earlier readings at the same times
    for name, (f, phi) in dict(fit_fractures=(fe["cases"]["fit_fractures"]["f"], 1.0), fit_volume=(1.0, fe["cases"]["fit_volume"]["phi"])).items():
        TD = red_model(tg, f, phi)[0]
        out[name] = dict(f=f, phi=phi, rms_K=float(np.sqrt(np.mean((np.interp(t[ok], tg[:28], TD[:28]) - TD_obs[ok]) ** 2)) * dT),
                         T_at={str(int(y)): float(plateau - np.interp(y, tg, TD) * dT) for y in (2, 3, 5, 10, 20, 30)}, T_C=(plateau - TD * dT).tolist())
    save("RQ_red", out)


def chain():
    """Full injector-reservoir-producer chains with uneven flow (one arrangement), against the screening model."""
    from study_open import confirm_one
    out = []
    for sid, c in ((141, dict(tag="Cape Station-type example, Moderate", N=150, m=100.0, vf=1, uneven=dict(SCEN["moderate"], seed=0))),
                   (10, dict(tag="Soultz-type example, Poor", N=20, m=30.0, vf=1, uneven=dict(SCEN["poor"], seed=0)))):
        r = confirm_one((sid, c)); out.append(r)
        print({k: (np.round(v, 2).tolist() if isinstance(v, list) and v and isinstance(v[0], float) else v) for k, v in r.items() if k in ("tag", "failed", "T_res", "screen_T_res", "gross", "screen_gross", "wall_s")}, flush=True)
        save("RQ_chain", out)


if __name__ == "__main__":
    dict(examples=examples, needs=needs, volume=volume, scale=scale, red=red, red_curve=red_curve, chain=chain)[sys.argv[1]]()
