"""Study parts 2 and 3: open systems (flow through fractured rock): fracture density, tortuosity, channelling, flow and
stimulated volume for the EGS examples (types 10, 14, 19 of dev_types.py), and their sensitivity to the uncertainties
that dominate the closed loops.

Screening model: Gringarten parallel fractures (aharrsim.reservoir.ParallelFractures) + a calibrated wellbore enthalpy
loss + the plant models; selected designs are confirmed with the full injector -> reservoir -> producer chain.

  stimulated block  W (along the wells) x L (between the wells) x B (height); N fractures, spacing S = W / N
  tortuosity tau    flow path tau*L in each fracture, wetted area x tau, rock volume kept (slab thickness S / tau)
  active fraction f only f*N fractures take flow (rock volume kept: spacing S / f)
  sweep phi         only phi of the block is in thermal contact (rock volume lost)
  impedance         I = I0 [ a + (1 - a) (N0 / (f N)) tau ],  a = near-well share (0.5)
  circulation       closed-loop balance: pump dp = I Q - thermosiphon head + wellbore friction + plant drop (0.3 MPa);
                    the static columns and friction are anchored on the full-chain examples of dev_types.py (hydraulics())

    python study_open.py screen     -> results/ST2_screen.json
    python study_open.py confirm    -> results/ST2_confirm.json   (full chains)
    python study_open.py part3      -> results/ST3_open_uncertainty.json
"""
from __future__ import annotations

import math
import sys
import time

import numpy as np
import CoolProp.CoolProp as CP

from common import save, load
from functools import lru_cache

from aharrsim.reservoir import ParallelFractures, UnevenFractures, cbar_w
from aharrsim.systems import doublet, YEAR
from aharrsim import plant
from study_cl import k_of_T, tri

IN = 0.0254
G = 9.81
YRS = np.array([0.25, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30], float)

SYS = {
    10: dict(name="Soultz-type EGS doublet (type 10)", label="10", T_rock=200.0, depth=5000.0, W=1000.0, L=600.0, B=400.0, N0=20, m=30.0, T_inj=60.0, k=3.0,
             rhoc=2650 * 850.0, I0=0.3e9, p_wh=2.0e6, p_bot=47.5e6, c_loss=13e3, lift_kW_per_kgs=0.0, quoted=(3.0, 5.0), rho_w=900.0, T_s=15.0),
    # Project Red: stimulated volume as published by Fervo (April 2026): 975 m of lateral x 110 m between wells x 143 m high; 88 fractures (75-100 estimated)
    14: dict(name="Project Red-type horizontal doublet (type 14)", label="14", T_rock=190.0, depth=2300.0, W=975.0, L=110.0, B=143.0, N0=88, m=57.96, T_inj=60.0, k=2.8,
             rhoc=2650 * 850.0, I0=0.15e9, p_wh=2.0e6, p_bot=18.1e6, c_loss=13e3, lift_kW_per_kgs=0.0, quoted=(3.5, 3.5), rho_w=920.0, T_s=15.0, stretch=4),
    # Cape Station (Phase I pattern, 2024 test): producer between two injectors (two fracture wings), 1,524 m laterals 137 m apart; I0 from the reported
    # wellhead pressures (15 MPa injection, 2.07 MPa production): dp 13.2 MPa at 100 kg/s
    141: dict(name="Cape Station-type producer between two injectors (type 14b)", label="14b", T_rock=200.0, depth=2590.0, W=1524.0, L=137.0, B=143.0, N0=150, m=100.0,
              T_inj=60.0, k=2.8, rhoc=2650 * 850.0, I0=0.0894e9, p_wh=2.07e6, p_bot=28.0e6, c_loss=13e3, lift_kW_per_kgs=0.0, quoted=(10.0, 16.0), rho_w=900.0,
              T_s=15.0, wings=2, stretch=2),
    19: dict(name="Superhot EGS doublet (type 19)", label="19", T_rock=375.0, depth=3500.0, W=800.0, L=300.0, B=400.0, N0=20, m=40.0, T_inj=60.0, k=2.5,
             rhoc=2700 * 950.0, I0=0.3e9, p_wh=5.0e6, p_bot=30.0e6, c_loss=7e3, lift_kW_per_kgs=0.0, quoted=(50.0, 50.0), rho_w=900.0, T_s=15.0),
}


for _k, _v in SYS.items():
    _v["sid"] = _k


@lru_cache(maxsize=4096)
def _cbar(T_inj, T_rock, p):
    return cbar_w(T_inj, T_rock, p)


def cbar(s, T_rock=None):
    """Mean specific heat of the water between injection and rock temperature [J/kg/K] (see aharrsim.reservoir.cbar_w)."""
    return _cbar(float(s["T_inj"]), round(float(s["T_rock"] if T_rock is None else T_rock), 1), float(s["p_bot"]))


DT_ID = {10: 10.0, 14: 14.0, 141: 14.5, 19: 19.0}
_hyd = {}


def hydraulics(sid):
    """Anchors for the closed-loop pressure balance of system sid, from the full-chain example in results/DT_types.json:
    static weight of the injection column, static weight of the production column against reservoir outlet temperature,
    and the single-phase friction of both wells at the example's flow."""
    if sid not in _hyd:
        r = {float(x["id"]): x for x in load("DT_types")}[DT_ID[sid]]
        a = r["assumptions"]; lp = a["loop"]; fr = a.get("friction_MPa", {})
        fi, fp = fr.get("injector", 0.0), fr.get("producer", 0.0)          # two-phase producers: friction stays inside the static table
        T = np.array([r["T_res"][y] for y in sorted(lp, key=float)]); cp_ = np.array([lp[y]["col_prod_MPa"] - fp for y in sorted(lp, key=float)])
        o = np.argsort(T); T, cp_ = T[o], cp_[o]
        keep = np.concatenate(([True], np.diff(T) > 0.5))
        h = dict(static_inj=lp["1"]["col_inj_MPa"] + fi, T=T[keep], static_prod=cp_[keep], fr_inj0=fi, fr_prod0=fp, m0=float(a["mdot"]), dm=[0.0], dcol=[0.0])
        if fp == 0.0:
            # flashing producer: its column (gravity + two-phase friction) does not scale like single-phase friction, so the change with
            # flow is taken from the full-chain confirmation runs (results/ST2_confirm.json), when they exist
            try:
                pts = {}
                for c in load("ST2_confirm"):
                    if c.get("sid") == sid and "failed" not in c:
                        m = c["m"]; big = m > 65.0
                        fri = fi * (m / h["m0"]) ** 1.8 * ((6.094 / 8.535) ** 4.8 if big else 1.0)
                        col = h["static_inj"] - fri - c["thermosiphon_MPa"][0]          # producer BHP - WHP at year 1 (rock temperature)
                        pts[m] = col - float(np.interp(c["T_res"][0], h["T"], h["static_prod"]))
                ms = sorted(pts)
                if len(ms) > 1:
                    dm = [0.0, 65.0] + [m for m in ms if m > 65.0]; dc = [0.0, 0.0] + [pts[m] for m in ms if m > 65.0]
                    if len(dm) > 3:                                                    # linear extrapolation beyond the largest confirmed flow
                        dm.append(2.0 * dm[-1]); dc.append(dc[-1] + (dc[-1] - dc[-2]) / (dm[-2] - dm[-3]) * (dm[-1] - dm[-2]))
                    h["dm"], h["dcol"] = dm, dc
            except Exception:
                pass
        _hyd[sid] = h
    return _hyd[sid]


def circulation_MW(sid, s, m, T_res, I, W=None, eta=0.7, dp_plant=0.3):
    """Circulation power [MW] at reservoir outlet temperatures T_res (array): closed-loop balance with the anchored columns."""
    h = hydraulics(sid)
    T = np.asarray(T_res, float)
    sp = np.interp(T, h["T"], h["static_prod"])
    hot = T > h["T"][-1]                     # hotter than the example ever was: scale the column weight with liquid density (liquid producers only)
    if np.any(hot) and h["fr_prod0"] > 0:
        p_mid = 0.5 * (s["p_bot"] + s["p_wh"])
        rho0 = CP.PropsSI("D", "T", h["T"][-1] + 273.15, "P", p_mid, "Water")
        sp = np.where(hot, h["static_prod"][-1] * np.array([CP.PropsSI("D", "T", t + 273.15, "P", p_mid, "Water") if hh else rho0 for t, hh in zip(T, hot)]) / rho0, sp)
    big = m > (110.0 if sid == 141 else 65.0)              # 9-5/8 in casing above this flow (as in the full-chain confirmation)
    scale = (m / h["m0"]) ** 1.8 * ((6.094 / 8.535) ** 4.8 if big else 1.0)
    if W is not None and "stretch" in s:                    # horizontal wells: friction length grows with the lateral
        lat0 = s["W"] / (3.0 if sid == 141 else 1.0); lat = W / (3.0 if sid == 141 else 1.0)
        scale *= (s["depth"] + lat) / (s["depth"] + lat0)
    sp = sp + float(np.interp(m, h["dm"], h["dcol"]))
    dp = I * (m / s["rho_w"]) / 1e6 - (h["static_inj"] - sp) + (h["fr_inj0"] + h["fr_prod0"]) * scale + dp_plant
    rho_s = 985.0
    return np.maximum(dp, 0.0) * 1e6 * m / rho_s / eta / 1e6, dp


_tab = {}


def _specific_work(h, p, T_amb):
    """Gross specific work [MJ/kg] of the best plant for a stream (h, p)."""
    T = CP.PropsSI("T", "H", h, "P", p, "Water") - 273.15
    w = [0.0]
    if T > 60:
        try:
            w.append(plant.orc_practical(1.0, h, p, T_amb)["W_net"] / 1e6)
        except Exception:
            pass
    if T > 150:
        for dbl in (False, True):
            try:
                w.append(plant.flash(1.0, h, p, T_cond_C=T_amb + 25.0, double=dbl)["W_net"] / 1e6)
            except Exception:
                pass
    return max(w)


def work_table(p, T_amb):
    """Interpolation table of specific work against enthalpy (plant output is proportional to mass flow)."""
    key = (p, T_amb)
    if key not in _tab:
        try:
            disk = load("ST_work_tables")
        except Exception:
            disk = {}
        sk = f"{p:.0f}|{T_amb:.0f}"
        if sk not in disk:
            hs = np.linspace(2.0e5, 2.9e6, 46)
            disk[sk] = dict(h=hs.tolist(), w=[_specific_work(h, p, T_amb) for h in hs])
            save("ST_work_tables", disk)
        _tab[key] = (np.array(disk[sk]["h"]), np.array(disk[sk]["w"]))
    return _tab[key]


def gross_MWe(m, T_res_C, s, T_amb=15.0):
    """Gross MWe from m kg/s leaving the reservoir at T_res: wellbore enthalpy loss (gravity + conduction), then the best plant."""
    T_res_C = max(T_res_C, s["T_inj"] + 1.0)
    h_res = CP.PropsSI("H", "T", T_res_C + 273.15, "P", s["p_bot"], "Water")
    T_mid = s["T_s"] + 0.5 * (s["T_rock"] - s["T_s"])
    q_loss = s["c_loss"] * max(T_res_C - T_mid, 0.0)
    h = h_res - G * s["depth"] - q_loss / m
    hs, ws = work_table(s["p_wh"], T_amb)
    return m * float(np.interp(h, hs, ws))


def design(s, N=None, tau=1.0, f=1.0, phi=1.0, m=None, W=None, L=None, B=None, T_rock=None, k=None, rhoc=None, T_amb=15.0, a_nw=0.5, I0=None,
           years=YRS, uneven=None):
    """Evaluate one fractured-reservoir design with the screening model. Returns histories and 30-year figures.
    uneven = dict(lorenz, corridor_m, seeds): unequal flows between the fractures (aharrsim.reservoir.UnevenFractures), the
    response averaged over the seeded arrangements; the impedance is left as it is (it is measured at the wellheads)."""
    N = N or s["N0"]; m = m or s["m"]; W = W or s["W"]; L = L or s["L"]; B = B or s["B"]
    T_rock = s["T_rock"] if T_rock is None else T_rock; k = k or s["k"]; rhoc = rhoc or s["rhoc"]; I0 = s["I0"] if I0 is None else I0
    Q = m / s["rho_w"]
    wings = s.get("wings", 1)                  # 2: each fracture flows towards an injector on either side of the producer
    cw = cbar(s, T_rock)
    Na = max(f * N, 1.0)                       # active fractures
    S_eff = W / Na / tau                       # slab thickness per active fracture with the rock volume kept
    pf = ParallelFractures(Q=Q, N=Na * wings, L_flow=tau * L, B=phi * B, S=S_eff, k_r=k, rhoc_r=rhoc, T_rock_C=T_rock, T_inj_C=s["T_inj"],
                           rhoc_w=cw * s["rho_w"])
    if uneven and uneven.get("lorenz", 0.0) > 0.0:
        n_pat = int(min(max(round(Na), 2), uneven.get("n_max", 160)))
        cor = max(int(round(uneven.get("corridor_m", 0.0) / S_eff)), 1)
        TDs = np.array([UnevenFractures(Q=Q, N=Na * wings, L_flow=tau * L, B=phi * B, S=S_eff, k_r=k, rhoc_r=rhoc, T_rock_C=T_rock, T_inj_C=s["T_inj"],
                                        rhoc_w=cw * s["rho_w"], lorenz=uneven["lorenz"], corridor=cor, seed=sd, n_pattern=n_pat).T_D(years * YEAR)
                        for sd in uneven.get("seeds", (0, 1, 2))])
        TD = TDs.mean(axis=0); TD_sd = TDs.std(axis=0)
    else:
        TD = pf.T_D(years * YEAR); TD_sd = np.zeros_like(TD)
    TD = np.maximum.accumulate(np.minimum.accumulate(TD[::-1])[::-1])     # the step response is monotone (guards against inversion noise)
    T_res = T_rock - TD * (T_rock - s["T_inj"])
    ss = dict(s); ss["T_rock"] = T_rock
    gross = np.array([gross_MWe(m, T, ss, T_amb) for T in T_res])
    I = I0 * (a_nw + (1 - a_nw) * (s["N0"] / Na) * tau)
    W_circ_t, dp_t = circulation_MW(s["sid"], s, m, T_res, I, W=W)
    W_circ = float(np.interp(1.0, years, W_circ_t))            # year-1 circulation power (it rises as the production column cools)
    W_lift = s["lift_kW_per_kgs"] * m / 1e3
    net = np.maximum(gross - W_circ_t - W_lift, 0.0)
    mean_net = float(np.trapezoid(net, years) / (years[-1] - years[0]))
    mean_gross = float(np.trapezoid(gross, years) / (years[-1] - years[0]))
    V = W * L * B * phi * wings
    t_sweep = V * rhoc / (m * cw) / YEAR       # front-crossing time of the block [yr]: V (rho c)_r / (mdot cbar_w)
    drop = np.where(TD > 0.1)[0]
    t10 = float(years[drop[0]]) if len(drop) else None
    return dict(N=N, tau=tau, f=f, phi=phi, m=m, W=W, L=L, B=B, T_rock=T_rock, V_km3=V / 1e9, S_m=W / N, T_res=T_res.tolist(), gross=gross.tolist(), net=net.tolist(),
                T30=float(T_res[-1]), T1=float(np.interp(1.0, years, T_res)), gross1=float(np.interp(1.0, years, gross)), gross30=float(gross[-1]),
                net1=float(np.interp(1.0, years, net)), net30=float(net[-1]), mean_net=mean_net, mean_gross=mean_gross, W_circ=W_circ, W_lift=W_lift,
                I_MPa_per_Ls=I / 1e9, dp_MPa=Q * I / 1e6, dp_pump_MPa=float(np.interp(1.0, years, dp_t)), W_circ30=float(W_circ_t[-1]), t_sweep_yr=t_sweep, t_10pct_yr=t10, area_km2=2 * Na * wings * tau * L * phi * B / 1e6,
                x_30=30.0 / t_sweep, y_30=0.5 * (W / Na / tau) / math.sqrt(k / rhoc * 30.0 * YEAR), TD=TD.tolist(),
                T_sd_K=(TD_sd * (T_rock - s["T_inj"])).tolist(), GWh=float(np.trapezoid(np.concatenate([[net[0]], net]), np.concatenate([[0.0], years])) * 8.766))


def strip(d):
    return {k: v for k, v in d.items() if k not in ("T_res", "gross", "net", "TD", "T_sd_K")}


def screen():
    out = dict(years=YRS.tolist(), systems={})
    Ns = [5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000]
    for sid, s in SYS.items():
        o = dict(name=s["name"], base=design(s))
        # A: density x tortuosity at base flow and volume
        o["density"] = [dict(tau=tau, rows=[strip(design(s, N=N, tau=tau)) for N in Ns]) for tau in (1.0, 2.0, 5.0)]
        o["density_curves"] = [design(s, N=N) for N in (5, 20, 100, 1000)] + [design(s, N=s["N0"], tau=5.0)]
        # B: flow at fixed volume, dense fractures
        ms = [10, 20, 30, 40, 60, 80, 100, 150, 200, 300]
        o["flow"] = [dict(label=lab, rows=[strip(design(s, N=N, m=float(m))) for m in ms]) for lab, N in (("base fractures", s["N0"]), ("dense (1 m ... 10 m spacing)", int(s["W"] / 5)))]
        # C: volume scale-up at dense fracturing (S = 10 m), flows up to the quoted capacity
        vol = []
        for vf in (1, 2, 4, 8, 16, 32, 64):
            Wv, Lv = geometry(sid, vf)
            rows = [strip(design(s, N=int(Wv / 10.0), m=float(m), W=Wv, L=Lv)) for m in ms]
            vol.append(dict(vf=vf, W=Wv, L=Lv, V_km3=Wv * Lv * s["B"] * s.get("wings", 1) / 1e9, rows=rows))
        o["volume"] = vol
        # D: channelling
        Nd = int(s["W"] / 10.0)
        o["channel"] = [strip(design(s, N=Nd, f=f, phi=phi)) for f in (1.0, 0.5, 0.2, 0.1, 0.05) for phi in (1.0, 0.5, 0.25)]
        # E: what the quoted capacity needs: smallest (flow, volume) with net year-30 >= 90 % of year 1 and mean net >= target
        need = []
        for target in sorted(set(s["quoted"])):
            best = None
            for v in vol:
                for r in v["rows"]:
                    if r["mean_net"] >= target and r["net30"] >= 0.9 * r["net1"]:
                        if best is None or r["m"] * v["V_km3"] < best["m"] * best["V_km3"]:
                            best = dict(r, vf=v["vf"])
            need.append(dict(target=target, design=best))
        o["need"] = need
        out["systems"][sid] = o
    save("ST2_screen", out)
    for sid, o in out["systems"].items():
        b = o["base"]; print(sid, "base: T30", round(b["T30"], 1), "gross1", round(b["gross1"], 2), "gross30", round(b["gross30"], 2), "net mean", round(b["mean_net"], 2), "t_sweep", round(b["t_sweep_yr"], 1))
        for d in o["density"]:
            print("  tau", d["tau"], [(r["N"], round(r["T30"]), round(r["mean_net"], 2)) for r in d["rows"]])
        for fl in o["flow"]:
            print("  flow", fl["label"], [(r["m"], round(r["T30"]), round(r["mean_net"], 2)) for r in fl["rows"]])
        for v in o["volume"]:
            print("  vol x", v["vf"], round(v["V_km3"], 3), [(r["m"], round(r["T30"]), round(r["mean_net"], 1)) for r in v["rows"]])
        print("  need", [(n["target"], n["design"] and (n["design"]["m"], n["design"]["V_km3"], n["design"]["N"], round(n["design"]["mean_net"], 1))) for n in o["need"]])
    return out


def need_fine():
    """What each quoted capacity needs, on a fine grid: the smallest flow whose 30-year mean net output reaches the target when the volume is
    just large enough for the net output of year 30 to stay within 10 % of year 1 (fractures 10 m apart). Added to results/ST2_screen.json."""
    out = load("ST2_screen")
    for sid, s in SYS.items():
        res = []
        for target in sorted(set(s["quoted"])):
            found = None
            for m in np.arange(20.0, 400.1, 2.5):
                def held(vf):
                    Wv, Lv = geometry(sid, vf)
                    d = design(s, N=max(int(Wv / 10.0), 2), m=float(m), W=Wv, L=Lv)
                    return d, d["net30"] >= 0.9 * d["net1"] and d["net1"] > 0
                d_hi, ok = held(256.0)
                if not ok or d_hi["mean_net"] < target:
                    continue
                lo, hi = 0.25, 256.0
                for _ in range(40):
                    mid = math.sqrt(lo * hi)
                    if held(mid)[1]:
                        hi = mid
                    else:
                        lo = mid
                d, _ = held(hi)
                if d["mean_net"] >= target:
                    found = dict(strip(d), vf=hi, gross1=d["gross1"])
                    break
            res.append(dict(target=target, design=found))
            print(sid, "target", target, found and dict(m=found["m"], V_km3=round(found["V_km3"], 3), vf=round(found["vf"], 2), mean_net=round(found["mean_net"], 2),
                                                           x30=round(found["x_30"], 2), y30=round(found["y_30"], 2), W_circ=round(found["W_circ"], 2), gross1=round(found["gross1"], 2)), flush=True)
        out["systems"][str(sid)]["need_fine"] = res
    save("ST2_screen", out)


# ------------------------------------------------------------------------------------------------ full-chain confirmation
CONFIRM = {
    10: [dict(tag="example", N=20, m=30.0, vf=1), dict(tag="dense, example flow", N=100, m=30.0, vf=1), dict(tag="dense, 80 kg/s", N=100, m=80.0, vf=1),
         dict(tag="dense, 100 kg/s, 2x volume", N=141, m=100.0, vf=2), dict(tag="dense, 150 kg/s, 4x volume", N=200, m=150.0, vf=4)],
    14: [dict(tag="example", N=88, m=57.96, vf=1), dict(tag="8x volume, 58 kg/s", N=390, m=57.96, vf=8), dict(tag="16x volume, 80 kg/s", N=390, m=80.0, vf=16),
         dict(tag="16x volume, 100 kg/s", N=390, m=100.0, vf=16)],
    141: [dict(tag="example", N=150, m=100.0, vf=1), dict(tag="2x volume, 60 kg/s", N=304, m=60.0, vf=2), dict(tag="4x volume, 100 kg/s", N=304, m=100.0, vf=4),
          dict(tag="8x volume, 200 kg/s", N=304, m=200.0, vf=8)],
    19: [dict(tag="example", N=20, m=40.0, vf=1), dict(tag="dense, example flow", N=80, m=40.0, vf=1), dict(tag="dense, 80 kg/s, 2x volume", N=113, m=80.0, vf=2),
         dict(tag="dense, 150 kg/s, 4x volume", N=160, m=150.0, vf=4), dict(tag="dense, 200 kg/s, 8x volume", N=226, m=200.0, vf=8)],
}


def geometry(sid, vf):
    """Block width and flow length for vf times the example's volume. Vertical doublets grow in both directions; the horizontal-well
    systems lengthen their laterals first (up to `stretch` times) and widen the well spacing after that."""
    s = SYS[sid]
    Wv = s["W"] * math.sqrt(vf) if "stretch" not in s else s["W"] * min(vf, s["stretch"])
    return Wv, s["L"] * vf / (Wv / s["W"])


def confirm_one(args):
    sid, c = args
    s = SYS[sid]
    t0 = time.time()
    W, L = geometry(sid, c["vf"])
    m = c["m"]; Q = m / s["rho_w"]
    rho, cp = (2650.0, 850.0) if sid != 19 else (2700.0, 950.0)
    rock = [dict(k=s["k"], rho=rho, cp=cp)]
    geo = dict(gradient_C_per_km=(s["T_rock"] - 15.0) / (s["depth"] / 1000.0))
    wings = s.get("wings", 1)
    fac = lambda T_bot: ParallelFractures(Q=Q, N=c["N"] * wings, L_flow=L, B=s["B"], S=W / c["N"], k_r=s["k"], rhoc_r=s["rhoc"], T_rock_C=s["T_rock"],
                                           T_inj_C=lambda t: T_bot(t) - 273.15, rhoc_w=cbar(s) * s["rho_w"])
    un = c.get("uneven")
    if un:                                         # unequal flows between the fractures (one arrangement)
        cor = max(int(round(un.get("corridor_m", 0.0) / (W / c["N"]))), 1)
        fac = lambda T_bot: UnevenFractures(Q=Q, N=c["N"] * wings, L_flow=L, B=s["B"], S=W / c["N"], k_r=s["k"], rhoc_r=s["rhoc"], T_rock_C=s["T_rock"],
                                             T_inj_C=lambda t: T_bot(t) - 273.15, rhoc_w=cbar(s) * s["rho_w"], lorenz=un["lorenz"], corridor=cor,
                                             seed=un.get("seed", 0), n_pattern=min(c["N"], 160))
    big = m > (110.0 if sid == 141 else 65.0)      # Cape Phase I produced about 100 kg/s through 7 in casing
    d_in = (8.535 if big else 6.094) * IN
    layers = [dict(od_in=9.625 if big else 7.0, material="steel"), dict(od_in=12.25 if big else 8.5, material="cement")]
    kw = dict(layers=layers)
    if sid == 10:
        kw.update(p_inj=1.0e6, p_wh=2.0e6, p_res=49.0e6, PI=20e-6 * max(1.0, m / 30.0))
    elif sid == 14:
        survey = dict(md_m=[0.0, s["depth"], s["depth"] + W], tvd_m=[0.0, s["depth"], s["depth"]])
        kw.update(p_inj=3.0e6, p_wh=2.0e6, p_res=22.0e6, PI=15e-6 * max(1.0, m / 57.96), prod_survey=survey, inj_survey=survey)
    elif sid == 141:
        survey = dict(md_m=[0.0, s["depth"], s["depth"] + W / 3.0], tvd_m=[0.0, s["depth"], s["depth"]])
        kw.update(p_inj=15.0e6, p_wh=2.07e6, p_res=30.0e6, PI=50e-6 * max(1.0, m / 100.0), prod_survey=survey, inj_survey=survey, mdot_inj=m / 2.0)
    else:
        kw.update(p_inj=5.0e6, p_wh=5.0e6, p_res=34.0e6, PI=10e-6 * max(1.0, m / 40.0))
    try:
        ch = doublet(f"confirm {sid}", fac, s["depth"], d_in, m, s["T_inj"], rock, 15.0, geo, 30.0, **kw)
        ts = ch["producer"][1].ts
        out = dict(sid=sid, tag=c["tag"], N=c["N"], m=m, vf=c["vf"], W=W, L=L, years=[], T_wh=[], T_res=[], x=[], gross=[], pump_kW=[],
                   unconverged=int(ch["producer"][1].meta.get("unconverged_steps", 0)), wall_s=None)
        from study_cl import power
        for y in (1, 5, 10, 20, 30):
            t = y * YEAR
            h = float(np.interp(t, ts["t"], ts["h_out"])); p = float(np.interp(t, ts["t"], ts["p_out"]))
            out["years"].append(y); out["T_wh"].append(float(np.interp(t, ts["t"], ts["T_out"])) - 273.15); out["T_res"].append(ch["T_res"](t) - 273.15)
            out["x"].append(float(np.interp(t, ts["t"], ts["x_out"]))); out["gross"].append(max(power(m, h, p, 15.0), 0.0)); out["pump_kW"].append(float(np.interp(t, ts["t"], ts["pump_kW"])))
        scr = design(s, N=c["N"], m=m, W=W, L=L, uneven=dict(un, seeds=(un.get("seed", 0),)) if un else None)
        out["screen_gross"] = [float(np.interp(y, YRS, scr["gross"])) for y in out["years"]]
        out["screen_T_res"] = [float(np.interp(y, YRS, scr["T_res"])) for y in out["years"]]
        out["screen_net"] = [float(np.interp(y, YRS, scr["net"])) for y in out["years"]]
        from dev_types import loop_pump
        lp = loop_pump(ch, m, scr["I_MPa_per_Ls"] * 1e9, s["T_inj"], out["years"], rho_res=s["rho_w"])
        out["loop_pump_MW"] = [lp[y]["kW"] / 1e3 for y in out["years"]]          # closed-loop balance on the simulated columns
        out["screen_pump_MW"] = [float(np.interp(y, YRS, np.array(scr["gross"]) - np.array(scr["net"]))) if scr["net"][0] > 0 else None for y in out["years"]]
        out["thermosiphon_MPa"] = [lp[y]["thermosiphon_MPa"] for y in out["years"]]
        out["W_circ"] = scr["W_circ"]; out["wall_s"] = time.time() - t0
        return out
    except Exception as e:
        import traceback
        return dict(sid=sid, tag=c["tag"], failed=f"{type(e).__name__}: {e}", tb=traceback.format_exc()[-600:])


def confirm():
    from concurrent.futures import ProcessPoolExecutor, as_completed
    jobs = [(sid, c) for sid in SYS for c in CONFIRM[sid]]
    out = []
    with ProcessPoolExecutor(2) as ex:
        for f in as_completed([ex.submit(confirm_one, j) for j in jobs]):
            r = f.result(); out.append(r)
            print(r["sid"], r["tag"], r.get("failed") or ("T_wh", [round(v) for v in r["T_wh"]], "gross", [round(v, 2) for v in r["gross"]], "screen", [round(v, 2) for v in r["screen_gross"]],
                                                          "pump", [round(v, 2) for v in r["loop_pump_MW"]], "screen W_circ yr1", round(r["W_circ"], 2), round(r["wall_s"])), flush=True)
            save("ST2_confirm", out)


def confirm_refresh():
    """Recompute the screening columns of results/ST2_confirm.json (after a change to the screening model) without rerunning the chains."""
    out = load("ST2_confirm")
    for r in out:
        if "failed" in r:
            continue
        sc = SYS[r["sid"]]; W, L = geometry(r["sid"], r["vf"])
        scr = design(sc, N=r["N"], m=r["m"], W=W, L=L)
        r["screen_gross"] = [float(np.interp(y, YRS, scr["gross"])) for y in r["years"]]
        r["screen_T_res"] = [float(np.interp(y, YRS, scr["T_res"])) for y in r["years"]]
        r["screen_net"] = [float(np.interp(y, YRS, scr["net"])) for y in r["years"]]
        r["screen_pump_MW"] = [float(np.interp(y, YRS, np.array(scr["gross"]) - np.array(scr["net"]))) for y in r["years"]]
        r["W_circ"] = scr["W_circ"]; r["V_km3"] = scr["V_km3"]; r["t_sweep_yr"] = scr["t_sweep_yr"]
    save("ST2_confirm", out)
    return out


# ------------------------------------------------------------------------------------------------ cyclic single fracture
def cyclic():
    from aharrsim.reservoir import CyclicFracture
    p = 0.5e6
    out = []
    for T_rock in (150.0, 200.0, 250.0):
        for af in (1, 2, 4, 8, 16):
            cf = CyclicFracture(L_f=300.0, B=100.0 * af, w=0.008 / af, k_r=3.0, rhoc_r=2.4e6, T_rock_C=T_rock, T_inj_C=40.0, n_x=60, n_z=60, z_max=25.0)
            cyc = cf.simulate(V_cycle=200.0, t_inj=4 * 3600, t_soak=12 * 3600, t_prod=8 * 3600, n_cycles=180)
            T = cyc[-1]["T_mean_C"]; E = cyc[-1]["E_J"]
            m_avg = 200.0 * 950.0 / 86400.0
            h = CP.PropsSI("H", "T", T + 273.15 - 5.0, "P", p, "Water")       # 5 K lost in the production leg (type 12 result)
            w = plant.orc_practical(m_avg, h, p, 15.0)["W_net"] / 1e6 if T > 65 else 0.0
            out.append(dict(T_rock=T_rock, area_factor=af, area_m2=2 * 300.0 * 100.0 * af, T_day1=cyc[0]["T_mean_C"], T_day180=T, MWth=E / 86400 / 1e6, MWe=max(w, 0.0),
                            MWe_20_wells=20 * max(w, 0.0), recovery=(T - 40.0) / (T_rock - 40.0)))
            print(out[-1], flush=True)
    save("ST2_cyclic", out)


# ------------------------------------------------------------------------------------------------ part 3
DESIGNS = {   # (label, kwargs for design())
    10: [("example (20 fractures, 30 kg/s)", dict()), ("100 fractures, 80 kg/s", dict(N=100, m=80.0)), ("2x volume, 100 kg/s", dict(N=141, m=100.0, vf=2))],
    14: [("example (88 fractures, 58 kg/s)", dict()), ("8x volume, 58 kg/s", dict(N=390, vf=8)), ("16x volume, 80 kg/s", dict(N=390, m=80.0, vf=16))],
    141: [("example (150 fractures, 100 kg/s)", dict()), ("4x volume, 100 kg/s", dict(N=304, vf=4)), ("8x volume, 200 kg/s", dict(N=304, m=200.0, vf=8))],
    19: [("example (20 fractures, 40 kg/s)", dict()), ("2x volume, 80 kg/s", dict(N=113, m=80.0, vf=2)), ("8x volume, 200 kg/s", dict(N=226, m=200.0, vf=8))],
}
P3 = [   # variable, low, high, label
    ("T_rock", 0.8, 1.2, "Gradient −20 / +20 % (rock temperature)"), ("k", 2.0, 4.0, "Rock conductivity 2.0 / 4.0 W/m/K (at 25 °C)"),
    ("kT", 1, 1, "k falls with temperature (Vosteen)"), ("rhoc", 0.85, 1.15, "Rock heat capacity −15 / +15 %"),
    ("f", 0.2, 1.0, "Active fractures 20 / 100 %"), ("phi", 0.5, 1.0, "Swept volume 50 / 100 %"), ("I0", 2.0, 0.5, "Impedance ×2 / ×0.5"),
    ("T_amb", 45.0, 15.0, "Ambient 45 / 15 °C"),
]


def eval_design(sid, kw, T_rock_f=1.0, k25=None, kT=False, rhoc_f=1.0, f=1.0, phi=1.0, I0_f=1.0, T_amb=15.0):
    s = SYS[sid]
    kw = dict(kw); vf = kw.pop("vf", 1)
    W, L = geometry(sid, vf)
    T_rock = s["T_s"] + (s["T_rock"] - s["T_s"]) * T_rock_f
    k = s["k"] if k25 is None else k25
    if kT:
        k = float(k_of_T(k, 0.5 * (T_rock + s["T_inj"])))     # conductivity at the mean of rock and injection temperature
    return design(s, W=W, L=L, T_rock=T_rock, k=k, rhoc=s["rhoc"] * rhoc_f, f=f, phi=phi, I0=s["I0"] * I0_f, T_amb=T_amb, **kw)


def part3():
    out = dict(years=YRS.tolist(), systems={})
    rng = np.random.default_rng(7)
    n = 3000
    for sid in SYS:
        o = []
        for label, kw in DESIGNS[sid]:
            base = eval_design(sid, kw)
            torn = []
            for var, lo, hi, lab in P3:
                rows = []
                for v in (lo, hi):
                    a = dict(T_rock_f=v) if var == "T_rock" else dict(k25=v) if var == "k" else dict(kT=True) if var == "kT" else dict(rhoc_f=v) if var == "rhoc" else \
                        dict(f=v) if var == "f" else dict(phi=v) if var == "phi" else dict(I0_f=v) if var == "I0" else dict(T_amb=v)
                    rows.append(strip(eval_design(sid, kw, **a)))
                torn.append(dict(var=var, label=lab, low=rows[0], high=rows[1]))
            mc = {}
            for T_amb in (15.0, 35.0):
                U = rng.random((n, 6))
                res = []
                for i in range(n):
                    d = eval_design(sid, kw, T_rock_f=float(tri(U[i, 0], 0.8, 1.0, 1.2)), k25=float(tri(U[i, 1], 2.0, SYS[sid]["k"], 4.0)), kT=True,
                                    rhoc_f=float(tri(U[i, 2], 0.85, 1.0, 1.15)), f=float(tri(U[i, 3], 0.1, 0.5, 1.0)), phi=float(tri(U[i, 4], 0.4, 0.8, 1.0)),
                                    I0_f=float(np.exp(tri(U[i, 5], np.log(0.5), 0.0, np.log(2.0)))), T_amb=T_amb)
                    res.append((d["mean_net"], d["net1"], d["net30"], d["T30"], d["mean_gross"]))
                A = np.array(res)
                mc[str(int(T_amb))] = dict(mean_net=np.percentile(A[:, 0], [10, 50, 90]).tolist(), net1=np.percentile(A[:, 1], [10, 50, 90]).tolist(),
                                           net30=np.percentile(A[:, 2], [10, 50, 90]).tolist(), T30=np.percentile(A[:, 3], [10, 50, 90]).tolist(),
                                           mean_gross=np.percentile(A[:, 4], [10, 50, 90]).tolist(), samples_mean_net=np.sort(A[:, 0])[:: n // 100].tolist(),
                                           p_zero30=float(np.mean(A[:, 2] <= 0.01)))
            # combined favourable / unfavourable ends of the closed-loop drivers (gradient, conductivity) with the full net-power history
            fav = eval_design(sid, kw, T_rock_f=1.2, k25=4.0, rhoc_f=1.15)
            unf = eval_design(sid, kw, T_rock_f=0.8, k25=2.0, kT=True, rhoc_f=0.85)
            hot = eval_design(sid, kw, T_amb=45.0)
            o.append(dict(label=label, kw=kw, base=base, tornado=torn, mc=mc, favourable=fav, unfavourable=unf, hot_ambient=hot))
            print(sid, label, "base mean net", round(base["mean_net"], 2), "P10/50/90 @15C", [round(v, 2) for v in mc["15"]["mean_net"]], "@35C", [round(v, 2) for v in mc["35"]["mean_net"]],
                  "fav", round(fav["mean_net"], 2), "unf", round(unf["mean_net"], 2), flush=True)
        out["systems"][sid] = dict(name=SYS[sid]["name"], quoted=SYS[sid]["quoted"], designs=o)
    save("ST3_open_uncertainty", out)


if __name__ == "__main__":
    dict(screen=screen, confirm=confirm, refresh=confirm_refresh, cyclic=cyclic, part3=part3, need=need_fine)[sys.argv[1] if len(sys.argv) > 1 else "screen"]()
