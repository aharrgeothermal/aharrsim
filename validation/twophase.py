"""Validation of the flashing (steam/water) up-flow hydraulics against measured
flowing pressure/temperature surveys of geothermal production wells.

Data: Ambastha, A.K. & Gudmundsson, J.S. (1986) SGP-TR-100, Appendix C (and
the workshop paper SGP-TR-93 pp. 183-188, Tables 1 and 3).  Pressures in the
appendix are gauge (bar-g); 1.013 bar is added.  The comparison marches the
same pressure-gradient model used inside aharrsim's riser (homogeneous or
drift-flux, Churchill friction, acceleration, IAPWS-95 properties) downward
from the shallowest measured point, with the flowing enthalpy fixed by the
adiabatic steady-flow energy equation (h + g z + v^2/2 = const).
"""
from __future__ import annotations

import math

import numpy as np

from common import save
from aharrsim.props import WaterIAPWS
from aharrsim.correlations import two_phase_gradient, friction_churchill, G_ACC

FT = 0.3048

WELLS = [
    dict(name="Los Azufres Az-18", W=26.7, h=1607e3,
         strings=[(0.0, 959.0, 0.7296 * FT), (959.0, 1324.0, 0.5153 * FT)],
         z=[0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1050, 1100, 1150, 1200, 1250, 1300, 1324],
         P=[30, 31.2, 32.2, 33.1, 34.2, 35.2, 36.2, 37.4, 38.4, 39.5, 40.8, 41.9, 42.9, 43.8, 45, 46.5, 50.1, 52.1],
         T=[238, 242.5, 244, 246, 248, 249, 251, 253, 254, 256, 257.5, 258.8, 260.2, 261.5, 263, 267, 271, 272],
         z_max=1250.0, orkiszewski=dict(mean_bar=-1.1, sd_bar=1.2, mean_pct=-2.65, sd_pct=2.2),
         note="feed zone near 1250 m (static liquid gradient below); compared down to 1250 m"),
    dict(name="HGP-A, 13.9 kg/s", W=13.9, h=1966e3,
         strings=[(0.0, 680.0, 0.802 * FT), (680.0, 1966.0, 0.5833 * FT)],
         z=[0, 160, 305, 401, 465, 642, 770, 914, 1067, 1219, 1316, 1380, 1524, 1684, 1765, 1829, 1925],
         P=[3.2, 4, 4.6, 4.9, 5.2, 5.9, 7.2, 8.6, 9.8, 11.5, 12.2, 12.6, 14.1, 15.5, 16.1, 16.7, 16.7],
         T=[146, 152.7, 158.4, 160, 162.4, 166.2, 174.3, 181.3, 186.5, 192.4, 196.2, 196.2, 200, 202.7, 205.4, 207, 207],
         z_max=1316.0, orkiszewski=dict(mean_bar=0.6, sd_bar=0.4, mean_pct=6.1, sd_pct=2.7),
         note="two-phase throughout; second feed zone at ~1311 m (4300 ft), so only 0-1316 m carries the full flow"),
    dict(name="HGP-A, 13.6 kg/s", W=13.6, h=1966e3,
         strings=[(0.0, 680.0, 0.802 * FT), (680.0, 1966.0, 0.5833 * FT)],
         z=[0, 127, 305, 401, 465, 642, 770, 914, 1067, 1219, 1316, 1380, 1524, 1684, 1765, 1829, 1925],
         P=[6, 6.3, 6.6, 6.9, 7, 7.9, 8.6, 10.6, 12.5, 14.1, 15.2, 15.4, 16.3, 17.8, 18.4, 19.4, 20.1],
         T=[166, 170.3, 175, 176, 178.4, 181.9, 181.9, 191.1, 197.8, 204, 208.1, 208.6, 212.4, 212.2, 213.8, 216.8, 219.2],
         z_max=1316.0, orkiszewski=None, note="compared above the second feed zone (0-1316 m)"),
    dict(name="HGP-A, 12.8 kg/s", W=12.8, h=1966e3,
         strings=[(0.0, 680.0, 0.802 * FT), (680.0, 1966.0, 0.5833 * FT)],
         z=[0, 160, 305, 401, 465, 642, 770, 914, 1067, 1219, 1316, 1380, 1524, 1684, 1765, 1829, 1925],
         P=[15, 16.1, 17.2, 17.5, 17.7, 18.8, 20.1, 21.6, 23.3, 24.7, 25.9, 26.5, 28.8, 31, 31.6, 32.9, 34],
         T=[201, 203.2, 206.8, 208.1, 208.6, 213, 216.8, 219.2, 221.8, 223, 225.7, 227.8, 232.4, 237.8, 238.5, 240.5, 243.2],
         z_max=1316.0, orkiszewski=None, note="compared above the second feed zone (0-1316 m)"),
    dict(name="Cerro Prieto M-90", W=45.0, h=1343e3,
         strings=[(0.0, 1299.0, 0.5808 * FT)],
         z=[25, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1225, 1250, 1299],
         P=[40.9, 42.5, 45.3, 48.1, 51, 54, 57.4, 60.5, 64.3, 68.3, 72.2, 76.8, 82, 83.6, 85, 88.5],
         T=[249, 250.5, 254, 257.5, 261, 264.5, 269, 272, 275.5, 279, 283, 287, 290.5, 291.5, 292, 292],
         z_max=1299.0, orkiszewski=dict(mean_bar=-0.3, sd_bar=0.8, mean_pct=-0.6, sd_pct=1.1),
         note="measured T 4-10 C below Tsat(P): dissolved CO2 likely (not modelled)"),
    dict(name="Utah State 14-2", W=40.9, h=1648e3,
         strings=[(0.0, 913.0, 0.7433 * FT)],
         z=[10, 41, 85, 169, 233, 296, 359, 388, 420, 445, 469, 495, 516, 539, 573, 616, 645, 675, 701, 725,
            747, 772, 786, 807, 828, 846, 863, 884, 901, 913],
         P=[27, 27.8, 29.4, 31.2, 32.8, 34.5, 36.1, 37, 37.9, 38.7, 39.5, 40.3, 41.2, 42, 43.2, 45, 46.2, 47.4,
            48.7, 49.9, 51.2, 52.5, 53.3, 54.5, 55.8, 57, 58.3, 59.6, 60.8, 61.6],
         T=None, z_max=913.0, orkiszewski=dict(mean_bar=-6.7, sd_bar=4.6, mean_pct=-13.6, sd_pct=6.9),
         note="pressure only; Orkiszewski also failed here (authors found no cause)"),
]


def march(well, model, rough=0.0006 * FT, W_=None, dz=5.0):
    """March pressure downward from the shallowest survey point."""
    f = WaterIAPWS()
    W = well["W"] if W_ is None else W_
    z0 = well["z"][0]
    p = (well["P"][0] + 1.01325) * 1e5
    h_ref = well["h"]
    zs = [z0]
    ps = [p]
    Ts = []
    xs = []

    def D_at(z):
        for a, b, d in well["strings"]:
            if a <= z < b + 1e-9:
                return d
        return well["strings"][-1][2]

    def terms(p, h, D):
        G = W / (math.pi * D * D / 4.0)
        tp, rho, mu, x, rl, rg, ml, mg, sig, T = f.point(p, h)
        if tp:
            rg_, fric, rm, al = two_phase_gradient(model, G, D, x, rl, rg, ml, mg, sig, rough, 1.0)
            return float(rg_), float(fric), float(rm), T, x, G
        Re = G * D / mu
        fr = friction_churchill(Re, rough / D)
        return rho, float(fr * G * G / (2 * D * rho)), rho, T, 0.0, G

    z = z0
    h = h_ref
    st = terms(p, h, D_at(z + 0.1))
    Ts.append(st[3]); xs.append(st[4])
    while z < well["z_max"] - 1e-6:
        step = min(dz, well["z_max"] - z)
        D = D_at(z + 0.5 * step)
        st0 = terms(p, h, D)
        # adiabatic: h + g*elev + v^2/2 = const  -> h increases downward
        h1 = h + G_ACC * step
        p1 = p + st0[0] * G_ACC * step
        for _ in range(20):
            st1 = terms(p1, h1, D)
            G = st0[5]
            dp = 0.5 * (st0[0] + st1[0]) * G_ACC * step + 0.5 * (st0[1] + st1[1]) * step \
                + G * G * (1.0 / st0[2] - 1.0 / st1[2])
            # kinetic-energy correction to the enthalpy
            h1n = h + G_ACC * step + 0.5 * ((G / st0[2]) ** 2 - (G / st1[2]) ** 2)
            p1n = p + dp
            if abs(p1n - p1) < 1.0 and abs(h1n - h1) < 0.1:
                p1, h1 = p1n, h1n
                break
            p1, h1 = 0.5 * (p1 + p1n), h1n
        z += step
        p, h = p1, h1
        st = terms(p, h, D)
        zs.append(z); ps.append(p); Ts.append(st[3]); xs.append(st[4])
    return np.array(zs), np.array(ps), np.array(Ts), np.array(xs)


def main():
    out = []
    for w in WELLS:
        row = dict(name=w["name"], W=w["W"], h_kJ=w["h"] / 1e3, note=w["note"], orkiszewski=w["orkiszewski"],
                   z_meas=w["z"], P_meas=[p + 1.01325 for p in w["P"]], T_meas=w["T"])
        zm = np.array(w["z"], float)
        sel = (zm >= w["z"][0]) & (zm <= w["z_max"] + 1e-6)
        Pm = np.array(w["P"], float)[sel] + 1.01325
        for model in ("hem", "drift_flux"):
            zs, ps, Ts, xs = march(w, model)
            Pc = np.interp(zm[sel], zs, ps) / 1e5
            err = Pc - Pm
            pct = 100 * (Pc - Pm) / Pm
            row[model] = dict(z=zs, P=ps / 1e5, T=Ts - 273.15, x=xs,
                              mean_bar=float(err[1:].mean()), sd_bar=float(err[1:].std(ddof=1)),
                              mean_pct=float(pct[1:].mean()), sd_pct=float(pct[1:].std(ddof=1)),
                              bottom_err_bar=float(err[-1]))
            if w["T"] is not None:
                Tm = np.array(w["T"], float)[sel]
                Tc = np.interp(zm[sel], zs, Ts - 273.15)
                row[model]["T_err_mean"] = float(np.mean(Tc - Tm))
            print(f"{w['name']:22s} {model:10s} mean {err[1:].mean():+6.2f} bar (SD {err[1:].std(ddof=1):4.2f}),"
                  f" {pct[1:].mean():+6.2f} % (SD {pct[1:].std(ddof=1):4.2f}); bottom {err[-1]:+6.2f} bar;"
                  f" x_top {xs[0]:.3f}")
        if w["orkiszewski"]:
            o = w["orkiszewski"]
            print(f"{'':22s} Orkiszewski (paper): {o['mean_bar']:+.1f} bar (SD {o['sd_bar']}), {o['mean_pct']:+.2f} % (SD {o['sd_pct']})")
        out.append(row)
    save("T1_twophase", out)


if __name__ == "__main__":
    main()
