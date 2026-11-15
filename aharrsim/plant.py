"""Surface plants: what the produced fluid is worth as electricity or heat.

orc()        subcritical organic Rankine cycle on a pressurised liquid (or two-phase) geothermal stream:
             best of six working fluids and of the evaporating temperature, optional recuperator,
             5 K pinch, turbine 0.85, feed pump 0.75, generator 0.96, condenser fans/pumps as a
             fraction of the heat rejected.  Air cooling: T_cond = T_amb + 15 K; wet cooling: T_wb + 10 K.
flash()      single or double flash of a two-phase / hot-water stream: separator pressure(s) optimised,
             steam turbine with Baumann wet-steam correction, condenser at T_cond, gas extraction and
             cooling parasitics.  dry_steam() is the same turbine on a steam feed.
heat_pump()  ground-source heat pump: COP from a second-law efficiency on the Carnot COP.
exergy()     thermal exergy of a stream relative to the dead state.
geophires_eta_u() / geophires_reinjection_T(): the GEOPHIRES v3 utilisation-efficiency correlations
             (Aspen-derived, NREL) used here only as a published reference for benchmarking.

Energy accounting: W_net = W_gross - feed/condensate pumps - cooling-system power - gas extraction;
circulation / downhole pumping of the wells is deducted by the caller (systems.py).
"""
from __future__ import annotations

import math

import numpy as np
import CoolProp.CoolProp as CP

ORC_FLUIDS = ["IsoButane", "n-Butane", "R245fa", "Isopentane", "n-Pentane", "Toluene"]
PINCH, ETA_T, ETA_P, ETA_G, EPS_REC = 5.0, 0.85, 0.75, 0.96, 0.85
FAN_AIR, FAN_WET = 0.02, 0.012          # kWe per kWth rejected (air-cooled condenser fans / cooling tower + pumps)
_W = "Water"


def _h(T_C, p):
    return CP.PropsSI("H", "T", T_C + 273.15, "P", p, _W)


def _T(h, p):
    return CP.PropsSI("T", "H", h, "P", p, _W) - 273.15


def exergy(m, h, p, T0_C=20.0):
    """Thermal exergy [W] of m kg/s at (h, p) relative to liquid water at T0 and the same pressure."""
    T0 = T0_C + 273.15
    s = CP.PropsSI("S", "H", h, "P", p, _W)
    h0 = CP.PropsSI("H", "T", T0, "P", p, _W)
    s0 = CP.PropsSI("S", "T", T0, "P", p, _W)
    return m * ((h - h0) - T0 * (s - s0))


# ----------------------------------------------------------------------------
def orc(m_b, h_in, p_b, T_cond_C=30.0, T_b_min_C=40.0, cooling="air", fluids=ORC_FLUIDS, dT_ev=2.0):
    """Best net ORC power for a brine stream m_b [kg/s] at enthalpy h_in [J/kg] and pressure p_b [Pa].
    The brine may be two-phase (it condenses in the evaporator).  Returns a dict (W_net, W_gross, Q_in, Q_rej,
    fluid, T_ev, recuperator, T_b_out, eta_th)."""
    fan = FAN_AIR if cooling == "air" else FAN_WET
    T_hot = _T(h_in, p_b)
    best = dict(W_net=0.0, W_gross=0.0, Q_in=0.0, Q_rej=0.0, fluid=None, T_ev=None, recuperator=False, T_b_out=T_hot, eta_th=0.0)
    for fl in fluids:
        Tc = CP.PropsSI("Tcrit", fl) - 273.15
        p1 = CP.PropsSI("P", "T", T_cond_C + 273.15, "Q", 0, fl)
        h1 = CP.PropsSI("H", "T", T_cond_C + 273.15, "Q", 0, fl)
        s1 = CP.PropsSI("S", "T", T_cond_C + 273.15, "Q", 0, fl)
        T_top = min(Tc - 3.0, T_hot - PINCH - 0.5)
        if T_top <= T_cond_C + 10.0:
            continue
        for T_ev in np.arange(T_cond_C + 10.0, T_top, dT_ev):
            p_ev = CP.PropsSI("P", "T", T_ev + 273.15, "Q", 1, fl)
            h2 = h1 + (CP.PropsSI("H", "P", p_ev, "S", s1, fl) - h1) / ETA_P
            T2 = CP.PropsSI("T", "P", p_ev, "H", h2, fl) - 273.15
            hf = CP.PropsSI("H", "T", T_ev + 273.15, "Q", 0, fl)
            h3 = CP.PropsSI("H", "T", T_ev + 273.15, "Q", 1, fl)
            s3 = CP.PropsSI("S", "T", T_ev + 273.15, "Q", 1, fl)
            h4 = h3 - ETA_T * (h3 - CP.PropsSI("H", "P", p1, "S", s3, fl))
            T4 = CP.PropsSI("T", "P", p1, "H", h4, fl) - 273.15
            for recup in (False, True):
                q_rec = 0.0
                if recup:
                    if T4 <= T2 + PINCH:
                        continue
                    h_hot_min = max(CP.PropsSI("H", "P", p1, "T", T2 + PINCH + 273.15, fl) if T2 + PINCH > T_cond_C + 0.5
                                    else h4, CP.PropsSI("H", "P", p1, "Q", 1, fl))
                    q_hot = max(0.0, h4 - h_hot_min)
                    q_cold = max(0.0, min(CP.PropsSI("H", "P", p_ev, "T", T4 - PINCH + 273.15, fl), hf) - h2)
                    q_rec = EPS_REC * min(q_hot, q_cold)
                    if q_rec <= 0:
                        continue
                h2r = h2 + q_rec
                T2r = CP.PropsSI("T", "P", p_ev, "H", h2r, fl) - 273.15
                # working-fluid flow limited by the pinch at the bubble point and by the brine outlet limit
                m_pinch = m_b * (h_in - _h(T_ev + PINCH, p_b)) / (h3 - hf)
                T_out_min = max(T_b_min_C, T2r + PINCH)
                m_out = m_b * (h_in - _h(T_out_min, p_b)) / (h3 - h2r)
                m_wf = min(m_pinch, m_out)
                if m_wf <= 0:
                    continue
                Q_in = m_wf * (h3 - h2r)
                Q_rej = m_wf * (h4 - q_rec - h1)
                W_gross = m_wf * (h3 - h4) * ETA_G
                W_net = W_gross - m_wf * (h2 - h1) - fan * Q_rej
                if W_net > best["W_net"]:
                    best = dict(W_net=W_net, W_gross=W_gross, Q_in=Q_in, Q_rej=Q_rej, fluid=fl, T_ev=float(T_ev),
                                recuperator=recup, T_b_out=_T(h_in - Q_in / m_b, p_b), eta_th=W_net / Q_in)
    return best


# ----------------------------------------------------------------------------
def _steam_turbine(h_in, p_in, p_out, x_in=1.0):
    """Expansion of steam from (h_in, p_in) to p_out: dry efficiency ETA_T with the Baumann rule for wetness."""
    s_in = CP.PropsSI("S", "H", h_in, "P", p_in, _W)
    h_s = CP.PropsSI("H", "P", p_out, "S", s_in, _W)
    x_s = CP.PropsSI("Q", "P", p_out, "H", h_s, _W)
    x_out_s = x_s if 0 <= x_s <= 1 else 1.0
    eta = ETA_T * 0.5 * (x_in + x_out_s) if x_out_s < 1.0 else ETA_T        # Baumann: 1 % per 1 % mean wetness
    for _ in range(3):
        h_out = h_in - eta * (h_in - h_s)
        x_out = CP.PropsSI("Q", "P", p_out, "H", h_out, _W)
        x_out = x_out if 0 <= x_out <= 1 else 1.0
        eta = ETA_T * 0.5 * (x_in + x_out) if x_out < 1.0 else ETA_T
    return h_out, x_out


def flash(m, h_in, p_in, T_cond_C=40.0, double=False, p_sep=None, p_sep2=None, cooling="wet", ncg_frac=0.01, p_min=1.2e5):
    """Single (or double) flash plant on m kg/s arriving at enthalpy h_in and pressure p_in.  Separator pressure(s)
    optimised for net power unless given (never above the arrival pressure, nor above saturation of a liquid
    feed).  Condenser at T_cond_C (direct-contact, wet cooling by default).  Parasitics: cooling system
    (FAN_* x heat rejected, pumps included), gas extraction ncg_frac of gross.
    Returns a dict (W_net, W_gross, m_steam, x_sep, p_sep, p_sep2, h_brine_out, T_brine_out, Q_rej)."""
    p_c = CP.PropsSI("P", "T", T_cond_C + 273.15, "Q", 0, _W)
    fan = FAN_AIR if cooling == "air" else FAN_WET
    T_in = _T(h_in, p_in)
    x_in = CP.PropsSI("Q", "H", h_in, "P", p_in, _W)
    p_hi = p_in if 0.0 <= x_in <= 1.0 else min(p_in, CP.PropsSI("P", "T", T_in + 273.15, "Q", 0, _W))

    def evaluate(p1, p2=None):
        if p1 <= p_c * 1.5:
            return None
        hf1, hg1 = CP.PropsSI("H", "P", p1, "Q", 0, _W), CP.PropsSI("H", "P", p1, "Q", 1, _W)
        x1 = (h_in - hf1) / (hg1 - hf1)
        if x1 <= 0.0:
            return None
        x1 = min(x1, 1.0)
        m_s1 = m * x1
        m_b1 = m - m_s1
        h_t1, x_e1 = _steam_turbine(hg1, p1, p_c)
        W = m_s1 * (hg1 - h_t1)
        Q_rej = m_s1 * (h_t1 - CP.PropsSI("H", "P", p_c, "Q", 0, _W))
        out = dict(p_sep=p1, x_sep=x1, m_steam=m_s1, h_brine_out=hf1, p_sep2=None)
        if p2 is not None and m_b1 > 0:
            if not (p_c * 1.5 < p2 < p1):
                return None
            hf2, hg2 = CP.PropsSI("H", "P", p2, "Q", 0, _W), CP.PropsSI("H", "P", p2, "Q", 1, _W)
            x2 = (hf1 - hf2) / (hg2 - hf2)
            if x2 <= 0:
                return None
            m_s2 = m_b1 * x2
            h_t2, x_e2 = _steam_turbine(hg2, p2, p_c)
            W += m_s2 * (hg2 - h_t2)
            Q_rej += m_s2 * (h_t2 - CP.PropsSI("H", "P", p_c, "Q", 0, _W))
            out.update(p_sep2=p2, m_steam=m_s1 + m_s2, h_brine_out=hf2, x_sep2=x2)
        W_gross = W * ETA_G
        W_net = W_gross * (1.0 - ncg_frac) - fan * Q_rej
        out.update(W_gross=W_gross, W_net=W_net, Q_rej=Q_rej, T_brine_out=_T(out["h_brine_out"], max(out.get("p_sep2") or p1, 1e5)))
        return out

    best = None
    grid1 = [p_sep] if p_sep else np.geomspace(max(p_min, 1.5 * p_c), max(p_hi * 0.98, 2 * p_min), 40)
    for p1 in grid1:
        if not double:
            r = evaluate(p1)
            if r and (best is None or r["W_net"] > best["W_net"]):
                best = r
        else:
            grid2 = [p_sep2] if p_sep2 else np.geomspace(max(p_min, 1.5 * p_c), p1 * 0.95, 25)
            for p2 in grid2:
                r = evaluate(p1, p2)
                if r and (best is None or r["W_net"] > best["W_net"]):
                    best = r
    if best is None:
        best = dict(W_net=0.0, W_gross=0.0, m_steam=0.0, x_sep=0.0, p_sep=None, p_sep2=None, h_brine_out=h_in,
                    T_brine_out=_T(h_in, 1e6), Q_rej=0.0)
    return best


def dry_steam(m, h_in, p_in, T_cond_C=40.0, cooling="wet", ncg_frac=0.015):
    """Steam feed (saturated or superheated) expanded directly to the condenser."""
    p_c = CP.PropsSI("P", "T", T_cond_C + 273.15, "Q", 0, _W)
    fan = FAN_AIR if cooling == "air" else FAN_WET
    h_out, x_out = _steam_turbine(h_in, p_in, p_c, x_in=min(1.0, max(0.0, CP.PropsSI("Q", "H", h_in, "P", p_in, _W))) if
                                  0 <= CP.PropsSI("Q", "H", h_in, "P", p_in, _W) <= 1 else 1.0)
    W_gross = m * (h_in - h_out) * ETA_G
    Q_rej = m * (h_out - CP.PropsSI("H", "P", p_c, "Q", 0, _W))
    return dict(W_gross=W_gross, W_net=W_gross * (1 - ncg_frac) - fan * Q_rej, Q_rej=Q_rej, x_exhaust=x_out)


# ----------------------------------------------------------------------------
def heat_pump(Q_source_W, T_source_C, T_sink_C, eta_II=0.45):
    """Heat pump lifting Q_source from T_source (evaporator) to T_sink (condenser): COP = eta_II * T_sink/(T_sink - T_source)."""
    Th, Tc = T_sink_C + 273.15, T_source_C + 273.15
    COP = eta_II * Th / max(Th - Tc, 1.0)
    W = Q_source_W / (COP - 1.0)
    return dict(COP=COP, W_e=W, Q_heat=Q_source_W + W)


def direct_use(m, h_in, p, T_return_C):
    return m * (h_in - _h(T_return_C, p))


# ----------------------------------------------------------------------------
# GEOPHIRES v3 reference correlations (SurfacePlant*.py, NREL; ambient-temperature band 15-25 C used here)
def geophires_eta_u(kind, T_in_C, T_amb_C=15.0):
    T = np.asarray(T_in_C, float)
    if T_amb_C < 15.0:
        frac = (T_amb_C - 5.0) / 10.0
        coef = {"orc": ((0.0, 2.746e-3, -8.3806e-2), (0.0, 2.713e-3, -9.1841e-2)),
                "sorc": ((-1.55e-5, 7.604e-3, -3.78e-1), (-1.499e-5, 7.4268e-3, -3.7915e-1)),
                "flash1": ((-4.27318e-7, 8.65629e-4, 1.78931e-1), (-5.85412e-7, 9.68352e-4, 1.58056e-1)),
                "flash2": ((-1.2e-6, 1.22731e-3, 2.26956e-1), (-1.42165e-6, 1.3705e-3, 1.99847e-1))}[kind]
    else:
        frac = (T_amb_C - 15.0) / 10.0
        coef = {"orc": ((0.0, 2.713e-3, -9.1841e-2), (0.0, 2.676e-3, -1.012e-1)),
                "sorc": ((-1.499e-5, 7.4268e-3, -3.7915e-1), (-1.55e-5, 7.55136e-3, -4.041e-1)),
                "flash1": ((-5.85412e-7, 9.68352e-4, 1.58056e-1), (-7.78996e-7, 1.0923e-3, 1.33708e-1)),
                "flash2": ((-1.42165e-6, 1.3705e-3, 1.99847e-1), (-1.66771e-6, 1.53079e-3, 1.69439e-1))}[kind]
    (a2, a1, a0), (b2, b1, b0) = coef
    return (1 - frac) * (a2 * T ** 2 + a1 * T + a0) + frac * (b2 * T ** 2 + b1 * T + b0)


def geophires_availability(T0_C, T1_C, T2_C):
    """GEOPHIRES 'availability_water' [J/kg] of liquid water between T1 and T2 with dead state T0."""
    A, B, Cc = 4.041650, -1.204e-2, 1.605e-5
    T0, T1, T2 = T0_C + 273.15, T1_C + 273.15, T2_C + 273.15
    av = ((A - B * T0) * (T1 - T2) + (B - Cc * T0) / 2.0 * (T1 ** 2 - T2 ** 2) + Cc / 3.0 * (T1 ** 3 - T2 ** 3)
          - A * T0 * math.log(T1 / T2)) * 2.2046 / 947.83
    return av * 1e6


def geophires_reinjection_T(kind, T_in_C, T_amb_C=15.0):
    T = np.asarray(T_in_C, float)
    frac = (T_amb_C - 15.0) / 10.0 if T_amb_C >= 15 else (T_amb_C - 5.0) / 10.0
    if T_amb_C < 15.0:
        coef = {"orc": ((0.0, 0.0894, 55.6), (0.0, 0.0894, 62.6)), "sorc": ((0.0, 0.02, 49.26), (0.0, 0.02, 56.26)),
                "flash1": ((-1.11519e-3, 7.79126e-1, -10.2242), (-1.10232e-3, 7.83893e-1, -5.17039)),
                "flash2": ((-7.70928e-4, 5.02466e-1, 5.22091), (-7.69455e-4, 5.09406e-1, 11.6859))}[kind]
    else:
        coef = {"orc": ((0.0, 0.0894, 62.6), (0.0, 0.0894, 69.6)), "sorc": ((0.0, 0.02, 56.26), (0.0, 0.02, 63.26)),
                "flash1": ((-1.10232e-3, 7.83893e-1, -5.17039), (-1.08914e-3, 7.88562e-1, -1.89707e-1)),
                "flash2": ((-7.69455e-4, 5.09406e-1, 11.6859), (-7.67751e-4, 5.16356e-1, 18.0798))}[kind]
    (a2, a1, a0), (b2, b1, b0) = coef
    return (1 - frac) * (a2 * T ** 2 + a1 * T + a0) + frac * (b2 * T ** 2 + b1 * T + b0)


def orc_practical(m, h_in, p_in, T_amb_C=15.0, kind="orc"):
    """Practical binary-plant estimate: GEOPHIRES utilisation efficiency x availability of the stream from its
    temperature down to ambient (liquid-water availability as in GEOPHIRES), brine leaving at the GEOPHIRES
    reinjection temperature.  kind 'orc' (subcritical) or 'sorc' (supercritical).  Returns dict (W_net, T_b_out, eta_u)."""
    T_in = _T(h_in, p_in)
    if T_in <= 70.0:
        return dict(W_net=0.0, T_b_out=T_in, eta_u=0.0)
    T_c = min(T_in, 250.0)              # the correlations were fitted up to about 250 C; hold the efficiency there
    eta = float(geophires_eta_u(kind, T_c, T_amb_C))
    av = geophires_availability(T_amb_C, T_in, T_amb_C)
    T_out = float(geophires_reinjection_T(kind, T_c, T_amb_C))
    return dict(W_net=max(0.0, eta * av * m), T_b_out=min(T_out, T_in), eta_u=max(eta, 0.0))
