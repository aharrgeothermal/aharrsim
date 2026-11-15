"""Case definition: geometry, materials, rock, geotherm, schedule, numerics.

A case is normally read from a YAML/JSON dictionary (see examples/).  Units in
the input are given by key suffixes: _m, _in (inches), _C, _MPa, _kg_s,
_days, _s.  Internally everything is SI with temperatures in kelvin.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, List, Optional
import math

import numpy as np

from .props import WaterIAPWS, ConstantFluid

INCH = 0.0254
DAY = 86400.0
YEAR = 365.25 * DAY


def _radius(d: dict, key: str) -> Optional[float]:
    """Return a radius [m] from a dict holding <key>_m (diameter in m), <key>_in or <key>_mm."""
    if f"{key}_m" in d:
        return 0.5 * float(d[f"{key}_m"])
    if f"{key}_in" in d:
        return 0.5 * float(d[f"{key}_in"]) * INCH
    if f"{key}_mm" in d:
        return 0.5 * float(d[f"{key}_mm"]) / 1000.0
    return None


@dataclass
class Material:
    name: str
    k: float          # W/m/K
    rho: float        # kg/m3
    cp: float         # J/kg/K

    @property
    def rhoc(self):
        return self.rho * self.cp


@dataclass
class CasingSection:
    """Wellbore construction over an MD interval, listed from the annulus outwards.

    r_inner : radius of the annulus outer wall (innermost casing ID, or the
              open-hole radius when `layers` is empty) [m]
    layers  : list of (r_out [m], Material); contiguous, the last r_out is the
              borehole radius
    """
    top: float
    bottom: float
    r_inner: float
    layers: List[tuple]

    @property
    def r_hole(self):
        return self.layers[-1][0] if self.layers else self.r_inner


@dataclass
class TubingSection:
    """Centre (return) pipe over an MD interval.

    r_in   : inner radius [m]
    layers : list of (r_out [m], k [W/m/K]) from inside out; for vacuum-insulated
             tubing give the insulation as an equivalent conductivity layer
    R_wall : optional override of the total wall resistance per unit length [K m/W]
    """
    top: float
    bottom: float
    r_in: float
    layers: List[tuple]
    R_wall: Optional[float] = None

    @property
    def r_out(self):
        return self.layers[-1][0] if self.layers else self.r_in

    def wall_resistance(self):
        if self.R_wall is not None:
            return self.R_wall
        R, r = 0.0, self.r_in
        for r_out, k in self.layers:
            R += math.log(r_out / r) / (2.0 * math.pi * k)
            r = r_out
        return R


@dataclass
class RockLayer:
    top: float        # TVD m
    bottom: float     # TVD m
    k: float
    rho: float
    cp: float


@dataclass
class Period:
    duration: float          # s
    mdot: float              # kg/s
    T_in: float              # K
    label: str = ""


@dataclass
class Survey:
    md: np.ndarray
    tvd: np.ndarray

    def tvd_at(self, s):
        s = np.asarray(s, float)
        out = np.interp(s, self.md, self.tvd)
        # extrapolate beyond the last station along the last segment
        if self.md.size >= 2:
            slope = (self.tvd[-1] - self.tvd[-2]) / (self.md[-1] - self.md[-2])
            beyond = s > self.md[-1]
            out = np.where(beyond, self.tvd[-1] + slope * (s - self.md[-1]), out)
        return out

    def md_at_tvd(self, z):
        """First MD at which TVD reaches z (used to place layer boundaries)."""
        if z <= self.tvd[0]:
            return float(self.md[0])
        for i in range(1, self.md.size):
            if self.tvd[i] >= z > self.tvd[i - 1] or (self.tvd[i] >= z and self.tvd[i] == self.tvd[i - 1]):
                if self.tvd[i] == self.tvd[i - 1]:
                    return float(self.md[i - 1])
                f = (z - self.tvd[i - 1]) / (self.tvd[i] - self.tvd[i - 1])
                return float(self.md[i - 1] + f * (self.md[i] - self.md[i - 1]))
        return None


@dataclass
class Case:
    name: str
    fluid: object
    survey: Survey
    md_total: float
    casing: List[CasingSection]
    tubing: List[TubingSection]
    rock_layers: List[RockLayer]
    T_surface: float                       # K
    geotherm: Callable                     # T0(tvd) -> K
    schedule: List[Period]
    pressure_bc: dict                      # {'outlet': Pa} or {'inlet': Pa}
    min_inlet_p: float = 0.2e6             # Pa; below this the loop is throttled at the wellhead
    flow_direction: str = "annulus_down"   # 'annulus_down' (return up the tubing) or 'tubing_down' (return up the annulus)
    mode: str = "coaxial"                  # 'coaxial' | 'production' (one channel, flow up from a feed zone at TD) |
                                           # 'injection' (one channel, flow down into the formation at TD) |
                                           # 'utube' (U-tube borehole heat exchanger: both legs in a grouted borehole)
    far_field: str = "dirichlet"           # outer radial boundary: 'dirichlet' (undisturbed T) or 'adiabatic' (unit cell)
    inflow: Optional[dict] = None          # production/injection: {'T': K or callable(t)->K, 'h': J/kg or callable,
                                           #   'p_res': Pa, 'PI': kg/s/Pa (productivity or injectivity index)}
    pump: Optional[dict] = None            # production: {'md': m, 'eta': -} downhole pump (ESP / line-shaft)
    utube: Optional[dict] = None           # utube: {'r_pipe_in', 'r_pipe_out', 'k_pipe', 'shank', 'k_grout'} (m, W/m/K)
    bottom_element: object = None          # coaxial: reservoir element in series at TD, with .T_out(T_in_K, t_s) -> K
    single_storage: bool = False           # production/injection: keep the fluid heat-storage term (off: steady-flow fluid energy
                                           # equation, transient rock only; the mixture density of a flashing column changes
                                           # too fast between iterations for the linearised storage term)
    roughness: float = 4.6e-5
    turn_loss_K: float = 1.0
    # numerics
    dz: float = 50.0
    dz_min: float = 2.0
    coupling_frac: float = 0.2
    r_far: Optional[float] = None
    radial_growth: float = 1.15
    ext_length: Optional[float] = None
    dt_initial: float = 300.0
    dt_growth: float = 1.3
    dt_max: float = 30 * DAY
    output_times: List[float] = field(default_factory=list)
    two_phase_model: str = "drift_flux"
    max_iter: int = 40
    tol_h: float = 5.0          # J/kg
    tol_T: float = 1e-3         # K
    tol_p: float = 200.0        # Pa
    scheme: str = "box"         # 'box' (2nd order) or 'upwind'
    time_scheme: str = "bdf2"   # 'bdf2' (2nd order, variable step) or 'euler'
    # options (verification switches)
    gravity_work: bool = True
    include_pressure: bool = True
    rock_model: str = "grid"    # 'grid' or 'resistance'
    R_formation: Optional[float] = None
    h_wall_override: Optional[float] = None
    fluid_storage: bool = True
    notes: str = ""

    # ------------------------------------------------------------------
    def rock_at(self, tvd):
        for L in self.rock_layers:
            if L.top <= tvd < L.bottom:
                return L
        return self.rock_layers[-1] if tvd >= self.rock_layers[-1].bottom else self.rock_layers[0]

    @property
    def t_end(self):
        return sum(p.duration for p in self.schedule)


# ----------------------------------------------------------------------------
# Parsing
# ----------------------------------------------------------------------------
def _material(name, mats):
    if name not in mats:
        raise KeyError(f"material '{name}' not defined under 'materials'")
    m = mats[name]
    return Material(name, float(m["k"]), float(m["rho"]), float(m["cp"]))


def _geotherm_from(d: dict, T_s: float, layers: List[RockLayer]):
    """Build T0(tvd) [K] from the geotherm block."""
    if "table" in d:
        z = np.array([r[0] for r in d["table"]], float)
        T = np.array([r[1] for r in d["table"]], float) + 273.15
        slope_end = (T[-1] - T[-2]) / (z[-1] - z[-2])

        def f(tvd):
            tvd = np.asarray(tvd, float)
            out = np.interp(tvd, z, T)
            return np.where(tvd > z[-1], T[-1] + slope_end * (tvd - z[-1]), out)
        return f
    if "segments" in d:          # [[top_m, gradient_C_per_km], ...]
        seg = sorted([(float(a), float(b) / 1000.0) for a, b in d["segments"]])
        tops = np.array([s[0] for s in seg])
        grads = np.array([s[1] for s in seg])
        Ttop = [T_s]
        for i in range(1, len(seg)):
            Ttop.append(Ttop[-1] + grads[i - 1] * (tops[i] - tops[i - 1]))
        Ttop = np.array(Ttop)

        def f(tvd):
            tvd = np.asarray(tvd, float)
            i = np.clip(np.searchsorted(tops, tvd, side="right") - 1, 0, len(tops) - 1)
            return Ttop[i] + grads[i] * (tvd - tops[i])
        return f
    if "gradient_C_per_km" in d:
        g = float(d["gradient_C_per_km"]) / 1000.0
        return lambda tvd: T_s + g * np.asarray(tvd, float)
    if "heat_flow_mW_m2" in d:   # conductive geotherm through the layered rock
        q = float(d["heat_flow_mW_m2"]) / 1000.0
        zs = np.linspace(0.0, 30000.0, 30001)
        kz = np.array([next((L.k for L in layers if L.top <= z < L.bottom), layers[-1].k) for z in zs])
        Tz = T_s + np.concatenate([[0.0], np.cumsum(q / kz[:-1] * np.diff(zs))])
        return lambda tvd: np.interp(np.asarray(tvd, float), zs, Tz)
    raise ValueError("geotherm: give one of table / segments / gradient_C_per_km / heat_flow_mW_m2")


def _inflow_from(d):
    """Feed-zone description for production / injection wells (temperatures in C, pressures in MPa in the file)."""
    if d is None:
        return None
    out = {}
    if "T_C" in d:
        out["T"] = (lambda t, f=d["T_C"]: float(f(t)) + 273.15) if callable(d["T_C"]) else float(d["T_C"]) + 273.15
    if "T_K" in d:
        out["T"] = d["T_K"]
    if "h_J_kg" in d:
        out["h"] = d["h_J_kg"]
    if "p_res_MPa" in d:
        out["p_res"] = float(d["p_res_MPa"]) * 1e6
    if "PI_kg_s_per_MPa" in d:
        out["PI"] = float(d["PI_kg_s_per_MPa"]) / 1e6
    return out


def _pump_from(d):
    if d is None:
        return None
    return dict(md=float(d["md_m"]), eta=float(d.get("efficiency", 0.7)))


def _utube_from(d):
    if d is None:
        return None
    return dict(r_pipe_in=0.5 * _len(d, "pipe_id"), r_pipe_out=0.5 * _len(d, "pipe_od"), k_pipe=float(d.get("k_pipe", 0.4)),
                shank=_len(d, "shank_spacing"), k_grout=float(d.get("k_grout", 1.5)))


def _len(d, key):
    """Length from '<key>_m' or '<key>_mm' or '<key>_in'."""
    if key + "_m" in d:
        return float(d[key + "_m"])
    if key + "_mm" in d:
        return float(d[key + "_mm"]) / 1000.0
    if key + "_in" in d:
        return float(d[key + "_in"]) * 0.0254
    raise KeyError(key)


def case_from_dict(d: dict) -> Case:
    name = d.get("name", "case")
    # fluid ------------------------------------------------------------------
    fl = d.get("fluid", "iapws")
    if isinstance(fl, dict) and "constant" in fl:
        c = fl["constant"]
        fluid = ConstantFluid(rho=c.get("rho", 1000.0), cp=c.get("cp", 4180.0),
                              mu=c.get("mu", 1e-3), k=c.get("k", 0.6))
    else:
        fluid = WaterIAPWS(backend=(fl.get("backend", "BICUBIC&HEOS") if isinstance(fl, dict) else "BICUBIC&HEOS"))
    # well -------------------------------------------------------------------
    w = d["well"]
    if "survey" in w:
        md = np.array(w["survey"]["md_m"], float)
        tvd = np.array(w["survey"]["tvd_m"], float)
    else:
        depth = float(w["depth_m"])
        lat = float(w.get("lateral_m", 0.0))
        md = np.array([0.0, depth] + ([depth + lat] if lat > 0 else []))
        tvd = np.array([0.0, depth] + ([depth] if lat > 0 else []))
    survey = Survey(md, tvd)
    md_total = float(md[-1])
    mats = d.get("materials", {})
    mats.setdefault("steel", {"k": 45.0, "rho": 7850.0, "cp": 490.0})
    mats.setdefault("cement", {"k": 1.0, "rho": 1900.0, "cp": 1000.0})
    casing = []
    for cs in w["casing"]:
        r_in = _radius(cs, "id")
        layers = []
        for L in cs.get("layers", []):
            r_out = _radius(L, "od")
            if "material" in L:
                mat = _material(L["material"], mats)
            else:
                mat = Material(L.get("name", "custom"), float(L["k"]), float(L["rho"]), float(L["cp"]))
            layers.append((r_out, mat))
        top = float(cs.get("top_m", 0.0))
        bot = float(cs.get("bottom_m", md_total))
        r = r_in
        for r_out, _ in layers:
            if r_out <= r:
                raise ValueError(f"casing layers must increase in radius ({top}-{bot} m)")
            r = r_out
        casing.append(CasingSection(top, bot, r_in, layers))
    tubing = []
    for ts in w.get("tubing", []):
        r_in = _radius(ts, "id")
        layers = [(_radius(L, "od"), float(L["k"])) for L in ts.get("layers", [])]
        tubing.append(TubingSection(float(ts.get("top_m", 0.0)), float(ts.get("bottom_m", md_total)),
                                    r_in, layers, ts.get("R_wall_K_m_W")))
    # rock -------------------------------------------------------------------
    r = d["rock"]
    rock_layers = [RockLayer(float(L.get("top_m", 0.0)), float(L.get("bottom_m", 1e6)),
                             float(L["k"]), float(L["rho"]), float(L["cp"])) for L in r["layers"]]
    T_s = float(r["surface_T_C"]) + 273.15
    geotherm = _geotherm_from(r["geotherm"], T_s, rock_layers)
    # operation --------------------------------------------------------------
    op = d["operation"]
    sched = []
    for p in op["schedule"]:
        dur = float(p["duration_days"]) * DAY if "duration_days" in p else float(p["duration_s"])
        sched.append(Period(dur, float(p["mdot_kg_s"]), float(p["T_in_C"]) + 273.15, p.get("label", "")))
    pbc = op.get("pressure_bc", {"outlet_MPa": 1.0})
    if "outlet_MPa" in pbc:
        pressure_bc = {"outlet": float(pbc["outlet_MPa"]) * 1e6}
    elif "bottom_MPa" in pbc:
        pressure_bc = {"bottom": float(pbc["bottom_MPa"]) * 1e6}      # production well marched up from TD
    else:
        pressure_bc = {"inlet": float(pbc["inlet_MPa"]) * 1e6}
    min_inlet_p = float(pbc.get("min_inlet_MPa", 0.2)) * 1e6
    # numerics ---------------------------------------------------------------
    n = d.get("numerics", {})
    o = d.get("options", {})
    case = Case(
        name=name, fluid=fluid, survey=survey, md_total=md_total, casing=casing, tubing=tubing,
        rock_layers=rock_layers, T_surface=T_s, geotherm=geotherm, schedule=sched,
        pressure_bc=pressure_bc, min_inlet_p=min_inlet_p, roughness=float(w.get("roughness_m", 4.6e-5)),
        turn_loss_K=float(w.get("turn_loss_K", 1.0)),
        dz=float(n.get("dz_m", 50.0)),
        dz_min=float(n.get("dz_min_m", 2.0)),
        coupling_frac=float(n.get("coupling_frac", 0.2)),
        r_far=(None if n.get("r_far_m", "auto") in (None, "auto") else float(n["r_far_m"])),
        radial_growth=float(n.get("radial_growth", 1.15)),
        ext_length=(None if n.get("ext_length_m") is None else float(n["ext_length_m"])),
        dt_initial=float(n.get("dt_initial_s", 300.0)),
        dt_growth=float(n.get("dt_growth", 1.3)),
        dt_max=float(n.get("dt_max_days", 30.0)) * DAY,
        output_times=[float(t) * DAY for t in n.get("output_times_days", [])],
        two_phase_model=n.get("two_phase_model", "drift_flux"),
        max_iter=int(n.get("max_iter", 40)),
        tol_h=float(n.get("tol_h", 5.0)),
        scheme=n.get("scheme", "box"),
        time_scheme=n.get("time_scheme", "bdf2"),
        gravity_work=bool(o.get("gravity_work", True)),
        include_pressure=bool(o.get("include_pressure", True)),
        flow_direction=str(w.get("flow_direction", op.get("flow_direction", "annulus_down"))),
        mode=str(w.get("mode", "coaxial")), far_field=str(o.get("far_field", "dirichlet")),
        inflow=_inflow_from(w.get("inflow")), pump=_pump_from(w.get("pump")), utube=_utube_from(w.get("utube")),
        bottom_element=w.get("bottom_element"), single_storage=bool(o.get("single_storage", False)),
        rock_model=o.get("rock_model", "grid"),
        R_formation=o.get("R_formation_K_m_W"),
        h_wall_override=o.get("h_wall_override_W_m2K"),
        fluid_storage=bool(o.get("fluid_storage", True)),
        notes=d.get("notes", ""),
    )
    _validate(case)
    return case


def _validate(c: Case):
    if c.flow_direction not in ("annulus_down", "tubing_down"):
        raise ValueError("flow_direction must be 'annulus_down' (return up the tubing) or 'tubing_down' "
                         "(return up the annulus)")
    if c.mode not in ("coaxial", "production", "injection", "utube"):
        raise ValueError("mode must be coaxial, production, injection or utube")
    if c.far_field not in ("dirichlet", "adiabatic"):
        raise ValueError("far_field must be dirichlet or adiabatic")
    if c.mode == "utube" and c.utube is None:
        raise ValueError("mode utube needs a 'utube' block (pipe_id, pipe_od, shank_spacing, k_pipe, k_grout)")
    if c.mode == "production" and (c.inflow is None or ("T" not in c.inflow and "h" not in c.inflow)):
        raise ValueError("mode production needs an 'inflow' block with the feed temperature (T_C) or enthalpy")
    if c.pump is not None and c.mode != "production":
        raise ValueError("a downhole pump is only used in production mode")
    if c.pump is not None and (c.inflow is None or "p_res" not in c.inflow or "PI" not in c.inflow):
        raise ValueError("a downhole pump needs inflow.p_res_MPa and inflow.PI_kg_s_per_MPa (to split the pressure march)")
    if c.bottom_element is not None and (c.mode != "coaxial" or c.flow_direction != "annulus_down"):
        raise ValueError("bottom_element is only used with the coaxial, annulus-down configuration")
    s = np.linspace(0, c.md_total, 2001)[:-1] + 0.5 * c.md_total / 2000
    for x in s:
        cs = [q for q in c.casing if q.top <= x < q.bottom]
        if not cs:
            raise ValueError(f"no casing/open-hole section covers MD {x:.1f} m")
        if c.mode == "coaxial":
            ts = [q for q in c.tubing if q.top <= x < q.bottom]
            if not ts:
                raise ValueError(f"no tubing section covers MD {x:.1f} m (tubing must run to TD)")
            if ts[0].r_out >= cs[0].r_inner:
                raise ValueError(f"tubing OD >= casing ID at MD {x:.1f} m")
    if c.mode == "utube":
        u = c.utube
        if 2 * u["r_pipe_out"] >= u["shank"] or u["shank"] / 2 + u["r_pipe_out"] >= min(q.r_inner for q in c.casing):
            raise ValueError("U-tube pipes overlap each other or the borehole wall")
