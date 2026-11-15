"""Prose and page skeleton for the uncertainty / open-system study report. Numbers come from the results via build(C, f)."""
from __future__ import annotations


def pct(a, b):
    return 100.0 * (a - b) / b if b else float("nan")


def build(C, f):
    t = {}
    sgn = lambda v, d=0: ("+" if v >= 0 else "−") + f(abs(v), d)

    # ------------------------------------------------------------------ Part 1
    def oat(c, var, e=False):
        r = C.get(("oate" if e else "oat") + str(c), {}).get(var)
        return (r["low"], r["high"]) if r else (float("nan"), float("nan"))

    def ph(c, var, v):
        return C.get(f"phys{c}_{var}_{v}", dict(T=float("nan"), MWth=float("nan"), MWe=float("nan")))

    def mc(c, k):
        return C.get(f"mc{c}", {}).get(k, [float("nan")] * 3)

    b7H, b8H, b20H = C.get("b7_H", 0), C.get("b8_H", 0), C.get("b20_H", 0)
    b7T, b8T, b20T = C.get("b7_T", 0), C.get("b8_T", 0), C.get("b20_T", 0)
    b8E, b20E = C.get("b8_E", 0), C.get("b20_E", 0)

    g, k_, cp, ins = oat(7, "grad"), oat(7, "k"), oat(7, "cp"), oat(7, "ins_k")
    kT7, nu3, nu10, nu100, lost7 = ph(7, "kT", True), ph(7, "nu", 3.0), ph(7, "nu", 10.0), ph(7, "nu", 100.0), ph(7, "ins_k", 0.3)
    m7 = mc(7, "MWth")
    t["p1_text_7"] = (
        f"Base case: 3 km, 80 °C/km, 10 kg/s in at 40 °C, giving {f(b7T, 0)} °C and {f(b7H, 2)} MWth at year 30 and no electricity. "
        f"Gradient (60–100 °C/km) moves the heat from {f(g[0], 2)} to {f(g[1], 2)} MWth and rock conductivity (2.0–3.5 W/m/K) from {f(k_[0], 2)} to {f(k_[1], 2)} MWth; "
        f"the two are the whole of the subsurface uncertainty. Heat capacity (±15 %) moves it by {f(abs(cp[1] - cp[0]), 2)} MWth and tubing insulation, anywhere "
        f"between 0.01 and 0.10 W/m/K, by {f(abs(ins[1] - ins[0]), 2)} MWth; even with the insulation lost altogether the well gives {f(lost7['MWth'], 2)} MWth, because at 10 kg/s "
        f"the water spends too little time in a 3 km string to lose much. Letting conductivity fall with temperature takes {f(-pct(kT7['MWth'], b7H), 0)} % off. "
        f"With all four sampled together the heat is {f(m7[0], 2)} / {f(m7[1], 2)} / {f(m7[2], 2)} MWth (P10 / P50 / P90): the base case sits above the median because it ignores the "
        f"fall of conductivity with temperature. None of {C.get('mc7_n', 0)} runs produces useful electricity. The one thing that changes the answer is the physics outside the casing: "
        f"if the lower half of the well sits in a convecting formation, the heat rises to {f(nu3['MWth'], 1)} MWth at a Nusselt number of 3, {f(nu10['MWth'], 1)} MWth at 10 "
        f"and {f(nu100['MWth'], 1)} MWth ({f(nu100['T'], 0)} °C) when the casing is held at reservoir temperature by condensing steam.")

    g, k_, cp = oat(8, "grad"), oat(8, "k"), oat(8, "cp")
    ge, ke = oat(8, "grad", True), oat(8, "k", True)
    sp, nl = oat(8, "spacing"), oat(8, "n_lat")
    kT8, n3, n10 = ph(8, "kT", True), ph(8, "nu", 3.0), ph(8, "nu", 10.0)
    m8, m8e = mc(8, "MWth"), mc(8, "MWe")
    t["p1_text_8"] = (
        f"Base case: two 4.5 km verticals joined by twelve 5 km laterals at 80 m spacing (69 km of hole), 55 °C/km, 100 kg/s, giving {f(b8T, 0)} °C, {f(b8H, 0)} MWth "
        f"and {f(b8E, 2)} MWe gross at year 30. Gradient (45–65 °C/km) and conductivity (2.0–4.0 W/m/K) each move the heat by about a fifth either way "
        f"({f(g[0], 0)}–{f(g[1], 0)} and {f(k_[0], 0)}–{f(k_[1], 0)} MWth) and the electricity by far more ({f(ge[0], 2)}–{f(ge[1], 2)} and {f(ke[0], 2)}–{f(ke[1], 2)} MWe), "
        f"because at 95–130 °C the binary plant's conversion efficiency climbs steeply with temperature. Heat capacity matters a little more here than in a single well "
        f"({f(cp[0], 1)}–{f(cp[1], 1)} MWth) because the laterals draw down a bounded block between them. The design levers are of the same size as the uncertainties: "
        f"halving the lateral spacing to 40 m costs {f(b8H - sp[0], 0)} MWth, and going from 8 to 16 laterals spans {f(nl[0], 0)}–{f(nl[1], 0)} MWth. "
        f"Temperature-dependent conductivity takes the base case to {f(kT8['MWth'], 0)} MWth and {f(kT8['MWe'], 2)} MWe. "
        f"Monte Carlo: {f(m8[0], 0)} / {f(m8[1], 0)} / {f(m8[2], 0)} MWth and {f(m8e[0], 2)} / {f(m8e[1], 2)} / {f(m8e[2], 2)} MWe gross (P10 / P50 / P90), "
        f"against the 8.2 MWe quoted for the commercial loop. A formation that convects around the laterals would close the gap: {f(n3['MWth'], 0)} MWth and {f(n3['MWe'], 1)} MWe at Nu 3, "
        f"{f(n10['MWth'], 0)} MWth and {f(n10['MWe'], 1)} MWe at Nu 10.")

    g, k_, cp, ins = oat(20, "grad"), oat(20, "k"), oat(20, "cp"), oat(20, "ins_k")
    ge, ke, ie = oat(20, "grad", True), oat(20, "k", True), oat(20, "ins_k", True)
    kT20, lost = ph(20, "kT", True), ph(20, "ins_k", 0.3)
    m20, m20e = mc(20, "MWth"), mc(20, "MWe")
    t["p1_text_20"] = (
        f"Base case: 20 km, 500 °C at total depth, 0.30 m bore, 20 kg/s, giving {f(b20T, 0)} °C, {f(b20H, 1)} MWth and {f(b20E, 2)} MWe gross at year 30. "
        f"Gradient (400–600 °C at depth) and conductivity (2.0–4.0 W/m/K) have almost identical effect: {f(g[0], 1)}–{f(g[1], 1)} MWth and {f(ge[0], 2)}–{f(ge[1], 2)} MWe for one, "
        f"{f(k_[0], 1)}–{f(k_[1], 1)} MWth and {f(ke[0], 2)}–{f(ke[1], 2)} MWe for the other. Unlike the 3 km well, insulation is a first-order variable: the water spends hours "
        f"rising past 20 km of colder descending water, and a tubing conductivity of 0.10 instead of 0.03 W/m/K cuts the electricity to {f(ie[0], 2)} MWe; with the vacuum lost "
        f"(0.3 W/m/K) the well returns {f(lost['T'], 0)} °C and {f(lost['MWe'], 2)} MWe. The largest single correction is physical, not statistical: at 300–500 °C granite conducts 35–45 % less "
        f"well than at room temperature, and with that included the base case falls to {f(kT20['T'], 0)} °C, {f(kT20['MWth'], 1)} MWth and {f(kT20['MWe'], 2)} MWe "
        f"({f(-pct(kT20['MWe'], b20E), 0)} % less electricity). Monte Carlo with that effect in every run: {f(m20[0], 1)} / {f(m20[1], 1)} / {f(m20[2], 1)} MWth and "
        f"{f(m20e[0], 2)} / {f(m20e[1], 2)} / {f(m20e[2], 2)} MWe (P10 / P50 / P90).")

    # ------------------------------------------------------------------ Part 2
    def s(sid):
        return C.get(f"s{sid}")

    def dens(o, tau, N):
        for d in o["density"]:
            if d["tau"] == tau:
                for r in d["rows"]:
                    if r["N"] == N:
                        return r
        return dict(mean_net=float("nan"), dp_MPa=float("nan"), T30=float("nan"))

    def vol(o, vf, m):
        for v in o["volume"]:
            if v["vf"] == vf:
                for r in v["rows"]:
                    if r["m"] == m:
                        return dict(r, V=v["V_km3"])
        return dict(mean_net=float("nan"), net30=float("nan"), W_circ=float("nan"), V=float("nan"))

    def chan(o, fr, phi):
        for r in o["channel"]:
            if r["f"] == fr and r["phi"] == phi:
                return r
        return dict(mean_net=float("nan"), net30=float("nan"), T30=float("nan"))

    if s(10):
        o = s(10); b = o["base"]
        t["p2_text_10"] = (
            f"Base case: 200 °C rock at 5 km, a 1,000 × 600 × 400 m block ({f(b['V_km3'], 2)} km³) cut by 20 fractures 50 m apart, 30 kg/s. It gives {f(b['gross1'], 2)} MWe gross and "
            f"{f(b['net1'], 2)} MWe net with no thermal decline in 30 years (the block holds {f(b['t_sweep_yr'], 0)} years of heat at that flow). The base case is therefore not short of "
            f"fracture area, and adding fractures does little thermally: 5,000 fractures instead of 20 raise the net only to {f(dens(o, 1.0, 5000)['mean_net'], 2)} MWe, all of it from lower "
            f"flow resistance. Tortuosity does not help at all: a path five times longer through the same rock multiplies the pressure drop "
            f"({f(dens(o, 5.0, 20)['dp_MPa'], 0)} MPa against {f(dens(o, 1.0, 20)['dp_MPa'], 0)} MPa at 20 fractures) and adds nothing to a temperature that is already at the rock's. "
            f"What dense fracturing buys is the ability to carry more flow without cooling: with fractures 10 m apart the same block sustains 80 kg/s and "
            f"{f(vol(o, 1, 80.0)['mean_net'], 1)} MWe net, which is inside the quoted 3–5 MWe; the base fracture set at 80 kg/s cools to nothing useful by year 30. "
            f"Beyond that the limit is hydraulic, not thermal: at 150 kg/s and four times the volume the gross is about {f(vol(o, 4, 150.0)['mean_net'] + vol(o, 4, 150.0)['W_circ'], 0)} MWe "
            f"but {f(vol(o, 4, 150.0)['W_circ'], 1)} MWe goes to the injection pumps, and the net never passes {f(max(r['mean_net'] for v in o['volume'] for r in v['rows']), 1)} MWe at any volume. "
            f"5 MWe net from one doublet at 200 °C needs lower impedance than the 0.3 MPa per L/s assumed, or two producers.")
    if s(14):
        o = s(14); b = o["base"]
        v8, v16 = vol(o, 8, 60.0), vol(o, 16, 100.0)
        t["p2_text_14"] = (
            f"Base case: 190 °C rock, a 1,000 m lateral pair 100 m apart with 250 m of fracture height ({f(b['V_km3'], 3)} km³), 64 fractures, 58 kg/s. Year 1 reproduces the reported "
            f"output: {f(b['gross1'], 1)} MWe gross, {f(b['net1'], 1)} MWe net. It does not last: the stimulated block holds only {f(b['t_sweep_yr'], 1)} years of heat at that flow, "
            f"the outlet has lost 10 % of its temperature rise by year {f(b['t_10pct_yr'], 0)}, and the 30-year mean net is {f(b['mean_net'], 2)} MWe. Fracture density and tortuosity cannot "
            f"fix that, because the fractures are already close enough (16 m) to drain the rock between them; 5,000 fractures give {f(dens(o, 1.0, 5000)['mean_net'], 2)} MWe. "
            f"The shortfall is volume. Eight times the volume ({f(v8['V_km3'] if 'V_km3' in v8 else v8['V'], 2)} km³) sustains 60 kg/s and {f(v8['mean_net'], 1)} MWe net for 30 years; "
            f"sixteen times ({f(v16['V'], 1)} km³, four 1 km lateral pairs' worth of length at 400 m well spacing) sustains 100 kg/s and {f(v16['mean_net'], 1)} MWe, "
            f"which meets the quoted 3.5 MWe as a 30-year figure. In practice a short-lived test and a 30-year plant are different claims: the published 3.5 MWe is a 30-day test result, "
            f"and the model agrees with it as such.")
    if s(19):
        o = s(19); b = o["base"]
        v2, v8a, v8b = vol(o, 2, 80.0), vol(o, 8, 200.0), vol(o, 8, 300.0)
        c1, c2 = chan(o, 1.0, 0.5), chan(o, 0.2, 1.0)
        t["p2_text_19"] = (
            f"Base case: 375 °C rock, an 800 × 300 × 400 m block ({f(b['V_km3'], 2)} km³), 20 fractures 40 m apart, 40 kg/s: {f(b['net1'], 1)} MWe net in year 1 falling to "
            f"{f(b['net30'], 1)} MWe by year 30 (mean {f(b['mean_net'], 1)} MWe) as the outlet cools to {f(b['T30'], 0)} °C. Here density does help, because the base spacing is too wide "
            f"for the rock between fractures to keep up: 100 fractures hold 375 °C for 30 years and give {f(dens(o, 1.0, 100)['mean_net'], 1)} MWe, and tortuosity 2 at 20 fractures gives "
            f"{f(dens(o, 2.0, 20)['mean_net'], 1)} MWe, a real gain from the longer contact path. That is the ceiling for the base flow: each kg/s of 375 °C water is worth about 0.29 MWe, "
            f"so the quoted 50 MWe is a statement about flow. Doubling the volume and the flow (80 kg/s) gives {f(v2['mean_net'], 0)} MWe; eight times the volume "
            f"({f(v8a['V'], 2)} km³) gives {f(v8a['mean_net'], 0)} MWe at 200 kg/s and {f(v8b['mean_net'], 0)} MWe at 300 kg/s, with {f(v8a['W_circ'], 0)}–{f(v8b['W_circ'], 0)} MWe of pumping. "
            f"So 50 MWe per well pair needs roughly 200 kg/s through about 0.8 km³ of uniformly fractured 375 °C rock: three to five times the flow of the best EGS producers to date. "
            f"This case is also the most exposed to poor sweep: if half the block is bypassed the dense design's year-30 output is {f(c1['net30'], 1)} MWe (mean {f(c1['mean_net'], 1)}), "
            f"and if only a fifth of the fractures take flow the mean is {f(c2['mean_net'], 1)} MWe.")

    # ------------------------------------------------------------------ Part 3
    def p3(sid, i):
        o = C.get(f"p3_{sid}")
        return o["designs"][i] if o else None

    def tor(d, var):
        for x in d["tornado"]:
            if x["var"] == var:
                return x["low"]["mean_net"], x["high"]["mean_net"]
        return float("nan"), float("nan")

    def p3text(sid, i, extra):
        d = p3(sid, i)
        if not d:
            return ""
        b = d["base"]["mean_net"]; T = tor(d, "T_rock"); k = tor(d, "k"); rc = tor(d, "rhoc"); ff = tor(d, "f"); ph_ = tor(d, "phi"); I = tor(d, "I0"); ta = tor(d, "T_amb")
        m15, m35 = d["mc"]["15"]["mean_net"], d["mc"]["35"]["mean_net"]
        return (
            f"Design carried forward: {d['label']}, {f(b, 1)} MWe mean net in the base case. A geothermal gradient 20 % either side of the estimate, the range Part 1 used, "
            f"takes it from {f(T[0], 1)} to {f(T[1], 1)} MWe. Conductivity across 2.0–4.0 W/m/K changes it by {f(abs(k[1] - k[0]), 2)} MWe and heat capacity (±15 %) by "
            f"{f(abs(rc[1] - rc[0]), 2)} MWe: with fractures 10 m apart the rock between them empties whatever its conductivity. A 45 °C ambient instead of 15 °C leaves {f(ta[0], 1)} MWe. "
            f"Of the flow-path unknowns, impedance (half to double) spans {f(I[0], 1)}–{f(I[1], 1)} MWe, flow confined to a fifth of the fractures leaves {f(ff[0], 1)} MWe and half the "
            f"volume swept {f(ph_[0], 1)} MWe. Sampling everything together gives {f(m15[0], 1)} / {f(m15[1], 1)} / {f(m15[2], 1)} MWe at a 15 °C site and "
            f"{f(m35[0], 1)} / {f(m35[1], 1)} / {f(m35[2], 1)} MWe at 35 °C (P10 / P50 / P90). " + extra)

    t["p3_text_10"] = p3text(10, 1, "The median is two thirds of the base case because every flow-path distribution is one-sided: channelling and bypass can only subtract. "
                             "The quoted 3–5 MWe sits between the median and the P90 at a cool site, and above the P90 at a hot one. "
                             "This design is also fragile late in life: on the base volume, 80 kg/s through half the fractures breaks through thermally, and in more than half the samples there is no net output left at year 30. "
                             "Doubling the volume cuts that to one in ten.")
    t["p3_text_14"] = p3text(14, 1, "With eight times the published volume the quoted 3.5 MWe is at about the P90 for a cool site; at the published volume the 30-year mean is below 0.2 MWe in nine cases in ten.")
    t["p3_text_19"] = p3text(19, 2, "The spread is the widest of the three in absolute terms, 28 to 88 MWe, and almost all of it is temperature: at 300 °C the same flow gives half the output of 375 °C, "
                             "at 450 °C more than twice. The quoted 50 MWe is near the median of this design, which is the design built to meet it.")

    # ------------------------------------------------------------------ tiles and summary
    d10, d14, d19 = p3(10, 1), p3(14, 1), p3(19, 2)
    t["tile_cl"] = f"{f(m7[0], 2)}–{f(m7[2], 2)}"
    t["tile_cl20"] = f"{f(m20e[0], 2)}–{f(m20e[2], 2)}"
    t["tile_cl8"] = f"{f(m8e[0], 1)}–{f(m8e[2], 1)}"
    if s(10):
        t["tile_dens"] = f"{sgn(pct(dens(s(10), 1.0, 5000)['mean_net'], s(10)['base']['mean_net']))} %"
        t["dens14"] = f(dens(s(14), 1.0, 5000)["mean_net"], 2); t["base14"] = f(s(14)["base"]["mean_net"], 2)
        t["dens19"] = f(dens(s(19), 1.0, 100)["mean_net"], 1); t["base19"] = f(s(19)["base"]["mean_net"], 1)
        t["dens10"] = f(dens(s(10), 1.0, 5000)["mean_net"], 2); t["base10"] = f(s(10)["base"]["mean_net"], 2)
        t["v19"] = f(vol(s(19), 8, 200.0)["mean_net"], 0); t["v14"] = f(vol(s(14), 16, 100.0)["mean_net"], 1); t["v10"] = f(vol(s(10), 1, 80.0)["mean_net"], 1)
        t["sweep14"] = f(s(14)["base"]["t_sweep_yr"], 0)
    if d19:
        t["p3_T19"] = f"{f(tor(d19, 'T_rock')[0], 0)}–{f(tor(d19, 'T_rock')[1], 0)}"
        t["p3_T10"] = f"{f(tor(d10, 'T_rock')[0], 1)}–{f(tor(d10, 'T_rock')[1], 1)}"
        t["p3_T14"] = f"{f(tor(d14, 'T_rock')[0], 1)}–{f(tor(d14, 'T_rock')[1], 1)}"
        t["p3_b10"], t["p3_b14"], t["p3_b19"] = f(d10["base"]["mean_net"], 1), f(d14["base"]["mean_net"], 1), f(d19["base"]["mean_net"], 0)
        t["p3_mc10"] = " / ".join(f(v, 1) for v in d10["mc"]["15"]["mean_net"]); t["p3_mc14"] = " / ".join(f(v, 1) for v in d14["mc"]["15"]["mean_net"])
        t["p3_mc19"] = " / ".join(f(v, 0) for v in d19["mc"]["15"]["mean_net"])
        t["p3_k"] = f(max(1.0, 1.0 * __import__("math").ceil(max(abs(tor(d, 'k')[1] - tor(d, 'k')[0]) / d["base"]["mean_net"] * 100 for d in (d10, d14, d19)))), 0)
    t["kT20_pct"] = f(-pct(kT20["MWe"], b20E), 0); t["kT7_pct"] = f(-pct(kT7["MWth"], b7H), 0)
    t["nu10_H"] = f(nu10["MWth"], 1); t["b7_H"] = f(b7H, 2); t["b20_E"] = f(b20E, 2); t["b8_E"] = f(b8E, 1); t["kT20_E"] = f(kT20["MWe"], 2)
    t["mc7"] = " / ".join(f(v, 2) for v in m7); t["mc8"] = " / ".join(f(v, 2) for v in m8e); t["mc20"] = " / ".join(f(v, 2) for v in m20e)
    cy = C.get("cy") or []
    if cy:
        hi = max(cy, key=lambda r: r["MWe_20_wells"]); lo = min((r for r in cy if r["T_rock"] == hi["T_rock"]), key=lambda r: r["area_m2"])
        t["cy_hi"] = f(hi["MWe_20_wells"], 1); t["cy_lo"] = f(lo["MWe_20_wells"], 1); t["cy_T"] = f(hi["T_rock"], 0)
        t["cy_Thi"] = f(hi["T_day180"], 0); t["cy_Tlo"] = f(lo["T_day180"], 0)
    return t


PAGE = r"""
<div class="wrap">
<header>
  <div class="eyebrow">aharrsim 2.1 · uncertainty and open-system study · {{date}}</div>
  <h1>How uncertain are the closed loops, and what would the open systems need?</h1>
  <p class="lede">Three questions put to the simulator. How far do the closed-loop results move when gradient, conductivity, insulation and the
  physics outside the casing are varied? Can the fractured open systems be brought up to their quoted capacities by making the fracture network
  far denser and more tortuous? And what happens to those open systems when the uncertainties that mattered for the closed loops are applied to them?</p>
  <nav class="toc" aria-label="Sections">
    <a href="#summary">Summary</a><a href="#method">Method</a><a href="#part1">1 · Closed loops</a><a href="#part2">2 · Open systems</a>
    <a href="#part3">3 · Open systems, uncertain</a><a href="#discussion">Discussion</a><a href="#limits">Limits</a><a href="#repro">Reproduce</a>
  </nav>
</header>

<section id="summary">
  <h2>Summary</h2>
  <div class="tiles">
    <div class="tile"><div class="k">3 km closed loop</div><div class="v">{{tile_cl}}<small>MWth</small></div><div class="s">P10–P90 heat at year 30; no electricity in any run</div></div>
    <div class="tile"><div class="k">20 km closed loop</div><div class="v">{{tile_cl20}}<small>MWe</small></div><div class="s">P10–P90 gross at year 30, against {{b20_E}} MWe at constant conductivity</div></div>
    <div class="tile"><div class="k">Denser fractures alone</div><div class="v">{{tile_dens}}</div><div class="s">net electricity of the 200 °C doublet, 20 → 5,000 fractures at fixed flow and volume</div></div>
    <div class="tile"><div class="k">Gradient ±20 %</div><div class="v">{{p3_T19}}<small>MWe</small></div><div class="s">superhot design built for {{p3_b19}} MWe; conductivity moves it by {{p3_k}} % or less</div></div>
  </div>
  <div class="prose">
    <p><b>1 · Closed loops.</b> Two variables carry the uncertainty: geothermal gradient and rock conductivity, each worth roughly ±20–30 % on heat
    and more on electricity. Rock heat capacity is negligible. Insulation is negligible for the 3 km well and first-order for the 20 km well.
    The largest single effect is one the earlier study left out: rock conductivity falls with temperature, which costs the 3 km well {{kT7_pct}} % of its heat
    and the 20 km well {{kT20_pct}} % of its electricity ({{b20_E}} → {{kT20_E}} MWe). No combination inside the plausible ranges turns a conduction-fed loop into a
    power producer of the size quoted for it. Only a different heat source does: a casing in a convecting or steam-bearing formation lifts the 3 km
    well from {{b7_H}} to {{nu10_H}} MWth and upwards.</p>
    <p><b>2 · Open systems.</b> Fracture density and tortuosity are the wrong levers once fractures are closer than a few tens of metres: the 200 °C doublet
    goes from {{base10}} to {{dens10}} MWe, the Fervo-type doublet from {{base14}} to {{dens14}} MWe, and only the superhot doublet, whose base spacing is too wide, gains
    materially ({{base19}} → {{dens19}} MWe). Tortuosity raises pumping power faster than it raises temperature. What reaches the quoted figures is
    flow through enough rock volume: {{v10}} MWe from 80 kg/s through the existing 0.24 km³ at 200 °C; {{v14}} MWe sustained from 100 kg/s through 0.4 km³ for the
    Fervo type (the published volume holds {{sweep14}} years of heat at the published flow); {{v19}} MWe from 200 kg/s through 0.8 km³ at 375 °C. Dense
    fracturing is what makes those flows sustainable without cooling. It is necessary and not sufficient.</p>
    <p><b>3 · Open systems with the closed-loop uncertainties.</b> The ranking inverts. Conductivity, decisive for closed loops, changes the
    dense open designs by {{p3_k}} % or less, and heat capacity by about as little, because advection through closely spaced fractures does not wait for conduction. Rock temperature dominates:
    a gradient ±20 % moves the three designs across {{p3_T10}}, {{p3_T14}} and {{p3_T19}} MWe (bases {{p3_b10}}, {{p3_b14}}, {{p3_b19}}). Ambient temperature is second for the two sub-200 °C
    systems, and the flow-path unknowns (impedance, channelling, swept fraction) come after. With everything sampled the P10 / P50 / P90 mean net outputs are
    {{p3_mc10}}, {{p3_mc14}} and {{p3_mc19}} MWe: the quoted capacities are reachable, at about the median for the superhot design and nearer the P90 for the other two.</p>
  </div>
</section>

<section id="method">
  <h2>Method</h2>
  <div class="prose">
    <p><b>Closed loops (Part 1).</b> Three cases from the twenty-type study were rerun with the full wellbore–rock solver: the 3 km coaxial retrofit (type 7),
    the multilateral loop (type 8) and the 20 km coaxial loop (type 20). Each uncertain variable was moved to a low and a high value with the rest at base
    (the tornado charts), the main design levers were moved the same way for scale, and flow was swept. Two pieces of physics were added to the rock model:
    conductivity falling with temperature (Vosteen &amp; Schellschmidt 2003, applied layer by layer at the undisturbed temperature) and a convecting zone,
    represented as an effective conductivity Nu × k over a depth interval. A Latin-hypercube Monte Carlo of full simulations then sampled the
    uncertain variables together, with temperature-dependent conductivity on in every run. Electricity is the better of a practical binary plant and a double flash,
    gross of the circulation pump unless stated.</p>
    <p><b>Open systems (Part 2).</b> The three fractured doublets (types 10, 14, 19) were described by a rectangular stimulated block of width W along the wells,
    length L between them and height B, cut by N parallel fractures (the Gringarten solution already in the code). Tortuosity τ lengthens the flow path to τL and
    thins the rock slabs to keep the volume constant. Flow resistance follows I = I<sub>0</sub>[a + (1 − a)(N<sub>0</sub>/N)τ] with half of the base impedance
    in the near-well region (a = 0.5) and the rest shared between fractures; pumping power is Q²I at 70 % efficiency plus any production-pump lift.
    Channelling is a fraction f of fractures carrying the flow and a fraction φ of the height swept. A fast screening version (tabulated plant work, linear
    wellbore heat loss) ran the sweeps; a subset was rerun through the full injector → fractures → producer chain. The metric is net electricity averaged over 30 years.</p>
    <p><b>Open systems under uncertainty (Part 3).</b> One design per system from Part 2 was taken forward and subjected to the Part 1 ranges
    (gradient ±20 %, conductivity 2.0–4.0 W/m/K falling with temperature, heat capacity ±15 %, ambient 15–45 °C) together with the flow-path unknowns
    (impedance ×0.5–×2, active fractures 10–100 %, swept volume 40–100 %), one at a time and in a 3,000-sample Monte Carlo.</p>
  </div>
  <h4>Part 1 ranges</h4>
  {{T_ranges}}
</section>

<section id="part1">
  <h2>1 · Closed-loop uncertainty</h2>
  {{p1_html}}
  <h3 id="p1-opt">Convecting zone: does more flow use it better?</h3>
  <div class="prose wide"><p>If the formation can deliver more heat, the obvious response is to circulate more. For the 3 km well it does not pay:
  doubling the flow to 20 kg/s adds 20–30 % to the heat and takes 25–35 °C off the outlet, so the electricity falls; halving it does the reverse.
  Across a Nusselt number of 10 to 100 and 5 to 20 kg/s the well gives 2–5 MWth and never more than a quarter of a megawatt of electricity. Once the rock is no longer the
  bottleneck the completion is: the cement sheath, the film in the annulus and the heat leaking from the rising water back to the descending water. Larger casing and tubing
  were tried and are worse (lower annular velocity, more leak area). The multilateral loop is different because it has 60 km of open hole in contact with the formation:
  in a mildly convecting formation it reaches the electricity quoted for it, and exceeds it at Nu 30. These are no longer closed-loop results in the usual sense:
  the heat arrives by fluid movement in the formation, and where that permeability exists an open well would produce more from it.</p></div>
  {{T_opt}}
  <div class="callout"><b>Carried into Part 3.</b> Gradient (±20 %, as rock temperature at depth), conductivity (2.0–4.0 W/m/K, falling with temperature),
  heat capacity (±15 %) and ambient temperature. Insulation has no counterpart in an open doublet, whose wells flow fast enough that wellbore loss is 1–3 °C.</div>
</section>

<section id="part2">
  <h2>2 · Open systems: fracture density, tortuosity, flow and volume</h2>
  {{p2_html}}
  <h3 id="p2-cyc">Cyclic single-fracture wells (types 12 and 13)</h3>
  <div class="prose wide"><p>The huff-and-puff well was given the same treatment: the fracture's wall area was raised sixteen-fold at constant stored volume
  (a thinner, wider or multi-wing fracture). In {{cy_T}} °C rock the produced water at day 180 rises from {{cy_Tlo}} to {{cy_Thi}} °C and the electricity from twenty wells from
  {{cy_lo}} to {{cy_hi}} MWe. Area helps, but each well still moves only 200 m³ a day, about 2.3 kg/s averaged; the limit is the volume cycled, and that is set by the
  fracture's compliance, not its area.</p></div>
  {{cyc_figs}}
</section>

<section id="part3">
  <h2>3 · Open systems with the closed-loop uncertainties</h2>
  {{p3_html}}
</section>

<section id="discussion">
  <h2>Discussion</h2>
  <div class="prose">
    <p><b>Closed loops are conduction-limited, so they inherit the rock's uncertainty in full.</b> Heat output scales almost linearly with gradient and a little
    less than linearly with conductivity, and nothing in the completion can compensate. Where the output is electricity the sensitivity is amplified two- to three-fold
    by the plant, because these loops deliver 95–150 °C water, on the steep part of the conversion curve. The spread found here (P10 to P90 about ±15 % of the median
    for heat and ±30–40 % for electricity, with the full tornado ranges about twice that) is what a developer carries until the first well has been logged.</p>
    <p><b>Constant conductivity is optimistic for deep, hot loops.</b> The correction is small at 3 km and large at 20 km, where it removes about two fifths of the
    electricity. It should be the default for anything above about 200 °C.</p>
    <p><b>Open systems are limited by flow and volume, then by temperature.</b> A fracture network needs to be dense enough that the rock between fractures drains
    in the project life, which means a spacing of the order of 10–20 m for 30 years. Beyond that, more fractures only lower the pressure drop. The capacity then is
    (flow) × (specific work at the rock temperature), sustained for as long as the volume lasts; the specific work is about 0.07 MWe per kg/s at 200 °C and 0.29 at 375 °C.
    Tortuosity is a cost: it is thermally equivalent to more fractures and hydraulically worse.</p>
    <p><b>The quoted capacities are flow claims.</b> 3–5 MWe at 200 °C is 80–150 kg/s; 50 MWe at 375 °C is 200 kg/s. Whether one well pair can pass that flow at a
    tolerable pressure drop, through several hundred fractures that all take their share, is the open question, and it is the one the Monte Carlo spreads on.</p>
    <p><b>What to measure first.</b> For a closed loop: temperature at depth and conductivity on core, because nothing else matters. For an open system: temperature at depth,
    then injectivity and the flow profile along the lateral (how many fractures take flow), then tracer-swept volume. Conductivity is not worth measuring for an open system.</p>
  </div>
</section>

<section id="limits">
  <h2>Limits of this study</h2>
  <div class="prose">
    <p>The fracture model is the parallel-plate idealisation: equal flow per active fracture, no heat from rock outside the stimulated block (which adds a little late in life
    and makes the small-volume cases slightly pessimistic), no thermo-elastic change of aperture (cooling tends to open the most-cooled fractures and worsen channelling, which is
    optimistic here). The impedance model is an assumption, anchored only at the base case; the half of it placed near the wells does not fall with fracture count. Flows above about
    65 kg/s assume 9⅝ in production casing, and no account is taken of induced seismicity or of water loss at the injection pressures implied (16 MPa and more for the
    high-flow cases). The convecting-zone cases use an effective conductivity, not a simulation of convection. Monte Carlo distributions are triangular judgements, not data.
    Plant models and reference conditions are those of the twenty-type study.</p>
  </div>
</section>

<section id="repro">
  <h2>Reproduce</h2>
  <div class="prose">
    <p>From the package root: <span class="mono">python3 validation/study_cl.py oat | mc | opt</span> (Part 1; about an hour on two cores),
    <span class="mono">python3 validation/study_open.py screen | confirm | cyclic | part3</span> (Parts 2 and 3; a few minutes each),
    then <span class="mono">python3 report/build_st_report.py</span>. Results are written to <span class="mono">validation/results/ST*.json</span>.</p>
  </div>
  <h3>Sources</h3>
  <ol class="sources">
    <li>Khodayar, M. &amp; Björnsson, S. (2024), Conventional geothermal systems and unconventional geothermal developments: an overview, Open Journal of Geology 14:196 (the twenty development types and their quoted capacities).</li>
    <li>Vosteen, H.-D. &amp; Schellschmidt, R. (2003), Influence of temperature on thermal conductivity, thermal capacity and thermal diffusivity for different types of rock, Physics and Chemistry of the Earth 28:499.</li>
    <li>Gringarten, A. C., Witherspoon, P. A. &amp; Ohnishi, Y. (1975), Theory of heat extraction from fractured hot dry rock, J. Geophys. Res. 80:1120. Stehfest, H. (1970), Comm. ACM 13:47.</li>
    <li>Norbeck, J. H. &amp; Latimer, T. (2023), Commercial-scale demonstration of a first-of-a-kind enhanced geothermal system (Project Red: 63 L/s, 191 °C, 3.5 MWe over a 30-day test).</li>
    <li>DiPippo, R. (2016), Geothermal Power Plants, 4th ed.; Beckers, K. F. &amp; McCabe, K. (2019), GEOPHIRES v2.0 (plant correlations).</li>
  </ol>
  <footer>aharrsim 2.1 · built {{date}} from validation/results</footer>
</section>
</div>

<script src="https://cdnjs.cloudflare.com/ajax/libs/echarts/5.4.3/echarts.min.js"></script>
<script>
const DATA = {{data}};
{{js}}
</script>
"""
