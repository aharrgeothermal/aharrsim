"""Prose and page skeleton of the laboratory-comparison report. Numbers come from the results through build(C, f, REP)."""
from __future__ import annotations

EXTRA_CSS = """
/* temperature field: one hue, stronger where the pack has cooled; the untouched pack recedes into the surface */
:root { --h0: #f3f6fb; --h1: #cde2fb; --h2: #86b6ef; --h3: #3987e5; --h4: #1c5cab; --h5: #0d366b; }
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) { --h0: #1f2730; --h1: #184f95; --h2: #256abf; --h3: #5598e7; --h4: #9ec5f4; --h5: #dcebfb; } }
:root[data-theme="dark"] { --h0: #1f2730; --h1: #184f95; --h2: #256abf; --h3: #5598e7; --h4: #9ec5f4; --h5: #dcebfb; }
.chart.field { height: 560px; max-width: 380px; margin-inline: auto; }
.figrow .chart { height: 340px; }
.figrow .chart.tall { height: 450px; }
.figrow .chart.field { height: 560px; }
.figrow.one .chart { height: 360px; }
.verdicts { display: grid; grid-template-columns: repeat(3, minmax(0, 1fr)); gap: 12px; margin-top: 1.2rem; }
@media (max-width: 900px) { .verdicts { grid-template-columns: minmax(0, 1fr); } }
.verdict { background: var(--surface); border: 1px solid var(--rule); border-radius: 8px; padding: 14px 16px; min-width: 0; }
.verdict h4 { margin: 0 0 4px; font-size: 1.02rem; }
.verdict p { margin: 0.35rem 0 0; font-size: 0.93rem; color: var(--ink-2); }
.eq { font-family: var(--f-mono); font-size: 0.9rem; background: var(--surface); border: 1px solid var(--rule); border-radius: 6px; padding: 10px 14px; margin: 0.8rem 0;
      max-width: 78ch; overflow-wrap: anywhere; }
dl.kv { max-width: 86ch; }
"""


def build(C, f, REP):
    t = {}
    t["x_rep"] = f(C["x_rep"], 2); t["x_pub"] = f(C["x_pub"], 2); t["x_fit"] = f(C["x_fit"], 2); t["x_unc"] = f(C["x_unc"], 2)
    t["on_pub"] = f(C["on_pub"], 0); t["on_fit"] = f(C["on_fit"], 0); t["on_unc"] = f(C["on_unc"], 0)
    t["K_fit"] = f(C["K_fit"] * 1e10, 1); t["a_fit"] = f(C["a_fit"], 1) if C.get("a_fit") else "–"
    t["maxdev"] = f(C["maxdev"], 1); t["maxdev_a"] = f(C["maxdev_a"], 1)
    t["rec_onset"] = f(C["rec_onset"], 0); t["s12_onset"] = f(C["s12_onset"], 0)
    t["D30"] = f(C["D30"], 0); t["D30_10"] = f(C["D30_10"], 0); t["D30_box"] = f(C["D30_box"], 0)
    t["V_th"] = f(C["V_th"], 2)
    t["E_plateau"] = f(4.2e6 * 0.5 / 3600.0 * 39.0 * 30.0 * 60.0 / 1e6, 0)
    t["tile_x"] = f(C["x_rep"], 2)
    t["rhoc_b"] = f(C["rhoc_b"], 2); t["lam"] = f(C["lam"], 1); t["ts_cyl"] = f(C["ts_cyl"], 0)
    t["avg_pub"] = f(C["avg_pub"], 1); t["avg_fit"] = f(C["avg_fit"], 1); t["avg_unc"] = f(C["avg_unc"], 1)
    t["T60_pub"] = f(C["T60_pub"], 0); t["T60_fit"] = f(C["T60_fit"], 0)
    t["on_visc"] = f(C["on_visc"], 0); t["on_buoy"] = f(C["on_buoy"], 0); t["on_a2"] = f(C["on_a2"], 0); t["on_a4"] = f(C["on_a4"], 0)
    t["on_glass"] = f(C["on_glass"], 0); t["on_gap"] = f(C["on_gap"], 0); t["d_gap"] = f(C["on_gap"] - C["on_pub"], 0)
    t["E_60"] = f(C["E_60"], 0); t["E_place"] = f(C["E_place"], 0); t["pct_out"] = f(C["pct_out"], 0)
    t["r_40"] = f(C["r_40"], 0); t["r_speed"] = f(C["r_speed"], 1)
    t["f_first"] = C["f_first"]; t["f_down"] = f(C["f_down"], 1); t["f_60"] = f(C["f_60"], 0); t["fp_60"] = f(C["fp_60"], 0); t["fp_down"] = f(C["fp_down"], 1)
    t["rec_E"] = f(C["rec_E"], 0); t["nom_E"] = f(C["nom_E"], 0); t["s12_E"] = f(C["s12_E"], 0); t["rec_kW"] = f(C["rec_kW"], 1); t["nom_kW"] = f(C["nom_kW"], 1)
    t["rec_theta60"] = f(100 * C["rec_theta60"], 0); t["nom_theta60"] = f(100 * C["nom_theta60"], 0)
    t["xbox99"] = f(C["xbox99"], 2); t["xbox90"] = f(C["xbox90"], 2); t["xbox50"] = f(C["xbox50"], 2)
    t["xun99"] = f(C["xun99"], 2); t["xun90"] = f(C["xun90"], 2); t["xun50"] = f(C["xun50"], 1); t["px"] = f(C["px"], 2)
    t["b_T1"] = f(C["b_T1"], 0); t["b_T30"] = f(C["b_T30"], 0); t["d_T1"] = f(C["d_T1"], 0); t["d_T10"] = f(C["d_T10"], 0); t["d_T30"] = f(C["d_T30"], 0)
    t["d_first"] = f(C["d_first"], 1); t["g_T30"] = f(C["g_T30"], 0); t["D_holds"] = f(C["D_holds"], 0); t["box_bt"] = f(C["box_bt"], 1); t["box_half"] = f(C["box_half"], 1)
    t["d_MWe1"] = f(C["d_MWe1"], 2); t["d_MWe30"] = f(C["d_MWe30"], 2); t["d_H1"] = f(C["d_H1"], 1); t["d_H30"] = f(C["d_H30"], 1); t["b_first"] = f(C["b_first"], 1)
    t["ib_on"] = f(C["ib_on"], 0); t["sb_on"] = f(C["sb_on"], 0); t["ib_avg"] = f(C["ib_avg"], 1); t["sb_avg"] = f(C["sb_avg"], 1)
    t["fb_a13"] = f(100 * C["fb_a13"], 0); t["fb_b13"] = f(100 * C["fb_b13"], 0); t["fb_a14"] = f(100 * C["fb_a14"], 0); t["fb_b14"] = f(100 * C["fb_b14"], 0)
    t["fb_nc"] = f(100 * C["fb_nc"], 0); t["fb_el"] = f(100 * C["fb_el"], 0); t["fb_a12_2"] = f(100 * C["fb_a12_2"], 0)
    t["fb_a13_first"] = f(C["fb_a13_first"], 1); t["fb_b13_first"] = f(C["fb_b13_first"], 1)
    t["rate_fit"] = f(C["rate_fit"], 2)
    return t


PAGE = r"""
<div class="wrap">
<header>
  <div class="eyebrow">aharrsim 2.3 · laboratory comparison · {{date}}</div>
  <h1>A single-well EGS in a sandbox, set against the simulator</h1>
  <p class="lede">Ma and co-authors (2026) built a one-well enhanced geothermal system in a 1.5 m sandbox and circulated water through it for an hour at a time.
  This page sets their measurements against what aharrsim gives for the same arrangement: the type 11 example of the benchmark, the swept-rock chart,
  and a model of the sandbox itself, which was added for the purpose.</p>
  <nav class="toc" aria-label="Sections">
    <a href="#summary">Summary</a><a href="#study">The experiment</a><a href="#ours">Where it sits in the benchmark</a><a href="#model">The sandbox in the model</a>
    <a href="#optimum">The recommended settings</a><a href="#field">From sandbox to field</a><a href="#updown">Injection above or below</a>
    <a href="#changes">Changes to the simulator</a><a href="#limits">Limits</a>
  </nav>
</header>

<section id="summary">
  <h2>Summary</h2>
  <div class="tiles">
    <div class="tile"><div class="k">Flushed share at breakthrough</div><div class="v">{{tile_x}}</div><div class="s">measured, the same at 0.3, 0.5 and 0.6 m³/h. Model {{x_pub}} with the published properties, {{x_fit}} with one property fitted</div></div>
    <div class="tile"><div class="k">Heat taken on the plateau</div><div class="v">{{E_plateau}}<small>MJ</small></div><div class="s">at 0.3, 0.5 and 0.6 m³/h alike: flow sets how fast the heat is taken, not how much</div></div>
    <div class="tile"><div class="k">Recommended settings</div><div class="v">{{rec_onset}}<small>min</small></div><div class="s">plateau in the model; the paper compares settings at 60 minutes, so the optimum is the test that had only just reached breakthrough</div></div>
    <div class="tile"><div class="k">Spacing for 30 years</div><div class="v">{{D30_10}}–{{D30_box}}<small>m</small></div><div class="s">between the two zones at 10 kg/s, from a large block to a bounded one; 0.9–1.2 m in the sandbox</div></div>
  </div>
  <div class="verdicts">
    <div class="verdict"><h4>Does it agree with what we simulated?</h4>
      <p>Yes. The sandbox is the benchmark's type 11 at about 1 : 150 and behaves as the swept-rock chart says: a plateau at rock temperature, then a decline, with the
      plateau set by the heat capacity flushed and not by time or flow separately. The paper's reading of its own results as "optimal settings" does not carry to another scale.</p></div>
    <div class="verdict"><h4>Did the simulation need updating?</h4>
      <p>No result changes. The simulator had no model of flow between two screens on one well, so two were added: a numerical block model (flow in a porous block with
      the weight and viscosity of the cold water) and a dipole element for the well solver. The type 11 example was run again with them.</p></div>
    <div class="verdict"><h4>What is new for us?</h4>
      <p>A measured sweep for a single well: breakthrough at {{x_rep}} of the heat capacity between the screens. And a point the paper does not raise: it injects above the
      producer, where cold water sinks towards it. At field scale the model makes that arrangement costly once the vertical permeability reaches 10⁻¹³ m².</p></div>
  </div>
  <div class="prose">
    <p><b>The experiment.</b> A closed, insulated steel cylinder 1.5 m high and 1 m across, packed with 1 mm glass beads and filled with water at 50–70 °C. One well on the axis:
    cold water leaves through an upper screen, sweeps the pack and returns through a lower screen at the base, 90 or 120 cm away, at 0.2–0.6 m³/h. The beads are small enough that
    water and solid are always at the same temperature, so this is swept rock in its sharpest form.</p>
    <p><b>Agreement.</b> The reported plateaus of 50, 30 and 25 minutes at 0.3, 0.5 and 0.6 m³/h are one and the same flushed volume, {{V_th}} m³ of bead pack, or {{x_rep}} of the
    heat capacity between the screens. A wider spacing delays the decline, and a higher flow brings it forward in proportion. All of this is the behaviour of the
    swept-rock chart of the benchmark.</p>
    <p><b>The model of the sandbox.</b> With the published properties the model leaves the plateau after {{on_pub}} minutes where the experiment holds for 30, the same sign as the
    difference the authors report for their own model. The cause is the weight of the cold water: at the published permeability it sinks towards the
    production screen faster than the measurements allow. The permeability was not measured. Lowering it to {{K_fit}} × 10⁻¹⁰ m², so as to match one number (the mean production
    temperature of the nominal case), gives a plateau of {{on_fit}} minutes and reproduces the other reported means within {{maxdev}} K. Keeping the published
    mean and making the horizontal permeability {{a_fit}} times the vertical does the same (within {{maxdev_a}} K).</p>
    <p><b>The recommended settings.</b> The paper recommends 0.4 m³/h, 31 °C injection and 120 cm spacing. In the model that setting holds its plateau for {{rec_onset}} minutes, so a
    comparison at 60 minutes sees almost no decline. Run for three hours it declines like the others. The rule that carries to any scale is one number: keep the flushed share below the
    breakthrough value for the life required. Flow trades power against life one for one, injection temperature scales the heat without moving the breakthrough, and spacing sets the volume.</p>
    <p><b>At field scale.</b> A well with zones 150 m apart at 10 kg/s is at first arrival after {{d_first}} years in a large fractured block and after {{box_bt}} years in a block
    bounded like the sandbox. Thirty years needs the zones {{D30_10}} m apart in a large block (10 % of the rise lost) and {{D30_box}} m in a bounded one (none lost). The benchmark's type 11 example (three fractures, 0.002 km³) cools within years for the
    same reason the sandbox cools within the hour, and its conclusion stands: the output of a single well is small, about {{d_H1}} MWth and {{d_MWe1}} MWe at 183 °C.</p>
  </div>
</section>

<section id="study">
  <h2>The experiment</h2>
  <div class="prose">
    <p>Ma, Li, Chen, Gao, Han, Lei and Wang, "Experimental investigation on heat transfer characteristics of single-well enhanced geothermal system", Geothermal Energy 14:9,
    published February 2026 (PowerChina Northwest Engineering, Xi'an, and the Guangzhou Institute of Energy Conversion). The news item that prompted this note summarises the same paper.</p>
  </div>
  <dl class="kv">
    <dt>Vessel</dt><dd>Stainless steel cylinder, 1.5 m high, 1 m in diameter, 3 mm plate, insulated outside; closed</dd>
    <dt>Medium</dt><dd>Glass beads 0.8–1.2 mm; porosity 0.35; 1,500 kg/m³, 800 J/kg/K, 2.5 W/m/K; permeability 8.92 × 10⁻¹⁰ m² (Kozeny–Carman estimate, not measured)</dd>
    <dt>Well</dt><dd>On the axis. Two screens 10 cm long and 16 mm across: injection through the upper one, on a coaxial pipe through the lid; production through the lower one, fixed at the base</dd>
    <dt>Variables</dt><dd>Flow 0.2–0.6 m³/h; pack 50, 60, 70 °C; injection 21, 31, 41 °C; screen spacing 90 and 120 cm. Each test ran two hours and was repeated three times;
    nothing later than 70 minutes is reported</dd>
    <dt>Sensors</dt><dd>44 PT100: two vertical lines with sensors 10 cm apart and one line across the diameter, 55 cm above the base, with sensors 8 cm apart; electromagnetic flowmeter</dd>
    <dt>Heat flow</dt><dd>C ρ q (T<sub>out</sub> − T<sub>in</sub>) with C ρ = 4.2 MJ/m³/K (the reported means reproduce to 0.1 % with this constant)</dd>
  </dl>
  {{T_rep}}
  <div class="callout"><b>How firm are the numbers?</b> The plateaus, the means over 60 minutes and the heat flows are consistent with each other. The temperature at 60 minutes
  of the nominal case is not: the text gives 37 °C, 42 °C and 46 °C for it in three places. It is treated here as 42 °C, known to about ±5 K.</div>
</section>

<section id="ours">
  <h2>Where it sits in the benchmark</h2>
  <div class="prose">
    <p><b>The type.</b> Type 11 of the benchmark is the single-well EGS: cold water down the annulus into an upper zone, through fractured rock, back through a lower zone and an
    insulated inner pipe. The example is a 4.8 km well in 183 °C rock at 10 kg/s, with the rock described as three fractures 30 m apart giving a flow path of 150 m
    (0.002 km³ swept uniformly). It gives {{b_T1}} °C at the wellhead in year 1 and {{b_T30}} °C in year 30. The sandbox is the same arrangement with a spacing of about one metre.</p>
    <p><b>The regime.</b> The benchmark divides the types by how the water meets the rock. Here it sweeps it. The swept-rock chart has two groups:</p>
  </div>
  <div class="eq">x = t / t<sub>s</sub>,&nbsp; t<sub>s</sub> = V (ρc)<sub>b</sub> / (Q (ρc)<sub>w</sub>)&nbsp;&nbsp;·&nbsp;&nbsp; y = (S / 2) / √(α t)</div>
  <div class="prose">
    <p>x is the share of the block's heat capacity that the water has had the chance to carry off; y compares half the spacing of the flow paths with the distance heat conducts in
    the time. In a pack of 1 mm beads a grain gives up its heat in under a second, so y is zero for every practical purpose and the chart predicts a sharp front: full
    temperature until the front arrives, with an ideal uniform sweep at x = 1.</p>
    <p><b>What the experiment adds.</b> The flow is not uniform. It leaves one short screen and converges on another, so the water on the shortest path arrives first.
    The three reported plateaus give the same x = {{x_rep}}, taking V as the cylinder of bead pack between the screen centres ({{ts_cyl}} minutes to flush at 0.5 m³/h with a
    heat capacity of {{rhoc_b}} MJ/m³/K). That is a measured sweep for a single well in a block about as wide as the screens are apart: a little over half of an ideal sweep.
    For comparison, the production record of Project Red examined in the manuscript put the product of active and swept fractions at 0.43–0.53.</p>
  </div>
</section>

<section id="model">
  <h2>The sandbox in the model</h2>
  <div class="prose wide">
    <p>The block model added to the simulator (<span class="mono">aharrsim.porous</span>) solves Darcy flow and heat transport in a cylinder round a well, in radius and height, with the
    water and the solid at one temperature. The flow is solved again as the pack cools, because cold water is both heavier and more viscous than hot. It was given the box, the
    screens and the bead pack as published (heat capacity {{rhoc_b}} MJ/m³/K with 1,500 kg/m³ read as the density of the pack; conductivity {{lam}} W/m/K; dispersivity one bead diameter)
    and run in three forms: with the published permeability, with the flow not affected by temperature, and with the permeability fitted to one reported number.</p>
  </div>
  {{F_typ}}
  <div class="prose wide">
    <p><b>Published properties.</b> The plateau ends at {{on_pub}} minutes against 30 reported, and the 60-minute mean is {{avg_pub}} °C against 55.3 °C. The plateau no longer scales with
    flow: at 0.2 m³/h the cold water sinks under its own weight and arrives in half the time that the scaling gives, where the experiment reports almost no change in the hour.</p>
    <p><b>Flow not affected by temperature.</b> Plateau {{on_unc}} minutes, the same flushed share ({{x_unc}}) at every flow, mean {{avg_unc}} °C. Closer, and still early.</p>
    <p><b>What the two effects do.</b> Taken alone, the viscosity of the cold water lengthens the plateau to {{on_visc}} minutes. Water at 21 °C is twice as viscous as water at 60 °C,
    a mobility ratio of about 0.5: a favourable displacement, which slows the shortest paths and evens out the sweep. Taken alone, its weight shortens the plateau to {{on_buoy}} minutes,
    because the injection screen is above the production screen and the cold water slumps towards it. The measurements sit close to the viscosity-only result, so the weight of the
    cold water acts far less in the experiment than the published permeability implies.</p>
    <p><b>One property fitted.</b> The strength of the buoyant flow is proportional to the vertical permeability, which the paper estimates from a formula and does not measure
    (the formula gives 5.6 × 10⁻¹⁰ m² for 1 mm beads at the stated porosity, not the 8.92 × 10⁻¹⁰ m² printed). An isotropic permeability of {{K_fit}} × 10⁻¹⁰ m² matches the mean of the nominal
    case; so does the published mean with the horizontal permeability {{a_fit}} times the vertical. Finer or tighter packing, layering of the pack as it was filled, or air released
    from the hot water and trapped in the pores would each lower it; the data cannot say which. With {{K_fit}} × 10⁻¹⁰ m² the plateau is {{on_fit}} minutes, the temperature at 60 minutes
    {{T60_fit}} °C, and every other series of the paper follows without further adjustment:</p>
  </div>
  {{T_cmp}}
  {{T_decay}}
  <h3>Inside the box</h3>
  <div class="prose wide">
    <p>The cold water spreads outwards from the upper screen as a body about 40 cm in radius after half an hour, with a tongue drawn down the well towards the production screen.
    It then fills the box downwards. Along the paper's line of sensors across the box, 55 cm above the base, model and report agree. After 40 minutes a drop of 20 K has reached
    {{r_40}} cm from the axis: past the sensors at about 36 cm and short of the outermost at about 44 cm, where the paper reports a 20 K drop over more than 80 % of the monitored
    distance and little change at the ends of the line. It moves outwards at {{r_speed}} cm/min (reported: about 0.75).</p>
    <p>Along the vertical line, 10 cm from the wall, they differ. The first cooling arrives at about the reported time and depth, 20–30 minutes and the level of the injection screen.
    From there the model cools the line down to {{f_60}} cm of depth within the hour, where the paper has the front moving down at 1.0 cm/min, to about 80 cm, and the pack below 110 cm
    untouched. No permeability tried changes that much.</p>
    <p>The paper's account of this line is hard to square with its own line across the box, which has the pack 20 K colder at 95 cm of depth after 40 minutes, and with its heat flows,
    by which {{pct_out}} % of the heat in the box has left after an hour ({{E_60}} of {{E_place}} MJ). It also gives the position of the line twice, differently. The production
    temperature integrates the whole flow and is the firmer test; the comparison along the vertical line is left open until the figures can be read.</p>
  </div>
  {{F_field}}
  {{T_fronts}}
  <h3>What else moves the model</h3>
  {{T_sens}}
  <div class="prose wide">
    <p>Heat capacity is not the explanation: glass at its usual 2,500 kg/m³ moves the plateau by a minute. Reading the 90 cm as the gap between the screens adds {{d_gap}} minutes.
    Conductivity and dispersion do not matter over an hour. The properties that move the result are those that govern the vertical flow.</p>
  </div>
</section>

<section id="optimum">
  <h2>The recommended settings</h2>
  <div class="prose wide">
    <p>The paper concludes that 0.4 m³/h, 31 °C injection and 120 cm spacing is the best compromise between heat flow and stability of the production temperature, and limits
    that conclusion to the conditions tested. Its tests ran for two hours; the settings are compared by the decline after 60 minutes and the means over that hour. The three choices
    do different things, and none of them is an optimum that survives a change of scale or of the time at which the comparison is made.</p>
  </div>
  {{F_rec}}
  {{T_rec}}
  <div class="prose wide">
    <ul class="tight">
      <li><b>Flow.</b> The heat taken before breakthrough is fixed by the heat capacity swept: {{E_plateau}} MJ in the nominal geometry at 0.3, 0.5 and 0.6 m³/h alike.
      A lower flow takes the same heat more slowly. Compared at 60 minutes, that looks like stability.</li>
      <li><b>Injection temperature.</b> The heat flow is proportional to the difference between pack and injection temperature. The paper's measure of decline, read as a percentage
      of the initial temperature in °C, falls when the injection is warmer simply because there is less difference to lose. As a share of the rise nothing is gained: by the reported
      means the produced water keeps 0.88 of its rise over the hour at 21 °C and 0.87 at 41 °C, and the model gives the same plateau at both. Raising the injection from 21 to 31 °C
      gives up a quarter of the heat for no gain in life.</li>
      <li><b>Spacing.</b> This is the one real lever: it adds swept volume. In the model 120 cm lengthens the plateau of the nominal case from {{on_fit}} to {{s12_onset}} minutes.
      It is limited in the sandbox by the height of the box and in the field by the length of well that can be placed in hot, permeable rock.</li>
    </ul>
    <p>The recommended combination has a plateau of {{rec_onset}} minutes in the model, at a heat flow of {{rec_kW}} kW against {{nom_kW}} kW for the nominal case, and keeps {{rec_theta60}} % of
    its rise at 60 minutes against {{nom_theta60}} %. What the paper found is the setting that had only just reached breakthrough at the hour of the comparison. The design rule that can be taken from the same data is</p>
  </div>
  <div class="eq">Q τ (ρc)<sub>w</sub> / (V (ρc)<sub>b</sub>) &lt; x at breakthrough</div>
  <div class="prose wide"><p>for a life τ at full temperature, with V the block between the screens and x at breakthrough {{x_rep}} as measured for this geometry. It is the volume rule of the
  benchmark, V ≥ ṁ c̄ τ / (0.8 (ρc)<sub>r</sub>), with the 0.8 of an ideal sweep replaced by the sweep a single well achieves.</p></div>
</section>

<section id="field">
  <h2>From sandbox to field</h2>
  <div class="prose wide">
    <p>With heat carried by the water alone, the decline depends on time only through x = t Q (ρc)<sub>w</sub> / (D³ (ρc)<sub>b</sub>), where D is the spacing of the two zones, and on
    the shape of the block. Two shapes bracket the field case. In the sandbox the block is a closed cylinder barely wider than the spacing. In a large fractured block the water
    is free to spread, which is the classical source and sink pair; its streamlines can be integrated once and for all, and that is the element added to the well solver
    (<span class="mono">reservoir.SingleWellDipole</span>).</p>
  </div>
  {{F_curves}}
  <div class="prose wide">
    <p><b>Two shapes, two declines.</b> The bounded block reaches first arrival sooner (x = {{xbox99}} against {{xun99}}) because the wall keeps all of the flow close to the well, and then
    collapses: half the rise is gone at x = {{xbox50}}, when the cylinder between the screens has been flushed once. The open block arrives later and declines slowly, as (0.73 / x)<sup>¼</sup>:
    half the rise is still there at x = {{xun50}}, because the flow keeps reaching fresh rock farther out. The reported decline of the sandbox, about 0.8 °C/min after the plateau
    ({{rate_fit}} in the fitted model), is of the bounded kind. Which kind a field well shows depends on whether the permeable rock ends near the well.</p>
    <p><b>The type 11 example.</b> The same well was run with four descriptions of the rock at its toe. The benchmark's three fractures hold {{b_first}} years of heat and give {{b_T1}} → {{b_T30}} °C
    at the wellhead. Zones 150 m apart in a large block arrive after {{d_first}} years and give {{d_T1}} → {{d_T30}} °C over 30 years, {{d_H1}} → {{d_H30}} MWth and {{d_MWe1}} → {{d_MWe30}} MWe.
    A block bounded like the sandbox would arrive after {{box_bt}} years and lose half its rise by year {{box_half}}. Holding 90 % of the rise for 30 years needs the zones {{D_holds}} m apart
    ({{g_T30}} °C at the wellhead in year 30). In every case the well starts at 4–5 MWth and under half a megawatt of electricity: the benchmark's reading of the type as small is unchanged.</p>
  </div>
  {{T_t11}}
  <h3>Spacing a single well needs</h3>
  {{T_space}}
</section>

<section id="updown">
  <h2>Injection above or below</h2>
  <div class="prose wide">
    <p>The sandbox injects through the upper screen and produces through the lower, as the single-well EGS concept and the type 11 example do. Cold water is heavier than hot. Injected
    above the producer it sinks towards it; injected below, it stays down and the producer draws on the hot rock above. How much this matters depends on the ratio of the buoyant flow,
    K Δρ g / μ, to the flow that is pumped.</p>
    <p><b>In the sandbox</b> the ratio is about one at the published permeability and a fifth of that at the fitted one. The model gives a plateau of {{sb_on}} minutes with injection
    through the lower screen against {{on_pub}} with the published permeability, and {{ib_on}} against {{on_fit}} with the fitted one. The experiment did not try it; swapping the two
    screens would test both the effect and the permeability.</p>
    <p><b>At field scale</b> the density difference is six times larger (60 °C water against 183 °C rock), so the same question turns on the vertical permeability of the fractured zone.
    For zones 150 m apart at 10 kg/s the ratio is about 0.1 at 10⁻¹⁴ m², one at 10⁻¹³ m² and ten at 10⁻¹² m², with the pumped flow taken half-way between the zones. After 30 years the produced water keeps {{fb_a14}} % of its rise with injection above and {{fb_b14}} % with injection below at 10⁻¹⁴ m²;
    at 10⁻¹³ m² the figures are {{fb_a13}} % and {{fb_b13}} %. At 10⁻¹² m² injection above short-circuits within two years. By a line-source estimate, driving 10 kg/s between two 20 m zones takes roughly 25 MPa at 10⁻¹⁴ m² and 2.5 MPa at 10⁻¹³ m².
    Rock permeable enough to take the flow at a workable pressure is therefore rock in which the arrangement of the sandbox is the unfavourable one.</p>
  </div>
  {{F_updown}}
  {{T_updown}}
  <div class="callout"><b>Read with care.</b> The field block is uniform, isotropic and closed, and is at one temperature before injection. A real stimulated zone is none of these, and its vertical
  permeability is the unknown that decides the result. The direction of the effect is robust; its size is not established.</div>
</section>

<section id="changes">
  <h2>Changes to the simulator</h2>
  <div class="prose wide">
    <ul class="tight">
      <li><span class="mono">aharrsim.porous.AxisymmetricBlock</span>: Darcy flow and heat transport in a porous or densely fractured cylinder round a well, with screens, buoyancy,
      temperature-dependent viscosity, anisotropic permeability and closed boundaries. Second-order upwind in space and time, so a thermal front stays sharp.</li>
      <li><span class="mono">aharrsim.reservoir.SingleWellDipole</span>: two zones on one well in a large block, as an element for the well solver. Sharp front on every streamline,
      or the parallel-fracture response of Gringarten along every streamline when the flowing fractures are a given distance apart.</li>
      <li><span class="mono">validation/segs_sandbox.py</span>: the checks, the experiment, the fit and the field cases of this page.</li>
      <li>No existing result changes. The manuscript and the benchmark tables stand as they are.</li>
    </ul>
  </div>
  {{T_checks}}
</section>

<section id="limits">
  <h2>Limits</h2>
  <div class="prose">
    <p>The comparison rests on the numbers in the text of the paper. Its table of operating conditions and its figures could not be read, so the default values of each series are
    inferred and no measured curve was digitised. Several reported numbers disagree with each other, as listed above.</p>
    <p>The fitted permeability is one number chosen to match one mean. It shows that a weaker buoyant flow reconciles model and experiment; it does not show why the flow is weaker.
    The model is axisymmetric, so any sinking plume is forced onto the axis.</p>
    <p>The sandbox is a porous medium with water and solid at one temperature. Fractured rock with flow paths metres apart cools earlier along each path; the dipole element has that
    option, and the benchmark's chart carries it as y. The sandbox says nothing about injectivity, about how much of a stimulated zone takes flow, or about the wellbore.</p>
    <p>The field cases assume the permeable rock is uniform. The spacings in the table are for full sweep of that rock by the dipole pattern and scale with the cube root of flow times life.</p>
  </div>
  <h3>Reproduce</h3>
  <div class="prose">
    <p>From the package root: <span class="mono">python3 validation/segs_sandbox.py checks | lab | calib | rec | field</span> (about 25 minutes in all), then
    <span class="mono">python3 report/build_sg_report.py</span>. Results are written to <span class="mono">validation/results/SG_*.json</span>.</p>
  </div>
  <h3>Sources</h3>
  <ol class="sources">
    <li>Ma, B., Li, A., Chen, W., Gao, X., Han, W., Lei, Y. &amp; Wang, L. (2026), Experimental investigation on heat transfer characteristics of single-well enhanced geothermal system,
    Geothermal Energy 14:9, <a href="https://link.springer.com/article/10.1186/s40517-026-00379-4">doi:10.1186/s40517-026-00379-4</a>.</li>
    <li>"Lab study reveals optimal settings for single-well geothermal heat extraction", <a href="https://bioengineer.org/lab-study-reveals-optimal-settings-for-single-well-geothermal-heat-extraction/">Bioengineer.org</a>, 8 November 2026 (news summary of the paper).</li>
    <li>Gringarten, A. C., Witherspoon, P. A. &amp; Ohnishi, Y. (1975), Theory of heat extraction from fractured hot dry rock, J. Geophys. Res. 80:1120.</li>
    <li>Khodayar, M. &amp; Björnsson, S. (2024), Conventional geothermal systems and unconventional geothermal developments: an overview, Open Journal of Geology 14:196 (type 11).</li>
  </ol>
  <footer>aharrsim 2.3 · built {{date}} from validation/results</footer>
</section>
</div>

<script src="https://cdnjs.cloudflare.com/ajax/libs/echarts/5.4.3/echarts.min.js"></script>
<script>
const DATA = {{data}};
{{js}}
</script>
"""
