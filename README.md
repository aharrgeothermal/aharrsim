# aharrsim — transient wellbore–reservoir–plant simulator for geothermal wells

Version 2.4 (November 2026). The core is a transient heat-transfer and hydraulics model of a vertical
(or deviated / L-shaped) coaxial closed loop, 1–12+ km deep:

* cold water is pumped **down the annulus** and returns hot **up the centre tubing** (or the reverse,
  `flow_direction: tubing_down`, which the report shows delivers about a quarter less heat);
* no fluid exchange with the rock — heat arrives by **conduction only**;
* working fluid is **water/steam** (IAPWS-95), including flashing in the riser when the
  wellhead pressure is below saturation.

The main output is the **heat delivered to surface**, `Q = m (h_out − h_in)`, with the outlet
temperature, steam quality, pressures, pumping requirement, exergy and depth profiles.

Since v2.0 the same solver also runs production, injection and U-tube wells, chains of wells and reservoir elements, and
surface plants, so that one example of every geothermal development type can be run on a common basis (below).

The simulator accompanies the paper "Swept rock or touched rock: a common-basis benchmark of twenty geothermal development types and an
audit of their quoted capacities" (Wrobel, McGeer and Sellar); `paper/make_figs.py` draws its figures from `validation/results/`.
The change history is in `CHANGELOG.md`.

## Heat delivered (reference design)

10 kg/s injected at 40 °C, 30 °C/km from 15 °C, rock k 2.8 W/m/K, 3½″ × 4½″ vacuum-insulated tubing to TD
(insulation k 0.03 W/m/K), 7″ liner below 3 km, 5 MPa at the wellhead:

| Depth | Rock T at TD (°C) | Heat, 1 yr (MWth) | 10 yr | 30 yr | Outlet T, 10 yr (°C) | 30-yr total (GWh) |
|---|---|---|---|---|---|---|
| 4 km | 135 | 0.43 | 0.36 | 0.33 | 49 | 95 |
| 6 km | 195 | 1.13 | 0.95 | 0.88 | 64 | 249 |
| 8 km | 255 | 2.02 | 1.72 | 1.61 | 82 | 452 |
| 10 km | 315 | 3.07 | 2.65 | 2.48 | 105 | 694 |
| 12 km | 375 | 4.22 | 3.67 | 3.46 | 129 | 963 |

Sensitivities (gradient 20–60 °C/km, insulation, rock conductivity, flow rate, larger completion, steam flashing) and the
number of wells needed for a 5 GW field (heat and power, including a best-case design and a check with temperature-dependent
rock conductivity) are in the report.

## Beyond the coaxial loop (v2.0)

The same solver now runs the other well configurations, and two further modules supply the reservoir
and the surface plant, so every development type in Khodayar & Björnsson's (2024) classification can be
simulated (`validation/dev_types.py` runs one example of each of the 20 types; results in the report):

| Module | What it adds | Benchmark |
|---|---|---|
| `mode: production` / `injection` | one flow channel (casing or tubing with a static annulus), feed zone at TD by temperature or enthalpy, flashing riser, downhole pump (ESP / LSP), productivity index and wellhead choke margin, bottom-pressure or wellhead-pressure boundary | exact steady solutions (2nd order), Ramey / Hasan–Kabir transient (≤1.5 K), measured flowing two-phase wells (Az-18 −1.9 %, M-90 +0.3 %) |
| `mode: utube` | single U-tube borehole heat exchanger, both legs coupled to the borehole wall through the multipole delta circuit | pygfunction g-function + multipole: wall temperature within 0.1 K, fluid within 0.2 K |
| `far_field: adiabatic` | unit-cell outer boundary for well arrays, lateral stacks and finite heat stores | – |
| `bottom_element` | reservoir element in series at TD of a coaxial well (single-well EGS) | energy balance |
| `aharrsim.systems` | chains of segments (injector → laterals → producer, injector → reservoir → producer) with transient hand-over; circulation power from the pressure balance round the loop (v2.2) | GeoCLUSTER `utube` database (SBT): 0.1–0.7 K over 1–40 years, 2–20 km laterals; GEOPHIRES-X on ten whole systems (v2.2) |
| `aharrsim.reservoir` | convective reservoir; Gringarten–Sauty doublet breakthrough; Gringarten 1975 parallel-fracture EGS (as GEOPHIRES; de Hoog inversion since v2.2); cyclic single-well fracture (Lauwerier-checked); lumped heat store | 30-digit mpmath inversion to 1e−9 (50 m spacing) and 4e−6 (sharp front); erfc limit to 1e−9; thin-slab limit to 0.003; breakthrough time exact; Lauwerier within 0.8 K |
| `aharrsim.plant` | ORC cycle (six fluids, recuperator), single/double flash, dry steam, heat pump, exergy; GEOPHIRES practical correlations | GEOPHIRES utilisation efficiencies: binary within ±11 %, flash within 1–3 % (35 °C condenser) |

```bash
cd validation
python verify_ext.py        # V8-V11 and B4 above (~25 min; needs pygfunction and mpmath)
python dev_types.py         # the 20 development-type examples (~15 min on 2 cores)
python dev_types_table.py   # summary table with the electricity assessment
```

Results of the examples as of v2.2 (15 °C ambient; GEOPHIRES practical binary for liquid water, flash or steam turbine for two-phase
and steam, 40 °C condenser; heat relative to each case's return temperature; heat pumps consume electricity; `python
validation/dev_types_table.py md` prints the table):

| # | Type (example) | Horizon | Water T first → last (°C) | MWth | Gross MWe first → last | Pumps MWe first → last | Net MWe | Quoted |
|---|---|---|---|---|---|---|---|---|
| 1 | Low-temperature / hot sedimentary aquifer doublet (Paris Basin type) | year 30 | 69 → 63 | 6.50 | 0.00 | 0.10 | 0.00 | <2 MWe |
| 2 | Medium-temperature pumped well + binary plant (Nevada type) | year 30 | 161 → 147 | 24.36 | 3.00 → 2.24 | 0.38 → 0.40 | 1.83 | 2–6 MWe (<7) |
| 3 | High-temperature self-flowing well, flash plant (Wairakei type) | year 30 | 175 → 175 | 30.82 | 4.53 | 0.00 | 4.53 | 3–25 MWe (up to 30–50) |
| 3b | High-temperature vapour-dominated (dry steam, Geysers type) | year 30 | 170 → 170 | 36.74 | 7.31 | 0.00 | 7.31 | 3–25 MWe (dry steam) |
| 4 | Modular ORC on surface waste heat (Climeon type) | steady | 110 → 110 | 8.41 | 0.56 | 0.00 | 0.56 | 355 kWe per module; 1–2 MWe stacked |
| 5 | Ground-source heat pump, 150 m U-tube borehole | year 10 | 5 → 5 | 0.01 | 0.00 | 0.00 | 0.00 | — |
| 6 | Shallow geothermal: 400 m coaxial borehole heat exchanger + heat pump | year 10 | 7 → 6 | 0.01 | 0.00 | 0.00 | 0.00 | — |
| 7 | Single-well deep coaxial retrofit in a hot field (GreenLoop type) | year 30 | 60 → 56 | 0.66 | 0.00 | 0.00 | 0.00 | 2–9 MWe claimed; 2 MWe retrofit |
| 8 | Multilateral deep closed loop (Eavor-Loop 1 type) | year 30 | 136 → 111 | 29.85 | 2.29 → 1.16 | 0.00 | 1.16 | 2.2 MWe (independent); 8.2 MWe (company, four loops) |
| 8 | the same loop at 30 kg/s, the flow that gives most electricity at year 30 | year 30 | 212 → 191 | 19.23 | 2.93 → 2.10 | 0.00 | 2.10 | as above |
| 9 | Doublet with horizontal closed-loop connections (Chevron-type, parameters assumed) | year 30 | 56 → 52 | 1.82 | 0.00 | 0.00 | 0.00 | undisclosed |
| 10 | Doublet EGS, hot dry rock (Soultz / UDDGP type) | year 30 | 182 → 185 | 15.97 | 1.89 | 0.28 | 1.61 | 3–5 MWe |
| 11 | Hybrid single-well EGS with surface heat exchanger (Eden Project type) | year 30 | 135 → 57 | 0.73 | 0.23 → 0.00 | 0.00 | 0.00 | small |
| 12 | Single-well cyclic fracture loop (Sage HeatRoot type) | day 365 | 84 → 72 | 0.29 | 0.01 → 0.00 | 0.00 | 0.00 | not stated |
| 13 | Multi-well cyclic EGS with subsurface thermal storage (Sage HeatCycle / Battery+ type) | day 365 | 84 → 72 | 5.88 | 0.18 → 0.10 | 0.00 | 0.10 | up to 50 MWe |
| 14 | Horizontal doublet EGS (Fervo Project Red type) | year 30 | 184 → 64 | 0.92 | 3.60 → 0.00 | 0.92 → 1.10 | 0.00 | 3.5 MWe pilot (63 l/s) |
| 14b | Horizontal EGS, one producer between two injectors (Fervo Cape Station type) | year 30 | 195 → 64 | 1.66 | 7.48 → 0.00 | 1.92 → 2.26 | 0.00 | >10 MWe per production well (2024 test); about 16 MWe gross (current design) |
| 15 | Fluid co-production from producing oil / gas wells (8 wells) | year 30 | 92 → 93 | 16.54 | 0.72 | 0.00 | 0.72 | 20–350 kWe |
| 16 | Closed loop in a single abandoned hydrocarbon well (CeraPhi / wells2watts type) | year 30 | 42 → 39 | 0.19 | 0.00 | 0.00 | 0.00 | — |
| 17 | Recovery of man-made heat after thermal EOR (SAGD chamber, HT ESP circulation) | year 5 | 112 → 58 | 0.67 | 0.24 → 0.00 | 0.00 | 0.00 | to be proven |
| 18 | Supercritical fluid produced from below a volcanic field (IDDP-2 type) | year 30 | 311 → 311 | 48.53 | 13.39 | 0.00 | 13.39 | >30 MWe expected |
| 19 | Superhot-rock EGS doublet (Newberry / AltaRock type) | year 30 | 264 → 264 | 35.47 | 11.93 → 5.67 | 0.00 → 0.46 | 5.21 | up to 50 MWe (??) |
| 20 | Ultradeep millimetre-wave drilled closed loop (Quaise type) | year 30 | 171 → 146 | 8.89 | 0.99 → 0.59 | 0.00 | 0.59 | 100 MWe |

The Project Red-type and Cape Station-type examples reproduce the tests in their first year and then exhaust their published
or assumed blocks (front-crossing time 4.5 and 10 years at the test flows): no net output after year 4 and year 10.

A report on the development-type study as of v2.0 (code changes, benchmark charts, a schematic and parameter card per type, time
series, depth profiles and discussion) is `aharrsim_dev_types_report.html`, built by `report/build_dt_report.py`; types 10, 11, 14 and
19 and the basis of the type 8 figure have changed since (see above).

Reading: below about 80 °C the water is heat only; 80–120 °C gives a binary plant 10–20 % of a small exergy; 120–200 °C is the
binary range; above 200 °C flash and steam plants take over. Anything that sweeps a convecting or fractured reservoir reaches megawatts
per well; anything that relies on conduction through a borehole wall stays at about a megawatt or less per well whatever the depth, so
the closed-loop capacities quoted by developers (2–9 MWe per retrofitted well, 100 MWe per 20 km plant) are well above what the
validated conduction model gives; the multilateral loop's nameplate of about 2 MWe per loop is matched only in 263 °C rock, and is two
to three times what conduction gives at Geretsried's 150–170 °C. The conventional figures are consistent with the model; the EGS
figures are test values, held only with enough flow through enough swept rock (next section).

## Uncertainty and open-system study (v2.1; figures as of v2.2)

`validation/study_cl.py` and `validation/study_open.py` ask three questions. `report/build_st_report.py` wrote the v2.1 report
(`report/aharrsim_uncertainty_report.html`), whose fractured-system figures are superseded by the ones below.

1. **Closed-loop uncertainty** (types 7, 8, 20; tornado, flow sweep, 160-run Latin-hypercube Monte Carlo of full simulations; type 8
   at 30 kg/s since v2.2). Gradient and rock conductivity carry it (each 10–30 % on heat, two to three times that on electricity);
   heat capacity is negligible; insulation is negligible at 3 km and first-order at 20 km. Conductivity falling with temperature
   (Vosteen & Schellschmidt, applied in 250 m layers at the undisturbed temperature) removes 11–22 % of the heat and 26–44 % of the
   electricity of the hot loops. Year-30 P10 / P50 / P90: 0.49 / 0.59 / 0.70 MWth (3 km), 1.22 / 1.61 / 2.04 MWe (multilateral loop
   at 30 kg/s; base 2.10), 0.23 / 0.33 / 0.46 MWe (20 km). Only a convecting formation (effective conductivity Nu × k) changes the scale.
2. **Open systems** (types 10, 14, 14b, 19, cyclic 12–13): fracture count 5–5,000, tortuosity 1–5, channelling, flow 10–400 kg/s,
   volume ×0.25–×64, in a screening model checked against the full injector → fractures → producer chain (18 designs at five times:
   within 5 %, median 1 %). Density and tortuosity saturate once the spacing is 10–20 m; the quoted capacities need flow through
   swept volume, held for 30 years: 3 MWe at 200 °C from 62 kg/s through 0.13 km³; the 3.5 MWe of the Project Red geometry from
   70 kg/s through 0.15 km³ (ten times the published volume); 10 MWe from a Cape Station-type producer from 195 kg/s through 0.43 km³;
   50 MWe at 375 °C from 190 kg/s through 0.43 km³. All four designs sit at x ≈ 0.8, y ≈ 0.15 on the swept-rock chart.
3. **Open systems with the closed-loop uncertainties** (tornado and 3,000-sample Monte Carlo): conductivity and heat capacity stop
   mattering (under 2 % for a well-fractured design); a gradient 20 % either way moves the four designs across 1.4–6.6, 1.5–7.0,
   2.7–15.9 and 24–119 MWe; ambient temperature and the flow-path unknowns (impedance, active fraction of fractures, swept fraction
   of the volume) follow. 30-year mean net P10 / P50 / P90 at 15 °C: 1.3 / 3.2 / 5.3 MWe (type 10, 100 fractures, 80 kg/s),
   1.8 / 3.4 / 5.2 (type 14, 16× volume, 80 kg/s), 1.3 / 5.5 / 10.8 (type 14b, 8× volume, 200 kg/s), 29 / 48 / 79 (type 19,
   8× volume, 200 kg/s). The medians lie 10–45 % below the base cases because channelling and bypass can only subtract.

```
python3 validation/study_cl.py oat|mc|opt        # Part 1 (about an hour on two cores; resumes; --fresh to restart); oat30 8, mc30 8
python3 validation/study_open.py screen|confirm|cyclic|part3|need
```

Multilateral-loop cases above 100 kg/s or with 6⅛ in laterals did not finish within the 240 s per-case limit and are omitted.

## Rock quality and lifetime output (v2.4)

`aharrsim.reservoir.UnevenFractures` is the parallel-fracture element without the assumption of equal flow: fractures carry the flows given by a
pattern of multipliers (or by a Lorenz coefficient, a corridor width and a seed) and exchange heat through the rock between them.
`validation/study_rock.py` runs three scenarios of rock quality on the four fractured systems; `validation/zoned_verify.py` verifies the element;
`validation/lifetime.py` gives the net electricity and heat of every example per kilometre of hole. Results and findings are in `CHANGELOG.md`.

```python
from aharrsim.reservoir import UnevenFractures
r = UnevenFractures(Q=0.11, N=300, L_flow=137.0, B=143.0, S=10.2, k_r=2.8, rhoc_r=2650 * 850.0, T_rock_C=200.0, T_inj_C=60.0,
                    rhoc_w=4160.0 * 900.0, lorenz=0.5, corridor=5, seed=0)
r.T_D([3.156e7, 3.156e8])      # dimensionless drop of the mixed outlet after 1 and 10 years
```

## Physics and numerics

| Part | Model |
|---|---|
| Fluid energy | Steady-flow energy equation in enthalpy for annulus and tubing, including gravity work (`m g dz`), kinetic energy and fluid heat storage; `T = T(p, h)` so compression/expansion, friction heating and flashing come out of the equation of state. Box scheme (2nd order) with exponential weighting where a cell's transfer number exceeds 1; automatic axial refinement where the tubing is poorly insulated (short counter-flow coupling length). |
| Rock + completion | Axisymmetric r–z finite volumes (casing steel, cement and layered rock as separate materials), radial and vertical conduction, far-field at the undisturbed geotherm. An equilibrium source makes the undisturbed geotherm an exact discrete steady state. |
| Time | Variable-step BDF2 (2nd order), geometric step growth, piecewise-constant operating schedule (flow, inlet temperature, shut-ins). |
| Coupling | Fluid enthalpies and all solid temperatures solved in one sparse system each iteration; `T(h)` linearised (Newton) and iterated with the pressure profile. |
| Heat transfer coefficients | Gnielinski (pipe), Gnielinski (2009) concentric annulus (inner/outer wall), laminar/transition blending; Chen factor in two-phase. |
| Flow direction | Annulus down / tubing up (default) or tubing down / annulus up; the solver works with a feed stream and a riser. |
| Hydraulics | Pressure marching each iteration: gravity, Churchill friction (all regimes), acceleration; two-phase homogeneous (HEM) or drift-flux (Zuber–Findlay with Hasan & Kabir 2010 parameters). |
| Properties | IAPWS-95 (CoolProp HEOS); bicubic p–h tables for speed (≤0.011 K in liquid, ≤0.05 K in steam); exact IAPWS-2008/2011 viscosity and conductivity; saturation table (clustered up to the critical point) for two-phase states. Valid to 1000 MPa, so 12 km bottom-hole pressures are in range. |

Not modelled: groundwater advection (pure conduction), temperature-dependent rock conductivity within a run
(the uncertainty study applies it layer by layer at the undisturbed temperature), natural convection in a shut-in annulus,
non-condensable gases, and in the fracture element unequal flow between fractures, heat from the rock outside the stimulated
block and thermo-elastic changes of aperture. Critical (choked) two-phase flow in the riser
is flagged with a warning but not limited. When the hot riser drives the loop harder than needed (thermosiphon),
the inlet is held at a minimum pressure and a wellhead choke takes the excess; the choke pressure drop is reported.

## Install and run

```bash
pip install numpy scipy CoolProp pyyaml matplotlib      # Python >= 3.9 (pygfunction, mpmath for verify_ext.py and screening.py)
python -m aharrsim run examples/deep_8km.yaml --plot --progress
```

Outputs (in `examples/deep_8km_results/`): `timeseries.csv`, `profiles.csv`, `summary.json`,
`summary.png`. All inputs are documented in `examples/deep_8km.yaml`.

From Python:

```python
from aharrsim import load_case, Simulator
sim = Simulator(load_case("examples/deep_8km.yaml"))
res = sim.run()
print(res.ts["T_out"][-1] - 273.15, res.ts["Q_wellhead"][-1] / 1e6)   # degC, MWth
```

## Verification and validation

Scripts in `validation/` regenerate every number in the report (`aharrsim_report.html`); results are
written to `validation/results/*.json`.

```bash
cd validation
python verify.py          # analytical checks, convergence, properties, reversed flow   (~25 min)
python field.py           # Qingdao, Weggis, Xi'an, Songyuan + conduction bounds
python twophase.py        # flowing geothermal wells (two-phase riser hydraulics)
git clone https://github.com/pnnl/GeoCLUSTER external/GeoCLUSTER   # public benchmark code + database
python benchmarks.py      # GeoCLUSTER database, SBT, published WellTH/SBT/TETRAD results
python deep_study.py      # 4-12 km base cases and sensitivities, 20-60 C/km, best case (~2.5 h on 2 cores)
python field_sizing.py    # wells needed for 5 GWth / 5 GWe from the deep-study results
python direction_study.py # tubing-down vs annulus-down over 17 casing x tubing combinations (~30 min)
python direction_analysis.py
python ../report/build_report.py
```

| Check | Result |
|---|---|
| Exact steady reverse-circulation solution | second-order convergence; 4×10⁻⁴ K max error at 50 m cells |
| Carslaw & Jaeger constant-temperature cylinder, 1 h–30 yr | within 0.1 % of the exact heat flow |
| Energy balance of the coupled solution | closes to round-off (~10⁻¹³) |
| Reversed flow: exact forward-circulation solution; mirror symmetry of the two directions | second-order convergence; directions agree to 10⁻¹² K with constant properties |
| Grid / time-step refinement, real 8 km design | all variants within 0.1 K from 30 days on |
| IAPWS-95 check values; property tables | 3×10⁻⁹; ≤0.011 K liquid, ≤0.05 K steam |
| Qingdao 2605 m field test (138 d) | outlet +0.2 K; mean heat +6 % |
| Weggis 2302 m (2 years, rock data assumed) | annual heat +9 % / +4 %, inside input uncertainty |
| GeoCLUSTER database (Sandia), 4–5 km, 126 points | RMS 0.10 K, max 0.38 K; outlet pressure within 0.2 MPa |
| SBT (GeoCLUSTER code), vertical 4–10 km | max 0.57 K |
| WellTH (Kubacka & Karayiannis 2024), 1–5 km | max 1.9 K |
| SBT (Beckers et al. 2022), 2–4 km | 20-yr mean heat within 1.1 % |
| TETRAD (Nalla et al. 2004), 5.6 km | −19 to +13 K (TETRAD shows almost no insulation effect) |
| Flowing geothermal wells (two-phase) | within 2 % for quality ≤ 0.35; 7–31 % at quality ≈ 0.6, < 20 bar |
| Xi'an 2000 m, Songyuan 2044 m | not usable: reported heat exceeds the conduction bound for the published rock |

Field tests only reach 2.6 km; the 4–12 km range rests on the analytical checks (no depth limit) and the
code-to-code comparisons.
