# Change history

## 2.4.1 (16 November 2026)

Repository release only: the simulator and every result are as in 2.4. The 2.4.0 release on GitHub and Zenodo (https://doi.org/10.5281/zenodo.24270031) lacked part of the `validation` folder: `study_rock.py`, `zoned_verify.py`, the 2.4 version of `study_open.py`, `results/RQ_volume.json`, `results/ZF_verify.json` and `results/superseded/`. They are in this release.

## 2.4 (November 2026)

**Rock quality and a lifetime measure.** Added after a co-author's review of the manuscript. The benchmark tables stand as they are.

| Addition | What it does |
|---|---|
| `aharrsim.reservoir.UnevenFractures` | parallel fractures that carry unequal flows and exchange heat through the rock between them: the Gringarten solution without the assumption of equal flow. In Laplace space the fracture temperatures obey a tridiagonal system along the flow path; the outlet is a matrix exponential and the produced water the flow-weighted mean |
| `aharrsim.reservoir.ZonedFractures`, `zoned_dense_limit` | the limit of wide zones that do not exchange heat (flow-capacity and storage-capacity form), with its closed form for closely spaced fractures |
| `flow_pattern`, `lorenz_of`, `lorenz_curve`, `lorenz_sigma` | flows per fracture with a given Lorenz coefficient (lognormal quantiles, optional corridors, seeded order) and the coefficient of a measured allocation |
| `validation/zoned_verify.py` | verification: equal flows, closed form, heat balance, convergence, and an independent finite-volume solution (`results/ZF_verify.json`) |
| `validation/study_rock.py` | three scenarios of rock quality (Good, Moderate, Poor: Lorenz coefficients 0.2, 0.5 and 0.8) on the four fractured systems, the volume each needs, the scale of the unevenness, the Project Red record read as uneven flow, and two full-chain confirmations (`results/RQ_*.json`) |
| `validation/lifetime.py` | net electricity and heat over the simulated life of every example per kilometre of hole (`results/LT_lifetime.json`) |
| `study_open.design(..., uneven=...)`, `confirm_one` | the screening model and the full chain with uneven flow |

Findings:

* **The element is verified three ways.** Equal flows return the parallel-fracture solution to 3 × 10⁻¹⁴; alternating zones agree with a finite-volume solution within 0.007–0.020 in T_D (the
  discretisation error of the latter); and wide zones tend to the mix of separate blocks (0.33 at 20 m, 0.007 at 1300 m).
* **Designs sized for even flow have no margin.** The four designs that just hold their quoted capacities keep 95 % of their 30-year net output in Good rock, 42–55 % in Moderate and
  17–32 % in Poor. Holding the outlet for 30 years takes 1.1, 2.7–2.8 and 5.4–6.0 times the swept volume.
* **Scale matters as much as the coefficient.** Conduction evens out differences between fractures closer than the conduction length (34 m after 30 years). With the same six flows at a
  Lorenz coefficient of 0.5 the year-30 outlet of the Cape Station-type design is 153 °C between single fractures 10 m apart and 112 °C between corridors 160 m wide (even flow: 196 °C).
* **The Project Red record does not distinguish its readings.** Fewer active fractures (0.45), half the volume swept, and uneven flow among all the fractures (Lorenz coefficient 0.35)
  fit the 614 days within the digitising accuracy and diverge afterwards.
* **Lifetime net electricity per kilometre of hole** spans two orders of magnitude: 290–790 GWh/km for wells in convecting reservoirs, 11–42 GWh/km for the fractured examples at
  190–200 °C (76 and 156 for the horizontal-well designs sized to hold their quoted capacities), and 6–9 GWh/km for the sealed loops that generate at all.

Also: the reference list of the manuscript was checked against its sources; the Qingdao field test is from Bu, Ran and Zhang (2019) and the GeoCLUSTER database from Beckers et al. (2023).

Figures: `paper/make_figs.py` draws the thirteen figures at the journal's widths (188 and 140 mm) with no lettering below 7 points, and `AUDIT=1 python paper/make_figs.py` checks each one for lettering that overlaps other lettering, lies on or within a point of a line or marker, or is cut off. No number in any figure changed; type 12 in Fig. 3 is now drawn at its value (5 kWe) instead of at the axis floor.

## 2.3 (November 2026)

**A laboratory comparison and two elements for a single well.** No existing result changes; the benchmark tables stand as they are.

Ma et al. (2026, Geothermal Energy 14:9, doi:10.1186/s40517-026-00379-4) circulated water through a single-well EGS built in a 1.5 m sandbox of 1 mm glass beads:
injection through an upper screen, production through a lower one 90 or 120 cm away, 0.2–0.6 m³/h, pack at 50–70 °C. It is type 11 of the benchmark at about 1 : 150,
in the swept-rock regime with water and solid at one temperature.

| Addition | What it does |
|---|---|
| `aharrsim.porous.AxisymmetricBlock` | Darcy flow and heat transport in a porous or densely fractured cylinder round a well (radius and height), with screens, buoyancy, temperature-dependent viscosity, anisotropic permeability and closed boundaries; MUSCL in space, Heun in time; the flow is solved again as the block cools |
| `aharrsim.reservoir.SingleWellDipole` | two zones on one well in a large block, as a reservoir element: sharp front on every streamline of the source–sink pair (first arrival at 0.746 D³(ρc)_b/(Q(ρc)_w)), or the Gringarten parallel-fracture response along every streamline when the flowing fractures are a given distance apart |
| `validation/segs_sandbox.py` | `checks`, `lab`, `calib`, `rec`, `field` (about 25 minutes in all; `curves` refreshes the dimensionless curves only); results in `validation/results/SG_*.json` |
| `report/build_sg_report.py` | builds `report/aharrsim_segs_report.html` |

Findings:

* **The experiment agrees with the swept-rock chart.** The reported plateaus (50, 30 and 25 min at 0.3, 0.5 and 0.6 m³/h) are one flushed share, x = 0.55 of the heat capacity of the
  bead pack between the screens: 41 MJ taken on the plateau at every flow. Flow sets how fast the heat is taken, not how much.
* **Model of the sandbox.** With the published properties the plateau of the nominal case ends at 21 min against 30 reported (60-minute mean 50.8 against 55.3 °C), because the
  cold water sinks to the production screen faster than measured. The permeability was estimated, not measured. Fitted to one number (the nominal mean), an isotropic permeability of
  1.8 × 10⁻¹⁰ m² (published: 8.92 × 10⁻¹⁰ m²) gives a plateau of 28 min and reproduces the other reported means within 1.5 K; the published mean with the horizontal
  permeability 2.8 times the vertical does the same within 1.0 K. Taken alone, the viscosity of the cold water lengthens the plateau (30 min) and its weight shortens it (17 min).
  The temperatures on the line of sensors across the box agree with the model; those reported for the vertical line do not, and are not consistent with the paper's own heat flows.
* **The recommended settings** (0.4 m³/h, 31 °C, 120 cm) have a plateau of 56 min in the model: the paper compares settings at 60 minutes, so its optimum is the setting that had only
  just reached breakthrough. A warmer injection gives up heat for no gain in life; spacing is the one real lever. The rule that carries is Qτ(ρc)_w / (V(ρc)_b) below the
  breakthrough share.
* **Field scale.** Zones 150 m apart at 10 kg/s reach first arrival after 4.3 years in a large fractured block and after 2.4 years in a block bounded like the sandbox;
  thirty years needs the zones 249 m apart in a large block (10 % of the rise lost) and 348 m in a bounded one (first arrival). The type 11 well with the dipole element gives 160 → 123 °C at the wellhead over 30 years,
  5.1 → 3.5 MWth and 0.39 → 0.16 MWe (benchmark element: 135 → 57 °C). The reading of the type as small is unchanged.
* **Injection above or below.** The sandbox and the type 11 example inject above the producer, where cold water sinks towards it. For the 150 m dipole the produced water keeps,
  after 30 years, 68 % of its rise with injection above and 76 % with injection below at a vertical permeability of 10⁻¹⁴ m²; 27 % and 97 % at 10⁻¹³ m². Uniform closed block; the direction is robust, the size is not established.

Checks: energy balance to 10⁻¹³; piston front at 0.997 of the sweep time; streamline table within 3 × 10⁻⁵; block model against the dipole element within 0.009 of the share lost;
sandbox grid within 0.04 K of the finest grid.

## 2.2 (November 2026)

**Three corrections.** They change the fractured-system results of v2.0 and v2.1; the coaxial and closed-loop results are not affected.

| Correction | Was | Now | Effect |
|---|---|---|---|
| Laplace inversion in `reservoir.ParallelFractures` | Gaver–Stehfest, 14 terms | de Hoog, Knight & Stokes (1982), `reservoir.dehoog` | Stehfest smears a sharp thermal front: wrong by up to 0.28 in T_D on the sharp-front test and by 0.12 (15 K) in the dense designs of the v2.1 study, where it brought the decline forward. Fractures 10–50 m apart at the example flows barely change. `T_D(t, method="stehfest")` reproduces the old results. |
| Heat capacity of water in the fracture element | a fixed volumetric value, equivalent to 4.64 kJ/kg/K | enthalpy mean between injection and rock temperature, `reservoir.cbar_w` (4.16 kJ/kg/K for 60 → 200 °C, 4.72 for 60 → 375 °C) | the fixed value was 12 % above the mean at 200 °C and made every fractured block cool a tenth too soon |
| Circulation power of fractured doublets | Q²·I/η plus a separate downhole pump | pressure balance round the loop with the thermosiphon credited: Δp_pump = I·Q + (p_bh − p_wh)_prod − (p_bh − p_wh)_inj + 0.3 MPa | type 10 pumps 0.48 → 0.28 MWe; the superhot doublet circulates unaided for more than 20 years |

Two corrections of basis: the type 14 example now uses the stimulated volume published for Project Red (975 × 110 × 143 m, 0.015 km³),
and the 8.2 MWe quoted for type 8 is the nameplate of the four-loop Geretsried project (about 2 MWe per loop), not of one loop.
The type 8 example is also run at 30 kg/s, the flow that gives most electricity (2.93 → 2.10 MWe against 2.29 → 1.16 MWe at 100 kg/s).

**Three additions.**

| Script | What it does | Result |
|---|---|---|
| `validation/geophires_compare.py` | runs the ten examples that GEOPHIRES-X 3.20.1 can represent (types 1, 2, 3, 8, 9, 10, 14, 14b, 15, 19) with the same inputs | liquid cases agree within 0.4–4.5 K and 1–8 % in output; thermal breakthrough (half the rise lost) at 4.42 against 4.37 years (type 14) and 10.12 against 9.99 years (type 14b). Differences that remain: Ramey has no expansion cooling, no flashing riser, c_p at a mean temperature, lifting of co-produced water |
| `validation/screening.py` | two dimensionless charts. Swept rock: outlet against x = t/t_s (flushed share of the block's heat capacity) and y = (S/2)/√(αt). Touched rock: line-source conductance U′(t) = 4πk/[ln(4αt/r_b²) − 0.5772] in an exchanger formula with NTU = U′L/(ṁc) | swept chart against the full injector → fractures → producer chains: 90 points, RMS 1.0 K, max 3.8 K. Touched chart against the simulated loops (laterals; coaxial wells to 3 km): 21 points, RMS 2.9 K, max 6.9 K, always high; 5–35 K high for the 20 km coaxial well |
| `validation/field_cases.py` | Project Red (614 days of production, Fervo Energy, April 2026), the 2024 Cape Station test and the first Geretsried loop, read with the same model | below |

Rules read off the charts. Fractures closer than 0.4√(αt) (14 m for 30 years) add nothing thermally; with that spacing the outlet
holds until x = 0.8, so a design life τ needs V ≥ ṁ c̄_w τ / [0.8 (ρc)_r]: one cubic kilometre per 450 kg/s for 30 years. A sealed
hole takes about 3 W per metre per kelvin (3.9 / 3.1 / 2.9 at 1 / 10 / 30 years, k = 2.8 W/m/K, 8½ in); NTU = 1 is 0.7 kg/s per
kilometre; the most electricity a kilometre of hole gives after 30 years is 8 / 20 / 38 / 66 kWe at 150 / 200 / 250 / 300 °C.
Per megawatt of gross electricity for 30 years:

| Rock temperature (°C) | Flow (kg/s) | Swept rock (km³) | or sealed hole (km) |
|---|---|---|---|
| 150 | 35 | 0.075 | 124 |
| 200 | 13 | 0.029 | 51 |
| 250 | 6.5 | 0.014 | 26 |
| 300 | 4.4 | 0.010 | 15 |
| 375 | 2.7 | 0.007 | – |

Field cases. Project Red (36 kg/s as operated, 175 °C, published volume): the block would be flushed in 7.2 years; by day 614
x = 0.23 and the heat produced is 20 % of the heat in place, as the developer reports. An ideally swept block would hold its
temperature to year 3; the record turns down at about day 500, which fixes the product of the active fraction of fractures and
the swept fraction of the volume at 0.43–0.53. Sustained net output (1.4 MW) is 40 % of the 3.5 MWe test figure. Cape Station:
the 2024 test is reproduced (195 °C; 8.0–10.6 MWe at 107 kg/s, depending on the correlation and the season); 16 MW gross at 221 °C
needs 146–155 kg/s and 0.12 km³ of swept rock per decade. Geretsried: with the flow fitted to the reported 8.5 MWth from four
lateral pairs (37 kg/s) the model gives 106 °C at the outlet and 0.41 MWe at 4.8 % conversion. The developer states 4.3 % on an
average duty of 11.8 MWth since start-up (about 0.5 MWe); no constant flow gives both heat figures, and the 71 kg/s whose mean is
11.8 MWth gives 90 °C, 3.9 % and 0.46 MWe. The
thermal nameplate (16 MWth per loop) is the conduction figure, the electrical nameplate (2.05 MWe per loop) is two to three times
what the chart gives at 150–170 °C. The Project Red temperature record is digitised from a published figure
(`validation/results/FE_project_red_digitised.json`); the Cape Station fracture height and
count are assumed, and the Geretsried flow is fitted to the reported heat.

```
python3 validation/geophires_compare.py run     # needs GEOPHIRES-X (git clone https://github.com/NREL/GEOPHIRES-X && pip install ./GEOPHIRES-X)
python3 validation/geophires_compare.py table
python3 validation/screening.py                 # charts and checks (needs mpmath); `loopflow` reruns the type 8 flow sweep
python3 validation/field_cases.py
python3 validation/study_open.py need           # smallest flow and swept volume that hold each quoted capacity
python3 validation/study_cl.py oat30 8          # type 8 uncertainty at 30 kg/s (also mc30 8)
```

The figures of the accompanying paper ("Swept rock or touched rock") are drawn by `paper/make_figs.py` from `validation/results/`.
The two HTML reports in `report/` on the development types (`aharrsim_dev_types_report.html`, v2.0) and on the uncertainty study
(`aharrsim_uncertainty_report.html`, v2.1) are kept as snapshots: their fractured-system figures predate the three corrections of v2.2 and
are superseded by `validation/results/`. The result files that the corrections replaced are kept in `validation/results/superseded/`.
