# Reproducing the tables and figures of the manuscript

"Thermal performance of geothermal development concepts and assessment of reported capacities" (Wrobel, McGeer and Sellar), release 2.5.0.

Every number in the manuscript is read from a file in `validation/results/`, and every such file is written by a script in `validation/`.
The result files are part of the archive, so the tables and figures can be regenerated without rerunning the simulations.

## Environment

Python 3.9 or later. `pip install -r requirements.txt` installs numpy, scipy, CoolProp, pyyaml and matplotlib, and for the validation scripts
h5py, pytest, pygfunction and mpmath. `pip install -e .` installs the package itself. `pytest -q` runs the seven tests (about 15 s).
The GEOPHIRES-X runs of Table 5 need GEOPHIRES-X 3.20.1 installed from source; their inputs and outputs are stored in `results/GX_runs.json`,
so the table can be rebuilt without it. Run the commands below from `validation/`. Times are for two cores.

The calculations added in release 2.5.0 were run with Python 3.13.18, numpy 2.5.3, scipy 1.18.1, CoolProp 8.0.0, matplotlib 3.11.2 and mpmath 1.4.1. A rerun of
all 22 benchmark cases in that environment reproduced the archived summary of release 2.4.1 without any difference, before the two changes of 2.5.0 to the
benchmark were made (the productivity index of type 3b and the circulation load of type 8; see `CHANGELOG.md`).

## Tables

| Table | Command | Result file |
|---|---|---|
| 1 Verification | `python verify.py`, `python verify_ext.py`, `python twophase.py`, `python field.py`, `python benchmarks.py`, `python zoned_verify.py` | `V1…V11_*.json`, `T1_twophase.json`, `F_*.json`, `B1…B4_*.json`, `ZF_verify.json` |
| 2 Examples | the case definitions in `dev_types.py` (one function per type) | – |
| 3 Benchmark | `python dev_types.py` (22 cases, about 15 min; `python dev_types.py 8` runs one type), then `python dev_types_table.py` (`md` prints the table) | `DT_types.json`, `DT_summary.json`; the row of type 8 at 30 kg/s: `python screening.py loopflow` → `SC_loop_flow.json` |
| 4 Lifetime output per kilometre | `python lifetime.py` | `LT_lifetime.json` |
| 5 GEOPHIRES-X comparison | `python geophires_compare.py run` (needs GEOPHIRES-X), `python geophires_compare.py table` | `GX_runs.json`, `GX_compare.json` |
| 6 Closed loops, sensitivities and percentiles | `python study_cl.py oat`, `mc`, `oat30 8`, `mc30 8`; `python study_review.py boot` | `ST1_oat.json`, `ST1_mc.json`, `ST1_oat_8_m30.json`, `ST1_mc_8_m30.json`, `RV6_percentiles.json` |
| 7 Fractured designs, net and gross targets | `python study_open.py screen`, `python study_open.py need`; `python study_review.py gross`; `python study_review.py share` | `ST2_screen.json` (`need_fine`), `RV1_gross_need.json`, `RV2_impedance_share.json` |
| 8 Fractured designs under uncertainty | `python study_open.py part3` (seed 7, 3000 draws per design and ambient temperature) | `ST3_open_uncertainty.json` |
| 9 Flow-heterogeneity scenarios | `python study_rock.py examples`, `needs`, `volume`, `scale`, `chain`; `python study_review.py seeds` | `RQ_examples.json`, `RQ_needs.json`, `RQ_volume.json`, `RQ_scale.json`, `RQ_chain.json`, `RV3_arrangements.json` |
| 10, 11 Screening charts | `python screening.py` | `SC_screening.json` |
| 12 Field records | `python field_cases.py`; `python study_rock.py red`, `red_curve`; `python study_review.py red` | `FE_field_cases.json`, `FE_project_red_digitised.json` (the digitised record), `RQ_red.json`, `RV4_red_flow_basis.json` |
| C1 Basis of the capacity figures | from the cited sources; no calculation | – |

## Calculations quoted in the text

| Where | Command | Result file |
|---|---|---|
| Section 5.3, circulation load of the sealed loops | `python dev_types.py 5 6 7 8 9 11 16 20` (`loop_pump`, `well_pump_kW`; `loop_circulation` for types 5 and 6) | `DT_types.json` |
| Section 5.3, the two multilateral-loop cases without a result | `python study_cl_timeouts.py limit 2400`, `python study_cl_timeouts.py pressure`, then `loads` | `RV7_type8_timeouts.json`, `RV7_type8_inlet_pressure.json` |
| Section 7.2, loop balance over the flows of Table 6 | `python study_cl_timeouts.py flows`; the 40 kg/s variant of type 20 in `dev_types.py` | `RV7_type8_flows.json`, `DT_types.json` |
| Section 7.1, full chain against the screening form | `python study_open.py confirm` | `ST2_confirm.json` |
| Section 8.1 and Appendix D, conduction along the flow path and from the surrounding rock | `python study_review.py paths`, `python study_review.py paths_grid` | `RV8_conduction_paths.json`, `RV8_conduction_grid.json` |
| Appendix B, feed-pressure floor and deliverability of the unaided wells | `python dev_types.py` (`flags` of each case) | `DT_types.json`, `DT_summary.json` |

`study_review.py outside` (`RV5_outside_conduction.json`) is the earlier two-dimensional form of the surrounding-rock calculation, in the dense-fracture
limit; it is kept for reference and superseded by `paths`.

## Figures

`AHARRSIM_ROOT=. python paper/make_figs.py` draws all thirteen figures into `paper/figs/` from the result files (`python paper/make_figs.py 4 8` draws
Figs 4 and 8 only). `AUDIT=1` also checks the lettering of each figure. The figures read: Fig. 1 `V9*`, `B4`, `V10`, `V11`, `ST2_confirm`; Figs 2–4 `DT_types`,
`DT_summary`, `SC_loop_flow`; Fig. 5 `LT_lifetime`; Fig. 6 `GX_compare`; Fig. 7 `ST1_*`; Fig. 8 `ST2_screen`, `RV1_gross_need`; Fig. 9
`ST3_open_uncertainty`; Fig. 10 `RQ_*`; Figs 11 and 12 `SC_screening`, `FE_field_cases`, `SC_loop_flow`; Fig. 13 `FE_field_cases`, `RQ_red`.

## Random numbers

`study_cl.py` (Latin hypercube): seed 20261113. `study_open.py part3` (simple random sampling): seed 7. `study_rock.py` and `study_review.py seeds`
(arrangements of the unequal flows along the well): seeds 0 to 4 and 0 to 19. `study_review.py boot`: seed 1, 4000 resamples.
