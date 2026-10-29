# Condensed-fuel pyrolysis: sourced inputs and runnable cases

The new `fs-thermal::enthalpy::pyrolysis` owner evolves actual virgin and char
inventories on the existing heat-conduction graph. The examples below exercise
heat-driven fuel release and compare the shared kinetic update to published
MaCFP model output. They do not establish experimental validation or a complete
PMMA/wood compartment-fire prediction.

## Retained primary sources

The three imported files are unchanged from
[MaCFP/matl-db, revision 9ae21bebd5dbc7478b8637fbc4df1c81291ec6ea](https://github.com/MaCFP/matl-db/tree/9ae21bebd5dbc7478b8637fbc4df1c81291ec6ea).

| Local file | Original path | Original Git blob |
| --- | --- | --- |
| `macfp_nist_pmma_2021.json` | `PMMA/Material_Properties/2021/MaCFP_PMMA_NIST.json` | `36bc54ad1a7e98bf67dc5c542562094efebb3d53` |
| `macfp_nist_tga_fds.csv` | `Verification/Model_predictions/MaCFP_PMMA_NIST_dynamic_tga_10K_FDS.csv` | `5a3f215a1b9da88525b84d7b50d0d5d7f31e9844` |
| `MACFP_LICENSE.txt` | `LICENSE` | `ce0d732f78242a2a802888a6cd652f87d9e2268d` |

Copyright (c) 2019 MaCFP; the supplied MIT license is retained. Cite the material
database as Batiot et al., *Measurement and Computation of Fire Phenomena
(MaCFP) Condensed Phase Material Database*, DOI
[10.18434/mds2-2586](https://doi.org/10.18434/mds2-2586), with the revision above.

The NIST property submission identifies first-order kinetics calibrated from
TGA: A=2.85e13 s⁻¹ and E=191000 J/mol, with zero solid yield. Its scalar thermal
inputs are Cp=2200 J/(kg K), density=1100 kg/m³, conductivity=0.20 W/(m K),
heat of pyrolysis=870000 J/kg and emissivity=0.9. The file distinguishes
calibration, literature values and assumptions. This is one published property
set, not an assertion that it is the most accurate PMMA model.

The governing reaction/storage equations are separately sourced to the
[FDS solid-phase technical chapter](https://github.com/firemodels/fds/blob/a26081f3b24915552fab7ebf078185af817d9043/Manuals/FDS_Technical_Reference_Guide/Solid_Chapter.tex),
revision `a26081f3b24915552fab7ebf078185af817d9043`, blob
`2a2a51c948ac036b8db584de25fdc25ef1084175`. This is an independent Rust
implementation using FrankenSim's existing thermal/material operators.

## Heat-driven one-dimensional slab

From the repository root:

```bash
cargo run -p fs-thermal --example heated_pyrolysis_slab -- \
  noncharring 60 0.01 48 50000
cargo run -p fs-thermal --example heated_pyrolysis_slab -- \
  authored-char 60 0.01 48 50000
```

Arguments are material scenario, end time [s], maximum boundary-update interval
[s], through-thickness cells, incident irradiance [W/m²], and an optional sixth
geometry selector (`fixed` or `receding`). Defaults retain `fixed` and match the
first command. CSV goes to stdout; configuration and conservation checks go to
stderr. The example partitions the end time directly, avoiding a spurious tiny
final interval from repeated clock addition.

The authored slab is 0.01 m² by 6 mm and initially 300.15 K, containing 66 g.
Its back and sides are insulated. The front absorbs 90% of the specified
irradiance and exchanges thermal radiation and h=10 W/(m² K) convection with
a 300.15 K environment. A bounded surface balance retains the front half-cell
conduction resistance; exactly its net watts enter the original thermal graph.
The net face power is frozen over each requested boundary interval. Interior
conduction uses the original stability bound; temperature and local reaction
extent are solved together. No independent heat-of-pyrolysis debit is added.

The non-charring case uses the scalar NIST inputs above. The second case is an
**authored char sensitivity**, changing residue yield to 20% and char heat
capacity to 1200 J/(kg K). It is not a wood material fit. Conductivity remains
explicitly fixed in both cases, with no inferred porosity/char mixing rule.

To retain the supplied constant reaction enthalpy, this demonstrator declares
an effective outgoing enthalpy:

```text
h_out(T) = [h_virgin(T) + 870000 - y_char h_char(T)] / (1 - y_char)
```

This is an effective solid-model energy convention, not an MMA NASA caloric
model. A real gas coupling must instead supply the original gas enthalpy and
compatible solid formation references; `fs-airflow::compartment::pyrolysis`
provides that source-bound adapter. The demo's scalar value must not be silently
reinterpreted as a gas with a different caloric law.

### Actual 60-second results

| Scenario | Cells | Maximum dt (s) | Released volatile (g) | Retained char (g) | Front face (K) |
| --- | ---: | ---: | ---: | ---: | ---: |
| Non-charring | 48 | 0.01 | 1.167932203 | 0 | 685.234667080 |
| Non-charring | 48 | 0.005 | 1.167725006 | 0 | 685.233974619 |
| Authored char | 48 | 0.01 | 0.932940286 | 0.233235071 | 685.178393363 |
| Authored char | 48 | 0.005 | 0.932774833 | 0.233193708 | 685.177705976 |

The explicit spatial refinement matters more than the timestep change in this
case. At dt=0.005 s, the non-charring emitted masses are:

| Cells | Cell width (mm) | Released volatile (g) | Front face (K) |
| ---: | ---: | ---: | ---: |
| 12 | 0.5000 | 0.915060651 | 690.413431200 |
| 24 | 0.2500 | 1.115985106 | 686.083189034 |
| 48 | 0.1250 | 1.167725006 | 685.233974619 |
| 96 | 0.0625 | 1.180131873 | 685.073350653 |

The 48→96-cell difference is about 1.05% of the finer emitted mass; it remains
a reported spatial sensitivity, not an error certificate. The 48-cell default
avoids the coarse 12-cell release bias while keeping the example small.
Across the ten retained slab runs, the largest absolute mass imbalance is
1.39e-17 kg and the largest energy imbalance is 1.81e-5 J. Numerical balance
does not remove the fixed-geometry, constitutive or resolution approximations.

## Receding noncharring slab

The optional geometry mode uses the same retained NIST/MaCFP scalar inputs and
effective outgoing enthalpy, with explicit constant-density shrinkage:

```bash
cargo run --release -p fs-thermal --example heated_pyrolysis_slab -- \
  noncharring 120 0.02 24 50000 receding
cargo run --release -p fs-thermal --example heated_pyrolysis_slab -- \
  noncharring 120 0.02 48 50000 receding
```

The original implicit heat/reaction operator advances each interval. Retained
mass determines the slab length at the supplied density; the back face stays
fixed. A monotone limited reconstruction conservatively remaps specific
enthalpy in retained-mass coordinate, then the original calorics recover
temperatures and original series links use the new lengths. The front surface
balance uses its actual current half-cell resistance. No small residual cell
mass is retained to keep the solver running.

Both runs completed 6,000 intervals with zero geometry retries:

| Quantity at 120 s | 24 cells | 48 cells |
|---|---:|---:|
| Emitted volatile mass (g) | 7.759298058 | 7.904790867 |
| Remaining thickness (mm) | 5.294609267 | 5.281382648 |
| Surface recession (mm) | 0.705390733 | 0.718617352 |
| Front face temperature (K) | 688.246669 | 686.013949 |
| Back cell temperature (K) | 419.152825 | 418.892591 |
| Final volatile rate (g/s) | 0.137994753 | 0.138924125 |
| Whole mass residual (kg) | +4.44e-16 | −2.78e-16 |
| Whole energy residual (J) | −1.16e-6 | −2.23e-6 |

Doubling the cells changes released mass and recession by **1.875% of the
coarse result** and front face temperature by **−2.233 K**. Front-cell centers
sample different physical locations across meshes, so the common face is the
appropriate boundary-temperature comparison. This is a two-mesh sensitivity
check, not a convergence-order or mesh-independence claim. The accumulated
absolute remap energy discrepancy is about 1.14e-8 J in each run; the separate
enthalpy-redistribution diagnostic measures profile redistribution, not lost
heat or a certified error bound.

The [retained receding results](receding_slab_example_results.json) contain
commands, source identity, actual endpoint accounts and 30-second checkpoints.
The original five-argument fixed invocation remains operational.
`authored-char` with `receding` refuses: this homogeneous shrinking owner does
not infer char density, pores, swelling or residue mechanics. Both the boundary
composition and geometry update remain first order; no experimental validation
is established by these numerical checks.

## Published FDS kinetic comparison

```bash
cargo run -p fs-thermal --example macfp_tga -- \
  data/reference/compartment_fire/pyrolysis/macfp_nist_tga_fds.csv 0.1
cargo run -p fs-thermal --example macfp_tga -- \
  data/reference/compartment_fire/pyrolysis/macfp_nist_tga_fds.csv 0.05
```

This uses the **same** `Arrhenius::isothermal_step` as the coupled solid.
The retained FDS material-temperature history is interpolated linearly, and
the backward-Euler rate is evaluated at each substep endpoint. This comparison
prescribes temperature; it does not solve or validate the solid heat equation.
The corresponding upstream FDS deck has different thermal inputs from the
scalar NIST JSON, so comparing the common kinetic law at its reported
temperature avoids falsely treating them as the same thermal simulation.

The CSV labels mass dimensionless but starts at 5. Both histories are
normalized to initial inventory 1. All 374 strictly positive source rows are
included in the main error statistic, including the final near-burnout row.
The secondary population of 373 rows has at least 1% of initial source mass.

| Maximum dt (s) | Steps | All-positive RMSE | All-positive max error | ≥1% RMSE | ≥1% max error |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0.1 | 22443 | 4.445455e-4 | 7.829503e-3 | 1.838596e-4 | 6.524859e-4 |
| 0.05 | 44884 | 4.376982e-4 | 7.815755e-3 | 1.682967e-4 | 5.954270e-4 |

The first source zero is at 2244.1 s. The model still retains mass fractions
0.007176276 and 0.007163719 respectively; no burnout threshold is applied to
force agreement. The later 210 source rows are excluded because their material
temperature freezes after depletion. This is published **model-to-model
verification**, not experimental validation or an independent property fit.

## Scope and retained results

The owner currently supports an endothermic, first-order single reaction with
distinct sensible virgin/char calorics. Gas leaves at its local production
temperature; its time-integrated enthalpy is retained. The base fixed-geometry
explicit path can still exhaust its conduction budget near depletion. The new
implicit and homogeneous receding owners address that numerical/geometric
limit without introducing residual mass. Changing pore structure, char
oxidation, multicomponent wood mechanisms, mechanical work and gas-volume
displacement remain separate capabilities requiring appropriate data.

The [machine-readable results](pyrolysis_example_results.json) retain the
actual configurations, ten-second trajectory checkpoints, conservation outputs and
kinetic comparisons. Verification used native Cargo on the original thermal
source and focused supporting crates because DSR/RCH were unavailable. The
original fixed-geometry increment passed 164 thermal tests; with implicit and
receding extensions the native suite passes 178, with one existing performance
test ignored. This is not a full workspace/DSR pass. Calibration data are not relabeled as
validation; Dalmarnock data are not used.
