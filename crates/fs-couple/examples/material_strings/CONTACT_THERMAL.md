# Finite contact heat in the material-string renderer

`material_strings PERFORMANCE NEW_OUTPUT.wav --contact-thermal PROFILE`
connects the original finite contact/enthalpy solve to the existing scheduled
pressure and PCM writer. Add the usual `--body-profile`, `--stop-profile`,
observer options and `--bow-trace` as needed. A selected shell pressure observer
still receives the same accepted mechanical step. The Instrument Lab launcher
`apps/instrument-lab/tools/render_material_strings.py` accepts the same
`--contact-thermal PROFILE` option, stages that profile byte for byte and
forwards it to the native renderer. Existing `--pack` substitutions and
analysis calibration stay on their original path.

Plastic shear work heats its named film cell. An explicitly routed series
drive deposits its actual dashpot loss into its named cell. Conductive links
exchange heat using the two original materials' conductivities. At every
nonlinear endpoint the film re-queries its original card for its temperature,
then the existing simultaneous contact solve uses that elasticity, yield
stress and plastic viscosity. Temperature-dependent modulus also retains the
original reversible heat/entropy account. The resulting traction acts on
the actual string/body/cavity state and reaches the normal pressure output.

The profile supplies geometry, source temperature intervals, heat destinations
and numerical controls. It supplies **no substitute heat capacity, density,
conductivity, friction coefficient, heat split or temperature correction to
the audio**. Every cell names a material pack already loaded by the same
assembly; its cached material card is used throughout.

## Required source data

Each referenced pack must supply `density` [kg/m³],
`specific_heat_capacity` [J/(kg·K)] and `thermal_conductivity` [W/(m·K)] over the
declared temperature interval at the performance's original `omega` query.
Existing source selection and dimensional checks apply. Missing or ambiguous
claims refuse. Finite heat capacity is integrated by the existing
`ResolvedHeatCapacityCurve` owner; the initial temperature is retained as an
exact caloric node with specific enthalpy zero. Only the enthalpy datum is
shifted. Original source charts and identities are retained and logged.

A film also retains its mechanically selected `shear_modulus`,
`shear_yield_stress` and `plastic_shear_viscosity` claims. The thermal option
adds density to the original film source requirements during assembly.
The film cell must use exactly that card, its initial temperature and its
physical volume `area × thickness`. A finite ribbon therefore needs one
separate film cell for **each physical contact port**, using that port's
distributed area.

Density is constant over the admitted interval. The existing thermal film
owner admits constant or affine shear modulus; it rejects nonlinear modulus
paths and derivative jumps because their additional thermodynamics are not
part of this model. Constitutive source applicability is checked on every
accepted endpoint; leaving the caloric/material domain stops the render.

## Profile format

One record per line, in this order. All quantities use SI units. Comments,
extra fields and trailing records are refused.

```text
frankensim-contact-thermal-v1
boundary adiabatic
solver TRIALS ITERATIONS DIFFERENCE_J CELL_TOL_J THERMAL_TOL_J TOTAL_TOL_J NEWTON_ABS NEWTON_REL KRYLOV_RESTART CALORIC_ERROR_J_KG
cells N
cell ID PACK T_LO T_INITIAL T_HI X0 Y0 Z0 X1 Y1 Z1 VOLUME SURFACE_AREA
films N
film PORT CELL_ID
series N
series-heat PORT CELL_ID
links N
link ID CELL_A CELL_B AREA DISTANCE_A DISTANCE_B
```

Repeat each `cell`, `film`, `series-heat` or `link` record the stated number
of times. IDs are unsigned integers; ports use the original bank's zero-based
physical contact order. `PACK` is the exact relative path used in the
performance, not a separate source lookup.

| Record | Meaning |
| --- | --- |
| `solver` | Maximum mechanical trials per sample, maximum nonlinear iterations, finite-difference perturbation in **cell energy** [J], per-cell/thermal/combined energy tolerances [J], Newton absolute/relative tolerances, Krylov restart size, and source caloric integration tolerance [J/kg]. These are numerical resolution controls, not material parameters. |
| `cell` | Original source pack; lower/initial/upper temperature [K] with the initial temperature strictly inside the interval; explicit cell endpoints [m], volume [m³] and surface area [m²]. Mass comes from source density times this volume. |
| `film` | One retained material shear port and its unique finite cell. Bind every existing film exactly once. |
| `series-heat` | One existing series-compliant drive and its receiving cell. Other drive losses retain their original external-loss account. A drive cannot be routed twice. |
| `link` | Two distinct existing cells, physical interface area [m²], and positive distances from their temperature locations to that interface [m]. Conductance is `A / (d_a/k_a + d_b/k_b)`, re-evaluated from both original source cards at every endpoint. |

Limits are 64 KiB, 256 cells, 1,024 links, 64 film bindings and 64 series heat
routes. A profile needs at least one physical heat route. Film and series
counts may individually be zero. Zero source conductivity gives zero
conductance; no replacement conductivity is introduced.

Neighbor-cell temperatures drive conduction; film temperatures also drive
the contact constitutive law. Thermal reassembly of strings, body and hair
requires their corresponding thermomechanical owners.

The explicit outer boundary is adiabatic. There is no inferred infinite
heat bath or convective coefficient. A neighboring finite cell must have its
own declared geometry and source material. Its stored heat remains in the
same conservation account. Profile cells must refer to already loaded
assembly packs; introducing an unrelated reservoir material requires an
explicit material source owner before this route can use it.

## Scope and output

This route retains fixed film geometry and tangent frames, positive
compression, solid material, the original sample clock, canonical bow/force
controls, and any selected body/cavity pressure observer. Independent finger or
fret stops can press, release and migrate along their original strings through
the existing `motion` and `station-motion` blocks in `--stop-profile`.
Their exact normal/attachment work enters the same mechanical and thermal
transaction; rejected trials advance neither geometry nor heat. Changing
callback size or enabling the trace does not switch to an isothermal step.

The bow/film footprint and thermal grid remain fixed. Moving those footprints,
thermal expansion and melt/recontact history still require corresponding
transport/work terms in the physical owner. Stop dissipation remains in its
original mechanical loss account; it is not relabeled as film heat.
Geometric stretching and retained longitudinal waves use their original
nonlinear mechanical/contact solve at every proposed thermal endpoint. Each
trial starts from accepted string, film and drive histories, and film plastic
heat, reversible affine-modulus exchange and explicitly routed drive loss enter
the same finite enthalpy network. String/body constitutive coefficients stay at
their declared source state; this contact-only model does not infer bulk string
heating or thermal expansion. Solid caloric state, affine film `G(T)` and the
existing residual-compliance exclusion remain required.

The render log records material/caloric identities, actual cell masses and
initial temperatures, then final temperatures and stored enthalpies.
`--bow-trace` records the same accepted traction and slip history as the
ordinary renderer. Tracing changes callback size only; it does not run a
second simulation. Existing WAV headroom and clipping reporting still apply.

`--decimate` can be combined with this profile and `--bow-trace`. Heat feedback
and trace rows retain every original mechanical step while the shared pressure
observer emits 48 kHz PCM. The declared mechanical horizon must contain complete
output intervals; admission occurs before either output file is created. The
observer keeps its causal filter delay and never restarts thermal history at
the audio clock.

The existing `shear-layer-bow.material-strings` pack is an **isothermal
manufactured integration fixture** and does not supply the additional thermal
source data. The native tests supply explicitly manufactured thermal claims
to verify this entire route. They are not measured rosin properties or a
calibrated violin/cello/bass preset. To simulate a specimen, provide supported
thermal and constitutive claims for that same specimen; changing those
claims then propagates through heat, traction, motion and pressure.

The library equivalent is
`MultiBowPerformance::with_contact_thermal_network(...)`. It borrows the
existing performance, owns its supplied finite thermal network, and implements
the same `PressureRenderer` trait. `into_thermal()` recovers accepted thermal
history for continuation; no mechanical or thermal state is reset.
