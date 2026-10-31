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

Version 2 also connects source-bound dry contacts to finite thermal strips on
both original participants. The actual bow/string and fiber crossings move
those footprints through material coordinates. Original cells retain their
heat when the contact leaves and revisits them. Normal release turns off both
friction heating and contact conduction through the existing mechanical owner.

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
owner admits constant or affine equilibrium shear modulus and Maxwell-arm
moduli. The ordinary `shear-memory` declaration selects each arm's modulus and
relaxation-time claims; no separate thermal-memory input is needed. At each
trial temperature the same retained sources determine the arm stiffnesses and
relaxation times. The arm's internal slip remains with its original mechanical
history when those properties change.

All retained elastic strains contribute to the film's mechanical entropy:
`s = -0.5 * k_inf'(T) * e² - sum(0.5 * k_j'(T) * r_j²)` [J/K], where `e`
is the equilibrium elastic slip and `r_j` is the relative slip stored in arm
`j`. Mechanical internal energy is `U = psi + T*s`. The finite heat solve
includes coefficient work and changes of `T*s`, together with the original
plastic and relaxation dissipation. Changing a branch stiffness therefore
cannot silently create or remove stored energy.

Nonlinear modulus paths and derivative jumps remain unsupported because their
additional thermodynamics are not part of this model. Constitutive source
applicability is checked on every accepted endpoint; leaving the
caloric/material domain stops the render. These are constitutive capabilities,
not a claim that an available material pack identifies measured rosin behavior.

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
the contact constitutive law. Version 4 can additionally route the original fiber material's hereditary
loss to its finite body. Thermal reassembly of strings, body and hair still
requires their corresponding thermomechanical owners.

## Moving dry-contact profiles

Use `frankensim-contact-thermal-v2` and insert the following block after the
last `cell`, before `films`. All other v1 records and meanings are retained.

```text
dry-contacts N
dry-contact PORT STRING_COORDINATE_ORIGIN fixed
footprint A_LO A_HI A_WIDTH B_LO B_HI B_WIDTH
surface N_A
strip CELL_ID LO HI WIDTH HALF_DEPTH
surface N_B
strip CELL_ID LO HI WIDTH HALF_DEPTH
```

Repeat each `strip` according to its surface count and repeat the complete
contact block according to `dry-contacts`. The first surface is the original
interface's string participant A; the second is participant B. For a material
fiber ribbon, replace the contact header with:

```text
dry-contact PORT STRING_COORDINATE_ORIGIN fiber FIBER_COORDINATE_ORIGIN
```

Origins map original string/fiber coordinates into the strips. Initial footprint
centres must match the actual contacts. Fiber receivers retain the original
fiber material and crossing. Bow travel, support acceleration and pressing use
the same sample/work accounts and midpoint contact maps as mechanics.

Each side's footprint area `(HI - LO) * WIDTH` must equal that physical port's
source-declared nominal contact area. A ribbon therefore needs a route for each
heated physical port, using its original distributed area. Strips are ordered,
contiguous and nonoverlapping; a cell cannot belong to both sides of one route.
The strip volume `(HI - LO) * WIDTH * 2 * HALF_DEPTH` must equal its original
cell volume. Every cell's material identity must match its ordered source
interface participant. Its density, heat capacity and conductivity come from
the same assembly material card. No contact-plane thermal mass is introduced.

The massless contact-plane temperature follows friction power and both surfaces'
finite `k*A/d` fluxes. Signed heat can cool a hotter body even without slip.
Cross-surface links refuse duplicated contact conduction. Geometry, source
identity, caloric mass and coverage are admitted before output and at every step.

Dry routes may share the network with explicit series-drive heat. Retained
shear-film bindings and moving dry routes cannot be combined by this profile;
set `films 0` for a dry profile. V1 film/series profiles keep their fixed-bow
restriction. Missing hair thermal properties are not supplied from rosin or
another material. Constant dry laws retain their original analytic interface
temperature path. A `coulomb` bow declaration may instead select source
coefficients that depend on `T`. Those original claims are pinned and evaluated
at the interface temperature solved jointly with traction and finite-cell heat.
The same slip work heats the cells, and their temperatures feed back into the
next coupled endpoint. No coefficients are inferred from a recording.

Temperature-dependent dry declarations require version 2 heat routes for every
affected physical port, including each distributed ribbon contact. Missing routes
refuse before output. The original `bow-interface` syntax and explicit T/p/v
applicability are retained. Pressure and incoming, average and endpoint slip
speed are checked against the real mechanics. Coefficient response to pressure,
speed or additional history is not supplied by this temperature-only coupling;
such response axes refuse rather than remaining frozen. Stribeck and hyperbolic
coefficient sources retain their constant-only restriction.

The interface temperature has no invented heat capacity. Its additional Newton
coordinate is scaled by the actual original half-cell conductances and sample
duration, solely for numerical units. Each accepted interface must close its
own flux/friction-work balance as well as the existing cell and total energy
accounts. Release disconnects heat exchange and preserves the source history;
failed trials publish neither new coefficients nor mechanical or thermal state.
This does not model a transported rosin layer or create missing rosin measurements.

The explicit outer boundary is adiabatic. There is no inferred infinite
heat bath or convective coefficient. A neighboring finite cell must have its
own declared geometry and source material. Its stored heat remains in the
same conservation account. Profile cells must refer to already loaded
assembly packs; introducing an unrelated reservoir material requires an
explicit material source owner before this route can use it.

## Version 3: a retained film against moving thermal material

`frankensim-contact-thermal-v3` retains the version 2 records and adds a
`film-surfaces` block immediately after `films`. The film remains attached to
its original string patch, with the same finite mass and elastic/plastic history.
The opposing material slides through its thermal footprint. Each receiving cell
retains its heat when contact leaves and revisits it.

For example, after declaring an existing `film 0 10`, a new block can read:

```text
film-surfaces 1
film-surface 0 +1
film-footprint 0.00025 0.00035 0.01 0.00005
surface 4
strip 20 0.0000 0.0001 0.01 0.00001
strip 21 0.0001 0.0002 0.01 0.00001
strip 22 0.0002 0.0003 0.01 0.00001
strip 23 0.0003 0.0004 0.01 0.00001
```

These are illustrative geometry records, not measured instrument dimensions.
All named cells must already exist. The film footprint fields are its initial
counterface-coordinate endpoints, width and film centre-to-interface distance
in metres. Its area and full thickness must match the original mechanical film.
Each receiver strip uses the same finite-volume geometry and sourced transport
requirements as a dry thermal strip. This route uses `dry-contacts 0` and the
original retained-film mechanical map.

With `bow_fibers`, a complete `film-surfaces` binding (one for every film)
explicitly attaches those films to their original string material patches.
`bow_fiber_motion` or `bow_fiber_program` may then accelerate, reverse and
translate the opposing hair through the existing source-derived fiber model.
Its changing material crossing changes the counterface velocity maps while
the string-side tangent, film mass and elastic/plastic/Maxwell history remain
attached. Every receiving strip must use the same original material card as
the physical fibers. Incomplete bindings refuse moving-fiber playback.

Canonical solved normal contacts and `bow_fiber_normal_program` use that same
moving hair geometry, including the original normal attachment work and support
inertia. The string-side film station must stay fixed: `bow_station_motion`
still requires film mass/history transport and is not admitted on this path.
Independent finger/fret motion remains supported through its original owner.
An invalid thermal endpoint or geometry update advances neither support
program, film history, finite cell state nor the accepted pressure sample.

The orientation is exactly `+1` or `-1`. At each nonlinear trial the receiver
interval translates by minus this orientation times the **actual candidate
elastic plus plastic relative displacement**. The commanded bow velocity does
not substitute for that displacement. Changing traction can therefore change
both the resulting motion and which receiving material is heated in that same
coupled endpoint. No separate stroke cursor advances during a rejected trial.

Heat crosses the original half-film and half-body Fourier resistances in series,
using their respective candidate temperatures and pinned conductivity claims.
Equal and opposite signed rates enter the existing cell energy equations;
plastic and viscoelastic relaxation heat remain deposited once by the original
layer owner. Fixed film-to-substrate links may remain, but a second link from the
film to a declared receiving strip refuses because it would duplicate conduction.
Finite coverage is checked before output and at every trial. Leaving the strip
refuses the step while preserving accepted motion, heat, trace and output state.

This is an attached-film, sliding-counterface approximation with adiabatic outer
boundaries. It does not transport rosin mass, renew or reset shear history, change
the string-side film attachment, model wear, or supply missing interface measurements.
Thermal cells and time steps still require physical resolution/refinement.

## Version 4: hereditary heat in the physical fiber bodies

`frankensim-contact-thermal-v4` retains all version 3 records and adds this
block after `series`, before `links`:

```text
material-heat-budget MAX_DIMENSION MAX_DENSE_BYTES MAX_MULTIPLY_TERMS
fiber-heat N
fiber-body STRIP LAYER CELL_ID
```

Repeat `fiber-body` for each selected physical material region. `STRIP` is
the original ribbon's zero-based across-width cell, and `LAYER` is the
zero-based radial material layer. A region cannot be selected twice, two
regions cannot share a receiving cell, and a fiber body cannot use a film's
cell. A body cell may also receive contact heat through the existing dry
contact or retained-film conduction route, provided that route admits the
same source and physical geometry.

The original `bow_fiber_memory` declaration supplies every modulus increment
and relaxation time. The heat profile supplies no additional damping value.
The geometry producer selects the complete original bending and axial memory
for the requested strip/layer. With `bow_fiber_strain`, axial heat comes from
the nonlinear full-strain material update in the same contact trial; bending
heat comes from the original complete held-force trajectory. Those sources
are retained separately and deposited once, rather than assigning a fraction
of total mechanical loss to a body.

The receiving volume must equal the original annular layer area times the
fiber length and number of fibers represented by the strip:
`pi * (r_outer^2 - r_inner^2) * length * fiber_count / strip_count`.
Density and this volume determine the finite mass. The cell must use the
same material card, source temperature and supported caloric law as its
mechanical region. One cell represents a uniform temperature over that
whole region; it cannot be replaced by a small contact footprint while
retaining the whole region's dissipated heat. Finer longitudinal temperature
variation requires its corresponding spatial material/heat discretization.

The dimension bound is at most 256 for the existing exact linear heat
observation. Byte and multiply-work limits are aggregate over the declared
receiving bodies and are divided conservatively among them before the
original operator admission. Source, geometry, mass, applicability and budget
checks run before the renderer opens its output; the physical owner repeats
material admission during each thermal trial. Missing source properties and
unsupported external acoustic trajectories still refuse.

A material-only profile may set `dry-contacts 0`, `films 0`,
`film-surfaces 0` and `series 0`. Combining fiber heat with admitted
contact routes and source-derived conductive links lets that heat change
contact temperature, traction and the same accepted pressure history.
The existing command, trace, decimation and PCM paths consume this owner.

This bulk hereditary route currently requires fixed density, geometry and
source mechanical coefficients throughout the admitted temperature interval.
It does not freeze a temperature-dependent relaxation time or modulus and
call that thermal feedback. Temperature-dependent body laws require their
own thermomechanical continuation. Film temperature-dependent laws continue
to use the existing simultaneous thermal/contact solve. Unselected physical
losses keep their original external-loss account.

## Scope and output

This route retains fixed film geometry and tangent frames, positive
compression, solid material, the original sample clock, canonical bow/force
controls, and any selected body/cavity pressure observer. Independent finger or
fret stops can press, release and migrate along their original strings through
the existing `motion` and `station-motion` blocks in `--stop-profile`.
Their exact normal/attachment work enters the same mechanical and thermal
transaction; rejected trials advance neither geometry nor heat. Changing
callback size or enabling the trace does not switch to an isothermal step.

The film patch and thermal grid remain fixed in their material coordinates;
the counterface region under the film follows the solved relative displacement.
Moving the film itself, thermal expansion and melt/recontact history still require corresponding
transport/work terms in the physical owner. Stop dissipation remains in its
original mechanical loss account; it is not relabeled as film heat.
Geometric stretching and retained longitudinal waves use their original
nonlinear mechanical/contact solve at every proposed thermal endpoint. Each
trial starts from accepted string, film and drive histories, and film plastic
heat, reversible affine-modulus exchange and explicitly routed drive loss enter
the same finite enthalpy network. String/body constitutive coefficients stay at
their declared source state. Version 4 explicitly routes source fiber hereditary
heat into its finite bodies; bulk string heating and thermal expansion still
require their corresponding physical owners. Solid caloric state, affine film `G(T)` and the
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
Moving dry surfaces select `with_dry_contact_thermal_network(...)` through the
same performance loop, finite thermal owner, observer, trace and PCM writer.
