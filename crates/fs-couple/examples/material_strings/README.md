# Material packs to a plucked or bowed instrument

This executable connects the existing normalized material-pack decoder, solid
property resolver, string/body/cavity assembly, optional torsion, canonical bow
schedule, physical load compiler and streamed PCM writer. It contains no modal
frequency presets, separate pluck oscillator, per-string audio mixer or second
physics integrator. All strings see one shared plate and optional cavity.

The supplied coefficients are **synthetic integration fixtures**, not measured
steel, gold, silver, wood, rosin or instrument data. The output is a reduced
physical-model demonstration, not a validated violin/cello/bass recording.

## Run without editing Rust

From the repository root, with its existing Rust toolchain and dependencies:

```bash
demo="$(mktemp -d)"
cp crates/fs-couple/examples/material_strings/*.tsv "$demo/"
cp crates/fs-couple/examples/material_strings/*.material-strings "$demo/"
cargo run -p xtask -- matdb-pack \
  --manifest "$demo/wire-manifest.tsv" --out "$demo/wire.fsmatpk"
cargo run -p xtask -- matdb-pack \
  --manifest "$demo/body-manifest.tsv" --out "$demo/body.fsmatpk"
cargo run --release -p fs-couple --example material_strings -- \
  "$demo/pluck.material-strings" "$demo/pluck.wav"
cargo run --release -p fs-couple --example material_strings -- \
  "$demo/pluck-and-bow.material-strings" "$demo/pluck-and-bow.wav"
```

Outputs must not exist. The application completes material, geometry, frequency,
control and preload admission before creating a WAV. A later physics/I/O failure
can leave an incomplete output prefix and is reported as failure; it is not a
complete performance. The original PCM encoder performs the declared Pa-to-PCM
conversion and reports clips. There is no automatic amplitude normalization.

Change a material TSV, compile a new pack to a new filename, and replace only its
pack path in the performance. Density, Young modulus and bending viscosity then
reach the same original string mass, dispersion, loss and attachment maps. With
torsion selected, shear modulus and shear viscosity also reach surface slip and
torque. Body changes alter the original plate reduction and its reciprocal load.
Use the same performance geometry/tension to compare materials without hidden
pitch normalization. The original packs, normalization records, rights and
physical-source identities remain available throughout rendering.

### Material-derived string/body joints

Every input version can replace the numeric bridge stiffness and damping at
the end of a `string` record with an ordered interface source:

```text
string wire.fsmatpk 0.5 0.0003 20 4 3 0.8 interface bridge.fsintpk 0.0001 normal_stiffness normal_viscosity
```

The fields after `interface` are pack path, physical area in m², stiffness-per-area
property and resistance-per-area property. The same `ResolvedLinearInterface`
owner used by material shell connections integrates these claims once into N/m
and N·s/m. The performance's `state` supplies the frozen T and angular frequency;
explicit source applicability must cover the retained physical modes. The
same admission also checks the coupled stiffness bound, including later gas
and rotational-joint changes. A cavity neck requires coverage down to zero
because its isolated inertia receives its restoring force through body links.
These conservative linear bounds do not qualify a nonlinear trajectory's full spectrum. The
original bilateral string/body solver supplies equal and opposite forces,
storage, loss and body pressure. Changing either interface property in the data
therefore changes reciprocal loading and the sound without editing a solver.

The ordered surfaces must match the **outermost string layer** and the actual
body material under its original physical force footprint. Uniform plates,
regional plates, shells and shell assemblies use this same check. A scalar
interface spanning different resolved materials refuses; it does not choose a
material from a name, a modal weight, a tied node or the interface's expectation.
Bulk substitution also requires a compatible interface source for that pair.
This is a fixed-state linear bilateral joint, not a model of a shaped bridge,
separation, slip, adhesive failure or changing contact area.

After building the wire and body packs above, the synthetic complete example is:

```bash
cargo run -p xtask -- matdb-pack \
  --manifest "$demo/bridge-manifest.tsv" --out "$demo/bridge.fsintpk"
cargo run --release -p fs-couple --example material_strings -- \
  "$demo/source-bridge.material-strings" "$demo/source-bridge.wav"
```

`bridge.tsv` carries the two source properties and exact ordered surfaces. The
example combines that joint with a gas-derived physical cavity neck. Its
coefficients are deliberately synthetic and carry no measured-instrument claim.

## Input records

The fixed-order header uses ASCII whitespace. File paths are relative to the
performance file's directory, contain no spaces, and cannot contain `..`.
Use LF line endings. After `schedule`, provide the **exact canonical bytes** from
`fs_scenario::gesture::GestureSchedule::to_canonical_bytes()`; the fixtures show
that existing tab-separated format. This executable does not invent another
ramp syntax or interpolator.

```text
frankensim-material-strings-v1
audio RATE_HZ SAMPLES FULL_SCALE_PA MAX_BLOCK
state SOLID_T_K SOURCE_OMEGA_RAD_S
gas T_K P_PA RELATIVE_HUMIDITY
body PACK LENGTH_M WIDTH_M THICKNESS_M MODES LISTENER_M
strings COUNT
string PACK FULL_SPAN_M RADIUS_M TENSION_N BENDING_MODES TWIST_MODES BRIDGE_STATION K_N_M D_N_S_M
... one string record per member ...
cavity none
friction coulomb MU_STATIC MU_KINETIC
forces retain PORT_COUNT EVENT_COUNT
port STRING_INDEX STATION centreline INITIAL_N
... one port record per load ...
force SAMPLE PORT_INDEX HELD_N
... one record per force assignment ...
schedule
... original canonical gesture bytes ...
```

`gas mixture T P mole COUNT SPECIES FRACTION ...` (or `mass`) delegates to the
existing named-gas parser used by wind instruments. Its species/state refusals
remain binding. That **same** frozen gas supplies all string drag, cavity
mechanics and compact body/neck pressure observations. Gas changes are a new
realization, not a live atmosphere mutation or an output equalizer.

For source-defined species and temperature-dependent heat capacity, replace that
row with `gas thermochemical T P mole COUNT SPECIES FRACTION ...` (or `mass`).
Supply each species' normalized species, NASA-9 and pure transport claim packs
through `--gas-component SPECIES.pack NASA9.pack TRANSPORT.pack`; repeat the
option for every component. These are the same source loader and constitutive
owners used by `music_render`. Declare `--gas-water none --gas-phase frozen`
explicitly, or identify the water species and choose the existing `liquid-water`
admission policy. Source identities, composition, validity and phase admission
stay attached to the resolved gas and are included in the render log. Unsupported
or incomplete source data refuse before rendering; a species name never supplies
missing thermal or transport properties. The atmosphere remains frozen during a
performance, with the same resolved state reaching all of the consumers above.

For a cavity whose neck loss follows the actual atmosphere, use:

```text
cavity physical-neck VOLUME_M3 RADIUS_M LENGTH_M INNER_END OUTER_END LISTENER_M
```

Each end is explicitly `flanged` or `unflanged`. The inner end opens into the
cavity; the outer end opens into the exterior. For example:

```text
cavity physical-neck 0.003 0.01 0.02 flanged unflanged 1
```

This selects the existing `CavityNeckLoss::RadiationAtResonance` owner. The
shared gas sets neck inertia and cavity compressibility, Bessel viscothermal
neck resistance, and compact radiation resistance at the declared outer end.
The two declared ends also supply their existing inertial end corrections,
counted once. Changing gas composition, temperature, pressure, or source-backed
thermal/transport data therefore changes the physical cavity load carried by
the same string/body mechanics. It does not require a new damping coefficient.
The selection works with the uniform plate, regional plate, curved shell and
assembled material-body profiles; their original signed volume maps remain
the common pressure coupling.

The resistance is evaluated at the **isolated Helmholtz frequency** and held
through the performance. This remains a uniform-pressure, circular compact-neck
approximation, not a broadband neck impedance, distributed cavity, f-hole model,
or exterior radiation solution. The original compact-radiation `ka <= 1`
admission applies. Pressure output retains the existing compact baffled,
zero-retardation body/neck observer; selecting an unflanged load does not add a
new spatial microphone/directivity model.

`physical-cavity.material-strings` uses this load with the original pluck-and-bow
gesture and synthetic string/body packs. After the pack commands above:

```bash
cargo run --release -p fs-couple --example material_strings -- \
  "$demo/physical-cavity.material-strings" "$demo/physical-cavity.wav"
```

`cavity neck VOLUME_M3 RADIUS_M LENGTH_M RESISTANCE_PA_S_M3 LISTENER_M` selects the
existing Helmholtz model with an explicitly specified neck resistance. It replaces
`cavity none`. Compression, neck motion, body motion and bridge reactions remain
in the same mechanical solve. This legacy research-fixture form keeps its
original two flanged ends and constant authored resistance; gas changes do not
recompute that supplied resistance. Existing files keep their original behavior.

`friction stribeck MU_STATIC MU_KINETIC CHARACTERISTIC_SPEED_M_S VISCOUS_PER_SPEED`
is also accepted. These coefficients retain **caller authority**. Bulk material
data do not certify a rosin/string friction law. The render log labels these
numeric inputs as a research override.

Alternatively, replace the friction row with source-bound ordered interfaces:

```text
friction interfaces 2
bow-interface first string-a-bow.fsintpk 293.15 0.0001 2 0 coulomb mu_static mu_kinetic
bow-interface second string-b-bow.fsintpk 293.15 0.0001 2 0 stribeck mu_static mu_kinetic stribeck_speed viscous_per_speed
```

The fields after the track and normalized interface-pack path are **interface
temperature [K], nominal contact area [m²], query contact pressure [Pa], and query
slip speed [m/s]**. These are followed by the declared law and its property names.
Coulomb requires static and kinetic dimensionless coefficients; Stribeck also
requires characteristic speed [m/s] and a viscous coefficient per speed [s/m].
Each bow track needs exactly one declaration, in any order. Different contacts
may use different packs and laws; repeated pack paths are decoded once.

The original `NormalizedInterfacePack` must declare `medium = dry`. Its ordered
surface A must be the selected string's **actual outer material state**, including
the outermost sleeve for a bonded string. Surface B identifies the prescribed
bow surface; texture, third body and named history remain part of the original
source identity. Changing a string's outer material while retaining an interface
for the old surface refuses instead of reusing its friction coefficients.

The selected claims must provide constant coefficients with explicit finite
bare-SI validity on `T`, `p` and `v`. The query chooses a supported immutable
source; its pressure and speed are not imposed on the physical solver. The
common applicability intersection is checked against the actual normal load
divided by nominal area and the actual relative surface velocity, including
torsional motion, at the original solver's initial, step-average and endpoint
checks. Missing, conflicting, nonconstant or out-of-domain sources refuse. A
released bow does not require positive-pressure friction data at zero load.

Interface temperature, nominal area and source history are frozen inputs. This
binding does not infer a Hertz patch, evolving hair compliance, rosin heating or
wear. The original packs and property-use receipts survive with the performance;
loading a source does not promote it to experimentally validated instrument data.

An optional record immediately after the complete friction declaration selects
a finite-width rigid ribbon instead of the default point bow:

```text
bow_ribbon WIDTH_M ODD_CELLS
```

The width is an interval along the actual material string, centred at the
original gesture station; the complete interval must remain strictly inside
that string. The odd cell count is an explicit spatial resolution from 1 to 63.
This first file profile accepts exactly one original bow track: its speed acts
at every cell, while its **total normal load** is divided evenly. Each cell uses
its own material-derived transverse/torsional surface map and solves its own
slip and traction in the existing simultaneous string/body/cavity mechanics.
The original schedule, including ramps and release, survives unchanged.

With source-bound friction, the declared nominal area belongs to the **whole
ribbon**. Area and normal load are divided by the same cell count, retaining
the same nominal pressure and original ordered interface card at each cell.
Runtime applicability uses each cell's actual slip. Numerical setup, control
storage and compilation work are charged for every cell under the same global
limits. Independent finger/obstacle stops remain available; a `bow-contact`
binding that derives normal load from a stop refuses this prescribed uniform
load profile. This rigid, uniformly loaded limit does not infer compliant bow
hair or evolving rosin. Refine spatial cells, structural modes and time together.

The included synthetic `source-bow.material-strings` uses the same string, body,
gas and gesture as `pluck-and-bow.material-strings`, with its friction carried by
an ordered interface pack and `bow_ribbon 0.01 5`. Removing that optional record
restores the point contact. After the initial wire/body build above:

```bash
cargo run -p xtask -- matdb-pack \
  --manifest "$demo/bow-manifest.tsv" --out "$demo/bow.fsintpk"
cargo run --release -p fs-couple --example material_strings -- \
  "$demo/source-bow.material-strings" "$demo/source-bow.wav"
```

Changing `bow.tsv` and compiling a new interface pack changes excitation through
the retained physical contact; no string/body solver or per-note sound parameter
needs editing. This fixture is deliberately synthetic and remains unqualified
as a measured bow/string model.

`forces preload ...` settles the entire initially unexcited instrument before
the window, with bows released. A force event at sample zero can release this
initial load. Its initial stored energy is reported separately, not relabeled as
work during the performance. `retain` starts from rest and applies a load-on
transient. `forces retain 0 0` supplies no external load program. A `surface` port
adds the original bow-facing torque when twist is present; `centreline` applies
only transverse force. Same-sample force and bow changes precede one joint step.
Distinct force ports add; repeated assignments to one port use the last value.

Each canonical bow track identifies its actual material string and fixed station.
One track may stay released for a purely plucked instrument. Multiple tracks use
the shared simultaneous-contact solver, not independent voices. Speeds and normal
loads may ramp or reverse; zero normal load lifts that bow without silencing the
instrument. Commands outside the observed control window and station changes
refuse rather than disappearing or silently rebuilding geometry.

## Rotational bridges from physical footprints

All three input versions accept an optional block after ALL string/section/polarization records
and optional `thermal-supports`, before `cavity`:

```text
rotational-bridges COUNT
rotational-bridge STRING STATION POSITIVE_BODY_BRIDGE_STRING NEGATIVE_BODY_BRIDGE_STRING K_NM_PER_RAD D_NMS_PER_RAD REST_ANGLE_RAD
```

Each row joins the selected string's sourced twist to the body's scalar rocking
coordinate with one reciprocal angular spring/dashpot. The two body indices
select the original point footprints in a **v2 `--body-profile`**, in string
order. The compiler derives their separation from that profile's actual mesh
nodes and the performance's body dimensions, then projects their signed
displacement difference on the same material body eigenbasis. No modal angular
shape or invented string radius is accepted in the file. Changing the string's
shear material or the body's material/geometry therefore changes the joint's
physical response through the original source owners.

The ordered footprint pair declares the scalar rocking direction assumed to
align with this string torque joint. This is a small-angle chord rotation, exact
for rigid rocking and a secant approximation for a deforming plate. The input
does not supply or infer a full 3D string/body orientation or certify a local
bridge slope. Joint stiffness and damping remain explicitly supplied mechanical
data; bulk material identity does not determine a bridge contact law.

String and body receive opposite torques in their existing coupled solve, with
angular spring energy and dashpot loss included in its work ledger. Body motion
and pressure can consequently respond to twist after every bow is lifted. The
original transverse bridges, full body viscosity, shared gas and cavity remain
active. Omit the block, or select zero joints, to retain the original path.
Missing twist modes, absent/v1 profiles, coincident footprints, bad indices and
nonpassive or unrepresentable coefficients refuse before output creation.

The runnable `torsional-bridge.material-strings` and `torsional-bridge.body`
fixtures use a surface pluck/release and two strings on one material body. Build
them with the directional body pack as shown in
[BODY_PROFILES.md](BODY_PROFILES.md#reciprocal-rotational-bridges). They remain
synthetic integration data, not a measured bridge or instrument.

## V2: independent material layers and support constraints

V1 files remain accepted unchanged. V2 changes only the version header and string
records; audio, gas, body, force, cavity and canonical bow records are identical.
Each V2 string explicitly selects a support law and then a section construction:

```text
frankensim-material-strings-v2
... same audio, state, gas, body and strings-count records ...
string core.fsmatpk 0.5 0.0003 tension 20 4 3 0.8 3000 0.1
section bonded 0.0002 1
layer sleeve.fsmatpk 0.0003
... remaining original records ...
```

`section homogeneous` selects the original uniform specimen. A bonded section
uses the string record's pack for its solid core, then 1..=15 concentric sleeves
in strictly increasing outer-radius order. The final layer must end at the
string's declared outer radius. These are **perfectly bonded annuli**, not helical
windings, frictional sleeves or an approximation certified for wrapped strings.
Every layer supplies its own density, Young modulus and bending viscosity. When
twist modes are selected, every layer also supplies shear modulus and viscosity.
Missing sleeve data never falls back to the core. Each original source receipt,
loss band and spatial assignment remains attached to the shared string bank.
Repeated references load one immutable pack once, not one copy per layer.

Replace `tension 20` with `extension 0.4998 0.001` to hold the 0.5 m loaded span
against a 0.4998 m stress-free length with a declared 0.001 strain limit. The
original support owner computes tension from the section's total EA and imposed
strain. Replacing a layer's modulus therefore changes tension on this path; the
`tension` path holds the explicitly applied force instead. Slack, compression,
and out-of-limit strain refuse. No frequency target or automatic retuning is
introduced. The section geometry describes the resolved loaded state; Poisson
contraction, differential thermal strain and live axial creep are not supplied.

The `bonded-pluck.material-strings` fixture uses the existing wire pack and a
separate synthetic sleeve pack. After the original build commands, run:

```bash
cargo run -p xtask -- matdb-pack \
  --manifest "$demo/sleeve-manifest.tsv" --out "$demo/sleeve.fsmatpk"
cargo run --release -p fs-couple --example material_strings -- \
  "$demo/bonded-pluck.material-strings" "$demo/bonded-pluck.wav"
```

Change only the sleeve source or geometry to compare constructions. Plucks,
surface torque, simultaneous bows, shared body/cavity feedback, gas loading and
PCM all use their original owners. Torsion has no invented direct pressure gain:
it affects sound through physical surface/contact coupling and, when explicitly
selected, reciprocal rotational bridges. None of these fixtures is a measurement
of gold, silver or a real string.

## Material thermal expansion at fixed supports

A homogeneous V2/V3 string with `extension L0 STRAIN_LIMIT` can additionally
declare the temperature at which `L0` was stress-free. Put this optional block
after all string/section/polarization records and before other optional records:

```text
thermal-supports COUNT
thermal-support STRING_INDEX REFERENCE_T_K
```

The current solid temperature remains the performance's `state` temperature.
The compiler integrates the same material card's
`linear_thermal_expansion_coefficient` [1/K] from the reference temperature to
that current state, using the existing source owner. Scalar alpha and
piecewise-linear alpha(T) curves are supported; crossed knots are integrated,
not replaced by the endpoint coefficient times the temperature change.
The existing string binder computes `TENSION = EA * ((L-L0)/L0 - integral(alpha dT))`
with current material stiffness and prescribed current radius. Expansion receipts
remain attached through torsion, both bending planes, bows, body coupling and
the existing pressure renderer. There is no frequency normalization.

`thermal-bow.material-strings` uses the supplied synthetic wire expansion data:

```bash
cargo run --release -p fs-couple --example material_strings -- \
  "$demo/thermal-bow.material-strings" "$demo/thermal-bow.wav"
```

Change its `state` temperature or only the wire's alpha(T) data to change the
fixed-support tension and resulting sound. The gas has its independently
declared temperature. Omitting the block preserves ordinary extension/tension
semantics and does not query expansion data. Missing or ambiguous alpha,
out-of-domain temperature paths, slack/compression and excessive strain refuse;
negative expansion is admitted by the existing material law. Each string may
be named once, and force-controlled or layered strings refuse this selection.
These are frozen homogeneous thermoelastic comparisons. They do not solve heat
transport during playback, support/body expansion, radial contraction, layered
thermal mismatch, yielding or melting. The fixture is synthetic integration
data, not a calibrated instrument-temperature prediction.

## Causal gas loading of the moving strings

All three input versions accept one optional numerical-resolution record after
the complete string/section/polarization records, optional thermal supports and any rotational-bridge block,
before an optional `prepared-linear` row and then `cavity`:

```text
cylinder-gas-loading 0.001 30 12
```

The fields are `X_MIN X_MAX INTERVALS`. This selects the shared unsteady-Stokes
cylinder operator for every transverse mode in every string. The original
homogeneous or bonded section supplies its **actual outer radius**, and the
performance's existing gas supplies density and dynamic viscosity. No second
radius, gas record, material name switch or per-mode air-loss coefficient is
accepted here. Changing pressure, the species fractions or a layer pack therefore
rebuilds the mechanical response through those original physical sources.

The operator adds the surrounding gas's potential-flow inertia and a positive
spectral approximation of its viscous memory. Each memory branch is a real
inertial coordinate connected to the string by a dashpot. Both transverse planes,
the string/body reactions, fluid kinetic energy and viscous work participate in
the existing coupled solve. The old sampled string-air resistance is replaced,
so it is not counted twice; the original material bending viscosity and sourced
torsion remain attached. This selection precedes additional bows, optional
`--stop-profile` contacts, static preload and the canonical schedule. Memory
positions are unforced static gauges; a released preload accelerates the fluid
through the actual string coupling.

The dimensionless radial diffusion coordinate is
`x = radius * sqrt(gas_density * relaxation_rate / gas_viscosity)`.
`X_MIN` must lie in `(0, 0.5]`, `X_MAX` must be larger and at most `3000`,
and `INTERVALS` must be in `1..=4096`. The midpoint grid is uniform in
`log(x)`. The analytically integrated slow tail below `X_MIN` becomes a
direct resistance; the fast tail above `X_MAX` is truncated. These are three
independent approximation controls. Refine both cutoffs, interval count, the
string/body basis and time/control discretization before making accuracy claims.
A finer memory grid alone does not establish instrument fidelity.

The supplied `cylinder-flow.material-strings` fixture selects a V3 bonded string
with one bending mode in each of two planes, one twist mode, a two-mode body and
a cavity neck. Twelve memory branches per transverse mode add 24 coordinates,
for **30 complete-network coordinates**. Build its independent sleeve pack using
the earlier copied TSV files, then render it:

```bash
cargo run -p xtask -- matdb-pack \
  --manifest "$demo/sleeve-manifest.tsv" --out "$demo/sleeve.fsmatpk"
cargo run --release -p fs-couple --example material_strings -- \
  "$demo/cylinder-flow.material-strings" "$demo/cylinder-flow.wav"
```

Use a new sleeve-pack filename if that pack was already built. The fixture
coefficients and this deliberately modest resolution are synthetic integration
data. The native tests additionally combine two layered strings, separate regional
body footprints, full body viscosity, torque transmission, a cavity, and a
compliant receiver-bound stop in the same source-preserving path.

V3 retains V2's support and section records and requires one additional row after
each complete section: `polarizations single` or
`polarizations two BOW_X BOW_Y BRIDGE_X BRIDGE_Y`. The two direction vectors
must be unit vectors in the string's transverse plane. The fixture's
`polarizations two 0.6 0.8 1 0` excites both physical planes while its scalar
bridge footprint receives motion along the first axis. The circular material
operator supplies the same inertia and dispersion in both planes; their physical
coupling directions determine which motions reach contacts and the shared body.

The executable counts **all** requested material, body, neck and fluid
coordinates against its 128-mode bound. It computes a connection allowance for
every string bridge, angular joint, memory dashpot, body tensor column and cavity
compression/neck-loss channel. The original setup estimate
`N * (L + 1)^2` remains bounded by 16,777,216 work units and the shared connection
ceiling is 4096. A zero body-viscosity tensor may consume fewer columns than the
declared allowance. The ordinary exact propagator additionally admits at most
**42 modes in each connected linear island**; this stricter computational bound
can refuse a model that fits the file-level counts. No requested mode or memory
branch is silently removed. The render report includes the declared resolution
and complete numerical allocation.

The continuum assumes small transverse motion of a long, isolated circular
cylinder in an unconfined incompressible Newtonian gas, using local two-dimensional
flow along the string. It does not resolve finite-amplitude wakes, neighboring
strings, end/axial flow, rotary drag, rarefaction or compressible radiation.
The atmosphere remains fixed during one realization. Body/neck observation keeps
its existing compact-radiation approximation.

### Explicit preparation for larger coupled gas models

An optional record after `cylinder-gas-loading` (when present) and before
`cavity` selects the shared prepared exact propagator:

```text
prepared-linear MAX_MODES MAX_DENSE_BYTES MAX_MULTIPLY_TERMS
```

All three values must be positive. The shared owner accepts at most 512 modes;
this executable retains its stricter **128 complete-coordinate limit**.
`MAX_DENSE_BYTES` bounds additional dense preparation/cache storage, and
`MAX_MULTIPLY_TERMS` bounds cold scalar multiply-add work. Both are explicit
refusal budgets, not accuracy parameters or permissions to drop weak modes.

The selection is activated before torque joints and gas memory are attached.
It caches the complete fixed-period transition and separate physical resistance
integrals while retaining every original position, velocity, held-force input,
material state and signed connection. Fluid kinetic energy, intrinsic material
loss, body viscosity, contact compliance, pressure and per-connection work keep
their existing owners. The ordinary path remains selected when the row is absent.

`prepared-cylinder-flow.material-strings` has one layered two-plane string with
both bridge-direction components nonzero, one twist mode, two body modes and a
cavity. It requests 24 fluid branches per transverse mode: **54 coordinates**
and 51 linear connections. Its declaration is:

```text
cylinder-gas-loading 0.001 30 24
prepared-linear 64 67108864 200000000000
```

For the supplied positive-loss materials, the final operator's conservative
preflight estimate is about **37.6 MB of additional dense storage and
161.2 billion scalar multiply-add terms**, within its explicit 64 MiB and
200-billion-term budgets. The work estimate reserves worst-case Taylor products
and scaling levels; it is not elapsed runtime or a real-time performance claim.
Preparation can be expensive even though the accepted fixed-period transition is
cached. The native regression therefore performs one such construction and only
32 release samples.

With the existing wire, body and sleeve packs built, render to a new output:

```bash
cargo run --release -p fs-couple --example material_strings -- \
  "$demo/prepared-cylinder-flow.material-strings" "$demo/prepared-cylinder-flow.wav"
```

Omitting only `prepared-linear` from this fixture produces the ordinary
connected-block capacity refusal. Undersized mode, byte or work allowances also
refuse. Increasing the prepared capacity does not establish convergence of the
fluid spectrum or the instrument basis.

The initial string/body bank must still fit its ordinary constructor before this
cold selection can be applied. The option admits larger subsequent torque/fluid
connections; it does not rescue an initial bank already refused by that
constructor. Generic callers needing an initially larger network can explicitly
use `CoupledModalSystem::new_prepared` or `new_with_columns_prepared` with their
own complete models and budgets.

## Nonlinear stretching with the same bow, body and gas memory

An optional record after `prepared-linear` (or the cylinder-loading record when
preparation is omitted), before `cavity`, selects source-derived geometric
stretching:

```text
geometric-stretching 0.2 64 100000 1e-11 1e-10
```

The fields are the maximum continuous-span slope, maximum constitutive
iterations, maximum generalized reaction [N/sqrt(kg)], absolute force residual
tolerance [N/sqrt(kg)] and relative force residual tolerance. The original
material or bonded section supplies axial rigidity EA; there is no independent
nonlinear stiffness or pitch parameter. Omitting the record retains the
original linear route.

For retained longitudinal inertia, select `longitudinal-waves` in that same
position instead of `geometric-stretching`:

```text
longitudinal-waves 0.2 0.01 64 100000 1e-11 1e-10 2 1 2
```

The fields are maximum transverse slope, maximum absolute dynamic axial strain,
the four constitutive iteration/force tolerances used above, then a count and
strictly increasing one-based axial sine indices. This record selects the same
indices for each original string; the Rust API accepts a separate selection per
string. Both transverse planes contribute to the same material strain
`u_x + (y_x² + z_x²)/2`. The existing `fs-nlmodal::LongitudinalString` supplies
the strain energy and gradient; no separate instrument solver is introduced.

Solid material density and EA determine the longitudinal frequencies
`omega_j = j*pi/L*sqrt(EA/mu_solid)`. Transverse motion retains its actual
gas-loaded mass normalization. Swapping a homogeneous material or a bonded
layer changes these coefficients automatically; changing the gas does not
assign transverse added fluid mass to axial motion. Every requested axial mode
must fit the original material elasticity/density applicability, clock and
complete-network budget. Unknown/out-of-band physics refuses without retuning
or silently deleting modes.

The axial oscillators join the original exact linear network. Their complete
positive strain energy replaces their linear elastic storage once, through the
same bounded discrete-gradient/friction solve. All original body/cavity,
torsion, bridge, gas-memory and pressure coordinates remain present. Both
reaction convergence and the complete work/energy balance must close. The
reported nonlinear correction can be signed because the exact linear owner
already includes axial stiffness; the complete physical energy remains positive.

`longitudinal-cylinder-bow.material-strings` is a runnable bonded-string example
with two transverse planes, torsion, gas memory, body, cavity and retained axial
modes. Build the existing wire/body/sleeve packs and render:

```bash
cargo run --release -p fs-couple --example material_strings -- \
  "$demo/longitudinal-cylinder-bow.material-strings" "$demo/longitudinal-cylinder-bow.wav"
```

The axial extension is explicitly elastic and lossless. Original bending,
torsional, bridge, body and gas losses remain active; bending viscosity is not
reused as an unsupported axial loss law. Unretained axial channels equilibrate
statically, so a two-mode selection is not a fully resolved longitudinal wave
field. Fixed endpoints, moderate slopes, small constitutive strain and positive
local tension remain required. There is no winding slip, temperature evolution,
nonlinear static preload or independent longitudinal end-force port into the
body. Axial strain reaches the receiver through its reciprocal effect on
transverse motion and the original bridge connections. This is not a claim of
measured violin/cello/bass equivalence.

`nonlinear-cylinder-bow.material-strings` combines a bonded string's two
transverse planes and material torsion with causal cylinder-gas memory, the
shared material plate, its cavity and an actual bow gesture. With the existing
wire, body and sleeve packs built, render to a new output:

```bash
cargo run --release -p fs-couple --example material_strings -- \
  "$demo/nonlinear-cylinder-bow.material-strings" "$demo/nonlinear-cylinder-bow.wav"
```

Both transverse planes contribute to one shared strain and tension. The same
loaded string mass normalizes the geometric storage, including selected gas
added inertia. The original complete linear network and scalar or simultaneous
friction solver participate in each `fs_phs::discrete_gradient` constitutive
trial; both reaction residual and full storage/work balance must close before
state is published.
The acceleration observer uses the endpoint geometric force. Existing gesture
scheduling and PCM writing consume this same state.

`nonlinear-double-bow.material-strings` uses two existing canonical bow tracks
to drive two bonded, two-plane, torsional strings through the same body, cavity
and fluid memory. Each bow retains its own traction, actual slip, work and loss;
all bows see the same constitutively accepted mechanical endpoint. No new
performance format or independent nonlinear voice is needed:

```bash
cargo run --release -p fs-couple --example material_strings -- \
  "$demo/nonlinear-double-bow.material-strings" "$demo/nonlinear-double-bow.wav"
```

`nonlinear-contact-bow.material-strings` combines the same source-derived
stretching, bonded string, two transverse planes, torsion, body, cavity and
gas memory with the existing normal-contact profile. The profile presses the
support into the string and withdraws it; the jointly solved normal reaction
sets the bow's friction load, including automatic release on separation:

```bash
cargo run --release -p fs-couple --example material_strings -- \
  "$demo/nonlinear-contact-bow.material-strings" "$demo/nonlinear-contact-bow.wav" \
  --stop-profile "$demo/contact-bow.stops"
```

Independent compliant stops and contact-loaded bows use the original joint
normal/tangential solver inside each stretching trial. Both endpoints of contact
storage, normal-contact loss and signed support work `-reaction * offset_change`
join the geometric, network and bow work balance before publication. Normal
translation uses the existing `motion`/`knot` block. A separate `station-motion`
block can additionally move independent finger/stop contacts along the original
string span, as described below. Existing stop profiles and gesture tracks remain
the input format.

This is the moderate-slope, fixed-end Kirchhoff–Carrier model. Nonlinear
harmonics and strong coupling require clock/mode refinement; residual closure
alone does not establish trajectory accuracy. This selection supports multiple
prescribed or contact-loaded bows, fixed, normally translating or migrating independent supports and
`forces retain`. Static preload still refuses because it requires a joint
nonlinear equilibrium solve. Select gas loading and rotational bridges before
attaching stretching. Axial wave
inertia, winding slip, evolving temperature and measured instrument qualification
remain outside this composition.

## Moving finger/stop stations in the same material instrument

`migrating-finger.material-strings` and `migrating-finger.stops` combine a fixed
contact-loaded bow with an independently pressing, migrating and withdrawing
compliant finger. The string retains four bending modes in each transverse
plane, material torsion, its bonded section's axial rigidity, shared material
body/cavity and causal gas memory. With the existing wire, body and sleeve packs:

```bash
cargo run --release -p fs-couple --example material_strings -- \
  "$demo/migrating-finger.material-strings" "$demo/migrating-finger.wav" \
  --stop-profile "$demo/migrating-finger.stops"
```

After the stop and optional `bow-contact` declarations, a v1 stop profile may
contain each of the normal and longitudinal motion blocks once, in either order:

```text
station-motion 1 0.00003 1048576 3
station-knot 480 0.2 0.43
station-knot 1920 0.2 0.48
station-knot 2880 0.2 0.46
```

The header fields are maximum longitudinal speed [m/s], maximum distance per
mechanical sample [m], maximum geometry/response work terms per update, and knot
count. Each knot gives the elapsed **endpoint sample** followed by one interior
station fraction per stop in declaration order. Sample zero implicitly retains
each authored station. The original endpoint sampler linearly interpolates
between knots and holds the final position. The `motion` block still supplies
independent normal offsets on that same clock. Both paths admit at most 128
knots and 4096 scalar coordinates; the complete stop profile remains at most
16 KiB. Native geometry updates admit at most 1,048,576 work terms.

The runtime obtains every new displacement/force map from the existing material
string source, including its gas-loaded mass and selected transverse directions.
It retains the full supported span and every string, torsion, body, cavity and
fluid coordinate. A declared receiver footprint and the contact normal stay
fixed. No speaking length, tension, modal frequency or observer gain is edited.

For normal closure `x = Nᵀq - s`, the endpoint identity is
`Δx = mean(N)ᵀΔq + ΔNᵀmean(q) - Δs`. The original contact reaction therefore
acts through `-mean(N)`, and its prescribed geometry work is
`reaction * (ΔNᵀmean(q) - Δs)`. The original unilateral potential and damping
use the complete closure change. Both geometry supplies, the actual normal
reaction, bow friction and source-derived stretching must close one common
energy balance before publication. Acceleration pressure uses the endpoint
attachment map and endpoint geometric force. Reports expose attachment work
separately from total geometry work; rejected motion changes neither accepted
station, normal offset, modal state, contact report nor sample clock.

A stop that supplies a bow's normal reaction must retain that bow's station;
its tangential migration/slip is not part of this model. Other independent
stops can migrate alongside it. Excessive speed or displacement, endpoints
outside the span, nonpositive endpoint/mean-map contact compliance and existing
material, force, penetration, energy or pressure limits refuse. Constant
station paths preserve the original fixed-stop arithmetic. This is a bounded
moving compliant contact model, not a solved finger, distributed tissue patch,
rigid stopping constraint or measured violin/cello qualification. Resolve the
contact footprint and refine spatial modes and time before treating a resulting
pitch or timbre as converged; four retained modes per plane are an example
configuration, not an accuracy claim.

## Material profile and physical limits

The input profile queries **bare SI axes `T` and `omega`**, at the explicitly
declared solid state; it does not guess other conditioning coordinates or strip
incompatible quantity metadata. Packs with incompatible typed axes or additional
required contexts need the existing general material-state API, not coercion.
All claims use the original `SingleClaimOnly` policy. No ambiguous property,
missing loss, out-of-domain mode, or absent shear data becomes a fallback value.

Every homogeneous string and every bonded layer needs `density` [kg/m3], `young_modulus` [Pa], and
`kelvin_voigt_bending_viscosity` [Pa s]. Positive `TWIST_MODES` additionally requires
`shear_modulus` [Pa] and `kelvin_voigt_shear_viscosity` [Pa s]; zero does not query or
infer those properties. The body additionally needs `poisson_ratio` [1]. Bending
and shear loss must have an explicit valid angular-frequency band; all actual
retained modes, not only the source query frequency, undergo the existing checks.

This executable selects homogeneous or bonded circular strings, fixed tension or
fixed-extension supports and fixed rotational ends. Its shared body is the legacy
simply supported isotropic flat plate or the explicit material chart described in
[BODY_PROFILES.md](BODY_PROFILES.md), with an optional uniform-pressure
cavity/circular neck. Bridges are interior transverse connections on the full
span **including afterlength**, plus explicitly selected scalar rotational
spring/dashpot joints. It does not model carved wood shells, a full 3D unilateral
bridge-contact law, f-hole geometry, distributed cavity
modes, helical winding/slip, solved finger disengagement,
bow-hair normal-load distribution, or exterior radiation back-loading. The body
and neck observations retain the existing compact baffled, zero-retardation
approximation. Source-state validity is not measured instrument-fidelity evidence.

Input caps: 1 MiB performance; 8 MiB per distinct material pack and 64 MiB total
raw packs; 8 strings, at most 16 layers per string, 16 bow contacts, 16 rotational
bridges, 64 force ports, 128 complete-network modes including fluid memory,
16 plate modes; 16,384 events per input family; 262,144 compiled controls per
family; 16,777,216 source/projection work units per compiler; at most 192 kHz,
600 seconds or 28,800,000 samples. Existing numerical owners retain their stricter
limits. The example additionally caps total mechanical storage at 100 J,
connection forces at 1000 N, bow traction at 100 N, slip at 10 m/s and observer
pressure at 1 MPa. These are refusal ceilings, not assertions of model accuracy.

## Native checks

```bash
cargo test -p fs-couple --example material_strings -- --test-threads=1
```

The tests exercise actual normalized-pack decoding, material-only substitutions,
shared gas response, source-derived cavity-neck inertia/loss and end conditions,
causal cylinder memory, joint pluck/double-bow control, callback equivalence, exact
existing PCM output, and source/input refusals. The fixture TSVs use the existing
`xtask matdb-pack` compiler, not a new importer.
