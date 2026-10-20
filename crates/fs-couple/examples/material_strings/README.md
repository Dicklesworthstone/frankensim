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

`cavity neck VOLUME_M3 RADIUS_M LENGTH_M RESISTANCE_PA_S_M3 LISTENER_M` selects the
existing Helmholtz model with an explicitly specified neck resistance. It replaces
`cavity none`. Compression, neck motion, body motion and bridge reactions remain
in the same mechanical solve.

`friction stribeck MU_STATIC MU_KINETIC CHARACTERISTIC_SPEED_M_S VISCOUS_PER_SPEED`
is also accepted. These coefficients retain **caller authority**. Bulk material
data do not certify a rosin/string friction law. This input profile does not yet
load external dry-interface cards; the existing Rust API supports that binding.

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
and before `cavity`:

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

## Causal gas loading of the moving strings

All three input versions accept one optional numerical-resolution record after
the complete string/section/polarization records and any rotational-bridge block,
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
added inertia. The original complete linear network and scalar friction solver
participate in each `fs_phs::discrete_gradient` constitutive trial; both reaction
residual and full storage/work balance must close before state is published.
The acceleration observer uses the endpoint geometric force. Existing gesture
scheduling and PCM writing consume this same state.

This is the moderate-slope, fixed-end Kirchhoff–Carrier model. Nonlinear
harmonics and strong coupling require clock/mode refinement; residual closure
alone does not establish trajectory accuracy. This selection currently supports
one prescribed bow and `forces retain`. Additional simultaneous bows, independent
stops (including moving stops and contact-loaded bows) and linear static preload
refuse because they require a joint nonlinear contact/equilibrium solve. Select
gas loading and rotational bridges before attaching stretching. Axial wave
inertia, winding slip, evolving temperature and measured instrument qualification
remain outside this composition.

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
shared gas response and causal cylinder memory, joint pluck/double-bow control, callback equivalence, exact
existing PCM output, and source/input refusals. The fixture TSVs use the existing
`xtask matdb-pack` compiler, not a new importer.
