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

`material_strings` is also an installable executable. Its binary target and the
example below compile the same frontend and use the same material, contact,
body and atmosphere owners:

```bash
cargo install --locked --path crates/fs-couple --bin material_strings
material_strings PERFORMANCE NEW_OUTPUT.wav --body-profile PROFILE
```

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

### Material shear layers in the original contact solve

`shear-layer-bow.material-strings` selects an isothermal elastic-plastic layer
through source data and physical geometry:

```text
friction layers 1
shear-layer bow manufactured-shear-layer.fsmatpk 0.000001 0.0001 0.05
bow_ribbon 0.01 3
```

Each `shear-layer` row contains **track ID, bulk material-pack path, area [m²],
thickness [m], maximum elastic engineering shear**. Supply exactly one row per
authored bow track. The area describes the whole ribbon; each cell receives
one share of area, stiffness, yield force and plastic resistance. Thickness
and the original material source remain unchanged.

After the wire/body pack commands above:

```bash
cargo run -p xtask -- matdb-pack \
  --manifest "$demo/manufactured-shear-layer-manifest.tsv" \
  --out "$demo/manufactured-shear-layer.fsmatpk"
cargo run --release -p fs-couple --example material_strings -- \
  "$demo/shear-layer-bow.material-strings" "$demo/shear-layer-bow.wav"
```

The pack must supply `shear_modulus` [Pa], `shear_yield_stress` [Pa] and
`plastic_shear_viscosity` [Pa·s]. The shared material owner derives
`k = A G/h`, `Y = A tau_y` and `c = A eta_plastic/h`. Elastic strain stores
energy; plastic flow dissipates yield work plus viscous work. Elastic
unloading can return energy to the string/body system. The same simultaneous
contact solve, optional fiber compliance, body/cavity motion and pressure
observer carry that response into PCM. Changing the pack or layer geometry
does not require changing Rust.

The `state` row supplies the frozen temperature and material query. Every
selected property needs explicit temperature support and applicability from
static storage through the sample clock's Nyquist frequency. Constants and
temperature-only curves are supported. Pressure-, rate-, strain- or
frequency-dependent coefficients require a more general constitutive owner
and refuse here. Oscillatory viscosity and friction coefficients cannot
substitute for plastic shear viscosity.

This first layer model requires retained positive compression, fixed layer
geometry/tangent frame, and an elastic shear limit in `(0, 0.1]`. Prescribed
bow loads must stay positive; contact-loaded bows instead author zero extra
load and use the simultaneously solved normal reaction. Opening refuses
atomically, preserving stored layer energy and accepted history. Layer
opening/recontact, wear, thermal evolution, migration and nonlinear string
geometric storage are not yet composed. Prepared-linear propagation is
supported within its existing independent restrictions.

The supplied layer is **manufactured integration data**. Its constants are
not measurements of rosin or evidence of realistic violin/cello sound.
Real rosin needs a supported constitutive source, including the relevant
temperature and rate behavior. The executable retains the complete pack,
query receipts, geometry and derived layer coefficients in the same source
catalog and diagnostics.

### Material fibers in the same bow/string/body solve

`fiber-ribbon.material-strings` adds an explicitly described fiber bundle to
the original ribbon contact. After the copy, wire-pack and body-pack commands
above, compile its additional material source and run the same executable:

```bash
cargo run -p xtask -- matdb-pack \
  --manifest "$demo/fiber-manifest.tsv" --out "$demo/fiber.fsmatpk"
cargo run --release -p fs-couple --example material_strings -- \
  "$demo/fiber-ribbon.material-strings" "$demo/fiber-ribbon.wav"
```

The optional record follows `bow_ribbon WIDTH_M CELLS`:

```text
bow_fibers PACK RADIUS_M LENGTH_M FIBER_COUNT TOTAL_TENSION_N CONTACT_FRACTION AXIAL_MODES NORMAL_MODES
```

Each cell receives its quadrature share of the **original fiber count**, mass,
axial stiffness, bending stiffness, tension and material losses. Independent
axial and transverse wave modes follow the source section and length. Contact
traction stretches the fibers, and their response changes slip, string motion
and body pressure in the original joint solver. `FIBER_COUNT` describes the
modeled interacting bundle; a complete multilayer bow hank does not imply that
all its hairs touch the string.

The executable also admits the fibers' **normal** receiver maps. With the same
material packs, the second fixture closes contact compliance, fiber deflection,
the resulting normal force, tangential stick/slip and physical pressure in one
joint solve:

```bash
cargo run --release -p fs-couple --example material_strings -- \
  "$demo/fiber-contact.material-strings" "$demo/fiber-contact.wav" \
  --stop-profile "$demo/fiber-contact.stops"
```

Its stop profile uses the original ribbon cell identity for both location and
receiver, for example:

```text
stop centre 0 ribbon:0 normal 0 -1 fiber:0 -0.000001 5000 1 0 10 0.001
bow-contact 0 centre
```

`ribbon:CELL` uses the runtime's exact midpoint geometry; `fiber:CELL` binds
that same cell's material-derived normal motion. Cell zero is the centre, then
the remaining cells run from left to right. Each cell needs its own normal law
and `bow-contact` binding. The string must retain two transverse polarizations;
the explicit contact normal must be perpendicular to its bow tangent and have
the correct torsional orientation. The schedule supplies only the common
tangential speed and **zero normal force**. Each actual normal reaction feeds
its corresponding friction law directly; it is not divided by the cell count
again. Separation removes that cell's traction automatically.

The negative gap in this manufactured fixture is an initial elastic
interference of one micrometre. Its stored contact energy is included from the
first step; it is not a claimed measured bow force or a settled preload.
Changing fiber tension changes transverse compliance, solved load and sound.
Changing the gap to a positive clearance leaves initially stationary separated
surfaces unforced. Contact stiffness is authored **per cell**, so its physical
quadrature must also be preserved when refining contact count.

Both fixtures use fixed interior contact stations and rigid fiber endpoints
in a uniformly translating tangential support frame. Normal offset trajectories
on fiber receivers refuse: moving the endpoints would also require their
frame-acceleration forces, which this frozen model does not supply. Ordinary
fixed/bridge stop motion remains available. The original `fiber-ribbon`
fixture retains prescribed uniform normal load.
It refuses speed changes and nonlinear-string options until those combinations
have the corresponding moving-support/assembly implementation. It does not
model flexible bow stick, contact migration, inter-hair rubbing or thermal
rosin. Axial and normal mode counts, ribbon cells and physical sample rate are
independent refinement axes. A selected mode outside source validity or the
clock's Nyquist guard is refused, never discarded.

The source pack must explicitly provide `density`, `young_modulus`,
`kelvin_voigt_young_viscosity` and `kelvin_voigt_bending_viscosity` with frequency
support. Zero viscosity must be declared. The included `fiber.tsv` is
**manufactured validation data**, not horsehair measurements or an asserted
violin improvement. The primary-source horsehair pack under
`data/matdb/seed-v1/bow-hair-horse-pitteroff-woodhouse-acustica1998` contains
static/dynamic measurements and geometry, but its 10–100 Hz loss range does
not establish a broadband Kelvin–Voigt law for these modes. It is not silently
extrapolated or converted into an arbitrary damping preset.

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

### A winding has its own wire geometry and material

`wound-pluck.material-strings` uses the explicit record

```text
section close-wound CORE_RADIUS WRAP_PACK WIRE_RADIUS AXIAL_PITCH
```

The string row's radius is the finished radius: core radius plus two winding
wire radii. The example uses a 0.2 mm core radius, a 0.05 mm wire radius and
0.101 mm pitch. Set torsional mode count to zero for this reduction. After the
usual fixture-pack commands above, build a separate winding pack and render:

```bash
cargo run -p xtask -- matdb-pack \
  --manifest "$demo/wire-manifest.tsv" --out "$demo/wrap.fsmatpk"
cargo run --release -p fs-couple --example material_strings -- \
  "$demo/wound-pluck.material-strings" "$demo/wound-pluck.wav"
```

The wrap resolves only its density and shear modulus. Actual helical wire
length determines its mass; close-coiled spring mechanics determine its
axial stiffness. The core and wrap share end displacement, with sliding allowed
along their interface. A free-ended wrap needs a different axial boundary law.
The core supplies bending stiffness and its explicitly selected Kelvin–Voigt
viscosity or generalized-Maxwell relaxation spectrum. The latter uses the same
`bending-loss prony` record as a homogeneous string and requires no viscosity
claim. Each original core arm remains an internal strain in the joint
string/body/cavity solver, projected with the **total core-plus-wrap mass**.
The wrap remains elastic in this approximation; it never inherits the core's
relaxation or counts as a bonded bending layer. Replacing the winding's pack therefore changes mass, tension under
fixed extension, frequencies, gas loading and the shared body's physical
pressure. Explicit fixed tension and fixed extension keep their different
experimental meanings. The exposed wrap material also supplies the identity
for any selected bow or bridge interface.

This is the **sliding, core-dominated, close-coiled approximation** described by
[Kemp (2017)](https://doi.org/10.1371/journal.pone.0184803), with exact helix
mass geometry. It does not treat the winding as a bonded solid sleeve. The
admission limits are wire/coil radius ≤ 1/4, pitch/circumference ≤ 0.15 and
maximum modal wavenumber × pitch ≤ 0.2; neighbouring turns must not intersect.
Those limits are declared scale assumptions, not measured error guarantees.
Resolved winding friction, wrapping prestress, winding loss, thermal mismatch,
torsional response and longitudinal winding waves remain unavailable. Requests
for those extensions refuse instead of borrowing homogeneous-string laws.
This fixture demonstrates the constitutive path; it does not validate a real
instrument or establish acoustic equivalence to a resolved wound string.

For a sourced bulk-silver wrapping material, compile the bounded NIST-based pack
and use the same geometry, core and force program:

```bash
cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/silver-pure-wrap-nist1995/manifest.tsv \
  --out "$demo/silver-wrap.fsmatpk"
cargo run --release -p fs-couple --example material_strings -- \
  "$demo/silver-wound-pluck.material-strings" "$demo/silver-wound-pluck.wav"
```

This example retains the manufactured core and body fixtures. Only the wrap
changes. The pack combines the NIST review's bulk-silver reference density with
an explicitly derived isotropic shear modulus at 293.15 K; it supplies no silver
damping. The winding remains lossless under the declared reduction, and this
example does not represent a measured silver-wound instrument.

### Material bending memory on coupled strings

`relaxing-strings.material-strings` selects one source-defined Maxwell arm on
each of two homogeneous strings. Each retains four bending modes and two
torsional modes at its declared tension, 20 N or 22 N. The original pluck-and-bow
schedule drives the first string; both load the same body and physical cavity
neck in the same gas mixture. Build its synthetic wire pack and render after
the copy and body-pack commands above:

```bash
cargo run -p xtask -- matdb-pack \
  --manifest "$demo/relaxing-wire-manifest.tsv" --out "$demo/relaxing-wire.fsmatpk"
cargo run --release -p fs-couple --example material_strings -- \
  "$demo/relaxing-strings.material-strings" "$demo/relaxing-strings.wav"
```

The optional per-string constitutive selector is accepted in every input version:

```text
bending-loss prony COUNT EXCESS_MODULUS_PROPERTY RELAXATION_TIME_PROPERTY ...
```

Place it immediately after the `string` row in v1, after the complete `section`
declaration in v2, or after `polarizations` in v3. Supply 1–8 property pairs;
moduli are nonnegative excess stiffness in Pa and times are positive seconds.
Omission or `bending-loss kelvin-voigt` retains the existing viscosity selection.
Prony selection requires an explicit positive `equilibrium_young_modulus` that
equals `young_modulus` exactly. The latter describes the relaxed backbone.
The selected path does not require `kelvin_voigt_bending_viscosity`;
`relaxing-wire.tsv` deliberately has none. Torsion retains its independently
selected shear modulus and shear viscosity.

The existing `bending-memory` / `string-memory` block before `cavity` remains
available for additional per-layer memory alongside explicitly selected
Kelvin–Voigt loss, including perfectly bonded sections and close-wound cores. See
[the shared string-memory interface](../../../../docs/MATERIAL_STRING_BENDING_MEMORY.md).
The two selections may coexist on different strings and share body histories,
gas loading, bow contact and PCM output. Selecting the same string through both
records refuses before material loading, preventing duplicate physical loss.
All declared string and body histories share one setup-work reservation.

The fixture declares density 7800 kg/m³, equilibrium modulus 200 GPa, excess
bending modulus 80 GPa and relaxation time 0.001 s. Every property covers
280–320 K and 0–200000 rad/s. Source admission requires coverage from zero
through the complete coupled instantaneous frequency bound, including later gas
and body changes. Every positive arm also requires `dt / tau <= 2`; at 48 kHz
this example has `dt / tau = 1/48`. Unsupported coefficients or bands refuse
before output creation.

The original shared solver carries eight internal material histories alongside
15 mechanical modes, including torsion, body and cavity. Changing a selected
modulus or time in the TSV, rebuilding its pack and changing only the pack path
in the performance changes the constitutive history loading the common body.
Its pressure is what the existing PCM writer receives. Static preload
initializes relaxed histories, and neighboring-string motion uses the same
reciprocal mechanics. The exclusive per-string selection supports homogeneous
strings and close-wound cores on the exact solver; bonded sections use the
existing additional-memory interface and `prepared-linear` refuses constitutive
histories. Core relaxation in a wound string keeps the original helix mass,
finished gas-loading radius and wrapping material identity. All participating
source bands must cover the complete instantaneous coupled frequency bound. These models
supply bending relaxation; the torsional selection below supplies independent
shear histories. Axial relaxation, interlayer slip and measured wound-string
behavior remain outside their scope. All supplied values
are synthetic integration data, not measured instrument properties.

### Material shear memory on torsional strings

Torsional modes can select their own source-backed relaxation spectrum for each
physical layer. Put this record after the optional `bending-loss` row (or in its
place when bending loss is omitted), after the string's complete section and
polarization declaration:

```text
torsion-memory 2
torsion-layer shear_modulus 0
torsion-layer shear_modulus 2 coating_delta_1 coating_tau_1 coating_delta_2 coating_tau_2
```

The count must equal the number of physical layers: the core first, then sleeves
from inside to outside. A homogeneous string supplies one layer. Each layer names
its positive equilibrium shear modulus and 0–8 pairs of excess shear modulus in
Pa and positive relaxation time in seconds. The equilibrium value must equal
that layer's original `shear_modulus` exactly. The example above leaves the core
elastic and gives the coating two independent arms; the core pack needs none of
the coating's properties. Each layer may use different names, times and arm
counts. Zero-arm layers still admit the equilibrium property's frequency band.

The material section owner weights every arm by that layer's actual polar second
moment. Its original solid rotary inertia, torsional mode shapes and parallel
shear viscosity remain attached. The existing reciprocal solver carries one
angular history per positive arm and retained torsional mode; adding a second
transverse polarization does not duplicate these histories. They coexist with
bending and body memory, gas loading, longitudinal coordinates, bow contact and
the shared body/cavity pressure output.

A performance change can therefore consist only of replacing a core or sleeve
pack, or changing its property selectors. The pack resolver asks only the
corresponding layers for those claims, and the existing performer and PCM writer
receive the resulting mechanics. Omission preserves the existing torsional law;
an explicit zero-arm selection preserves its arithmetic without dummy material
properties. Missing data, conflicting dimensions, repeated source pairs,
unsupported source bands and nonpositive times refuse instead of substituting
an elastic response.

Torsional modes must be present. All declared bending, torsional and body
histories share the original 256-state and setup-work limits, reserved before
material loading. Torsional memory requires the exact history owner;
`prepared-linear` refuses it. This is small-strain torsion of perfectly bonded
concentric sections with explicit linear shear relaxation. It does not model
winding slip, finite-strain axial relaxation or measured instrument fidelity.

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

`cavity exterior-neck` takes the same fields and hands the outer end to a
computed exterior load:
- **What the neck keeps:** its mass to the inner end only, L + δ_in, and the
  same compressibility.
  - Its Bessel viscothermal resistance is frozen at the isolated frequency of
    the two declared ends.
- **What the neck drops:** the compact outer end correction and the compact
  radiation resistance. The outer end's mass and radiation come from the
  `body-radiation skin` load, through the skin's aperture caps (below).
- Without that load the neck has no outer mass. So `exterior-neck` refuses
  without `body-radiation skin`, and a compact `neck` or `physical-neck`
  refuses beside it, because its outer end would be counted twice.
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

The material-layer form above selects bulk constitutive properties and layer geometry.
A separate supported choice replaces the friction row with source-bound ordered interfaces:

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
The original dry-friction schedule, including ramps and release, survives unchanged.
Material shear layers retain ramps but require positive compression throughout.

With source-bound friction, the declared nominal area belongs to the **whole
ribbon**. Area and normal load are divided by the same cell count, retaining
the same nominal pressure and original ordered interface card at each cell.
Runtime applicability uses each cell's actual slip. Numerical setup, control
storage and compilation work are charged for every cell under the same global
limits. Independent finger/obstacle stops remain available. A `bow-contact`
binding may instead derive a cell's load from a co-located stop; the original
gesture must then prescribe zero normal force. Its source area is still that
cell's quadrature share, and applicability uses the actual solved load.
The rigid limit does not infer compliant bow hair or evolving rosin; use
`bow_fibers` for retained material-fiber motion. Refine spatial cells,
structural modes and time together.

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

### Tangential contact with an existing material body

A stop profile can return a bow contact's tangential reaction to a physical
receiver already present in the material instrument:

```text
frankensim-string-stops-v1
bow-surface 0 bridge:0
```

Pass that file with `--stop-profile`. The first index selects the zero-based
bow track in the original schedule. `bridge:STRING` selects the actual body
footprint used by that string's sole bridge connection; it does not select
the string's displacement or create another body. The footprint is projected
from the receiver's original material basis. Changing the receiver material
therefore updates the contact's inertia, compliance, loss and force/displacement
map in the same solve that produces string motion and body/cavity pressure.

The selected bridge axis must be collinear with the contacted string's bow
tangent in the profile's common transverse frame. Opposite axes reverse the
receiver map's sign; oblique axes refuse because one scalar footprint cannot
describe motion in an unrepresented direction. Single-plane sources use their
existing first transverse axis. Two-plane sources retain their explicitly
declared bow and bridge directions. Source-defined friction and normal-load
laws keep their original applicability requirements.

These records can be the entire profile, or accompany independent `stop` and
`bow-contact` records before the motion blocks. Surface binding precedes normal
contact and static preload. Each bow may bind once; distinct bows may share a
receiver and retain their coupled reactions. Absent indices, repeated bow
bindings and uniform ribbon profiles refuse before material loading. Ribbons
do not declare individual receiver footprints in this format. If a
`bow-contact` also supplies this bow's normal load, its normal receiver must
belong to the same physical body, with a supported perpendicular footprint.
One scalar bridge cannot supply both perpendicular rows; use distinct authored
footprints on that body when both directions are dynamically resolved. A fixed
normal stop cannot supply this dynamically resolved receiver's normal contact.

The gesture speed remains an additional prescribed relative actuator speed;
use zero when the relative surface motion should come entirely from the solved
mechanics and physical applied forces. This profile describes a string rubbing
an existing material body. It does not infer a bow-hair structure, rotated
receiver frame or measured instrument qualification. The existing observation,
gesture scheduling, contact energy accounting and PCM writer use the same
accepted mechanical state.

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

All three values must be positive. The shared owner accepts at most 1024 modes,
and so does this executable's complete-coordinate limit.
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

## Computed exterior reaction on the shared body

The optional `body-radiation` record before `cavity` selects a geometry- and
gas-derived load on the original material body. An optional `shell_observer`
record may follow it. Regional flat plates use the
Rayleigh half-space model; closed curved shells and assemblies use the shared
Helmholtz boundary-element solver:

```text
body-radiation rayleigh SIDES MIN_HZ MAX_HZ FREQUENCIES MAX_CELLS MAX_SWEEP_TERMS MAX_POLES
body-radiation bem MIN_HZ MAX_HZ FREQUENCIES MAX_PANELS MAX_SWEEP_TERMS MAX_POLES
shell_angular_refinement LEVEL
body-radiation skin MIN_HZ MAX_HZ FREQUENCIES MAX_PANELS MAX_SWEEP_TERMS MAX_POLES LEVEL MAX_ANCHOR_M
```

`rayleigh-pluck.material-strings` and `orthotropic.body` provide a small runnable
example. After the initial copy and wire-pack commands above:

```bash
cargo run -p xtask -- matdb-pack \
  --manifest "$demo/directional-body-manifest.tsv" --out "$demo/directional-body.fsmatpk"
cargo run --release -p fs-couple --example material_strings -- \
  "$demo/rayleigh-pluck.material-strings" "$demo/rayleigh-pluck.wav" \
  --body-profile crates/fs-couple/examples/material_strings/orthotropic.body
```

The original retained, mass-normalized DKT displacement fields supply three
area points per triangle. Every signed body mode enters the full Rayleigh
impedance matrix, including modes whose integrated displacement is zero.
The shared positive-real fitter preserves mutual modal reactions; its acoustic
storage and loss join the same exact string/body/cavity/bow solve. Radiation
therefore alters actual string and body motion before the pressure observer.
Material, grain, thickness, gas pressure and gas-composition changes rebuild
this load from their physical data. No damping ratio or output gain substitutes
for the air reaction. The native report retains the original specimen, gas,
complete mode/cell/pole counts and sampled complex/resistance fit errors.

`SIDES` is one or two independent infinite-baffle half-spaces; a body with a
modeled interior cavity requires one exterior side. Supply an odd 33–257
frequency count. The band must cover all retained free-body frequencies, and
the fitter's highest 8x off-band storage poles must resolve on the mechanical
clock. All declared cells, fitted poles and existing material histories count
against the original work and exact-state limits; none are silently dropped.
The complete fit keeps the shared 5% peak/2% RMS complex and 10% peak/5% RMS
resistance gates. Those sampled estimates do not bound mesh or out-of-band
error. Finer source meshes and wider bands can require fewer retained modes
or a higher physical clock to fit the explicit resource envelope.

The `bem` route requires a v3 shell or v4 assembly profile with a closed,
outward-oriented exterior skin. It uses the original three-dimensional
displacement fields, triangle normals and areas to compute the full mutual
modal load. Pressure, temperature and mixture changes enter through the same
gas state used by the string bank and cavity; the fitted load reacts on the
same string/body/contact solve. The report records the basis identity, gas,
panel and pole counts, actual work, backward residual, outgoing-flux balance,
and sampled fit errors. The source is the existing symmetric P1/P0 Galerkin
owner: the original transpose work map retains reciprocity before projection.
The unchanged strict matrix and passive-fit gates still apply; no collocation
result is symmetrized. `MAX_PANELS` bounds the source mesh and its cache, without
decimation. `MAX_SWEEP_TERMS` bounds the owner's complete per-frequency work
times the declared frequency count, including singular moments, angular Gram
assembly and dense factorization. The 80-panel two-source 65-frequency fixture
requires 10,535,395,000 units and declares 12 billion.

The optional `shell_angular_refinement` record immediately after `body-radiation
bem` selects the original angular level 0, 1 or 2; omission retains level 0.
Its geometric slab-resonance exclusion and angular-resolution limits refuse
unresolved bands, without falling back to a different solver. These numerical
limits and discrete passivity do not certify acoustic convergence.

`closed-shell-pluck.material-strings` and `closed-shell.body` supply a complete
data-only route through this same shell/BEM owner. After the initial copy and
wire-pack commands above:

```bash
cargo run -p xtask -- matdb-pack \
  --manifest "$demo/closed-shell-manifest.tsv" --out "$demo/closed-shell.fsmatpk"
cargo run --release -p fs-couple --bin material_strings -- \
  "$demo/closed-shell-pluck.material-strings" "$demo/closed-shell-pluck.wav" \
  --body-profile crates/fs-couple/examples/material_strings/closed-shell.body
```

The profile is a manufactured 42-node, 80-triangle closed shell with one
retained free coordinate. Its soft isotropic material and supports keep the
complete coupled fixture within a small spectral envelope. It is an integration
example, not a violin shape or measured wood. Its material pack, physical shape,
string pack, tension and gas mixture remain ordinary input data. The CLI admits
the shell parser's complete 4 MiB profile budget; flat profiles retain their
existing 64 KiB limit.

Open shells, uncapped sound holes, inconsistent face orientations and
multiply shared edges are refused. Coincident assembly boundary nodes are welded
by exact coordinates for the acoustic topology; no gaps are filled or holes
silently capped. This does not certify vertex-manifold quality, absence of
self-intersection or geometric/acoustic convergence. In particular, a violin with f-holes still
needs an explicit coupled interior/exterior opening model; merely selecting
`bem` does not create one. The `skin` route below is one.

The `skin` route computes the same symmetric P1/P0 Galerkin load on the body
profile's declared `acoustic-skin-v1` (below) instead of the shell's own facets:
- **Skin motion:** the skin moves as its anchored points on the parts, checked
  within `MAX_ANCHOR_M`.
- **Angular level:** `LEVEL` is the owner's angular level, 0 to 3.
- **Mouth port:** with `cavity exterior-neck`, the aperture caps carry the
  declared shares of the one neck volume flow. That flow is one more port, on
  the neck's own coordinate.
- **Coincident skin:** a skin that is the shell's own surface, each vertex
  anchored at its node, gives the `bem` load bit for bit.
- **Band:** the owner's slab chart (k times the smallest box width below π)
  caps the band.
  - The load covers the ascending prefix of body modes from `MIN_HZ` to
    `MAX_HZ`. The modes above it carry no exterior load.
  - The report states how many modes are loaded (`body_modes`) and how many
    are left one-way (`one_way_body_modes`). It refuses if a loaded mode lies
    below `MIN_HZ`.
  - With a mouth, `MIN_HZ` must not exceed the neck's resonance with its
    declared ends.
- **Observer:** a `one_way_observer` on the same skin hears this loaded
  mechanics, and its report names the load (`two_way_load`, `two_way_modes`).

The `rayleigh` and `bem` routes require the ordinary exact owner;
`prepared-linear` is refused for them. The `skin` route also runs on
`prepared-linear`, whose air bank carries the load. Without the `shell_observer` selection below, the receiver remains the
existing compact, zero-retardation pressure observation of the **loaded**
acceleration. The Rayleigh example retains one body mode to keep the complete
passive acoustic image small; it is an integration fixture, not a claim of
converged instrument timbre.

## Finite-distance pressure from the loaded material shell

The optional `shell_observer` record follows `body-radiation bem` and precedes
`cavity none`. It sends the original accepted body modal velocities through the
shared causal BEM observation transfer. The same closed shell, signed
mass-normalized displacement fields, explicit mouth-flow traces and physical
gas supply both mechanical reaction and the listener. The existing reciprocal
P1/P0 source supplies its actual mixed layer densities to its adaptive receiver;
panel pressure means are used only for mechanical work. Pressure and all three
gradient components use the same Green field and retained angular resolution.
Material, geometry and gas changes rebuild both from their original data.

```text
shell_observer ID X Y Z CLEARANCE_M PASS_HZ STOP_HZ LATENCY N_FFT FFT... TAIL_TOL DTFT_TOL MAX_FREQ_ENTRIES MAX_SWEEP_TERMS MAX_REPLAY_TERMS
shell_quadrature REL_TOL MAX_DEPTH MAX_KERNEL_EVALUATIONS
```

The optional `shell_quadrature` record follows `shell_observer`. It exposes the
existing adaptive receiver controls; omission uses relative tolerance `1e-7`,
maximum depth 14 and two million kernel evaluations per frequency. This budget
is separate from `MAX_SWEEP_TERMS`, which bounds complete boundary-source work.
Exhausted quadrature refuses without substituting a centroid receiver.

`X Y Z` is the listener position in metres in the physical body frame.
`CLEARANCE_M` must be positive and every surface panel must remain farther from
the listener than that distance. `PASS_HZ` and `STOP_HZ` declare the useful
pass band and the end of the smooth reconstruction taper; the original BEM load
must cover every positive-frequency bin through the taper's upper frequency.
The finest FFT's first positive bin must not lie below the load's lower bound.
`LATENCY` is the explicit causal
reconstruction delay in mechanical samples. Supply at least four successive
doubling power-of-two FFT lengths and bounded transfer and sweep work budgets.
`MAX_REPLAY_TERMS` bounds retained convolution work per accepted mechanical
sample, including pressure and all three gradient outputs. The shared
observation owner checks causal tail energy and complex
DTFT reconstruction error against the declared tolerances.

After the wire and closed-shell packs above are built, the corresponding
integration fixture uses the same installed executable:

```bash
cargo run --release -p fs-couple --bin material_strings -- \
  "$demo/closed-shell-observed.material-strings" "$demo/closed-shell-observed.wav" \
  --body-profile crates/fs-couple/examples/material_strings/closed-shell.body
```

This fixture declares an 8 kHz mechanical/output clock, a 120 Hz pass band,
a taper ending at 400 Hz and 256 samples of reconstruction latency. Its shell,
material and string remain the manufactured integration source described above.
It is a deliberately narrow-band example, not a full-band violin recording.

When selected, finite-distance pressure replaces the compact pressure channel
in the existing scheduled renderer and PCM writer. The observer receives each
accepted body velocity once; callback boundaries do not restart its history.
It adds no second force solver and does not debit radiation energy again.
Physical travel delay, declared reconstruction latency and the authored output
window are retained. The renderer does not normalize pressure, shift the signal
earlier or append an undeclared tail.

Attachment requires a cold performance with the source-derived shell BEM load.
An arbitrary fitted load, Rayleigh plate, unrelated lumped cavity or undeclared
opening cannot substitute for that matching exterior problem. Explicit compact
mouths from the original distributed cavity contribute accepted physical volume
flow alongside structural velocities. An already advanced
performance refuses a new observer because its past boundary motion is unknown.
The report exposes basis and loading identities, gas, panel/mode counts,
transfer controls and work estimates.

This is homogeneous exterior Helmholtz propagation. Gas density and sound speed
affect radiation and travel, but frequency-dependent molecular absorption is
not included. The FIR realization tolerances describe reconstruction of the
computed transfer; they do not certify mesh convergence, material accuracy,
room acoustics or agreement with a measured instrument. Sound holes still need
an explicit coupled interior/exterior opening model.

## Finite-distance pressure through a declared skin, one way

The optional `one_way_observer` record stands where `shell_observer` would and
precedes `cavity`. It sends the accepted body modal velocities, and a lumped
cavity's neck volume flow, through the same shared causal BEM observation
transfer to one receiver, but feeds **no** air reaction back. The mechanics may
therefore be the exact or the `prepared-linear` owner. It refuses alongside
`shell_observer` or a `rayleigh` or `bem` air load. Beside `body-radiation skin`
it hears the mechanics that the same skin's load reacts on.

```text
one_way_observer ID X Y Z CLEARANCE_M PASS_HZ STOP_HZ LATENCY N_FFT FFT... TAIL_TOL DTFT_TOL MAX_FREQ_ENTRIES MAX_REPLAY_TERMS RULE [MAX_PAIRS] MAX_SKIN_PANELS MAX_ANCHOR_M MAX_BOUNDARY_TERMS WORKERS LOSS_FRACTION SHIFT_FRACTION
```

The receiver, band, latency, FFT ladder, tolerances and replay allowance mean
what they mean for `shell_observer`. `RULE` is the panel quadrature,
`triangle-near MAX_PAIRS` or `triangle MAX_PAIRS`. The transfer is baked on one
ruled Burton-Miller operator at every frequency: the geometry-routed owner's
switch from plain CBIE at the skin's box bound is a step in every transfer, and
a bake across it fails its DFT ladder. The mode receipt keeps the routed owner,
plain CBIE below the bound, at each mode's own frequency.
`MAX_BOUNDARY_TERMS` bounds the complete dense work, panels^2 (panels + sources
+ 1) per frequency, over every bin through the taper plus one solve per mode.
`WORKERS` threads the dense assembly, the LU and the per-field solves; no bit
depends on the count. Every bin below the stop must also keep six panels per
wavelength on the skin's largest panel.

The skin is data in the v3 or v4 body profile, after the body and before any
`interior-fe-v1` volume:

```text
acoustic-skin-v1
skin-vertices N
skin-vertex X Y Z K                           (x/y fractions of the body length and width, z in metres, as shell nodes)
skin-anchor W COMPONENT TRIANGLE B0 B1 B2     (K rows: a weighted point in one part's triangle)
skin-triangles M
skin-triangle A B C                           (closed and outward)
skin-apertures P
skin-aperture ID FLOW_SHARE T
skin-cap TRIANGLE                             (T rows)
```

Each vertex moves as the weighted sum of its anchored points on the parts, and
the observer checks every anchored point against the vertex, within
`MAX_ANCHOR_M`, before any boundary work. A panel's normal velocity per mode is
the mean of its vertices'. Aperture caps carry `FLOW_SHARE` of the lumped
cavity's outward neck flow, uniform over the cap, and no structural motion. A
performance with a lumped cavity needs apertures, and one without needs none.
The skin may be the shell's own surface (the integration test anchors the
closed shell's triangles at their own nodes) or a separate surface of the same
body; the generator that writes it states how it was made.

Because nothing is fed back, the report states for every retained mode, at its
own frequency and on the same skin and route, what this owner leaves out:
- the radiation damping ratio and the added-mass frequency shift;
- the mode's wood damping, half-power bandwidth and spacing to the nearest mode;
- the largest first-order air mixing into another mode.

A mode stays one-way when its radiation damping is at most `LOSS_FRACTION` of
its wood damping and its shift at most `SHIFT_FRACTION` of both its bandwidth
and its spacing. The command prints one line per mode and the list of modes
outside that rule, which would have to be fed back. These are estimates on the
declared skin, not a convergence certificate, and they do not change the render.

The transfer and the receipt are the expensive part (hours on a 3000-panel skin).
They depend only on the transfer declaration and the body-mode frequencies. The
declaration's hash binds the skin operator (panels, every source trace, gas,
clearance and rule), the receiver, the clock, the band and the FFT ladder.
`--observer-store DIR` (with `--body-profile`, and only for a performance with
`one_way_observer`) keeps them in an existing directory, in one file per setup
key. The first render bakes and writes the entry, and every later render with
the same key loads it. So notes that differ only in strings, bow, stops or
length share one entry.
- The entry holds the kernels, the ladder, the impedance rows and the plain-CBIE
  limit as exact f64 bits, and a hash over all of them.
- A loaded setup renders the computed one bit for bit.
- A damaged, truncated or foreign entry refuses. Nothing is recomputed behind
  it, and an existing entry is never replaced.
- The key cannot see a change to the boundary solver or the bake that keeps
  their hash domains, so each file name also carries the hash of the executable
  that wrote it: a rebuilt binary bakes afresh.
- The command prints `setup_key`, `setup_loaded` and the entry's path.
- Fill a store with one render before starting others on the same key.

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
infer those properties. A pack that also states `shear_viscosity_constant_loss_reference_omega` [rad/s]
wrote its shear viscosity as `eta_t G / omega_ref` from a measured torsional loss factor `eta_t`. Every twist
mode then gets `eta_t` at its own frequency (the shell's constant-loss spectrum fitted over 50 Hz-40 kHz; a
twist mode outside that band refuses) instead of the Kelvin-Voigt `eta_t omega / omega_ref`. The bending
viscosity stays Kelvin-Voigt. The body additionally needs `poisson_ratio` [1]. Bending
and shear loss must have an explicit valid angular-frequency band; all actual
retained modes, not only the source query frequency, undergo the existing checks.

`packs/violin-string-dominant-{e,a,d,g}-pickering1985-section.fsmatpk` are the four Thomastik Dominant violin
strings measured by Pickering (1985). `xtask matdb-pack` builds them from the sources that
`gen_dominant_string_packs.py` writes into `data/matdb/seed-v1`. Each is an equivalent homogeneous section of the
string's own mean radius (E 0.1555 mm, A 0.3405 mm, D 0.399 mm, G 0.40275 mm; use that radius in the string record):
- its density gives the measured linear density;
- its `young_modulus` gives the least-squares inharmonicity of the printed partials;
- its `shear_modulus` gives the measured torsional stiffness.

Each pack states its torsional loss factor with the constant-loss anchor: 0.0032 (E), 0.071 (A), 0.036 (D) and
0.119 (G). It is the pendulum-and-string decrement with the pendulum's own gravity stiffness taken out. The
bending viscosity is written as 0 because no transverse internal loss is measured; such a string loses
transverse energy only to its air drag and its supports.

This executable selects homogeneous, bonded or close-wound circular strings, fixed tension or
fixed-extension supports and fixed rotational ends. Its shared body is the legacy
simply supported isotropic flat plate or an explicit material plate, shell or
assembly described in [BODY_PROFILES.md](BODY_PROFILES.md), with an optional
uniform-pressure cavity/circular neck. Bridges are interior transverse
connections on the full span **including afterlength**, plus explicitly selected
scalar rotational spring/dashpot joints. A full 3D unilateral bridge-contact
law, coupled f-hole/interior/exterior geometry, distributed cavity modes,
resolved winding contact/slip and solved finger disengagement remain outside
this composition. The computed planar and closed-shell exterior reactions and
the optional finite-distance closed-shell observer are described above.
Without that observer, body and neck pressure retains the existing compact
baffled, zero-retardation approximation. Source-state validity is not measured
instrument-fidelity evidence.

Input caps: 1 MiB performance; 8 MiB per distinct material pack and 64 MiB total
raw packs; 8 strings, at most 16 layers per string, 16 bow contacts, 16 rotational
bridges, 64 force ports, 128 complete-network modes including fluid memory,
16 plate modes; 16,384 events per input family; 262,144 compiled controls per
family; 16,777,216 source/projection work units per compiler; at most 384 kHz,
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
existing PCM output, source/input refusals, and finite-distance shell pressure
from the original accepted body velocities. The fixture TSVs use the existing
`xtask matdb-pack` compiler, not a new importer.
