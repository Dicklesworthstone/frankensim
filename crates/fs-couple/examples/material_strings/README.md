# Material packs to a plucked or bowed instrument

This executable connects the existing normalized material-pack decoder, solid
property resolver, string/body/cavity assembly, optional torsion, canonical bow
schedule, physical load compiler and streamed PCM writer. It contains no modal
frequency presets, separate pluck oscillator, per-string audio mixer or second
physics integrator. All strings see one shared plate and optional cavity.

The basic wire/body coefficients are **synthetic integration fixtures**, not measured
steel, gold, silver, wood, rosin or instrument data. The output is a reduced
physical-model demonstration, not a validated violin/cello/bass recording.

`section measured` also accepts source-bound whole-string line mass, flexural
rigidity, diameter, span and fixed tension through this same performer. The
[measured-string workflow](MEASURED_STRINGS.md) compiles the existing published
Lampis cello-string pack and couples it to a manufactured laboratory plate.
It retains the source's unknown measurement conditions and effective-coefficient
limitations; it does not infer hidden winding, axial, torsional or intrinsic-loss
properties, or claim a complete cello model.

### Shared material cards

Every bulk material path accepts the canonical **material-card pack**
(`FSMCDPK`) used by material reeds, duct walls and other material consumers.
Its original chemistry, phase, process, claim evidence and associated model pack
remain intact throughout string, body, region, fiber, shear-layer, normal-stop and
thermal construction. Ordered bow, bridge, holding and normal-stop interfaces therefore name those
same original material identities. Replacing a material file re-resolves the
original structural laws and their coupled pressure response; it does not
select a different instrument solver.

Legacy claim packs (`FSMATPK`) retain their historical synthetic
`source-pack:<pack_id>` association and exact numerical behavior. The binary
header chooses the decoder, independently of the filename suffix. A malformed
canonical card refuses rather than being reinterpreted as legacy data. The
complete artifact, including any associated models, counts against the existing
source-byte limits.

The `state` row still declares absolute temperature in kelvin and angular
frequency in radians per second. Compatible typed axes retain their quantity
metadata; temperature-difference and cyclic-frequency axes refuse instead of
silently changing those conventions. Missing state coordinates and material
properties still refuse. Associated models are retained as source evidence;
this ingress does not make arbitrary constitutive laws executable or extend
the existing solver's physical domain.
The existing source-derived normal-contact owner also requires matching axis
schemas across both bulk materials and the interface: mixing a bare legacy
axis with a typed axis refuses. Canonical loading preserves that distinction;
it does not rewrite source cards to make their schemas appear equal.

### Held material environment

An optional block immediately after `state` supplies additional source
coordinates for individual bulk packs:

```text
state 293.15 1000
material-states 2
material-state wire.fsmatpk 1 moisture_content 0.12
material-state body.fsmatpk 1 moisture_content 0.08
gas mixture 293.15 101325 mole 2 N2 0.8 O2 0.2
```

Each row gives the original relative pack path, coordinate count, then pairs
of source axis name and value. The example requires packs that explicitly
declare those moisture axes; the supplied constant wire/body fixtures do not,
and refuse that extra coordinate. Up to 64 packs and eight extra coordinates
per pack are supported. Repeated pack paths or axes, unused pack paths, missing
source semantics, incompatible declarations across claims and out-of-domain
states refuse. `T` and `omega` remain reserved for the original `state` row.

A held coordinate must be source-qualified static pressure, composition,
moisture or concentration. Its exact quantity and basis come from the source
card. In particular, the solid's declared moisture content is a supplied
material state; the renderer does not infer it from atmospheric relative
humidity. Changing that coordinate re-resolves the original density, moduli
and supported loss coefficients. Source-declared curves and surfaces therefore
alter the existing string, section, fiber and isotropic plate laws, and their
joint pressure output, without replacing the immutable material card or
choosing a different instrument model. The resolved state identity changes
with the coordinates.

The existing support prescription remains authoritative: maintained force stays
fixed, while fixed extension changes its resulting force with the resolved
axial stiffness. There is no implicit pitch correction. Temperature re-queries
for explicit contact-thermal profiles retain the same extra coordinates;
thermal expansion still needs its existing qualified source law and reference
state. This input does not supply a sorption, diffusion, swelling or evolving
composition law. A dynamic or otherwise unsupported constitutive dependency
still refuses, and every retained frequency must remain in the source's
applicability band.

`section measured` keeps its own explicit observation axes, including unknown
measurement conditions. A bulk `material-state` binding for such a string
refuses instead of overriding those observations. Existing performances with
no `material-states` block keep their original query and behavior.

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

`--bow-trace NEW_FILE` also writes the bow's slip history. The file has one
tab-separated row per accepted sample and bow port, with these columns:
- the physical sample;
- the port;
- its regime (`Sticking`, `Sliding`, `Released` or a layer's);
- the slip velocity in m/s, which is the step's slip distance times the rate;
- the traction and the normal force in N.

With the flag, the render runs one sample per callback, and the WAV bytes are
unchanged. Like the WAV, the file must not exist already, and it is created
only after admission.

### Move the bow position along the string

`migrating-bow.material-strings` moves a finite-width bow across the original
vibrating string/body/cavity system and holds the last position. Run it with
the same wire/body packs and command used above for `pluck-and-bow`.
Insert this optional block after any `bow_ribbon` / `bow_fibers` declarations
and before `forces`:

```text
bow_station_motion 2 0.00003 1048576 3
bow_station_knot 480 0.2
bow_station_knot 2400 0.28
bow_station_knot 3600 0.23
```

The first row declares maximum travel speed along the string [m/s], maximum
travel per mechanical sample [m], geometry-update work limit, and knot count.
Each knot gives an elapsed mechanical sample and one position fraction per
original bow track. With a ribbon, the single position is its centre; the
existing ribbon geometry places every cell and checks the full footprint.
Positions interpolate linearly from the initial schedule position and hold
after the final knot. Bow speed and load remain in the canonical gesture
schedule; the position path does not retune or restart the string.

The same force/velocity contact maps carry transverse and torsional motion,
body response, material memory and microphone pressure. Independent moving
finger stops may share the physical clock. A `bow-contact` normal follows its
bow automatically, including material fiber receivers: the same solved normal
reaction supplies friction while the normal patch's displacement-map motion
contributes its physical actuator work to the shared energy balance. Bow
speed remains the tangential surface speed; longitudinal station travel is a
separate geometry control. The `bow_fiber_motion` control described below also
moves the crossing along the hair and retains its separate normal-patch
restriction. Retained shear/thermal films still refuse because
their moving material histories are not transported. This is a midpoint geometry
discretization: refine time, retained modes and ribbon cells before making
accuracy claims. It does not establish measured instrument realism.

### Native application source packets and output clock

The same performer can run from an explicit set of attached source files:

```bash
material_strings performance.material-strings output.wav \
  --self-contained \
  --attached-source materials/wire.fsmatpk \
  --attached-source materials/body.fsmatpk
```

Names inside the performance and optional body, stop, contact-thermal, interface,
and gas sources must resolve to those exact attachment paths. Add each referenced
file with another `--attached-source`. Restricted execution snapshots the
admitted files before decoding them; an unlisted reference refuses rather than
reading another file. It also refuses `--observer-store`, so computed observer
setups come from the attached physical sources.

Without `--self-contained`, existing CLI path resolution is unchanged: material
references are relative to the performance file, while profile and explicit gas
source options use paths supplied by the caller.

For mechanics above 48 kHz, `--decimate` writes 48 kHz PCM using the existing shared
pressure observer. The mechanical state, material histories, gestures, and
simulation time step remain at the declared input rate. The rate must be a
supported integer multiple of 48 kHz, and the mechanical horizon must contain
complete output intervals. The observer retains its causal filter delay and
does not append or discard mechanical samples; its observation policy is printed
with the render diagnostics.

`--decimate` composes with `--contact-thermal` and `--bow-trace`. Contact heat
continues through the original coupled mechanical/enthalpy owner at every
mechanical step. The trace observes each accepted mechanical sample before
pressure decimation, so its indices and row count retain the mechanical clock.
Only the microphone pressure is resampled. Trace I/O failure refuses the output
without pretending that already accepted mechanical steps were rolled back.

`--workers N` (1 to 256; default 1) declares threads for the dense BEM air
loads computed before the render. A `body-radiation skin` load solves its
frequencies on them, and dense factorizations share them. It is a scheduling
budget only: every load, and so every output sample, is bit-identical for every
count. Each skin frequency in flight holds its own dense operator and factor,
about 32 n^2 bytes for n reciprocal unknowns (0.37 GB at n = 3389).

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

Material shear layers compose with either `geometric-stretching` or
`longitudinal-waves`. For example, insert
`geometric-stretching 0.2 64 100000 1e-11 1e-10` before `cavity` in
`shear-layer-bow.material-strings`; alternatively select
`longitudinal-waves 0.2 0.01 64 100000 1e-11 1e-10 2 1 2`.
The original strain energy, layer elastic/plastic history and any series-drive
spring/dashpot remain in one staged contact solve and complete energy balance.
Changing only layer source data changes the shared pressure output. These two
strain selections remain mutually exclusive. Contact-thermal continuation uses
that same nonlinear solve at each proposed thermal endpoint, retaining accepted
string, film and drive histories until the complete balance is accepted.
The film and declared receiving cells evolve thermally; the string/body law
remains at its declared source state. This contact-only continuation does not
infer bulk string heating or thermal expansion. Fixed layer geometry, retained
compression and the existing residual-compliance exclusion remain required.

The `state` row supplies the initial temperature and material query. Without
an explicit contact-thermal profile it stays frozen. Every
selected property needs explicit temperature support and applicability from
static storage through the sample clock's Nyquist frequency. Constants and
temperature-only curves are supported. Pressure-, rate-, strain- or
frequency-dependent coefficients require a more general constitutive owner
and refuse here. Oscillatory viscosity and friction coefficients cannot
substitute for plastic shear viscosity.

An optional source-memory row immediately after each `shear-layer` adds the
material's retained Maxwell shear arms to that same contact:

```text
shear-layer bow coating.fsmatpk 0.000001 0.0001 0.05
shear-memory 2 shear_delta_g_1 shear_tau_1 shear_delta_g_2 shear_tau_2
```

Each pair names measured or explicitly manufactured `DeltaG` [Pa] and `tau`
[s] properties in the original pack. `shear_modulus` is the equilibrium
modulus, and each arm stiffness is `area * DeltaG / thickness`. Ribbon cells
partition the area once and retain separate shear histories. Arms relax even
below plastic yield and at held total slip; their spring energy and physical
dashpot heat enter the original simultaneous string/fiber/body solve and
pressure observer. No modal damping ratio or audio filter supplies this loss.
The source must provide positive values and DC-through-Nyquist applicability.
The physical clock must resolve every arm with `dt <= 2*tau`; increasing
sample rate supplies the temporal refinement. There are at most 64 arms per
contact, subject also to the original resolved-property budget.

The existing `--contact-thermal` route deposits both plastic and Maxwell heat
into each original film cell. It retains every arm on temperature refresh,
including temperature-dependent relaxation times. Arm moduli must be constant
over the full thermal interval; affine equilibrium modulus retains the existing
mechanical-entropy account. A temperature-dependent arm modulus refuses because
its additional entropy is not supplied by this rung. These source declarations
describe a constitutive material model; they do not identify real rosin data or
remove the fixed-film, retained-compression and geometry limitations below.

This first layer model requires retained positive compression, fixed layer
geometry/tangent frame, and an elastic shear limit in `(0, 0.1]`. Prescribed
bow loads must stay positive; contact-loaded bows instead author zero extra
load and use the simultaneously solved normal reaction. Opening refuses
atomically, preserving stored layer energy and accepted history. Layer
opening/recontact, wear, migration and nonlinear string geometric storage
are not yet composed. Explicit finite contact heating is available through
`--contact-thermal PROFILE`; see [CONTACT_THERMAL.md](CONTACT_THERMAL.md).
That route uses the same source cards, canonical schedule and pressure/PCM
owner, with source heat capacity, density and conductivity, retained strain,
and material temperature feedback. No thermal constants are inferred. Prepared-linear propagation is
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

The numeric tension field specifies the maintained force of the **entire physical
bundle**, independently of its material stiffness. To keep the same tightened
support spacing while changing fiber material, replace that field with one of
these explicit boundary conditions:

```text
bow_fibers PACK RADIUS_M LENGTH_M FIBER_COUNT fixed-extension L0_M STRAIN_LIMIT CONTACT_FRACTION AXIAL_MODES NORMAL_MODES
bow_fibers PACK RADIUS_M LENGTH_M FIBER_COUNT thermal-extension L0_M STRAIN_LIMIT REFERENCE_K CONTACT_FRACTION AXIAL_MODES NORMAL_MODES
```

`LENGTH_M` is the current distance between the rigid fiber supports. `L0_M`
is each fiber's stress-free length at the current material state for
`fixed-extension`, or at `REFERENCE_K` for `thermal-extension`. The original
section/support owner computes each fiber's tensile force from the same
resolved `young_modulus` and cross-section used for its waves, then multiplies
by the physical fiber count exactly once. Refining ribbon cells only partitions
that bundle. No nominal force field is silently overridden and no pitch is
retuned. The original source pack, query, support prescription and reference
temperature remain retained by the native input and material catalog.

The thermal form requires the fiber pack's
`linear_thermal_expansion_coefficient` over the complete path from
`REFERENCE_K` to the material `state` temperature. The shared material owner
integrates every crossed temperature knot; endpoint alpha times temperature
change is not substituted. It checks total, thermal and elastic strains against
`STRAIN_LIMIT` in `(0,1)` and refuses missing data, unsupported paths and a
slack or compressed resultant. This is frozen-state, small-strain initialization
with explicitly supplied current geometry: it does not predict bow-stick
expansion, changing support separation, live thermal evolution or Poisson
contraction. `STRAIN_LIMIT` is an admission assumption, not a measured yield
limit. Optional `bow_fiber_strain` retains its separate additional-motion limits.

For example, the existing manufactured `fiber-contact.material-strings`
fixture can replace `180 40 0.4` with
`180 fixed-extension 0.594 0.05 0.4`, using its unchanged fiber pack and
`fiber-contact.stops`. With the supplied 0.6 m span, 40 micrometre radius and
3 GPa modulus, that boundary condition gives about 27.42 N for 180 fibers.
Changing the fiber modulus now changes installed tension as well as EA/EI.
With the normal-contact profile, the resulting fiber deflection changes the
solved bow load, traction and the existing string/body pressure output.
The plain prescribed-load ribbon does not infer this normal-contact geometry.

An independently characterized isotropic Young relaxation spectrum can be
attached immediately after `bow_fibers`:

```text
bow_fiber_memory MAX_STATES EQUILIBRIUM_MODULUS_PROPERTY ARMS DELTA_E_PROPERTY TAU_PROPERTY ...
```

Every property is resolved from the **same fiber pack and state**. Moduli have
units Pa and times have units s. The equilibrium modulus must equal the pack's
`young_modulus`. The original section owner retains each modulus increment and
relaxation time; the fiber owner projects that one law into axial extension and
bending, and the original simultaneous contact solver advances its histories.
Changing these source data changes traction, slip, string/body response and the
emitted pressure. No spectrum is estimated from the requested sound.

`MAX_STATES` bounds the complete ribbon history (cells × retained modes ×
positive material arms), with a ceiling of 256; it does not cap or discard modes.
The source must explicitly support the transient frequency range and the
physical clock must resolve the selected relaxation times. Fiber memory currently
requires the numeric maintained-force support form; fixed-extension and thermal
preloads refuse because the original owner does not initialize their hereditary
prestrain histories. Kelvin–Voigt
viscosities remain explicitly sourced parallel mechanisms; they are not inferred
from the relaxation spectrum. This linear-strain law cannot yet be combined
with `bow_fiber_strain`, whose nonlinear strain history requires a different
constitutive coupling. Missing or incompatible source data refuse before WAV
creation. The measured horsehair pack cited below does not supply this spectrum;
the executable tests use plainly identified manufactured material data.

For an accelerating or reversing stroke, explicitly add this record after
`bow_fibers` and optional `bow_fiber_memory` or `bow_fiber_strain`:

```text
bow_fiber_motion MAX_SPEED_M_S MAX_ACCELERATION_M_S2 MAX_STEP_M MAX_UPDATE_TERMS
```

The original bow track now supplies endpoint support velocities. A control
tick activates its target at the unchanged schedule boundary; the support
reaches that target at the **following mechanical endpoint**, one physical
interval later. Initial velocity is already co-moving. On each mechanical
interval the renderer derives the uniform-support acceleration,
applies its mass-projected inertial load to the actual axial fiber modes,
integrates support travel, and moves the contact along the original fiber
material coordinates. The existing joint solver still owns string/fiber/body
waves, traction, loss, pressure and PCM output. Stroke endpoints and reversals
must remain inside the original hair span; excessive acceleration or map work
refuses before a physical step commits. Low-rate velocity steps represent
finite-interval accelerations, so refine the control rate together with the
mechanical rate. This is not a reset or a change to string frequency.

The original dry-interface material domain is checked against both actual
endpoint support speeds and both endpoint maps, as well as the work-conjugate
average slip. String and fiber motion both contribute; a mean speed cannot
hide an out-of-domain endpoint. These sampled checks do not certify unsampled
continuous-time extrema. Endpoint columns reserve `2 * modes * ports` additional
terms within the same explicit map-work budget.

This motion rung accepts prescribed loads or the solved reactions of canonical
normal fiber contacts, with dry friction. The normal receiver follows the same
moving material crossing as the tangent. The original unilateral contact law
resolves indentation and release, so separation releases friction automatically;
there is no independent force envelope for a contact-loaded bow. Arbitrary
authored normal receiver rows and retained rosin/thermal patches still refuse.
The native mechanical ledger
is in the accelerating support frame, with inertial-force work included in
external work. `FiberSupportFrame` also reports laboratory mechanical energy
and physical support work. It adds the rigid-translation kinetic energy and
the cross term between support velocity and the original fibers' modal momentum;
coordinate inertia work is never classified as frictional heat. Normal work
from migration along the fibers is included in physical support work and exposed
as `normal_attachment_work_j`. That share is already in native `obstacle_work_j`:
subtract it when combining the native accounts with laboratory support work.
String-side contact migration and normal indentation keep their separate work.

For a continuous acceleration or reversal, use this alternative record followed
by velocity endpoints on the original mechanical sample clock:

```text
bow_fiber_program INITIAL_SPEED_M_S MAX_SPEED_M_S MAX_ACCELERATION_M_S2 MAX_ABS_TRANSLATION_M MAX_UPDATE_TERMS KNOTS
fiber_speed SAMPLE VELOCITY_M_S
```

For example, a 48 kHz performance of at least 384 samples can use:

```text
bow_fiber_program 0 0.1 20 0.01 1048576 3
fiber_speed 128 0.02
fiber_speed 256 -0.02
fiber_speed 384 0
```

Velocity varies linearly between declared endpoints, and holds its final value
after the last knot. Sample 1 is the end of the first physical interval. The
original bow track must have zero tangential speeds. It prescribes normal load
only for unbound contacts; contact-loaded bows require zero authored load and
use their solved normal reaction. The program supplies endpoint velocities to the **same** support
inertia and moving-contact owner as `bow_fiber_motion`; choosing both refuses.
Velocity, acceleration, translation at reversal extrema, and held final travel
are checked before playback. The finite material hair span is also enforced by
the original contact owner. Callback block size does not change this motion.

To press and withdraw the same fiber supports, add a normal velocity program
after `bow_fiber_motion` or `bow_fiber_program` and its endpoints:

```text
bow_fiber_normal_program INITIAL_SPEED_M_S MAX_SPEED_M_S MAX_ACCELERATION_M_S2 MAX_ABS_TRANSLATION_M MAX_UPDATE_TERMS KNOTS
fiber_normal_speed SAMPLE VELOCITY_M_S
```

For example, this 48 kHz program presses toward an initially separated string
and returns to its original clearance over 256 mechanical intervals:

```text
bow_fiber_motion 1 1000 0.001 1048576
bow_fiber_normal_program 0 0.1 20 0.01 1048576 4
fiber_normal_speed 64 -0.02
fiber_normal_speed 128 0
fiber_normal_speed 192 0.02
fiber_normal_speed 256 0
```

The tangential support velocity still comes from the original bow track in
this example. Normal velocity uses the same linear interpolation, finite
horizon and acceleration checks as `fiber_speed`. Its integral translates the
original contact clearance, while its acceleration supplies the mass-projected
inertial load to each cell's original transverse modes. Consequently changing
fiber material, count, tension or geometry changes the resulting normal force
and traction in the existing string/body/acoustic solve. The program supplies
no normal-force envelope. Contact load remains zero during separation, and
pressing work remains mechanical work rather than frictional heat. Independent
stop-motion records must keep zero offset placeholders at these fiber contacts;
other stops can retain their own motion.

This is uniform translation of rigid supports along one common contact normal;
it does not model bow-stick flexibility or a rotating normal frame. The material
and dry-contact restrictions of the axial motion rung still apply. The report's
`normal` field exposes normal travel, velocity, inertia and contact translation
work. The top-level `FiberSupportFrame` energy, support work and inertial work
include both translation directions. When replacing native support-frame work
with laboratory support work, subtract `normal.contact_translation_work_j` from
native obstacle work as well as `normal_attachment_work_j`; both shares are
already included in the complete physical support work. No optional trajectory
changes the physical sample clock or PCM observation path.

Omitting the axial motion records retains the explicit constant-speed,
frozen-station model used by the fixtures below. These motion selections do not
establish measured-instrument sound quality; material, geometry and resolution
still require validation.

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

Both fixtures begin with fixed interior contact stations and rigid fiber endpoints
in a uniformly translating tangential support frame. Adding the optional
`bow_station_motion` block above to `fiber-contact.material-strings` moves each
string-side normal and tangential contact together at the original ribbon width;
the same `fiber-contact.stops` profile and zero-load gesture schedule apply.
The contact coordinate along each fiber and its endpoint frame remain fixed,
so this is longitudinal travel on the string, not travel along the bow hair.
Normal offset trajectories
on fiber receivers refuse: moving the endpoints would also require their
frame-acceleration forces, which this frozen model does not supply. Ordinary
fixed/bridge stop motion remains available. The original `fiber-ribbon`
fixture retains prescribed uniform normal load.
The same fiber continuum now composes with `geometric-stretching` or
`longitudinal-waves`. Its modes append after the original string/body/fluid
coordinates without replacing their nonlinear strain storage, layered material
histories or pressure observer. Equal and opposite fiber/string tractions and
any resolved normal reactions participate in that same constitutive solve.
For example, add `geometric-stretching 0.2 64 100000 1e-11 1e-10` before the
`cavity` row in `fiber-contact.material-strings` to retain the source-derived
string stretching energy with its existing fiber contact profile. The
longitudinal selection instead retains the declared axial wave coordinates;
the two string strain images remain mutually exclusive.

To also couple the fibers' own normal deflection to axial strain and tension,
add this optional record immediately after `bow_fibers`:

```text
bow_fiber_strain MAX_SLOPE MAX_STRAIN MAX_ITERS MAX_FORCE ABS_TOL REL_TOL
```

The first two values declare the motion envelope (slope at most 0.3 and strain
at most 0.1); the remaining four declare the joint numerical solve. They must
match any `geometric-stretching` or `longitudinal-waves` solve configuration.
These are admission and numerical limits, not fitted material measurements.
The source pack still supplies all stiffness, mass and loss coefficients.
Normal deflection now changes the same retained fibers' axial waves and tension
through `u_x + y_x²/2`. The unchanged contact solve propagates their reactions
to string motion, body pressure and WAV output. No second fiber oscillator is
introduced. Omission keeps the explicitly linear fiber rung.

`fiber-strain-contact.material-strings` enables both fiber and string strain
with the existing manufactured materials, geometry and normal-contact profile:

```bash
cargo run --release -p fs-couple --example material_strings -- \
  "$demo/fiber-strain-contact.material-strings" "$demo/fiber-strain-contact.wav" \
  --stop-profile "$demo/fiber-contact.stops"
```

This remains a physical integration fixture, not a measured violin or horsehair
model. Retained axial harmonics must be refined alongside normal modes; omitted
axial harmonics equilibrate statically. Excess slope/strain, local loss of
tension, unsupported resolution, incompatible solve settings and a failed
joint solve refuse. Force-loaded initialization without normal obstacles or
retained shear layers uses the same nonlinear preload owner described below.
Viscosity remains the original
small-strain operator; finite-slope viscous coupling is not supplied.

Without `bow_fiber_motion` or `bow_fiber_program`, changing support speed still
refuses. The optional motion controls supply frame-acceleration forces and the
moving material crossing described above. Neither path models a flexible bow
stick, inter-hair rubbing or transported thermal rosin histories. Axial and normal mode
counts, ribbon cells and physical sample rate are
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

Bonded concentric sections use this same exclusive selection. Each physical
layer resolves the declared property roles on its own material pack: it supplies
its own equilibrium modulus, excess moduli and time constants. The original
section integrates each `DeltaE * I_layer`; each time survives in the exact
solver. Neither a common time constant nor a scalar viscosity is fitted.
Two transverse planes use the same physical section with separate strain
histories; torsion remains an independent sourced channel. All layers together
admit at most 128 source arms, and all declared string/body histories share the
existing 256-state limit and setup budget.

`relaxing-layered-strings.material-strings` uses that same synthetic pack in
the core and sleeve of two bonded strings, with two transverse planes per
string. It can be rendered with the command above after changing the input and
output filenames. Assigning a different pack to a `layer` row changes only
that layer's density, stiffness and relaxation spectrum. The finished radius,
gas loading, support prescription and shared bridge/body/cavity mechanics stay
on their original construction path. This fixture is a bonded-section example;
it does not assert that a sliding helical winding is a bonded annulus.

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
strings, bonded sections and close-wound cores on the exact solver;
`prepared-linear` refuses constitutive histories. Core relaxation in a wound string keeps the original helix mass,
finished gas-loading radius and wrapping material identity. All participating
source bands must cover the complete instantaneous coupled frequency bound. These models
supply bending relaxation; the torsional selection below supplies independent
shear histories. Axial relaxation, interlayer slip and measured wound-string
behavior remain outside their scope. All supplied values
are synthetic integration data, not measured instrument properties.

### Tension loss

A string's construction (its core, winding and their contacts) loses a fraction
eta_F of the tension's share of each mode's stiffness. This is Valette and
Cuesta's eta_F term: Q_F^-1 = eta_F T k² / (T k² + EI k⁴). Each mode carries it
as viscous damping, zeta_F = eta_F T k² / (2 mu omega²), beside the pack's
bending law and the gas drag, and it is booked as material heat. It changes no
mode frequency.

The optional per-string record selects one row of
[`tension-loss.tsv`](tension-loss.tsv) by its id:

```text
tension-loss lampis-2024-dominant-cello-g
```

Place it after the optional `bending-loss` row, or in that row's place.

The table holds cited constructions, not the strings in these examples:
- Lampis's Dominant cello G, 9.0e-4. It is inverted from that string's mean Q
  of 1064 on its own data, so it is an upper bound.
- Mansour, Woodhouse and Scavone's Kaplan cello strings, 1.2e-4 to 2.3e-4.
- Galluzzo, Woodhouse and Mansour's steel violin E, 5e-5.

Each row's source text says what was measured and what was left out. No violin
string's eta_F is measured.

The receipt prints:
- the selected row with its source, and lists the others;
- with no row selected, `eta_F=stated unknown, entered as 0`, and no tension
  loss term is formed.

An unknown id refuses before material loading, as does a repeated record.

In the library:
- an authored modal ratio or Rayleigh pair refuses beside any eta_F, because
  each is already a complete reduced loss;
- a layered thermoelastic string refuses a nonzero eta_F, because none of its
  layers owns that heat.

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
remain binding. By default that same frozen gas supplies string drag, cavity
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
performance, with each original regional source reaching its own consumers.

### Independent cavity and exterior atmospheres

An optional `cavity_gas` row immediately after `gas` and before `body` declares
the interior medium with the same native gas grammar. Omission inherits `gas`.
For example, a helium-filled body in the original exterior atmosphere is:

```text
gas mixture 293.15 101325 mole 2 N2 0.8 O2 0.2
cavity_gas mixture 293.15 101325 mole 1 He 1
```

The existing spatial cavity FE, volume-compression springs, neck inertias and
optional thermal-diffusion histories use `cavity_gas`. String drag, the complete
exterior acoustic load and microphone pressure still use `gas`. Internal gas
changes therefore alter the body's actual reciprocal pressure loading before
sound observation. Both records retain their own species/source receipts;
`cavity_gas thermochemical ...` uses the same attached source registry.

The two reference pressures must agree to roundoff. A mean pressure difference
requires a static structural load/prestress solve, which this acoustic path
does not perform. Temperature and composition can differ as frozen regional
states; mean heat exchange and species diffusion are outside this model.
Open-mouth regions additionally retain the original compact acoustic limits.

A distinct lumped interior is supported only by `cavity exterior-neck` with
its actual `body-radiation skin` load and aperture geometry. The internal neck
keeps its original inner-end mass and frozen Bessel loss; the exterior load
owns the outer inertia and radiation in the exterior gas. The existing neck
resistance still uses its interior-gas reference frequency, not a broadband
neck law. Mixed-media `neck` and `physical-neck` declarations refuse because
their compact outer-end laws include inertia/loss from the one supplied gas.
An explicit `cavity_gas` without any lumped or spatial cavity also refuses.

`regional-cavity.material-strings` uses the original manufactured
`closed-shell-cavity.body`, `closed-shell.fsmatpk` and `wire.fsmatpk`:

```bash
cargo run --release -p fs-couple --example material_strings -- \
  "$demo/regional-cavity.material-strings" "$demo/regional-cavity.wav" \
  --body-profile "$demo/closed-shell-cavity.body"
```

Build the body pack from `closed-shell-manifest.tsv` with the existing
`xtask matdb-pack` command and copy that original body profile alongside the
performance. This closed-shell specimen and compact body pressure remain
integration examples, not measured violin/cello/bass or free-field validation.

### Molecular caloric response inside the material body

An optional `cavity_molecular` row after `gas` and optional `cavity_gas`, before
`body`, binds molecular reservoirs to the existing spatial thermal cavity:

```text
cavity_molecular MAX_CALORIC_COORDINATES MAX_MEMORY_STATES LAW
```

The body profile must declare its original `interior-fe-v1` geometry and
`air-thermal isothermal-walls ...` heat eigenwindow, as in
`closed-shell-thermal.body`. `LAW` uses the shared material relaxation grammar:
`iso9613-terrestrial-air`, subject to its composition/domain restrictions, or
`declared N` followed by explicit species, reservoirs and collision rates.
The law is evaluated against `cavity_gas`; omission of that record inherits the
original `gas` record. Thermochemical interiors use the same attached species,
NASA-9 and transport registry. Exterior observer absorption is a separate
physical source and never supplies the cavity's molecular data.

For example, these **synthetic integration coefficients** give nitrogen one
rotational reservoir with a 25 ms reference time at 101325 Pa and explicitly
declare oxygen non-relaxing:

```text
gas mixture 293.15 101325 mole 1 He 1
cavity_gas mixture 293.15 101325 mole 2 N2 0.8 O2 0.2
cavity_molecular 2 96 declared 2 N2 1 rotation 2 power 2533.125 293.15 0 280 320 O2 0
```

Changing the interior composition, pressure, temperature or authored collision
coefficients changes the joint thermal/molecular bulk response, reciprocal
body/string force and pressure through the original simulation owner. Every
present species needs an explicit law entry; zero reservoirs must be declared,
and legacy air records without species identities refuse. `rotation DOF` and
`vibration MODES THETA_K DEGENERACY ...` use the existing reservoir laws;
`power P_TAU_REF T_REF EXPONENT T_MIN T_MAX` and
`partners N ID LOG_P_TAU SLOPE T_MIN T_MAX ...` use their existing rate laws.

`MAX_CALORIC_COORDINATES` reserves gas temperature plus active reservoirs per
heat mode (1–128). `MAX_MEMORY_STATES` reserves the complete joint heat and
unresolved molecular histories (1–4096, also subject to the frontend's existing
complete-state and setup-work ceilings). This replaces the heat-only memory
allowance, so thermal energy is not counted twice. The original thermal FE,
projection, cancellation and native propagation limits remain binding.
Missing thermal geometry, a separate lumped cavity and `prepared-linear`
propagation refuse. Both mean gases remain frozen during the performance;
these linear acoustic histories do not evolve mean heat or species inventories.

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

Selected claims require explicit finite bare-SI validity on `T`, `p` and `v`.
Constant coefficient sources retain the original path. A `coulomb` declaration
may select coefficients that depend on `T` when every affected contact has an
explicit finite-cell dry route in a version 2 `--contact-thermal` profile; see
[CONTACT_THERMAL.md](CONTACT_THERMAL.md). The existing heat solver then solves
interface temperature, source-refreshed traction and the original cell energies
together. Coefficient response to pressure, speed or additional history is not
supplied by this temperature-only route. Stribeck and hyperbolic coefficient
sources remain constant throughout their applicability.

The query chooses a supported immutable source; its pressure and speed are not
imposed on the physical solver. Actual normal load divided by nominal area and
actual relative surface velocity, including torsion, enter the original initial,
step-average and endpoint checks. Missing, conflicting, unsupported or
out-of-domain sources refuse. A released bow does not require positive-pressure
friction data at zero load.

Nominal area and declared surface history remain explicit inputs. The original
packs, selected claims and actual coefficient-query receipts survive with the
performance. This binding does not infer a Hertz patch, missing hair or rosin
properties, a transported rosin layer, or wear; loading a source does not promote
it to experimentally validated instrument data.

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

### Material-derived normal finger/support contacts

The same `--stop-profile` also accepts `material-stop` records. These resolve the
actual outer string material and an explicitly named receiver material through
the existing Hertz/Hunt--Crossley point-contact owner. Changing either bulk pack
changes its reduced modulus and normal reaction in the original coupled solve.
There are no independent stiffness, exponent or contact-loss overrides:

```text
material-stop ID STRING STATION SIDE fixed GAP_M RECEIVER_BULK_PACK INTERFACE_PACK MODEL_PACK KAPPA_MAX_PER_M KAPPA_MIN_PER_M MAX_RATE_M_S HALFSPACE_DEPTH_M LAYER_THICKNESS_M RADIUS_REL_UNCERT MODULUS_REL_UNCERT LOAD_REL_UNCERT MAX_FORCE_N MAX_INDENT_M
```

The record belongs under the existing `frankensim-string-stops-v1` header.
`SIDE` is `positive`, `negative`, or `normal NX NY`; station and optional
normal/station-motion programs use the existing physical controls. `GAP_M`
is signed unloaded surface clearance, so negative values explicitly preload
the contact. `fixed` selects a prescribed receiver. A moving receiver can instead
be an actual node of the retained material shell or shell assembly:

```text
shell-receiver ID SHELL_PART NODE NX NY NZ RX_M RY_M RZ_M
material-stop ID STRING STATION SIDE shell:ID GAP_M @receiver INTERFACE_PACK MODEL_PACK KAPPA_MAX_PER_M KAPPA_MIN_PER_M MAX_RATE_M_S HALFSPACE_DEPTH_M LAYER_THICKNESS_M RADIUS_REL_UNCERT MODULUS_REL_UNCERT LOAD_REL_UNCERT MAX_FORCE_N MAX_INDENT_M
```

Declare each named receiver before its contacts and before the motion blocks.
`SHELL_PART` is the original assembly component index, or zero for a single
shell. `NODE` addresses that part's original physical mesh. The finite unit
vector `N` and node-to-contact arm `R` are in the body's global XYZ frame.
The exact map is `N · (U + Theta × R)`, so the reaction includes both force
and the moment `R × (force N)`. The string-side normal and this body direction
must describe the same physical direction in their respective declared frames;
the adapter does not infer a string-to-body orientation. `negative` reverses
both sides together. An authored `stop` can also use `shell:ID` with its existing
explicit constitutive coefficients.

For a `material-stop`, `@receiver` is mandatory with `shell:ID`: the actual
incident body material supplies the original immutable bulk card, physical
state and elastic source receipts. A separate receiver pack cannot override
them. All elements touching the node must have the same resolved material
state; interfaces between different materials require their own physical patch
description and refuse here. Assembly parts share the original assembled body
basis while retaining their own material identities. Body motion reacts in
the same simultaneous normal/friction solve as string, bridge and cavity
motion; no new body mass, delayed force or prescribed-motion surrogate is added.

`INTERFACE_PACK` is an existing `NormalizedInterfacePack`, ordered string
surface A then receiver surface B. `MODEL_PACK` is an existing
`NormalizedModelPack` with its complete normalization receipts and exactly
one supported normal-law selection. The record explicitly associates these
artifacts; the interface, model and source domains must agree. The common
material state uses the performance's temperature and angular frequency.
The original point owner requires density, Young's modulus, Poisson ratio and
yield stress from both bulk cards, explicit adhesion data from the dry
interface, and the selected model's applicability parameters. Missing data,
unsupported laws, incompatible ordered materials and adhesive domains refuse.

Both positive principal curvatures describe the **relative local gap**; they
are geometry data, with the larger first. Equal values describe a circular
patch; a zero curvature/line-contact geometry is unsupported by this input.
Depth, layer thickness, uncertainty and rate bounds must be explicitly supplied.
The declared depth and layer thickness cannot exceed the actual contacting
string constituent: radius for a homogeneous string, winding-wire radius for
a close-wound string, or outer-sleeve thickness for a bonded string. A shell
receiver additionally bounds both depths by the smallest actual thickness of
its incident sections. The current point law requires the receiver's actual
isotropic elastic body state; an orthotropic-only body does not become an
isotropic contact merely because its pack also contains scalar modulus data. The
original finite-patch owner additionally checks patch size, strain, pressure,
yield and rate applicability. This does not infer a tissue model, treat a whole
winding as homogeneous bulk, or certify a measured instrument.

Complete bulk/interface/model receipts and the admitted operating box stay
attached to the runtime normal contact. Each completed mechanical/thermal
candidate checks both contact indentations, the actual discrete closure rate
(including changing station maps and support offsets), and the mechanical
sample duration before any state, work report or clock is published. Solver
scratch iterates remain unpublished. The source state and relative curvature
are frozen during playback; changing material temperature/history or contact
curvature requires a newly admitted source, not a hidden coefficient update.
All source paths use the existing bounded native attachment/file loader.

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
and optional `thermal-supports`, before any `translational-bridges` and `cavity` records:

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

## Independent translational bridges from physical shell nodes

All three input versions accept an optional block after all string/section/
polarization records, optional `thermal-supports` and optional
`rotational-bridges`, before `cylinder-gas-loading`, `prepared-linear` and
`cavity`:

```text
translational-bridges COUNT
translational-bridge STRING STATION SX SY SHELL_COMPONENT NODE GX GY GZ K_N_PER_M C_NS_PER_M REST_M
```

Each row adds a centreline spring/dashpot between one interior string station
and one physical node of a **v3 shell or v4 assembled-shell `--body-profile`**.
`STRING`, `SHELL_COMPONENT` and `NODE` are zero-based indices; a standalone
shell has only component zero. `STATION` is a fraction strictly between zero
and one of the full supported string, including afterlength. `[SX, SY]` is a
signed unit direction in that string's transverse frame; `[GX, GY, GZ]` is a
signed unit direction in the body's global frame. A single-plane string
requires `SY = 0`. Stiffness and viscous resistance must be finite,
nonnegative and not both zero; the stress-free extension is signed metres.

The body displacement row comes from the selected node's original global XYZ
translations in every retained shell `physical_modes` vector, preserving mass
normalization and mode signs. Assembled components share the same assembled
eigenbasis: `SHELL_COMPONENT` selects physical geometry, while the modal
attachment addresses the combined body at index zero. Input never substitutes
authored modal coefficients, a normal-force footprint or a separately solved
lateral basis for these translations.

The two directions declare the same physical signed axis in different frames.
The file currently supplies no string-to-world frame, so the geometry author
must ensure that these declarations are consistent, including across rows.
The parser checks unit directions and available coordinates; it cannot
independently verify their relative orientation. Joint stiffness and damping
remain explicit mechanical data, rather than properties inferred from the
bulk string or body material.

The relative coordinate is
`[SX, SY] dot y_string - [GX, GY, GZ] dot u_body - REST_M`.
Spring energy, dashpot loss and equal-and-opposite reactions enter the
existing shared solve. The original bridges remain active. Distinct
directions can therefore load both transverse planes through different
body displacements; these centreline ports apply no direct twist torque.
Material loss, torsion, gas, cavity and pressure observation remain in the
same system.

Changing a body material pack rebuilds the shell basis and nodal projections.
Changing string material or surrounding gas rebuilds the string projection
with its actual mass and material operator. Cold gas, cavity and rotational
reconstruction retains the declarations, original body basis, material memory
and observers.

The native API accepts one complete nonempty binding before playback or
string motion and before nonlinear storage is attached. Scalar
`residual-flexibility exact-static-green` refuses additional directed bridges.
Missing nodes/components, absent second-plane coordinates, zero projected
body motion, invalid laws and complete-network budget violations also refuse.
Flat-plate profiles do not provide the required shell nodal translations.
Omitting the block or selecting zero rows preserves the existing path;
at most 16 additional translational bridges are accepted.

`directional-bridge.material-strings` adds a second-plane string connection
to the global X displacement of node 11 in `arched.body`, alongside the
original bridge, layered string, torsion, cavity and causal gas loading.
Using the wire, sleeve and directional body packs prepared by the arched-body
recipe in [BODY_PROFILES.md](BODY_PROFILES.md):

```bash
cp crates/fs-couple/examples/material_strings/directional-bridge.material-strings "$demo/"
cargo run --release -p fs-couple --example material_strings -- \
  "$demo/directional-bridge.material-strings" "$demo/directional-bridge.wav" \
  --body-profile crates/fs-couple/examples/material_strings/arched.body
```

This is synthetic integration data. It demonstrates an additional physical
reaction path and material-dependent pressure response; it does not establish
measured bridge properties or violin/cello/bass fidelity.

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
the complete string/section/polarization records, optional thermal supports and any rotational or translational bridge blocks,
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
every string bridge, additional directed translation, angular joint, memory dashpot, body tensor column and cavity
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

## Omitted bending modes at the supports

An optional record after `prepared-linear` (or where it would be), before
`geometric-stretching` and `cavity`:

```text
residual-flexibility exact-static-green
```

A point support seen through a string's N retained sines cannot form its kink.
The retained basis alone is therefore a Rayleigh-Ritz upper bound, and open and
stopped notes come out sharp by an error that falls as 1/N. The modes past N
still yield, almost statically, under a support force. With this record every
string's bridge spring and every fixed stop is seen through that compliance:
the exact static Green's function of the same tension, EI and pinned span,
minus the retained modes' static sum. Nothing is fitted; no tension, length or
pitch changes.

- The bridge spring condenses exactly to k/(1 + k G_bb). Its dashpot condenses
  to first order, c/(1 + k G_bb)^2.
- A fixed stop's closure gains the bridge cross term. The stops share the
  condensed compliance G_ss' - k_eff G_sb G_s'b as one matrix. The stop owner
  stores f G f / 2 and books its step-held response's numerical loss
  dR G dR / 2, reported per step as `obstacle_residual_energy_j` and
  `obstacle_residual_loss_j`. A pinching pair at one station is solved
  together on its net force.
- The receipt prints each condensed bridge, the stops' compliance and pairs.
  Every run prints `linear_frequency_upper_bound_rad_s`, the coupled stiffness
  bound the Nyquist guard checks.
- Refused: a released static preload with stops (use `forces retain`), receiver
  stops, moving or migrating stops, bound nonlinear storage, material bridge
  interfaces, additional directed translational bridges, rotational and constant-loss bridges, and a string with other than
  one bridge connection.
- The bow port carries no residual yet. That is the next step.

## Nonlinear stretching with the same bow, body and gas memory

An optional record after `prepared-linear` and `residual-flexibility` (or the
cylinder-loading record when both are omitted), before `cavity`, selects
source-derived geometric stretching:

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

The default axial extension is explicitly elastic and lossless. Original bending,
torsional, bridge, body and gas losses remain active; bending viscosity is not
reused as an unsupported axial loss law. Unretained axial channels equilibrate
statically, so a two-mode selection is not a fully resolved longitudinal wave
field. Fixed endpoints, moderate slopes, small constitutive strain and positive
local tension remain required. There is no winding slip, temperature evolution
or independent longitudinal end-force port into the
body. Axial strain reaches the receiver through its reciprocal effect on
transverse motion and the original bridge connections. This is not a claim of
measured violin/cello/bass equivalence.

### Complete material strain-rate loss

Immediately after `longitudinal-waves`, the existing axial viscosity property
can supply the complete retained moderate-slope constitutive loss:

```text
longitudinal-loss finite-strain-kelvin-voigt
```

Every original solid layer must supply `kelvin_voigt_young_viscosity` [Pa s].
The actual annular areas give `C_a = sum(eta_i A_i)` [N s]. The selected law acts
on `u_xt + y_x*y_xt + z_x*z_xt`, so transverse motion can dissipate energy even
when the axial coordinates have zero instantaneous velocity. Both transverse
polarizations, the mean strain and all retained axial harmonics participate.
Material substitution changes that same loss without a new damping ratio.

The original nonlinear/friction solve evaluates the strain derivative at
midpoint coordinates and uses the accepted displacement increment / dt.
For this quadratic strain, that is exactly the finite strain increment / dt.
The work-conjugate reaction dissipates
`C_a L dt (edot_0^2 + sum(edot_j^2)/2)` and acts reciprocally on axial and
transverse motion. Axial diagonal damping is zero for this selection, preventing
double counting. The loss is included in the existing network energy balance
and `nonlinear_dissipation_j`; endpoint acceleration pressure uses the actual
endpoint strain and velocity. No extra time step or output filter is added.

`longitudinal-loss kelvin-voigt` retains the explicitly linear `u_xt` option;
`elastic` (or omission) retains no axial viscosity. Missing, negative,
frequency-dependent or out-of-band data refuses. The complete small-signal
coupled frequency bound is checked again when later material histories are
attached; it does not certify every nonlinear harmonic. The separate moving-end
support adapter explicitly refuses this option. This does not add temperature
evolution, self-heating, finite-slope hereditary stress, winding slip or a fully
geometrically exact string law, and does not establish measured audio fidelity.

### Source-derived longitudinal material memory

After `longitudinal-waves` and the optional `longitudinal-loss` record, an
explicit uniaxial relaxation spectrum can be selected from the original string
material packs:

```text
longitudinal-memory MAX_STATES EQUILIBRIUM_MODULUS_PROPERTY ARMS DELTA_E_PROPERTY TAU_PROPERTY ...
```

The native input applies these property roles to every original homogeneous or
bonded layer. Each layer retains its own values and relaxation times; the Rust
API `with_material_longitudinal_relaxation` also supports independent per-string
and per-layer selections. Every equilibrium claim must equal that layer's
original `young_modulus`. Modulus increments are nonnegative pascals and times
are positive seconds, with explicit source frequency support. Missing data,
inconsistent quantities, source-band escape, an insufficient complete history
budget or a prepared backend without this input's exact state owner refuses.

For axial sine index `n`, each layer and arm contributes
`Delta_E * A_layer * (n*pi/L)^2 / mu_solid` to the instantaneous stiffness.
The same coupled Maxwell owner integrates its independent history and physical
loss. The projection uses the actual layer area and original solid inertia;
transverse gas-added mass and bending area moments do not enter it. A material
change therefore changes axial force and its reciprocal interaction with the
transverse strings, bearings, body, cavity and observed pressure. No axial loss
factor or sound filter is inferred. Existing body, bending, torsional and fluid
histories survive, and `MAX_STATES` is counted in the shared setup budget.

This is linear hereditary stress in `u_x` about the maintained preload. The
existing finite-slope **elastic** strain coupling remains; this addition does
not supply finite-slope hereditary stress, preload creep, self-heating, evolving
temperature or a resolved winding. Unretained axial modes still equilibrate
elastically. A separately selected axial Kelvin–Voigt viscosity is a distinct
parallel physical mechanism and must not count the same relaxation twice.
Numerical and measured-instrument convergence remain separate obligations.

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
`forces retain`. A released force-loaded pluck can instead use `forces preload`
as described below. Select gas loading and rotational bridges before
attaching stretching. Axial wave
inertia, winding slip, evolving temperature and measured instrument qualification
remain outside this composition.

### Force-loaded nonlinear plucks

`nonlinear-preload.material-strings` uses the same wire/body/sleeve packs as the
other bonded examples. It holds a declared transverse load against the complete
string, shared body and cavity before sample zero, then releases that load:

```text
forces preload 1 1
port 0 0.2 centreline 0.5
force 0 0 0
```

Both `geometric-stretching` and `longitudinal-waves` use their original
positive strain energy and exact static gradient during this initialization.
Selected conservative bow-fiber strain participates in that same solve.
The latter also equilibrates its retained axial coordinates; selected linear
axial stiffness is replaced exactly once. Both transverse planes, material
torsion, body/cavity motion and selected gas memory remain in the same network.
Changing the material EA changes the loaded shape, tension and released sound.
No pitch target or mode-by-mode pluck amplitude is introduced.

One factorization of the existing bilateral static response preconditions a
bounded energy-descent solve with safeguarded secant steps. The existing
constitutive iteration cap counts trial evaluations including backtracking;
its force tolerances and source strain bounds apply to the actual equilibrium.
The accepted report is available from `nonlinear_preload_diagnostics()`. No
sample advances during initialization. Its stored energy is an initial condition,
not an inferred history of actuator work. Selected Maxwell arms start relaxed at
that shape, and unloaded gas-memory position gauges retain zero displacement.

Compliant finger/fret and other normal contacts now participate in that same
equilibrium. Their original unilateral potentials determine contact activity
and reciprocal reactions; stationary contact has no dashpot force. Changing a
material-bound receiver changes its source-derived contact law, initial loaded
shape and released sound without changing the retained string/body solver.
Every normal force, indentation and source envelope is checked at the accepted
shape before the original mechanical state and relaxed histories publish.

`nonlinear-stopped-preload.material-strings` combines the bonded two-plane,
torsional string, causal gas memory and shared body/cavity with the existing
compliant stop. Use the same wire/body/sleeve packs and copy
`stopped-pluck.stops` into the performance directory:

```bash
cargo run --release -p fs-couple --example material_strings -- \
  "$demo/nonlinear-stopped-preload.material-strings" "$demo/nonlinear-stopped-preload.wav" \
  --stop-profile "$demo/stopped-pluck.stops"
```

The declared 0.5 N load is released at sample zero. Its support is a manufactured
law, not measured finger tissue. Existing `material-stop` records can select
the original bulk/interface/model sources instead. The 128-trial cap includes
normal-contact evaluations and all backtracking; no separate settling clock or
initial transient removal is used.

Retained shear layers, obstacle residual compliance and nonlinear-only support
of otherwise free coordinates still require a more general static solve and
refuse this initialization. Their existing dynamic `forces retain` path is unchanged.
Failed convergence, cancellation or any source/energy/force limit publishes no
accepted mechanical or material state. Refine modes, time and the physical
model before claiming instrument accuracy; the example is manufactured
integration data, not a measured violin/cello/bass specimen.

## Source-bound axial bearings on the full string span

The retained longitudinal waves can now drive the actual material shell body
through an interior axial bearing. Select `longitudinal-waves`, then place the
optional bearing records after `longitudinal-loss` (or directly after
`longitudinal-waves` when no axial-loss record is present), before bending
memory or cavity records:

```text
axial-bearings COUNT
axial-bearing STRING STATION SHELL_PART NODE TX TY TZ RX_M RY_M RZ_M REST_M interface PACK.fsintpk AREA_M2 STIFFNESS_PROPERTY VISCOSITY_PROPERTY
```

`STATION` is strictly inside the complete supported string span, including any
afterlength. `TX TY TZ` is its positive unit axial direction in the body's global
coordinates. `RX_M RY_M RZ_M` is the arm from the actual shell node to the
bearing point. The body displacement is `t dot (U + Theta cross r)`; the same
geometry applies the equal-and-opposite force and its `r cross t` moment. Shell
assembly parts address their original physical nodes while sharing one common
modal body. This requires a sourced shell or assembled-shell `--body-profile`.

The interface pack supplies stiffness and viscosity per area through the same
resolver as ordinary material bridges. Its ordered surfaces must match the
actual outer string material and every incident material state at the selected
body node. Original source identities, frozen temperature, applicability bands,
area, property choices and geometry remain attached to the rendered instrument.
A declared constant-loss anchor uses the existing prepared propagation owner.
Changing the section material changes the axial inertia, stiffness and
displacement projection; changing the interface changes the reciprocal bearing
reaction. Both affect the original body pressure through the same coupled solve.

Every bearing is counted in the complete connection/work budget. The new and
existing bridge laws, plus axial elasticity and explicit axial loss, must cover
the conservative retained linear frequency interval after coupling. This can
refuse a source whose band is too narrow; it does not certify nonlinear
harmonics or omitted modes. These are interior bearings on a fixed-end model,
not moving-end inertia, inferred break-angle prestress or a validated instrument
specimen. Supply actual geometry and interface measurements for the intended
instrument; no violin, cello or bass coefficients are invented.

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

A stop that supplies a bow's normal reaction follows that bow's
`bow_station_motion` trajectory. Keep its entries at their original source
stations in the independent `station-motion` block; the bow program supplies
the actual co-located normal and tangential endpoints on every sample. Bind
independent stop programs before bow motion in the Rust API, as the native
loader does. Other independent stops can migrate alongside it. Excessive speed
or displacement, endpoints
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
`prepared-linear`, whose air bank carries the load. Its poles count against the prepared owner's byte and work
budget (the `prepared-linear` record), not the exact owner's state and setup ceiling, so `MAX_POLES` may reach
the shared pole ceiling (16384). Without a matching finite-distance observer below, the receiver remains the
existing compact, zero-retardation pressure observation of the **loaded**
acceleration. The Rayleigh example retains one body mode to keep the complete
passive acoustic image small; it is an integration fixture, not a claim of
converged instrument timbre.

## Finite-distance pressure from the loaded material plate

`rayleigh-observed.material-strings` sends the **same** original retained plate
velocities through the shared piano/wind Rayleigh receiver. It preserves each
signed modal surface field, its physical quadrature area, distance, propagation
time and material-gas law. A cancelling mode is not replaced by an absolute
area or a positive monopole. The physical air reaction continues in the original
string/body solve; observing it introduces no second mechanical load.

The two records follow `body-radiation rayleigh` and precede `cavity none`:

```text
plate_observer X_M Y_M Z_M MAX_HZ MAX_SURFACE_CELLS MAX_REPLAY_TERMS
observer_gas stokes-kirchhoff MAX_PATHS MAX_FILTER_TAPS MAX_KERNEL_TERMS
```

The gas record uses the shared native propagation grammar. Classical
Stokes–Kirchhoff transport uses the performance's complete gas state, including
viscosity, thermal conductivity and heat capacity as well as density and sound
speed. `lossless` is an explicit comparison. Neither choice supplies missing
species relaxation data; other shared selections retain their original source
and composition eligibility checks. To declare the gas's molecular law, use
`observer_gas molecular MAX_PATHS MAX_FILTER_TAPS MAX_KERNEL_TERMS LAW`,
where `LAW` is the shared material relaxation declaration used by wind ducts.
The species and evaluated relaxation state bind to this performance's original
gas record, so an unsupported mixture or incomplete collision-partner table
refuses. Each **actual surface-to-receiver path**
is prepared separately (only bit-identical distances share a filter), before
the signed pressures sum. No average-distance atmospheric filter is applied
to the final signal.

The original three-point triangle quadrature serves both reaction and
observation. `MAX_SURFACE_CELLS` bounds all its cells, without dropping any;
`MAX_REPLAY_TERMS` bounds the actual retained convolution and derivative work.
`MAX_HZ` covers every retained body reference mode and sets the reported
equivalent-disc resolution diagnostic. It is a declared resolved band, not an
output filter or proof that nonlinear excitation remains in-band. Mesh and
structural-mode refinement still determine spatial accuracy.

The plate occupies z=0 in its original chart; the microphone must be at least
5 cm into z>0. One observed half-space produces one pressure signal even if
the reaction load describes two independent exterior faces. Cavity mouths and
other observers are refused here because their sound needs complete spatial
observation. Closed shells use their own routes below. This is an infinite
baffle problem, not an unbaffled violin exterior or a room model.

All numerical derivative, interpolation and common gas padding delays are
reported and retained, alongside the physical flight time. The native pressure
writer receives these accepted observations on the same mechanical timeline;
it does not normalize the signal, remove latency or append an unrequested
ring-out. With the packs and regional profile prepared as above:

```sh
cargo run --release -p fs-couple --example material_strings -- \
  "$demo/rayleigh-observed.material-strings" "$demo/rayleigh-observed.wav" \
  --body-profile "$demo/orthotropic.body"
```

Changing the body material, its grain or the gas therefore rebuilds the
original mechanical response and its finite-distance microphone from data.
This fixture verifies the connection; its synthetic sources and one retained
body mode do not qualify a measured instrument.

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
connections on the full span **including afterlength**, with optional independent
translational directions from shell nodes and scalar rotational spring/dashpot joints. A full 3D unilateral bridge-contact
law, coupled f-hole/interior/exterior geometry, distributed cavity modes,
resolved winding contact/slip and solved finger disengagement remain outside
this composition. The computed planar and closed-shell exterior reactions and
the optional finite-distance closed-shell observer are described above.
Without that observer, body and neck pressure retains the existing compact
baffled, zero-retardation approximation. Source-state validity is not measured
instrument-fidelity evidence.

Input caps: 1 MiB performance; 8 MiB per distinct material pack and 64 MiB total
raw packs; 8 strings, at most 16 layers per string, 16 bow contacts, 16 rotational
and 16 additional translational bridges, 64 force ports, 128 complete-network modes including fluid memory,
16 plate modes; 16,384 events per input family; 262,144 compiled controls per
family; 16,777,216 source/projection work units per compiler; at most 384 kHz,
600 seconds or 28,800,000 samples. Existing numerical owners retain their stricter
limits. The example additionally caps total mechanical storage at 100 J,
connection forces at 1000 N, bow traction at 100 N, slip at 10 m/s and observer
pressure at 1 MPa. These are refusal ceilings, not assertions of model accuracy.

## Native checks

The complete three-dimensional tensor body path is runnable with
`tensor-body.material-strings` and selected `tensor-body.body`. It shares the
same source property binding and hereditary normal-traction elimination as
material reeds, and retains body feedback onto the strings through the
original bridges. See [complete tensor memory](BODY_PROFILES.md#complete-three-dimensional-tensor-memory)
for the source-pack recipe, native grammar and explicit model/capacity limits.
Selecting `tensor-body-shared.body` with the same performance and packs
refines the original cantilever to eight triangles while retaining all 144
original material histories through exact equal-time coordinate sharing.
It keeps the original 36-state ceiling and reconstructs every point's source
strain before admission. See [shared history coordinates](BODY_PROFILES.md#exact-shared-material-history-coordinates)
for the optional native selector and its initial-history domain.

```bash
cargo test -p fs-couple --example material_strings -- --test-threads=1
```

The tests exercise actual normalized-pack decoding, material-only substitutions,
shared gas response, source-derived cavity-neck inertia/loss and end conditions,
causal cylinder memory, joint pluck/double-bow control, callback equivalence, exact
existing PCM output, source/input refusals, and finite-distance shell pressure
from the original accepted body velocities. The fixture TSVs use the existing
`xtask matdb-pack` compiler, not a new importer.
