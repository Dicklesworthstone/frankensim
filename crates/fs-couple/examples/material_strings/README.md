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
it affects sound through the physical surface/contact coupling when that coupling
is active. None of these fixtures is a measurement of gold, silver or a real string.

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
fixed-extension supports, fixed rotational ends, one simply supported **isotropic flat plate**, and one
optional uniform-pressure cavity/circular neck. Bridges are interior transverse
connections on the full span **including afterlength**. It does not model carved
wood shells, torque-bearing bridge contact, f-hole geometry, distributed cavity
modes, helical winding/slip, nonlinear tension, solved finger disengagement,
bow-hair normal-load distribution, or exterior radiation back-loading. The body
and neck observations retain the existing compact baffled, zero-retardation
approximation. Source-state validity is not measured instrument-fidelity evidence.

Input caps: 1 MiB performance; 8 MiB per distinct material pack and 64 MiB total
raw packs; 8 strings, at most 16 layers per string, 16 bow contacts, 64 force ports, 128 complete-network modes,
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
shared gas response, joint pluck/double-bow control, callback equivalence, exact
existing PCM output, and source/input refusals. The fixture TSVs use the existing
`xtask matdb-pack` compiler, not a new importer.
