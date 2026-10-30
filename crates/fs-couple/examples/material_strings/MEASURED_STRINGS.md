# Measured whole-string coefficients on a laboratory plate

`measured-cello-g.material-strings` uses the published nominal section
coefficients for a Thomastik-Infeld Dominant cello G string from the existing
[Lampis source pack](../../../../data/matdb/seed-v1/cello-g-strings-whole-lampis-2024/manifest.tsv).
The receiver is the **manufactured laboratory plate** from `body-manifest.tsv`,
with the fixture's explicit interior spring/dashpot connection. This is an
unqualified assembly-coefficient experiment, not a complete cello, a reproduced
monochord measurement, or a validated cello recording. The original measured
test span is 0.70 m; the laboratory support and receiver are additional model
choices, so no nominal pitch is imposed.

The source pack cites Lampis, Mayer and Chatziioannou,
[JASA Express Letters 4, 113201 (2024)](https://doi.org/10.1121/10.0034330),
Table 1, and [Acta Acustica 8, 44 (2024)](https://doi.org/10.1051/aacus/2024034),
Table 1. The fixture selects the JASA Express Letters Dominant values; it does
not mix the two sessions. Attribution and the source observations remain in
the normalized pack.

## Compile the original data and render

Run from the repository root with the existing Rust toolchain:

```bash
case_dir="$(mktemp -d)"
cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/cello-g-strings-whole-lampis-2024/manifest.tsv \
  --out "$case_dir/cello-g-whole.fsmatpk"
cargo run -p xtask -- matdb-pack \
  --manifest crates/fs-couple/examples/material_strings/body-manifest.tsv \
  --out "$case_dir/body.fsmatpk"
cp crates/fs-couple/examples/material_strings/measured-cello-g.material-strings \
  "$case_dir/performance.material-strings"
cargo run --release -p fs-couple --example material_strings -- \
  "$case_dir/performance.material-strings" "$case_dir/measured-string-lab.wav"
```

These commands use the normal material compiler and the existing native
performer. The fixture starts from its declared transverse preload and releases
it into the same reciprocal string/body solve. It retains two transverse
polarizations and emits physical pressure through the original observation and
PCM path. The bow track is released; this example does not validate a bow-contact
law. No pitch correction, recorded waveform, or audio filter supplies the
string response.

The fixture selects `prepared-linear 64 67108864 200000000000`. The existing
fixed-system owner precomputes the coupled physical propagator and its work
integrals once, making playback practical without repeatedly integrating the
same operator. The original material coefficients, mechanical clock, state
variables, reciprocal couplings and work model remain in force. The three
explicit allowances bound modes, preparation bytes and preparation work.

To open the same sources in Instrument Lab:

```bash
python3 apps/instrument-lab/tools/render_material_strings.py \
  --prepare-only --performance "$case_dir/performance.material-strings" \
  --out "$case_dir/lab-experiment"
```

Open `$case_dir/lab-experiment/experiment.native.json` in the app's native
performer. Run the app with the native `material_strings` executable configured,
as described in [the native comparison guide](../../../../apps/instrument-lab/NATIVE_COMPARISON.md).
The packet retains the original whole-string and plate packs. A material or
support edit is admitted by the same Rust owner when rendering.

### Prepare directly from original source manifests

The launcher can invoke the existing Rust material compiler for each declared
source. Build the ordinary `xtask` executable once (`cargo build -p xtask`), then
prepare this experiment without manually producing intermediate packs:

```bash
python3 apps/instrument-lab/tools/render_material_strings.py \
  --prepare-only \
  --performance crates/fs-couple/examples/material_strings/measured-cello-g.material-strings \
  --pack-compiler target/debug/xtask \
  --pack-manifest cello-g-whole.fsmatpk=data/matdb/seed-v1/cello-g-strings-whole-lampis-2024/manifest.tsv \
  --pack-manifest body.fsmatpk=crates/fs-couple/examples/material_strings/body-manifest.tsv \
  --out "$case_dir/from-original-sources"
```

Each binding invokes `xtask matdb-pack --manifest ... --out ...`. The canonical
compiler resolves the manifest's source files, units, uncertainty, validity,
identity and attribution; Python does not interpret those properties. Replacing
a manifest or editing its declared source data compiles a fresh pack at the
same referenced path. The resulting pack bytes enter both the portable
experiment and the native render. Use `--renderer PATH --lab-experiment` in
place of `--prepare-only` to compile, render and save the same experiment.

For a source that must retain a canonical material identity, add
`--pack-material-state REFERENCED_PATH CHEMISTRY PHASE PROCESS REVISION`
for that manifest target. The launcher forwards these explicit identity fields
to the same compiler's `--material-state` option. Its verified material-card
bytes then reach the native performer and the saved Lab experiment unchanged.
This matters when an ordered bow/string or finger/string contact names the
original material identities. It does not infer a contact pair or change any
material property. Identity fields must describe the actual source state;
revision is an unsigned 32-bit integer. Without this option the compiler keeps
its original property, interface, model or species envelope. The compiler
refuses applying a material-card identity to a non-material source.

Precompiled `--pack` bindings and manifest bindings can be combined, but each
referenced path has exactly one explicit binding. Compilation failure stops
the experiment; a previously compiled pack is never substituted for a failed
source. This path changes no string setup, selectors, body, gas or gesture.
Different measured spans, diameters or tensions still require the corresponding
explicit setup and are checked by the native measured-string owner.

## Explicit property and state selection

The section record is one line:

```text
section measured AXIS_COUNT AXIS VALUE ... LENGTH_PROPERTY DIAMETER_PROPERTY TENSION_PROPERTY LINEAR_DENSITY_PROPERTY BENDING_STIFFNESS_PROPERTY no-intrinsic-loss
```

For this source the coordinate is `temperature_known 0`. It means the
measurement temperature is unstated. It is not zero kelvin, room temperature,
or the temperature from the performance's `state` or `gas` row. Body material
state and laboratory atmosphere remain independent inputs; neither fills in
missing measurement conditions.

| Selected property | Physical role | Nominal source value |
| --- | --- | --- |
| `playing_length` | Full source test span | 0.70 m |
| `diameter_dominant` | Finished outside diameter | 0.00115 m |
| `tension_dominant` | Applied source test force | 116.3 N |
| `linear_density_dominant` | Published effective line-mass estimate | 0.0062 kg/m |
| `bending_stiffness_dominant` | Published whole-string flexural-rigidity estimate | 0.000303 N m² |

Each selector must name a distinct positive scalar of the required dimensions
with constant-within-validity source support. The original card, query,
observations, uncertainty and usage receipts remain attached. A declared
frequency applicability band constrains the retained modes; an unstated band
remains unstated. The string row's length, radius and fixed tension must agree
with the selected source values. Changing those setup values alone refuses.
A different complete measured source can supply its own coefficients and setup
without modifying the physics implementation.

The existing prestressed Euler–Bernoulli owner uses line mass and flexural
rigidity directly. It does not turn them into fictitious bulk density,
Young's modulus, a core/wrap construction, shear modulus or axial stiffness.
Torsion, fixed-extension/thermal prestress, geometric stretching, longitudinal
waves, layer relaxation and source-bound material contact need additional
physical source information and are refused for this section.

## What the experiment does not establish

The source's own observation explains that line mass was **inferred from**
`T / (2 L f0)^2`, and flexural rigidity was fitted from modal-frequency
spacing. The line mass is not an independently weighed dry-solid mass.
The reference measurement atmosphere is not supplied by this pack, and its
fluid contribution has not been de-embedded. No reference gas, correction or
inertial subtraction is inferred here.

Consequently the published-data fixture does **not** select
`cylinder-gas-loading`, which would add a new fluid-inertia continuum to those
effective coefficients. Its declared laboratory gas still supplies the
existing sampled cylinder resistance; removing the optional continuum does
not remove all gas effects. A physically qualified prediction of this measured
string under a different atmosphere needs the reference conditions and an
appropriate coefficient identification. The stated source uncertainty alone
does not establish that qualification.

The generic measured-section API can use the shared causal cylinder-fluid
operator when the supplied coefficients are independently established
mechanical/dry coefficients. The synthetic bank integration tests exercise
that path, including actual radius, fluid inertia, reciprocal force and closed
work, without presenting it as a validated Lampis atmosphere substitution.

`no-intrinsic-loss` is an explicit modeling limitation. The source includes a
mean quality factor, but a mean over modes does not identify an intrinsic
bending, tension, torsional or frequency-dependent loss law. The performer
does not convert that mean into viscosity or silently apply a tension-loss
table. The laboratory model still dissipates energy through its declared gas
and support/body laws. Its decay is not a validated reproduction of the
measured string's decay.

A full cello additionally needs sourced body and bridge geometry/materials,
cavity and radiation data, support and contact laws, appropriate resolution,
and comparison with measurements. The laboratory plate does not stand in as a
claimed cello body.
