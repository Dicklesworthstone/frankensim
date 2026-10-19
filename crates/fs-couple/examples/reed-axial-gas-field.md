# Authored gas composition and temperature along a reed bore

[reed-axial-gas-field.performance](reed-axial-gas-field.performance) drives the
existing massive reed, regional duct and baffled pressure receiver from one
authored axial gas field. The inlet contains pure nitrogen. Carbon dioxide is
exactly absent there and reaches a mole fraction of 0.5 at the outlet. Temperature
also varies along the bore. Geometry, source material packs and primitive states
are data; the same thermochemistry, transport and acoustic owners handle each
admitted mixture.

This is a prescribed stationary field with a common static pressure. It does
not solve gas transport, diffusion, heating or chemical reactions. The fixture
exercises material propagation; it is not a measured clarinet or oboe specimen.

## File grammar

Place one optional block after the complete `segments` list and before
`compile_limits`:

```text
gas_field continuous STEPS_PER_SPAN STATIONS mole SPECIES_COUNT ID...
gas_field midpoint CELLS_PER_SEGMENT STATIONS mole SPECIES_COUNT ID...
gas_station X_M T_K FRACTION...
```

The example selects the existing continuous Magnus propagation owner; the older
numeric spelling `gas_field 8 ...` remains shorthand for `gas_field midpoint 8 ...`.
The mode changes only the field's spatial integration through the shared duct
owner, with the same source-backed primitive states and pressure renderer.

The example supplies:

```text
gas_field continuous 16 3 mole 2 N2 CO2
gas_station 0 300 1 0
gas_station 0.285 320 0.75 0.25
gas_station 0.57 310 0.5 0.5
```

Coordinates are metres, temperatures are kelvin, and fractions follow the
declared species-column order. IDs must match the canonical IDs in explicitly
loaded species packs. Changing their order together with their columns preserves
the field. Each station must supply all fractions; authored fractions are never
repaired or normalized. Exact zeros stay zero. A declared species can be zero at
every station and later be enabled solely by changing its fractions.

The ambient record supplies the one common static pressure, 100000 Pa here.
The first station must be exactly zero, subsequent coordinates strictly increase,
and the last station must exactly equal the sum of the authored bore lengths.
An explicit gas suffix on any segment conflicts with the field. The injected
`ReedPerformance::from_bytes_with_gas_field` API also refuses a second field
in the file. An optional `exterior` record retains its existing independent
meaning and must satisfy the existing common-pressure rule.

The field reevaluates the existing NASA-9 mixture, ideal-gas EOS and pure-source
Wilke/WMS transport at authored stations and the selected propagation owner's
material sample positions. It interpolates primitive temperature and mole
fractions using `AxialGasField`, without blending precomputed density or sound
speed. The exact inlet state sets reed density and characteristic impedance.
Local states drive bore propagation; the exact outlet state supplies the default
exterior loading and microphone medium. An explicit exterior replaces that
default in both consumers. Downstream composition can therefore alter the reed
motion through acoustic feedback as well as the rendered pressure.

Continuous mode evaluates the original field at the two quadrature positions
of every fourth-order Magnus step, splitting exactly at material/geometry
boundaries. The step count applies to each smooth span, and the existing
`fs-duct` owner limits the complete binding to 65536 steps. Midpoint mode
samples the field at cell midpoints and permits at most 64 generated cells.
Both modes preserve the original geometry and operation order for a uniform
field. The file admits cylinder/cone geometry; nonuniform tone holes are
refused. The underlying continuous duct API also supports smooth profiles.

Cold field input is bounded to 2..=4096 stations, 1..=128 species and 16384
station/species pairs. Midpoint resolution must be 1..=64 cells per segment;
continuous resolution must be 1..=65536 steps per smooth span and also satisfy
the owner's aggregate work cap. Every active source must admit the sampled
temperature and pressure under its own validity and phase policy.

## Build the supplied source packs and render

Run from the repository root. These commands compile the repository's existing
human-readable source data; generated binary packs are not checked in.

```sh
mkdir -p target/reed-field-packs
for source in nitrogen nitrogen-nasa9-mcbride2002 nitrogen-transport-lemmon2004 carbon-dioxide carbon-dioxide-nasa9-mcbride2002 carbon-dioxide-transport-laesecke2017-huber2016
do
  cargo run -p xtask -- matdb-pack \
    --manifest "data/matdb/seed-v1/$source/manifest.tsv" \
    --out "target/reed-field-packs/$source.pack"
done

cargo run -p fs-couple --bin music_render -- wind \
  crates/fs-couple/examples/reed-axial-gas-field.performance \
  target/reed-axial-gas-field.wav --block 113 \
  --gas-component target/reed-field-packs/nitrogen.pack \
    target/reed-field-packs/nitrogen-nasa9-mcbride2002.pack \
    target/reed-field-packs/nitrogen-transport-lemmon2004.pack \
  --gas-component target/reed-field-packs/carbon-dioxide.pack \
    target/reed-field-packs/carbon-dioxide-nasa9-mcbride2002.pack \
    target/reed-field-packs/carbon-dioxide-transport-laesecke2017-huber2016.pack \
  --gas-water none --gas-phase frozen
```

Choose a new output filename when rerunning: the existing WAV writer refuses to
overwrite outputs. The same file enters an ensemble through
`music_render ensemble OUT.wav --full-scale-pa 2 --reed INPUT.performance`
with the same `--gas-component`, `--gas-water` and `--gas-phase` arguments.
Both loaders use the same library parser. Source packs need no implicit lookup,
and arbitrary gas names do not assign physical constants.

For a pure-nitrogen comparison, retain the two declared source IDs and replace
the three station rows with `0 300 1 0`, `0.285 300 1 0` and
`0.57 300 1 0` (each prefixed by `gas_station`). To hold the external microphone
medium fixed across comparisons, add `exterior 300 100000 0` immediately before
`termination flanged` in both files. This isolates the changed bore feedback
from a changed exterior medium.

## Sources and retained evidence

The commands bind [nitrogen](../../../data/matdb/seed-v1/nitrogen/manifest.tsv)
and [carbon dioxide](../../../data/matdb/seed-v1/carbon-dioxide/manifest.tsv)
species associations, their explicit NASA-9 packs, and the indicated pure
transport packs. Their manifests retain the original source references,
normalization evidence, redistribution terms and bounded validity. The field
does not extend those domains.

`ReedPerformance::gas_field()` retains the complete declared source catalogue
and original stations. `gas_records()` also exposes each authored knot's
resolved state and active source-pack hashes, which the existing wind/ensemble
sidecar writer records. A locally inactive component has no physical evaluation
receipt for that sample. It remains in the field catalogue, allowing its authored
fraction to become nonzero elsewhere.

The pack adapter needs a positive initial composition to bind that catalogue.
The parser uses an exactly summing dyadic seed solely for this construction.
It never evaluates that seed at any physical coordinate; all physical samples
use the authored field. Focused tests compare different seed compositions through
both resolved states and actual rendered pressure in continuous and midpoint
modes, check zero-inlet and globally
inactive species, source/refusal paths, downstream feedback and callback parity,
and preserve the original pressure bytes when the optional block is absent.

The resulting model retains the existing regional transfer fit, compact reed
law and infinite-baffle receiver assumptions. It adds authored access to the
shared gas field, without claiming complete fluid motion or measured instrument
accuracy.
