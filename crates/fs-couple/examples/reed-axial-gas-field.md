# Gas composition and temperature along a reed bore

The same performance loader accepts prescribed `gas_field` stations or generates
them from a conservative `gas_steady` heat/species problem. Both feed the original
reed, duct, exterior load and microphone. The steady example and grammar are
described below; they use the same source packs as the prescribed example.

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
field. The file admits cylinder/cone geometry and the shared smooth radius
profiles described below. Nonuniform midpoint fields still refuse profiles;
select continuous propagation for a varying medium in a smooth bore. Tone-hole
field combinations retain the original continuous duct owner's restrictions.

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

## Generate a steady field from flow, wall heat and species injection

[reed-steady-gas.performance](reed-steady-gas.performance) supplies nitrogen at
330 K, wall cooling, and a downstream carbon-dioxide injection. The existing
`SteadyAxialGas` owner solves stationary mass/species transport and source-bound
NASA-9 enthalpy balance, then publishes the same `AxialGasField` consumed above.
Changing the source species, inlet fractions, flow or wall conditions changes
the generated thermodynamic and transport states and the actual reed feedback.

Replace the prescribed block with:

```text
gas_steady continuous STEPS_PER_SPAN CELLS INJECTIONS mole SPECIES_COUNT ID...
gas_inlet MASS_FLOW_KG_S T_K FRACTION...
gas_cell END_M WALL_CONDUCTANCE_W_K WALL_T_K INITIAL_T_K
gas_injection CELL MASS_FLOW_KG_S T_K FRACTION...
gas_solve MIN_T_K MAX_T_K H_TOL_J_KG MAX_ITER DT_S MAX_RELAX MAX_SUBSTEPS P_TOL_W E_TOL_J MAX_MACH
```

`midpoint CELLS_PER_SEGMENT` is also available in the first record. Supply one
`gas_cell` row per thermal cell and then one `gas_injection` row per injection;
zero injections requires no injection rows. Injection indices start at zero.
Fractions on every inlet follow the declared species-column order. They are
canonicalized with the source IDs once, retaining exact zero components and the
original species, NASA-9 and transport pack identities.

Cell ends must increase from zero to the exact complete bore length. The loader
integrates each cell's volume and lateral wall area from the original cylinder,
cone and smooth-profile records, including cells that cross geometry boundaries. The format
does not duplicate bore radii or volumes. `WALL_CONDUCTANCE_W_K` is the explicitly
supplied total wall-to-gas conductance of that cell; it is not a guessed heat
transfer coefficient. `WALL_T_K` is a prescribed wall temperature and
`INITIAL_T_K` only seeds the bounded steady relaxation. The wall's own thermal
storage is not advanced by this stationary problem.

The inlet mass flow is a prescribed stationary throughflow, independent of the
oscillatory reed flow. It sets enthalpy/species transport and the admitted mean
Mach number; this does not solve a steady pressure drop from the blowing
gesture. Acoustic propagation uses the generated frozen field at the ambient
static pressure. This route requires an unbranched bore: its
stationary mass balance does not infer leakage through tone holes or their
chimneys. A separately solved branch field can still use the prescribed route
where admitted by the original duct owner.

The caloric bracket must lie inside the loaded source domains. The enthalpy
residual tolerance is absolute J/kg, cell power tolerance is W, and the accepted
step energy tolerance is J. `MAX_ITER` is at most 128, thermal cells and
injections are each at most 64, and relaxation/substep limits are each at most
1,000,000. The parser also bounds retained state/species products. `MAX_MACH`
must be positive and no greater than 0.1. Nonconvergence, unresolved geometry,
unsupported source or phase, excess mean flow and conflicting segment/field
declarations refuse before a renderer is returned.

Build the packs using the commands above, then render the steady example by
changing the performance and output paths in that same `music_render wind`
command to `crates/fs-couple/examples/reed-steady-gas.performance` and
`target/reed-steady-gas.wav`. The same loader also serves reed ensemble parts.
No extra physics flags or compiled gas names are required.

`ReedPerformance::steady_gas()` retains the original problem, declared solver
controls and accepted `SteadyGasSolution`. The wind/ensemble sidecar includes
the generated field plus cell geometry, inlet/injection compositions, mass
flows, actual power and energy residuals, work used and maximum mean Mach.
The same solution field is used by the acoustic operator and receiver; the
sidecar does not reconstruct it from the rendered audio.

This is a stationary donor-cell heat/species model with primitive interpolation
between its emitted faces. The existing thermal tests establish first-order
spatial convergence. Diffusion, reaction, condensation, transient atmosphere,
mean-flow acoustic correction and measured instrument qualification require
their own existing-owner extensions and validation.

## Smooth bore geometry in the same performance

Each of these records occupies one entry in the existing `segments` list:

```text
profile exponential STEPS R0_M AREA_FLARE_PER_M LENGTH_M
profile power STEPS R0_M X0_M EXPONENT LENGTH_M
profile spline STEPS natural COUNT X_M R_M ...
profile spline STEPS clamped LEFT_SLOPE RIGHT_SLOPE COUNT X_M R_M ...
```

They use the same `fs-duct::RadiusProfile` constructors as the acoustic API:
`r(x) = R0 exp(AREA_FLARE*x/2)`, `r(x) = R0 ((X0+x)/X0)^EXPONENT`, or a C2
cubic through the supplied stations. Spline abscissae retain their original
values; the first is the local origin and the last minus the first is the
segment length. Slopes are `dr/dx`. Nonpositive/unresolved profiles, invalid
station order, failed positivity admission and invalid step budgets return the
original constructor refusal; no endpoint cone is substituted.

For example, replace `cylinder 0.0015 0.57` in the steady example with:

```text
profile exponential 16 0.0015 0.4 0.57
```

The existing `music_render wind` and ensemble commands then generate the
heat/species field and microphone pressure using that same smooth bore.
`STEPS` is the declared profile integration resolution. Thermal geometry uses
fourth-order two-node Gauss panels, split at measured spline knots and thermal
cell boundaries, to integrate volume and the sloped lateral wall area. The
complete thermal binding admits at most 65536 panels. Cylinder/cone measures
retain their exact formulas. Refine profile `STEPS`, thermal cells and
`gas_steady continuous STEPS_PER_SPAN` separately: they control geometry,
stationary transport and acoustic field integration, respectively. No geometric
error enclosure, condensation law or additional material data is implied.
