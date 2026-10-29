# A live heated string in a finite thermochemical gas

`finite_gas_string` advances one vibrating material string while electrical
heating changes its temperature, resistance and material loss. A fixed-support
thermoelastic scenario also changes tension, stiffness and pitch through sourced
material properties, without resetting the vibration. The same string
exchanges heat with a finite gas parcel, whose temperature, pressure, caloric
state and transport properties evolve during the sound. The result is a pressure
WAV and a synchronized CSV of motion, temperatures, electrical response and
whole-system energy.

The initial string is in static equilibrium under an explicitly specified point
load. The load is held for the declared number of samples and then removed.
Modal displacement and velocity persist throughout the run; changing temperature
does not restart the vibration or synthesize a fresh decaying tone.

Each accepted sample composes the existing string, conductor, caloric, gas and
conduction owners:

- Source-resolved electrical work becomes solid heat once. Current and voltage
  drives retain their sign; their dissipated power is nonnegative.
- Kelvin–Voigt bending dissipation heats the solid. Aerodynamic drag heats the
  surrounding gas.
- A single signed conductive transfer leaves the solid and enters the gas.
  Its series resistance uses the source conductivities and explicit layer
  lengths and area.
- Accepted temperatures update the material loss, electrical resistance and
  gas state used by subsequent motion and pressure observations.
- Fixed-support thermoelastic strings also update sourced thermal strain and
  elastic modulus. Axial elastic storage and axial/transverse elastic entropy
  supply signed reversible heat to the same solid/gas exchange transaction.

The coupled library publishes the string and parcel together only after all
source, temperature, energy and cancellation checks pass. The executable waits
for an accepted complete trajectory before publishing its outputs.

## Run with source packs

The solid source in this directory is a **manufactured CC0 integration fixture**.
Its deliberately small heat capacity, 0.1 J/(kg K), makes feedback visible in a
short run. Its coefficients are not measurements of a real wire or instrument.
The nitrogen inputs below use the repository's existing species association,
NASA-9 caloric source and paired dilute-gas transport source. Species names do
not supply missing physical data.

From the repository root, compile all four source packs and run the native
binary:

```bash
run_dir="$(mktemp -d)"
demo=crates/fs-couple/examples/finite_gas_string
cargo run -p xtask -- matdb-pack \
  --manifest "$demo/manifest.tsv" \
  --out "$run_dir/string.pack"
cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/nitrogen/manifest.tsv \
  --out "$run_dir/nitrogen.species"
cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/nitrogen-nasa9-mcbride2002/manifest.tsv \
  --out "$run_dir/nitrogen.nasa9"
cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/nitrogen-transport-svehla1995/manifest.tsv \
  --out "$run_dir/nitrogen.transport"

cargo run --release -p fs-couple --bin finite_gas_string -- \
  "$demo/heated.string" "$run_dir/string.pack" "$run_dir/heated.wav" \
  --gas-record crates/fs-couple/examples/self_heating_vibration/nitrogen.gas \
  --gas-component "$run_dir/nitrogen.species" "$run_dir/nitrogen.nasa9" "$run_dir/nitrogen.transport" \
  --gas-water none --gas-phase frozen \
  > "$run_dir/heated.csv"
```

The supplied scenario renders 3,200 samples at 16 kHz, or 0.2 seconds, with four
string modes and a 0.1 A drive. The initial gas is 300 K and 85,000 Pa. The WAV
path must be new, so a later experiment cannot overwrite a previous render.
Source identities, initial gas provenance, output scale and clipping count go
to stderr; the CSV goes to stdout.

The binary can also be run as `cargo run -p fs-couple --example
finite_gas_string -- ...` with the same arguments. Use the workspace's native
build execution policy when building either target.

## Heating a string with fixed supports

`expanding.string` uses the same electrical, gas and conduction path, with an
explicit stress-free reference length of 0.648 m at 300 K and loaded length
0.65 m. Its `fixed_support` prescription derives tension from
`E(T) A [(L-Lref)/Lref - integral(alpha dT)]`. It supplies no target pitch or
prescribed tensile force.

Compile `expanding-manifest.tsv` in place of the first solid manifest, then use
`expanding.string` with the same gas arguments:

```bash
cargo run -p xtask -- matdb-pack \
  --manifest "$demo/expanding-manifest.tsv" \
  --out "$run_dir/expanding.pack"
cargo run --release -p fs-couple --bin finite_gas_string -- \
  "$demo/expanding.string" "$run_dir/expanding.pack" "$run_dir/expanding.wav" \
  --gas-record crates/fs-couple/examples/self_heating_vibration/nitrogen.gas \
  --gas-component "$run_dir/nitrogen.species" "$run_dir/nitrogen.nasa9" "$run_dir/nitrogen.transport" \
  --gas-water none --gas-phase frozen > "$run_dir/expanding.csv"
```

This second pack is also a **manufactured integration fixture**. It declares
`alpha = 2e-5 /K`, affine `E(T) = 2e9 - 1e6 (T-300)` Pa, and caloric heat
capacity 10 J/(kg K). The larger heat capacity admits the full thermoelastic
storage contribution. These are synthetic laws, not measurements of an
instrument string. To isolate expansion, copy the source, change the expansion
coefficient to zero, and rebuild the pack. To isolate stiffness softening,
replace the modulus curve with a constant source. Geometry, motion, electrical
heating, heat transfer and observation use the same implementation in each run.

The CSV exposes the actual evolving tension and fundamental frequency alongside
the signed reversible heat and stored-energy changes. The waveform remains the
observer of the continuing physical state; there is no fitted pitch envelope.

## Bonded layers with independent materials and heat

`bonded.string` uses the same executable for a core and one or more circular
sleeves. Each layer retains its own card, mass, caloric cell, temperature,
electrical response, thermal strain, elastic storage and bending loss. The
native layered section integrates mass, axial stiffness and bending stiffness;
the native conductor divides the total drive current between the layers.
Changing a sleeve pack changes these coupled properties without changing code.

Compile the existing expanding pack for the core and the new manufactured
sleeve pack, and place the scenario beside its relative sleeve pack path:

```bash
cargo run -p xtask -- matdb-pack \
  --manifest "$demo/expanding-manifest.tsv" --out "$run_dir/core.pack"
cargo run -p xtask -- matdb-pack \
  --manifest "$demo/sleeve-manifest.tsv" --out "$run_dir/sleeve.pack"
cp "$demo/bonded.string" "$run_dir/bonded.string"
cargo run --release -p fs-couple --bin finite_gas_string -- \
  "$run_dir/bonded.string" "$run_dir/core.pack" "$run_dir/bonded.wav" \
  --gas-record crates/fs-couple/examples/self_heating_vibration/nitrogen.gas \
  --gas-component "$run_dir/nitrogen.species" "$run_dir/nitrogen.nasa9" "$run_dir/nitrogen.transport" \
  --gas-water none --gas-phase frozen > "$run_dir/bonded.csv"
```

Both solid packs are manufactured integration fixtures. The sleeve doubles the
core fixture's density, expansion coefficient, resistivity and bending viscosity.
These choices demonstrate material substitution, not measured instrument
construction. The geometry is a perfectly bonded core and annular sleeve;
it does not represent a sliding helical winding or a wound piano string.

In bonded mode, `heat` applies its solid power to the outer cell. Neighboring
cells exchange equal and opposite Fourier heat through their original source
conductivities. The two authored resistance geometry factors are distance/area
[inverse metres]; the interface conductance is
`G = 1 / (shape_inner/k_inner + shape_outer/k_outer)`.
The fixture uses logarithmic radial resistances between mid-area representative
radii, with full cylindrical interface coverage and no extra contact resistance.
Those factors and the interface area are geometric discretization data. Add
annuli and refine them when a single temperature per material cannot resolve the
radial gradient; this executable does not automatically generate a radial mesh.

Only the outer lateral surface exchanges heat with the finite gas. Source
conductivity pins, exposed area and the native lumped Biot bound are checked;
there is no additional ambient convection reservoir. Gas drag heats that same
parcel. Thermal strain, tension, modulus, layer loss and gas properties feed
the continuing motion and physical pressure. Axial and transverse reversible
heat returns to the specific layer that supplied it.

The first 27 CSV columns retain their order. `solid_temperature_k` is the
mass-weighted layer temperature in bonded mode; appended
`layer_0_temperature_k`, `layer_1_temperature_k`, ... expose the individual
cells, from core outward. Initial and endpoint electrical columns retain their
time meanings. Bonded Joule heating uses the accepted endpoint current and
resistance, so its per-sample work is `endpoint_current_a² * endpoint_resistance_ohm * dt`.
The original uniform path keeps its initial-state electrical quadrature.
The joint accepted-state energy balance includes every layer and the finite gas.

## Pressure-driven gas exchange from data

Uniform fixed-tension, uniform fixed-support thermoelastic and bonded strings
all accept the same finite-reservoir `nozzle` input. Their existing electrical
drive and reciprocal conductive exchange remain active during species flow.
`heated-flow.string` adds a reservoir to the first uniform heated-string example,
using the original `string.pack` and nitrogen packs compiled above:

```bash
cp "$demo/heated-flow.string" "$demo/nitrogen-reservoir.gas" "$run_dir/"
cargo run --release -p fs-couple --bin finite_gas_string -- \
  "$run_dir/heated-flow.string" "$run_dir/string.pack" "$run_dir/heated-flow.wav" \
  --gas-record crates/fs-couple/examples/self_heating_vibration/nitrogen.gas \
  --gas-component "$run_dir/nitrogen.species" "$run_dir/nitrogen.nasa9" "$run_dir/nitrogen.transport" \
  --gas-water none --gas-phase frozen > "$run_dir/heated-flow.csv"
```

The same `nozzle` row can be placed after `gas` in `expanding.string` to retain
its fixed-support thermal expansion and elastic feedback while the atmosphere
changes. Select `drive none` for an unpowered string or retain either the
explicit signed current or voltage drive. A material string does not need a
second solid layer to exchange gas with a reservoir.

`bonded-flow.string` adds a finite second reservoir through the existing
source-resolved converging-nozzle owner. It uses the same core, sleeve and
nitrogen packs compiled above. Copy its explicit reservoir record beside the
scenario and run the same binary:

```bash
cp "$demo/bonded-flow.string" "$demo/nitrogen-reservoir.gas" "$run_dir/"
cargo run --release -p fs-couple --bin finite_gas_string -- \
  "$run_dir/bonded-flow.string" "$run_dir/core.pack" "$run_dir/bonded-flow.wav" \
  --gas-record crates/fs-couple/examples/self_heating_vibration/nitrogen.gas \
  --gas-component "$run_dir/nitrogen.species" "$run_dir/nitrogen.nasa9" "$run_dir/nitrogen.transport" \
  --gas-water none --gas-phase frozen > "$run_dir/bonded-flow.csv"
```

The authored reservoir starts at 300 K and 95,000 Pa in 2 × 10⁻⁶ m³. The
restriction area is 10⁻¹⁰ m² and its explicit discharge coefficient is one,
the ideal-nozzle assumption. These values define an integration example,
not a measured valve. The reservoir is finite: its mass, energy, temperature
and pressure change as it supplies the string's surroundings.

The optional `nozzle` row follows `gas`, before `conduction`. The first field
names a gas record relative to the scenario directory; both records use the
same explicit `--gas-component` catalog and frozen-phase policy. To inject
another species or mixture, change the reservoir record's IDs/fractions and
provide its species, caloric and transport packs through the existing repeated
flags. Material names supply no missing properties. The initial local mixture
need not contain the injected species. Both volumes retain the original source
laws and formation-energy references after mixing.

Each sample first advances the original solid/electrical/heat transaction,
then transfers upstream species and enthalpy through the native nozzle. Actual
post-heating pressures choose the direction; reversing the pressure difference
reverses transfer without a different simulation. The final surrounding
composition, density, viscosity and sound speed reach subsequent drag and the
same continuing pressure observer. Modal states and solid temperatures persist.
Nozzle closure is explicit zero area or zero discharge coefficient.

The nozzle and both parcel inversions reuse the authored `gas` temperature
bracket, specific-energy tolerance and iteration limit. The `nozzle` row adds
extensive mass/energy tolerances and a maximum donor-inventory fraction per
sample. Both initial temperatures must be in that bracket. Source-domain exit,
excess withdrawal, pressure-order overshoot or a failed combined energy check
refuses the trajectory; neither a hidden reservoir nor automatic clipping repairs
it. These limits do not estimate time-discretization error.

With flow enabled, columns after the existing 27 fields and any layer temperatures expose both
masses, both extensive internal energies `m*u`, reservoir temperature/pressure,
their combined internal energy, signed mass rate into the surroundings,
cumulative signed mass/enthalpy transfer, and each parcel's storage change.
Pairs of columns `gas_mole_fraction_ID` and `reservoir_mole_fraction_ID` report
the actual final mixture for every initial species in either volume. An absent
species has fraction zero. Positive transfer is into the string's surroundings;
advected enthalpy can itself be signed under a valid source reference.

The existing `storage_change_j` and cumulative energy defect include **both**
variable-mass gas inventories. Transferred enthalpy is internal to that pair,
so it is not added as external heat. The energy ledger uses the combined
thermal/flow allowance. The first 27 columns and existing layer-temperature
columns retain their order; omitting `nozzle` preserves the old output schema
and heat-only route.

The string, local parcel, reservoir and observer publish together only after
thermal, flow, electrical and energy admission. A late refusal returns no
partial trajectory. Electrical source work enters solid heat once; conductive
heat and transported enthalpy cancel internally in the combined energy account.
Uniform strings retain their original initial-state Joule quadrature; bonded
strings retain endpoint Joule quadrature. This remains first-order heat-then-flow
composition of well-mixed, fixed-volume inventories. It adds no resolved jet
momentum, spatial mixing, chamber resonance or radiation load.

## Explicit scenario

`heated.string` is a bounded, whitespace-delimited input. Records must appear in
this order, beginning with `frankensim-finite-gas-string-v1`; extra fields and
records are refused. All physical values use SI units.

| Record | Fields, in order |
| --- | --- |
| `audio` | sample rate [Hz], sample count, pressure at PCM full scale [Pa] |
| `string` | length [m], radius [m], fixed tension [N] **or** `fixed_support`, mode count |
| `release` | point-load station as a fraction of length, signed force [N], first sample without the load |
| `listener` | compact-observer distance [m] |
| `solid` | initial temperature [K], minimum and maximum caloric temperatures [K], source query angular frequency [rad/s], enthalpy interpolation allowance [J/kg], maximum temperature change per sample [K] |
| `gas` | fixed volume [m³], minimum and maximum energy-inversion temperatures [K], specific internal-energy tolerance [J/kg], maximum iterations, maximum temperature change per sample [K] |
| `nozzle` (optional, uniform or bonded) | reservoir gas-record path, reservoir volume [m³], flow area [m²], discharge coefficient [0,1], mass tolerance [kg], extensive energy tolerance [J], maximum donor fraction (0,1) |
| `conduction` | exchange area [m²], gas layer distance [m], solid layer distance [m] |
| `exchange` | maximum absolute transfer per sample [J], source residual tolerance [J], energy tolerance [J], maximum iterations |
| `drive` | `none`, `current A`, or `voltage V` |
| `heat` | signed external solid power [W], signed external gas power [W] |
| `thermoelastic` (optional alternative) | reference stress-free length [m], reference temperature [K], linear strain limit, maximum absolute reversible heat per sample [J], energy residual tolerance [J], maximum iterations |

After `heat`, `fixed_support` requires exactly one of the uniform
`thermoelastic` record or this bonded block:

| Record | Fields, in order |
| --- | --- |
| `bonded` | core radius [m], sleeve count (1–15), common stress-free reference length [m], reference temperature [K], linear strain limit |
| `layer` (one per sleeve, core outward) | source pack path relative to this scenario's directory, outer radius [m], contact area with preceding layer [m²], inner and outer half-path resistance geometry [1/m] |
| `layered_solve` | maximum thermal evaluations per inner solve, maximum Newton iterations, finite-difference energy step [J], maximum cell energy residual [J], thermal balance tolerance [J], combined mechanical/thermal tolerance [J], Newton absolute residual tolerance [J] |

The primary solid CLI argument is always the core. Sleeve paths cannot be
absolute or contain parent-directory components. Radii must increase, and the
last must equal the `string` radius; interface and gas areas must fit their
physical cylindrical surfaces. Each layer uses the declared solid temperature
range and maximum temperature increment. `exchange` still controls the outer
solid/gas solve; `layered_solve` bounds each native nonlinear thermal solve
within it. Up to `exchange.maximum_iterations + 2` inner proposals may be
needed. No uniform reversible-heat bracket is silently ignored: the bonded
block replaces the entire uniform `thermoelastic` record.

Numeric fixed tension forbids both thermal-support blocks. The original
numeric-tension and uniform fixed-support scenarios remain valid.

The parser admits 8,000–192,000 Hz, 1–262,144 samples and 1–32 modes. Input is
limited to 8 KiB. The exchange area may be zero for an insulated interface and
may not exceed the cylinder surface. Negative prescribed heat removes energy;
negative force reverses excitation. A signed electrical drive reverses current,
while resistance remains positive.

The gas record uses `thermochemical T P mole COUNT ID FRACTION ...`, or `mass`.
Repeat `--gas-component` for each explicitly paired species, NASA-9 and transport
pack triple in a mixture. The evolving parcel requires `--gas-phase frozen`
and an explicit water identity policy. It does not perform evolving phase
equilibrium.

Source support applies to every evaluated proposal. The supplied solid covers
280–400 K; the supplied nitrogen transport has a declared pressure interval of
75,000–101,325 Pa. Widening an inversion or exchange bracket does not extend a
source's domain. An oversized bracket can therefore refuse during an endpoint
probe even when the initial state is supported.

## Useful comparisons

Copy the scenario and change `drive current 0.1` to `drive none` to compare the
same mechanical release without Joule heating. `drive voltage 0.3` instead
exercises the changing current of a fixed-voltage source as resistance evolves.
Drive values are held over the run. The uniform path evaluates electrical work
from initial resistance; the bonded path uses accepted endpoint resistance.
Both report initial and endpoint electrical response.

Increase the first `gas` value from `0.000001` to `0.0001` m³ to enlarge the
thermal inventory while retaining initial temperature, pressure and composition.
This changes the gas temperature response and subsequent heat exchange with the
same vibrating string. The volume does not introduce a chamber resonance.

Changing the viscosity curve in `properties.tsv` and rebuilding the solid pack
changes material loss. Changing the gas record and all corresponding source
packs changes the surrounding species. No executable recompilation is needed
for these material, gas, geometry, drive or excitation changes. Existing source
domains and the declared step controls still govern admission.

## Pressure and energy output

Each CSV row describes one accepted sample. The first 22 columns retain their order:
time, physical observer pressure, point displacement and velocity, both
temperatures, gas pressure, cumulative material heat, drag heat, signed
solid-to-gas exchange, Joule heat, prescribed solid and gas heat, mechanical
work, actual stored-energy change, cumulative energy defect and tolerance,
initial and endpoint current, initial and endpoint resistance, and first-mode
damping ratio. Five appended columns report actual tension, first-mode frequency,
cumulative signed reversible heat, axial free-energy change and mechanical
entropy-energy change. With `drive none`, the four electrical response columns
and Joule heat are zero. The three reversible-storage columns are zero for the
original fixed-tension model.

Stored energy is measured directly from accepted modal states, source-resolved
axial storage and mechanical entropy, solid specific enthalpy, and gas specific
internal energy. Thermoelastic mechanical internal energy is kinetic energy
plus axial/transverse free energy plus `T*S`, with `S = -d(Psi)/dT` at fixed
modal coordinates and support length. Its change is relative to the
initial loaded string and initial thermal states. Material dissipation, drag
and conductive transfer are internal exchanges; they cancel from the combined
balance. The cumulative defect compares the three actual storage changes plus
declared other losses with mechanical work, Joule work and prescribed heat.
Numerical residuals are not deposited as physical heat.

The WAV maps the specified `full_scale_pa` to PCM16 full scale and reports any
clipping. It performs no peak normalization. The CSV retains the physical
pressure before PCM scaling or clipping, so amplitude comparisons across runs
remain meaningful.

This is a first-order, frozen-law-per-sample, small-displacement string model
with fixed geometry/density and either one uniform temperature or one temperature
per declared bonded annulus. The original
fixed-tension path keeps Young's modulus constant. The fixed-support path
admits scalar thermal expansion and globally affine Young's modulus with the
corresponding reversible heat. Slack/compression, large strain, moving supports,
radial contraction and melting refuse or require a different model. The
well-mixed, fixed-volume gas is a thermal inventory. The optional nozzle on either string path
uses the existing upstream-enthalpy flow work account; moving-volume pressure
work, resolved jets, chamber acoustics and phase change remain outside this example.
The compact pressure observer supplies no reciprocal radiation load.

Focused tests in `tests.rs` cover strict signed inputs, continuous motion
and deterministic pressure, an independent three-storage energy integral,
material/species/inventory/drive changes reaching the same output, sourced
thermal expansion reaching tension, pitch and continuing pressure, and source,
excess-heat and cancellation refusal. Bonded tests additionally check the same source-driven
trajectory, independent section tension/frequency and parallel resistance,
endpoint Joule work, per-layer temperatures, and data-only sleeve/gas swaps.
The nozzle regressions decode canonical source-pack bytes and both gas records,
check species inventories and an independent NASA caloric integral, compare
mixture data reaching pressure/drag, and cover closed-path parity, outward flow,
invalid inputs and missing/mismatched finite reservoirs.
Run them through the example target:

```bash
cargo test -p fs-couple --example finite_gas_string
```
