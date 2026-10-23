# Material vibration exchanging heat with a finite gas

The existing `self_heating_vibration` executable evolves one source-bound axial
material specimen, including inertia, Maxwell branch memory, physical dissipation
and temperature-dependent relaxation. It now accepts a finite gas reservoir.
Heat exchanged with that reservoir changes its actual internal energy,
temperature, pressure and transport state. Those changes feed back into the
same material/thermal solve within each proposed step.

The original command without `--finite-gas` retains its infinite 300 K
environment and original CSV columns. The finite mode uses
`fs_couple::finite_gas_thermal::FiniteGasThermalEnvironment` inside the original
`Oscillator::step` feedback callback. It does not freeze the gas at its initial
temperature or run a separate material oscillator.

## Run from existing source data

The specimen table is a **manufactured integration fixture**, with deliberately
small heat capacity to expose feedback. It is not measured metal or polymer
data. The following gas uses the repository's nitrogen species association,
NASA-9 caloric coefficients and explicitly paired dilute-gas transport source.
The native loader is shared with `music_render` and `material_strings`.

From the repository root:

```bash
run_dir="$(mktemp -d)"
cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/nitrogen/manifest.tsv \
  --out "$run_dir/nitrogen.species"
cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/nitrogen-nasa9-mcbride2002/manifest.tsv \
  --out "$run_dir/nitrogen.nasa9"
cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/nitrogen-transport-svehla1995/manifest.tsv \
  --out "$run_dir/nitrogen.transport"

cargo run --release -p fs-couple --example self_heating_vibration -- \
  --table crates/fs-couple/examples/self_heating_vibration/material.txt \
  0.5 0 0.04 800 100 \
  --finite-gas 1e-8 200 600 1e-5 \
  --gas-record crates/fs-couple/examples/self_heating_vibration/nitrogen.gas \
  --gas-component "$run_dir/nitrogen.species" "$run_dir/nitrogen.nasa9" "$run_dir/nitrogen.transport" \
  --gas-water none --gas-phase frozen \
  > "$run_dir/finite-gas.csv"
```

Use the workspace's DSR/RCH execution policy for native builds and tests.
These commands describe the actual executable interface; source review alone
does not establish a completed native run.

The positional values are initial end velocity [m/s], payload [kg], duration
[s], step count and declared film coefficient [W/(m² K)]. The original
`--pack FILE BRANCH_COUNT` material ingress is also supported.

`--finite-gas VOLUME_M3 TMIN_K TMAX_K MAX_EXCHANGE_J` declares the fixed gas
volume, its source-energy inversion bracket and the symmetric signed exchange
bracket for each step. Every evaluated proposal must be supported by the
original material/gas sources. A large exchange bracket is not a harmless
accuracy knob: an unsupported endpoint probe refuses. The temperature interval
also bounds the permitted gas temperature increment.

Two optional inputs are available:

- `--gas-power W`: signed external thermal power supplied to the gas,
  default zero. This is separate from heat exchanged with the solid.
- `--solid-temperature K`: initial solid temperature, default 300 K.
  It must lie in the original material chart.

The gas declaration uses the shared grammar:
`thermochemical T P mole COUNT ID FRACTION ...`, or `mass`.
Repeat `--gas-component` for each explicitly paired species, NASA-9 and
transport source triple. Species names never select missing coefficients.

The evolving parcel requires `--gas-phase frozen`. An initial
`liquid-water` saturation check would not establish phase validity throughout
later heat evolution, so it is refused rather than silently downgraded.
The supplied nitrogen transport has an authored dilute application rectangle
of 75,000–101,325 Pa; the gas pressure must remain within it. The temperature
bracket does not broaden that source support.

## Physical comparisons

Change only the gas volume from `1e-8` to `1e-4` m³ to increase thermal
inventory while retaining the same initial temperature, pressure, composition,
specimen, excitation and declared film coefficient. The larger reservoir
approaches the original fixed-ambient experiment; it is still finite.

Add `--gas-power 0.001` and set the initial velocity to zero to heat the solid
through the gas. The reported solid-to-gas transfer becomes negative. Change
only the gas composition or its caloric source to investigate finite heat
capacity effects without changing the mechanical source.

The film coefficient remains explicitly declared. Changing gas transport does
not automatically infer flow, a new convection correlation or a new surface
coefficient. A correlation-driven interface requires its own admitted geometry
and flow model.

## Accepted history and energy

The finite mode reports physical time, displacement, velocity, both
temperatures, gas pressure, cumulative material heat, numerical loss, signed
solid-to-gas transfer, external gas heat, radiation heat, combined stored
energy, and independently measured cumulative energy defect/tolerance.

One Q leaves the solid and enters the gas. Physical Maxwell loss becomes
solid heat exactly once. Backward-Euler losses remain numerical losses and
are never deposited as heat. The cumulative balance measures storage directly
from accepted mechanical, solid-enthalpy and gas-internal-energy states;
internal transfers cancel. Radiation is disabled by this example's original
source adapter and is zero in its output.

Source identities and initial gas composition are printed to stderr. CSV
publication starts only after the entire trajectory has been accepted, so a
late source, convergence or cancellation refusal emits no partial trajectory.

This is a fixed-geometry, small-strain, single-solid-phase axial reduction
with a uniform body and well-mixed fixed-volume gas. It does not add acoustic
radiation, pressure work, gas motion, a chamber resonance, a resolved thermal
boundary layer or phase change. Step refinement and independent physical
validation remain necessary.

Focused regression definitions are in `finite_gas_tests.rs`: finite-volume
and caloric-source changes reaching actual motion; an independent NASA heat
capacity integral and whole-system energy; reverse transfer from external gas
heating; deterministic replay; source/temperature/cancellation refusal.
