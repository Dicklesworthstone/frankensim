# Source-driven layered sections: electricity, heat and sound

`fs_material::section_temperature::ConcentricSectionTemperatureLaw` retains
one original material card and source-claim selection per layer. A temperature
update re-evaluates every requirement in the original bundle, not only Young's
modulus. Mechanical, thermal, electrical and auxiliary properties keep the same
source selection and non-temperature coordinates. Unsupported states refuse.

`section_temperature::heating::HeatedConcentricSection` connects those updates
to the existing `fs-thermal` enthalpy timeline. Each annulus retains its original
mass and its source-backed Cp/density chart. Every thermal substep uses one
material evaluation for per-layer Joule heat and link conductance, then refreshes
the endpoint mechanical section. A failed advance returns no candidate. Generated
heat does not become a persistent authored source, and the thermal audit already
includes it: never add the electrical-energy report to that audit again.

## Run the actual acoustic consumer

From the configured FrankenSim workspace:

```sh
cargo run -p fs-couple --example heated_section_sound -- \
  --demo 2.0 0.2 1.0 0.01 comparison
```

The arguments are sleeve electrical-resistivity scale, terminal voltage [V],
heating duration [s], PCM full scale [Pa], and output prefix. The core is unchanged.
The demo is synthetic data, not a claim about gold, silver, steel or a measured
specimen. Repeat with another sleeve scale to change data without editing physics.

For real input, use two existing canonical `NormalizedMaterialCardPack` files:

```sh
cargo run -p fs-couple --example heated_section_sound -- \
  --packs core.pack sleeve.pack 0.2 1.0 0.01 comparison_from_packs
```

The example's declared source context is `T` in kelvin over [300,400] K; it
requires positive density, Young's modulus, specific heat, thermal conductivity
and electrical resistivity. Density must be constant on that interval. Missing,
ambiguous, differently conditioned or unsupported data refuse; no nominal-metal
fallback is used. The library API additionally supports full caller-supplied
nonthermal context and per-layer absolute-temperature axis names.

The same section geometry, fixed 10 N tension, pluck, ambient gas, listener and
acoustic renderer produce the two signals. `comparison.csv` contains physical
pressure samples for cold and heated states. `comparison.wav` is a stereo
comparison, **cold on the left and heated on the right**, using the existing
PCM encoder with one shared explicit full scale. No separate peak normalization
hides material effects. Clipping is reported. Existing output files are not
overwritten. The example also reports layer temperatures, axial stiffness,
integrated electrical heat, stored-energy change and the thermal residual.

These are two frozen-state 50 ms sound observations, not continuous sound during
heating. The model is a perfectly bonded concentric slender member, not a helical
winding model. Each annulus has one temperature; interface area and half-path
thermal resistance geometry are explicit modelling inputs. This rung does not
include thermodynamic expansion/work feedback, density-driven geometry changes, phase
transitions, skin effects, or a resolved intralayer temperature field. The existing
conduction CFL does not certify nonlinear electrical-source stability; independently
refine the requested maximum step. Source/model uncertainty remains upstream.

Focused native checks (not executed in the connector-only implementation session):

```sh
cargo test -p fs-material section_temperature
cargo test -p fs-couple --example heated_section_sound
```

## Constrained heating: retain the physical supports, not the old pitch

The optional `HeatedConcentricSection::with_axial_constraint` binding adds
source-pinned thermal expansion and per-layer stress on the SAME internal
thermal steps. Each original layer bundle must explicitly include
`linear_thermal_expansion_coefficient` in 1/K, with finite signed admission.
The shared source integrator evaluates the full alpha(T) path, not endpoint
alpha times temperature change. This expansion rung uses canonical `T` and a
common original stress-free reference temperature. Missing, unsupported or
incompatible expansion data refuse rather than becoming zero expansion.

```rust
use fs_material::section::thermal::AxialConstraint;
use fs_material::conductor::OhmicDrive;
use fs_couple::material_assembly::concentric::thermal::heating::
    compile_heated_concentric_material_assembly;

// `body` is a configured HeatedConcentricSection; alpha is already sourced.
let constrained = body.with_axial_constraint(
    AxialConstraint::FixedStrain(1e-4), 0.01,
)?;
let (heated, audit) = constrained.advance_to(
    target_s, OhmicDrive::Voltage(voltage_v), step_config, checkpoint,
)?;
let state = heated.axial_state().expect("constraint was attached");
let layer_stresses = state.response().layers();
let actual_support_force_n = state.response().force_n();
let sound = compile_heated_concentric_material_assembly(
    &template, &heated, &other_material_bindings,
)?.realize()?;
```

The acoustic adapter borrows the accepted thermal specimen and uses its actual
support resultant with the same mass and stiffness, without retuning or restoring
the initial tension. It requires matching physical length and stationary
fixed-strain supports. A slack or compressed result remains a valid generic
material response, but refuses at the taut-string boundary. Existing damping,
source and acoustic admission checks remain active.

Generic rods also support `FixedForce` (zero means free expansion) and
`LinearSpring`. Differential expansion can leave opposing layer stresses even
at zero resultant; every layer's total, thermal and elastic strains is checked.
A failure at an intermediate thermal substep discards the complete proposed
thermal and axial advance. Attaching supports after heating retains the original
stress-free reference rather than resetting it at the new temperature.

This is explicitly ONE-WAY thermal-to-quasistatic mechanics. Elastic energy is
reported as a diagnostic, not added to heat or presented as a closed total-energy
balance. Thermal geometry/reference mass stay fixed; force/spring-controlled
strain does not move that mesh and those controls are not admitted by this
fixed-end acoustic observer. Buckling, thermodynamic backreaction, finite strain
and continuous sound during heating remain outside the model. The CLI example
above retains its original fixed-force comparison; the constrained path is the
library API shown here.

Focused additional native check:

```sh
cargo test -p fs-couple --test heated_constrained_section
```
