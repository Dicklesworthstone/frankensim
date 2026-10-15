# Finite electrical storage, material heat and mechanical response

`fs_couple::material_circuit::MaterialCircuit` connects existing source-resolved
conductors to the original `fs_phs::circuit` Kirchhoff descriptor solver. It
retains capacitor charges and inductor fluxes while every resistor is refreshed
from a `ResolvedConductor`. The original graph owns C/L, ports, transformers and
topology. No new RC/RLC integrator, circuit solver or resistivity law is used.

Every resistor branch must have exactly one material binding. Numeric resistor
values on the input graph are replaced, not used as unsourced fallback values.
A step returns each branch's signed voltage/current and nonnegative material
heat. These use the actual algebraic midpoint of the native descriptor step;
the static Hamiltonian gradient is zero on its constraint coordinates and would
incorrectly report no current or heat there. Algebraic multipliers are not
stored-energy history. Only charge and flux survive material re-evaluation.

Branch electrical work, total material heat, native dissipation and electrical
storage/supply balance are independently compared before returning a candidate.
Closed graphs take no external inputs: every joule of heat then depletes finite
stored electrical energy. Driven graphs retain external port work separately.
Original branch indices, complete source bundles and native step reports survive.

## Coupled finite supply

`GasPiston::circuit_thermal_candidate` passes a circuit candidate through the
existing state-dependent power and finite-body heat-exchange solve. The body
resistivity is resolved at the SAME enthalpy/phase midpoint as heat transport.
All trial circuits start from the same accepted charge/flux: iteration cannot
drain the supply repeatedly. The final audit includes electrical storage, body
enthalpy, gas internal energy, piston kinetic energy, exterior pressure work,
applied force work, both prescribed heats and external electrical port work.

`MaterialOhmicBody` supplies the source/geometry binding without an independent
ideal current/voltage. The existing `MaterialJouleHeating` ideal-drive API uses
that same conductor resolution, including same-card identity, thermal-coordinate
updates, static-law admission and whole-body volume checks.

```rust,ignore
use fs_couple::material_circuit::{MaterialCircuit, MaterialResistor};
use fs_couple::material_circuit::thermal::CircuitThermalSource;
use fs_couple::cavity_phs::piston::thermal::power::joule::MaterialOhmicBody;
use fs_phs::circuit::{Branch, CircuitGraph};

// ohmic is a MaterialOhmicBody bound to the same original material card,
// source conditions, length and area as the finite thermal body.
let conductor = ohmic.resolve_conductor(body.initial, body.volume_m3)?;
let circuit = MaterialCircuit::new(CircuitGraph {
    node_count: 2,
    branches: vec![
        (1, 0, Branch::Capacitor { farads: 0.002 }),
        (1, 0, Branch::Resistor { ohms: 1.0 }), // replaced by material binding
    ],
    transformers: vec![],
}, vec![MaterialResistor { branch: 1, conductor }],
   vec![], vec![0.2], &[])?; // authored 0.2 C gives 10 J, not infinite supply

let candidate = piston.circuit_thermal_candidate(
    input, body, piston_controls, exchange_controls, source_controls,
    CircuitThermalSource {
        circuit: &circuit, resistor_branch: 1, body: ohmic, inputs: &[],
    }, circuit_controls, conductance_callback, || gate.is_requested(),
)?;
let (next_piston, next_body, next_circuit, report) = candidate.into_parts();
// Publish ALL endpoints together in the application's physical state owner.
// Resistor heat and interface heat are already included. Do not apply again.
```

The generic circuit path supports multiple resistors and exposes their separate
heat. The gas/piston convenience path supports exactly one whole-body resistor;
additional resistive sinks refuse rather than being silently deposited into the
wrong body. A multi-body circuit needs a simultaneous thermal mapping of those
branches. The circuit source may have multiple capacitors, inductors, sources
and transformers within the explicit dense-state bounds.

## Physical and numerical limits

C/L and circuit topology stay fixed; their values remain explicitly authored
circuit data. They are not yet temperature-dependent material dielectric or
magnetic laws. Resistance can depend on the shared temperature and declared
liquid mass fraction. A latent plateau may change phase, resistance and supply
depletion without changing temperature. Missing data are not a switch model.

The body has fixed geometry and a uniform enthalpy state. No spatial current
profile, moving induction geometry, electromagnetic mechanical force, contact
arc, switching topology, dielectric loss or flowing melt is implied. Midpoint
integration does not prove monotone charge decay or thermal equilibrium approach
at arbitrarily large timesteps. Bounds and energy residual checks are not a
constitutive qualification or a time-discretization error certificate.

Six circuit integration tests and six coupled finite-supply tests exercise the
original owners, charge/flux history, resistor partitions, transformer/source
work, resistivity-only trajectories, latent-phase feedback and atomic refusal.

```sh
cargo test -p fs-couple --test material_circuit
cargo test -p fs-couple --lib cavity_phs::piston::thermal::power::joule::tests::circuit
```

Independent synthetic NumPy/SciPy equations also check linear circuit refinement
and combined electrical/thermal/gas/mechanical energy. They are not native Rust
execution or measured-material validation.
