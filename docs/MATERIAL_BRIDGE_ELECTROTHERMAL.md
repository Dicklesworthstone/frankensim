# Finite electrical supply, material heat and spatial deformation

`fs_couple::strip_enthalpy::bridge::electrical::step_candidate` connects the
existing source-pinned material network, charge/flux circuit, deposited beam and
enthalpy network in one nonlinear thermal residual. A material's resistivity can
change current and Joule heat; that heat changes its relaxation and sag. The
resulting mechanical dissipation heats its original cells as well. Conductivity
and resistivity are evaluated at the same proposed endpoint, not staggered into
independent accepted thermal steps.

## Original owners and one transaction

The adapter invokes the original `circuit_enthalpy::step_candidate`. Its material
callback evaluates the original `Bridge::advance_with_cell_energy`, deposits
physical Maxwell heat into the mapped thermal cells, and returns the original
material-derived electrical and thermal coefficients. The circuit owner performs
its ordinary Kirchhoff/storage update and branch heat partition. The same
enthalpy residual includes both generated heats, conduction and authored external
boundaries. No second circuit integrator, beam solver, constitutive law, Joule
partition or Newton solver exists in this adapter.

Every trial starts from the same accepted charge/flux and Maxwell histories. The
final fresh trial carries its own beam candidate through the circuit/thermal
acceptance, so the returned mechanical state is not a replay at another guessed
temperature. Cancellation, source-domain failure, either work budget, or a final
combined-energy refusal returns none of the three endpoints.

The prepared network exposes `resistor_destinations()` from its original cell
bindings. This avoids a second caller-authored heat map that could send a branch's
loss to a different material. Both thermal-conductivity snapshots and
**electrical-only conductor snapshots** are checked against the mechanical source
card and fixed coordinates. A missing conductivity requirement on an isolated
cell does not bypass electrical/mechanical source compatibility.

## Usage

Prepare the normal `MaterialNetwork`, including its explicit `CellResistor` and
`ConductionLink` bindings, and the original `MaterialCircuit`. Capacitor charges
and inductor fluxes belong to that circuit; the adapter does not manufacture an
ideal supply to replace them.

```rust,ignore
use fs_couple::strip_enthalpy::bridge::electrical::{self, ElectricalConfig};

let transport = network.prepare(&thermal, cancelled)?;
let candidate = electrical::step_candidate(
    &thermal, &beam, &bridge_bindings, &circuit, &transport, &held_ports,
    target_s, gravity_m_s2, &additional_downward_n_m,
    ElectricalConfig {
        bridge: beam_feedback_controls,
        circuit: circuit_work_controls,
        maximum_circuit_trials,
    },
    &energy_guesses_j, cancelled,
)?;
let (heat_next, beam_next, circuit_next, report) = candidate.into_parts();
```

Keep the prepared network across accepted steps to preserve selected property
claims. The original `BridgeEnthalpyConfig` bounds mechanical trials and beam
regime. `maximum_circuit_trials` separately bounds every electrical evaluation,
including discarded Jacobian stencils and final admission. Energy guesses are
joules in sorted thermal CellId order; line loads are in bridge axial order.
Authored boundary heaters, cooling and reservoirs remain active.

`report.electrical` is the original circuit step, including branch currents,
material resistors, storage and supplied work. `report.mechanical` is the original
beam step, including reactions, displacement bounds and independent load work.
Original material samples, per-layer energy and the numerical report remain
available. `thermal.external_heat_j` includes both internally generated heats;
do not add `deposited_joule_j` or `deposited_mechanical_heat_j` again.

`other_external_heat_j` is evaluated independently from the original boundary
rate owner before adding mechanical or electrical heat. The final gate checks

```
delta electrical storage + delta elastic storage + delta thermal storage
+ mechanical numerical dissipation
- electrical supply work - beam load work - other external heat = 0.
```

Mechanical numerical damping is not physical heat. Gravity is included once in
beam load work, not additionally as a thermal source or a second potential-energy
subtraction. This total gate is independent of the original component gates and
also checks the sum of both deposited heats against actual thermal input.

## Evidence and model limits

Six G0/G1/G3 Rust regressions cover analytic finite-capacitor heat with either
charge sign, electrical-data-only changes to sag, endpoint resistance with
load/unload memory and independent owner replay, unloaded circuit-path parity,
electrical-only source mismatches/domain refusal, and final cancellation/budgets.
They have been authored, not executed, in the editing environment:

```sh
rch exec -- cargo test -p fs-couple --test bridge_electrothermal \
  --test bridge_enthalpy --test bridge_enthalpy_growth --test circuit_enthalpy
```

DSR, its documented fallback, RCH, Cargo and rustc are unavailable locally.
Independent NumPy/SciPy equations passed 216 parameter cases plus four accepted
load/unload steps, with maximum energy discrepancy 3.84e-15 J. For an explicitly
synthetic specimen, changing only resistivity from 0.1 to 0.2 ohm m changes
endpoint temperature from 303.921209 to 301.980171 K and downward tip displacement
from 0.880247 to 0.867685 mm. These are equation checks, not Rust output or
calibrated material predictions. Whitespace and reviewed-blob checks do not
establish Rust compilation or experimental validity.

The original quasistatic slender/small-strain beam and caloric restrictions
remain: constant reference density, temperature-independent spring energy and
free strain, one solid phase, fixed thermal geometry. Electrical C/L, conductor
path geometry and topology remain fixed during a step. Temperature-dependent
R, k and relaxation are admitted only within their source support. The electrical
midpoint and thermal/mechanical endpoint combination is first-order in general;
nonlinear residual convergence is not a time-accuracy or unconditional-passivity
claim. No Lorentz force, field-gap motion, piezoelectric response, thermoelectric
heat, melting or automatic circuit-topology migration is inferred. Series
contacts require their separate heat-routing owner and are rejected by this
bulk-only circuit entry point. Split at events and load/source discontinuities;
target-time events stay pending until joint acceptance.
