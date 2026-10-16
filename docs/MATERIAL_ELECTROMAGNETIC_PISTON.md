# Reciprocal circuit force, finite material heat and shared gas motion

`cavity_phs::piston::double_acting::electrical` connects the existing
`MaterialCircuit` directly to the original two-gas piston. Its `network` module
also connects the same circuit current and motion to the existing finite
material/enthalpy network. This extends `MATERIAL_PISTON_CIRCUIT.md`: electrical
energy can now do mechanical work directly, not only through heating the gas.

## One reciprocal connection, original physical owners

A dedicated existing voltage-source branch represents a frozen linear
translational transducer. With positive branch current from first to second
terminal and positive velocity expanding the primary gas:

```text
back_emf = k * midpoint_velocity
force    = k * branch_current
native_source_flow = -branch_current
```

Thus the native electrical port work plus mechanical force work cancels when
`dx = dt * midpoint_velocity`. The single signed factor `k` has units
N/A = V s/m. Changing winding orientation reverses it. This is the linear
moving-conductor/voice-coil approximation (uniform transverse field gives
k = B * active_wire_length). `k` is explicit field/geometry or measured assembly
data, not inferred from a metal name or replaced by separate arbitrary force
and voltage constants. It is frozen for this step; a changing magnetic field or
position-dependent stored magnetic energy is not implemented.

`TranslationalCircuitPort.branch` names the ORIGINAL graph branch index. That
branch must be a voltage source with its own exclusive source-port index.
`ElectricalMotionInput.inputs` supplies every circuit input but must contain
zero at this reserved port. Applied electrical supplies belong to other ports.
The adapter refuses shared source-port indices and nonzero placeholders rather
than overwriting a user's input. Zero `k` explicitly disconnects mechanics but
leaves the voltage-source branch as an electrical short, not an open circuit.

The circuit still owns charge, flux, KCL/KVL and material resistor heat. Each
gas still owns its complete caloric/EOS/transport/acoustic state. The original
piston owns pressure work and one momentum. The existing fallible Newton solver
closes velocity feedback; it does not replace either physical time integrator.
Every numerical probe starts from unchanged accepted histories. Fresh velocity,
kinematic and work checks run before returning either endpoint.

## Finite material network

`MaterialDoubleActingGasPiston::electromagnetic_network_candidate` retains the
original source-resolved rigid body's canonical inertia while also returning
all thermal cells and the circuit. There is no second mass parameter or source
resampling when the complete mechanism is returned.

At each proposed thermal endpoint, `MaterialNetwork` supplies actual resistance,
conductivity and phase data. The electrical-motion solve uses that resistance
for the SAME current that supplies both electromagnetic force and Joule heat.
Its mechanical callback solves both gas heat contacts and common motion. The
original `circuit_enthalpy` residual deposits the resulting resistor heat into
its existing destinations; no second circuit solve follows the force solve.
The ordinary circuit/enthalpy API delegates to that same private implementation.

The two gas contacts name finite cells through `DoubleNetworkDrive`. When both
name one cell, heats are grouped into one debit and total contact area must fit
its exposure. Both original same-card conductivity evaluations are retained.
Missing properties still refuse a zero-area contact; ambient convection on an
occupied gas-contact cell is refused rather than duplicated. Other cells keep
their existing source, reservoir and conductive-link boundary laws.

The independent combined balance is:

```text
Delta E_circuit + sum(Delta H_cells) + Delta U_a + Delta U_b + Delta K
  = work of OTHER electrical ports + external network heat
    + prescribed gas heats + additional applied force work
```

Joule heat, gas/cell contact heat, and electromagnetic transfer are internal.
External thermal supply is evaluated by the original rate law BEFORE adding
internal transfers. The implementation checks actual source-map debits, native
port work, both gas heat residuals, local cell residuals and combined energy.
Budget exhaustion or final cancellation returns no partially updated owner.

## Calling pattern

After constructing the original material cards, circuit graph, active thermal
cells, shared mechanism and explicit controls:

```rust
let proposal = mechanism.electromagnetic_network_candidate(
    &materials, &thermal, target_s, &circuit, &circuit_inputs,
    drive, network_controls, &initial_cell_energy_increments_j,
    || false,
)?;
let (thermal, circuit, mechanism, report) = proposal.into_parts();
```

`drive` is `electrical::network::ElectromagneticNetworkDrive`. It contains the
existing paired gas/thermal controls and geometry, the reserved transducer
branch and factor, and `ElectricalMotionSolve`. Its electrical tolerances must
match `network_controls.network.circuit`. Initial cell guesses are in sorted
CellId order; the velocity guess is explicitly supplied and never overwritten
with a discarded trial's history. All nested circuit and boundary solves count
against global work limits, including numerical Jacobian probes.

The standalone `DoubleActingGasPiston::electrical_candidate` exposes the direct
connection without thermal cells. In that method resistor heat is a reported
output, NOT an evolved finite thermal store. Use the network method to retain it.

## Evidence and limits

Synthetic reference comparisons test motoring, generation/reaction force,
resistance changes, winding reversal, other electrical supplies and refinement.
Native regression sources use the original Rust circuit, material and gas owners,
including finite heat destinations, same-cell coverage and final cancellation.
Reference equations do not establish native Rust execution or measured accuracy.

The network couples endpoint thermal coefficients with midpoint electrical and
gas mechanics and is generally first order; refine time separately from solver
tolerances. All media are uniform closed frozen-composition gases. Circuit C/L
and topology, thermal geometry and caloric charts remain fixed. The mechanical
body is rigid with frozen material state; heating network cells does not silently
soften, melt or expand it. There is no magnetic field solve, saturation,
variable-inductance work, eddy/skin loss, switching, combustion, mass exchange or
emitted-audio model. Updated acoustic material states are exposed, not rendered.
