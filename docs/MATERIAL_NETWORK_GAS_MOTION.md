# Finite material networks drive gas and mechanical motion

`GasPiston::material_network_candidate` connects the existing finite circuit /
thermal network to the existing closed gas and translating inertial wall.
Multiple resistor bodies retain their separate caloric and transport sources.
One explicitly selected cell exchanges heat with the gas; all cells can still
conduct through the network. Electrical storage, material enthalpy and gas /
mechanical history are accepted together.

The complete feedback is: source-bound resistance -> circuit current -> heat in
the correct cell -> cell temperature or phase -> heat transfer into the gas ->
gas pressure and boundary motion -> new gas midpoint and heat rate. Capacitors
and inductors are finite energy stores, not extra prescribed heating terms.

## Existing owners, one coupled interval

The circuit/thermal network's original Newton solve still uses cell-energy
increments in sorted `CellId` order. At each proposed network endpoint the
material owner resolves the original resistors and thermal links. The new gas
boundary solves `Q = dt G (T_cell_endpoint - T_gas_midpoint)`, starting from the
same accepted gas/velocity for every trial. It delegates caloric inversion,
pressure work and motion to the existing `GasPiston::step`, and scalar heat
feedback to the existing bounded solver. The thermal source receives **-Q/dt**
at exactly the selected cell. Circuit branch heating remains deposited once by
`circuit_enthalpy`.

A final fresh network evaluation retains the matching gas candidate and exact
interface material record. Original local thermal, total thermal, electrical,
impulse, gas-energy and source residual checks remain. A separate all-owner
balance includes electrical storage, cell enthalpy, gas internal energy,
piston kinetic energy, exterior pressure work, external force work, electrical
port supply and independent prescribed heat. Internal Joule and gas exchange
are removed from external supply, not counted twice.

Cancellation, unsupported source queries, an inadmissible finite-difference
probe, exhausted work or a late energy failure returns no physical endpoint.
No nonlinear trial consumes accepted capacitor charge or changes gas volume.
The global boundary budget includes numerical probes and final admission;
each boundary solve also retains its bounded heat / motion / material work.

## Data-only use

The original `MaterialNetwork` supplies each body's caloric curve, card,
resistor geometry and link geometry. The selected gas-contact cell must have
an explicit conductivity binding, even for zero transfer area. No additional
property callback or metal / atmosphere preset is required.

```rust,ignore
use fs_couple::cavity_phs::piston::thermal::boundary::network::{
    NetworkGasDrive, NetworkGasSolve,
    material::GasConductionGeometry,
};

// Existing network, thermal history, circuit history, gas piston and numerical
// controls. Use the same physical duration on both sides; split at events.
let drive = NetworkGasDrive {
    cell: gas_wall_cell,
    input: piston_input,          // contains separate prescribed gas heat
    piston: piston_controls,
    exchange: heat_controls,
};
let controls = NetworkGasSolve {
    network: circuit_thermal_controls,
    maximum_boundary_solves: 300,
    total_energy_tolerance_j: 1e-6,
};
let geometry = GasConductionGeometry {
    area_m2: 0.01,
    gas_distance_m: 0.001,
    cell_distance_m: 0.02,
};
let proposal = gas.material_network_candidate(
    &network, &thermal, target_s, &circuit, &ports,
    drive, controls, geometry, &energy_guess_j, || false,
)?;
let (next_thermal, next_circuit, next_gas, report) = proposal.into_parts();
// Publish all three. Heat is already in each endpoint: never add it again.
let pressure_pa = next_gas.gas().state().state().pressure;
let sound_speed_m_s = next_gas.gas().state().state().sound_speed;
```

`MaterialConductionLayers` remains the conduction-law owner. Its two coefficients
come from the gas's actual midpoint and the body's actual network endpoint.
The body source must match the network's source bundle exactly, preserving its
claim selection, other coordinates and receipts. Liquid fraction changes use
the existing phase chart, including constant-temperature latent plateaus.
There is no second thermal store attached to the gas contact.

The material adapter admits a patch no larger than the cell's exposed area
after existing link and reservoir coverage, and requires zero ambient convection
on that cell. It does not silently expose the gas patch to ambient a second
time. Other cells retain their ordinary boundaries. The generic
`GasPiston::network_candidate` permits a caller-owned source and interface law;
surface coverage and any additional source accounting then remain explicit
caller obligations. `thermal_boundary_candidate` is also usable separately for
a prescribed reservoir or within another finite-source residual.

## Physical and numerical scope

This is a fixed thermal grid, fixed circuit topology / capacitance / inductance,
one uniform closed frozen-composition gas and one translating inertial wall.
The thermal cell is a chamber wall or other heat store, not necessarily the
translating wall. It does not gain mechanical deformation from its thermal or
phase state. There is no combustion, mass flow, moving melt, spatial convection,
multicylinder linkage or automatic reconfiguration of an existing acoustic
renderer. The updated common gas state exposes pressure, transport and acoustic
properties through its existing interface.

Network rates use backward Euler and endpoint material properties, while gas
energy and motion retain their midpoint method. The combined method is generally
**first order**, not second order merely because the gas owner is midpoint.
Energy balance and converged algebraic residuals are not time-error bounds,
uniqueness proofs or unconditional passivity claims. Refine time independently.

Focused native regressions (already included by the piston test selection):

```sh
cargo test -p fs-couple --lib cavity_phs::piston::thermal::boundary:: -- --test-threads=1
```

Seven authored Rust tests cover exact isolation parity, original source records,
heating/cooling motion, data-only changes through two bodies and finite supply,
external-source separation, all-owner energy, late cancellation, work exhaustion,
clock mismatch and surface admission. Local native execution was unavailable
during authoring; independent Python simultaneous-equation checks are not a
Rust test pass or measured-material qualification.

## Open flow through the same finite network

`GasPiston::nozzle_network_candidate` now composes the original conservative
nozzle transfer with this existing circuit/enthalpy/heat/motion solve.
`nozzle_material_network_candidate` supplies the data-only material ingress,
including an explicit conductivity or resistivity basis for every electrical
branch. `MaterialGasPiston` exposes the same operation while retaining its
original regional body, source receipts and material-derived inertia.

The nozzle exchanges species and upstream total enthalpy **once**, at the
accepted initial pressures. Every subsequent network trial starts from that
same transferred chamber. The reservoir is never withdrawn again during
Newton or heat-feedback iteration. Its new mixture therefore reaches the
actual gas pressure, conductivity, heat capacity and acoustic state used by
the wall calculation. Resistor and wall conductivity still share the prepared
network's original material bundle at each actual thermal endpoint.

```rust,ignore
let candidate = piston.nozzle_material_network_candidate(
    &reservoir, nozzle, &binding, &thermal, target_s, &circuit, &ports,
    &electrical_bases, drive, nozzle_controls, network_controls,
    wall_geometry, &energy_guess_j, || false,
)?;
let (next_thermal, next_circuit, next_piston, next_reservoir, report) =
    candidate.into_parts();
// Publish all four endpoints together; every internal transfer is already applied.
```

`report.network` preserves the original circuit, cell, wall and motion accounts
for the post-transfer stage. The outer report independently checks both gas
masses and the full energy change in electrical storage, thermal cells,
reservoir, chamber and piston motion against external electrical work,
prescribed heat, outside-pressure work and force work. Nozzle enthalpy, Joule
heat and reciprocal wall heat cancel as internal transfers. A nonzero mass
flow retains the chamber's explicit missing-stream-port audit boundary; the
outer paired report provides its complete finite-system balance.

The reservoir is a uniform fixed-volume ideal gas, and the original nozzle
source, donor-fraction and pressure-crossing bounds remain. This is a
**first-order transfer-then-network split over one duration**, requiring time
refinement. The thermal geometry and piston inertia remain fixed. No reaction
chemistry, jet thrust, moving thermal mesh or simultaneous nozzle/network
solve is inferred. Failure or late cancellation returns no tentative physical
endpoint. A shut nozzle reproduces the original closed-network states and
port history exactly.

The regression sources exercise both flow directions, one withdrawal across
many nonlinear trials, conductivity-only electrical data, material and incoming
mixture changes through actual heat/motion, an independent all-store trajectory
balance, shut-path parity and complete rollback. Native Rust tests could not
run in the authoring environment; their presence is not a native pass.

