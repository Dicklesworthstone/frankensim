# Many thermal contacts and gas chambers on one material network

`cavity_phs::piston::thermal::boundary::ports` gives one gas multiple thermal
contacts. `ports::network` connects several independent gas/piston states to
one existing material/circuit/enthalpy network. One chamber may contact several
cells; several chambers may contact the same cell. Each retains its own gas
composition, mass, volume, momentum, material sources and heat records.

## One gas state, distinct signed contacts

Each chamber performs one scalar total-heat solve per network trial. All its
contacts see the same gas energy midpoint from the existing gas/momentum owner.
The individual heat is `dt * G_contact * (T_cell_endpoint - T_gas_midpoint)`.
The sum supplies the gas; each receiving cell is debited by its own contact's
heat. A hot contact and a cold contact can carry opposite nonzero heat while
the net gas heat is zero. Neither averaging temperatures nor dividing net heat
by conductance can recover that exchange. Every original contact source is
retained, including at equilibrium or explicit zero conductance.

The difference between evaluated contact heat and total heat applied to the gas
is explicitly checked against source and energy tolerances. Contact values are
not rescaled to hide it. The old single-contact API shares the same source-solve
implementation and keeps its original types and arithmetic. Empty contacts
reproduce isolated gas motion; empty chambers reproduce the circuit/thermal
network without a new physical model.

## Data-only network coupling

```rust,ignore
use fs_couple::cavity_phs::piston::thermal::boundary::ports::network::{
    Chamber, Contact, NetworkSolve,
    material::{MaterialContact, step_candidate},
};
use fs_couple::cavity_phs::piston::thermal::boundary::network::material::GasConductionGeometry;

// network, thermal, circuit, gas_a/gas_b and original source/solver controls
// are caller-owned existing objects. All inputs use this exact same interval.
let chambers = [
    Chamber { id: 10, piston: &gas_a, input: input_a,
        piston_solve: motion_a, exchange: heat_a },
    Chamber { id: 20, piston: &gas_b, input: input_b,
        piston_solve: motion_b, exchange: heat_b },
];
let geometry = GasConductionGeometry {
    area_m2: 0.008, gas_distance_m: 0.001, cell_distance_m: 0.02,
};
let contacts = [
    MaterialContact { contact: Contact { id: 1, chamber: 10, cell: hot_cell }, geometry },
    MaterialContact { contact: Contact { id: 2, chamber: 10, cell: cold_cell }, geometry },
    MaterialContact { contact: Contact { id: 3, chamber: 20, cell: hot_cell }, geometry },
];
let controls = NetworkSolve {
    network: network_controls, maximum_chambers: 2, maximum_contacts: 3,
    maximum_boundary_solves: 600, total_energy_tolerance_j: 1e-6,
};
let candidate = step_candidate(&network, &thermal, target_s, &circuit, &ports,
    &chambers, &contacts, controls, &energy_guess_j, || cancellation_requested())?;
let (next_thermal, next_circuit, report, balance) = candidate.into_parts();
for chamber in report.source.chambers {
    let (next_gas, contact_records, gas_balance) = chamber.boundary.into_parts();
    // Install under chamber.id with ALL other returned endpoints.
}
// Electrical/contact heat is already in these endpoints. Do not add it again.
```

`MaterialNetwork` remains the original caloric, electrical and thermal binding.
The adapter reuses the existing `MaterialConductionLayers` source and geometry
checks, comparing each body-side source with the network's actual endpoint
conductivity bundle. There is no replacement material registry, default gas,
resistivity law, circuit solver or temperature update. Explicit phase coordinates
can change conductivity during a latent plateau through the existing chart.

Patch areas are summed per cell across *all* chambers. The sum must fit the
exposed area after authored thermal links/reservoirs, and ambient convection on
a gas-contact cell must be zero to avoid duplicate surface exchange. This is an
area budget, not a geometric proof that two positioned patches do not overlap.
Missing conductivity refuses even for a zero-area patch. Original claim choices
and other physical conditions are retained. Custom source laws can use the
generic `network::step_candidate`; they own their physical coverage assumptions.

## Conservation, history and limits

All network proposals start from the same electrical and thermal histories.
Every chamber proposal starts from its own same accepted gas and wall velocity.
Rejected trials cannot drain a capacitor multiple times or repeatedly advance a
gas. Stable IDs control computation and report ordering, not input-vector order.
Duplicate contact/chamber IDs, foreign cells, mismatched clocks and repeated
references to the same gas state refuse before physical stepping.

Each original solver retains its own physical gates. The final combined balance
contains electrical storage, all cell enthalpies, every gas internal energy and
wall kinetic energy, exterior-pressure work, applied-force work, external
circuit supply and independent prescribed heat. Internal electrical heat and
all gas-contact transfers are excluded from external supply. Global chamber
work includes discarded numerical trials and final admission. Source failures,
insufficient work, unsupported numerical probes, cancellation and late balance
failures return no joint endpoint.

This model has fixed thermal geometry and circuit topology/C/L, uniform closed
frozen-composition gases, and independently translating walls. It does not model
mass exchange, valves, combustion, shared crankshafts, acoustic propagation
between chambers, spatial convection, melt motion or live acoustic-renderer
reconfiguration. The common updated gas states expose their own pressure,
transport and acoustic properties. Temperature-endpoint network stepping plus
gas midpoint stepping is generally first order. Energy conservation is neither
an error bound nor a uniqueness, entropy-production or unconditional-stability
claim; source support and temporal refinement remain necessary.

Focused native tests, included by the existing piston test selection:

```sh
cargo test -p fs-couple --lib cavity_phs::piston::thermal::boundary:: -- --test-threads=1
```

Nine new tests cover single/empty parity, opposite signed heat, material-source
changes through both chambers and finite electrical history, original source
records, canonical ordering, independent gas heating, aggregate surface limits,
source failures, final-work exhaustion and cancellation without partial state.
Local cargo/rustc were unavailable during authoring. Independent simultaneous
Python equations and refinement checks are not native Rust test passes or
experimental material validation.
