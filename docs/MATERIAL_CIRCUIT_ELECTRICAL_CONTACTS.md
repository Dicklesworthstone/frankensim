# Electrical contact heat on the shared material circuit and thermal network

An electrical connection is not only a bulk resistor. The interface's specific
resistance belongs to its ordered surfaces, finish, pressure, environment and
history. `material_circuit::contacts::SeriesElectricalContact` consumes the
existing `fs_material::electrical_interface` result and actual contact area:

```
R_branch = R_bulk + sum(R_specific_contact / contact_area)
P_bulk = I^2 R_bulk
P_contact_i = I^2 R_contact_i
```

The original Kirchhoff descriptor still solves node potentials, current, stored
charge and flux. There is no second RC/RLC integrator or fabricated bulk material
for a contact. The complete source and geometry for each voltage drop and heat
record remain attached. `MaterialCircuitStep::heat_j` includes bulk AND contact
losses; `resistors` contains only bulk heat, and `contacts` contains interface
heat. A `ResistorHeat::voltage_v` is the bulk drop, which equals the branch
terminal voltage only when no series contact is present.

## Electrical storage and source changes

Use `MaterialCircuit::new_with_contacts` with numeric, explicitly authored C/L,
or call `rebind_contacts_candidate` on a circuit from `new_with_storage` to
combine source-backed R/C/L with source-backed interfaces. Contact rebinding is
an explicit zero-time change that retains charge, flux and storage sources and
checks algebraic compatibility. Storage rebinding likewise retains contacts
and its existing parameter-work account.

Advance contact-bound circuits with `step_with_contacts_candidate`, supplying
fresh bulk AND interface sources. The ordinary contact-free `step_candidate`
refuses a contact-bound circuit instead of dropping its resistance or heat.
Contact IDs and branch assignments remain fixed within one step. Explicitly
sourced zero resistance is accepted; absent source data is not perfect contact.
Zero area is not approximated by a large finite resistance: circuit opening,
switching impulses and arcs need their own models and are not supported here.

## Finite-body heat deposition

`circuit_enthalpy::contacts::step_candidate` adds the missing two-sided heat
route to the original circuit/enthalpy nonlinear solve. Each `ContactHeatCells`
provides its contact ID, two receiving `CellId`s in source A/B order and an
explicit fraction entering A. The remainder enters B. Bulk heat destinations
remain separate. No default 50/50 or inferred heat-partition law is supplied.

The public data-only `ElectricalContactNetwork` wrapper binds contacts to
existing thermal links. The link supplies A/B order and physical area; the same
cell caloric charts supply temperatures and declared phase fractions. The
shared `ThermalInterfacePoint` retains both material identities, independent
coordinates, original pressure and other external conditions. Resistance is
resolved at the actual proposed thermal endpoint, not at an independently
advanced temperature. Phase-dependent contact resistance can therefore change
current during a constant-temperature latent plateau.

```rust,ignore
use fs_couple::circuit_enthalpy::contacts::material::ElectricalContactNetwork;

let network = ElectricalContactNetwork {
    bulk: &material_network,
    contacts: &contact_bindings,  // MaterialElectricalContact records
    electrical_bases: &bases,    // explicit conductivity/resistivity per branch
};
let initial = network.resolve(&thermal, || false)?;
let circuit = circuit.rebind_contacts_candidate(initial.contacts, &ports)?;
let proposal = network.step_candidate(
    &thermal, target_s, &circuit, &ports, controls, &energy_guess_j, || false,
)?;
let (next_thermal, next_circuit, report) = proposal.into_parts();
// Publish both. All bulk/contact heat is already included in the thermal state.
let interface_heat_j = report.electrical.contacts[0].heat_j;
```

Every thermal trial starts from the SAME accepted circuit and enthalpy history.
Individual contact heat is partitioned and deposited once; the thermal solver
still owns conduction, reservoirs, independent heat and phase storage. The
final accepted proposal checks per-cell residuals, total thermal energy,
electrical terminal work, branch loss and the combined storage/supply balance.
Newton stencils and fresh final admission spend the original circuit budget.
Cancellation or failed source/energy checks publish neither endpoint.

Initially selected bulk and accepted contact card/claim identities must remain
unchanged throughout the solve. Unsupported or ambiguous source queries refuse;
explicit claim pins remain authoritative. The convenience wrapper supports the
existing conductivity/resistivity basis choices and complete source snapshots.
A state that changes the selected claim is refused, not silently accepted as a
different material law. Use the generic contact callback to compose additional
material thermal-contact resistance or radiation; keep all source records in
its payload and do not add electrical heat a second time.

Legacy contact-free circuit, storage, thermal and piston entry points retain
their existing APIs. They do not automatically gain two-sided contact heat
routes: use the contact-aware entry point for contact-bound circuits. A legacy
coupled call refuses absent routes rather than silently discarding contact heat.

## Scope and checks

This is a fixed-topology, fixed-geometry Ohmic network with zero interfacial
energy storage and uniform receiving thermal cells. Contact area and heat split
are physical inputs, not solved contact mechanics or spreading temperature.
Static scalar laws or one-dimensional curves on declared thermal coordinates
are admitted by the material wrapper; dynamic-only source laws refuse. It does
not model thermoelectric effects, arcing, degradation, an open contact or a
spatial current field. Interface thermal conductance and electrical resistance
are different properties; neither is inferred from the other.

Thermal rates/materials use endpoint values and current uses the original
circuit midpoint. The coupled method is generally first order. Residual and
energy tests are not a timestep-error bound or measured-material qualification.

```
cargo test -p fs-couple --test material_circuit_contacts --test circuit_electrical_contacts
```

Eight authored Rust regressions cover analytic circuit discharge and separate
losses, area/zero/reversed-current cases, storage-source preservation, explicit
heat splits, data-only changes through two bodies and finite electrical supply,
latent-phase feedback, independent heat and refusal/cancellation. Native Rust
compilation and tests were unavailable during authoring. Independent Python
checks of the equations are numerical references, not native test passes.
