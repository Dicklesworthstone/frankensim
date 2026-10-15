# Contact heating, material interfaces and finite thermal networks

`HeatedImpactSystem::step_with_thermal_network` connects the existing collision
and caloric-reservoir runtime to `fs_thermal::enthalpy::DepositionSimulation`.
All receiving cells and links advance on one simultaneous thermal network.
Several contacts may heat the same cell; a cell shared by several links is not
advanced separately for each pair. The original material observer resolves
additional sources and conductances on each actual thermal substep.

## One accepted state

`ReservoirCell { reservoir, cell }` binds every caloric reservoir to exactly one
active network cell, and every active cell must be represented. The adapter
checks the original borrowed caloric object, invariant mass, reference volume,
whole surface area and accepted enthalpy. A merely similar temperature or a
cloned replacement chart is not sufficient. The network clock must equal the
mechanical sample clock. Initial network events must already be processed.

Contact dissipation comes from the existing coupled mechanical solve. The
adapter distributes this known work as a held power over the mechanical sample;
the original thermal solver integrates it over its actual CFL substeps. It is
not replayed as a persistent boundary source or added again by the observer.
Conduction, boundary cooling, phase fractions and thermal time integration remain
with the existing network and original material charts.

Source refresh, mechanical state and history, all caloric reservoirs, the network
and the thermal observer are accepted together. A late material, source,
numerical, budget, topology or cancellation failure leaves them all unchanged.
Observer Clone must isolate mutable physical state; external effects or panics
cannot be rolled back. Pending events, including endpoint events, are refused:
births and scheduled boundary/topology changes need an explicit coupled boundary.

The network audit includes contact work as thermal input. The combined account
subtracts that internal transfer and preserves the independently accumulated
outside heat. Per-cell net exchange is reconstructed from the original solver's
storage change minus its integrated contact source; it is not an independently
resolved flux measurement. Numerical energy error stays visible in the combined
residual rather than becoming an invented heating/cooling term.

## Mechanical contact controls conduction and cooling coverage

`step_with_contact_network` uses the same transaction and adds
`ContactNetworkLink { contact, link_id, cells, selection }` declarations.
The ordered thermal partners must carry caloric curves belonging to the original
ordered mechanical bulk cards. Each contact's area-specific thermal resistance
is queried from its ORIGINAL interface card at the retained interface-law point,
using the existing `fs_conduction::interface::InterfaceResistance` evaluator.
There is no duplicate interface table, material-name inference or fallback.

The original circle/ellipse contact area supplies R''/area through the existing
series-resistance operator. At exactly zero area the network link is absent.
That separation exposes the covered areas to the network's original convection
law. Recontact covers the areas and restores conduction without resetting heat.
All selected links are replaced together through the network's existing contact
update operation; unrelated support links survive. A link ID cannot be reused
for different ordered partners, asymmetric coverage cannot be silently converted,
and the thermal observer cannot overwrite these source-derived conductances.

The interface source is required even while separated. Pressure-equivalent
finite-gap footprints are not accepted as physical thermal transfer geometry.
The existing thermal-network budgets also enforce cumulative surface coverage
and active-link limits across multiple contacts.

```rust,ignore
use fs_couple::render::plate::impact::material::heating::network::ReservoirCell;
use fs_couple::render::plate::impact::material::heating::network::contact::ContactNetworkLink;
use fs_matdb::SelectionPolicy;
use fs_thermal::enthalpy::{CellId, StepConfig};

let bindings = vec![
    ReservoirCell { reservoir: "left".into(), cell: CellId(9) },
    ReservoirCell { reservoir: "right".into(), cell: CellId(3) },
    ReservoirCell { reservoir: "support".into(), cell: CellId(77) },
];
let links = vec![ContactNetworkLink {
    contact: "touch".into(), link_id: 101,
    cells: [CellId(9), CellId(3)], selection: SelectionPolicy::SingleClaimOnly,
}];
let frame = heated.step_with_contact_network(
    &mut thermal, &mut material_observer, &bindings, &links,
    StepConfig { max_dt_s: thermal_dt, energy_tolerance_j: thermal_allowance },
    &forces, &gate,
    |contact, stores| resolve_original_contact(contact, stores),
)?;
```

The source callback is the existing unloaded material-feedback seam described in
`MATERIAL_CONTACT_REFRESH.md`; it may return `Ok(None)` to keep fixed contact
coefficients. Changing source material properties affects the same original
mechanical contact, interface queries, heat stores and network, not a separate
instrument-specific simulation.

## Physical limits and verification

Reference geometry and uniform thermal cells remain assumptions. The contact
area and selected interface coefficients are held at the proposed mechanical
endpoint for the thermal interval. This is first-order sample-held coupling,
not continuous event-resolved contact conduction. Thermal subcycling does not
remove that approximation; refine the mechanical clock as well. The existing
conduction CFL does not bound arbitrary nonlinear source stiffness.

This does not supply spreading resistance, local flash temperature, deformation,
flooding, phase-changing contact mechanics, or loaded thermoelastic work. Material
source refresh remains unloaded-only. Heat partition, thermal lumping, physical
partner assignment, and supported source domains remain explicit inputs.

The existing focused target now includes fixed-network parity, analytic
three-node refinement, shared contact heat and other sources, complete observer
rollback, source-selected multi-contact conduction, ambient-area coverage and
invalid binding/interface/topology cases:

```
cargo test -p fs-couple --test material_contact_heat -- --test-threads=1
```

Synthetic coefficients in these tests do not claim measured-material accuracy.
