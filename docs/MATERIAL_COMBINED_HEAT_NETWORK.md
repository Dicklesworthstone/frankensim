# One material assembly for electrical and thermal exchange

`fs_couple::circuit_enthalpy::contacts::heat::MaterialHeatNetwork` composes an
existing `ElectricalContactNetwork`, optional `MaterialContactLink` records,
and an optional `MaterialRadiation` enclosure. It introduces no material law,
Kirchhoff solver, radiosity solver or time integrator.

```rust,ignore
let assembly = MaterialHeatNetwork {
    electrical: &electrical_contacts,
    thermal_contacts: &thermal_contacts,
    radiation: Some(&enclosure_binding),
};
let selected = assembly.prepare(&cx, &thermal)?;
let proposal = selected.step_candidate(
    &cx, &thermal, target_s, &circuit, &source_ports, controls, &energy_guess_j,
)?;
let (next_thermal, next_circuit, report) = proposal.into_parts();
```

Every evaluation uses the same trial enthalpies. Fresh bulk conductances receive
thermal contact resistance in series exactly once. Signed radiative powers are
added exactly once. The existing circuit owner deposits bulk and electrical
contact losses into their original cells. Authored independent heat sources,
convection, reservoirs, event handling and latent storage stay active in the
original thermal owner. Neither a rejected trial nor a source query advances
accepted charge, flux, cell enthalpy or the event cursor.

Preparation retains the original bulk, electrical-interface, thermal-interface
and emissivity claim selections. It does not cache physical properties. Each
subsequent evaluation queries those claims at the actual temperature and phase.
An unsupported selected source refuses rather than switching to another model.
Preparing again is an explicit new source-selection operation.

Radiating cells must use their existing bulk cards, caloric laws and complete
physical conditions. A radiation binding cannot quietly give the same body a
different pressure, temperature-axis meaning or liquid-fraction convention.
The original radiation and contact owners retain geometry, source and surface
coverage admission. All original source receipts and per-mechanism reports are
available in `report.source.material`; the complete bulk material snapshots
remain in its `material` member for downstream consumers.

## Energy and scope

The standard circuit/enthalpy report preserves the actual source arithmetic,
including tiny net-enclosure roundoff. Radiation in this closed enclosure is
internal exchange, not external supply. The additional combined gate checks
`report.total_energy_residual_j + dt * sample.radiation_net_power_w()` against
the existing total-energy tolerance. True non-electrical external heat is
`report.non_electrical_heat_j - dt * sample.radiation_net_power_w()`. No surface
heat is clipped or redistributed to manufacture an exactly zero total.

This is fixed-geometry, fixed-topology Ohmic transport with uniform thermal
cells, supplied contact areas/heat splits, and a closed opaque gray-diffuse
radiation enclosure. Thermal and electrical contact resistances are independent
properties, not inferred from one another. View factors and capacitance/
inductance remain fixed. Radiation may coexist with conduction or convection
on an exposed surface because these are distinct mechanisms. No participating
gas radiation, contact mechanics, spreading-temperature law, arc or flowing
melt is implied. Overall time accuracy remains that of the existing coupled
endpoint/midpoint method, generally first order. Energy balance is not an
experimental validation or a time-error certificate.
