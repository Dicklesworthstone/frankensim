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

## Directional heat conduction

For an orthotropic material, declare `DirectionalConductionLink` records from
`circuit_enthalpy::material` and use
`assembly.prepare_with_directional_conduction(&cx, &thermal, &directions)`.
Each original link side supplies three principal conductivity property keys,
their explicit selection policy, a nonzero material-frame identity, and a unit
temperature-gradient direction expressed in those material axes. All links
touching one cell must agree on its principal roles, selected claims and frame;
their directions can differ. Link endpoint order, area and half-cell distances
remain with the original thermal geometry.

The projected coefficient is `k_n(T) = sum_i n_i^2 k_i(T)`. Both sides enter the
existing series-resistance law, followed by any declared thermal contact.
Principal claims are selected once, retained through source continuation, and
queried at every candidate cell state alongside its electrical properties.
The normal material sample carries all three principal receipts and the applied
link coefficients. An orthotropic-only cell can set its legacy scalar
`MaterialCell::conductivity` binding to `None`; no substitute scalar is invented.
Changing orientation or principal conductivity data therefore changes heat,
temperature-dependent resistance and the original coupled electrical history.

This also composes with electrical interfaces and enclosure radiation through
the same prepared heat network. For additional links supplied by another
physical owner, use
`prepare_with_directional_conduction_and_link_conductances` and continue to
provide their live values through `resolve_with_link_conductances`. Such links
must be disjoint from the oriented Fourier bindings. A source-integrated
isotropic layer cannot replace an explicitly directional link; it can still
serve other links in the same assembly.

This is a one-dimensional link/Galerkin reduction with a declared gradient
direction. It does not solve arbitrary three-dimensional anisotropic flow or
the transverse gradients of a generally oriented insulated rod. Those require
the spatial tensor-conduction owner. Source intervals, nonnegative principal
conductivities, unit directions and unchanged material endpoints remain explicit
admission conditions.

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
