# State-dependent material contacts in finite thermal networks

`fs_conduction::transient::deposition::contact::apply_material_contacts` adds
explicit contact resistance to the original bulk-only thermal link coefficients.
Use it inside the existing implicit thermal, circuit/enthalpy, or generic
network/gas coefficient callback. It introduces no new time integrator.

The physical path is bulk conduction A -> interface resistance -> bulk
conduction B. The original `InterfaceResistance` resolves the contact source;
`SeriesThermalResistance` composes `1/G_bulk + R''/area`. The existing thermal
residual exchanges the resulting signed heat once between the two cells.
Electrical losses, radiative exchange and independent sources remain with their
existing owners. Changing contact data can therefore change temperatures,
phase fractions, electrical resistance, and finite capacitor/inductor history.

## Source and state binding

`fs_material::thermal::interface::ThermalInterfacePoint` uses the existing
`ThermalMaterialPoint` to update two independent temperatures and optionally two
liquid mass fractions. Both actual bulk-card content identities must match the
phase states. The interface card's ordered surface material IDs must match the
same bulk cards. Texture, medium, environment, history and all other supplied
query coordinates remain part of the original interface source.

An optional mean-temperature axis explicitly means `(T_A+T_B)/2`. It is not an
inferred equilibrium or solved interfacial temperature. All evolving axes must
be distinct; pressure/gap/lot and other undeclared coordinates remain unchanged.
Only static dimensional scalar contact resistance and one-dimensional curves
on declared thermal coordinates are admitted. Missing, ambiguous, out-of-domain
or unsupported source laws refuse; exact claim pins remain exact pins.

Each `MaterialContactLink` supplies the existing link ID, two `ContactCell`
bindings, the original interface condition point and claim selection. The
caloric curves must be the same borrowed instances used by the thermal timeline.
Link A/B order and area remain authoritative. No material names choose a contact
coefficient and a material swap must supply the corresponding interface data.

```rust,ignore
// In the original circuit_enthalpy::step_candidate coefficient callback:
let mut coefficients = material_network.resolve(trial_thermal, cancelled)?;
let contacts = apply_material_contacts(
    trial_thermal, &contact_bindings,
    &mut coefficients.conductances_w_k, cancelled,
)?;
// Return the updated coefficients with BOTH original source payloads retained:
let source = (coefficients.source, contacts);
```

Every evaluation starts with fresh **bulk-only** coefficients. Do not feed a
previously contact-corrected map back into this function: that counts a contact
twice. Every selected link must have a supplied bulk coefficient; unselected
links are unchanged. The returned sample retains original contact receipts,
endpoint phases, contact area, bulk conductance and the series-resistance report.
Bulk provenance stays in the original material-network payload; the series
report does not invent a replacement bulk-property receipt.

Overrides are staged before changing the coefficient map. Late source refusal,
invalid geometry or cancellation preserves the entire input map. A caller's
accepted physical state is never advanced by this coefficient evaluation.
Zero contact area or zero bulk conductivity isolates a link exactly, while
still checking the original contact source and both thermal identities.

## Numerical and physical scope

These are fixed-geometry uniform cells with negligible interfacial heat storage.
The contact law sees the declared representative cell states, not separately
solved face temperatures after internal layer drops. Its data must support that
lumped closure. This is not a contact-area, spreading-resistance, adhesion,
pressure-history, damage, interface phase-mixture or moving-melt solver.

The calling thermal scheme and its source/residual checks remain authoritative.
For the existing backward-Euler circuit/enthalpy composition the overall method
is generally first order. Nonlinear convergence and energy balance are not
experimental validation or a timestep-error bound.

Focused native tests:

```sh
cargo test -p fs-conduction --test contact_material_state --test deposition_contact
cargo test -p fs-couple --test circuit_material_contact
```

Tests cover two-sided phase/source binding, implicit latent-contact feedback,
source-only changes through both body temperatures and capacitor discharge,
series composition, isolation, and cancellation/refusal without partial updates.
Local Rust execution was unavailable during authoring; the independent numerical
reference checks are not native test passes or measured-material qualification.
