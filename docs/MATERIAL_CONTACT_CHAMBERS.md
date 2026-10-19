# Electrical contact heat, finite material bodies and multiple gas chambers

`fs_couple::cavity_phs::piston::thermal::boundary::ports::network::step_with_contacts_candidate`
connects the existing electrical-contact/enthalpy solve to the existing
multi-chamber gas and independent-wall solve. It accepts the same chamber and
gas-surface descriptions as `step_candidate`, plus explicit `ContactHeatCells`
and a `ContactThermalCoefficients` callback.

A circuit contact contributes its source-derived series resistance to the
Kirchhoff solve. Its actual Joule heat is split between its two thermal bodies
by the existing contact-deposition implementation. Each body's endpoint
enthalpy determines its material state, electrical properties and conductivity.
Each chamber then uses those same thermal endpoints and its own source-resolved
gas midpoint for all its heat contacts and pressure/momentum work.

Electrical contact heat routes and gas-wall surface routes are different
physical interfaces. Their IDs and areas are not interchangeable. A cell can
receive electrical heat while exchanging signed thermal energy with more than
one gas. Opposite gas-surface transfers remain individually reported even when
their sum is zero. No contact heat is entered again as external supply.

The ordinary and contact-aware gas APIs use one shared network residual. The
original circuit/enthalpy owner alone deposits both bulk and interface losses.
The original electrical-contact source check is reused, not copied into the gas
solver. All trials start from the same accepted circuit, thermal timeline and
gases. Local thermal checks, circuit work/storage checks and the final combined
energy check must pass before any endpoint is returned. Nonlinear trial work,
including the final evaluation, remains under the existing budgets.

The returned `NetworkCandidate` retains `ContactThermalSample` as its material
payload, including canonical electrical heat destinations. `into_parts()`
returns the thermal/circuit endpoints and all gas candidates together. Do not
advance a returned chamber again or apply any reported heat a second time.

## Composition

The callback can use `ElectricalContactNetwork::resolve` for material inputs.
It may also apply the existing thermal-contact and radiation adapters to the
same trial coefficients, retaining those source results in its payload. The gas
contact callback receives that payload and the actual gas midpoint; use the
existing `MaterialConductionLayers::resolve_from_material` to consume a complete
thermal/electrical snapshot, not a second conductivity-only query.

The contact-free `step_candidate` retains its API and source shape. Empty
contact bindings use the original contact-free circuit entry point. No gas
chambers reduces to the original electrical-contact/enthalpy evolution; isolated
gases retain their own original pressure/momentum evolution.

## Data-only material wiring

`network::material::step_with_electrical_contacts_candidate` accepts the existing
`ElectricalContactNetwork`, the thermal/circuit histories, all `Chamber`s and
all gas `MaterialContact` patches. No handwritten property callback is needed:

```rust,ignore
let proposal = network::material::step_with_electrical_contacts_candidate(
    &electrical_contact_network, &thermal, target_s, &circuit, &source_ports,
    &chambers, &gas_patches, controls, &energy_guess_j, || cancelled(),
)?;
```

Electrical heat destinations come from the original contact binding. All gas
patch areas are summed across chambers against each cell's remaining exposed
surface, using the same admission as the ordinary material path. Electrical
conductivity/resistivity choices remain explicit per branch. The ordinary
multi-chamber entry point now also uses prepared transport sources; its new
`step_candidate_with_bases` variant exposes the same explicit basis choices.
Legacy `step_candidate` remains resistivity-only.

`ElectricalContactNetwork::prepare` returns a borrowed
`PreparedElectricalContactNetwork`. It composes the existing bulk preparation
with pinned interface claims, never cached temperatures or coefficients. The
ordinary contact-only thermal step and the new chamber adapter both use it.
Fresh trial values still come from the same bulk/interface evaluators, including
phase coordinates and complete source conditions. A competing claim appearing
later does not make the retained, valid pin ambiguous; leaving the pinned
source's supported domain refuses rather than substituting the competing law.
Preparing at a new starting state explicitly permits new policy selections.
The actual contact source must still agree with the accepted circuit; contact
reconfiguration is not hidden in this preparation.

For custom radiation or thermal-contact composition, prepare once, call
`prepared.resolve(trial_thermal, cancelled)` in the existing callback, and retain
all source payloads. `heat_destinations` exposes the original contact/cell routes
without inferring their partition. Neither preparation nor evaluation advances
physical history.

## Physical scope and checks

The model has fixed circuit/thermal geometry and capacitance/inductance,
explicit contact area and heat partition, uniform thermal cells, closed
frozen-composition gas chambers, and independent translating walls. It does not
infer contact mechanics, arcs, spatial spreading resistance, shared crank
motion, combustion or gas mass exchange. The thermal endpoint/gas midpoint
composition remains generally first order. Energy balance is not a timestep
error estimate or experimental qualification.

Focused native command:

```sh
cargo test -p fs-couple --lib cavity_phs::piston::thermal::boundary::ports::network::
cargo test -p fs-couple --test circuit_electrical_contacts
```

New regressions use real source cards and the existing physical solvers for
contact-only material changes reaching both bodies and both gas motions,
isolation and contact-free reductions, and cancellation/refusal without partial
history. They are synthetic numerical fixtures, not measured-material data.

Further regressions compare the data-only path to the same prepared callbacks,
check aggregate area and cancellation, and distinguish overlapping competing
claims from an unsupported pinned source. The contact-only material-change
trajectory now runs through the data-only multi-chamber entry point.
