# Finite contact heat exchange on the collision clock

`HeatedImpactSystem::step_with_coupled_transport` connects the contact-heating
runtime to interacting finite thermal stores. The callback receives all old
thermal states together, each store's actual contact dissipation and a read-only
unpublished mechanical endpoint. Both partners can therefore exchange heat
without either callback mutating or advancing the other partner prematurely.

The preceding `step_with_transport` API remains available for independent
thermal environments and delegates to the same acceptance path. See
`MATERIAL_CONTACT_HEAT.md` for attachment, heat routing and existing limitations.

## One source and geometry path

`ImpactThermalTrial::contact_observation(key)` evaluates the original material
normal-contact owner at the proposed mechanical endpoint. Its pressure and
footprint can drive the selected thermal interface. `contact_patch_area_m2(key)`
uses the circular radius or both true elliptic axes. Separation has exactly zero
area. A finite-gap pressure-equivalent footprint is not silently promoted to a
physical thermal contact area: that case refuses.

`fs_conduction::interface::lumped_pair::advance_contact_pair` consumes two
`UniformEnthalpyStepInput` values and the original `InterfaceResistance`. It
forms R''/area through the existing source-bearing interface owner and calls
`solve_lumped_contact_step`. Both contact-heating sources enter that same
implicit finite-body equation. A heat-only intermediate state need not be
admissible when simultaneous outward conduction keeps the actual endpoint valid.

At zero area the thermal link is exactly open, not an arbitrarily large finite
resistance. Each body still receives its irreversible contact work, if any.
Neither temperature nor latent fraction resets on separation or recontact.
The original interface receipt, resistance network and finite-body solve report
remain in the returned proposals. A zero-area report retains the interface but
has no fabricated resistance or contact-solve report.

## Composition

The two stores below were attached to `heated` as `left` and `right`; their
material charts and thermal masses remain the runtime's original inputs.
`interface` is an already resolved area-specific resistance; `left_transport`
and `right_transport` are the original admitted bulk transport objects.
The caller must bind the selected interface to the correct physical partners
and choose its source coordinates; material-name guessing is not performed.

```rust,ignore
use fs_conduction::interface::lumped_pair::advance_contact_pair;
use fs_couple::render::plate::impact::ImpactError;
use fs_couple::render::plate::impact::material::heating::transport::CoupledThermalStep;

let frame = heated.step_with_coupled_transport(&forces, &gate, |trial| {
    let stores = trial.reservoirs();
    if stores.len() != 2 || stores[0].key != "left" || stores[1].key != "right" {
        return Err(ImpactError::Invalid("expected the two bound thermal partners"));
    }
    let area = trial.contact_patch_area_m2("touch")?;
    let transfers = advance_contact_pair(
        &cx,
        [stores[0].thermal, stores[1].thermal],
        [&left_transport, &right_transport],
        &interface,
        area,
        thermal_energy_allowance_j,
    ).map_err(|error| ImpactError::Owner(error.to_string()))?;
    Ok(CoupledThermalStep {
        transfers: Vec::from(transfers),
        boundary_heat_j: 0.0, // This entire two-body thermal collection is isolated.
    })
})?;
```

A fixed physical link can instead supply its declared nonzero area even during
mechanical separation. The general joint callback also accepts other thermal
network owners; do not apply pair solves sequentially to a shared node and call
that a simultaneous network solve.

## Energy and acceptance

The finite solver removes Q from A and supplies the same Q to B. Those signs
appear in the canonical per-store `external_heat_j`; their sum is zero for an
internal link. Contact dissipation is a separate, already counted internal heat
source. It is never redeposited from endpoint force or the instantaneous receipt.

`CoupledThermalStep::boundary_heat_j` independently declares heat from outside
the entire collection. Returning a recipient credit without its donor debit
refuses for an isolated collection even when the recipient's local balance is
satisfied. Wrong proposal counts, foreign charts, material-domain excursions,
local/cumulative/combined energy failures and cancellation publish neither
mechanical state/history/time nor any thermal store. Callbacks must return
proposals, not mutate external live state that the adapter cannot roll back.

## Explicit physical limits

Mechanics still uses frozen constitutive parameters: this is not backward
thermal feedback into stiffness, shape, yield or melting mechanics. The thermal
step uses endpoint contact geometry and holds the selected R'' during the step;
it is not continuous collision/event detection or a resolved integral of contact
area over a tick. Resolve time-dependent source coordinates explicitly and use
timestep refinement to assess that approximation.

Each thermal store is uniform with fixed reference geometry and negligible pV
work. The original reduced Biot gate spreads total contact conductance over its
whole surface; it does not establish local spreading or flash-temperature
validity under a small patch. No ambient, radiation or gap-gas conduction is
invented. The example/source fixtures are synthetic, not measured materials.

Focused verification uses the existing target:

```
cargo test -p fs-couple --test material_contact_heat -- --test-threads=1
```

The added regressions compare the original finite solver with an independent
constant-Cp identity; compare unchanged mechanical trajectories and endpoint
contact observations; vary only thermal-interface data and latent heat; exercise
opening, recontact and a joint-valid/source-only-invalid endpoint; and verify
cancellation, missing donor debit, source, geometry and clock refusal behavior.
