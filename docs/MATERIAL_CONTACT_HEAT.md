# Contact work, material enthalpy and thermal exchange

`fs_couple::render::plate::impact::material::heating::HeatedImpactSystem`
connects the existing source-bound impact runtime to finite thermal reservoirs.
It uses the original contact force law, time integrator and material phase chart.

## Physical data path

A `MaterialContact` retains both bulk states, its ordered interface and the
selected normal law. The existing `fs-dcontact` owner now exposes irreversible
power per obstacle, sharing its original `loss_force_at`, quadrature and gap
calculation. The heated runtime samples that power at the **same midpoint and
Gonzalez port velocity** used by `fs-phs`. Endpoint force times endpoint velocity
is not a substitute for the accepted step's irreversible work.

Construct a `CaloricReservoir` for each thermal destination from its original
`EquilibriumEnthalpyPhaseCurve`, admitted initial state, invariant mass and
explicit volume/area. This also accepts `ResolvedHeatCapacityCurve::curve()`;
keep that source object alive to retain its complete property-query receipts.
The equilibrium chart handles temperature-dependent heat capacity and latent
heat without a second caloric formula. Geometry and thermal lumping are caller
assumptions; thermal mass is not inferred from the contacting surface.

A `ContactHeatRoute` names an existing contact and two reservoir keys. Its
explicit fraction goes to A; B receives the remainder. No equal split, effusivity
correlation or material-name inference is silently selected. Multiple contacts
can share a reservoir. A contact may be routed only once. Unrouted contact,
modal damping, felt crushing, radiation and fluid losses remain outside these
thermal stores, rather than being incorrectly assigned to contact heating.

```rust,ignore
use fs_couple::render::plate::impact::material::heating::{
    CaloricReservoir, ContactHeatRoute, HeatedImpactSystem,
};

// `mechanics` is an existing source-bound ImpactSystem. The curves and initial
// states below are already admitted by their original fs-material owners.
let stores = vec![
    CaloricReservoir::new("a".into(), &curve_a, state_a, mass_a, volume_a, area_a)?,
    CaloricReservoir::new("b".into(), &curve_b, state_b, mass_b, volume_b, area_b)?,
];
let routes = vec![ContactHeatRoute {
    contact: "interface".into(), to_a: "a".into(), to_b: "b".into(),
    fraction_to_a: declared_heat_partition,
}];
let mut heated = HeatedImpactSystem::new(mechanics.prepare_analytic()?, stores,
    routes, absolute_thermal_energy_allowance_j)?;
let frame = heated.step(&generalized_forces, &gate)?;
```

Preparation/attachment preserves the current mechanical state and clock; old
losses are never replayed as new heat. Heat inventories use compensated summation
and an absolute energy budget, checked locally and cumulatively. A material-domain
failure, thermal imbalance, mechanical refusal or cancellation publishes neither
new mechanics nor any new thermal reservoir. The read-only mechanics accessor
prevents advancing behind the thermal clock. `into_parts` returns both owners.

## Existing conduction, convection and finite-reservoir owners

`step_with_transport` passes each store through the canonical
`UniformEnthalpyStepInput` / `UniformEnthalpyStep<R>` carrier. Its internal heat
is the store's actual contact work; initial state, mass, volume, surface area
and duration all come from the same heated runtime. Every store is visited even
when dry, so thermal relaxation continues after separation.

```rust,ignore
let frame = heated.step_with_transport(&generalized_forces, &gate, |key, input| {
    environments[key].advance(&cx, input)
        .map_err(|error| ImpactError::Owner(error.to_string()))
})?;
```

For `fs_conduction::lumped::LumpedThermalEnvironment`, this invokes the existing
implicit enthalpy solver, conductivity/emissivity models and Biot admission.
No conduction law is copied into the collision runtime. A different owner can
return the same carrier, including a finite thermal contact proposal. Such a
callback must **propose**: it must not mutate an external live reservoir before
the whole transaction succeeds. Its original concrete reports are returned for
joint publication. This adapter cannot undo arbitrary callback side effects.
Returned states must belong to the exact original chart. The solver's declared
energy allowance cannot exceed the caller's thermal allowance. State, boundary
heat and cumulative inventories are independently checked before publication.
An adiabatic continuation starts from the actual accepted thermal endpoint;
it neither erases cooling nor repairs a thermal solver's disclosed residual.

## Energy ownership and limits

The combined balance is

```
mechanical energy change + thermal enthalpy change
+ mechanical dissipation not deposited in the stores
- applied mechanical work - signed boundary heat = numerical residual.
```

`mechanical.dissipated_energy_j` remains the original total mechanical loss.
Do not add `deposited_heat_j` to it again when checking the combined balance.
`remaining_mechanical_loss_j` discloses signed subtraction roundoff rather than
clamping it into invented physical heating/cooling.

This is **one-way thermalization with frozen mechanical constitutive data**.
It does not implement thermal expansion, temperature-dependent contact response,
plasticity, fracture, flash-temperature fields or moving melt. A thermal absorber
may cross a phase plateau without implying that a rigid colliding body has become
liquid. Fixed contact maps, prescribed thermal geometry and a declared partition
remain assumptions. Internal adaptive substepping is not exposed by this wrapper;
use an explicitly chosen fixed timestep and the existing Newton work limits.
The adiabatic path reuses prepared buffers; the report-returning transport path
allocates at most 64 store proposals. Neither is a hard-real-time certification.

Focused tests:

```
cargo test -p fs-couple --test material_contact_heat -- --test-threads=1
cargo test -p fs-couple --test material_impact -- --test-threads=1
```

Fixtures are explicitly synthetic, not measured properties of real metals.
