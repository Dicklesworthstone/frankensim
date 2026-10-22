# Solved support pressure, material heat transfer, and finite electrical supply

`ThermoplasticSolid::advance_electrical_with_supports_candidate` adds actual
whole-face support pressure to the existing electrical/3D-solid/enthalpy solve.
The input `ThermalSupportContact` is the original contact owner's source binding,
not a second pressure or contact-geometry model.

```rust,ignore
let next = solid.advance_electrical_with_supports_candidate(
    &thermal, &circuit, &prepared_heat_network, &support_contacts,
    &ports, target_s, &prescribed_m, &nodal_forces_n,
    mechanical_settings, electrical_solid_controls, &energy_guesses_j, cx,
)?;
let actual_contacts = &next.report().support_contacts;
let (thermal, solid, circuit, report) = next.into_parts();
```

At every thermal trial, native tetrahedral equilibrium supplies compression
`p = -n . sigma . n`, using the original reference-face outward normal. The
original interface evaluator queries the explicitly pinned resistance at that
pressure and both endpoint temperatures. The existing series-resistance owner
combines it with the material network's freshly resolved bulk conductance:

```
1 / G_total = 1 / G_bulk + R_interface_per_area / A
```

The same existing circuit/enthalpy residual uses this conductance along with
endpoint resistivity, radiation, physical plastic heat and admitted reversible
heat. Thus a change in load or interface material data reaches thermal exchange,
current, finite charge, yield response and compatible displacement together.
Every trial restarts from the accepted histories. No interface heat is added
again as a source: the thermal link already transfers it between its two cells.
Original electrical, material and energy reports remain intact. The new
`support_contacts` account gives pressure, original interface receipt, bulk
coefficient, combined conductance and actual signed heat.

Support links must NOT also appear in `MaterialHeatNetwork::thermal_contacts`;
that would count their interface resistance twice. Electrical contacts on these
same links are also refused: their fixed source conditions cannot silently stand
in for a coupled pressure-dependent electrical law. Unrelated thermal/electrical
interfaces and authored radiation remain available. The prepared network must
still include the support link's bulk material conduction. Missing data refuses;
explicit zero bulk conductivity insulates but still validates the contact source.

Unpowered consumers use `advance_with_series_thermal_contacts` with the same
support bindings and a callback supplying fresh bulk-only coefficients. Both
paths share the original face admission, pressure evaluation and series law.
The original interface-only API remains available and still rejects a second
conductance owner. Empty support lists preserve the original powered arithmetic.

## Model boundaries and verification

This is quasistatic small-strain/small-rotation J2 mechanics with fully solid
caloric charts and the existing energetic-hardening restrictions. Contacts are
complete prescribed exterior faces to finite thermal supports, not two deformable
meshes. Footprint, bulk path geometry, electrical topology and radiative view
factors stay fixed within an interval. Tension refuses: opening, friction,
evolving area, electrical contact pressure feedback, field forces, melting and
inertial dynamics are not inferred. Endpoint thermal/mechanical response with
midpoint circuit current is generally first order; split at events and refine.

Three unpowered and five powered regressions reuse the loaded-solid fixtures.
They cover the analytical two-cell BE limit, material-source retention, bulk
isolation, compression-to-current/displacement feedback, both charge signs,
original-component replay, unloading, source expiration, duplicate corrections,
work exhaustion, cancellation and empty-support equivalence. The existing
whole-cell electrical fixture's path area is corrected to `tet_volume / length`;
its original independent resistance formula and physical assertions are retained.

```sh
rch exec -- cargo test -p fs-couple --test solid_plastic_material_transport
```

Native Rust compilation/tests remain unverified in the editing environment:
DSR, its documented fallback, RCH and Cargo could not start. Independent
Python/SciPy equations passed 108 synthetic pressure/interface/resistivity/charge/timestep
cases and four linear series cases. These are equation checks, not execution of
Rust, experimental calibration, or proof of rollback or temporal convergence.
