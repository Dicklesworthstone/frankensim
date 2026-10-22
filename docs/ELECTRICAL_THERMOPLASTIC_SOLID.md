# Finite electrical supply coupled to loaded 3D material response

`ThermoplasticSolid::advance_electrical_candidate` connects the native loaded
P1/J2 tetrahedral solid to the existing contact-capable circuit/enthalpy solve.
One prepared `MaterialHeatNetwork` supplies bulk resistivity/conductivity,
thermal-contact resistance, electrical contacts, and closed-enclosure radiation.
No dummy heater, second thermal timeline, material law or nonlinear solver is
introduced. The old explicit observer path remains unchanged.

```rust,ignore
let prepared = heat_network.prepare(cx, &thermal)?;
let candidate = solid.advance_electrical_candidate(
    &thermal, &circuit, &prepared, &ports, target_s,
    &prescribed_displacements_m, &nodal_forces_n, mechanical_settings,
    electrical_solid_controls, &cell_energy_guesses_j, cx,
)?;
let (thermal_next, solid_next, circuit_next, report) = candidate.into_parts();
```

`ElectricalSolidSolve` is in `solid_enthalpy::plastic::electrical`. It contains
existing `CircuitEnthalpySolve` controls, a separate all-trial mechanical work
budget, and the native yield-history tolerance. Guesses use joules in sorted
thermal-cell-ID order. Structural material and energy arrays retain tetrahedron
connectivity order; resistor heat destinations come from the prepared network.

## One accepted endpoint

Each residual or Jacobian trial starts from the identical accepted charge/flux,
plastic history, displacement and enthalpy. The original native mechanical
update contributes its physical plastic heat and signed reversible caloric
exchange. The original circuit owner deposits bulk and two-sided interface
Joule heat once. The final accepted callback carries its actual native solid
solution; there is no post-solve mechanical reconstruction at another state.

Temperature-dependent resistivity changes current, available electrical energy,
local temperatures, yield response, and the loaded solid's free displacement.
Material cards, original caloric objects and fixed conditions must agree before
solving. This uses the SAME admission function as the unpowered heat/solid
entry point. The material owner's rate-independent source check additionally
rejects transport that would freeze a mechanical or frequency coordinate while
that state changes. Dependency checks run once on pinned definitions; values
and source validity are still evaluated at every proposed thermal endpoint.

Circuit/thermal reports and complete endpoint material/radiation/contact reports
are retained in `report.circuit_thermal`. That subsystem report includes
mechanical heat in its callback-supply account. The separately reported CLOSED
balance is:

```
electrical storage change + elastic/hardening storage change
+ thermoelastic internal-energy correction + caloric storage change
+ numerical mechanical dissipation
- electrical supply work - mechanical boundary work - true external heat
= 0 (within the original combined joule tolerance)
```

Hardening storage and numerical loss do not become heat. Closed radiation's
roundoff-net transfer is reported but excluded from external supply. Separate
closed thermal and source-deposition checks prevent one cancellation in the
combined total from hiding a second missing heat term. Dead-force and prescribed
reaction work both remain with the original mechanical owner. Refusal, source
expiration, work-budget exhaustion, or cancellation returns no candidate.

## Boundaries and verification

This is the existing small-strain AND small-rotation, quasistatic, fully solid
energetic-hardening reduction. E, Poisson ratio, H and reference density stay
constant; the original supported isotropic expansion/reversible-heat model may
be used. C/L, electrical paths, thermal areas, view factors and fixed source
conditions are unchanged within an interval. Contact pressures and heat splits
remain explicit inputs; this does not supply contact mechanics, field forces,
piezoresistivity, melt flow, circuit-topology migration, or inertial dynamics.

Thermal transport is backward Euler with endpoint material/mechanical response;
electrical current retains the existing midpoint update. Coupled temporal
accuracy is generally first order, not an unconditional stability claim. Split
intervals at events and load/source discontinuities. Target events remain pending.
No explicit-observer substeps or observer history are silently imported into this
implicit entry point.

Five integration regressions under `solid_plastic_material_transport::heat`
exercise an independent two-tet/J2/RC/BE reduction, resistivity-only substitutions,
positive/negative charge with actual interface heating, original native circuit
and mechanical replay, continued unloading, final-evaluation budget rollback,
cancellation and source-domain failure. The existing thermal fixture is reused.

```sh
rch exec -- cargo test -p fs-couple --test solid_plastic_material_transport
```

Native Rust compilation/tests were not executed in the editing environment:
RCH, DSR, Cargo and rustc were unavailable. Independent Python/SciPy equations
passed 54 synthetic parameter/sign/contact cases with maximum combined-energy
discrepancy 8.33e-17 J. Those are equation checks, not Rust execution, temporal
convergence evidence, experimental calibration, or validation of rollback.
