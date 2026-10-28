# Source-integrated planar layer conduction

`fs_material::thermal::conduction` supplies a constitutive flux closure, not a
second thermal integrator. For a homogeneous source-free planar layer,
`j*d = integral(T_cold..T_hot) k(T) dT`. The implementation integrates the
original piecewise-linear conductivity chart exactly on each crossed segment.
It does not replace the layer by its endpoint k or average two material names.

`TwoLayerConduction` solves two such layers and an optional area-specific
contact resistance `R` together. The contact jump is `j*R`, while total heat
is `Q=A*j`. A bounded safeguarded Newton solve checks both independently
integrated layer powers and the actual contact jump. Kelvin representation
roundoff is reported separately: dividing it by a vanishing resistance is
not a physical heat-balance defect and does not relax the layer power gate. Swapping layers and
endpoints reverses Q. Positive, finite k and fixed phase/composition are
required; an insulating layer belongs to the existing explicit insulation lane.
Both charts must support the full boundary-temperature interval.

The response includes both interface temperatures, actual heat residuals,
iteration count, a positive secant `G=Q/(Ta-Tb)`, and the physical endpoint
derivatives. At equilibrium the exact series-conductance limit supplies G.
For interface conductivities `kia,kib`, define
`D=da/kia+db/kib+R`. Then
`dQ/dTa=A*k(Ta)/(kia*D)` and `dQ/dTb=-A*k(Tb)/(kib*D)`.
These derivatives differentiate the physical equations, not solver iterations.
The secant is a flux lowering and MUST NOT be substituted for this tangent.

## Source preparation and original enthalpy storage

`conduction::source::ResolvedConductionLayer::new` binds an original
`EquilibriumEnthalpyPhaseCurve`, `ThermalMaterialPoint`, conductivity key,
explicit `ClaimSelection`, pure phase and temperature span. The existing
`query_line_typed` owner checks continuous source support, fixes all other
conditions, and retains every original breakpoint. The existing static-property
admission rejects undeclared operating axes and frequency-only sources.
The requested source span must remain in the caloric chart's solid or liquid
branch. A melting plateau is not admitted as solid merely because T is unchanged.

The prepared chart performs no source selection, receipt creation or caloric
state hashing during nonlinear flux trials. Complete original source lines are
shared with accepted payloads through `Arc`; no global cache is introduced.
Source uncertainty is retained, not converted into a numerical accuracy claim.

`conduction::network::PreparedLayerConduction` binds any explicit subset of an
existing `DepositionSimulation`'s links to these prepared material layers.
The original link area, two-sided surface coverage, cell geometry, material IDs
and actual borrowed caloric objects are checked. Reprepare after changing those
bindings. The original h, mass and clock may evolve normally.

Use `prepared.resolve(candidate, config, checkpoint)` as the coefficient
callback to `BackwardEulerStep::evaluate`. Its `StepCoefficients` replace only
the selected link conductances and contain no additional source watts:

```rust
let trial = step.evaluate(
    &cell_energy_increments_j,
    |candidate| prepared.resolve(candidate, layer_config, || Ok(())),
    || Ok(()),
)?;
let (next, layer_samples, energy_audit) = trial.accept(
    local_energy_tolerance_j, total_energy_tolerance_j, || Ok(()),
)?;
```

Application code should supply its cancellation checkpoint in place of the
example no-op closures. The outer nonlinear owner chooses and solves for the
energy increments; a small iteration update is not acceptance. Every trial must
resolve its own endpoint coefficients. Existing coupled solvers can consume the
same returned conductance map while retaining their original source map and
payload. Do not also add Q as a source: the original heat-rate owner already
transfers it exactly once on the two sides of the link. Untargeted links and
all authored source, convection and reservoir boundaries remain unchanged.

## Model and verification boundaries

This is a quasi-steady planar subcell closure. It adds no transient storage
inside a layer, resolved multidimensional conduction, anisotropy, radiation,
distributed layer generation, phase-front tracking or evolving geometry.
The thermal integrator's time error and the layer solve's watt tolerance remain
separate. A fixed-pressure/composition source line is not a reacting material law.

Six `layer_conduction` tests cover analytic slab limits, nonconstant k,
contact jumps and the near-perfect-contact limit, small integrals, derivatives
and refusal/cancellation. Four
`layer_conduction_network` tests use the original material-source and thermal
residual owners, including an independent two-cell backward-Euler endpoint.
Independent Python checks cover 256 randomized configurations and endpoint
finite differences, plus 64 independently solved 70-digit affine-layer systems.
Native Rust compilation and execution were unavailable during implementation;
these tests are not evidence of a completed workspace/DSR gate.
