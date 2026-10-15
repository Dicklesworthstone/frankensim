# Accepted material heat and stress -> the existing vibration pencil

`ThermalBody::initial_stress_pencil` connects an accepted source-bound thermal
shape to the original solid and modal owners. It requires only the execution
context: no new density, elastic coefficients, thermal strains, stress array or
support list can be supplied independently and accidentally describe another
physical state. Enable `fs-mbd-material`'s `thermoelastic-coupling` feature.

```rust,ignore
let body = mesh.solve_thermal(&cell_temperatures, &fixed_dofs, limits, &cx)?;
let vibration = body.initial_stress_pencil(&cx)?;
let km = vibration.pencil();
// km.stiffness, km.mass and km.free_dofs are the original fs-solid carrier.
// Feed them directly to the existing fs-modal spectral APIs.
```

The same method is available on the latest accepted shape of a solved thermal
history. No second heat timeline, shape solve or material query is needed:

```rust,ignore
let (hot, heat_report) = heated.advance_to(time_s, step_config, max_shapes, &cx)?;
let vibration = hot.configuration().initial_stress_pencil(&cx)?;
let actual_stress = vibration.stress();
```

Changing only heat capacity changes the temperature reached under a prescribed
heat input. The original material owner then supplies thermal strain and current
stiffness; compatible deformation determines stress, and stress changes the
incremental vibration operator. A restrained heated specimen need not behave
like an independently enlarged, unstressed specimen. The same input semantics
apply to structures and machine components; no instrument or material name
chooses any physical equation.

## One physical linearization, with an explicit scope

The solid owner assembles the conservative initial-stress approximation

    K_t = K_e + K_g
    K_g[(a,i),(b,j)] = delta_ij integral grad(N_a)^T sigma_0 grad(N_b) dV.

`K_e` is the existing elastic matrix with current source-resolved elastic data.
`sigma_0` is recovered from the accepted constrained thermal equilibrium. The
same source/mesh/thermal/constraint identity checks used by stress recovery must
pass before any matrix is returned. The full symmetric Mandel stress tensor,
including shear components, participates. Shape gradients use the existing
batched linear algebra owner; mass and elastic assembly stay in `linear3`.

This is a **small-displacement reference-mesh initial-stress approximation**.
It keeps the original consistent mass matrix, mesh gradients and support DOFs;
it does not mix stress on a reference chart with gradients on an independently
scaled chart. The separately retained `body()`/`geometry()` still describe actual
deformed volume, center of mass and inertia. Their updated coordinates are not
silently substituted into this spectral approximation. This is not the full
finite-strain tangent, a follower-load linearization, or an estimate of every
geometric effect on resonance. The caller's strain, source and geometry limits
continue to bound the underlying static state.

The weak form is the classical initial-stress form in Bleyer/Pierre,
*Computational Mechanics Numerical Tours*, “Linear Buckling Analysis of a 3D
solid”, equation (32). There is no new material law or frequency correction.

## Signed stability, not invented audible frequencies

`geometric_stiffness()` and `elastic_stiffness()` expose the separate matrices.
Tension and compression can contribute opposite signs. Neither the geometric
matrix nor the combined matrix is required to be positive definite. The spectral
owner must see the actual signed eigenvalues: taking an absolute value or
clamping a nonpositive eigenvalue into an audible frequency would conceal an
instability. Only a positive lambda admits frequency sqrt(lambda)/(2*pi).
Positive eigenvalues alone do not establish continuum stability or accuracy.

All-fixed specimens retain their stresses but have no vibration degrees of
freedom; this method returns the existing `NoFreeDofs` error rather than releasing
supports. Changing a support or equilibrium requires a new physical solve, not
editing a saved pencil. Cancellation/refusal returns no partial pencil and does
not change an accepted thermal history. Old snapshots retain their own operator
inputs after later heating or cooling.

This API does not integrate vibrating displacement, map mode amplitudes across
changing bases, account for work in a live stiffness change, or generate pressure
or audio by itself. Damping, radiation and force projections remain explicit
consumers. The heat report is still a heat-only report, not a total thermoelastic
energy certificate.

Focused native regressions:

```sh
cargo test -p fs-solid --test thermal_prestress
cargo test -p fs-mbd-material --features thermoelastic-coupling --test thermoelastic_body
```

These use synthetic source claims, analytical restraint/gradient identities and
the actual existing generalized eigensolver. They are not measured-timbre or
named-material accuracy assertions.
