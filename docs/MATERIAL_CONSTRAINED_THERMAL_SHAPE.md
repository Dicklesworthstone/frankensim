# Material temperature -> constrained shape, stress and conserved body inertia

`fs-mbd-material::deformation::TetrahedralMass` binds a material mass inventory
to a shared-node tetrahedral mesh. Its `solve_thermal` method, enabled by
`thermoelastic-coupling`, composes the existing material source resolver,
`fs-solid::linear3` equilibrium/stress/geometry operations and the existing
rigid-body inertia owner. No new thermal expansion, elasticity or dynamics law
is introduced.

## Physical input and use

Resolve each reference tetrahedron as one `RegionBinding` using
`VolumeMoments::tetrahedron` on its actual vertices. Use those SAME vertex and
connectivity arrays to bind the mesh:

```rust,ignore
let mesh = TetrahedralMass::bind(
    &reference_body, &nodes_m, &tetrahedra, &element_region_keys,
    max_nodes, || cancel_requested(),
)?;
let result = mesh.solve_thermal(
    &thermal_cells, &fixed_displacement_dofs, limits, &cx,
)?;
let current_inertia = result.body().inertia();
let current_nodes = result.geometry().nodes_m();
let stress_and_support_forces = result.stress();
```

Each `ThermalCell` names its original region and material card, current absolute
temperature and explicit source selections for stiffness and expansion. Input
order does not change the reference region order or the FEM connectivity order.
The complete original query, including pressure, composition and quantity-axis
metadata, is restored. Only temperature changes. Missing, repeated or foreign
assignments are errors, never permission to retain a stale cell.

The material owner resolves current Young's modulus and Poisson ratio and
integrates alpha from that cell's own reference temperature. The existing FEM
owner solves compatible displacement under the supplied zero-displacement
supports, checks its true equilibrium residual, and recovers actual elastic
stress, stored elastic energy and nodal internal/support forces. It then checks
strain bounds, nodal motion and element geometry before publishing the mesh.

Adjacent elements share vertices. A heated cell bonded to a cold cell therefore
cannot simply scale independently: their common face must have one displacement.
A fully constrained heated body stays the same shape and develops stress instead
of being fictitiously enlarged. Modifying only a source stiffness or expansion
law changes that physical response through the same solver.

## Same matter, different geometry

Reference density defines reference mass. The static FEM uses that Lagrangian
reference density, not a newly queried density on an unchanged reference volume.
Only current E and nu are requested for stiffness, so a density source supported
at the reference temperature alone does not block this closed-mass calculation.

Each P1 tetrahedron has an affine deformation map. Its actual updated vertex
coordinates supply analytic geometric moments through `VolumeMoments`; its
original mass is copied exactly into `SpatialInertia`. Current density is mass
per current volume. Complete tensors and centers are combined using the original
parallel-axis owner. This handles nonuniform deformation and shear, not just
uniform expansion. `body().inertia()` can be supplied to existing rigid/articulated
response calculations; `stress()` remains available to structural consumers.

`deform(current_nodes, ...)` also accepts compatible prescribed geometry from
other solvers. That route is a mass/geometry operation, not proof that the
supplied deformation satisfies any constitutive or equilibrium equation.

## Scope and checks

This delivery uses one reference region per tetrahedron, bounded by the existing
4096-region cap, fixed topology and per-cell uniform temperatures. It solves
static small-strain isotropic thermoelasticity. The source/model domain and
explicit strain/displacement/quality limits bound applicability; no yield or
finite-strain extrapolation is inferred. Temperatures can come from a thermal
solver, but this method does not itself evolve the heat timeline.

The resulting stress field is on the reference FEM chart. It is not silently
inserted into a geometric-stiffness or modal operator. A subsequent rigid-body
probe is not a live momentum transfer or a coupled thermomechanical time step.
This does not certify overlap freedom, solve topology changes, change mechanical
joint anchors, implement melting, or establish total thermomechanical energy
conservation. The previous uniform free-expansion and heating APIs are unchanged.

Focused tests:

```sh
cargo test -p fs-mbd-material --test deformation
cargo test -p fs-mbd-material --features thermoelastic-coupling --test thermoelastic_body
```

Tests use synthetic material data, independent volume quadrature and analytical
free/fully-constrained/bonded-cell references, and exercise the actual articulated
dynamics owner. They are not empirical accuracy claims for a named material.
