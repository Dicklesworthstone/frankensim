# Material-driven hydrostatics and rigid motion

`fs_mbd_material::hydrostatic` connects the existing source-backed body and
fluid material owners to partial-submersion forces, moments and rigid motion.
It does not introduce a second density registry or rigid-body integrator.

## One material-to-motion path

1. `resolve_body` resolves the solid regions' density claims and compiles their
   mass, center of mass and full inertia. `principal_body` adapts that SAME
   specimen to the existing free-body API.
2. `HydrostaticMedium::new` takes an existing `FluidTransport` owner, uniform
   gravity and a point on the free surface. It captures that checked transport
   snapshot once and retains a borrow of its source. `ResolvedLiquidState`,
   source-resolved gases and explicit transport-only data use the same path.
3. `DisplacementMesh` describes the fluid-excluding envelope in the original
   body frame. A tetrahedral partition is accepted, or `cuboid` supplies six
   cells. It is separate from the solid mass partition: a sealed hollow vessel
   displaces its enclosed envelope, not only its material walls. Flooded voids
   must be omitted by the geometry owner.
4. `load` clips that envelope against the current surface and returns buoyant
   force, torque in the actual principal frame, center of buoyancy and pressure
   potential. Material-induced balance changes do not move the authored mesh.
5. `step` reuses `RigidBodyIntegrator::step_kick_drift` and recomputes the load at
   the candidate endpoint. It returns a state only after the energy and
   cancellation checks succeed. Commit `step.dynamics.state_after` on success.

For an upright uniform cuboid of height H and body/fluid densities rho_b/rho_f,
its partially immersed equilibrium draft is H*rho_b/rho_f when 0<rho_b<rho_f.
Changing either density changes that result without changing a force equation.
Changing the distribution of solid material can also move the center of mass
and alter torque without changing the displaced volume or buoyancy center.

## Actual source-resolved example

```
cargo run -p fs-mbd-material --example floating_materials -- 500 1000
cargo run -p fs-mbd-material --example floating_materials -- 750 1000
cargo run -p fs-mbd-material --example floating_materials -- 500 1200
cargo test -p fs-mbd-material --test hydrostatic
cargo test -p fs-mbd-material --example floating_materials
```

The numbers are authored synthetic densities in kg/m3, not measured materials.
Both solid and liquid data enter immutable material cards and their existing
resolvers. Geometry (2 x 2 x 1 m), initial pose, gravity and timestep are held
fixed between runs. Output is CSV: height, vertical speed, immersed volume,
buoyancy and accumulated signed NUMERICAL energy change. There is no hidden
water preset, independent mass override or artificial fluid damping.

## Physical and numerical boundaries

The medium is a prescribed, homogeneous, constant-density hydrostatic reservoir
with a planar surface perpendicular to gravity. Fluid density above the plane
is neglected. Uniform exterior pressure cancels on the closed envelope; below
it the gauge pressure is rho_f*|g|*depth. Buoyancy excludes body weight, which
the existing integrator applies exactly once. No waves, added mass, drag,
slamming/water-entry impacts, capillarity, flooding, compressible stratification
or finite-container back-reaction is included. A gas snapshot is only a local
constant-density approximation here, not a gas-column EOS integration.

All tetrahedron cut cases integrate volume and first moments analytically.
Positive-volume slivers are not discarded by a wet-fraction cutoff. This is
floating-point geometry, not a certified enclosure or proof of a nonoverlapping,
watertight partition. The caller owns geometry coverage and rigidity.

The pressure potential is rho_f*|g|*integral(max(depth,0), envelope volume).
Its negative pose derivative is the buoyancy force/torque. `step` audits kinetic
energy + body gravity potential + this fluid-displacement potential with an
absolute joule allowance. Its first-order Lie kick-drift scheme is NOT
energy-exact. Passing a local defect gate does not establish global stability,
experimental accuracy or a general fluid-structure energy certificate. Timestep
refinement and the appropriate physical regime remain necessary.

A load pass traverses at most 65,536 cells and polls cancellation between cells.
A time step performs two such passes and one existing rigid-body step; failure
leaves caller state untouched. Fluid/material/envelope changes construct a new
model: the API does not invent a mass/energy-conserving live replacement event.

The Rust regressions cover analytic cut moments, all simplex permutations,
source-only substitutions, heterogeneous torque, potential derivatives,
endpoint feedback, heave refinement, long-run energy and failed transactions.
Reference: R. Fitzpatrick, *Fluid Mechanics*, Hydrostatics, University of Texas,
https://farside.ph.utexas.edu/teaching/336L/Fluidhtml/node27.html .
