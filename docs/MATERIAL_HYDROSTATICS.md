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
5. `step` reuses `RigidBodyIntegrator::step_verlet_with` and recomputes the load at
   the actual drifted pose for the second half-kick. It returns a state only after the energy and
   cancellation checks succeed. Commit `step.dynamics.state_after` on success.

For an upright uniform cuboid of height H and body/fluid densities rho_b/rho_f,
its partially immersed equilibrium draft is H*rho_b/rho_f when 0<rho_b<rho_f.
Changing either density changes that result without changing a force equation.
Changing the distribution of solid material can also move the center of mass
and alter torque without changing the displaced volume or buoyancy center.

## Initialize a floating or submerged equilibrium

`solve_heave_equilibrium` finds the vertical placement where the original
reservoir's buoyancy balances the original body's weight. Its bracket is a
translation from the supplied reference pose, positive opposite gravity.
Attitude and lateral position stay fixed. Both momenta in the returned state
are zero: this is an initialization, not a damping step on a live trajectory.

```rust
use fs_mbd_material::hydrostatic::{
    HeaveEquilibriumControls, solve_heave_equilibrium,
};

let equilibrium = solve_heave_equilibrium(
    &medium,
    &mesh,
    body.frame(),
    reference,
    HeaveEquilibriumControls {
        lower_offset_m: -2.0,
        upper_offset_m: 2.0,
        maximum_force_defect_n: 1e-8,
        maximum_evaluations: 64,
    },
    || false,
)?;
let initial_state = equilibrium.state;
```

The result includes the actual load, signed residual in newtons, number of
load evaluations, and the external torque needed to hold the chosen attitude.
A nonzero `holding_torque_body_nm` means a freely rotating body will not remain
in that attitude. The solve does not establish roll/pitch balance or stability.
No-root brackets, exhausted work budgets, cancellation and native reservoir
refusals return errors without publishing an approximate state.

The same `HydrostaticReservoir` interface covers these existing owners:

| Reservoir | Density and fluid inventory |
| --- | --- |
| `HydrostaticMedium` | One source-derived density below a prescribed plane |
| `LayeredHydrostaticMedium` | Stable discrete layers, including an optional overlying fluid |
| `LinearDensityHydrostaticMedium` | Continuous stable linear density gradient |
| `DensityProfileHydrostaticMedium` | Source-bound piecewise linear density profile |
| `FiniteHydrostaticTank` | One conserved liquid inventory and a solved waterline |
| `FiniteLayeredHydrostaticTank` | Conserved immiscible inventories and solved interfaces |

The equilibrium solver reuses each owner's existing load evaluation and frozen
material snapshots. In finite tanks, every candidate body placement solves the
conserved inventory again; the waterline therefore rises with displacement.
An inner tank-volume tolerance must resolve the requested outer force tolerance.
A uniform body denser than every available fluid can have no floating root;
a fully submerged neutrally buoyant body can have an interval of equilibria.
The caller's bracket selects the region to investigate without claiming uniqueness.

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

The basic `HydrostaticMedium` is a prescribed, homogeneous, constant-density hydrostatic reservoir
with a planar surface perpendicular to gravity. Fluid density above the plane
is neglected. Uniform exterior pressure cancels on the closed envelope; below
it the gauge pressure is rho_f*|g|*depth. Buoyancy excludes body weight, which
the existing integrator applies exactly once. No waves, added mass, drag,
slamming/water-entry impacts, capillarity, flooding, compressible stratification
are included. Finite-container back-reaction is supplied by the separate finite
tank owners above. A gas snapshot is only a local
constant-density approximation here, not a gas-column EOS integration.

All tetrahedron cut cases integrate volume and first moments analytically.
Positive-volume slivers are not discarded by a wet-fraction cutoff. This is
floating-point geometry, not a certified enclosure or proof of a nonoverlapping,
watertight partition. The caller owns geometry coverage and rigidity.

The pressure potential is rho_f*|g|*integral(max(depth,0), envelope volume).
Its negative pose derivative is the buoyancy force/torque. `step` audits kinetic
energy + body gravity potential + this fluid-displacement potential with an
absolute joule allowance. Its second-order symmetric endpoint-load scheme is NOT
energy-exact. Passing a local defect gate does not establish global stability,
experimental accuracy or a general fluid-structure energy certificate. Timestep
refinement and the appropriate physical regime remain necessary.

A load pass traverses at most 65,536 cells and polls cancellation between cells.
A time step performs two such passes and one existing rigid-body step; failure
leaves caller state untouched. Fluid/material/envelope changes construct a new
model: the API does not invent a mass/energy-conserving live replacement event.

An equilibrium solve performs at most 128 load queries. Each query retains
its reservoir's own cell and waterline budgets and cancellation checkpoints.
The outer scalar solve allocates no additional heap storage. It can supply a
balanced starting draft without simulating an artificial settling transient.

The Rust regressions cover analytic cut moments, all simplex permutations,
source-only substitutions, heterogeneous torque, potential derivatives,
endpoint feedback, heave refinement, long-run energy and failed transactions.
Reference: R. Fitzpatrick, *Fluid Mechanics*, Hydrostatics, University of Texas,
https://farside.ph.utexas.edu/teaching/336L/Fluidhtml/node27.html .

## Hydrostatic initial state for resolved free-surface flow

For simulations using the actual sparse D3Q19 flow owner,
`MaterialFreeSurfaceScaling::free_surface3_hydrostatic` constructs a planar
liquid column in balance with the material solver's density-dependent
`rho * gravity` force. Select a supporting Cartesian wall, a number of complete
wet layers, an explicit absolute gas-pressure reference, and a maximum
dimensionless density excursion. At least one gas layer is required;
capillarity and wall wetting must be absent. Existing permanent solid tiles
remain stationary obstacles.

The initializer uses the native collision and gas-link boundary equations.
Writing `h = |g_lattice|/(2*cs^2)`, the wet layer nearest the surface has density
`1/(1-h)`, and each deeper layer multiplies that density by `(1+h)/(1-h)`.
This includes the surface population correction needed for discrete rest.
The ordinary `step()` retains the declared gas pressure, so there is no
artificial settling run or pressure callback to maintain this initial state.

The profile is weakly compressible numerical hydrostatics, with mass
`fill*rho_lattice*rho_source*dx^3`. The density perturbation is not a physical
liquid equation of state. At fixed physical geometry and gravity, decreasing
`dt/dx` reduces the pressure and buoyancy error relative to incompressible
hydrostatics. Spatial refinement alone at fixed `dt/dx` does not eliminate that
error. The existing fully submerged stationary-solid load owner measures
pressure traction from actual reflected populations; this initializer does
not add partially wetted or moving-solid traction.

```
cargo test -p fs-lbm --test material_free_surface3 planar_hydrostatics
```

These tests specify 64 native rest steps in all six gravity directions,
pressure-datum and zero-gravity invariance, mass retention, and pressure,
submerged force and torque convergence as compressibility decreases.
