# fs-em contract

## Purpose and layer

L3 material-resolved electromagnetic fields. The current owner is a linear,
isotropic, fixed-state planar magnetic vector-potential problem:
`-div(mu^-1 grad(A_z)) = J_z`, `B = (d_y A_z, -d_x A_z)`,
`H = B/mu`. Coordinates and extrusion thickness are SI metres, potential is
Wb/m, impressed current is A/m2, and permeability is absolute H/m.

## Owners, public API and invariants

`PlanarMagnetostatics` owns an immutable `TriComplex2`, one complete resolved
material state per triangle, and the existing FEEC stiffness. Permeability
uses `fs-material::field_storage::linear_magnetic_permeability`, the exact
admission shared with source-derived circuit inductance. Unsupported
field/frequency dependence or operating-field bounds refuse. A different
material or temperature snapshot requires a new prepared problem.

The existing weighted Whitney one-form Hodge and incidence assemble
`d0^T M(1/mu) d0`. Impressed cell current enters exact P1 load moments.
Prescribed potentials must cover complete exterior edges and anchor every
connected component. Remaining exterior edges have zero natural tangential H.
The existing checked CG owner supplies the recurrence and a recomputed true
Euclidean residual; exhausted convergence never publishes a field.

`PlanarField` borrows its exact prepared owner and retains nodal potential,
cell B/H, volume-integrated B.H/2 energy and the solver report. Fields remain
separate on each side of a material interface. Constant local gauge offsets
are subtracted before differentiation.

`PlanarFrame` admits a finite right-handed orthonormal 3-D placement.
`CellSample` names an exact source triangle and consumer position in metres.
`sample_cells` validates barycentric and extrusion support before expressing
B in the consumer's global frame. It never guesses a nearest cell. The output
feeds the existing LBM moving-conductor/current/heat/force owner.

## Errors, work and execution

`MagneticError` distinguishes shared material/assembly refusals, unsupported
inputs, cancellation and an actual nonconvergence report. Assembly limits,
explicit sample limits and solver iteration controls bound work. Polls occur
through source admission, assembly, load construction, component traversal,
Krylov boundaries, field recovery and transfer. Failure leaves immutable input
owners unchanged and publishes no partial output.

Traversal and arithmetic use fixed mesh/source storage order and the existing
numerical owners. No cross-platform bitwise identity or certified roundoff
enclosure is claimed. There is no unsafe code or feature flag.

## Conformance and boundaries

Synthetic regressions cover affine fields and energy, layered permeability and
tangential-H continuity, impressed-current scaling, manufactured-field
refinement, temperature-source refresh, component/boundary admission, source
model refusal, cancellation and explicit geometry/frame transfer. The actual
two-step LBM regression changes only one region's permeability and checks
computed B, current/Joule heating, Lorentz force and evolved momentum/thermal
populations. Sources parse with the pinned Rustfmt; native compilation and
regression execution remain separate verification requirements.

This is a discrete P1 estimate on a finite planar domain with fixed linear
material state. It does not model axisymmetry, open/infinite boundaries,
permanent magnets, nonlinear or hysteretic ferromagnetism, magnetic loss,
induction feedback, magnetization forces or thermal field-work coupling.
The existing flow connection treats the solved field as imposed in its
low-induction regime. Point transfer is not a conservative remap and does not
certify an arbitrary target grid's discrete magnetic divergence. Numerical
fixtures are not empirical material or device validation.
