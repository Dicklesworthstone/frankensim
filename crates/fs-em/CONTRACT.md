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
spatial eddy-current/flow induction feedback, magnetization forces or thermal
field-work coupling.
The existing flow connection treats the solved field as imposed in its
low-induction regime. Point transfer is not a conservative remap and does not
certify an arbitrary target grid's discrete magnetic divergence. Numerical
fixtures are not empirical material or device validation.

## Distributed winding and circuit connection

`winding::PlanarWindingBasis` admits signed trianglewise `J_z / terminal_current`
profiles [1/m2], including turn count and outgoing/return orientation. Profile
integrals use cross-sectional area; flux linkage and energy use the mesh's
actual extruded volume. Each profile must close separately on every connected
mesh component to the caller's current-balance tolerance. Closure residuals
remain available. Unit-current fields use homogeneous prescribed potential
and the existing planar boundary admission and checked solver.

The active-section inductance is the symmetric field-energy Gram matrix
`L_ij = integral(B_i.B_j/mu dV)`. Every entry must agree with the original P1
`integral(A_i J_unit,j dV)` pairing within an explicit absolute tolerance.
Original fields, material snapshots and solver reports remain reachable.
Callers separately declare nonnegative series end/leakage self-inductances;
zero explicitly neglects those ends. These entries contribute to total flux
and stored energy, but not to the represented cross-section field. End mutual
coupling is unsupported. Singular total inductance refuses; no numerical
leakage is invented. The existing Cholesky owner constructs the inverse, whose
measured `L L_inverse - I` row-sum residual must pass the declared threshold.

`winding::circuit::WindingCircuit` binds every original graph inductor exactly
once through explicit branch indices. Positive winding current follows that
branch's first-to-second terminal orientation. Total `lambda = L i` replaces
placeholder branch inductance; the full inverse replaces only the original
flux storage block. Capacitor charges, algebraic rows, Kirchhoff stamps and
source ports retain their original correspondence. The existing `CircuitDae`
owns consistent initialization and audited evolution. Circuit-state currents
reconstruct the same field through the existing triangle/position/frame
transfer; resistor branches continue to own dissipation.

Extraction admits at most 16 windings with explicit basis-storage and numerical
budgets. The dense circuit adapter admits at most 32 nodes, 64 branches,
16 transformers, 32 source ports and 128 extended states. Cancellation or
refusal publishes no partial owner. Regressions cover actual mutually coupled
circuit evolution, an interleaved capacitor and reversed branch mapping,
source-work/heat/storage accounting, material-only substitution, winding sign
and turn scaling, extrusion depth, end storage, gauge closure and refusals.
Native test execution remains a separate verification requirement. This lane
keeps geometry and material state fixed during evolution; it does not infer
nonlinear magnetics, motion, skin/proximity losses, changing-material work or
three-dimensional winding-end fields.
