# fs-em contract

## Three-dimensional distributed windings and native mutual circuits

`spatial::vector_potential::winding::SpatialWindingBasis` connects complete
closed 3D unit-current profiles to the original source-bound Whitney field.
Profiles include physical turns, orientation and return paths in the original
mesh frame. Fixed induced-linear materials and homogeneous tangential magnetic
traces admit one native unit-current solve per winding; original local current
continuity, complete relative curl-kernel compatibility, gauge and source
checks remain authoritative. Current cannot disappear through a prescribed
boundary. Preparation explicitly bounds winding count and retained field data.

Every self/mutual entry compares `integral(B_i.H_j dV)` against the original
Whitney pairing `integral(A_j.J_unit_i dV)` in both directions. The complete
reciprocal matrix includes source geometry and all meshed end effects; only
explicitly declared omitted series self-inductance is added. There is no
empirical coupling coefficient or invented leakage for dependent modes.

`SpatialWindingCircuit` replaces selected original inductors through
`assemble_circuit_with_mutual_inductance`. It preserves caller winding order,
original branch/source indices, unrelated circuit stores and transformers.
Consistent initialization, actual coupled currents and audited stepping use
the native DAE. Full positive-definiteness admission remains with that owner;
dependent modes need physical leakage or an explicit ideal-transformer model.

Ordinary circuit steps reuse the small native system. Optional field
reconstruction performs no new 3D solve: original curl, stiffness, source and
constitutive operators recover potential, face flux, B/H, current and cell
energy at the actual currents. Current compatibility, Ampere residual,
divergence/conforming-flux and work allowances are rechecked at that amplitude.
Declared external inductive energy is separate from meshed field energy;
both already belong to the native winding storage and are never added twice.

Geometry, material sources, winding profiles and boundary traces stay fixed.
This linear reduction does not imply saturation, remanence, skin/proximity
diffusion, thermal evolution, moving winding transport or an infinite exterior.
G0/G1/G3 fixtures in `tests/spatial_winding.rs` independently normalize a compact
closed current loop, compare direct and reconstructed fields, check gauge,
polarity and material response, and solve a native mutual RL circuit against
independent coupled equations with unrelated LC storage. Return, source, budget,
dependent-mode and cancellation refusals preserve the native physical gates.
These authored fixtures do not constitute measured-device qualification.



## Driven solid-winding motional EMF and reciprocal Lorentz work

The nonlinear moving field and native winding circuit now admit driven solid
sections through their original current-constrained regions. With
`q = delta(A)/dt + v.grad(A)` and `u = n I/integral(sigma dA)`, the same
projected conductivity operator solves the magnetic residual and reconstructs
`E_conductor = u - P_sigma(q)`. Actual total `J = sigma E_conductor` supplies
both exact P1 Joule heat and `J cross B` work. The conductor-frame heat,
laboratory electric field and motional electric field remain distinct receipts.

The native terminal law adds the source-derived back-EMF
`V_motion = depth*n*mean_sigma(v.grad(A))` to direct
`F_new^T delta(A)/dt`, explicit end-inductance voltage and original DC loss.
This voltage is computed directly even at zero current. Separate DC and signed
motional drive work reconcile the complete field source work with the circuit;
neither creates an extra heat deposit. Negative mechanical work is supplied
by the imposed-motion owner, while physical Joule heat stays nonnegative.

`step_with_eulerian_velocity_and_sources_candidate` refreshes only original
conductor/resistor sources, retaining field/circuit history, magnetic laws,
end storage and topology. Solid profiles follow the refreshed conductivity;
an unchanged magnetic potential does not acquire artificial transformer EMF.
The sparse nonlinear operator, nonsymmetric FGMRES solve, true residual,
Peclet limits and original energy/cancellation gates remain authoritative.

This is prescribed Eulerian motion with invariant material occupancy and
coefficient fields along the motion. It does not translate finite interfaces,
move permanent magnetization, or account for thermal advection. Stranded cells
and fixed-recoil magnets remain stationary. The conductor-frame convention and
transport restrictions follow the original
[quasistatic Lorentz formulation](https://doc.comsol.com/6.3/doc/com.comsol.help.acdc/acdc_ug_theory.05.09.html)
and [Eulerian velocity scope](https://doc.comsol.com/6.4/doc/com.comsol.help.acdc/acdc_ug_magnetic_fields.08.040.html).

Focused G0/G1/G3 fixtures in `tests/nonlinear_winding_solid/moving.rs` cover
zero-motion equivalence, open-current motional voltage, motor/generator work,
total spatial heat and force, original conductivity refresh, and unchanged
histories on source/budget/cancellation refusal. They are authored numerical
fixtures; native execution and empirical device qualification remain separate.



## Causal solid-winding diffusion in native nonlinear circuits

`winding::nonlinear::NonlinearWindingCircuit::new_with_solid_regions` and its
explicit-port variant drive whole connected solid regions through the original
native terminal solve and `PlanarNonlinearEddy` history. The accepted spatial
potential retains causal skin/proximity memory under finite-amplitude drive
and source-defined magnetic saturation. The existing sparse magnetic operator,
projected consistent P1 conductivity mass, winding/end storage and native
`CircuitDae` remain authoritative; no terminal-impedance realization is fitted.

Harmonic and transient paths share `planar::eddy::solid::SolidWindingRegion`
and the original source normalization. With `S = integral(sigma dA)`, each
whole signed section enforces `integral(J dA) = n I` using
`J_unit = n sigma/S`, `R_DC = depth*n^2/S`, and
`E = n I/S - (delta(A)/dt - mean_sigma(delta(A)/dt))`. One region belongs to
one winding coordinate; explicit return sections, floating conductors and
disjoint stranded copper retain their existing physical roles.

Conductivity refresh updates the source functional and DC drive from the same
original cells while preserving accepted potential, current, circuit and clock.
Actual induction is `F_new^T delta(A)/dt`, plus explicit end-inductance voltage.
An endpoint-linkage change from reweighting a stationary field creates no EMF.
The old potential paired with new profiles is prepared once per terminal solve.
`inductive_voltage_v` reports physical induction; endpoint linkage remains a
diagnostic that can change when the conductor profile changes.

The original concrete field receipt carries actual total local E, J and
physical Joule heat, including uniform/circulating-current cross terms.
`uniform_drive_work_j` is already included in supplied field work and physical
heat; it is never an additional deposit. Stranded heat remains disjoint, and
backward-Euler magnetic Bregman loss remains numerical dissipation. Independent
terminal, field, end-storage and global work gates precede publication. Source,
topology, budget and cancellation refusals preserve accepted histories.
The stationary entry points retain fixed geometry; the Eulerian moving
extension above retains explicit material-transport restrictions.

G0/G1/G3 tests in `tests/nonlinear_winding_solid.rs` and its `slabs` module
compare bipolar diffusion/saturation pulses, native loaded-circuit recurrences,
stationary biased-field conductivity refresh, and two-grid skin/proximity
profiles against independent scalar and continuous-slab oracles. These authored
fixtures are not experimental material qualification; native execution requires
the workspace dependency repositories.



## Bias-dependent winding impedance with resolved eddy currents

`planar::eddy::harmonic::PlanarHarmonicWindings` prepares frequency-domain
small-signal responses from an admitted `PlanarIncrementalWindingBasis`.
The original DC operating point, signed distributed winding profiles, complete
nonlinear differential magnetic Hessian, source temperature, geometry, remanence
and boundary trace remain authoritative. Both longitudinal `dH/dB` and transverse
`H/B` enter the field tangent. Phasors are peaks with `exp(-i omega t)`, matching
`fs-phs::harmonic`; inductive impedance therefore has negative imaginary part.

Every unit winding solves `(K_t - i omega C_sigma) a_j = f_j`. `C_sigma` is the
existing consistent P1 conductivity mass minus each isolated region's weighted
mean projection. It enforces zero integrated axial current independently in
each floating conductor, with original topology/material admission. The sparse
Hessian and matrix-free regional projection feed the existing real-equivalent
block and cancellable FGMRES owners. The original diffusion/storage Jacobi
preconditioner is reused on both real components, with explicit restart, cycle,
field-payload and Krylov-storage limits. Full-mesh operator buffers are reused;
no dense field matrix is formed. Conductor-free frequencies reuse DC derivative
fields only after their true residual passes the current AC allowances.

Original source-paired linkage gives the complete reciprocal terminal matrix
`Z = R_copper - i omega (F^T a + L_end)` and effective inductance `-Im(Z)/omega`.
Every active stranded triangle binds its complete original conductor source;
its resistance integrates `rho J_unit^2` over the actual extruded volume.
Stranded and resolved conductor cells must be disjoint. Explicit end/leakage
storage is retained; end resistance is not inferred. Resolved core/shield eddy
currents produce actual frequency-dependent loss and inductance. Stranded
profiles remain fixed; driven solid sections use the additive path below.
Fitted hysteresis loss and full-wave effects remain outside this reduction.

Acceptance audits complex reciprocity and the COMPLETE mutual power identity
`Z_ij = R_ij + omega^2 a_i* C_sigma a_j - i omega a_i* K_t a_j`, including
declared end storage. Arbitrary phase-shifted terminal currents reuse the unit
fields, then independently reconstruct original-cell B, projected E, J and exact
P1 cycle-mean Joule power. Real terminal power equals spatial copper plus eddy
heat; reactive power equals `-2 omega` times mean quadratic perturbation storage.
DC copper heat is separate. A caller-declared cell-diameter/skin-depth maximum
uses the larger local incremental permeability, including the transverse
direction. This spatial admission bound is not an error certificate.

`tests/winding_harmonic.rs` compares native responses with independently soluble
field modes across DC bias, frequency, conductivity and extrusion depth. It
covers signed mutual coupling, phase-dependent power, floating-current
neutrality, original heat sources, end storage, the no-eddy limit, affine bias
and gauge invariance, transverse tangent/skin depth, double-counting/topology
refusal, work budgets, nonconvergence and cancellation. Manufactured sources
exercise the implementation; they do not validate an empirical material/device.
Native test execution remains a separate requirement where workspace dependency
repositories are available. Finite-amplitude nonlinear cycles, material/geometry
evolution and causal transient realization of the terminal matrix are not claimed.

### Resolved skin and proximity currents in solid windings

`new_with_solid_regions` adds `SolidWindingRegion { region, winding, signed_turns }`
bindings to the same prepared owner. Each whole connected conductive section
has the signed constraint `integral(J dA) = n I_w`. With
`S = integral(sigma dA)`, eliminating its axial potential gradient gives
`J_unit = n sigma/S`, `R_DC = depth*n^2/S`, and actual peak field
`E = n I_w/S + i omega (A - mean_sigma(A))`. The existing sparse magnetic
pencil, projected conductivity operator, current-closure topology, original
source/flux pairing and complex mutual power identity are unchanged.

The retained DC winding profile must already match this source-normalized
distribution. Its volume-weighted absolute discrepancy, scaled by
`max(1, abs(I_DC))`, must satisfy both the caller's weak-load and integrated
current allowances. No operating point or profile is silently replaced.
One conductive region has one winding coordinate and a declared polarity;
every return section is explicit. Other regions remain floating. A winding
may combine solid and stranded sections without sharing source cells.

`Conduction::heat_with_region_drives` adds the uniform drive after the original
conductivity projection and reuses its exact P1 total-field Joule quadrature.
Every driven region checks its actual complex integrated current against `n I`.
Solid-cell heat includes the local uniform/circulating-current cross terms;
separate DC and eddy densities would give the wrong spatial thermal load.
Returned solid, stranded and floating source-cell partitions are disjoint.
`mean_dissipated_power_w()` sums their actual heat once for circuit/thermal
consumers. An unpowered solid winding can still develop resolved proximity
currents while its integrated current remains zero.

`tests/winding_solid.rs` exercises four separated slab sections with two
independent currents against `J'' = -i omega mu sigma J`, with interface data
from Ampere's law. It checks spatial-current refinement, independently
integrated heat, changing neighboring-current phase at fixed winding current,
the uniform low-frequency limit, local Joule cross terms, mixed solid/stranded
sections plus a floating shield, explicit returns, source-profile refusals and
unchanged accepted bias after cancellation. These are manufactured planar
MQS tests, not empirical winding validation. Axial end redistribution, individual
insulated strands inside an authored homogeneous region, full 3-D end turns,
and displacement-current regimes are not inferred.


### Native circuits drive and load the resolved winding fields

`HarmonicWindingFrequency::circuit_response` connects the same full terminal
impedance to original native R/L/C, mutual-inductance, ideal-transformer and
voltage/current-source constraints. Each field-owned inductor branch is replaced
temporarily by a positive reference resistor and a compensating native current
port. One `HarmonicModel` factorization supplies the open-source and unit-port
responses; a bounded terminal Schur solve and final native replay cancel the
reference resistance exactly. Its value is numerical conditioning, not physical
damping, storage or heat. Existing sparse unit field solutions are reused.

The result retains original branch/node/transformer identities, individual
source outputs, authored control-index gaps and complex phases. Repeated controls
expand to private native ports so their individual currents are not lost; one
control cannot mix voltage and current units. Original external mutual links
must not touch field-owned windings, whose complete magnetic coupling already
belongs to the field terminal matrix. Actual source constraints, transformer
work, Kirchhoff currents, native residuals and winding terminal consistency
are checked independently of the small interface solve.

Spatial solid, stranded and floating-conductor losses are recovered once through
the same field response. Original circuit resistor heat and unrelated native
storage remain separate. Supplied complex power closes against their complete
real and reactive account, using the distinct conjugation conventions for
voltage-driven and current-driven native ports. Reference resistors and private
compensating sources never appear as physical losses or external supplies.
Input graph, native dense payload, interface storage, solve work and cancellation
are bounded before publication; true poles and inadmissible constraints refuse.

Five authored-source regression groups cover loaded coupled windings and
reference-resistance invariance, mixed ideal current/voltage clamps, repeated
source controls, ideal transformers with external mutual RLC at series resonance,
and binding/budget/cancellation/true-pole refusal. The same fixed-bias planar MQS
and source-domain limits apply; native test execution is a separate requirement
where the workspace dependency repositories are available.

## Ionized-gas current and Joule heating

`PlanarConduction` accepts the sealed `PlanarOhmicResponse` material response,
including `fs-material::ionized_gas::IonizedGasConductivity`. The same weighted
FEEC assembly, boundary current solve, source-work audit, cell current and
`E.J` heat integration are used for source-property conductors and ionized gas
with actual frozen populations. Gas populations, separate electron/heavy
temperatures and collisional rates remain in the original
material snapshot. Contact construction and material refresh keep their
stronger source-material requirements; no gas/solid contact law is inferred.

Field recovery checks the material's admission against the actual electric
field before publishing current or heat, including circuit voltage recovery.
This prevents a weak-drift gas law from silently extending into breakdown.
Resistance and total-current preparation accept an explicit positive reference
voltage, so a valid small physical drive need not admit a one-volt computational
basis. The retained reference field is admitted at its actual voltage;
`G=P/Vref²=I/Vref` supplies resistance, and circuit recovery scales once by
`V/Vref` before independently checking the actual field, current and power.
The gas response is static, unmagnetized and local, with a prescribed massive
neutral collision background. It does not add a discharge, Hall, radiation,
space-charge or magnetohydrodynamic solver. Thermal, flow or finite-rate
chemistry feedback must carry its actual new populations and temperatures into
a new accepted snapshot and resolve collision rates there. A kinetic endpoint
is never projected back onto LTE merely to obtain conductivity.

The material's checked Joule partition consumes the actual recovered FEEC
electric field and supplies separate electron/heavy thermal powers. Each
power density uses that field cell's physical volume before entering the
conservative thermal-pool update. Their sum is the original field's `E.J`;
the native resistor heat is not added a second time. Updated temperatures feed
the same source-owned material and field owners.

G1/G3 tests in `tests/ionized_gas.rs` use analytic H/H+/electron thermodynamics
and explicit synthetic collision sources, then solve the actual native mesh.
Affine voltage, current density, Joule power, native circuit resistance and
voltage reversal agree with independent slab formulas at two temperatures.
G0 tests refuse actual PDE and circuit fields above the carrier drift domain.
Further G1/G3 tests carry an accepted finite-rate ionization endpoint through
actual neutral-density momentum drag, field current, Joule power and native
circuit resistance, and vary electron temperature with populations held fixed.

`tests/ionization_feedback.rs` carries separate-temperature monatomic atomic
ionization through that same path. Original NIST event rates and source-derived
atomic detailed balance change carrier amounts; native cell-volume Joule powers
are held explicitly during the chemistry interval. Independent monatomic
calorics audit electrical work, changing electron heat capacity and chemical
storage once. Refreshed populations and collision rates then change the native
current and circuit resistance. This is an ionization substep with a subsequent
field feedback evaluation, not a continuously fixed-voltage transient or an
empirically validated discharge; elastic exchange remains its separate owner.

`tests/coulomb_feedback.rs` uses the same native field/circuit path with
reciprocal Coulomb momentum transfer. The material resolves correlated carrier
velocities, retaining both its original background drag and every charged pair.
Original FEEC cell volumes convert its friction heating and signed
drift-dependent thermal correction to two held powers. The existing gas energy
owner consumes the identical population-bound Coulomb source together with the
original neutral momentum source, refreshes their zero-drift exchange at actual
temperatures, and closes both pool budgets against one electrical-work input.
Refreshed two-temperature material changes native current and resistance.
This remains an explicitly held-power energy interval followed by field
feedback, without a second PDE, a free plasma common-velocity gauge, or an
additional `E.J` charge to either thermal pool.

`tests/collisional_ionization_feedback.rs` combines the two paths in one
reacting energy evolution. Both an ionizing parcel and a recombining parcel
consume the original atomic impact/reverse law, neutral momentum source and
Coulomb screening policy while carrier amounts and both temperatures change.
Independent monatomic budgets include changing electron heat capacity,
original chemical storage, both internal exchange channels and the one held
native electrical-work input. The final populations invalidate the original
Coulomb transport binding; an explicit refreshed model retains its original
screening and domain limits. Tests check actual neutral-density drag, compare
the final coupled conductivity with an independent two-fluid inverse and a
counterfactual holding the old populations, then verify native field current,
Joule power and circuit resistance. No equilibrium projection or additional
material/PDE owner is introduced.

## Prescribed Eulerian moving-conductor induction

`PlanarEddy::step_with_eulerian_velocity_candidate` continues the original
accepted vector potential, impressed-current history, floating electrical
regions and clock with a prescribed continuous P1 planar velocity. It solves
`K A_mid + P_sigma(Adot + v.grad(A_mid)) = load(J_impressed_mid)` on the same
fixed mesh and source-bound linear material snapshots. The existing
conductivity-weighted P1 mass/mean projection enforces zero integrated axial
current separately in every isolated conducting region. Uniform motional EMF
is cancelled by the floating scalar-potential gradient, without fictitious
eddy current or Joule heating.

Velocity feedback makes the operator nonsymmetric. The existing cancellable
FGMRES owner supplies the recurrence, with a positive diffusion/storage Jacobi
preconditioner and a correct matrix-free transpose. Full-mesh work buffers and
cell geometry are reused during Arnoldi actions; neither the advection matrix
nor regional projection is densified. Explicit restart/cycle/Krylov-storage
budgets bound the solve. Acceptance recomputes the actual retained increment's
true Euclidean residual against the larger of caller-supplied absolute [A m]
and relative targets. A caller-declared cell magnetic Peclet limit
`mu*sigma*max_vertex_speed*diameter/2` bounds admission of this unstabilized
Galerkin spatial discretization; it is not an error certificate.

Returned original-cell fields distinguish laboratory E from conductor-frame
`E'=E+v cross B`, with `J=sigma E'`. The SAME stage B and P1 current supply
`J cross B`, exact quadratic Joule heat, and exact P1 Lorentz-mechanical work.
The accepted account is `delta W + Q_Joule + W_mechanical = W_impressed`;
positive mechanical work is delivered by the field to prescribed motion.
Thermal and mechanical receivers must close their own balances with these
reported transfers. These are free-conduction-current forces, not an additional
Maxwell-stress body load or a magnetization/interface-force model.

This is a nonrelativistic, fixed-geometry, fixed-coefficient Eulerian MQS lane.
Motion must preserve the represented material occupancy/coefficient field,
including homogeneous through-flow with separately supplied boundary data.
Translating finite bodies, moving material interfaces or magnetic sources,
ALE/sliding meshes, electrical terminals, nonlinear magnetic laws, remanence,
hysteresis and displacement current remain outside it. Prescribing velocity
does not advance mechanical state or infer its energy supply. The stationary
and temperature-refresh operations remain separate. Physical convention and
fixed-domain applicability follow the official COMSOL AC/DC guide,
[Velocity (Lorentz Term)](https://doc.comsol.com/6.4/doc/com.comsol.help.acdc/acdc_ug_magnetic_fields.08.040.html).

Focused G1/G3/G4 tests compare the spatially varying-velocity field with an
independently assembled nonsymmetric two-mode system and separate degree-two
quadrature for heat and mechanical work. They cover source/material/depth
substitutions, genuine field feedback, uniform-field charge cancellation,
zero-velocity stationary continuation and work/source/velocity/cancellation
refusals. Synthetic coefficients do not qualify an empirical material/device.

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
cell B/H, volume-integrated constitutive energy and the solver report. Induced
cells use B.H/2; recoil cells retain the fixed-remanence reference below. Fields remain
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
nonlinear or hysteretic ferromagnetism, irreversible demagnetization, magnetic loss,
spatial eddy-current/flow induction feedback, distributed magnetostriction or
moving-geometry field-work coupling. The bounded static net-force connection
is described below; it is not a dynamic electromechanical energy integrator.
The existing flow connection treats the solved field as imposed in its
low-induction regime. Point transfer is not a conservative remap and does not
certify an arbitrary target grid's discrete magnetic divergence. Numerical
fixtures are not empirical material or device validation.

## Static magnetic load to native mechanics

`PlanarField::maxwell_wrench` integrates the computed Maxwell stress
`T = B tensor H - (B.H) I/2` over the closed boundary of explicitly selected
mesh triangles. Those triangles include the target and an exterior collar.
`TriComplex2::trace_for_faces` remains the topology owner. Both cells incident
on every contour edge must share one complete resolved exterior material
snapshot, use an induced-linear law without remanence, and have zero original
impressed current. Enclosed recoil magnets remain supported by the original
field solve and its source-envelope operating checks. The contour cannot touch
the outer mesh boundary or run along a material interface; incident triangles
must also lie on opposite geometric sides rather than form a folded embedding.
The original static solve retains its actual cell-current field for this admission;
induction-history recovery alone cannot satisfy it. Incremental responses on
magnetized problems also refuse because they omit the fixed bias and carry no
actual recoil operating-point admission. Unbiased current responses remain usable.

The model assumes a homogeneous linear exterior medium without mechanical
pressure response or magnetostriction. Its permeability comes from the same
material query as the field. No material name selects a constant or force law.
Contrasting enclosed permeability changes the actual field and thus the load.
The contour's exterior-side P1 B/H values remain piecewise constant: edge
force and first moment are exact for those discrete values, while the physical
force remains a mesh- and contour-dependent numerical estimate. No gradient
recovery, interfacial averaging or error certificate is implied.

Actual extrusion thickness multiplies the surface load. Geometric triangle
orientation determines outward normals for either clockwise or
counterclockwise complexes. Force and torque use the caller's `PlanarFrame`;
the torque origin is explicit in that same frame. The returned `fs_ga::Wrench`
is the canonical native articulated-dynamics input, not another force type or
mechanical integrator. End-face contributions cancel in the invariant planar
model. End fringing and distributed magnetostriction remain unsupported; the
electric owner below uses the same geometric integration. The fictitious contour is not a physical body
surface, so its traction must not be used as an inferred local solid load.

This is an instantaneous static connection. Motion requires updated geometry
and field solves with explicit source, field and mechanical work accounting;
it cannot reuse a frozen field as a complete energy-conserving actuator. An
overlapping J-cross-B load must not also be applied to the enclosed body.
Original mesh-cell and edge limits bound contour traversal/allocation. Fixed
edge order and compensated summation preserve deterministic selection-order
behavior. Cancellation and admission failures leave the original field intact.

Focused G1/G3/G4 regressions exercise zero force in a uniform field, refinement
toward an independent Lorentz-force integral, thickness and source scaling,
source-only permeability changes propagated through native free-body
acceleration, frame/origin transformations, and unsupported contour/history
refusals. These synthetic calculations do not validate an empirical material
or device. The stress convention follows J. L. Kirtley, MIT 6.685, Class Notes 1,
section 4.2: https://web.mit.edu/6.685/www/chapter1.pdf.

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
branch's first-to-second terminal orientation. Native `delta_lambda = L i` replaces
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


## Differential winding inductance at nonlinear bias

`winding::incremental::PlanarIncrementalWindingBasis` first solves the actual
nonlinear DC field from signed spatial winding profiles and terminal currents,
with the original material laws, fixed recoil sources and prescribed potential
trace. The operating field must pass the native nonlinear residual gate.
Its source-derived constitutive Hessian is then assembled once and reused by
checked CG for every unit-current derivative; perturbations are zero on exactly
the original prescribed vertices. The mesh and tangent remain sparse.

The resulting symmetric terminal matrix is `L_diff = d(lambda)/di`, including
declared linear end self-inductances. It is not `lambda/i` or `2W/i^2` in a
saturated or biased core. Both radial `dH/dB` and transverse `H/B` tangent
eigenvalues are retained. Every entry of the Hessian-energy Gram matrix must
agree with the original `integral(dA/di_j * J_unit,i dV)` pairing within the
explicit henry allowance. The existing inverse-residual gate admits the small
terminal inverse. The original operating field, actual absolute flux, signed
profiles, boundary, source identities and true-residual reports remain reachable.
Voltage evaluates `L_diff di/dt` with all mutual terms and no voltage from a
fixed remanent flux offset.

Source H(B) corners refuse within a caller-declared tesla exclusion band;
maximum-flux and recoil-envelope edges refuse because a two-sided derivative
is not supported there. The band is not a field-error certificate. This is a
local, reversible, fixed-temperature and fixed-geometry DC derivative, not an
arbitrary-amplitude circuit law, changing-temperature work law, hysteresis
model, skin/proximity impedance or three-dimensional end-field calculation.

G1/G3/G4 regressions compare saturation with an independent two-mode analytic
solution, including the full transverse tangent under nonzero boundary bias;
check signed mutual/turn/depth scaling, explicit end storage and inverse action;
compare independent nonlinear current perturbations; retain permanent flux
without spurious induced voltage; and exercise source corners, field bounds,
work limits, singular modes, nonconvergence and cancellation. The formulation
follows [ANSYS Theory Reference §5.5](https://ansyshelp.ansys.com/public/Views/Secured/corp/v251/en/ans_thry/thy_emg6.html).
These manufactured sources do not qualify an empirical core or device.

## Spatial material winding resistance and native-circuit heat

`winding::ohmic::OhmicWindingCircuit` binds every active original magnetic
triangle to its complete `ResolvedConductorState`. The same signed unit-current
profile and actual extruded mesh volume give each winding's intrinsic series
resistance, `R_w = sum(rho_e * J_unit,w^2 * volume)`. Conductivity and resistivity
use the canonical scalar conductor owner; no metal names or preset coefficients
enter the solver. Current-carrying windings must occupy disjoint triangles;
unresolved mixtures/fill factors and skin, proximity, eddy and end-turn losses
are not inferred.

Resistance is stamped on the mapped native flux coordinates, conjugate to the
full mutually coupled `L^-1 delta_lambda` currents. All original magnetic/capacitive
storage, Dirac structure, ports, ideal constraints and numeric graph resistors
are retained. The existing descriptor midpoint step remains the only electrical
time integrator. Per-triangle `E.J` heat, every graph-resistor heat, integrated
series `I^2 R` work, native dissipation and total stored-energy/source-work balance
must agree within explicitly supplied absolute power/energy allowances. These
checks do not certify temporal or spatial approximation error.

Each candidate begins with the original accepted fluxes/charges and time; source
refresh, nonlinear thermal trials and cancellation cannot debit accepted history.
A refresh retains the complete original card, coordinate schema, requirements,
selected source claims and independent electrical source basis. Only the owned
absolute-temperature coordinate may change. Fixed external conditions and the
prepared magnetic permeability stay fixed. A temperature-dependent permeability
source refuses changed temperature even if its endpoint values coincide, since
fixed-L evolution cannot account for magnetic parameter work along that path.
Data-only material replacements are prepared as a new magnetic/electrical model.
`OhmicWindingCircuit::rebind_material_candidate` transfers accepted history to
that model at the same clock. It retains the exact mesh/metric, original
homogeneous magnetic boundary, signed current profiles, end-inductance model,
branch correspondence and all nonmagnetic graph data. Every capacitor charge is
preserved bit-for-bit. Total physical linkage is conserved by translating native
incremental linkage as `delta_lambda_new=delta_lambda_old+(lambda0_old-lambda0_new)`;
identical bias retains native linkage bit-for-bit without adding and subtracting
the larger permanent flux. The native consistency solver re-equilibrates
algebraic voltages/currents. Incompatible ideal constraints refuse instead of
creating an unaccounted flux or charge impulse. Both actual endpoint fields must
satisfy their source recoil envelopes before a replacement is returned.

Replacement adopts the new full mutual-inductance block and spatial conductor
sources, including a different admitted resistivity/conductivity source basis.
The signed external work is the change in `W0 + delta_lambda^T L^-1 delta_lambda / 2`
at fixed total linkage. Its incremental part is independently compared with the
native field/capacitor storage change before adding the bias reference, so large
W0 cannot mask electrical-accounting errors. The full signed work is checked
against the caller's absolute work bound. The report retains before/after bias
and total field/capacitor energy, current, resistance and source-state identities.
This is work in the declared ZeroFlux/FixedRemanence field potential; it does not
predict material formation energy, an intermediate replacement trajectory,
hysteresis or irreversible demagnetization. No time, source-port work or Joule heat is
advanced. Subsequent original midpoint steps and field sampling use the new
material basis; resistance-only replacements have zero instantaneous storage work.
Source, geometry, ideal-constraint, work and late-cancellation refusals publish
nothing. G1/G3 regressions in `tests/winding_ohmic/rebind.rs` cover analytic mutual
current/work, nonzero-time continuation with actual new spatial heat, reverse and
zero-work replacements, and incompatible geometry/constraints.

This is externally supplied material reconfiguration, not continuous thermal
evolution, motion, hysteresis, mass transfer or a caloric heat debit. Changing
geometry or automatically closing magnetic parameter work into thermal storage
still requires its own coupled owner.

## Temperature-refreshed spatial induction

`PlanarEddy::step_with_materials_candidate` retains the original mesh, A history, fixed boundary traces, electrical region identities, source choices and clock while assembling source-pinned endpoint permeability. It reuses the original checked CG with endpoint magnetic effort (K_old A_old+K_next A_next)/2; insulating/floating algebraic constraints are checked at the actual endpoint. The accepted candidate owns its rebuilt operator through a shared immutable owner; a field recovered from that candidate borrows the candidate lifetime.

Each triangle reports actual old/new curl A, magnetic free energy, volume and signed parameter work P=V B_old dot B_next (1/mu_next−1/mu_old)/2. The native account is delta psi + Joule heat = impressed-source work + P. Parameter work is not automatically heat: a thermal or external-work owner must close it. Fixed-coefficient steps retain P=0. No constitutive loss is added to resolved sigma E². Fixed geometry and linear isotropic nondispersive mu remain required; moving media, hysteresis and unconditional coupled entropy/passivity are not claimed. Synthetic tests cover analytical one-mode evolution, repeated operator/history continuation, algebraic insulating refresh, unchanged-source refusal, assembly budget and cancellation.

## Source-resolved reversible permanent magnets

`PlanarMagnetostatics::new_with_material_models` accepts induced-linear and `PlanarRecoilMaterial` cells on the same mesh. Each magnet supplies a complete `ResolvedRecoilMaterial` and an explicit unit direction in the original planar frame. The same FEEC stiffness and checked CG solve impressed current together with the exact P1 weak load `integral (Br/mu_recoil).rot(grad(phi)) dV`. Retaining this whole moment includes material-interface and natural-boundary source contributions even when Br is constant inside a cell. Unprescribed exterior edges continue to impose zero tangential physical H.

Recovery reports `H=(B-Br)/mu_recoil`, source-envelope operating margins, and cell storage `V*(B-Br).H/2`, with the `FixedRemanence` reference explicit. Induced cells retain `V*B.H/2`. A source-envelope exit refuses the complete field, including fully prescribed solutions. This is fixed-state constitutive storage, not a force, hysteresis-loss, demagnetization-history or live parameter-work result. Existing explicit frame/triangle sampling consumes the resulting B. Fixed-remanence winding response requires the explicit affine constructor described below. The stationary eddy-history owner continues to refuse permanent-magnet bias until its source/storage trajectory accounts for Br.

The recoil convention follows [COMSOL's remanent-flux constitutive relation](https://doc.comsol.com/6.3/doc/com.comsol.help.acdc/acdc_ug_magnetic_fields.08.032.html); the shifted fixed-state energy reference follows the linear-recoil discussion in [FEMM's permanent-magnet energy note](https://femm.info/doku/doku.php?id=pmenergy). `material_permanent` adds manufactured-source regressions for constant-Br natural boundaries, current-plus-remanence loading, magnet/air interface fields, direction reversal, material/temperature swaps, source/domain refusals and consumer guards. These are analytic model checks, not measured catalog-material validation.

## Fixed-remanence affine winding and circuit response

`PlanarWindingBasis::new_with_fixed_remanence` uses the same linear operator for a source-admitted zero-current field and explicit incremental unit-current responses. Its additional stationarity allowance is in Wb and checks `integral H0.B_unit dV`; the inductance-pairing allowance remains in H. Retained-basis capacity includes the bias field. The generic entry also accepts entirely induced-linear meshes as zero bias.

Total linkage is `lambda=lambda0+L*i`. Native `CircuitDae` coordinates remain `delta_lambda=L*i`, including declared end inductance; initialization and current recovery do not interpret permanent flux as current. Native magnetic/capacitor energy is relative to the zero-current state. `bias_energy_j` retains its constant field reference W0, and `stored_energy_with_bias_j` adds it for reporting without subtracting a large constant during step accounting. Changing fixed Br data changes bias flux/field/storage without inventing electrical power.

`field_for_currents` reconstructs `A0+sum(i_w*A_unit,w)` and uses original physical recovery to admit every recoil cell. It also superposes the original static current snapshots, allowing this actual field to use the existing instantaneous Maxwell-force owner; arbitrary induction recovery still lacks that provenance. Sampling, checked circuit initialization/`step_audited`, and Ohmic candidate endpoints share the recoil admission. Incremental responses are explicitly marked and carry no actual magnet operating-point claim. Envelope exits/cancellation return no accepted replacement. Endpoint checks do not bound the continuous trajectory between steps. The raw DAE exposes incremental coordinates and carries no material-domain policy itself; fixed-remanence callers use the checked wrapper.

An active conducting recoil triangle must reuse its exact original complete material bundle, including its electrical claim. This admits unchanged recoil/conductivity sources without generic permeability and refuses temperature, remanence, slope, limit or source refresh that would leave stale fixed bias. Ordinary induced-material electrical refresh retains its previous policy. Geometry, material conditions and remanence stay fixed; no motional EMF, hysteresis or thermal parameter work is supplied. Manufactured-card `winding_permanent` regressions cover affine linkage/storage and layered inductance oracles, native RL/Ohmic work, mutual branch mapping, data swaps, source-envelope refusal and fixed conducting-magnet source binding.

## Material-driven reversible nonlinear planar fields

`planar::nonlinear::PlanarNonlinearMagnetostatics` admits a linear or reversible
H(|B|) source per original triangle. It shares mesh/metric, exterior-edge and
component-gauge admission, shape gradients, weighted Whitney Hodge, incidence,
sparse assembly and checked CG with the existing planar owner. It assembles
the actual weak residual `integral H(curl A).curl(v)-J v`, and the symmetric
tangent with radial eigenvalue dH/dB and transverse eigenvalue H/B. The tensor
is rotated into the gradient frame before the existing degree-1 Hodge assembly.

The complete initial field must satisfy prescribed potentials and source
envelopes. Newton uses explicit iteration/linear/backtracking limits and an
energy-decreasing Armijo search; trial fields outside the original H(B) range
are rejected, never extrapolated. Acceptance recomputes the free-coordinate
nonlinear Euclidean residual [A m], independent of Krylov convergence. The
reported target is the maximum of the absolute allowance and the relative
allowance times max(initial residual, impressed-current load). Exhaustion or
cancellation returns no partial field. Material, geometry and source owners
are immutable throughout; traversal and tie-breaking use fixed storage order.

Returned B/H, per-triangle energy/coenergy, total storage and actual work counts
remain tied to the prepared owner. Nonlinear energy integrates the source H(B),
not B.H/2. `sample_cells` reuses the original explicit triangle/frame transfer,
so nonlinear material changes reach existing imposed-field flow consumers.
This is a finite-domain fixed-geometry P1 estimate, not a discretization error
certificate, nonlinear winding inductance, hysteresis, eddy-current trajectory,
force/coenergy motion derivative, thermal parameter-work or remanence model.
The same reference conventions are described in the
[FEMM manual](https://www.femm.info/doku/lib/exe/fetch.php?media=upload%3Afiles%3Amanual.pdf),
magnetic energy/coenergy discussion. Synthetic tests do not validate a catalog
material. Focused G1/G3/G4 tests cover exact discrete saturation, source-only
substitutions, oblique mixed-material interface continuity, energy derivatives,
the linear limit, source envelopes, bounded nonconvergence and cancellation.

## Spatial dielectric electrostatics and field-derived circuit capacitance

`electrostatic::PlanarElectrostatics` retains one complete material state per
original triangle. Its absolute SI permittivity is admitted by the same
`linear_electric_permittivity` law used by the uniform `ResolvedCapacitor`.
The existing weighted Whitney Hodge, incidence, sparse assembly and checked CG
owners solve `-div(epsilon grad(V)) = rho_free`; no magnetic material or fictitious
conductivity substitutes for a dielectric. Fixed external temperature, pressure
and composition remain attached to their original source queries. Frequency,
field, strain and stress-dependent or applicability-bounded laws require an
appropriate constitutive model and refuse here.

Exterior voltage traces cover complete mesh edges and anchor every connected
component. Other exterior edges carry zero normal displacement. Recovery reports
cell `E=-grad(V)`, `D=epsilon E` and `integral E.D/2`, without smoothing material
interfaces. Prescribed-node weak reactions `K*V-f` are electrode free charges;
their sum must cancel the prescribed volume charge. Free-node residuals and
`U - (V.volume_load + V.electrode_charge)/2` have independent caller-supplied
absolute charge/energy gates in addition to the checked linear residual.
These finite-element balance checks do not certify spatial approximation error.

`capacitance` solves two disjoint complete terminal traces at one and zero volts
with no volume charge. The source-bound unit-voltage field gives `C=2U`, checked
against its independently summed positive-terminal charge. Its `circuit_branch`
feeds the existing Kirchhoff/pHS capacitor with that geometry-derived value;
the original field object retains all source and geometry identities. Changing
only dielectric data or a supported fixed temperature changes both the spatial
field and native circuit storage/decay. Live coefficient changes still require
parameter-work accounting, and are not implemented by overwriting a branch.

The model is fixed-geometry planar, isotropic, linear, nondispersive and lossless.
Fringing within the specified domain is solved; missing exterior/end geometry,
leakage/conduction, dielectric loss, breakdown, electrostriction and evolving
thermal work are not inferred. `tests/electrostatic.rs` covers analytic layered
dielectrics, volume-charge/electrode Gauss balance, voltage-offset and energy
identities, native RC response, source-state refresh, and explicit refusals using
synthetic claims rather than measured material qualification.

## Multiple dielectric electrodes and native circuit storage

`PlanarElectrostatics::capacitance_matrix` extends the same source-resolved
field owner to a complete ordered set of disjoint electrode traces. The scalar
and multiple-electrode APIs share trace admission; every terminal covers whole
exterior edges, and unassigned exterior edges retain the natural insulated
boundary. Exactly one original checked field solve is performed per terminal,
with explicit terminal, total-solve and retained-field capacity bounds.
Every unit field retains the original mesh, material states and charge reactions.

The full Maxwell matrix is the symmetric electric-energy Gram matrix in farads.
Its coefficients must independently agree with every terminal charge reaction,
reciprocity and the common-voltage nullspace within the caller's coefficient
tolerance. No reference terminal is silently removed. An explicitly grounded
submatrix preserves the remaining original terminal order. A floating electrode
is not eliminated by deleting its row: it still needs its physical charge
constraint in the receiving model.

`response` forms voltage differences before superposition and integrates the
resulting original cell fields for stored energy. `circuit_at_voltages` retains
the source matrix and returns native Kirchhoff capacitor branches with
`C_pair=-C_ij`, together with consistent signed initial branch charges.
Its explicit one-to-one electrode/node map and i-before-j orientation determine
the exact native capacitor charge ordering. Actual terminal charge and network
storage must match the field at the supplied voltages within independent C/J
tolerances. Every terminal needs a capacitive connection. Zero mutual entries
produce no branch; positive off-diagonal entries refuse passive circuit lowering
without clipping, even though the original field matrix can remain useful.

The native circuit owner supplies subsequent KCL/KVL, charge dynamics, source
work and resistor losses. The extracted network is the dielectric store itself,
not an additional parallel energy account. Geometry and dielectric state remain
fixed; parameter work, dielectric loss and leakage remain outside this
reduction. The static field-force connection below is not a motion/work closure.
A balance gate is not a spatial-error certificate.
Focused synthetic tests compare a three-electrode mesh with independent P1
coefficients and direct field charge/energy, verify gauge and terminal ordering,
propagate a regional dielectric-data substitution through actual native coupled
RC dynamics, and refuse a positive-mutual coarse discretization without altering
its field matrix. Full native integration execution is still required in an
environment able to compile the pinned dependency graph.

### Charge-constrained floating electrodes

`electrostatic::network::floating` prepares mixed voltage/total-charge control
on the original Maxwell matrix. Each floating conductor remains equipotential
while its total free charge is held at the supplied value, including nonzero
trapped charge. The retained terminals keep original order; floating inputs use
their explicitly selected order. Topological gauge admission includes common
electrodes joining disconnected dielectric pieces and refuses an entirely
unreferenced component before numerical factorization.

The original `fs-la` Cholesky owner factors the diagonally equilibrated floating
block once. Subsequent queries reuse that factor and the retained-to-floating
voltage map, without new FE solves. Explicit terminal and dense-work allowances
bound preparation; cancellation is checked around the bounded dense factor and
through assembly and field recovery. The Schur complement retains the actual
floating path rather than deleting its matrix row. It is a differential
capacitance at fixed floating charge; nonzero charge also produces an affine
charge offset and stored energy, recovered by the mixed-control solve.

Every accepted solution recovers the original spatial field and independent
weak electrode reactions. Actual-amplitude total-charge and terminal/field
energy gates apply separately from reduction residuals. Voltage differences
remain separate from their common reference to preserve small physical fields.
The full electrode state initializes the existing native capacitor network and
its charge/energy history without discarding the floating nodes. Positive
mutual entries on a coarse field discretization do not prevent this SPD block
solve, although they still refuse pair-capacitor circuit lowering.

Focused `electrostatic::network::floating` regressions cover analytic neutral
and charged three-electrode fields, several simultaneous charge constraints,
capacitance scaling, common-potential and terminal permutation invariance,
native RC decay at conserved floating charge, physical resistor heat, gauge
refusal, independent actual-state gates and cancellation/retry. This is frozen
linear planar electrostatics: no changing geometry/material work, dielectric
dispersion, spatial-error certificate or empirical device validation is added.

## Spatial material DC current and native resistor heat

`conduction::PlanarConduction` binds each original planar triangle to the
canonical `ResolvedConductorState`, preserving its complete source bundle and
explicit independent resistivity/conductivity basis. The weighted Whitney Hodge
and incidence assemble `d0^T M(sigma) d0` in actual conductivity dimensions.
Dielectric and conductor fields share scalar Dirichlet elimination, checked CG,
and complete-electrode trace admission through one internal owner. They retain
their distinct physical laws: conductor reactions are current [A], and
`integral E.J dV` is irreversible power [W], with no dielectric `1/2` factor.

The stationary equation is `-div(sigma grad(V))=0`. Prescribed voltage traces
cover complete exterior edges and anchor each connected component; remaining
edges are insulated. Recovery keeps original cell `E=-grad(V)` and `J=sigma E`
without averaging across interfaces. Weak reactions are current entering the
conductor at prescribed nodes. Net terminal current, free-node current residual,
and the difference between integrated Joule and terminal power must satisfy
independent absolute physical allowances. These discrete balance gates do not
certify spatial approximation error.

`resistance` solves the actual heterogeneous geometry at one volt and checks
its integrated power against the independently summed terminal current before
returning a native `Branch::Resistor`. Its retained source/geometry field is the
same dissipative component as that branch. `PlanarResistance::field_at_voltage`
recovers current and per-original-triangle `cell_power_w` at the native circuit's
actual stage voltage without re-solving. For a midpoint circuit step, the
midpoint branch voltage and `dt*cell_power_w` partition that same native resistor
heat; the spatial report is not an additional loss. Complete material snapshots
and actual extruded volumes remain accessible for an explicit receiving thermal
cell map. This electrical owner neither advances thermal state nor invents a
heat partition between unrelated bodies.

Temperature, composition and fixed pressure act through original resolved
source queries. New supported material snapshots require a newly prepared
operator. Both independent electrical source bases use the common scalar
conductor law; no material-name branch or inferred unit-length wire enters.
Frequency-, field-, strain- or rate-bounded sources refuse the stationary scalar
model. The mesh and state remain fixed during each solve. Insulating regions
must be represented by the conductor domain's actual boundary, not artificial
small positive conductivities. Floating unknown-potential electrodes, nonlinear
conduction, displacement current, induction,
Hall effects and geometry/thermal evolution are outside this admission.

`PlanarConduction::new_with_contacts` adds source-backed zero-thickness Ohmic
interfaces to that SAME native stiffness. Each `ConductionContact` binds two
distinct original exterior edges with exactly coincident paired endpoints and
opposite adjacent cell interiors to a complete `ResolvedElectricalInterface`.
The ordered source material states must match those cells, including revision.
The source must support finite positive R'' at static external conditions;
pressure or material changes act through source data and a new preparation.
Explicit zero resistance requires shared mesh vertices rather than a penalty.

The contact law is `J_n = (V_a - V_b)/R''`. Its exact P1 trace matrix is
`A/(6 R'') * [[2,1],[1,2]]`, assembled with equal/opposite side blocks using
actual extruded edge area. Contacts connect component gauges; a voltage trace
can therefore anchor another body through the interface. Internal contact
edges cannot also be external voltage/current electrodes. Duplicate or
nonmatching interfaces refuse before publication.

`ConductionField::contact_response` retains endpoint voltage jumps/current
densities, signed integrated current and exact nonnegative surface Joule power.
`cell_power_w` remains bulk heat; `power_w` is bulk PLUS surface heat, checked
against the same independent terminal work. Native current drive, two-terminal
resistance and circuit-voltage reconstruction all use this operator and total
loss. Surface heat is not silently assigned to either neighboring bulk cell.
Moving-contact, dielectric charge evolution and the bulk-only enthalpy/circuit
adapters explicitly refuse these contacts until their surface work, charge and
heat-partition evolution is owned. Static field/circuit resistance extraction
is supported. `tests/conduction/contacts.rs` checks analytic series resistance,
pressure and bulk substitutions, current-controlled contact connectivity,
nonuniform surface quadrature, total Joule balance and unsupported interfaces.

`new_with_phase_contacts` retains that SAME operator and admits only explicitly
owned liquid-mass-fraction axes in addition to its original fixed-external
conditions. Each declared fraction is an existing dimensional scalar in [0,1];
aliases, missing/type-mismatched axes and unowned response/applicability refuse.
Axis rows follow their original ordered traces through canonical sorting and
remain inspectable on the returned problem. The unchanged `new_with_contacts`
entry does not acquire permission to freeze arbitrary phase or field axes.
The composing thermal owner supplies the exact phase from its original enthalpy
chart and refreshes pinned source snapshots. `fs-couple` now provides that
explicit spatial-contact/enthalpy composition, including two-sided surface heat;
the bulk-only adapters above still refuse contact losses they do not own.
No phase fraction is inferred from temperature, and no phase-interface law,
motion or automatic electrical mixing rule is invented by this constructor.

Existing mesh/assembly and Krylov budgets bound work; cancellation is polled
during material admission, assembly, elimination, field/heat recovery and branch
reconstruction. Refusal never modifies an accepted source or field. G1/G3/G4
regressions in `tests/conduction.rs` compare independent series and parallel
material geometries, source-basis and temperature changes, thickness, voltage
reversal and gauge shifts. Actual repeated native RC steps must lose precisely
the recovered spatial Joule heat. Unsupported source models, trace/assembly/
iteration limits and late cancellation refuse. Synthetic fixtures are not
experimental material validation. The continuum current and dissipation laws
follow MIT *Electromagnetic Fields and Energy*, sections
[7.2](https://web.mit.edu/6.013_book/www/chapter7/7.2.html) and
[11.3](https://web.mit.edu/6.013_book/www/chapter11/11.3.html).

## Electric field, circuit and native mechanical load

`PlanarCapacitanceMatrix::field_at_voltages` reconstructs an actual
`ElectricField` from its retained unit-electrode potentials. This invokes the
same native gradient, weak reaction and physical-energy recovery as a direct
solve, with independent actual-voltage charge and energy tolerances. It does
not perform another Krylov solve. Returned potentials use the first electrode
as their explicit zero reference; a common offset cannot create a field.
Circuit consumers supply their instantaneous endpoint voltages in the original
terminal ordering, preserving the exact dielectric mesh and source snapshots.

`ElectricField::maxwell_wrench` integrates `D tensor E - (D.E) I/2` through the
same contour, outward normals, extrusion thickness, frame/moment transform and
compensated reduction used by magnetic forces. It returns the existing
`fs_ga::Wrench` accepted by native articulated dynamics. Both incident cells of
every contour edge must have the same complete dielectric snapshot, zero
original volume-charge source, and no prescribed electrode vertex on the edge.
The field retains per-cell charge sources even when their global sum vanishes.
Material interfaces, mesh boundaries, folded neighbours and invalid budgets
refuse before any load can be published. Enclosed material contrast and free
charge are supported through their actual computed field.

The assumed homogeneous exterior has no electrostriction or mechanical fluid
pressure response. The result is a net instantaneous load on the enclosed body,
not a distributed traction on the artificial contour. Motion still requires
updated geometry and a source/field/mechanical work closure; no frozen-field
actuator or overlapping free-charge force may be counted as that closure.
Synthetic G1/G3/G4 tests compare against an independent volume-charge force,
check material-only changes through native acceleration, reconstruct the force
after a native multi-electrode RC discharge, and exercise gauge/frame,
source-domain, budget and cancellation behavior. These are discrete numerical
connections, not measured dielectric-device qualification. Force-density and
exterior-medium conventions follow the energy discussion in MIT 6.013,
chapter 11: https://web.mit.edu/6.013_book/www/chapter11/11.9.html.

## Nonlinear reversible induction and spatial Ohmic heat

`planar::eddy::nonlinear::PlanarNonlinearEddy` retains the actual source-bound
P1 field, electrical regions, fixed boundary trace and physical clock. Its
backward-Euler endpoint equation is `F_H(A_next) + C_sigma delta A/dt = load(J_next)`.
The same material H(B) response, FEEC tangent, bounded Newton/Armijo search and
checked CG used by the nonlinear static owner solve this positive variational
problem. The induction term contributes its quadratic objective, residual and
tangent together. No second magnetic solver or inferred permeability is used.

Linear midpoint and nonlinear backward-Euler histories share conductor source
admission, connected-region topology, isolated-current projection, complete
electrical roles, balanced impressed-current loads and algebraic constraints.
They also share the exact consistent-P1 `sigma E_z^2` integration. Every output
`EddyCellHeat` retains its original source, triangle, nodal E/J, volume, power
and heat. A conductivity refresh pins the complete original card, requirements,
coordinate schema, source selection and electrical basis. Nonlinear conductors
also retain the H(B) claim. A temperature change is admitted only when the full
original H(B) range is still supported and its constitutive law is unchanged;
the explicit linear-temperature refresh described below owns changes in mu(T).

For a fixed law, the source-integrated endpoint energies satisfy
`delta psi + Joule heat + numerical dissipation = endpoint impressed work`.
The separately reported numerical term is the nonnegative magnetic Bregman
loss `V*(H_next dot (B_next-B_old) - (psi_next-psi_old))`. It uses the original
positive source slopes and angular convexity to avoid cancellation; the linear
limit is `V*|B_next-B_old|^2/(2*mu)`. This is time-discretization loss and MUST NOT
heat receiving material. Only the existing spatial Joule account is deposited.
Per-cell and total numerical losses, actual source work and the full nonlinear
residual remain visible. Failure of source support, nonlinear budget, current
balance, energy accounting or cancellation returns no candidate replacement.

This is a first-order stationary, isotropic, reversible magnetoquasistatic
section with fixed geometry. It does not add hysteresis, remanence evolution or
conductor motion. A recovered transient field uses the existing explicit
position/frame transfer without re-solving static equilibrium. Synthetic
`nonlinear_eddy` tests cover independent P1 decay, material-data substitution,
Joule/numerical separation, time refinement, conductivity refresh, balanced
excitation through conducting/insulating cells, source refusal and rollback.
The std-only production conduction leaf also has direct analytical projection,
power, floating-current and gauge/reversal tests; those leaf results do not
establish execution of the complete material/mesh/Newton integration.

## Circuit-driven nonlinear material windings

`winding::nonlinear::NonlinearWindingCircuit` couples actual stranded terminal
currents to `PlanarNonlinearEddy` and the original `fs-phs` Kirchhoff descriptor.
The authored signed `J_z/I` profiles reuse the original winding component
balance and exact P1 linkage pairing. A bound original inductor branch becomes
a private native power port. The default current-input causality has input `-I`
and output `V_a-V_b`.
Its placeholder inductance contributes no storage. Unbound inductors,
capacitors, transformers, graph resistors and external ports retain their
original circuit owners and coordinate ordering.

The existing bounded Newton--Krylov and fallible directional-difference owners
solve `V = delta(lambda)/dt + R_winding*I`, with endpoint winding currents and
`lambda = integral(A_z*J_unit dV) + L_end*I`. Every trial restarts both accepted
histories. Declared nonnegative end/leakage self-inductances are optional;
unrepresented end geometry, resistance and mutual coupling are not invented.
Disjoint stranded triangles use the same source/resistance and spatial Joule
owners as linear windings. They cannot also own resolved eddy conductivity.
Complete material/card/claim/basis continuity is retained across conductivity
temperature refresh. The original full nonlinear H(B) law remains frozen;
linear mu(T) may refresh through the explicit material-work path below.

The private port work `-dt*I*V` cancels the field/winding terminal work.
`NonlinearWindingStep.supplied_work_j` contains only original external source
work; its native descriptor receipt still visibly includes private ports.
The combined audit checks field, native circuit and end-inductance storage;
spatial core and stranded Joule heat; explicit graph heat; and separately
reported magnetic/end-inductance backward-Euler loss. Numerical loss is never
physical heating. Terminal voltage residual, source/linkage work pairing and
each heat decomposition must pass physical allowances before publication.
The candidate publishes the circuit and field together, preserving original
histories on source, budget, native constraint, nonlinear or cancellation refusal.

`NonlinearWindingCircuit::new_with_ports` also admits explicit voltage-input
causality for current-driven windings, including ideal current sources. That
private native voltage-source stamp receives the actual
`delta(lambda)/dt + R_winding*I` and returns `-I_native`; the same Newton owner
solves `I_native-I_field = 0`. Current-input and voltage-input windings can share
one nonlinear field, including signed mutual linkage. A declared positive
`residual_scale_ohm` only converts the current residual to the shared numerical
voltage scale. It is not a physical resistor and adds no heat, voltage drop,
stored energy or parasitic regularization. The independent declared
`current_tolerance_a` checks the actual current mismatch; the existing voltage
and combined energy allowances remain mandatory. Step receipts expose both
physical voltage and current residuals, separate from scaled Newton residuals.

Voltage-input initialization requires an explicit initial terminal voltage,
compatible original native fluxes/charges and matching admitted field current.
The original rank-aware native initializer owns that consistency check. Its
zero-time algebraic midpoint receipt determines the current; ambiguous initial
current splits may refuse instead of silently replacing the field history.
Finite voltage-input trials retain the accepted native state and let the
original descriptor step evolve its storage. They do not reinitialize or
project capacitor charge against a trial step-average voltage. Source-bound
graph-resistance refresh retains the same terminal causality and port ordering.
Original external inputs are held over each finite interval, and endpoint field
currents pair with native step-average voltages in this first-order method.
The receipt describes finite-step work; it is not an instantaneous-source
impulse model or an automatic resolution of incompatible ideal-source constraints.

This first-order coupling admits at most 16 windings within the existing
32-node/64-branch/32-port/128-state allocating circuit bounds. The native graph
must admit its explicitly selected private-port trials; singular ideal-source
topologies can refuse this partitioned lane. Geometry, proposed magnetic law and
native graph coefficients stay fixed during each native trial. No hysteresis
or experimental material-qualification claim is introduced.

`nonlinear_winding_circuit::current` supplies independent manufactured oracles
for signed ideal current drive through linear/saturating fields, declared end
storage and continuation; smooth-ramp work and numerical-loss time refinement;
native capacitor/inductor history; mixed causal and reversed mutual windings
with source-bound graph resistance; core eddy shielding
and heat; and initialization/current-budget/cancellation refusals. These cases
exercise the actual field and native circuit owners, not a replacement lumped
inductance model; their material laws do not claim measured qualification.

`NonlinearWindingCircuit::step_with_eulerian_velocity_candidate` resolves
moving bulk conductors inside this same terminal Newton solve. Its actual
nonlinear moving-field candidate changes the winding flux linkage and solved
load current; no fitted speed voltage or independent force current is added.
The public controls, step and candidate retain stationary default field types
and carry the complete moving-field step when this operation is selected.
The combined balance is `delta(storage) + Q + D_BE + W_mechanical = W_external`;
signed mechanical work is positive out of the field, while stranded, graph and
bulk conductor-frame Joule heat remain separate nonnegative physical losses.
Every trial preserves the accepted circuit/field history, and terminal voltage,
source/linkage reciprocity and the combined energy account are checked before
publication. Complete finite P1 velocity is required. Stranded winding cells
and fixed-recoil magnet cells must be stationary; source laws, conductivity,
mesh, material occupancy and profiles stay fixed. The underlying moving-field
divergence, connected-current, source-envelope, magnetic Peclet and nonlinear
residual bounds apply. Moving winding geometry, transported magnetization and
sliding interfaces remain outside this operation.

### Oriented permanent magnets with saturating cores and induction history

`IsotropicMagneticMaterial::Recoil(PlanarRecoilMaterial)` admits the same
source-backed, specimen-oriented recoil law as the original linear owner.
Linear regions, single-valued reversible H(B) cores and permanent magnets share
the original nonlinear P1/FEEC solve. Magnet cells evaluate
`H=(B-Br)/mu_recoil`, shifted storage `psi=|B-Br|^2/(2 mu_recoil)`, and
Legendre coenergy `psi*=mu_recoil |H|^2/2+Br.H`. This coenergy can be negative.
Remanence enters the constitutive weak residual independently of impressed
current, including on natural boundaries; it is not a fitted current source.

Each initial and trial field must satisfy the original parallel and transverse
recoil envelope. Newton backtracking handles both nonlinear H(B) range and
recoil-range refusals. Accepted fields expose per-cell operating points and
remaining source margins beside their declared energy reference. These are
numerical estimates, not certified bounds or an irreversible-demagnetization
model. Construction checks the positive recoil tangent without requiring an
inadmissible zero-field initial state; callers can supply a valid biased seed.

The same fixed law participates in `PlanarNonlinearEddy` and
`NonlinearWindingCircuit`. Physical storage uses the shifted energy, while
the recoil backward-Euler Bregman loss is computed directly as
`|B_next-B_old|^2/(2 mu_recoil)`. Constant Br cancels from this increment.
Winding linkage retains the actual magnetic bias; its time difference supplies
Faraday voltage. Spatial Ohmic loss remains the only magnetic material heat
in this lane. Accepted induction fields can transfer an exterior Maxwell
wrench while retaining enclosed magnets and saturating material.

Conductivity refresh requires the original complete material source and
unchanged recoil quantities. Magnetic-law refresh preserves the exact recoil
law and orientation. Temperature-dependent remanence, permeability or envelope
changes are refused; the caloric nonlinear-recoil path is also refused until
remanence work and entropy data have an owner. Geometry and remanence remain
fixed. No hysteresis, irreversible loss, moving-magnet induction or 3D end-field
claim follows from this extension.

### Source-pinned magnetic temperature refresh on continuing circuits

`PlanarNonlinearEddy::step_with_materials_candidate` and
`NonlinearWindingCircuit::step_with_materials_candidate` reassemble the same
original mesh at proposed source-pinned linear mu(T) states. Original material
cards, complete claims, coordinate schema, fixed context and constitutive kinds
remain pinned. Nonlinear H(B) and oriented recoil entries must remain exactly
unchanged. A refreshed
candidate owns its assembled problem and retains original field, circuit and
linkage histories through clone and continuation.

Each magnetic cell reports signed `parameter_work_j` equal to
`P=psi_next(B_old)-psi_old(B_old)`. The backward-Euler Bregman loss uses the
endpoint law for both field values. Native energy closure is consequently
`delta psi + Joule heat + numerical loss = external work + P`, including native
circuit and declared end storage for windings. Source work and parameter work
remain distinct. A caloric owner must separately resolve field entropy and
reversible heat; neither P nor numerical loss is automatically material heating.
Original fixed-law stepping uses the same implementation with zero P.

`nonlinear_winding_circuit::thermal` supplies independent RL and conducting-P1
checks for data-only mu(T), signed local/global work, physical/numerical heat
separation, owned restart, source correspondence and late cancellation. Native
execution must be reported separately from parser and source review.

`NonlinearWindingCircuit::with_graph_resistors` binds every ordinary graph
resistor to a complete `GraphResistorMaterial`. It reassembles only the original
native resistor coefficients from the common `ResolvedConductor` owner. Original
branch endpoints, transformers, private winding ports, state-coordinate maps,
charges, fluxes, field history and winding linkage stay intact. Exact native
storage equality is admitted before returning an independent candidate at the
same clock. A static scalar Ohmic source and complete unique original-resistor
coverage are required; winding branches cannot be replaced through this API.
The full source bundles and geometry remain with the accepted circuit and its
step report beside the original graph-local heat records. Explicit resistance
reconfiguration stores no parameter energy. Thermal consumers must admit the
accepted sources and refresh them at their actual endpoint through the shared
material/enthalpy owner; this API supplies no temperature or contact law.
Synthetic `nonlinear_winding_circuit` regressions exercise independent nonlinear
RL and coupled-winding/capacitor equations, source-only H(B)/resistivity swaps,
core/strand loss separation, temperature refresh and transactional failure.

## Accepted induction and nonlinear circuit fields drive mechanical loads

`PlanarEddy::maxwell_wrench` and `PlanarNonlinearEddy::maxwell_wrench` recover
the actual accepted B/H from their retained potential and share one source
admission with the original geometric Maxwell-stress integral. They do not
re-equilibrate the field. Every contour edge needs the same complete linear
material snapshot on both sides, zero original impressed current, and no
resolved eddy conductor. Enclosed source-defined nonlinear magnetic material, fixed recoil magnets,
stranded current and eddy currents remain admissible. The history supplies the
source support that a detached recovered static field intentionally lacks.
Mesh work bounds apply before field recovery. Original edge orientation,
extrusion, exterior sampling, compensated reduction, mechanical frame and
moment origin remain owned by the common contour implementation.

`NonlinearWindingCircuit::field()` exposes this same accepted induction state,
so source saturation/resistivity and actual voltage-driven current determine
its native mechanical wrench. The result can feed the existing articulated
body dynamics without a fitted force coefficient or material-name branch.
This is an instantaneous net magnetoquasistatic load on a fixed planar section,
not a local nodal traction map or a moving-field work closure. No energy is
debited by this observation, and an overlapping Lorentz body load must not be
added. Nonlinear/conducting exterior stress, displacement-current electric
stress, end fields, striction and exterior-medium pressure response remain
outside this admission.

G1/G3/G4 tests in `tests/magnetic_force/induction.rs` cover the linear history
limit and independent volume-Lorentz refinement, source-only saturation and
voltage changes through the native winding circuit and body acceleration,
retained decaying-eddy force versus an incorrect static reset, source/contour
refusals and late cancellation. Full native integration execution remains
required; source parsing and metadata compilation are not physical validation.

## Distributed material-conductor current boundaries

`PlanarConduction::solve_with_current_traces` accepts signed uniform inward
current density on actual exterior edges beside prescribed voltage traces.
The original planar metric supplies extruded edge area; exact P1 integration
supplies half the edge current to each endpoint of the existing FEEC solve.
Current traces are not assumed equipotential. Duplicate/interior edges,
nonfinite loads and simultaneous voltage/current traces on the same edge refuse;
every connected component retains a real voltage trace and return-current path.

Free-node current residuals subtract the authored load. Terminal power uses
the prescribed current at current-only vertices and the actual combined weak
reaction at voltage vertices. `nodal_applied_current_a` separately preserves
the applied trace contribution at mixed corners, allowing the voltage-source
reaction to be recovered without double counting. The original cell E.J and
terminal work must agree before acceptance. Zero-load voltage and native
two-terminal reductions retain their existing owner. G1/G3/G4 regressions check
layered material fields, extrusion, reversal, gauge, mixed corners, input order
and transactional cancellation. This is fixed-geometry isotropic DC conduction,
not an ideal equipotential current electrode, contact law or electrical storage.

## Oriented material conductivity in the native spatial solve

`PlanarConduction<C>`, its field and its resistance reduction accept the shared
source-bearing `PlanarConductorMaterial`. The default remains the original
scalar conductor. Oriented principal data and mixed scalar/tensor cells use the
same FEEC assembly, voltage/current boundaries, solve and power acceptance.
The existing weighted owner admits tensor positivity; current is `J=sigma E`
and the original per-cell integral of `E.J` must match terminal work. Tensor
fields expose their full matrix instead of an effective scalar surrogate.
Source-fixed frames must name the same geometry coordinate frame throughout
one problem. Fixed geometry, symmetric DC and static external-condition limits
remain; no Hall response or out-of-plane closure is inferred. Four focused
tests cover independent rotated fields, current traces, heterogeneous resistance,
isotropic degeneration, missing principal sources and inconsistent frames.

The same static field owner also admits an explicitly retained liquid-mass-
fraction coordinate from the shared conductor material type. Its physical
owner must supply that actual coordinate; no arbitrary source axis becomes
fixed merely by appearing in a query. All other source conditions retain the
existing static admission. The field remains a stationary snapshot; latent
enthalpy evolution, fixed-volume phase admission and circuit feedback belong
to the existing `fs-couple` thermal owner. No phase geometry, phase-interface
resistance or automatic electrical mixture law is inferred here.

## Native spatial leaky-dielectric charge and energy history

### Material replacement at retained native charge

`PlanarLeakyDielectric::rebind_material_candidate` substitutes new complete
dielectric/conductor source bundles on the original mesh at the accepted clock.
The native scalar-potential solver recomputes voltage at conserved free-node
charge; those charge coordinates are retained bit-for-bit. Original electrode
voltages remain held. Their charge transfer and signed electrical work are
explicit, separately from signed material-parameter work. The shared native
coefficient-work calculation includes all tensor cross terms and independently
closes `delta_U = W_electrodes + W_material`. Neither work account is Joule heat.

New material cards, electrical basis and material orientations are data; declared
mesh frames must agree. New dielectric laws are validity-wide constant and
subsequent ordinary temperature/phase updates remain pinned to the replacement
sources. Conductivity-only changes retain voltage and charge exactly and alter
subsequent native current/heat. Work bounds, assembly budgets, field/charge
gates and cancellation refuse transactionally. This explicit instantaneous
constitutive replacement does not predict an insertion trajectory, material
formation energy, thermal entropy, mass transport or radiative transients.

G1/G3/G4 tests in `electrostatic/rebind.rs` cover conserved-charge and held-voltage
capacitor limits, heterogeneous redistribution and reversal, rotated dielectric
data, electrical-basis substitution through continued relaxation, gauge
invariance, work refusal and late-cancellation retry. These are synthetic
analytical model checks, not measured material validation.

### Charge relaxation

`electrostatic::relaxation::PlanarLeakyDielectric<C>` retains physical P1
potentials, weak nodal charges, cell dielectric storage and time on the original
planar mesh. It uses the existing dielectric and scalar/oriented Ohmic Hodge
operators and checked potential solve. The midpoint interval satisfies
`(K_epsilon + dt K_sigma/2) V1 = (K_epsilon - dt K_sigma/2) V0`
at free vertices. No lumped fitted capacitance or duplicate cell constitutive
law replaces the native operators.

Original complete exterior voltage traces remain fixed in topology and may
ramp in value. Other boundaries isolate total conduction plus displacement
current. Electrode reactions include both terms and pair with midpoint voltage
for source work. Shared native cell recovery computes physical local E.J heat,
E.D/2 storage and displacement current. The candidate checks free-node current,
net electrode current, weak-charge/cell-storage agreement and
`delta U_field + Q_Joule - W_electrodes` before publication.

The dielectric source must declare validity-wide constant absolute scalar
permittivity. Conductivity may refresh at the actual state, retaining original
material identity, full source claims, basis, frame, phase axis and coordinate
meanings. The original accepted charge is the starting point of every trial;
failure or cancellation publishes no state. Variable permittivity, dispersion,
magnetic induction, breakdown, carrier chemistry and moving interfaces are not
represented.

Held-coefficient midpoint accuracy is second order, but charge sign is not
preserved for arbitrary intervals; the physical relaxation time and spatial
field still require independent resolution. Native regression cases live in
the fs-couple `conduction_phase_enthalpy::dielectric` lane and compare analytic
epsilon/sigma decay, heterogeneous interface charge and voltage-ramp work.
Source/parser review completed for this change. These regression cases have not
been executed in the current environment because pinned sibling dependencies
and configured remote runners are unavailable.


### Directional dielectric storage in native fields and charge relaxation

`PlanarElectrostatics::new_with_laws` assembles scalar or explicitly oriented
per-cell source permittivities through the original weighted FEEC Hodge
operator. Field recovery uses `D = epsilon E` and storage uses the physical
integral of `E.D/2`, including cross terms. The native multi-electrode
Gram/reaction/circuit extraction and superposed response retain that same
tensor storage. Existing pair-capacitor reduction still refuses nonpassive
pair coefficients. No scalar average or conductivity-derived frame is used.
Scalar construction remains available. The scalar accessor
`permittivity_f_m()` now returns `Option<&[f64]>`, with `None` for any
tensor/mixed field; the complete coefficient and source-law views are explicit.

`PlanarLeakyDielectric::new_with_dielectric_laws` retains these laws through
the original continuing charge owner and complete conductivity/temperature/
phase material refresh. Dielectric and conductor material axes may differ;
their geometry frames must identify the same physical mesh. Every selected
principal permittivity must declare validity-wide constancy, and every
rotated matrix entry must remain unchanged during a trajectory. The original
charge, current, storage, source-work and rollback gates remain active.

Maxwell force contours retain their explicitly isotropic exterior convention.
An anisotropic interior may be enclosed by a supported scalar collar.
Declared tensor exterior laws refuse even when their numerical eigenvalues
coincide. This does not supply universal anisotropic mechanical stress,
normal-direction closure, dispersion, induction or moving geometry.

The focused tests specify oblique displacement/storage, specimen/frame
rotation, native capacitance extraction and independent RC decay, and
anisotropic interiors with scalar force collars. Reconstructed runtime
sources match the previously parsed source lengths; final reconstructed
sources/tests received source/API review but could not be parsed, compiled
or executed after the executor became unavailable.


### Retained charge under source-backed epsilon(T)

`PlanarLeakyDielectric::new_with_thermal_dielectric_laws` retains the original
material derivative laws in the SAME native spatial charge owner. Its step uses
`(Kepsilon1+dt Ksigma/2)V1=q0-dt Ksigma V0/2`: q0 is the original accepted
charge, never reconstructed with replacement epsilon. Original conductivity
recovery supplies midpoint E.J; the original dielectric recovery supplies both
endpoint D and Helmholtz storage. Displacement current uses actual D1-D0.

Each original triangle reports the exact two-endpoint coefficient work
`P=-V E0.(epsilon1-epsilon0).E1/2`, closing
`delta_psi+Q_Joule-W_electrode-P=0`. This signed transfer is separate from
Ohmic loss. A thermal owner must also close field entropy/internal energy;
merely applying epsilon(T) and relabeling the storage residual as heat is not
admitted. The ordinary constant-epsilon construction and phase-refresh gates
remain intact. Fixed mesh, voltage traces, independent conductivity/dielectric
frames, original source bundle and selected claims retain their native checks.

Focused coupled regressions cover scalar and oriented source laws, both slope
signs and zero slope, continued charge/heat, electrode work and the isolated-charge
limit. Changed native source passed Rust syntax and bounded independent
accounting/source review; native compilation/test execution was unavailable.

### Saturating and fixed-recoil moving induction

`PlanarNonlinearEddy::step_with_eulerian_velocity_candidate` advances
source-bound reversible H(B) and stationary recoil bias with prescribed P1
conductor velocity. Backward Euler uses the original analytic FEEC tangent plus
the nonsymmetric motion derivative, bounded Newton/FGMRES and residual
backtracking. The linear and nonlinear paths share their matrix-free motion
action/transpose, floating-conductor projection and exact P1 current/heat/work
recovery; Arnoldi actions reuse full-mesh scratch.

One projected endpoint current supplies conductor-frame Joule heat and signed
Lorentz work, closing delta psi + Q + W_mechanical + D_BE = W_impressed. The
original convex Bregman D_BE remains numerical loss, never a thermal credit.
Complete-source slope bounds govern the caller-admitted magnetic Peclet limit.
Recoil cells require zero velocity, whether conducting or insulating. Stationary
magnets may bias adjacent moving metal, but permanent magnetization is not
transported. Geometry, material sources, remanence axes and electrical regions
remain fixed; every failed candidate preserves accepted history.

The scheme is first order in time. The centered P1 Peclet allowance is a
formulation admission, not an error certificate. Moving interfaces, thermal
parameter evolution within the operation, hysteresis, irreversible
demagnetization and terminal-circuit motion coupling require separate owners.
Five physical regressions cover an independently assembled nonlinear P1
residual and heat/work quadrature, an analytic stationary-magnet/moving-metal
transient with time refinement, stationary continuation, saturated uniform-field
charge cancellation and source/domain/work/cancellation refusals.

### Anisotropic moving-conductor charge, heat and force

`PlanarConduction::solve_with_motion` reuses the original source-resolved
scalar/rotated SPD conductivity Hodge and checked charge solve for per-original-cell
planar velocity and an imposed normal magnetic field. The weak motional load,
total-current boundary reactions, `J = Sigma(E + u cross B)`, conductor-frame
Joule heat and `J cross B` force share one current/material/volume owner.
`MovingConductionField` retains the authored motion and field snapshots and
separately reports heat, signed mechanical work, volume `E dot J` and independently
recovered terminal work. Admission checks current continuity, the sum of absolute
local work-pair defects, and `heat + mechanical = terminal power`. Unloaded
exterior boundaries insulate the total current. Velocity or field changes reuse
the prepared conductivity matrix. Original stationary and zero-motional-field
recovery remain intact.

Scope is a prescribed fixed-geometry quasistatic snapshot; no self-consistent
induction, Hall law, moving mesh or fluid/solid time integrator is inferred.
Six analytical tests cover open-circuit charge separation without spurious drag,
shorted braking, loaded generation and motor signs, rotated anisotropic currents,
imposed-current EMF, heterogeneous counter-moving layers and stationary/refusal
behavior. Synthetic sources verify these equations, not an empirical material fit.

### Three-dimensional current-free magnetic fields

`spatial::SpatialMagnetostatics` solves `H = -grad(phi)`,
`B = mu H + Br` and weak `div(B) = 0` on actual tetrahedra. Each original
cell retains its resolved positive permeability or source-backed oriented recoil
law. The original Whitney 1-form material Hodge, incidence/stiffness assembly
and checked CG own the numerical operator. Prepared geometry and stiffness are
reused across boundary solves. `try_element_geometry` supplies a checked
companion to the original tetrahedral geometry owner and refuses malformed,
singular or unrepresentable cells before field publication.

Explicit exterior faces accept continuous nodal magnetic potential [A] or
outward normal B [T]; unspecified faces impose zero normal flux. Each disconnected
pure-flux component must satisfy its own integrated flux compatibility and
receives one deterministic gauge. The independent full weak residual includes
the artificial gauge node, so the gauge cannot supply missing flux. Relative
potentials preserve field and energy under representable common gauge shifts.

Recovery returns actual cell H/B, original recoil operating-point admission,
weak boundary fluxes and interior normal-flux jumps. Recovered cell energy is
checked against the assembled Hodge energy. Fixed remanence uses the existing
shifted reference `mu |H|^2 / 2`; this is not an independent force or variable-
magnetization energy law. Failed domain, solve, work or cancellation checks
return no partial field.

This is a single-valued scalar-potential model on current-free domains with
fixed linear/recoil sources and explicit finite boundaries. Linked free-current
circulation, general 3D current-driven curl-curl fields, induction, saturation,
hysteresis, motion and infinite-domain accuracy require other owners. The mesh
must be a conforming nonoverlapping tetrahedral complex; local adjacency checks
do not certify absence of remote cell intersections. Raw P1 B is cellwise
constant and its interface jumps are reported, not claimed to be an exactly
equilibrated H(div) flux.

Seven analytical tests cover vector fields with potential/flux boundaries,
large gauge offsets, oriented recoil/air demagnetizing layers and energy,
harmonic-field mesh convergence, disconnected compatibility and malformed
geometry/boundaries, original recoil limits, true solver-budget exhaustion and
late cancellation. These are synthetic physical oracles, not measured device
validation or a general 3D electromagnetic capability claim.

### Saturating three-dimensional current-free fields

`spatial::nonlinear` uses the original current-free H = -grad(phi) weak flux
problem and source-resolved linear, reversible H(B), or oriented fixed-recoil
materials on the same tetrahedral mesh. It inverts the original monotone H(B)
segments without fitting a second law. The exact radial differential
permeability supplies a tensor Whitney tangent; bounded Newton/checked-CG
steps use decreasing integrated coenergy with source-range-preserving
backtracking. Geometry, scalar boundary admission, component gauges and the
reduced checked-CG owner are shared with the original linear spatial solver.

Results retain distinct source-integrated energy and coenergy, a Legendre
identity audit, full weak flux residuals (including gauge vertices), raw
interface flux jumps and recoil operating points. Physical boundary work uses
unknown reactions only on physically prescribed potential traces; elsewhere
it uses the authored normal flux. A loose residual allowance cannot turn an
interior defect into boundary supply. Invalid source ranges, unbalanced
Neumann flux, exhaustion and cancellation return no accepted field.

The focused spatial_nonlinear target covers reduction to the original linear
solver, rotated saturating/air layers, recoil/core/air flux and exact energy,
a refined nonlinear radial manufactured field, differential material response,
and source, energy, budget and cancellation refusals. This is a fixed-state
single-valued scalar-potential model; currents, hysteresis, eddy currents,
moving work and thermal evolution require their respective owners. P1 flux
jumps and finite exterior truncation remain discretization errors; these tests
do not establish an experimental material calibration or a performance claim.

### Three-dimensional Maxwell force and torque

`SpatialMagneticField` and `SpatialNonlinearField` expose `maxwell_wrench`
through one private 3D surface integrator. Original tetrahedron selection
defines a closed outward surface. Exact triangle area/centroid quadrature
integrates outer-cell Maxwell stress and torque. Every contour neighbor must
share one complete source-identical induced-linear exterior snapshot;
nonlinear/recoil interfaces and mesh boundaries refuse. Saturated and remanent
targets may be enclosed inside a linear collar. Canonical Se3 placement and an
explicit receiving-frame torque origin return the native Wrench with a
canonical compensated reduction, work budgets and cancellation.

The five focused force tests cover independent nonzero cuboid refinement on
two surfaces, saturated-core symmetry and interface refusal, frame/origin
behavior, and source/region admission. The cuboid result has first-order raw
field/surface convergence. Test-only extrapolation checks its independent
continuum limit; the public API returns the raw discrete stress result, with
no extrapolation or inferred error bound. Field and surface refinement are
required for useful force accuracy. This is fixed-geometry static loading;
moving-geometry work exchange, material striction, equilibrated stress and a
virtual-work cross-check require separate physical work.

### Current-carrying three-dimensional vector-potential magnetostatics

`spatial::vector_potential` prepares the actual lowest-order 3D Whitney curl-curl operator d1^T M_(nu,2) d1 for source-bound isotropic linear and oriented fixed-recoil media. It admits exact piecewise-constant impressed-current loads, exterior tangential-potential edge traces, and natural n cross H=0 elsewhere. Current continuity and the entire relative curl-kernel compatibility are checked, including harmonic modes on handled domains, without silently changing the source. The existing checked CG solves the explicitly compatible invariant ranges of semidefinite operators; neither the full curl-curl operator nor the auxiliary face operator is claimed globally SPD. The bounded original cohomology-metric projection removes all curl-free gauge components and must preserve every magnetic face flux. Accepted fields retain truthful free-edge equations, divergence-free face flux, reconstructed B/H, original recoil-envelope admission, shifted storage, and independent current/trace constitutive pairing. The shared Maxwell integrator additionally requires one source-identical induced-linear current-free exterior collar. This is static lowest-order magnetostatics; higher-order spaces, Coulomb saddle/HX solvers, 3D eddy or nonlinear vector evolution, and TEAM qualification remain open.

The six focused spatial_vector_potential tests cover an independent
current-driven quadratic-potential manufactured solution under refinement,
exact uniform vector fields with gradient-shifted traces, a handled-domain
harmonic source that passes local continuity but must be refused globally,
recoil energy and limits, actual budgets/cancellation, and the current-free
Maxwell collar. The stored static source/trace pairings are not time-integrated
work or Ohmic heat.

### Current-driven reversible saturation in three dimensions

`spatial::vector_potential::nonlinear` extends the original Whitney edge
magnetostatics to source-bound reversible H(B), linear and oriented fixed-recoil
cells. Exact original source energy supplies H and its full symmetric radial
reluctivity tangent; only the changing material Hodge is rebuilt. Geometry,
current loading, boundary admission and complete relative-kernel compatibility
remain shared with the linear vector owner. Newton minimizes integrated source
energy minus impressed-current pairing at fixed tangential A, with bounded
source-preserving backtracking and the existing checked CG/coexact owners.

Newton's assembled dual residual is projected into the admitted coexact tangent
space before its relative linear solve. The measured compatibility defect must
stay inside the original physical current allowance. Impressed current is never
changed: every accepted field must satisfy the original full free-edge equations,
source Legendre identity, div-curl/reconstructed face flux and independent
current/physical-trace energy pairing. Energy and coenergy remain distinct.
Maxwell loads use the existing current-free, source-identical linear collar.

The six focused nonlinear-vector tests cover reduction to the original linear
current solver, gradient-shifted traces, first-order saturated current-driven
manufactured-field convergence, rotated source energy derivatives, recoil/core/
air tangential-field and dual-energy agreement, actual source/iteration/work/
cancellation refusals, and force admission. The six original vector regressions
also pass after the shared extraction. This is fixed-state static reversible
magnetism; hysteresis, nonlinear induction, thermomagnetic evolution and
validated engineering force bounds are separate extensions.


## Prepared distributed Maxwell force and torque

`spatial::force::eggshell::PreparedSpatialEggshell` caches an original-mesh
P1 virtual-displacement functional. Nodal weights are one over a complete
target plateau and zero on every exterior mesh vertex. All changing-weight
cells and intermediate constant plateaus require one identical original
induced-linear, current-free exterior material when the functional is applied.
Source/current-carrying targets cannot hide inside fractional plateaus.
Preparation binds the actual immutable mesh/positions, not merely matching
dimensions or topology, and is reusable across accepted excitations.

Translation integrates `-T grad(w)` using original cell B/H and volumes.
Rotation uses the original P1 interpolant
`sum_i N_i w_i (omega cross (x_i - origin))`, including its complete gradient
moment. It is exactly rigid over the target plateau and vanishes outside.
Both functionals integrate the discrete constant cell fields exactly without
smoothing interfaces, inferring a stress inside magnetic matter, or requiring
another field solve. Deterministic compensated reduction returns the canonical
`fs_ga::Wrench`, transformed through the caller's SE(3) frame and moment origin.

Scalar linear/recoil, scalar reversible nonlinear, current-driven linear/recoil
vector, and current-driven nonlinear vector owners share this same sealed
geometry/stress implementation. Original source identity, current-free exterior,
field provenance, geometry budgets and cancellation gate every application.
This is instantaneous mechanical load transfer. It does not supply a moving
domain energy exchange, material striction, a source/control shape adjoint,
equilibrated stress or a verified physical error bound. Collar-profile and
mesh refinement remain necessary.

The six added cases in `spatial_maxwell_force/eggshell` check raw cuboid
force/torque against analytic targets and two collar profiles, independent
re-solved coenergy translation finite differences, reuse/frame/origin/geometry
identity, complete closed weights and source/fractional-plateau refusal,
a saturated enclosed core, and vector/current-collar admission. The original
surface-force regressions remain in the same target.


## Moving-conductor 3D A–phi induction on the native Whitney owner

`spatial::vector_potential::induction` continues an admitted
`SpatialVectorField` with fixed tangential-A traces and original linear or
fixed-recoil magnetic materials. Prepared sparse curl-curl and conductivity
blocks couple the magnetic increment to electric scalar potential:
`J = sigma * (-A_dot - grad(phi) + v cross B_mid)`. The same P1 stage current
enters Ampere's law, charge conservation, Ohmic heat and Lorentz mechanical work.
Exact tetrahedral P1-product integration supplies the quadratic transfers.

A deterministic relative tree gauge removes gradient freedom while retaining
harmonic circulation. Dynamic A is not independently projected into a static
coexact representative. Face-connected conducting components own separate
floating potentials, including distinct potential DOFs at incidental point or
edge contacts. Conductor boundaries/interfaces are naturally insulating.
The sparse nonsymmetric motion operator uses the original bounded, cancellable
FGMRES owner with a recomputed true residual of the actually retained increment. Rounded
zero-physics A increments are admitted only when this residual and all physical
gates still pass; solver tolerances are not relaxed. Full physical Ampere equations,
including gauge-eliminated rows, every component's charge equations, endpoint
insulating equations, discrete div-curl, recoil operating bounds and complete
field/Ohmic/Lorentz/impressed-work balance independently gate acceptance.
Unknown/assembly/Krylov budgets and a caller-admitted cell magnetic Peclet
allowance bound the prepared operation.

Conductivity must retain the identical original magnetic material bundle and
fixed external conditions. Impressed cell current is piecewise linear in time;
nodal velocity is held over a step. Coefficients and Eulerian occupancy stay
fixed. This is a lowest-order, linear/fixed-recoil MQS rung: moving interfaces,
ALE, terminals, displacement current, saturation, hysteresis and thermal
feedback are not inferred. Heat and mechanical transfers require their own
receivers; returned candidates do not silently advance either receiver.

`spatial_induction` checks real floating charge cancellation of uniform motion,
gradient-gauge invariance, an independently integrated charge-eliminated cube
mode and second-order time refinement, analytic three-dimensional motional
feedback/Ohmic/Lorentz work, separated floating conductors, fixed-recoil energy
and source limits, and bounded/cancelled deterministic rollback. Synthetic
materials qualify equations and discretization, not empirical devices.

## Source-backed orthotropic conductivity in native 3D induction

`fs_material::conductor::spatial` supplies a sealed shared conductor boundary
for the original scalar source, three explicit principal sources and heterogeneous
scalar/orthotropic cell maps. `ResolvedOrthotropicConductor3` retains the complete
resolved material bundle, exact conductivity/resistivity basis and property keys,
and finite orthonormal principal directions in physical mesh coordinates.
Resistivity is inverted in the principal frame; the original tensor rotation
owner forms the symmetric positive coefficient. Missing directions, incompatible
units, nonstatic sources and unrepresentable positive tensors are refused.

Both native induction histories, their candidates, temperature refresh and
Maxwell/eggshell load adapters retain the same conductor type. Existing scalar
`new` constructors preserve all-insulator type inference; `new_with_conductors`
admits the shared source owner. One coefficient supplies Whitney/charge mass,
both motional operators, J=sigma E, exact P1 E.J heat, charge conservation and
J-cross-B mechanical work. The actual matrix infinity norm bounds the magnetic
Peclet gate. Scalar arithmetic keeps its original multiplication order.

Temperature refresh keeps every principal key, source basis and exact fixed
frame; it cannot replace orientation or fall back to an unrelated scalar source.
This does not add Hall response, frame evolution, ALE, electrical memory or a
phase-mixing law. Magnetic constitutive and finite-caloric ownership remain in
their existing subsystems.

The focused `spatial_induction::orthotropic` cases specify independently
integrated cube-mode decay and radial-motion heat/Lorentz work for aligned and
oppositely rotated tensors, explicit principal resistivity inversion, source and
frame refusal, scalar compatibility, heterogeneous cells and backward-Euler
numerical-loss separation. The analytical charge-eliminated mass is
tr(sigma)/40-(sigma_xy+sigma_xz+sigma_yz)/120, so a rotation changes the native
trajectory and heat without a change to the physics code.

## Reversible saturation in three-dimensional A–phi histories

`spatial::vector_potential::induction::nonlinear::SpatialNonlinearInduction`
continues an admitted original nonlinear vector field with coupled conductivity,
floating charge and prescribed nodal motion. Backward Euler uses H(B_next) and
J=sigma*(-delta_A/dt-grad(phi)+v cross B_next). The original source's consistent
differential reluctivity enters the original FEEC tensor curl-curl operator and
bounded Newton/FGMRES. The shared `induction::coupling` implementation owns the
relative tree gauge, conducting components, sparse conductivity/motion blocks,
source compatibility and exact P1 current/Ohmic/Lorentz recovery for both
transient discretizations. No independent copy of those equations is introduced.

Original linear, reversible isotropic saturation and fixed-recoil cells may
coexist. Fixed tangential traces and actual prescribed-edge masks are retained.
Every nonlinear trial starts from the same accepted A and source history.
The retained full Ampere and charge equations, div-curl, source/recoil envelopes
and total work independently gate acceptance. Iteration, assembly and Krylov
work bounds and cancellation return no partial state.

The source-exact 3D Bregman loss D=H_next.delta_B-delta_psi is computed from
positive original segment slopes and angular convexity. Reports separately
retain signed magnetic storage change, physical Joule heat, Lorentz work,
impressed work and nonnegative numerical dissipation. D is never thermal heat,
hysteresis or a fitted core-loss term. This is first-order fixed-geometry,
fixed-magnetic-law MQS evolution; vector hysteresis, ALE, terminals, resolved
displacement current and mechanical state advancement remain separate work.

Seven new tests check an independently charge-eliminated BE mode, analytic
saturated decay and first-order refinement, exact multi-knot work/heat/loss,
three-dimensional motion/charge/Lorentz feedback, rotation/gauge covariance,
fixed-recoil energy/envelopes and transactional source/work/cancellation
refusal. Six original linear midpoint induction regressions remain passing.

## Source-pinned temperature refresh of three-dimensional induction

`SpatialInduction::step_with_materials_candidate` evolves the retained linear
A–phi history with endpoint source-backed conductivity and permeability.
Original geometry, Whitney bases and incidence allocations are shared; only
material-weighted operators are rebuilt. The complete original magnetic and
electrical source selections, conductivity basis, temperature coordinate,
fixed traces, conducting support and accepted history remain authoritative.
Per-cell receipts separate signed permeability parameter work from native
Ohmic heat and prescribed Lorentz work. Magnetic free-energy change plus
Ohmic heat plus mechanical work equals impressed work plus parameter work.
The endpoint field and source bundle are retained together; a static re-solve
cannot discard induction memory. This lower owner does not choose temperature
or deposit reversible heat: the finite-caloric consumer must account for
magnetic entropy/internal energy and positive combined heat capacity.

`SpatialNonlinearInduction::step_with_conductors_candidate` provides the same
source-pinned conductivity refresh for retained saturating, linear and fixed-
recoil magnetic histories. The complete original conductor bundle, electrical
basis and selected claims stay fixed; only the source-owned absolute temperature
may change. Original magnetic admission checks the whole H(B) or recoil range,
so equal endpoint coefficients cannot conceal changed magnetic storage. These
magnetic laws must remain temperature independent. Original relative gauges,
conducting components, geometry and A history are retained while the shared
conductivity operator is rebuilt. Local Joule heat, signed Lorentz work and
positive backward-Euler numerical loss remain separate native transfers.
The finite-enthalpy consumer in `fs-couple` determines temperature and deposits
only physical heat; this lower owner never mutates a thermal receiver.

## Maxwell loads from retained three-dimensional induction

`SpatialInduction` and `SpatialNonlinearInduction` now expose
`maxwell_wrench` and `eggshell_wrench` through the same private load bridge.
Original constitutive B/H is recovered directly from accepted A, with no
static re-solve that could discard eddy-current memory after source removal.
The existing closed-surface integration and prepared P1 virtual-displacement
functional retain canonical SE(3), moment-origin, deterministic reduction and
geometry-identity behavior. No duplicate magnetic stress or force model is added.

Every surface neighbor or eggshell transition/fractional plateau must be one
complete source-identical induced-linear exterior, without resolved conductivity
or impressed current. A momentarily idle conductor is not accepted as insulating
air. Conducting, recoil or saturating targets may lie on the enclosed weight-one
plateau. Mesh budgets bound complete surface-field recovery; cancellation and
refusal mutate no accepted electromagnetic state. Returned instantaneous
wrenches do not advance mechanics or debit mechanical work, and overlapping
Lorentz body loads must not be added again.

Two substantive tests use an independently integrated compact P1 current loop.
For manufactured m_z=0.125 and imposed B_x=0.4, raw initial torque_y is
0.0499977666885 N m versus the independent 0.05 N m reference. Retained torque
after impressed-current removal is 0.0383014007176 N m and decreases to
0.0116541343558 N m over four subsequent steps; the fresh static source-off solve
instead has zero torque. Both transient owners and both force functionals,
analytically equivalent linear-limit endpoints, frames/origins, late cancellation
and strict conductor/current/source collar refusals are exercised.
This finite-mesh manufactured result is not a universal force-error bound or
experimental device validation; field/collar refinement remains required.

## Nonlinear scalar material response in the original spatial conduction solve

`conduction::response::PlanarConductionResponse` closes original pinned
conductivity/resistivity sources at the actual recovered |E| and |J| of each
triangle. Its bounded iterations use the SAME FEEC assembly, voltage/current
traces, field recovery and charge/terminal-work admission as
`PlanarConduction`. The ordinary constructor keeps its static-source refusal.
The frozen linearization is private outside the explicit coupled-solve callback;
nonlinear response is never obtained by scaling a converged unit-voltage result.

Accepted candidates retain both the applied material coefficients that produced
the native fields and their complete source evaluations at those fields.
Independent conductivity and current-density defects, plus the sum of absolute
cell power defects, must pass alongside native algebraic/current/work gates.
Numerical iteration defects are not Joule heat. Rejected sources, budgets,
nonconvergence and cancellation leave the accepted owner unchanged.

This fixed-geometry scalar memoryless rung accepts prescribed native voltage and
current-density traces. A nonlinear equipotential total-current terminal needs
its own constraint; tensor, Hall, moving-conductor, field-storage and hysteretic
paths retain their own owners. Positive source transport does not prove Picard
convergence, constitutive uniqueness or dynamic stability.

Independent regressions exercise two-region nonlinear series continuity,
current-driven resistivity, simultaneous E/J dependence, data-only source swaps,
local E.J power, rotation/reversal/gauge invariance, constant-source identity,
separate power admission, continuation and refusal rollback.


## Total-current electrodes with original source response

`PlanarConduction::solve_current` imposes signed total current on complete equipotential electrodes using the original frozen FEEC operator. Actual-amplitude recovery checks both terminal currents and authored I*V independently; source-dependent callers repeat this constraint inside the existing E/J closure. `PlanarConductionResponse` supplies ordinary and thermally refreshed current candidates without unit-field material queries or scaling accepted nonlinear fields. Actual local E.J remains heat; constitutive defects remain numerical residuals. Six manufactured regression groups cover series and parallel heterogeneous specimens, nonuniform electrode reaction, data-only source changes, both signs/zero, constant-law parity, source/work gates and cancellation/retry. Bounds remain fixed geometry and memoryless scalar response; convergence and uniqueness are not guaranteed.

## Original source response under coupled electrical drives

`PlanarConductionResponse::solve_coupled_with_materials_candidate` admits an
external owner's actual native field through the same source iteration used by
voltage/current boundaries. The returned field must borrow that exact frozen
operator; another operator's field refuses. Each trial retains its proposed
driver state, but only the state paired with the admitted constitutive field
is returned. Circuit history, source work and external work limits remain with
the caller; response iteration does not authorize advancing physical time.

The complete original laws and applied operator can be transferred together to
a continuing owner and rejoined without matrix assembly. Original source
admission rejects mismatched law/operator pairs. The original conductivity,
current-density and local-power defect gates and cancellation behavior remain
unchanged. The conduction regressions compare direct and coupled execution and
check source continuity and cancellation after a native trial.

## Original surface contacts inside nonlinear bulk response

`PlanarConductionResponse` now assembles its original scalar E/J-dependent bulk
sources and source-backed surface contacts through the same native FEEC operator.
Each actual-amplitude voltage, total-current or external-circuit trial preserves
the contact traces, voltage jumps, current continuity and exact surface heat.
Contact losses remain separate from bulk E.J in the report and are included
once in total terminal power. Ordinary static bulk admission is unchanged.

The coupled material/contact-context method retains original surface cards,
selected claims, ordered material assignments and coordinate meanings. A thermal
owner may refresh its explicitly owned temperatures and phase fractions; no
source replacement or electrical-field-dependent interface law is inferred.
The same original sources and contact geometry survive continuation. Focused
manufactured regressions compare independent nonlinear series equations under
voltage and current control, surface-source changes, both polarities, local and
total heat, and missing/replaced-source refusal. Interfaces remain stationary,
Ohmic at each supplied context and without storage; nonlinear or hysteretic
contact laws require their own physical state and closure.

## Exact planar dielectric reflection

The opt-in `plane-wave` feature exposes
`plane_wave::PlanarDielectricInterface`, an exact monochromatic boundary solve
for two homogeneous isotropic half-spaces. It retains both original material
states. Incident epsilon and both permeabilities use the existing lossless
field-storage admission. Transmitted epsilon is either the explicitly chosen
nondispersive law or the original `ResolvedDebyeCapacitor` response, converted
from terminal admittance to bulk permittivity without a second Debye law.

The convention is exp(+i omega t), with a forward transmitted wave
exp(-i k_z z) and Im(k_z) <= 0. TE and TM amplitudes refer to the same
tangential electric-field axis on both sides. The TM coefficient therefore
differs in sign from a rotating optical p-polarization basis. Returned power
fractions are normal Poynting fluxes; transmitted power is evaluated at the
interface, before subsequent absorption in the half-space. Total internal
reflection retains its complex phase and evanescent field. The exact critical
angle avoids division by vanishing normal wavenumber.

Synthetic source-bound tests cover normal incidence, magnetic impedance,
Brewster cancellation, independent Poynting balance, total internal reflection,
critical incidence, grazing transparency, Debye phase and capacitor-geometry
invariance, retained source identity and explicit refusals. This is not a mesh
Maxwell solver, transient wave model, layered-medium solver, or experimental
qualification of an otherwise static material at an optical frequency.

### Coherent finite dielectric films

`plane_wave::layers::PlanarDielectricLayers` prepares up to 64 ordered finite
films between an explicit lossless incident source and an explicit lossless
or Debye exit source. Every film retains its original `DielectricResponse`,
thickness and sourced permeability. Constitutive data are cached once at the
selected frequency; angle/polarization queries are allocation-free. Zero
thickness is an exact identity after source admission. The half-space owner
and this stack share material-contrast and outgoing-normal-wave calculations.

Composition retains bounded homogeneous E/H states in common real reference
ports, then reconstructs forward scattering gains. It preserves resonant
coupling without subtracting products of rounded nearly unit reflections.
Propagation uses only the decaying/oscillatory exponential; an analytic
zero-normal-wave limit includes exact internal critical angles. Returned
tangential-E phases refer to the first and last film faces. Total internal reflection, finite-gap
evanescent tunnelling, Fabry-Perot resonance and coherent antireflection are
retained. Exit power is the normal Poynting flux entering the last half-space;
finite-film absorption is reported separately. Dissipative even/odd slab modes
retain weak Debye loss that subtracting reflected/transmitted powers would
round away. Wholly lossless films have exactly zero dissipative work.

The existing plane-wave tests add half-space/zero-thickness reductions,
quarter-wave and resonant phases, grazing subdivision and distinct-material
resonances, independent E/H transfer matrices with magnetic media, original
Debye/geometry invariance, weak-film loss, critical films/exits and thick absorbing/evanescent limits. This remains a coherent
homogeneous isotropic planar model: finite coherence length, roughness,
anisotropy and mesh Maxwell evolution are outside its scope.

### Coherent vector beam diffraction through the same films

`plane_wave::angular_spectrum::PlanarAngularSpectrum` accepts complex peak Ex/Ey
samples at the incident face of an original `PlanarDielectricLayers`. The input
is one explicitly periodic, band-limited rectangular cell in physical metres;
x varies fastest. Both grid dimensions are powers of two, capped at 4096 to
bound a single uninterruptible FFT pencil. All incident Fourier bins must lie
strictly inside the incident light cone. No modes are silently removed,
windowed or renormalized. A caller approximating an isolated beam must resolve
its bandwidth and converge the cell/padding independently; an authored hard
aperture's unresolved high spatial frequencies are not synthesized.

The original `fs-fft` plans analyze both complex components once. Every mode
uses its actual conserved transverse wave number and its own TE/TM response
from the original source/material/film owner. Complex amplitudes combine before
forming intensity. Source epsilon and mu remain distinct; scalar refractive
index or normal-incidence transmission cannot replace their vector response.
For the existing exp(+i omega t) convention, normalized inverse FFTs analyze
the field and unnormalized forward FFTs reconstruct exp(-i k.r) modes.

Each subsequent plane reuses the cached modes and plans, applies exact kz
propagation, enforces k.E=0, and obtains H=(k cross E)/(omega mu). It returns
all six complex field components and physical cycle-mean Poynting flux.
Incident, reflected, their coherent sum, and transmitted fields have explicit
face/distance conventions. Film-induced evanescent fields decay without a
growing transfer matrix. An exactly grazing exit is refused because the
supplied tangential-electric parameterization becomes singular there; the
underlying film solver retains its internal critical-angle limits.

Modewise power integrates over the original cell: reflected outward power,
exit-face transmitted power and physical film absorption are separate. Loss
after propagation into an absorbing exit is computed from the actual field,
without relabeling it as film absorption. Local flux retains interference and
its sign. Numerical-storage and FFT-work budgets precede allocation; source
response work is bounded by the admitted grid and original 64-layer maximum.
Cancellation polls between rows, pencils and bounded mode batches, publishing
no partial field. Further planes require six separable FFTs and linear modal
work, with no repeated constitutive or film solve. This is an operation-count
description, not a measured throughput claim.

Native leaf tests exercise both Maxwell curls/divergences, polarization and
signed power, and evanescent decay. Source-backed integration regressions cover
original-film reduction, an independent direct Fourier diffraction series,
TIR, Debye film/exit loss, budgets and cancellation. Their Cargo execution
requires the workspace constellation dependencies. This feature does not add
transients, incident near-field excitation, rough/finite film edges, arbitrary
three-dimensional scattering, or empirical material-frequency qualification.

### Two-sided coherent illumination and layer heating

`PlanarDielectricLayers::prepare_illumination` retains the same source-bound
film operator and prepares simultaneous illumination from its two physical
faces. Inputs are peak tangential electric fields in V/m, with one frequency,
transverse wavevector, polarization and exp(+i omega t) clock. Both exterior
channels must propagate in lossless media. Internal films may still be lossy,
evanescent, exactly critical or opaque; the existing one-sided API retains
its separate admission of lossy and total-internal-reflection exit media.

Two bounded sweeps prepare the complete physical E-field scattering matrix
and each film's positive even/odd loss transfers. Evaluation combines complex
source fields before computing outgoing powers or heat, including phase
cancellation and enhancement. Returned boundary E and power-conjugate H are
in V/m and A/m; port powers and each authored film's heat are in W/m². The
different exterior modal admittances enter transmission and reciprocity.
Layer heating comes directly from dissipative modes, so weak absorption is
retained without subtracting nearly equal incident and outgoing fluxes.
Wholly lossless films have exactly zero heat.

Prepared coherent evaluations allocate nothing and cost O(layers). Explicit
source factors describe a positive second moment C=B B*, with one peak-phasor
time-average factor 1/2; an independent-source helper supplies only diagonal
factors. Statistical evaluations preserve interference within each factor,
sum mutually uncorrelated factors, and obey caller mode/layer-work budgets
and cancellation. These are monochromatic ensemble averages, not spectral
densities, an inferred bandwidth, or a spatial coherence model. The report
does not advance temperatures or debit heat into a thermal owner; a coupled
consumer must book the returned finite-film absorption exactly once.

Focused regressions compare simultaneous-source interface continuity,
admittance-weighted reciprocity, independent E/H boundary equations and
per-layer Poynting work; phase-controlled thin-film heating; full complex
source covariance; and opaque, critical and evanescent internal limits.
Material frequency qualification and non-planar Maxwell evolution remain
outside this analytic boundary model.

## Source-backed nonlinear electrostatics and differential capacitance

`electrostatic::nonlinear::PlanarNonlinearElectrostatics` retains the original
planar P1 mesh, physical extrusion depth and material sources. Linear scalar or
explicitly framed tensor cells can coexist with the shared reversible radial
`D(E)` material law. The voltage functional integrates `D.dE`; reported physical
storage integrates `E.dD`. Their sum, rather than twice either quantity, is
checked against the actual volume-charge and electrode work pairing.

Prepared complete exterior-edge electrodes accept voltage or total free charge.
A charged floating conductor has one unknown potential and its own total-charge
equation. Shared traces can join disconnected mesh pieces; each resulting
component must have an actual voltage reference. Relative potentials preserve
representable small fields under a large common voltage bias. Volume charge,
individual floating-charge constraints, free-node residuals and total charge
balance are checked independently before a field is returned.

The analytic radial material tangent distinguishes longitudinal differential
permittivity from transverse secant permittivity. Bounded Newton steps use the
original checked-CG owner, native interval coefficient admission and diagonal
equilibration of the equipotential quotient. Prepared cell geometry applies the
tangent without assembling a global matrix. Line search uses direct material
coenergy increments, including source work, rather than subtracting two large
stored objectives. Material support, cancellation and finite nonlinear/linear
work budgets remain hard refusal boundaries; rejected trials publish no field.

`NonlinearElectricField::differential_capacitance` differentiates that accepted
operating point with volume charge and every floating total charge held fixed.
One prepared tangent is reused across voltage-terminal columns. Recovered
charge, reciprocity and gauge residuals must pass in farads. The result retains
its field and source identity; it is a local `dQ/dV`, not a global `Q=C*V` law.
A genuine source-law kink or exhausted field support refuses a unique two-sided
differential response. These fixed-geometry, memoryless, reversible laws do not
infer dielectric relaxation, hysteresis, thermal evolution or electrostriction.

## Nonlinear electric forces and original discrete virtual work

`NonlinearElectricField::maxwell_wrench` uses the original shared geometric
contour owner with an explicitly frozen isotropic constitutive potential.
Its stress is `E tensor D - w* I`. Parallel and transverse principal stresses
use independently integrated physical energy `w` and coenergy `w*`, retaining
small longitudinal pressure without subtracting two large biased quantities.
The original linear magnetic/electric callers retain their existing dyadic
pressure arithmetic through the same contour geometry and first moments.

Every contour edge requires one identical complete isotropic material source
on both sides, zero original free-charge density, and no electrode vertex.
Actual extrusion depth, winding, frame, origin, work bounds and cancellation
remain with the original owner. The enclosed region may contain charge and
contrasting materials. This gives a net instantaneous native `Wrench`, not an
invented solid-surface traction map or an additional overlapping rho E load.

`shape_forces` differentiates the original accepted mixed electrical P1
functional with respect to vertex positions and physical extrusion thickness.
Voltage controls and floating total charges stay fixed. Each cell transports
its original total free charge, so its density changes inversely with volume;
the consistent q/3 source remains unchanged. Isotropic material contexts stay
fixed. Anisotropy refuses until its deformation/frame convention is supplied.
The result retains the field and its actual controls. `virtual_work` contracts
the same derivative directly over cells without allocating a nodal force
vector; local displacement differences preserve rigid-translation invariance.
Its virtual displacement is in metres and its first-order result is in joules.

These envelope derivatives require no extra field solve or tangent assembly,
but inherit the accepted primal residual and spatial error. Independent
fixed-voltage/coenergy and fixed-charge/energy perturbations test the declared
held variables. A source-free nonlinear radial field tests the pressure term
against a distinct half-pairing falsifier. Refinement and original linear force
parity are implementation evidence, not a force enclosure or measured material
validation. Electrostriction, thermal/material evolution, ambient mechanical
pressure, finite moving-mesh work and a mechanical time step remain with their
corresponding physical owners.
