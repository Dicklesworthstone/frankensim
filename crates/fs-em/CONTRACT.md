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
model. End fringing, electric stress and distributed magnetostriction remain
unsupported. The fictitious integration contour is not a physical body
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
