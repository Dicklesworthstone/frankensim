# CONTRACT: fs-plate

## Unfactored full tensor-material shell operators

`TensorPronyShell::prepare_harmonic_operator` evaluates all original
tensor-Prony/Kelvin sections at one admitted angular frequency and returns a
`PreparedTensorPronyShellOperator` without factoring the bare structure.
Coupled physical owners can retain the entire shell field through a bare-shell
resonance. Original supported translational/rotary mass is borrowed unchanged;
stiffness includes complete material storage and the original lossless drilling
term; damping contains the complete material loss divided by angular frequency.
No modal basis, material refit or load is introduced.

The operator retains its source identity and exact frequency. Its
`admit_displacement` requires every original finite free coordinate and then
performs the original facet/layer full-cycle strain admission. A coupled owner
must apply that admission before publishing a physical response and separately
check its own residual and work balance. These frequency-specific matrices do
not define a frozen causal time material or its frequency derivative.

Preparation admits a conservative logical operator payload and the original
spatial/material bounds before allocation. The existing bare and port-loaded
harmonic paths share this assembly, preserve their dense-factor admission and
discard temporary sparse duplicates before factoring. Their existing material
frequency-derivative and field-recovery paths remain in force. Two additions to
`tests/shell_prony_harmonic.rs` check original mass identity, analytical storage
and loss across frequencies, prescribed full-field recovery, and original
band/strain/dimension/budget/cancellation refusals.

## Physical viscosity on CST/DKT shells

`shell::viscosity` supplies complete passive membrane, membrane/bending and
bending loss on the existing flat-facet shell field. Sparse spatial assembly
uses the original facet rotations and global-axis support map. Direct
projection first evaluates the same CST strains and DKT curvature quadrature,
then the original constitutive Gram channels, retaining weak loss, exact null
channels and every modal cross term. Neither operation changes elastic or
inertial operators or dissipates the numerical drilling gauge, prestress or
unrelated reinforcement. Symmetric A/B/D passivity uses the existing physical
six-component resultant owner shared with MITC; no second constitutive solver
or fitted modal damping law is introduced.

The thickness strain is `epsilon(z)=epsilon0+z*kappa`; the native DKT
`beta=-grad(w)` field already supplies its signed first-moment curvature.
This preserves signed B work for offset layers while leaving all D-only
quadratic work unchanged. Supplied resultants are in each original facet frame
with explicit SI units. Material source applicability, frequency support and
thermal state remain caller obligations. Projection yields `Phi^T C Phi`, with
units 1/s only for a mass-normalized basis; it does not certify modes or a
modal truncation.

Checked element, full/free-coordinate, mode and work bounds precede output
allocation. Cancellation at element and projection-row boundaries publishes
no partial result. The public `ShellModel` does not retain its mesh: callers
must supply its original CST/DKT geometry. Dimension/map checks do not prove
identity, and MITC oblique supports or tied assembly coordinates require their
own physical expansion owners. Focused regressions cover independent affine
strain work and B reversal, rigid/drilling nulls, original DKT parity, facet
rotation and supports, weak semidefinite loss, bounds and cancellation.

## Physical ply shear on the native MITC field

`MitcShellModel::prepare_laminate_shear` binds an original equilibrated laminate
recovery to the exact oriented section, inertia and reference surface used in
one native facet. It caches the existing condensed strain observation once;
each physical nodal snapshot forms one native polynomial, then arbitrary ply
and barycentric point queries allocate nothing. Actual ply stresses and strains
come from the material owner's continuous equilibrium profile, retaining full
off-axis coupling and physical shear work. Responses expose section, ply and
facet coordinates explicitly. This observation changes no assembled operator.

`PreparedThermalShell::laminate_shear_field` instead reuses the same thermal
owner's strain field, including thermally induced interior director motion at
fixed physical nodes. Local work integrates to the native H contribution, not
to membrane/bending or numerical drilling energy. Source mismatch, invalid
positions and cancellation refuse without advancing any owner. The source's
two-probe cylindrical-bending approximation remains binding; this is not a
general 3-D/free-edge or failure/delamination model. Focused native tests cover
mixed shear work, a loaded nonuniform facet, thermal director motion and refusal.

## Original MITC thermoelastic factors

`shell::mitc::thermoelastic::project_thermoelastic_factors` projects the
original homogeneous symmetric slab's reversible membrane stiffness and
bending relaxation factors separately. It reuses the original physical
six-coordinate nodal basis, directors, elastic director-bubble condensation
and positive seven-point element quadrature. Every retained cross-mode term
is represented by the full factor Gram. The source owner supplies the
constitutive square factors already in the original facet's engineering
strain frame; no thermal property, shear loss or drilling loss is inferred.

The uniform through-thickness thermal mode stores membrane energy without
dissipation under insulated faces and absent lateral diffusion. Odd thermal
diffusion modes belong to the caller's positive bending-memory spectrum.
Regional row metadata permits one element pass followed by separate grouping
for distinct material rates. Exact zero channels remain zero; compensated
short physical-strain and factor contractions retain weak signed channels.
No dense nodal thermal matrix, eigensolve, spectrum or time integrator is added.

The complete supplied basis is admitted up to the existing 1,024-mode shell
ceiling. Checked work and combined factor-entry allowances, including
quadrature padding, precede allocation. Original geometry/section validation
and cancellation remain binding. Static elastic bubble condensation is the
retained trial-space approximation, not a dynamic elimination of thermal
director states. This local slab projection is not a general layered thermal
model, a source-applicability proof or a spatial/roundoff certificate.

Five authored regressions cover affine membrane energy, independent bending
energy, signed cross terms, weak channels and exact nulls, curved director
condensation and frame covariance, complete-basis admission and refusal.
Independent numerical equations were checked; native Rust execution remains
unverified in this environment.

## Bonded-beam stretch in the shared shell equilibrium

The single-shell and tied-assembly preload owners include the original bonded
Euler--Bernoulli beams' Green--Lagrange centerline strain
`e = u' + (u'^2 + v'^2 + w'^2)/2`. Their original linear axial energy is
replaced by `EA/2 * integral(e^2) dx` through polynomial energy/force/tangent
corrections, without subtracting nearly equal full operators. Axial stress is
generated by the actual joint equilibrium; it is neither prescribed afterward
nor added twice to the vibration pencil. Compression remains signed and the
existing tangent-inertia admission refuses an unstable accepted increment.

The beam's original centroidal frame, cubic displacement fields, signed bending
rotations and eccentric rigid-arm map are shared with linear assembly. Five
Gauss points integrate the degree-eight axial storage exactly. Prepared strain
rows and the original component mortar maps are reused by every load increment;
matrix-free Jacobian actions apply the local derivative directly without a
dense beam tangent or dense tie transform. Sparse tangent assembly mirrors
identical contributions. A Bernstein convex-hull bound admits the gradient
over the complete segment, including extrema outside integration points.

Zero displacement preserves the original elastic pencil, and reference mass
is unchanged. Aggregate preload work admits at most 40,000 bonded segments.
The original load continuation, cancellation, force-balance, source/reference,
moderate-gradient and physical-rotation limits still apply. Bending, torsion,
offsets and glue retain their reference-frame linear laws: this does not supply
finite-rotation beam kinematics, shear deformation, torsional buckling,
post-buckling continuation or a new material constitutive law.

Regressions cover independent signed axial and Hermite geometric work, exact
quartic reinforcement energy in a loaded shell, energy/force/tangent
differences, sparse overlapping tie rows, complete-segment admission, and a
three-rung fourth-order strain-energy ladder. An eccentric tied beam chain
matches the monolithic loaded shell, including eliminated-seam loads; one-part
equilibrium and mass retain exact parity with the original owner. The strain
convention follows the centerline Green-strain construction in Wu and Tiso
(2016), doi:10.1007/s11044-015-9476-5, extended to both transverse directions.

## Harmonic resonance of the accepted loaded assembly

`shell::preload::assembly::vibration::harmonic::AssemblyHarmonicModel` caches
the complete original `Phi^T M Phi` and `Phi^T K_t Phi` once, using the same
projection as `AssemblyResponse`. Finite mass-orthogonality and eigenvector
residual terms remain present. Every frequency factor uses the shared
`fs-modal::harmonic` owner with all projected glue damping; no artificial
diagonal loss, modal retuning, mass replacement or second static preload is
introduced. Factors are reusable for arbitrary phased source/receiver cases.

Peak phasors use `Re(x exp(-i omega t))`. Complex physical component loads
pass through the original support/tie transpose, and both displacement
perturbations and velocities expand through its original forward map. The
accepted loaded mean shape remains unchanged and separate from the harmonic
perturbation. Physical component cycle power is recomputed independently and
must agree with retained-coordinate input power. The reported viscous power
is glue loss only; material/radiation loss is not inferred.

`PreparedAssemblyHarmonicResponse::solve_with_frequency_derivative` reuses the
same full factor for the exact source/response frequency tangent. Optional
physical force/moment slopes pass through the same support/tie map; both motion
and its derivative expand to original component coordinates. The original
source power and its derivative are independently recomputed there and checked
against the modal owner. No frequency perturbation or new factorization is
needed. Geometry, loaded equilibrium and the retained basis stay fixed. This
gives physical shell receivers the amplitude/slope pair used by fs-modal's
local driven-resonance refinement. The one-mode loaded-shell test checks full
physical displacement/velocity slopes, including a phased source derivative,
at DC and on both sides of resonance; malformed slopes and cancellation refuse.

Mode, projection-work and projection-payload limits are explicit, followed by
the shared frequency solver's separate factor limits and residual threshold.
No mode is silently removed. Cancelled or refused calls publish no result and
preserve both the loaded body and the immutable prepared factor. Singular
undamped frequencies refuse through the shared solver without damping jitter.
The response is limited to the supplied spectral windows and unchanged
linearized connection branches; it is not a modal-tail or nonlinear-validity
certificate. G1 tests compare one-mode loaded resonance/phase reversal with an
analytic oscillator, reconstruct original FE equations and physical glue heat,
and demonstrate the error from dropping cross-mode damping. G0/G4 tests cover
nonmatching slave loads, caps, malformed data, cancellation and deterministic
factor reuse. These are synthetic numerical checks, not specimen validation.

`harmonic::ports::AssemblyHarmonicPorts` additionally binds original physical
unit-effort attachment shapes and caches their complete `Phi^T T^T` force map.
Each supplied complex impedance maps conjugate port flow to resisting effort;
the shared solver assembles its full reactive, mutual and nonreciprocal loading
before factorization. The source M/K and glue C remain unchanged. Original
component reaction forces, physical port flows and load work are reconstructed
independently from the physical velocities, including supported and eliminated
slave DOFs. Physical applied power must equal glue dissipation plus external
load work. A Maxwell-port oscillator test checks frequency-dependent stiffness,
phase and power around the original loaded resonance; a tied-body test retains
mutual loading and tests the error from dropping those entries.

The loaded port factor also supplies `solve_with_frequency_derivative`: the
original component source slopes and complete constitutive `dZ/domega` enter
the shared same-factor tangent. An explicit zero matrix declares a locally
constant load; missing/nonfinite entries refuse. Both original physical flow
and reaction derivatives are reconstructed, including supported/slave DOFs,
and physical source/glue/load work slopes are independently balanced. Output
and observation limits precede physical allocation. Fixed component receivers
use the shared compensated transfer projection; physical port flow and effort
also expose amplitude/slope pairs for local resonance/antiresonance refinement.
The receiver, geometry, accepted branch and retained basis stay fixed. These
local response extrema do not certify complete spectra or modal-tail error.
Focused tests compare an analytic dispersive Maxwell-loaded shell at DC and
around resonance, distinguish its displacement and flow peaks, and compare a
phased nonmatching-tie/mutual-port tangent against frequency differences.

The port owner consumes the caller's attachment vectors without copying them,
with separate cold projection/payload limits and common per-frequency factor
limits. Pointwise Hermitian PSD admission is a numerical check only. Port order,
units, impedance provenance, source band, and ownership of each physical loss
remain explicit caller responsibilities. No acoustic surface, material law,
external storage/heat decomposition, or broadband validity is inferred.

## Shared original source for coupled resonance transients

`vibration::response::AssemblyResponseSource` owns the original complete
projected mass/stiffness, inverse mass, glue damping and canonical pHS source.
The existing `AssemblyResponse::prepare` reuses this owner without changing its
standalone held-force propagation. `AssemblyHarmonicPorts` exposes its borrowed
accepted vibration and original physical attachment vectors so coupled solvers
reuse the same body, mode ordering, support/tie transforms and force map.

The fs-couple exact assembly/MIMO transient uses this source to evolve structural
motion and original rational-load memory together. It reconstructs full physical
displacements about the unchanged loaded mean, velocities, reactions and load
impulses, and independently checks supplied work against structural energy,
original glue dissipation and signed external work before publishing either
history. The original projection budgets and retained-mode/fixed-branch
limitations remain; this extraction introduces no new FE operator, preload,
mode window, material law or external storage/heat certificate.

## Purpose and layer
Orthotropic thin-plate bending for instrument bodies: DKT (discrete
Kirchhoff triangle) elements with an orthotropic bending-rigidity matrix
from `fs_material::OrthotropicElastic`, membrane-prestress geometric
stiffness, offset stiffener beams, and modal analysis through
`fs_modal::slice_window` (certified intervals, inertia-certified counts).
Layer: **L3 FLUX** (deps include fs-material L3, fs-modal L1, fs-sparse L1,
fs-rep-nurbs L2).
Bead frankensim-fsim-plates-shells-kj3s0 (musical-acoustics program).

`shell` already owns spatial membrane/bending mechanics (the legacy CST/DKT
faceted operator), independently specified A/D sections, preload and radiation
projection. `shell::mitc` adds the separate shear-deformable triangle for IRP
M1.2-W1a. The legacy `assemble`, DKT stiffness and drilling constant are
unchanged. This new path never selects DKT as a fallback.

The new element uses Lee, Lee & Bathe (2014),
doi:10.1016/j.compstruc.2014.02.005, Eqs.6–8 and 15–17/Table 1: three physical
displacement/rotation nodes, two internal director rotations, cubic rotation
interpolation, six interior covariant shear tying points and the published
nonzero twisting distance 1/10000. The bubble stiffness is statically
condensed; the consistent mass undergoes the same congruence. This is a
small-strain, faceted-midsurface, linear-through-thickness shear-shell rung,
not an exact curved CAD chart, finite-strain shell, or 3-D solid cross-check.
Geometric director-gradient terms enter bending, including the displacement
terms needed to preserve rigid rotation on a curved component.

`NodalDirectors::new` requires explicit finite unit directions, a positive
minimum facet-normal dot product, alignment with every incident facet, and
the same bound between incident facets. A kink outside that declared envelope
refuses; junctions require separate component surfaces and an external tie.
No director smoothing or material defaults are inferred. `ShearShellSection`
requires independent positive-definite membrane A, bending D and transverse
shear, an explicit positive shear correction, areal mass, rotary inertia and
material angle. Its solid constructor requires BOTH transverse shear moduli.
`new_coupled` and `from_laminate` retain the original material section's full
membrane/bending coupling B and first mass moment I1 about the declared mesh
reference surface. The coupled stiffness block and inertia moments must be
positive definite. Both strain cross terms enter the original MITC quadrature
and elastic bubble condensation. Consistent laminate inertia integrates the
same director motion `u + z*(theta cross director)`, including translation/
rotation cross work; nodal lumping retains the nodal cross terms. Drilling
inertia remains an explicit numerical gauge in the existing I2 scale. The
homogeneous and complete laminate constructors use the SAME physical director
interpolation. This corrects the former homogeneous angular-vector rotary
interpolation on varying directors, while preserving flat-layer equivalence.
Consequently a vanishing coating approaches the bare curved shell continuously;
a varnish comparison cannot switch its underlying inertia discretization.

The source laminate already includes its transverse shear correction: it is
not multiplied by 5/6 a second time. Layer order, reference surface, density,
and orientation remain caller data. Prony layers lower their retained reference
storage into this fixed elastic pencil; this does not turn a frequency-domain
reference into a broadband hereditary solver. The coupled Kelvin--Voigt
viscosity constructor retains viscous B, checks passivity of the complete
six-component A/B/D block, and projects it through the SAME elastic trial space.
Neither B nor its viscous counterpart is required to be positive by itself.

`tests/laminate_shell.rs` covers flat homogeneous equivalence, signed B work,
first-moment kinetic work, reference-surface invariance, varying-director
kinetic energy, coupled admission and axial-load-induced curvature on three
meshes against closed-form CLT. Inline
viscosity tests check full proportional work and semidefinite coupled loss.
These are G0/G1 numerical checks, not validation of a historical varnish film,
audible superiority, or a particular instrument specimen.

New material-frame trigonometry uses
`fs_math::det`; `SheetMatrix::rotated` is not used by this path.

Supports declare a component mask and either global axes or a proper local
orthonormal frame. Repeated support declarations on one node refuse rather
than silently combining frames. A `MitcShellModel` retains an ordinary
`ShellModel` pencil plus the nodal support-coordinate transforms. Its
`project_load` and `expand` are the work-conjugate physical-coordinate map;
the pencil's `dof_map` alone must not be read as global coordinates when a
local support is present. `modes_shell(model.pencil(), ...)` uses the existing
eigensolver, and physical mode observation must use `expand`. The separate
wrapper avoids changing the legacy public `ShellModel` layout or existing
CST/DKT consumers.

`shell::reduction::rigid::mass_normalized_modes` accepts both lumped and
consistent symmetric positive definite shell inertia for a complete free pencil
in global physical coordinates. It preserves all supplied translation/rotation
cross inertia during two-pass mass orthonormalization. Consistent inertia uses
one shared sparse AMD/LDLT factor for every candidate's true
`sqrt((K phi)^T M^-1 (K phi))` residual; diagonal inertia keeps its direct path.
The original signed stiffness is retained, so supported or unstable candidates
receive their actual residual rather than a fabricated zero-frequency claim.
`mass_normalized_modes_with_options` exposes the existing sparse factor budget
and cancellation without publishing a partial basis. This remains a linearized
rigid-pose preparation helper, with geometric nullity and rotation limits still
enforced by the reduction owner. `tests/shell_rocking.rs` checks six-mode
orthogonality and continuum affine kinetic work on asymmetric heterogeneous
facets, signed residual oracles, and physical/budget/cancellation refusals.

The third nodal rotation is a numerical drilling gauge. The physical element
depends on `theta × director` and is blind to a node's director-axis rotation.
One area-averaged relative membrane-spin equation per node fixes that gauge;
unlike multiple facet penalties, it does not impose incompatible spin targets
on a curved fan. Its positive strength is an explicit numerical input, not a
material modulus. Both the flat patch and every finest benchmark carry the
required alpha=1e-4/1e-3/1e-2 displacement-insensitivity check (<1e-4 relative).
Consistent mass integrates the physical reference-area/director velocity with
the supplied `I0`, signed `I1`, and `I2`, using the original enriched
interpolation and bubble condensation. Its numerical drilling inertia acts on
the same single area-averaged relative membrane-spin equation per node as
stiffness, with weight `sum(area * I2 / 3)` over the incident facets. It
contributes zero work to every rigid velocity and preserves physical
translation/rotation cross inertia on curved fans. `ShellMass::Lumped` retains
its original separate nodal matrix. The drilling inertia remains a coordinate
regularization, not a measured rotational material law.
`tests/mitc_inertia.rs` checks the complete rigid 6-by-6 mass against independent
barycentric moment integration for flat and varying-director shells, both signs
of `I1`, rotated/translated fixtures, gauge separation, supported positive mass,
and unchanged lumped inertia. These checks certify the declared
reference-area/director model; they do not add a curved three-dimensional
volume-Jacobian model.

`assemble_mitc3p_with_reference_inertia` retains three physical
`ShearShellInertia` values per facet in triangle-corner order. Each contains
original reference-area moments `I0`, signed `I1`, and `I2` about the same mesh
surface and directors as the elastic section. The positive-definite moment
matrix is admitted at every corner and interpolated linearly, preserving its
positive cone. The admission compares finite binary products exactly, avoiding
singular matrices hidden by rounded square roots and overflow/underflow of
floating determinant products. This field supplies inertia while the section
retains its original A/B/D/shear stiffness. The existing positive five-by-five
quadrature integrates its degree-seven enriched velocity products exactly;
consistent mass retains full translation/director cross inertia and the original
bubble congruence. Constant corner values use the uniform arithmetic path.

Explicit lumping and the numerical drilling weight use the exact positive
partition `integral N_i I_k dA = area*(2*I_k,i + I_k,j + I_k,l)/12`, rather than
assigning the facet mean to every node. `reference_inertia_moments` exposes the
retained physical profile. `EnrichedMitcShell` restores its complete physical
inertia from that same profile, and `with_material_model` refuses any change
to the original corner inventory, including relocation at unchanged totals.
Temperature-dependent stiffness never supplies replacement reference density.
`tests/mitc_reference_inertia.rs` covers independent rigid six-by-six moments
on flat and varying-director fans, transformed geometry and both signs of
`I1`, constant-path bits, equal-total-mass frequency changes, exact interior
bubble inertia, material refresh, positive lumping and typed admission errors.

`MitcBudget` checks node/triangle caps before assembly allocation. Assembly
polls its caller cancellation callback before work, per element and per
quadrature point and publishes no partial pencil. Static solves use only
`fs-sparse` factorization, refuse singular/indefinite stiffness and require
an independently recomputed true residual <=1e-8 relative. They project full
world force/moment loads and return full world displacement/rotation values.
Sequential COO assembly and deterministic math specify same-host replay;
no new cross-ISA bitwise claim is made.

The new source tests are `shell::mitc::tests` and
`tests/shell_obstacle_course.rs`: G0 rigid motions, membrane/bending/shear
patches, mass, directors, support frames and cap/cap+1 refusals; G1 thin-limit
square-plate convergence; G3 node relabeling; and the Scordelis–Lo, pinched
cylinder, cut-out hemisphere and full-hemisphere G2 pointwise envelopes.
The obstacle target retains regular and counter-seeded distorted ladders,
observed orders and drilling checks. A separate report logs the legacy DKT
model-form rung without making it the acceptance path. Source presence and
formatting do not establish those numerical results: central native execution,
the existing CST/DKT and piano/percussion tests, both WASM cones and the
unchanged string golden remain required before .47/DEF-11 qualification.

For the hemisphere rows, convergence means strictly decreasing absolute
successive displacement changes over the last three rungs, with a signed
Richardson/Aitken limit (negative step ratio admits alternating convergence).
The cut-out limit must agree within the declared band with both MacNeal–Harder
1985 (0.094) and Ko et al. 2017 Tables 15/16 (0.093); the finest displacement
uses the latter normalization. The full hemisphere uses Table 17 (0.0924).
Regular rows retain 2% bands and distorted rows 5%; the distorted cut-out
retains N=4,8,16,32 and adds N=64 to resolve its non-contracting N=32 step.
Both regular hemisphere tests reassemble and solve the same ladders with D
multiplied by 1.05, requiring the extrapolated-limit reference check to reject
that planted defect. Membrane/shear matrices and geometry stay fixed in this
falsifier. Roof/cylinder reference-error monotonicity is unchanged. These
limits are pointwise estimates, not certified continuum-error enclosures.

## Public types and semantics

`shell::mitc::viscosity::ShearShellViscosity` supplies independent symmetric
membrane, bending and transverse-shear viscous resultants in physical facet
coordinates. Each tensor is positive semidefinite, including directional-only loss and zero; missing
couplings are not inferred. The shared strain kernel evaluates elastic and
viscous bilinear forms with identical MITC kinematics and quadrature.
`project_viscosity` uses the ELASTIC bubble transformation for the full
viscous congruence, followed by projection on caller-supplied expanded physical
modes. It preserves every modal cross term and adds no drilling loss, modal
decay fit or second viscosity condensation. This is a fixed condensed trial
space; dynamic bubble-coordinate elimination is not claimed.

`projection_terms` charges conservative checked condensation and complete
modal projection work before projection allocation. Existing mesh budgets and
director/section admission apply. Cancellation is polled per element, quadrature
point and mode row, with no partial returned matrix. Focused regression sources
check proportional elastic/viscous material work, full projection including cross
terms, explicit zero, indefinite tensors, budgets and cancellation; source
presence alone is not a native pass.

`shell::survey::{validate_midsurface, material_axis_angle}` expose the SAME
existing manifold/fan and global material-direction owners used by `MeshShell`.
A new material consumer can validate its actual midsurface without fictitious
thickness/material fields, and project grain independently of triangle-edge
ordering. This retains the existing bounds on duplicate vertices/facets,
connectedness, winding, vertex fans and unresolved tangent projections; it does
not add a global self-intersection certificate. `ShearShellSection::with_material_angle_rad`
preserves the physical section while selecting the resulting facet angle.

`shell::pressure::ShellPressurePort` selects original oriented midsurface
triangles and supplies one P1 unit-pressure load/normal-volume-flow operator.
The legacy DKT `ShellReduction::volume_port` and material MITC consumers share
that integral. Empty selections retain the complete mechanical coordinate
extent with zero fluid participation; omitted internal facets do not create
fictitious enclosed volume. Duplicate selection, out-of-range geometry,
nonfinite modes, projection budgets and cancellation refuse. G0/G3 manufactured
tests check closed-surface rigid motion, affine volume, pressure work, signed
orientation and exclusion of internal structure. This owner does not extrude
skins, infer fluid ownership, certify closure or supply radiation impedance.

`ShellReduction::open_screen_radiation` supplies an explicit zero-thickness
midsurface chart in the original mass-normalized source basis. Every facet
corner retains `n_f dot translation_mode,node`; adjacent folded facets keep
their distinct normals and prescribed free-edge velocities. Physical rotations
have zero acoustic offset in this chart. Thickness still determines source
mass, stiffness and damping; the separate extruded-skin chart is unchanged.
Optional conforming midpoint refinement shares vertices topologically and
interpolates each original facet's P1 trace exactly. It adds acoustic resolution
without adding structural modes, smoothing the geometry, grading the edge or
certifying source convergence. Final vertex/facet/complete-corner limits and
bounded cancellation polls apply before publication. Open topology and frequency
admission remain with the pressure-jump owner. `hessian_matrix_with_cancel`
shares the original complete Hessian arithmetic, with facet/row cancellation;
the existing unchecked scratch API retains its numerical order. Two source
tests cover folded traces, exact P1 work on three refinement rungs, topology
sharing, complete budgets and cancellation/replay; a tangent test pins arithmetic
parity. These are authored regression sources, not executed numerical acceptance.


`shell::arching` (IRP M1.2-W6) admits explicit cubic-Hermite section curves or
a sampled bicubic height grid, with a required unloaded/declared-loaded
reference state. Four-point nodal differentiation plus tensor Hermite
interpolation gives C1 height/gradient interpolation and reproduces bicubic
polynomials. Lofted cross sections must meet the long arch at x=0; no missing
section slopes or geometry parameters are defaulted. Graph lifting preserves
vertex order, triangles and holes and supplies analytic unit graph directors;
it does not unload a scan or certify arbitrary self-intersections. Graduation
is an independent `ThicknessField`; nonpositive samples refuse. Lifts and
graduation sampling have explicit work caps and node-boundary cancellation.
Meshing stays in L2 `fs-mesh` and is a dev dependency here, outside the runtime
WASM graph. The new `arched_outline_shell` target composes outline admission,
lift, physical graduation, `MeshShell` admission, DKT/MITC modal ladders and a
Bessel-shape-paired annular-film oracle. These tests require native execution
before any convergence or instrument-geometry qualification claim.

`ArchingParameters::NurbsHeightMap` (M6.1-W2) shares an existing immutable
`fs-rep-nurbs::NurbsSurface<f64>` with the same lift and graduation paths.
Physical x/y graph parameters are required: query positions and the planar
Jacobian must reproduce x/y and identity (16 eps positional, 1e-12 derivative
tolerances). General parametric charts refuse rather than silently changing
coordinates. Derivatives come from the representation owner; this adapter
does not fit, resample, smooth or duplicate the measured surface. The lift
charges twenty control-net visits per query (four per evaluation; a seam
corner takes four one-sided partials and the exact point) and retains its
existing caps and cancellation boundaries. Degree-one graphs have C0 cell
seams, where the ordinary first partial is undefined. A node on a seam keeps
the exact surface height, and its director uses the mean of the one-sided
slopes, each evaluated 2^-26 of the shorter adjacent span inside that span.
The bicubic C1 claim above applies only to the bicubic path.
`fs-scenario::assembly_generators::ag_011_arching_second_order_and_existing_mesh_lift`
tests the measured outline → existing mesher → shared NURBS → existing lift.

With the optional `scenario` feature, `ArchingSurface::from_measured` and
`ThicknessField::from_measured` lower fs-scenario's measured `ArchingSurface`
generator, the one owner of measured arching and graduation. They take its
height and thickness graphs and its declared reference state: `Unloaded`, or
`Loaded` described by the source label and state. The caller declares nothing,
and a source without a thickness table refuses. The feature is optional, so
fs-plate alone does not build fs-scenario's graph. The section and bicubic charts above carry no
source dimensions, bounds or ledger rows. They are declared-parameter charts
(M1.2-W6), so measured geometry with provenance arrives only through this
lowering.

`shell::junction` (IRP M1.2-W2a) couples smooth `ShellComponent` charts along
1-D edge polylines (`EdgeTie`) or 2-D overlapping patches (`FaceTie`) without
averaging nodal directors across fold seams. On each 1-D seam, the Wohlmuth
(2000) biorthogonal linear dual basis `psi_0(xi) = 2 - 3*xi`,
`psi_1(xi) = 3*xi - 1` produces a diagonal slave mortar mass `D` and local
transfer operator `C = D^{-1} M` with a minimum-norm local affine consistency
projection onto the master polyline's affine hull so that
`sum_j C_{ij} = 1` and `sum_j C_{ij} p_{m,j} = p_{s,i}` hold to machine
precision on both straight and curved seams. `TieOffset::Midsurface` ties
midsurfaces directly; `TieOffset::Skin { side }` (`side` in `{-1, +1}`)
evaluates the rigid attachment arm `d = side * 0.5 * h_master * n_master` from
explicit component nodal thickness and area-weighted nodal normals, enforcing
`u_slave = u_master + theta_master x d`. `RotationCoupling::Full` eliminates
all 6 slave DOFs per node; `RotationCoupling::TranslationsOnly` eliminates only
translations, leaving slave rotations free. `RotationCoupling::SeamTangent`
eliminates the translations and only the rotation about the seam tangent
(taken from the master polyline through `C`). That is the fold kinematics: at
a fold, `Full` would tie one side's bending slope to the other side's drilling
spin, so the drilling penalty would enter the junction. On a T-junction across
alpha 1e-4..1e-2, `Full` moves the first ten frequencies by 6.6e-4 (CST/DKT)
and 2.4e-3 (MITC3+), and `SeamTangent` by 3.2e-6 and 1.4e-6. `Full` stays
the coplanar-seam coupling and reproduces a conforming monolithic mesh, which
shares all six DOFs. A `FaceTie`, or a compliant glue tie, refuses
`SeamTangent`. Congruence reduction
`K_r = T^T K T`, `M_r = T^T M T` preserves exact symmetry, positive-semidefinite
strain energy, and zero net virtual interface power `P_m + P_s = 0`.
`TieBudget` enforces explicit caps on `max_segments`, `max_slave_dofs`,
`gap_tolerance_m`, and `angle_tolerance_rad`, returning an auditable
`TieReceipt`; conflicting cross-point constraints and duplicate slave edge spans
refuse with `PlateError::BadConstraint`.

`shell::junction::{GlueLine, TieKind::Compliant}` and `shell::ownership`
(IRP M1.2-W2b) extend multi-component shell assemblies to closed 3-D boxes and
Kelvin–Voigt compliant adhesive seams. At 3-way closed-box corners where three
`EdgeTie` polylines meet, Wohlmuth (2000, §2.3) wirebasket corner pinning and
transitive slave resolution eliminate corner slave DOFs against the shared
master corner root without double-constraining the third edge. `TieKind::Compliant(GlueLine)`
retains both sides' DOFs in the reduced basis and adds a symmetric positive-semidefinite
dual-mortar interface stiffness `glue_stiffness` (`N/m^3 * width_m` integrated
against the Wohlmuth linear dual basis along 1-D `EdgeTie`s or exact 2-D triangle
patches on `FaceTie`s) and optional viscous damping `glue_damping` (`N*s/m^3 * width_m`),
preserving all 6 rigid-body modes (`1/2 q_r^T K_glue q_r <= 1e-14`) and converging
at `O(1/k)` to the `TieKind::Rigid` spectrum as glue stiffness scales by `10^k`.
`FluidBoundaryOwnership::audit` verifies complete single-owner fluid-domain
attribution across all `(component, facet, side)` skin faces, classifies
`FaceTie` and tied-edge glue-overlap strips as `SolidInterior`, checks manifold
orientation compatibility across every tied edge, refuses untied cavity boundary
edges that are not covered by a declared `ApertureLoop`, and computes enclosed
cavity volume `1/6 sum_f (x_0 . (x_1 x x_2))` on the inner skin surface closed by
tied-edge bevel and corner triangles.

`shell::mitc::{MitcShellModel::strain_energy, interpolate_nested_shell, s_norm_error}`
(IRP M1.2-W1b) evaluate the elastic strain energy `1/2 u^T K u` and the Hiller &
Bathe (2003) / Lee, Lee & Bathe (2014, Eq. 20) relative s-norm error
`E_h = ||u_ref - u_h||_s^2 / ||u_ref||_s^2` on nested refinements (`N_ref > N_max`,
`N_ref % N_c == 0`, with coincident coarse-to-reference vertex verification).
`tests/shell_vibration_locking.rs` validates the MITC3+ shell against Leissa
(1973, *J. Sound Vib.* 31:257) thin rectangular plates (SSSS, CCCC, CSCS, CFFF,
SSFF), closed-form Mindlin (1951) thick plates (`t/a in {0.01, 0.1, 0.2}` with
`kappa = 5/6`), closed-form Flügge (Leissa 1973, NASA SP-288 ch. 2) shear-diaphragm
circular cylinders (`R/t = 100`), the Lee & Bathe (2014) / Chapelle & Bathe (2011)
asymptotic locking suite (`t/L in {1e-2, 1e-3, 1e-4}` on cylinder, hyperboloid,
and hyperbolic paraboloid in membrane- and bending-dominated regimes), and reports
orthotropic wood-like (`E_L / G_LT in {15, 20}`) transverse-shear frequency shifts
against CST/DKT.

- `PlateSection::{orthotropic, isotropic}` — bending rigidity D from
  engineering constants: `D11 = E1h³/12(1−ν12ν21)` etc., `D33 = G12h³/12`.
  Material axis 1 is the grain (L) axis — the matdb axis-convention
  contract. Raw constructors take coherent SI (Pa, m, kg/m³); the `_qty`
  variants, `Stiffener::qty`, and `AssemblyOptions::qty` are the
  dimensioned fs-qty front doors (Length/Density/Pressure/Area/
  SecondMomentOfArea/SurfaceTension), bit-identical to the raw paths and
  pinned against a hand-computed spruce D-matrix (D11 = 27.2503 N·m,
  D22 = 2.04378, D12 = 0.715322, D33 = 1.6875 exactly). ν stays a bare
  ratio (dimensionless); the elastic law keeps fs-material's raw-Pa
  contract.
- `PlateSection::orthotropic_plane_stress` consumes only E1, E2, nu12 and
  G12 plus thickness/density. It uses reciprocal nu21 = nu12 E2/E1 and
  requires positive reduced compliance, not an isotropic Poisson bound.
  No out-of-plane constants are synthesized. `orthotropic` delegates to
  the same calculation; nonfinite inputs and unrepresentable bending
  stiffness refuse. G1 compares the reduced stiffness to an independent
  2-by-2 compliance inversion and exchanges principal axes, including an
  admissible synthetic nu12 > 1/2. The constitutive reduction follows
  [Roylance, Laminated Composite Plates](https://web.mit.edu/course/3/3.11/www/modules/laminates.pdf),
  equation (4) and through-thickness moment integration.
- `PlateSection::orthotropic_plane_stress_at_angle` rotates a uniform material
  section counterclockwise from mesh x toward y using the energy-conjugate
  engineering-curvature transform `D_mesh = T^T D_material T`. The full D16/D26
  bending/twisting terms reach the existing DKT assembly. Zero angle preserves
  the aligned section exactly; isotropic constants are angle-independent.
  Angles must be finite radians. G1 compares a rotated section and an irregular
  DKT element to independent Hessian projection/compliance energy; G3 checks
  half-turn periodicity, isotropic invariance and reciprocal quarter-turn
  equivalence. This is a homogeneous single-layer model, without laminate
  extension/bending coupling. Piecewise grain is represented by independently
  rotated element sections, not by interpolation within an element.
- `PlateMesh::{rectangle, rectangle_boundary, from_unstructured, structured_equivalent, boundary_nodes, boundary_edges, total_area}` —
  structured and unstructured right-triangle and general 2D triangular meshes;
  `from_unstructured` validates quality, coordinate finiteness, connectivity, and positive
  signed area (`2A > 1e-15`). `structured_equivalent` reproduces `rectangle` bit-identically.
- `PlateChart` and `PlateRegion` — mid-surface soundboard and plate chart representations
  binding geometry, thickness, material section, sub-regions, and boundary supports.
  `with_element_sections` replaces the uniform section with one validated section
  per triangle. `section_field` reports the active uniform/per-element field.
  Region labels alone do not alter constitutive behavior; their triangle indices
  can select entries of the section field. Changing mesh topology requires a
  matching field; stale counts refuse at assembly.
- `PlateSectionField` and `assemble_with_sections` use each triangle's full D,
  thickness and density for both stiffness and translational/rotary mass in the
  existing DKT assembler. Sections require finite positive thickness/density and
  symmetric positive-definite D. The uniform assembler delegates to the same
  kernel. Mesh admission is rechecked before indexing mutable chart data.
  Interfaces coincide with triangle edges and share displacement/slope DOFs:
  perfect bonding is assumed. Thickness is centered on the common mid-surface;
  offsets, laminate coupling, interface slip and delamination are not modeled.
- `triangulate_soundboard` — parametric guitar/instrument soundboard mesh generator.
- `dkt_stiffness` — the 9×9 DKT bending stiffness. DOF convention per node:
  `(w, wx, wy)` with wx = ∂w/∂x, wy = ∂w/∂y (SLOPES; Batoz's published
  θ-tables are translated internally, producing βx = −w,x, βy = −w,y whose
  sign cancels in the curvature quadratic form). Exact three-mid-side-point
  integration of the quadratic integrand.
- `assemble(mesh, section, boundary, stiffeners, opts) -> PlateModel` —
  reduced (K, M) pencil with Dirichlet elimination
  (`EdgeSupport::{SimplySupported, Clamped}`), lumped mass (translational
  ρhA/3, rotary ρh³/12·A/3 — SPD by construction, as the fs-modal pencil
  contract requires), optional uniform membrane pre-tension T entering as
  the P1 geometric stiffness on w, and `Stiffener` Hermite beams: bending
  with the parallel-axis effective rigidity `EI + EAe²` on (w, slope-along),
  torsion GJ on the cross-slope, lumped translational beam mass.
- `PlateChart::assemble_with_mass(stiffeners, opts, transverse_mass, stiffener_mass)`
  selects the existing `TransverseMass::{Lumped, Linear, EdgeCubic}` panel law
  independently of `StiffenerMass::{Lumped, ConsistentHermite, ConsistentEccentric}`.
  `ConsistentHermite`
  exactly integrates `rho*A*N(s)^T*N(s)` for the same cubic Hermite transverse
  displacement used by beam bending. Its coordinates `(w1, dw/ds1, w2, dw/ds2)`
  map through the existing beam tangent into plate `(w, wx, wy)` DOFs, including
  support elimination. It replaces endpoint lumping; it never adds a second
  copy of beam mass. Geometry, total physical mass, beam stiffness and all
  existing assembly methods remain unchanged. Zero density admits an explicitly
  massless brace; negative/nonfinite density and overflowing consistent inertia
  refuse. This is translational Euler–Bernoulli inertia only: no eccentric
  axial motion, rotary inertia, cross-slope torsional inertia or shear deformation
  is inferred. The existing parallel-axis rigidity remains `EI + EAe²`.
- `StiffenerMass::ConsistentEccentric` adds the exact beam kinetic terms
  `rho*(I+A*e²)*integral((d/ds w_dot)²)` and
  `rho*A*e²*integral(beta_cross_dot²)` to that Hermite translation law.
  The first uses the derivative of the same cubic Hermite displacement;
  the cross slope uses the existing linear torsion interpolation and normal
  `(-ty,tx)`. With physical axial rotations `(wy,-wx,0)` and centroid arm
  `(0,0,e)`, both `A*e²` terms are translation of the offset beam centroid.
  `I` supplies only the given centroidal bending rotary inertia. The stiffness,
  total material mass and existing mass variants are unchanged; constant
  translation still has exactly `rho*A*L` inertia. Saint-Venant `J` is not a
  polar mass moment: no torsional section rotary inertia, independent axial
  DOFs or beam shear deformation is inferred. Missing/nonfinite geometric
  data and overflowing derived inertia refuse rather than dropping a term.
- `modes(model, window, opts)` — thin front over `fs_modal::slice_window`:
  every frequency arrives as a certified eigenvalue interval and the
  in-window count is inertia-certified.
- `PlateError` — typed refusals with stable `FS-PLATE-*` codes: bad
  section, section count (expected and actual), degenerate element (with the offending element id and 2A), bad
  boundary (with out-of-bounds node index and node count), bad stiffener,
  forwarded modal refusals.
- `loading::PlatePointStencil` locates a point on the actual triangular mesh
  within node/triangle admission bounds. It uses the lowest-index containing
  triangle at shared edges and P1 weights for the `(w, wx, wy)` interface
  fields. This interpolation is explicit and distinct from the DKT element's
  interior rotation polynomial. A physical transverse force and Cartesian
  moments `[Fz, Mx, My]` map to `[Fz, -My, Mx]` generalized efforts: physical
  rotations are `[wy, -wx, 0]`. Applying the same stencil and its transpose
  preserves virtual work and power. Point probes reconstruct displacement,
  slopes, transverse velocity and angular velocity; `load_vjp` differentiates
  force/moment amplitudes at a fixed location. Repeated additive transfers
  support caller-bounded distributed loads. Eliminated-DOF applied loads are
  returned separately; they are not the complete dynamic support reactions.
  Invalid/outside locations, malformed support maps and nonfinite arithmetic
  refuse, and a refused addition leaves the force vector unchanged. Moving
  paths relocate on the undeformed mesh; changing supports requires a new
  stencil. No path-coordinate derivative, contact model or moving-support
  kinematics is inferred. `tests/plate_loading.rs` checks independent rigid
  motion, force/moment resultants, virtual work, moving-edge continuity,
  support elimination and bounded atomic refusal.
- In the `fs-plate-transient` crate, `PlateDynamics` binds the existing sparse DKT
  pencil to `fs-time`'s structural residual, discrete adjoint and loaded
  trajectory interfaces. It borrows a fixed mesh/support model and uses
  `M = mass_scale*M0`, `K = stiffness_scale*K0`, and
  `C = mass_damping_per_s*M + stiffness_damping_s*K + C_material`.
  The optional material matrix is assembled by `viscosity::PlateBendingViscosity`
  and remains fixed under the four mechanical parameters. The four parameters
  precede the load's parameters in every VJP. Explicit stored transposes
  and all mass/stiffness/damping derivatives are implemented without dense
  matrices or per-action scratch. `PlateLoad` supplies a pure time-dependent
  reduced load and its parameter VJP, permitting moving mesh-point stencils.
  The supplied mass may be lumped or consistent, including displacement/slope
  and stiffener cross inertia. Diagonal mass keeps linear-time admission;
  consistent mass is checked for symmetry to assembly roundoff and positive
  sparse LDLT inertia once, without changing or lumping the physical matrix.
  `new_with_mass_options` exposes the original factor's numeric byte/pivot
  controls; cancellation is checked around that bounded synchronous work.
  Those controls also bound sparse passivity admission of the actual Rayleigh
  combination when stiffness-proportional damping is selected. Prestress may
  make K indefinite while mass damping still makes C passive; exact zero-loss
  directions remain unchanged. See `fs-plate-transient/CONTRACT.md` for the
  dimensionless arithmetic allowance and refusal boundary.
  The mass Jacobi preconditioner is exact only for lumped mass; consistent
  acceleration and its initialization pullback solve the full original matrix
  under explicit Krylov restart/cycle limits. An
  effective-operator Jacobi preconditioner supports the explicit adjoint
  solve. Forward Newton/FGMRES automatically uses the same physical diagonal
  through the model hook, using the integrator's actual effective mass,
  damping and tangent weights. This action allocates no scratch and refuses
  nonpositive, nonfinite or noninvertible diagonals before a step is accepted.
  DOF, full support-map, stored-entry and load-parameter caps precede scanning
  or allocation. Admission/preconditioner construction polls cancellation;
  operator/load callback work is bounded by the admitted model and caller.
  Recording, solve scratch, replay, checkpoint and observation budgets remain
  owned by the shared time integrator. Geometry, supports and reference
  matrices are fixed; stiffness scaling also scales any prestress/stiffeners.
  The adapter does not independently certify positive-semidefinite K or
  identify Rayleigh coefficients from physical damping. It retains the supplied DKT/inertia spatial model,
  linear small-deflection kinematics, and the supplied initial-state chain.
  `examples/moving_plate.rs` advances a force across a clamped DKT mesh and
  differentiates sampled deflection through all four mechanical parameters
  and force amplitude, including consistent initial acceleration. The focused
  transient tests compare complete q/v/a and gradients with independent dense
  stepping, check common-scale invariance, undamped energy conservation and
  damped decay, and exercise admission/cancellation/refusal behavior.

### Consistent CST/DKT shell inertia

`shell::CstDktMass::Consistent` integrates the existing P1 midsurface velocity
and P1 physical director rotations over each flat triangular facet and its
physical thickness. It retains inter-node mass terms and rotates the
through-thickness transverse inertia into global axes. This is a declared
kinetic interpolation on the existing CST/DKT stiffness, not an exact dynamic
interpolation of the DKT interior or a replacement for the MITC shell model.

All six rigid velocities reproduce the exact mass and spatial inertia of the
represented triangular prism or variable-thickness wedge. For uniform thickness,
the original vertex-lumped triangle has four times its continuum centroidal
midsurface orbital inertia. Physical
normal spin is already carried by the midsurface translations; it receives no
additional fictitious through-thickness rotational inertia in the consistent
model. The drilling coordinate remains numerically invertible through a
kinetic gauge on `theta_normal_dot - membrane_spin_dot`, using the original
`rho*h^3/12` scale. That term vanishes for every rigid motion. It is numerical
stabilization and is not a material micro-rotation law.

`MeshShell::assemble` consumes consistent inertia by default;
`MeshShell::assemble_with_mass` admits the original lumped comparison. The
survey's measured nodal thickness is P1 within each facet, symmetric about its
reference midsurface, with the supplied piecewise constant material density.
A positive seven-point degree-five triangle rule integrates `rho*h*Ni*Nj`
and `rho*h^3*Ni*Nj/12`, preserving the thickness-dependent center of mass,
orbital inertia and transverse director storage. The relative drilling gauge
uses that same local `rho*h^3/12`. Constant-thickness facets retain the original
homogeneous arithmetic exactly; stiffness retains its facet-average section.
The `assemble_shell[_sections|_materials]_with_mass` entry points expose the same
mass choice for homogeneous section data. Existing assembly entry points, generated
profile defaults, and the frozen bonded piano pencil retain their original
arithmetic. Stiffness, support maps and original signed shell-tie/offset
congruences are shared. No laminate first mass moment is guessed from a
homogeneous section. Tests exercise the actual survey consumer against an
independent triangular-prism and variable-thickness spatial inertia,
factorial-moment affine velocity/director integration through degree five,
frame changes, three subdivision levels, sparse positive-definite supported
mass, and explicit legacy profile equivalence.

### Graded reference density in CST/DKT shells

`assemble_shell[_sections|_materials]_with_reference_density` and
`MeshShell::assemble_with_reference_density` accept three positive densities
[kg/m³] in each original facet's corner order. The P1 field replaces only the
section's scalar reference density; neighboring facets may disagree at a shared
vertex, so material interfaces are not averaged together. The existing elastic
sections, support map and consistent translational/director mass remain the
owners of stiffness, dynamics and modal normalization. These entry points are
explicitly consistent-only; legacy lumped APIs retain their original behavior.

Measured shells integrate P1 density and P1 thickness together. Sixteen positive
Gauss–Duffy points integrate `rho*h*Ni*Nj` and `rho*h^3*Ni*Nj/12` through degree
six, including the existing relative drilling-rate gauge, which still vanishes
for all six rigid motions. Uniform density delegates to the original mass
kernel and retains its arithmetic, including measured nonuniform thickness.
`MeshShell::mass_with_reference_density_kg` integrates the same reference
inventory analytically, including density/thickness correlation. The survey's
original `mass_kg` continues to describe its original piecewise constant
materials; the explicit field does not mutate that specimen or invent a source
measurement, live thermal density update or mass-transfer law.

`tests/shell_reference_density.rs` checks independent complete spatial inertia,
degree-six director storage, equal-mass redistribution, sharp interface and
uniform-limit equivalence, malformed fields, density-scaled native elastic
eigenvalues and the existing mass-normalized acoustic pressure port. No new
eigensolver, acoustic model or time integrator is introduced; those consumers
receive the resulting physical mass matrix.

### Graded panels with bonded beam reinforcement

`assemble_stiffened_shell_with_reference_density` and its `materials` form
bind one original density triple per panel facet and one endpoint pair per
bonded `ShellBeam`. Both fields are admitted before shell assembly. Adjacent
panels or reinforcements can retain distinct densities at a shared node;
neither field is replaced by its mean or inferred from the other component.

The original consistent panel reference-mass kernel and existing selected
`BeamMass` kernel integrate their respective fields. The same eccentric
rigid-arm congruence transforms beam stiffness and inertia before sparse
accumulation and original support restriction. Beam translation, sectional
spin and offset coupling therefore retain the same density distribution.
The resulting native pencil supplies modal normalization and frequency;
no additional integration formula, mass owner or elastic law is introduced.
Every original bonded assembly entry retains its frozen legacy pencil.

Five regression definitions in `shell/stiffened/reference_density_tests.rs`
check independent complete spatial moments and frame covariance, polynomial
density--velocity correlation for both beam mass images, supported native
modes and density scaling, original uniform/empty/material-section paths,
and complete-field admission. Changed Rust was parsed with pinned Rustfmt
and independent moment arithmetic was checked; native execution of these
new regressions remains unperformed.

### Layered CST/DKT reference density and first mass moment

`ReferenceDensityMoments` and the `*_with_reference_density_moments` assembly
entry points extend that same mass owner to asymmetric coatings and layered
inventories. The immutable carrier stores `qk = integral rho*s^k ds`, all in
kg/m³, where `s = z/h` follows the mesh normal and its declared reference plane.
`from_laminate` consumes the existing material owner's `I0/I1/I2` without
reintegrating or replacing its layers. Local thickness scales the source's
layer fractions and reference offset. A positive-definite moment block is
required: `q0 > 0`, `q2 > 0`, and `q0*q2 > q1²`; `q1` remains signed.

The physical moments `h*q0`, `h²*q1`, and `h³*q2` are integrated with the same
positive degree-six rule and measured P1 thickness. Kinetic velocity is
`u_dot + z*(theta_dot cross normal)`, so the first moment supplies the signed
translation/director mass blocks. It is not a stiffness, damping, or added
normal-spin inertia. The original relative drilling gauge retains the second
moment scale and vanishes for rigid motion. Symmetric uniform-through-thickness
inputs delegate to the preceding scalar/P1-density kernel exactly.

These APIs change inertia only. Independent CST/DKT elastic sections remain
explicit; full laminate A/B/D constitutive coupling and transverse shear use
the existing MITC `ShearShellSection::from_laminate` path. The mesh and moments
must describe the same reference plane and normal; reversing a normal reverses
the first moment for an unchanged physical inventory. Tests in
`shell_layered_density.rs` cover independent full spatial moments, graded
density/thickness correlation, layer reversal, reference-plane shifts, normal
orientation, homogeneous limits, and the existing modal/acoustic consumers.

### Three-dimensional beam inertia

`shell::stiffened::timoshenko::{segment_matrices, assemble_beam_chain}` and
their existing `tie_beams_to_shells` consumer now use consistent beam mass.
One shared kinetic integral covers axial motion, torsion, both bending planes,
and section spin. It integrates the exact static Timoshenko displacement and
independent rotation fields with four Gauss points (degree-six exact); zero
shear compliance gives cubic Hermite displacement and its Rayleigh rotary
inertia. Axial/torsional fields are linear. The polar torsional inertia is
`rho*(Iy+Iz)`, independent of the Saint-Venant stiffness constant `J`.

Rigid motion has the continuum mass and spatial inertia for any segment
length: transverse orbital inertia about the center is `rho*A*L^3/12`,
plus `rho*I*L` for section spin. Endpoint lumping instead gives
`rho*A*L^3/4` for the orbital term. The default chain and shell-tie path no
longer incurs that single-segment error. General deformational frequencies
still require spatial refinement; the static shape functions do not make a
finite element exact for every dynamic wavelength.

`BeamMass` makes the discretization explicit on the `*_with_mass` variants.
`assemble_stiffened_shell` and its material-section form retain their frozen
legacy lumped piano pencil and bit-semantics version. Their `*_with_mass`
forms admit consistent mass; the same signed eccentric rigid-arm congruence
transforms stiffness and inertia, including translation/rotation coupling.
The shell's own mass choice is not altered. Prismatic centered sections,
small motions, fixed section orientation and straight segments remain the
model assumptions; no warping inertia or finite-rotation gyroscopic law is
implied. Tests cover independent polynomial kinetic integrals, arbitrary
spatial rigid motion with eccentric bonds, the shear-rigid limit, the
continuum Timoshenko spectrum, bonded/tied congruence, and the unchanged
legacy piano golden.

### Graded reference density in three-dimensional beams and shell ties

The `timoshenko::*_with_reference_density` segment and chain entry points
retain two original densities [kg/m³] per segment, in its node order. Density
is P1 along the original prismatic member and uniform across each section.
Each segment owns both endpoint values; shared-node material jumps are not
averaged. The field replaces the section's scalar density in axial, torsional,
both transverse and both section-rotation kinetic terms. Geometry, elastic
coefficients and shear areas retain their original owners and admission.

Consistent axial and torsional blocks use exact P1-density/P1-velocity moments.
The existing four positive Gauss points integrate the complete degree-seven
graded Timoshenko kinetic polynomial without extra quadrature. Section rotation
remains independent of centerline slope; the zero-shear-compliance limit retains
the Rayleigh/Euler--Bernoulli field. Polar torsional inertia still uses `Iy+Iz`,
never the Saint-Venant stiffness constant. Explicit lumping integrates each
`rho*N_i` partition, preserving total inventory and its first axial moment.
Uniform fields delegate to the scalar kernel with the supplied density and
preserve its arithmetic, including the original default and legacy paths.

`junction::tie_beams_to_shells_with_reference_density` transfers this complete
mass through the original nonmatching mortar and signed eccentric rigid-arm
map. Every beam and segment needs an explicit field. Its prestressed companion
and the corresponding chain/segment functions retain the independently supplied
axial forces in stiffness only. Shell inertia, existing supports and shell ties
are unchanged. The returned matrices retain the original inventory; resampling
temperature-dependent current density on frozen geometry is not mass transfer.

`tests/beam_reference_density.rs` checks independently integrated full spatial
mass moments, reversed members, static-polynomial kinetic energy in every
velocity field, finite shear and Rayleigh limits, positive lumped partitions,
uniform-path bits, discontinuous native chain assembly, and off-node shell-tie
mass/centroid/rotary work with prestress. The model remains a small-motion,
prismatic centered beam; this field does not introduce transverse eccentric
density, graded elastic coefficients, warping inertia or moving-frame transport.

### Signed axial prestress in beams and shell ties

`timoshenko::{prestressed_segment_matrices, assemble_prestressed_beam_chain}`
accept the caller's equilibrated axial section forces in newtons, positive in
tension and negative in compression. The original elastic element gains the
conservative total-slope initial-stress potential
`N/2 * integral(v'^2 + w'^2) dx`. Three Gauss points integrate its quartic
strain products exactly on the SAME Euler--Bernoulli or Timoshenko displacement
fields. The independent section rotation is not substituted for displacement
slope. Both bending planes and their signed rotation conventions are retained.
This is the total-slope beam-column model of Cheng and Pantelides (1988),
doi:10.1061/(ASCE)0733-9445(1988)114:5(1152).

`junction::tie_prestressed_beams_to_shells` transfers this complete tangent by
the original nonmatching mortar and eccentric rigid-arm map, including
translation/rotation cross terms and supported/slave coordinates. Every beam
and segment needs an explicit finite force. Mass, shell stiffness and existing
ties are unchanged. Zero-force fields reproduce the original segment, chain
and beam-tie pencils bit for bit. Compression-induced indefinite tangents are
retained for the existing stability owner, never clipped or regularized.
Nonfinite element or accumulated chain/tie entries refuse transactionally;
the tie retains its original work budgets and cancellation boundaries.

G0/G1/G3 tests cover independent continuum polynomial work, the closed-form
Euler geometric matrix, rotated two-plane work, zero-force/mass parity, signed
eccentric mortar virtual work, and finite-input/overflow/cancellation refusals.
Three mesh rungs bracket the pinned critical load against
`Pcr = EI*q^2/(1 + EI*q^2/(G*As))` (Euler limit `EI*q^2`), within 3e-5 relative
on the finest rung. Four rungs check the continuum bending/shear frequency
Rayleigh quotient with original rotary inertia to 1e-5 relative. The static
trial fields still require spatial refinement for a loaded dynamic problem.

This extension requires an already balanced, straight-segment reference and
physical supports. It does not solve prestress equilibrium, infer external
forces or shell membrane stress, change the nonlinear bonded-beam preload
model, or claim follower-load, torsional-buckling, finite-rotation or
post-buckling dynamics. The broader assembled-equilibrium bead remains open.

## Invariants
1. Element certificates (tested on an irregular triangle): stiffness
   symmetry; all three rigid-body motions store zero energy to roundoff;
   a constant-curvature field reproduces `A·κᵀDκ` EXACTLY — the patch test
   that pins the shape-function tables and the DOF convention (it caught a
   θ↔slope translation error during development).
2. Convergence to independent analytic references at two mesh densities
   with an O(h²)-or-better trend and ≤2% fine-mesh error: simply-supported
   Navier (measured 0.72% → 0.18%), clamped Leissa 35.992 (1.02% → 0.28%),
   orthotropic Navier rectangle (0.24%), continuum membrane in the
   prestress-dominated drumhead limit (0.41% → 0.10%).
3. Stiffener terms act at first-order Rayleigh scale (tested by term
   isolation — the gate that caught an asymmetric Hermite sign error which
   had inflated the fundamental 64%): weak-brace bending ≤ +1%, torsion
   ≤ +1%, lumped mass ≤ −1% on the quarter-line fixture; eccentricity
   multiplies a weak brace's rigidity visibly (measured 60.8 → 76.6 rad/s).
4. Orthotropy is load-bearing: swapping E_L/E_R moves the Navier
   fundamental by 89% on the spruce-like fixture (the wrong-grain-axis
   mutation gate).
4b. Stiffened-panel literature case: the Olson & Hazell (1977) clamped
   203×203×1.37 mm plate with one central 6.35×12.7 mm integral rib
   (E = 68.7 GPa, ν = 0.3, ρ = 2820 kg/m³; reference rows reproduced in
   Srivastava/Datta/Sheikh 2004 Table 4 and Thinh/Binh/Tu 2013 Table 1).
   The first five modes MAC-pair order-preservingly across two mesh
   densities (measured MAC ≥ 0.9997, floor 0.90) and sit within an
   authored 8% of the Olson–Hazell theory column (measured worst 1.7%,
   fine mesh: 712.2/762.9/981.7/997.6/1395.4 Hz vs
   718.1/751.4/997.4/1007.1/1419.8; experiment column recorded in the
   same JSON row, worst 5.2%). The envelope is authored against the
   model boundary: bending-only DKT + full-composite parallel-axis rib.
5. All modal results flow through fs-modal: counts are inertia-certified
   and every eigenvalue carries its M⁻¹-norm residual interval.

## Error model
Structural misuse panics only inside test helpers; the public API refuses
through `PlateError`. Degenerate triangles refuse with the element id;
stiffener paths refuse on out-of-mesh nodes, zero-length segments, or
non-positive constants. Modal-stage refusals (singular window endpoints,
unresolved windows, non-SPD mass) forward the fs-modal typed error.

## Determinism class
Bit-deterministic across repeat runs on one host: sequential assembly
through `fs_sparse::Coo` canonical accumulation; modal path inherits
fs-modal's determinism (tested there). Cross-ISA goldens are not recorded.

## Cancellation behavior
Spatial assembly and modal calls run to completion under their existing
controls. Moving-load location is bounded by admitted node/triangle counts.
The transient adapter polls admission and preconditioner construction; its
operator callbacks are bounded by the admitted sparse model. The shared time
integrator owns cancellation and accepted-prefix behavior during stepping
and reverse replay.

## Unsafe boundary
None. Workspace `unsafe_code = "deny"`.

## Feature flags
None. The shared generalized-alpha dynamics and adjoint adapter moved to the
`fs-plate-transient` crate on 2026-09-27. As an optional feature here, its
`fs-time`/`fs-solver` edge closed the package cycle
`fs-feec -> fs-couple -> fs-plate -> fs-solver -> fs-feec`, and every workspace
cargo command failed. Its source, tests and examples moved there on 2026-09-28.
Moving-load interpolation remains here.

## Conformance tests

`tests/plate_conformance.rs` (bead 3ez8g.13.2 and 3ez8g.3.9; the reimplementation
contract through the public surface): pt-001 the SS isotropic plate
lands on the exact Navier frequency with mesh convergence (0.12%
fine); pt-002 the clamped square lands on Leissa 35.992 (0.28%);
pt-003 the fs-qty front door is BIT-IDENTICAL to the plain door
through assemble + modes; pt-004 stiffener term isolation (a
vanishing rib is a no-op at 3e-10, a real rib stiffens 56%); pt-005
the Olson-Hazell fundamental-region pin (709.6 vs 718.1 Hz theory,
1.2%; the full five-mode gate lives inline; corpus row
acoustic-olson-hazell-1977-mode3); pt-006 refusals by name — this
battery FOUND that an out-of-range boundary node PANICKED inside
assemble; it now refuses as `PlateError::BadBoundary`; pt-007
unstructured mesh admission (`PlateMesh::from_unstructured` /
`structured_equivalent`) reproduces rectangle meshes, boundaries, and
modal eigenvalues BIT-IDENTICALLY; pt-008 soundboard plate chart modes
on a parametric guitar-top geometry. Inline unit modules untouched.
In-crate: element patch tests (symmetry, rigid-body, constant-curvature
exactness); SS Navier and clamped Leissa two-density convergence ladders
with trend assertions; orthotropic Navier + E_L/E_R swap mutation;
drumhead membrane limit; stiffener eccentricity discrimination in the
deliberately weak-brace regime (a stiff brace saturates into a line
support, where extra rigidity is invisible — measured and documented);
stiffener term isolation against Rayleigh scale; the Olson–Hazell
stiffened-panel literature case with a MAC mode-pairing table (two mesh
densities, order-preserving pairing gate); named refusals (degenerate
element, bad section, bad boundary, bad stiffener). All JSON-line evidence rows; all
modal counts inertia-certified through fs-modal.
Regional G1 integrates a manufactured constant-curvature field across two
sections with different density, thickness and grain on three refinements;
the global quadratic energy and translational/rotary mass match independent
piecewise continuum integrals. G3 requires bit-identical pencils when a uniform
section is repartitioned into equal element sections. G0 rejects wrong counts,
nonphysical D/mass inputs, stale topology and invalid connectivity.
Consistent Hermite beam tests compare every cubic polynomial velocity product
to its independent continuum integral, preserve rigid-translation mass and
subdivided-beam kinetic energy, and reverse segment direction. An assembled
oblique-beam test checks tangential versus cross-slope mapping, rigid tilt,
cubic bending, support elimination and unchanged stiffness/default pencils
with all three panel inertia choices. Eccentric-inertia tests integrate cubic
bending derivatives and linear cross slopes against independent polynomial
integrals; rigid tilt checks the supplied centroidal and parallel-axis moments.
They retain constant-translation mass, reverse directions/offset signs,
subdivide spans, exercise the centered/massless limits and refuse bad inputs.
The assembled oblique-beam case checks both physical slope maps, supports and
unchanged stiffness; changing Saint-Venant J cannot change the mass matrix.
These are discretization checks, not a
measured piano bridge-mobility or high-band convergence claim.

## No-claim boundaries
- The planar `PlateChart` is bending-only; the spatial `shell` owner has
  membrane/in-plane DOFs and relative-spin stabilization. Its CST/DKT model
  is a Kirchhoff model-form rung; the separate MITC3+ path supplies transverse
  shear. Neither source presence nor a pointwise obstacle-course band proves
  arbitrary shell geometry, instrument acoustics or measured material fidelity.
- MITC4 quads are NOT implemented (the new shear path is triangular); trigger-deferred
  at the kj3s0 close (lands with the first quad-mesh consumer or a DKT
  validation gap). Arbitrary 2D planar triangulations are admitted via
  `PlateMesh::from_unstructured` (bead 3ez8g.3.9); 3D CAD surface/solid meshes
  are flattened/extracted into planar mid-surfaces before plate admission.
- The stiffener is a straight Hermite beam on existing plate nodes with
  uniform cross-section: no variable cross-section (vibraphone undercuts),
  no curved braces, no beam shear deformation (Euler-Bernoulli, not
  Timoshenko in v1 — the eccentricity term dominates instrument bracing).
- The eccentricity model is the parallel-axis rigidity `EI + EAe²` acting
  through plate bending curvature; it does not model the membrane force
  the offset beam induces in the plate (requires membrane DOFs). Its optional
  consistent eccentric inertia follows the same existing beam kinematics;
  it does not add the missing in-plane structural load path.
- Lumped transverse and rotary mass remains the default. An explicit
  `PlateChart::assemble_consistent_transverse_mass` option uses the exact P1
  triangle integral for transverse panel inertia while retaining lumped
  slope and stiffener inertia. It preserves total rigid-translation mass;
  torsional beam rotary inertia remains separate. A second
  opt-in, `PlateChart::assemble_edge_cubic_transverse_mass`, integrates a
  cubic Bernstein deflection field whose edges reproduce nodal Hermite values
  and tangential slopes. Its symmetric interior control reproduces quadratic
  fields but is an explicit reconstruction, since DKT does not specify an
  interior deflection field. The existing rotary and stiffener inertia laws
  remain unless `assemble_with_mass` explicitly selects consistent Hermite
  translation or consistent eccentric stiffener inertia.
  `edge_cubic_transverse_shape` evaluates that same field at a supplied
  barycentric point; `edge_cubic_transverse_mean_shape` gives its exact area
  mean. They allow reciprocal point effort/motion and exact volume projection
  without silently returning to P1 after assembly.
  `edge_cubic_transverse_gradient_shape` differentiates that same field
  analytically in physical x/y coordinates. Flat-plate consumers recover
  physical axial rotations as `(w_y,-w_x,0)` and attachment motion as
  `u + theta x arm`; this does not change the DKT stiffness or claim a globally
  C1 field across facets. The field, derivatives and exact mass integral share
  the same Bernstein owner. Analytic tests cover quadratic fields on oblique
  triangles, vertex/edge limits and a nonquadratic cubic interior.
  For this chosen interior
  control, the positive three-point triangle rule samples exactly that area
  mean, although it is not a general cubic integration rule. Cropped panels
  or a quadratic rotation times a varying attachment arm require degree-three
  quadrature within each source facet. Neither option
  alone certifies Model D high-band bridge convergence.
- `PlateChart::assemble_with_kinetic_mass` independently selects transverse,
  `RotaryMass::{Lumped, Dkt}` and stiffener inertia. `Dkt` integrates
  `rho*h^3/12 * integral(beta_x^2 + beta_y^2) dA` using the original quadratic
  Batoz rotations that generate bending stiffness. A fixed positive six-point
  degree-four rule is exact for their products and retains inter-node and
  displacement/slope coupling. This term replaces panel slope lumping directly
  in assembly, without subtracting nearly equal assembled matrices. The result
  is the actual `PlateModel.m` used by modal, harmonic and transient consumers;
  original stiffness, section fields, supports and beam inertia remain intact.
  Existing entry points and `RotaryMass::Lumped` preserve legacy arithmetic.
  DKT rotations differ from gradients of the optional edge-cubic displacement;
  these remain explicitly different kinetic models. This is Rayleigh/Kirchhoff
  rotary inertia, not transverse shear, a solid-wave correction or a physical
  validity certificate at wavelengths comparable with thickness. Focused tests
  check independent edge constraints and integration, rigid/quadratic kinetic
  work, positive mass, exact support restriction and beam retention, and a
  three-rung bending-wavelength ladder for both P1 and cubic transverse fields
  against the declared continuum model.
  Nonfinite or underflowed consistent rotary inertia refuses without mass floors.
- Uniform isotropic pre-tension only (scalar T); tensor/nonuniform
  prestress fields join the soundboard-downbearing consumer.
- Elastic pencils remain (K, M). The explicit MITC viscosity projection above
  supplies a full material C on the same condensed physical trial space;
  time integration and its dissipation ledger remain downstream.
- Validation is against analytic references (Navier, Leissa 35.992,
  membrane) plus the Olson–Hazell stiffened-panel literature case, whose
  experiment column (holography-measured) is recorded but gated only
  through the 8% theory-column envelope. No measured-INSTRUMENT
  (guitar-top holography) comparison yet — that case joins the
  instrument-matdb bead's consuming demo.

## Graded reference density in the native flat-plate pencil

`PlateChart::with_reference_density` retains positive P1 volumetric density
at each triangle's original corners. It overrides section scalar density in
every chart mass assembly method, including the analytic edge-cubic kinetic
image; the section's stiffness/thickness and independent beam inertia retain
their own meanings. Per-triangle corner values can represent material jumps
without averaging neighboring specimens. The field is immutable through its
accessor and coverage is rechecked when assembling the publicly editable mesh.

Density is integrated jointly with the selected velocity fields. Exact
barycentric moments give consistent P1 translation and the degree-seven
edge-cubic translation integral. The latter forms the Bernstein `G*B` once,
then `B^T*(G*B)`, using bounded stack storage and no runtime quadrature.
The original DKT quadratic rotations use a positive degree-five rule; analytic
cubic-gradient rotations use exact weighted polynomial moments. Positive nodal
lumping preserves total mass and first in-plane moments, while retaining its
declared approximation to the quadratic kinetic form. The cubic kinetic owner
removes exactly these same graded panel slope row sums before inserting its
consistent rotation operator. Uniform profiles dispatch the original arithmetic.

The resulting full `PlateModel.m` is the native modal/harmonic/transient mass,
so inertia, inertial forces and mass normalization use one operator. No mean
density, modal retuning or added diagonal floor replaces that reference mass.
This fixed-reference Rayleigh/Kirchhoff plate retains elementwise constant
thickness and symmetric through-thickness density. It does not supply eccentric
laminate moments, evolving density transport, finite-strain thickness change,
or three-dimensional short-wavelength material physics.

`tests/plate_reference_density.rs` checks exact mass/centroid and inertial-load
moments in all kinetic images, quadratic continuum kinetic products, equal-mass
profiles producing different native supported modes, uniform arithmetic parity,
independent beam retention, vertex/frame changes, exact partition invariance,
source retention and malformed coverage. These are G0/G1/G3 numerical checks;
they do not establish measured specimen or audible realism.

## Physical compliant shell attachments

`shell::mitc::connections` adds finite bilateral attachment resultants to the
existing MITC pencil before eigensolution. Global translation and physical moment
rows use `MitcShellModel::project_load`, preserving `k*b*b^T` and
force/displacement work in the existing support frame. The original mass matrix
is preserved. Viscous rows project in the same physical basis without diagonal
fitting, including semidefinite loss and zero retained modal participation.

Explicit connection, footprint and work limits and cancellation apply.
Detached, cancelling, negative or nonfinite laws refuse transactionally.
Regressions cover independent static compliance, oblique support coordinates,
physical moment arms, common rigid motion and viscous work. These additions do
not implement continuum mortar, exact ties, connection mass, finite strain or a
prestressed equilibrium.

### Finite-area load and observation ports

`loading::PlatePatchStencil` composes the existing point maps over explicit
positive physical area quadrature, bounded by `PlatePatchBudget`. It transfers
a total uniformly distributed force/free couple or nonuniform pressure/couple
density, retains the full origin force/moment resultant, and aggregates the
applied load portions on eliminated DOFs. Loading is atomic across the whole
patch. Area-averaged motion and its transpose preserve virtual work and
instantaneous power for a uniform total load; nonuniform traction uses each
quadrature area exactly once. Material velocity is distinct from a moving
observer's convective derivative. The caller owns footprint geometry and
quadrature accuracy; an outside point refuses the patch rather than clipping
its area. There is no inferred contact law, dynamic support reaction, or path
derivative. The loading tests cover footprints spanning triangles/supports,
independent world resultants, work/power, subdivision and late-point rollback.

## Directional shell viscosity projection

MITC membrane, bending and transverse-shear loss admit positive-semidefinite
tensors through the shared material resultant owner, preserving physically
undamped strain directions without regularization. `project_viscosity`
expands the elastic bubble coordinates, evaluates each original quadrature
strain map once, and applies the original constitutive square factors to each
mode's strain before forming the modal Gram matrix. This retains resolved weak
loss that dense nodal assembly followed by projection can erase through
cancellation. Every modal cross term remains; exactly zero factor columns are
skipped without a magnitude cutoff. Temporary channel storage is linear in the
mode count, and projection remains quadratic in that count. The conservative
per-facet work expression is `1_000_000 + 4096*n + 84*n*n`, plus input coverage,
including all seven quadrature points and twelve factor slots. Affine
membrane-energy checks verify active directions, exact null modes, full modal
congruence and relative accuracy of representable weak material loss.

## Shared DKT curvature quadrature

`viscosity::PlateBendingViscosity` assembles complete anisotropic material
bending loss on the original spatial DKT support map, using this same
three-point curvature integration: `C_e = integral B^T V_b B dA`.
`ViscousPlateSection` admits finite symmetric positive-semidefinite resultants,
including exact zero and directional null channels. Its material constructor
retains the selected `PlaneStressViscosity`, applies the original orientation
and thickness-cubed reduction, and checks both endpoints of the caller's
declared resolved frequency band against the source. Cloned element selections
share their immutable source law. This band declaration does not certify a
trajectory's frequency content or extrapolate the material law.

The operator retains all displacement/slope cross terms; it is independent of
the assembled elastic tensor, prestress, stiffeners and inertia. A caller must
supply the original mesh geometry/connectivity used for the retained model.
Explicit triangle, full/free coordinate and pre-elimination entry ceilings
precede allocation; complete element coverage, valid support mapping, finite
assembly and cancellation are checked before publishing. There is no loss
floor, modal approximation, numerical dissipation or history-dependent Prony
model in this Kelvin--Voigt operator. `tests/plate_viscosity.rs` checks exact
constant-curvature power, rigid/null directions, source retention and band
refusals, unchanged damping under prestress/inertia changes, and bounded
admission/cancellation.

`dkt_curvature_quadrature` exposes the same three mid-side curvature matrices
used by `dkt_stiffness`, with their signed doubled area. Each matrix maps the
original transverse displacement and slope DOFs to engineering bending strains;
its integration weight is `twice_area/6`. Stiffness assembly retains its original
arithmetic and ordering. Constitutive consumers can therefore retain actual
material histories without reconstructing curvature from transverse nodal
motion alone or introducing another plate basis. Physical material quadrature
requires positively oriented triangles and retains the original geometry
refusals. This is the existing flat DKT discretization, not a shell or nonlinear
strain model. The coupled full-tensor reed tests exercise the shared curvature
and energy path; native execution of the new tests remains unverified.

## Harmonic local layer-strain observation

Native strain interpolation, elastic director condensation and cached real/
harmonic observations subtract one common nodal translation before weighted
sums. Uniform travel therefore produces exactly zero strain; represented
deformation and independent director bubbles survive unchanged. Full physical
coordinates remain in the owning state and inertia is unchanged. This avoids
amplifying roundoff in nominally zero basis-row sums into strain, stress or
material loss. It does not remove rotations or turn this linear shell into a
finite-rotation model. `rigid_translation_has_exactly_zero_native_and_cached_shell_strain`
and `shell_strain_preserves_deformation_and_independent_bubbles_under_translation`
cover constant/varying directors, real/complex fields and independent bubbles.
Both fail on the original source and pass with translation removed; all ten
existing strain/element tests also pass in the focused native test seam below.

`shell::mitc::strain::PreparedElementStrain` observes the original MITC trial
space with the same director interpolation, elastic bubble condensation and
assumed transverse shear as stiffness and viscosity. It prepares a fixed 48x18
real Bernstein observation map through the original point evaluator. The map
uses 6,912 retained bytes; no physical-mode-sized map or duplicate geometry is
retained. Coherent real/imaginary 18-DOF nodal motion is mapped once per query,
and `HarmonicElementStrain` can serve multiple signed layer intervals. Published
work ceilings are 1,048,576 setup, 8,192 per harmonic field and 4,096 per layer;
these are conservative admission units, not benchmark timings. Only the selected
facet is inspected, avoiding repeated whole-mesh audits.

Six quadratic triangular Bernstein controls represent the complete membrane,
curvature and assumed-shear field. Convexity bounds every point in the triangle
and between the two layer faces. Each control's exact harmonic-cycle Frobenius
peak is obtained from the real/imaginary Mandel Gram matrix, preserving phase.
An isotropic constant-Poisson-ratio plane-stress film includes
`epsilon_zz=-nu/(1-nu)*(epsilon_xx+epsilon_yy)`; engineering shears are divided by
two and tensor shear has weight two in the squared norm. Tensor component order
is [xx,yy,zz,xy,yz,zx]. Shear-energy correction is not applied to strain.

Nonfinite motion, invalid layer faces or Poisson ratio, invalid geometry,
unrepresentable arithmetic and cancellation refuse. Prepared maps are immutable
and remain reusable after refusal. The controlling Bernstein coefficient is not
an observed maximum location; independent component bounds can have different
locations/phases. These are retained linear incremental-field bounds with the
original floating-point/condensation roundoff, not outward-rounded interval,
continuum or 3D interlaminar certificates. They exclude absent mean/prestrain,
transients and omitted modes. The consumer supplies an authored applicability
envelope; no damage or historical-material threshold is inferred.

Six focused native tests pass using the current original geometry, sections,
directors, condensation, field and fs_math owners with typed error-container
facades. They cover physical norms/plane stress, coherent phase, both layer
faces, tensor order, curved varying directors with nonzero coupled A/B/D,
polynomial reproduction/enclosure, cached-map storage and cold/warm refusal.
This does not establish a complete Cargo workspace or whole-instrument result.

`PreparedElementStrain::real_field` also observes an actual real displacement
snapshot through the same cached map. `RealElementStrain::isotropic_layer_bound`
uses the same plane-stress tensor, signed film faces and spatial convex hull;
one real product and scaled Frobenius norm replace the two complex components
and cycle maximization. Its work ceilings are 4,096 field units and 2,048 layer
units. This bounds space at the supplied instant, including any retained mean
displacement in the input, without asserting a bound between time endpoints or
reconstructing missing prestress/cure strain. The original six focused tests and
two new real-snapshot equivalence/refusal tests pass against the original
numerical owners; a complete instrument render remains a separate check.

## Integrated thermal expansion on the original MITC shell

`shell::mitc::thermal` accepts already integrated engineering thermal strains
in each facet's declared tangent frame. The thickness profile is affine about
the original reference surface; its membrane and gradient coefficients may be
constant or P1-interpolated over a facet. `ShellThermalStrain::from_isotropic`
and `from_orthotropic_in_plane` lower the shared `fs-material` source-integrated
small-strain paths. The orthotropic constructor explicitly aligns principal
axis 3 with the shell normal and rotates axes 1/2 in plane. Callers retain the
original source receipts. `from_faces` requires an explicitly affine profile
and a midsurface reference; it never evaluates a coefficient at the endpoint
and multiplies it by a temperature difference.

`MitcShellModel::prepare_thermal` borrows its immutable original mesh, sections,
directors, support frames and stiffness/mass pencil. Thermal resultants use
the complete existing A/B/D block, including signed membrane/bending coupling.
The original enriched kinematics and quadrature assemble the 20-coordinate
thermal load, followed by the SAME elastic static condensation. Recovery
retains the inhomogeneous bubble displacement `K_bb^-1 f_thermal,b`; a zero
nodal displacement does not silently freeze the internal director response.
The original direct solver resolves compatible displacements with explicit
additional nodal forces and its unchanged true-residual gate.

`ThermalShellResponse` returns physical nodal motion, area-averaged membrane,
bending and transverse-shear resultants, actual thermally shifted strain
energy and material nodal forces. Material forces/storage exclude numerical
drilling contributions and are labeled accordingly, not presented as complete
support reactions. Reference inertia is unchanged. Cold preparation is linear
in facet count and retained fields support repeated recovery without repeating
elastic condensation. Cancellation is polled before work, at quadrature points
and before publishing; the native direct factorization remains uninterrupted.
All results are deterministic ordinary floating point, not error certificates.

The focused `tests/thermal_shell.rs` cases cover heating/cooling free expansion,
biaxial restraint stress and energy, compatible thermal-gradient bending,
signed coupled-section moments, spatially varying thermal bubble forces,
independent energy derivatives, rigid-frame covariance and atomic refusal.
An explicit three-knot source CTE curve drives the actual shell displacement
and oriented orthotropic lowering, discriminating exact path integration from
an endpoint-coefficient shortcut.
This is prescribed-temperature small-strain mechanics, with no heat transport,
thermal feedback, finite-strain geometry update, buckling tangent, interlaminar
recovery or measured-material validation. Nonlinear thickness profiles require
explicit thickness refinement rather than a face-only fit.

`LaminateThermalField` additionally binds the original
`fs-material::laminate::thermal::PreparedLaminateThermal` to each facet through
the same whole-section angle as `ShearShellSection::from_laminate`.
`prepare_laminate_thermal` checks the original A/B/D and shear operators and all
three mass moments before preparing a load. Shared `Arc` states retain the
original layers and source section without copying them per facet. Each layer
may have its own integrated expansion and affine thickness profile; no averaged
CTE, replacement homogeneous material, or inferred temperature profile enters.

Both thermal representations use the same native load assembly, affine bubble
condensation, constrained solve and force recovery. Layered recovery evaluates
physical elastic strain in each original layer, preserving residual stress
energy even when net membrane force and bending moment vanish. The engineering
strain transform and its work-conjugate transpose map section/facet strains and
stress resultants respectively. Transverse shear remains the original shell
law. `LaminateThermalField::layer_stress` exposes the two physical in-plane face
stresses at a supplied compatible strain; it makes no interlaminar claim.

`tests/laminate_thermal_shell.rs` adds independent equal-stiffness bimetal
curvature and irreducible-energy formulas, restraint loads, native energy/force
derivatives, rigid rotation, homogeneous reduction, original mass and atomic
refusal checks. These are authored numerical checks. Thermal histories of Prony
layers, temperature-dependent section updates, finite rotations, and coupled
heat/entropy evolution remain outside this shell adapter.

`PreparedThermalShell::strain_field` exposes the original full-facet strain
polynomial including the thermal director-bubble offset. It shares the existing
Bernstein conversion with cold/harmonic observation. The reduced elastic layer
bound covers all facet coordinates and both layer faces, with engineering-shear
weights gamma^2/2; no unmeasured normal thickness strain is reconstructed.
`LaminateThermalField::elastic_strain_bound` subtracts the original cached layer
free strains using the inverse engineering rotation. Its `section()` getter
lets temperature-refresh consumers reuse the exact already oriented section.
The spatial-gradient regression observes a nonzero interior thermal strain
field even with every physical nodal displacement held at zero.

### Thermal prestress and repeated equilibrium solves

`thermal::prestress` adds the signed P1 transverse membrane initial-stress
form `integral Nxx*w_x^2 + 2*Nxy*w_x*w_y + Nyy*w_y^2 dA` on planar panels with
constant aligned directors. N comes from the actual native thermal recovery,
including its affine bubble. This geometric form has no internal/director
coordinates, so adding it preserves the original elastic Schur complement,
support-frame maps and condensed mass for this stated approximation. Complete
internal director motion, physical rotations, displacement gradients and
support compatibility are admitted before publication. Tension stiffens and
compression may make the returned pencil indefinite; negative eigenvalues are
not repaired or converted to real frequencies. The existing spectral owner
interprets stability. Geometric bending/shear-stress terms, nonlinear base
geometry, follower loads and post-buckling remain outside this extension.

`PreparedMitcStatic` retains an immutable shared native model and its scaled
sparse LDLT factors. Repeated physical loads reuse the factors while preserving
the original ordering, inertia admission, compensated refinement and true
residual gate. The legacy static entry delegates to these same kernels without
copying its model. Shared ownership prevents a different operator from being
paired with the factor. Bounded checkpoints surround the existing sparse kernels.

`RealElementStrain::at_barycentric`, its Bernstein controls and the shared native
integration rule let thermal/work consumers integrate the same complete strain
polynomial without duplicating kinematics or quadrature. This includes the
accepted thermal bubble; controls bound the field and are not physical samples.

`tests/thermal_shell_prestress.rs` checks independent frequency-squared shifts,
compression past a signed stability threshold, free expansion, energy Hessians,
world/support-frame covariance, hidden bubble motion and atomic refusals.
`tests/prepared_mitc_static.rs` checks repeated independent manufactured motions,
legacy equivalence, shared ownership and refusal recovery. Targeted metadata
checks pass; these are not claims of executed numerical or measured validation.

### Retained thermal director inertia

`thermal::transient::EnrichedMitcShell` restores both independent internal
director coordinates per facet in the invertible basis `y=Tq+Br`. Original
nodal stiffness, physical reference inertia, support maps and drilling gauge are retained,
with cross inertia `T^T M B`, bubble inertia and the original bubble stiffness.
Thermal forces enter both momentum equations; internal motion is evolved,
never prescribed by a moving thermal lift or reset to its static equilibrium.

Fixed-step implicit midpoint reuses the existing scaled sparse LDLT and
compensated refinement kernels, with an owned effective factor shared by all
trials and reusable across accepted equal-dt intervals without copying the
sparse operators. The native report separates endpoint total-load work,
polarized quadratic storage change, stiffness-parameter work and complete kinetic-energy change. Actual
strain polynomials, transformed nodal and internal material forces, gauge
storage and whole-facet director bounds are recovered from the accepted full
motion. This transient owner assembles physical I0/I1/I2 translation/director
mass directly from the original element fields, including the nodal and
internal cross blocks. It does not inherit the modal pencil's numerical
relative-spin inertia regularization. The midpoint effective
operator retains the original algebraic drilling stiffness. Fully free
massless director-axis directions require compatible initial displacement and
velocity; unsupported efforts on those numerical directions refuse.
State/operator identity,
stationary support compatibility, finite arithmetic and final cancellation
are admitted before returning immutable proposals.

`EnrichedMitcShell::kinetic_energy_change_j` evaluates actual endpoint velocities
on the original full physical mass. The shared `fs-math::quadratic` reduction
retains endpoint difference/mean low words before the original compensated
sparse action and final contraction. The midpoint owner compares retained
kinetic and work reductions before rounding its report; absolute kinetic
totals remain observations. This keeps small director exchange beside an
unchanged dominant mode and supplies the same increment to thermal consumers.
The operator's symmetry, numerical accuracy and fixed-reference assumptions
remain required; this arithmetic is not a trajectory or interval certificate.

`with_material_model` changes the physical constitutive operator while retaining
the original FULL `T0` basis, physical mass, support map and numerical drilling
stiffness. It assembles `K(T)=Kref+T0^T[Kphysical(T)-Kphysical(ref)]T0` on each
original enriched facet, including the new nodal/internal cross terms. Actual
thermal loads and material forces are conjugate to this same fixed basis;
recovery never applies the new static condensation map to accepted motion.
The current model must retain mesh, directors, constraints and reference mass
moments. Rebinding is immutable and does not move the physical history.

`PreparedMitcMidpoint::advance_from` uses BOTH endpoint operators in the
momentum equation. With `h=dt/2`, its memory-free factor is `M+h C+h²K1`, and its
right-hand side is `M v0+h[f_mean-(K0+K1)q0/2]`. Exact Maxwell average-force
tangents and original-history forces extend these same operators below.
The native first-law account separates
`q0^T(K1-K0)q1/2` from load work and quadratic storage. This diagnostic parameter
work is not itself a material heat source; the coupled thermodynamic owner must
use its original free energy and entropy. The constant-operator entry shares
the same stepping kernel and retains exact zero parameter work without an
extra stiffness application. Numerical gauge forces come from the fixed
reference nodal gauge, with exactly zero independent-director components.

`with_viscosity` attaches one admitted `ShearShellViscosity` per original facet.
The SAME material A/B/D/transverse-shear law, native strain field, positive
quadrature and support-frame scatter assemble sparse C in the complete fixed
`T0q+Br` coordinates. Independent director rates and all nodal/director and
membrane/bending couplings are retained. There is no mass-proportional drag,
modal fit, numerical drilling loss or second constitutive law. The existing
viscosity assembly work bound is checked before allocation; the midpoint factor
is reused over fixed-dt intervals. All-zero viscosity retains the original
lossless path exactly.

Physical heat is independently evaluated from the original constitutive Gram
factors and midpoint strain rates on each facet's positive quadrature. The
native first law is
`delta(kinetic+quadratic+Maxwell storage)+heat-W_parameter-W_load`;
its signed defect is never relabeled as heat. `facet_dissipated_energy_j` names
the original thermal destinations but does not perform a thermal credit.
These per-facet amounts sum to the report's material heat and must not be
debited twice. Source validity and any external thermal deposition remain with
the caller. The viscosity is frozen over one interval; elastic refresh retains
its exact owner, and differing endpoint viscosity owners refuse. An explicit
between-step rebind preserves accepted displacement, velocity and reference
inertia without creating energy or resetting a history. Resolving a physical
decay still requires timestep convergence; midpoint is not L-stable.

`with_relaxation` adds explicit per-facet `ShearShellRelaxationBranch` spectra
above the EQUILIBRIUM sections. The same full A/B/D/shear work and original
positive quadrature retain all nodal/director couplings. Each facet branch
stores twenty local elastic coordinates: eighteen global nodal coordinates
and two residual director coordinates. No global-DOF-by-spectrum history,
modal truncation, drilling loss or static director substitution is introduced.
Source branch count, assembly work and retained branch/ramp/history payload
are bounded before allocation.

For each exact affine numerical displacement increment, the original
`PronyRamp` supplies the branch endpoint, average force and integrated heat.
Writing its average history as `a*r0+c*delta_q`, the existing prepared factor
adds `2*h²*c*K_branch`; its right-hand side adds `-h*a*K_branch*r0`.
The factor and prepared ramps are reused for every equal-step interval.
Material storage changes use the accepted history difference; heat is evaluated
independently in the original passive strain channels. The native report and
facet destinations expose Kelvin and Maxwell heat separately as well as their
sum, with Maxwell storage included in the mechanical first law. Small accepted
coordinate/velocity rounding defects remain numerical residuals.

`with_relaxation_history` requires explicit initial elastic branch coordinates;
zero means relaxed history even at a nonzero mechanical shape. Accepted states
retain the exact source owner. Both interval endpoints and prepared ramps must
share that owner, while elastic material refresh preserves it in the original
fixed basis. Changing a spectrum cannot silently transfer or reset its stored
energy. Failed and cancelled proposals leave the original mechanical and
material histories intact. The branch path is integrated exactly, while the
coupled motion still requires timestep convergence.

The scope is temperature/parameter-dependent LINEAR elasticity in a fixed
reference geometry with constant physical mass, optional frozen linear
viscosity and explicit frozen Maxwell spectra. Nonlinear prestress, remesh,
adjoint and experimental accuracy remain outside this lane.
`tests/thermal_shell_inertia.rs`
checks analytic bubble/cross inertia, independent dynamic internal motion,
static-limit recovery, repeated midpoint conservation and atomic refusals.
`tests/thermal_shell_material_basis.rs` and `tests/thermal_shell_variable_step.rs`
add fixed-basis physical observation and inertia, scalar softened-oscillator
and two-director momentum references, operator-parameter work, compatibility
and cancellation. Compiler checks do not establish numerical execution.
`tests/thermal_shell_viscosity.rs` adds full-operator parity with the existing
modal projection, independent damped-director convergence over three time
steps, nonproportional heat/work closure, rigid-motion and gauge-loss exclusion,
material refresh, exact zero-loss parity and transactional refusals. It also
checks two-branch Maxwell motion, history, storage and independently integrated
Kelvin/Maxwell heat against a continuous ODE reference over three time steps;
spatial A/B/D/shear memory, nonzero prehistory, null rigid motion, source identity,
budget bounds and cancellation exercise the same native sparse owner.

## Equilibrium of tied shell assemblies before modal reduction

`shell::preload::assembly::equilibrate_shell_assembly` solves the full CST/DKT
body in retained physical coordinates before producing its small-vibration
pencil. Components keep their independently specified membrane A/bending D,
reference inertia, essential supports and bonded beams. The existing junction
owner supplies rigid/compliant edge and face ties, including nonmatching and
eccentric interfaces. Its exact sparse map expands component deformation,
transfers dead loads and nonlinear internal forces by `T^T`, and assembles the
consistent equilibrium tangent by `T^T K T`. One shared map per component
avoids a separate dense assembly or per-element transfer matrix allocation.

The original preload membrane energy, analytic force/tangent, load-increment
driver, preconditioner and physical force acceptance are reused. Every accepted
increment must satisfy the original unpreconditioned equivalent-force residual,
moderate-gradient/physical-rotation limits (including eliminated slave nodes),
and positive tangent inertia. Symmetric tangent contributions use identical
mirrored addends. The returned `TiedShellModel` retains its original mass, maps,
glue stiffness/damping and mode expansion, with its stiffness replaced by the
accepted tangent. Glue is already included; it must not be added again. A
single part without ties preserves the original single-shell pencil arithmetic.

All components explicitly declare the existing `ReferenceLoadState::Unloaded`;
already loaded scans refuse instead of receiving a second load from an assumed
stress-free state. Physical supports are required, and disconnected rigid modes
or instability refuse without diagonal regularization. Essential supports on
rigid slave nodes refuse until the junction owner composes those constraints;
supports on masters propagate through the actual tie. Geometry is immutable.
Aggregate limits are 64 parts, 10,000 nodes and 40,000 facets, with the original
tie and nonlinear work budgets. Cancellation/refusal returns no partial body.

Focused regressions compare one-part equilibrium bitwise and a tied nonlinear
cantilever against the same monolithic body, including a slave-seam load, full
tangent and mass. Independent energy/residual differences exercise a
nonmatching eccentric seam. The scope remains moderate-deformation CST/DKT
membrane and bonded-beam centerline prestress with reference-frame linear
bending, torsion and glue. It does not claim nonlinear MITC bending, footprint
contact, torsional buckling, post-buckling or closure of M1.4-W1.

## Complete condensed-film strain coordinates

`MembraneReduction::stretching_strain_gram` and `strain_facets` expose the
original statically equilibrated quadratic strain columns and area-weighted
plane-stress metrics. For the complete symmetric modal-pair vector p, the
same incremental stretching energy is p^T G p/2. Installed tension and linear
bending remain outside this metric. Projection work is bounded before dense
allocation, and no weak direction or mixed term is truncated.

A uniform proportional hereditary spectrum can use an energy factor of this
field without repeating the in-plane solve: every condensed strain column
is equilibrated against the original in-plane trial space, and a common
scalar relaxation law with compatible initial history preserves its range.
This does not authorize heterogeneous spectra, independently evolving tensor
directions, free in-plane initial strains or time-dependent prestress. The
coupled impact owner supplies constitutive states, clock/domain admission and
heat; fs-plate still owns geometry and equilibrium strain energy only.

The inline strain test checks the complete Gram against original positive
facet energy, including mixed pairs and budget refusal. A bounded native run
also exercised the unchanged tangent kernels; it does not constitute a new
full mesh/modal construction qualification.

## Reciprocal tensile connections in assembled shell equilibrium

`shell::preload::assembly::equilibrate_tethered_shell_assembly` extends the
preceding dead-load assembly entry point with massless tensile connections
between physical shell bearings or a bearing and a fixed world anchor. Each
bearing identifies an existing facet, barycentrics and a reference offset;
the supplied reference tension and axial stiffness retain their original N
and N/m meaning. Zero stiffness explicitly prescribes tension. The original
slack/taut radial potential, cancellation-resistant extension, force, material
stiffness and geometric stiffness remain the sole constitutive owner.

Both original point maps compose with the existing component ties. Their
difference is coalesced into deterministic retained columns before force,
Jacobian action or sparse tangent assembly. This retains offset moments,
transferred slave forces, cross-body Hessian blocks and cancellation of shared
coordinates, including endpoints on the same component. Matrix-free actions
use the three-dimensional radial Hessian without forming dense coupling blocks.
The unchanged one-part/world-anchor path keeps its original 18-column arithmetic.

The original continuation scales the complete tether potential together with
dead loads, with the same physical force, deformation and stable-inertia gates.
The returned tangent already contains the connection stiffness and the original
mass is unchanged. Original tether responses report tension, length, relative
potential and equal/opposite endpoint forces, including world-anchor reactions.
Structural stored energy and tether potential changes are explicitly separate.
Final cancellation is polled after both displacements and responses are ready.

The reference body must already be supported and tied sufficiently to have a
stable stiffness; tethers do not silently regularize rigid modes. Work is capped
at 4,096 connections, 262,144 aggregate retained columns and 8,388,608 scalar
tangent contributions. Loaded reference geometries still refuse. This is the
original moderate-deformation shell and linear offset-rotation model with
massless scalar tensile laws, not a distributed vibrating string, finite-rotation
attachment, compressive contact, soundpost friction or complete M1.4-W1 closure.

## Compression-only elastic members in assembled equilibrium

`shell::preload::assembly::equilibrate_connected_shell_assembly` adds physical
compression-only axial members to the same solve as existing tensile tethers,
dead forces and shell/tie internal forces. `AssemblyStrut` supplies an explicit
free length [m] and positive axial rigidity EA [N]. Its stiffness is EA divided
by free length; no contact penalty or missing material property is inferred.
An overlength member models elastic interference preload, while a short member
remains separated until its bearings approach the declared free length.

The original radial connection owner evaluates compression
`C = (EA / L_free) * max(L_free - L, 0)` and the elastic potential
`U = (EA / (2 L_free)) * max(L_free - L, 0)^2`. Its signed radial gradient
and full Hessian retain the negative transverse geometric stiffness `-C/L`.
The same coalesced relative point map carries both endpoint forces, offset
moments, matrix-free actions and symmetric cross-body tangent blocks through
rigid/nonmatching ties. Existing tensile arithmetic and one-part APIs remain
unchanged. The original continuation scales whole connection potentials,
and its physical force and stable-inertia gates accept the final body.

`StrutResponse` reports compression, signed free-length gap, current length,
equal/opposite endpoint forces, absolute nonnegative elastic storage, and a
cancellation-resistant potential change relative to the reference. Connection
energy is separate from structural storage. Open members have exactly zero
force and tangent. At zero gap the open-side tangent is selected explicitly;
a two-sided linearization does not exist there. Small-vibration use is valid
only while the chosen active branch persists. Final cancellation still occurs
after all physical responses and component displacements are ready.

Tethers and struts share the existing aggregate count, mapped-column and
tangent-work caps. The unloaded shell assembly still needs physical supports
and ties; compression does not regularize reference rigid modes. This scope is
a massless elastic axial member with the original moderate-deformation shell
and linear bearing-offset map. It does not model member flexure/buckling,
distributed post inertia, finite contact footprints or friction, and does not
claim complete M1.4-W1 or mixed shell/solid assembly coverage.

## Certified vibrations about an assembled loaded equilibrium

`shell::preload::assembly::vibration::eigen_windows` consumes the accepted
`PreloadedShellAssembly` tangent and reference mass directly through the
existing fs-modal multi-window slicer. It retains each original inertia count,
mode, residual interval and shared-symbolic work receipt. Tethers, elastic
struts and glue stiffness already belong to the tangent and are never added
again. No geometry is rebuilt as stress free. A borrow keeps the source body
unchanged while the returned `AssemblyVibration` is in use.

Original junction maps expand mass-normalized modes and carry incremental
full-component loads by `Phi^T T^T`, including forces and moments at eliminated
nonmatching slave nodes. `perturbation_displacements` returns `T Phi q` and
also serves the conjugate velocity map. `displacements` includes the accepted
mean `T u_equilibrium`, so zero perturbation preserves the loaded body. Both
projections and reconstruction require Cx, poll between modes, and poll around
the original indivisible junction operations before publishing results.

The complete row-major `Phi^T C_glue Phi` projection retains cross-mode damping;
no diagonal approximation or omitted-mode coupling certificate is inferred.
The junction's original directional jump factors evaluate physical glue heat
as a sum of squares and modal damping as their full Gram matrix. Compensated
constraint/jump actions retain small relative motion beside common rigid motion;
the heat query needs no global scratch vector and never clips a negative matrix
quadratic. The preloaded vibration consumer uses that same borrowed-mode owner,
and the exact-response resistance witness uses the same source actions. Native
regressions compare source and assembled sparse power and check common-motion
invariance, weak relative modes, signed cross terms and nonfinite-heat refusal.
Only the existing glue loss is included. Other material, radiation and
frequency-dependent losses retain their own owners. Mode vectors are borrowed
during projection instead of allocating a second full eigenbasis. Explicit
mode-count and dense-projection-entry caps supplement the slicer's own budgets;
these output caps are checked after its bounded solve. Empty windows remain
valid, and nonpositive/zero-crossing mode intervals refuse instead of being
converted into audible frequencies.

This is the original monolithic equilibrium pencil with linear physical port
images, not a new eigensolver, time integrator or nonlinear modal reduction.
Perturbations require the original deformation limits and unchanged active
connection branches. Full-body response convergence and damping certification
outside the supplied windows are separate requirements; this scope does not
claim the complete M1.4-W2 measurement or hundred-thousand-DOF performance gate.

## Exact held-force response about an assembled loaded equilibrium

`shell::preload::assembly::vibration::response::AssemblyResponse` prepares a
fixed-clock forced response from the accepted body and its declared vibration
windows. It forms the full `M_r = Phi^T M Phi` and `K_r = Phi^T K_t Phi` from
the original sparse pencil. Finite modal residuals and orthogonality errors are
retained; identity mass and diagonal eigenvalue stiffness are not substituted.
The original sparse owner checks both projected matrices for positive
inertia, and the actual symmetric inverse of `M_r` must pass an explicit
floating-point residual check. No regularization changes the body.

Canonical state is `[q, p]`, with `p = M_r v`. The existing
`fs-phs::linear_exact::PreparedExactZoh` owns the complete coupled exponential,
held-input evolution, positive glue-heat integral and independent storage/work
accounts. The full `Phi^T C_glue Phi` matrix enters once, including its
cross-mode terms. Preparation is cached; warm steps perform no exponential.
No new integration scheme, diagonal damping approximation, Rayleigh loss or
material/radiation loss is introduced.

The junction retains the original directional viscosity square factors and
jump rows, including rotational width scaling and eccentric constraint terms.
Those rows move from the interface assembly patches; a second global sparse
operator is not stored. A bounded lazy projection through the same tie map and
mode basis supplies `ResistanceWitness::try_from_source_rows`. Its original
rectangular QR and componentwise reconstruction gate certify the unchanged
assembled resistance; no eigenvalue clipping or acceptance relaxation occurs.

Inputs are incremental full-component forces and moments about the already
balanced static preload. Original junction and modal maps supply the force
projection and physical outputs. Zero incremental input and zero canonical
state preserve the loaded equilibrium exactly. Frames include the loaded mean
displacement, endpoint velocity and the conjugate interval flow `Delta u / dt`.
Held-force work uses the actual displacement increment; its agreement with the
original exact owner's independent input integral and the energy balance must
both satisfy the unchanged `1e-12` relative gate. The reported storage is the
quadratic perturbation energy, excluding constant static preload energy.

Explicit projection-work and projection-buffer limits sum mass/tangent work,
retained physical glue-source payload, sparse projection staging, rectangular
witness buffers and QR work. They are checked before corresponding allocation
or numerical loops; numeric factors receive the same separate payload ceiling. The
existing exact owner's scalar-state and augmented-dimension caps remain. With
one independent input per retained mode this path admits at most 42 modes,
refusing larger requests without dropping any modes. The projection ceiling
excludes borrowed eigenvectors and the original bounded preparation cache; it
is not a process-memory or sparse large-body throughput claim. Cx is polled
between projection work units and around existing indivisible sparse/physical
maps and the bounded exact owner. State publication follows the final poll.
Invalid data, failed preparation/energy gates and cancellation leave the
physical state unchanged. Complete canonical snapshots restore bitwise on the
same body, mode order and clock.

G0/G1 regressions cover loaded equilibrium rest and refusals; an independent
forced-oscillator solution at angular steps 2.5, 1.25 and 0.625; full coupled
glue dynamics and heat against a three-rung physical-equation RK4 reference;
the failure of diagonalized glue damping; original FE perturbation storage;
nonmatching slave force/work reciprocity; original directional loss power;
source-projection byte/work refusal and late cancellation; and final-cancel/
checkpoint retry.
This is linear response in the declared windows. Nonlinear deformation,
connection/contact branch changes, frequency-dependent loss, omitted-mode
response and full-body truncation certification remain separate requirements.

## Continuing loads on the original assembled equilibrium

`shell::preload::assembly::continuation::ShellAssemblyLoadPath` prepares the
original unloaded CST/DKT body once and retains its accepted loaded equilibrium.
The original elastic reference pencil, reference mass, sparse mortar maps,
element laws, axial connection maps and reference preconditioner are reused
across subsequent dead-load changes. Accepted tangents never replace the
elastic reference, and displaced geometry is never reinterpreted as stress free.
This is the existing native equilibrium owner, not a modal or parallel solver.

Construction has the original one-shot assembly's arithmetic and load ramp:
dead forces and complete axial connection potentials grow together from zero.
Subsequent `advance` calls start at the last accepted physical coordinates and
interpolate complete old-to-new dead loads through the existing equal-increment
driver. Tensile and compression-only connection potentials remain fully active
throughout that continuation. Their original force, material and geometric
Hessians, slack/separated branches and reciprocal endpoint maps are unchanged.
No second prestress or connection stiffness is added to the accepted pencil.

Each increment must pass the original finite/moderate-deformation admission,
unpreconditioned equivalent-force balance and positive tangent-inertia gate.
The existing residual/Jacobian budget covers the requested load change. Initial
accepted coordinates must satisfy the supplied new deformation and force gates.
Complete full-component loads, including loads on eliminated/support rows,
remain attached to the accepted state; the original work-conjugate tie map owns
their projection. Failed solves, invalid loads, work exhaustion and cancellation
return without replacing any accepted load, displacement, tangent or connection
response. The final cancellation poll occurs after the complete candidate and
owned load arrays are ready. Identical full load bits skip the solve and retain
the entire accepted body; its counts remain those of the last actual solve,
and this no-op does not recertify under newly supplied tolerances.
The accepted equilibrium's stored energy remains the original shell/beam/glue
storage; axial connection potential and endpoint forces remain in their original
tether/strut response records. No load work or connection energy is relabeled
as structural storage or counted a second time.

Five focused G0/G1/G4 regressions exercise an independent isotropic cubic
membrane-energy law through loading, reversal and unloading; nonmatching tied
path partition versus cold equilibrium with original tensile/compressive laws;
preservation of the full accepted connection tangent on the first continued
trial; bitwise first-endpoint and unchanged-load equivalence; failure/domain/
budget isolation; and last-poll cancellation with deterministic retry after
caller source arrays are changed. Reference mass and maps remain unchanged.

These are quasi-static equilibria under the original CST/DKT moderate-rotation
approximation. No elapsed-time dynamics, damping loss, quasi-static path-work
integral, changing material history, changing supports/ties, follower force,
finite-rotation bending, global stability theorem or unstable-branch continuation
is supplied. A subsequent vibration calculation remains local to its accepted
equilibrium and must preserve the chosen axial connection branches. This slice
does not close the broader assembled-prestress bead or claim a calibrated
instrument specimen.

## Continuing axial settings on the original assembled equilibrium

`ShellAssemblyLoadPath::advance_with_axial_settings` continues the complete
component dead loads, each original tether's reference tension [N], and each
original compression-only strut's free length [m] from the last accepted state.
The setting arrays follow the original tether and strut input orders. Accepted
settings are available through `tether_reference_tensions_n` and
`strut_free_lengths_m`; actual endpoint tensions/compressions remain in the
accepted physical response records. Ordinary `advance` retains these accepted
settings while changing dead loads.

The existing equal-increment driver interpolates the scalar settings as well
as dead loads. At each increment the original tension-only and compression-only
potentials supply the direction-updating force and complete material/geometric
Hessian. Strut axial rigidity EA stays fixed and stiffness is recomputed as
`EA / interpolated_free_length`; a blend of endpoint stiffnesses or potentials
would represent a different physical member. Tether axial stiffness stays fixed.
Original endpoints, unloaded geometry, supports, tie maps, element laws, elastic
reference pencil, mass and preconditioner stay fixed. Immutable connection map
columns are shared when preparing candidate laws, and no new geometry or
reference matrix is assembled for a settings change.

The first continued balance check uses the old accepted setting exactly; each
candidate increment must pass the existing deformation, physical force balance
and positive tangent-inertia gates before its endpoint can be accepted. Setting
count mismatch, negative/nonfinite tension, free length below the original
1e-9 m floor, unrepresentable EA/L, solver or work-budget failure, and cancellation
leave the complete accepted loads, settings, equilibrium and response records
unchanged. The final cancellation poll follows all candidate construction.
Exact repeated settings and full load bits take the original no-op path and
retain the last actual solve's reports without recertifying under new options.

Connection potential changes are measured from the original unloaded geometry
using the new setting's law. Their difference between endpoints with different
settings is not an actuator-work audit. Structural storage remains the original
shell/beam/glue energy, separate from the axial connection responses. This
extension supplies no physical actuation time, damping, material history,
prescribed support motion, follower force, unstable-branch continuation or global
stability certificate. The original CST/DKT moderate-rotation and local
small-vibration limits continue to apply; changing a connection branch requires
a new equilibrium and its corresponding tangent.

Four focused G0/G1/G4 regressions cover an independently integrated cubic panel
equilibrium and tangent under simultaneous retensioning, free-length changes
and dead loads, including slack/open transitions; nonmatching tied path
partition/reversal versus cold equilibrium; original intermediate laws and full
sparse/matrix-free connection tangents against force finite differences;
invalid-setting and evaluation-budget atomicity; exact repeated-input no-op;
and last-poll cancellation/retry after caller source mutation. The broader
assembled-prestress bead remains open.

Native validation: `cargo check --offline --locked -p fs-plate --tests` passed,
including all four new regressions. Runtime execution remains unverified: the
shared native test build was killed by the 8 GiB memory limit while compiling
upstream asupersync, before any test ran. No executed-test pass is claimed.

### Graded reference density in tapered round shafts

`shell::stiffened::beam::round_element_with_reference_density` integrates a
positive linear reference-density field with the original linearly tapered
circular geometry. Translational and section-rotary inertia use the same
local density. The existing five-point Gauss rule integrates every kinetic
entry through degree nine; stiffness and its original interpolation stay fixed.
`RoundBeamModes::with_reference_density` admits one endpoint pair per original
geometry interval, including density jumps at material interfaces. It preserves
that field through mesh subdivision and pivot insertion, adds explicitly
attached bodies to the original mass matrix, then computes rigid inertia,
elastic frequencies, mass-normalized modes and reciprocal force rows from
that matrix. The field replaces `RoundBeamSpec::density_kg_m3`; no mean-density
or effective-tip-mass substitution is made. Uniform fields retain the original
scalar-path arithmetic. Density is reference data, not an evolving thermal or
finite-strain density; the existing cold planar Rayleigh-beam and bounded
frequency-slice assumptions remain unchanged.

Focused regressions cover independent polynomial kinetic integrals, constant
profile bit parity, discontinuous-density rigid inertia across three mesh
levels, equal-mass redistribution and density scaling of the actual spectrum,
attached inertia, and invalid-profile/numeric refusals.

### Spatial axial prestress in tapered round shafts

`shell::stiffened::beam::{round_element_with_axial_force,
assemble_prestressed_round_beam}` add the signed conservative initial-stress
potential `1/2 * integral(N(x) * w'(x)^2 dx)` to the original graded Rayleigh
element. Force is in newtons, positive in tension; independent linear endpoint
pairs preserve distributed force variation and jumps at original interfaces.
Three Gauss points integrate the geometric polynomial through degree five.
Original tapered elasticity, translational/rotary mass and attached-body maps
are retained. Zero force preserves the unloaded element bits.

`RoundBeamPencil` exposes the full sparse physical-coordinate tangent and mass
with two coordinates `[w, theta]` per node. Existing support, modal and stability
owners consume this pencil. Compression-induced negative eigenvalues are never
clipped, and a rigid rotation is not falsely projected out under direction-fixed
dead axial loads. This bounded cold owner requires an already equilibrated,
straight, slender reference and caller-supplied supports/balancing loads; it
does not infer thermal stress, solve prestress equilibrium, introduce shear or
follower-load physics, or claim finite-rotation/post-buckling dynamics.

Focused G0/G1/G2/G3 tests cover independent signed polynomial work, interface
refinement and inertia parity, the original unloaded pin-supported spectrum,
three-rung loaded Rayleigh frequencies and Euler buckling, signed stability
on both sides of the critical load, and incomplete/nonfinite/overflow refusals.

### Retained plastic bending and quasistatic springback

`plastic::PlasticBendingPlate<L>` integrates an existing `SmallStrainLaw` through
the thickness of the original DKT elements. Every original area/thickness point
retains its complete native three-dimensional strain and internal variables.
The existing checked plane-stress reduction eliminates transverse tractions and
condenses the same algorithmic tangent; engineering twisting curvature is
converted to tensor shear exactly once. Section moments, constrained membrane
resultants, DKT internal forces and their tangent use those same point returns.
Composite three-point Gauss–Lobatto thickness integration includes both faces
and exactly reproduces the existing elastic DKT stiffness. Yield-front accuracy
depends on the caller's thickness resolution; no plastic-front error certificate
is inferred from elastic quadrature exactness.

`evaluate` returns a candidate for imposed bending. `equilibrate` follows one
caller-selected dead-load/prescribed-coordinate increment using shared sparse
LDLT, a reused symbolic pattern, diagonal equilibration, a recomputed linear
residual and independent force/moment tolerances. Every Newton and backtracking
trial starts from the same accepted history. Unloading retains plastic strains
and can produce residual curvature and springback. State is sealed to its owner;
foreign states, domain violations, failed local/global solves, memory limits and
cancellation cannot mutate accepted history. Cancellation is polled at material
points and around sparse calls; an individual sparse factor/solve is synchronous.

This owner assumes a flat, homogeneous thin plate with zero midsurface extension,
small material strains and a caller-bounded DKT rotation field. A Bernstein bound
checks that quadratic rotation field throughout each triangle. Reported membrane
resultants are constraint reactions, not a solved free-extension condition.
Material axes are mesh axes unless the native law explicitly transforms them.
Finite rotation, membrane stretching, thickness change, forming limits, damage,
fracture, instability continuation and plastic-to-heat conversion are not claimed.
The returned force/tangent is numerical model output, not a material calibration.

G0/G1/G4 regressions cover equality with original elastic DKT stiffness including
twist, an independently derived J2 rectangular-section bending/springback law
under equal-and-opposite principal curvatures, actual incremental nodal-load
equilibrium and release, frozen-history tangent finite differences, and refused
or cancelled attempts with unchanged accepted histories.


## Native shell reference inertia independent of stiffness

`shell::inertia::prepare_reference_mass` exposes the original local P1
midsurface/director kinetic integral independently of elastic coefficients.
Positive P1 thickness and normalized reference density moments retain actual
material inventory and signed first-moment translation/rotation coupling.
`ShellReferenceMass` separates the physical matrix from the existing numerical
relative-drilling matrix. Normal director spin acquires no physical mass.
Existing homogeneous assembly retains its original arithmetic; general moment
assembly and the new preparation share the original point/moment kernels.

`reference_kinetic_bilinear` evaluates both kinetic products on the same
positive reference quadrature without building matrices or subtracting nearly
equal total/gauge energies. Pivoted square completion preserves nonnegative
physical quadratic storage for eccentric inventories. The mesh owner supplies
frame rotation, work-preserving coordinate scaling, aggregate budgets and
cancellation between facets. Retained dynamics must keep the original matter
inventory; no current-density update, finite-rotation kinematics or thermal
transfer is supplied by this kinetic API.

Five G0/G3 regressions in `tests/shell_retained_inertia.rs` define agreement
with native flat/oblique pencils, independent volume integration for tapered
graded layers, signed reference-surface shifts, rigid-spin/pure-drilling energy
separation and domain refusal. Rustfmt and scoped whitespace checks passed.
Native compilation and regression execution remain unverified because the
shared runner exhausted its available disk and memory.


The retained motion API also exposes `reference_kinetic_increment`: all five
before/after kinetic, numerical-loss and endpoint-work products are formed in
one traversal from the same original point amplitudes. This preserves the
relative-drilling work identity near rigid motion without a loose absolute
floor or subtraction from large physical energy. The native implementation is
allocation-free for one facet. A sixth focused regression covers a low-density
nearly rigid transition with gauge energy over 25 orders of magnitude below
physical energy and checks every product against the original mass matrices.

## Source-bound strength on native thermal shells

`thermal::strength` connects the original laminate maximum-stress and coupled
Tsai-Wu criteria to `PreparedThermalShell` at an actual nodal snapshot. The facet
must retain the same shared `PreparedLaminateThermal` as the selected strength
object; equal aggregate stiffness cannot substitute for source, layer,
temperature, reference-surface or free-strain identity. A cloned source snapshot
remains compatible; a refreshed snapshot requires native thermal preparation.

`laminate_strength_bound` recovers the original condensed strain polynomial,
rotates engineering membrane and director-gradient strain to section axes, and
uses the material owner for both faces of every layer in its elastic material
axes. It includes bending and thermal self-stress even when section force and
moment resultants cancel. The six Bernstein controls and two layer faces bound
the complete represented triangle/thickness because stress is affine and both
admitted utilization functions are convex. This is a conservative spatial bound,
not a peak location: an edge control is not an edge midpoint and a bound reaching
one does not prove local failure. `laminate_strength_at` assesses a physical
barycentric point separately.

The observation adds no solve, source query, hot allocation, stiffness change
or damage evolution. Explicit face-evaluation budgets bound 12 evaluations per
layer for a whole facet and two per layer at a point; cancellation is polled
before publication. Inputs and accepted mechanics remain immutable. The scope is
fixed-condition infinitesimal-strain elastic in-plane strength, excluding
interlaminar failure, delamination, time-interval bounds and numerical/model-error
certificates. No missing strength or interaction is inferred.

`laminate_strength_shell` contains analytical rotated-membrane/interaction and
zero-resultant thermal-bilayer oracles, spatial envelope coverage, and original
source/temperature, budget, nonfinite and cancellation refusal regressions.

## Source-bound DKT stress-life fatigue

`fatigue::FatiguePlate` assembles the original native DKT pencil from retained
homogeneous isotropic material states and per-element thicknesses. Its immutable
gauges recover signed normal stress at original midside quadrature points and
signed through-thickness depths. The native convention is epsilon(z) = z B q:
positive transverse-displacement curvature compresses the positive-z face.
Original plane-stress stiffness converts engineering twist exactly once, and
uniform membrane pretension contributes T/h to both normal stresses, including
at zero transverse displacement. Piecewise section properties remain local.

Gauge admission binds chemistry, process and absolute temperature; source-backed
S-N curves also bind the original card, named material state and every fixed
physical query coordinate. Elastic properties must retain fixed external source
conditions. The same section data creates stiffness, inertia and prepared
nine-entry stress rows, so an unrelated assembled pencil cannot substitute for
the observed material. There are no source queries or constitutive solves during
stress recovery.

Complete quadratic rotation-field Bernstein controls enforce the declared
small-slope domain. Linear vertex curvature and both faces bound the engineering
strain throughout each triangle and thickness, including pretension. Explicit
geometry, gauge, sample and aggregate-residue limits bound preparation and
observation. `PlateFatigueObservation` stages all gauges together, retains open
rainflow histories across clones, previews endpoint half cycles without closing
the live history, and finalizes atomically and idempotently.

`tests/fatigue.rs` checks the signed bending/twist relation, both faces and the
midsurface, pretension and rigid motion, original matrix equivalence, heterogeneous
sections, analytic Goodman/Miner cycles, source mismatch and transactional refusal.
This fixed-condition elastic observer excludes laminate-effective calibration,
plastic reconstruction, thermal cycling and unobserved material viscous stress.
Miner damage neither changes the native pencil nor supplies heat, cracks or a
validated service-life certificate.

## Static offset plus harmonic shell strain

`HarmonicElementStrain::isotropic_layer_bound_with_offset` observes a complete
harmonic cycle about an explicitly supplied real strain polynomial from the same
native facet and reference configuration. The complete spatial/thickness
Bernstein hull is retained. A scaled two-dimensional trust-region dual, with at
most 64 secular bisections per control tensor, bounds the affine cycle while
preserving component phase; zero-offset and static-only cases reuse the original
owners. `BIASED_LAYER_BOUND_WORK` admits this bounded query work. G0/G3 regressions
cover scalar/circular/hard-case solutions, phase rotation, representable extreme
scales, combined offset/vibration limit crossing, spatial containment and
cancellation. This is ordinary floating-point linear-shell observation, not an
interval certificate, nonlinear prestress solve or inferred eigenstrain.

## Multiaxial fatigue-limit windows on native DKT plates

`FatiguePlate::assemble_with_findley` prepares normal-Miner gauges, Findley
gauges or both on the same original source sections and native pencil. Findley
gauges recover complete signed plane-stress tensors in [xx, yy, zz, xy, yz, zx]
order, with engineering twist converted once and membrane pretension retained.
Their original calibration binds material/process/absolute temperature and,
when sourced, exact card, named state and every fixed typed query coordinate.

`PlateFindleyBudget` bounds aggregate retained tensor samples, actual plane
counts and declared projection/circle work. The shared material `FindleyCycle`
owns two-dimensional shear-path enclosing circles, complete-tensor closure and
finite-plane selection. Tensor history is retained once per gauge; prepared
projections and source objects are shared by fresh windows. Appends roll back by
saved sample counts, without copying retained history at every timestep.

`PlateFatigueObservation` publishes normal and Findley samples together.
`preview_findley`, `finish_findley` and `findley_results` expose fatigue-limit
usage separately from accumulated Miner damage. An explicit `restart_findley`
seeds the new window at the current accepted tensor and preserves lifetime
normal history. A finished Findley window must be restarted before further
samples. Overall `finish` remains atomic across both models, including a late
closure failure. No repetition count, finite life or stabilized material response
is inferred from a representative stress block.

Four focused tests in `tests/fatigue_findley.rs` cover signed twist/pretension,
native matrix equality, true rotating shear traction on a tilted plane, explicit
window continuation, source context, aggregate budgets and late rollback.

## Coherent harmonic laminate shear and complete-ply strain bounds

`PreparedLaminateShearFacet::evaluate_harmonic` extends the original
equilibrated cylindrical two-probe shear recovery to coherent complex motion.
Physical point responses retain section, material and facet phase and report
cycle-mean recoverable shear work. The original laminate owner caches three
quadratic Bernstein strain maps per ply (twelve additional scalars); the middle
coefficient is a polynomial control, not a physical midpoint sample.

Real and harmonic `isotropic_layer_bound` queries combine the original MITC
in-plane field with those recovered transverse shear maps over the complete
physical ply. Six triangle controls times three thickness controls cover shear
maxima inside the layer even when exterior-face shear vanishes. Bounds retain
common harmonic phase and tensor shear weights. Preparation retains the
original section, reference, inertia and orientation admission; observation
uses fixed local arrays and does not reinterpret the default coating shear or
loss operators.

The supplied Poisson ratio selects an explicit isotropic, zero-eigenstrain
plane-stress normal-strain closure. It is not a transverse-normal law inferred
from an orthotropic source. Known thermal-origin fields refuse that full-tensor
closure because normal free strain is not retained by the shear observer;
their native in-plane and physical ply-shear recovery remains available.
These observations retain the cylindrical approximation and do not claim a
general 3-D/free-edge solution, material failure, discretization-error bounds
or interval-certified floating-point bounds.

Ten focused native tests passed, using the original laminate integration,
MITC assembly/condensation, thermal response and sparse static solve. Coverage
includes coefficient reconstruction, phase/snapshot equivalence, mean stored
work, interface traction, interior maxima, complete-domain containment,
cancellation/reuse and the thermal closure refusal. The existing total-energy
regression now includes the condensed director's native membrane/bending
contribution when comparing recovered shear work with the full shell operator.

## CST/DKT hereditary membrane and bending resultants

`shell::relaxation::ShellRelaxationBranch` admits a passive excess A/B/D
stiffness block and its explicit positive relaxation time. A [N/m], B [N]
and D [N m] are increments above the unchanged equilibrium elastic pencil.
The existing coupled-resultant factor owns passivity, including signed B and
resolved weak directions; no viscosity units or frequency-sampled total
storage tensor are inferred. Material source support remains caller data.

`project_relaxation` uses the original support map, facet transformations,
CST membrane rows and signed three-point DKT curvature rule shared with
instantaneous shell viscosity. Geometry and physical strains are computed
once per facet across its spectrum. Drilling stabilization, prestress,
stiffeners and inertia acquire no material memory. Each complete modal Gram
matrix is accumulated in the existing double-double arithmetic owner.
Exactly equal positive f64 times aggregate; distinct times and every
cross-mode coefficient remain separate. Rank reduction, history and temporal
evolution belong to the existing consuming material-memory owner.

Spatial dimensions, total supplied branches, distinct times and projection
work are bounded before output allocation. Cancellation returns no partial
image. The public shell pencil does not carry mesh identity: matching mesh,
bare CST/DKT support map and basis remain a caller obligation. MITC, tied and
obliquely constrained pencils retain their own projection owners.

`tests/shell_relaxation.rs` checks independent affine A/B/D work with both B
signs and heterogeneous/equal facet times, correlated weak stiffness beside
rigid/drilling nulls, and source, budget and cancellation refusals. The three
new focused native tests and all six existing shell-viscosity tests passed
with the original numerical owners after sharing the strain contraction.
## Source-bound tensor-Prony harmonic CST/DKT shells

`shell::viscosity::harmonic::TensorPronyShell` assembles an original native shell
from retained per-facet `TensorPronyShellSection` laws, thicknesses and reference
densities. Its reference pencil uses the existing CST/DKT assembly, strong
supports, physical axial-rotation convention and selected consistent or lumped
inertia. One explicit reference frequency determines the unchanged numerical
drilling stiffness. The complete source geometry, material and support map are
retained together; matching dimensions of unrelated public models is not
accepted as source identity.

At each positive admitted frequency, the original material owner prepares the
full complex plane-stress constraint, including transverse normal/shear
coupling. The existing physical shell strain owner integrates the complete
anisotropic membrane and bending storage and loss. Current positive physical
storage is assembled directly beside a separately retained original drilling
operator. Subtracting a stiff reference physical matrix to recover soft
low-frequency storage is avoided. Drilling acquires neither material heat nor
frequency dispersion. Original mass is unchanged and never diagonalized.

The shared `fs-modal::harmonic` factor solves
`[K'(omega) - omega^2 M - i K''(omega)] q = f` with peak phasors under
`exp(-i omega t)`. Its optional full physical `B/Z` boundary load enters that
same factorization; all mutual and reactive entries survive. The shared result's
`average_viscous_power_w` is intrinsic material heat from the original Prony and
Kelvin source, counted once. Boundary absorbed work is reported separately.
The equivalent `C_eff = K''/omega` is a frequency-domain identity only, not a
memoryless replacement for hereditary material time integration.

Prepared material reductions recover all six Mandel strain phasors before a
response is published. The original strain ceiling applies to the exact maximum
norm over the complete phase cycle. Both section surfaces and all triangle
vertices are admitted: CST membrane strain is constant, DKT curvature is affine
over the facet, and the homogeneous strain is affine through thickness, so
convexity bounds this finite-element strain field. The original DKT midpoint
ordering supplies vertex extrapolation without a second curvature operator.

`with_frequency_derivative` reuses the same structural factor and the original
complex plane-stress material frequency derivative. Signed storage and
effective-viscosity derivatives retain every spatial cross term; neither is
forced to be a passive operator. Double-double arithmetic evaluates the
loss/omega quotient derivative from the original material values. Authored force
and complete attachment/impedance slopes enter the shared total-parameter solve,
along with explicit inertia and velocity frequency terms. Motion, intrinsic
heat and boundary-work slopes retain independent residual/power checks.

Explicit element/full/free-coordinate and per-assembly work ceilings bound the
native spatial operations. Separate dense numerical payload ceilings include
source matrices and material/resultant arrays before allocation; port and
tangent work retain the shared harmonic owner's own budgets. Cancellation or
source/numerical refusal publishes no partial response. Repeated excitations at
one frequency reuse the material reductions and factor.

The homogeneous constructor supplies a centered CST/DKT thin-shell model. It
does not supply MITC transverse-shear dispersion, nonlinear or prestressed/tied
assembly response, moving supports, DC, startup
transients, spatial convergence or experimental validation. The six focused
regressions in `tests/shell_prony_harmonic.rs` define analytic motion/heat/slopes,
full anisotropic loaded frequency changes, undamped drilling, strongly separated
reference/target moduli, whole-facet strain admission and source/budget/cancel
refusals. Native execution remains pending in an environment that can build the
original dependency graph.

## Independent hereditary MITC resultants

`shell::mitc::relaxation` admits a separate passive excess A/B/D and full
transverse-shear tensor with each positive material relaxation time on each
facet. Signed membrane/bending coupling is admitted through the complete
original-unit passive work matrix. Zero and directional increments are legal;
material tensors do not include numerical drilling. The supplied increments
are above the relaxed equilibrium section, with source frequency/temperature
support remaining explicit caller data.

Projection is bound to the retained `MitcShellModel`: its mesh, directors,
elastic sections, support frames and original equilibrium director-bubble map
are reused. Reduced vectors expand through that model's work-conjugate support
owner. Geometry and equilibrium condensation are prepared once per active
facet; physical strains are evaluated once at each original positive
seven-point quadrature site across the complete spectrum. Instantaneous
viscosity and hereditary projection share the same strain contraction.
Compensated Gram sums aggregate only exactly equal positive f64 times and keep
all cross-mode work. The original mass and elastic pencil are unchanged.

This is hereditary material on the original fixed condensed trial space. It
does not perform frequency-dependent elimination of independent bubbles or
add a second mechanical state owner. Explicit mesh, mode, source-branch,
distinct-time and scalar-work ceilings precede physical expansion and modal
square allocation; arithmetic overflow and cancellation refuse the image.
Focused tests compare original dimensional tensor work, varied directors,
oblique supports, signed A/B/D and independent shear spectra, rigid/drilling
nulls, resolved weak shear, zero spectra and admission bounds.

## Layered tensor-Prony harmonic sections

`TensorPronyShellSection::from_laminate` retains the original immutable
`fs-material::laminate::LaminateSection`, including physical layer faces,
orientations, source laws and declared reference surface. Full elastic layers
use an exact lossless tensor lowering with their original strain limits;
original Prony branches and frequency bands remain intact. Reduced-only
plane-stress layers are refused because they do not declare full transverse
strain constants. `law()` returns the optional homogeneous law; `laminate()`
returns the optional original stack, without an invented equivalent material.

Each frequency integrates storage, loss and signed frequency derivatives through
the existing laminate A/B/D owner and the existing coupled CST/DKT strain owner.
Both reference `model().k` and prepared harmonic stiffness include physical
membrane/bending coupling. Reference drilling remains separate and unchanged.
The existing laminate density-moment adapter supplies exact I0/I1/I2 to the
original consistent mass owner, including signed translation/director coupling.
Lumped mass is explicitly refused for layered sections. Mesh position and normal
must represent the declared reference surface and thickness orientation.

Prepared layer transforms are the same engineering transforms used by the
material owner's `T^T Q T` integration. All triangle vertices at both faces of
every layer recover the original six complex strains and admit the complete
harmonic cycle. The bound covers this affine thin-shell strain field only.
Layer count, source branches and retained recovery arrays enter existing work
and payload ceilings; preparation and admission remain cancellable per layer.
The shared physical-port response and analytic frequency tangent need no
separate laminate solve or finite differences.

This extension does not supply MITC transverse-shear dispersion, transverse
interface traction recovery, delamination, thermal evolution or spatial/model
error certification. The four native regressions in
`tests/shell_prony_laminate_harmonic.rs` cover offset rotated A/B/D, signed physical
inertia, mixed-layer motion/heat/frequency slopes and original layer-domain
refusals. Formatting/parsing checks passed; native execution remains pending
because the original dependency graph exceeds this session's available memory.

## Positive scaling of admitted shell material channels

`ShellViscosity::scaled` applies a finite positive multiplier to the original
A/B/D resultants and its square root to their cached passive strain channels.
The original factor is retained: reconstructing it from rounded scaled matrix
entries can distort a resolved weak correlated direction. Exact zeros, signed
couplings and the factor-one identity survive; nonzero coefficient/channel
underflow or overflow refuses. This fixed-size material operation performs no
new factorization and supplies no source-temperature or frequency authority.
The source-bound shell temperature bridge exercises actual affine CST strain
work across warm/cold factors, including a case where re-factorization lost
4.04% of the original weak-channel scaling.

## Original MITC material projection data

`MitcShellModel::original_geometry` exposes immutable references to the exact
mesh, directors and equilibrium sections retained by assembly.
`ShearShellSection::coupled_matrices` exposes the same facet-frame A/B/D/shear
coefficients used by the native operator, including the original rotation and
shear correction once. These accessors permit source-bound full-tensor material
projection without a duplicate specimen or a changed elastic condensation.

## Original layer harmonic stress, spring storage and material heat

`PreparedTensorPronyShell::prepare_material_observer` binds a complete material
observer to the same prepared CST/DKT frequency, original facets, global support
map, layer orientations and full tensor laws. Geometry and original strain rows
are cached once and shared by repeated harmonic solves and observation. The
observer adds original constitutive spectra; it performs no additional
mechanical factorization or pointwise material inverse.

Each original layer is observed at the existing three triangle midpoint
locations and two positive thickness Gauss points. The signed reference-surface
coordinate and original layer volume are retained. Full six-component Mandel
strain and stress stay in the original layer axes. Equilibrium and Maxwell
spring energy, Kelvin and Maxwell heat, and independent stress power are
integrated with compensated sums. The signed work-minus-heat residual remains
visible. Numerical drilling energy, mechanical kinetic energy and boundary
work are excluded from physical material storage and heat.

Shared whole-cycle admission checks every layer face at every triangle vertex
before integration. The original plane-stress recovery and material strain
ceiling remain binding throughout the affine CST/DKT field. For these kinematics,
the positive quadrature integrates the quadratic material energy exactly in real
arithmetic; it is not a spatial-discretization or continuum-accuracy certificate.
A supplied displacement field is observed without certifying its force equation,
transient startup history, temperature evolution or interface-traction closure.

An explicit additional budget covers original layer points, retained spectra,
one complete returned field and conservative work. Preparation and observation
honor cancellation and return no partial field on refusal. The parent harmonic
payload budget now includes its shared original facet cache.

Three native regression definitions compare rotated, offset mixed-layer fields
with independent original A/B/D quadratic work, assembled harmonic heat, original
volume and phase invariance; exercise homogeneous recovery and factor reuse; and
check payload/work/point caps, source-domain and cancellation refusals, and the
absence of material heat from pure drilling motion. Source review and exact
GitHub blob readback completed. Native compilation, formatting and execution
remain unverified while the workspace is unavailable.

### Graded inventory in loaded reinforced-shell vibration pencils

`shell::preload::{equilibrate_stiffened_shell_with_reference_density,
equilibrate_shell_materials_with_reference_density}` carry the original graded
panel/beam mass through the same nonlinear dead-load and tether equilibrium.
The corresponding shared stiffened assembly owns all density integration and
eccentric inertia. Continuation replaces only stiffness by the accepted native
tangent; the original mass matrix is retained exactly for downstream modes.
Independent membrane/bending sections and both beam inertia choices remain
explicit. Dead forces are caller-owned; density changes do not silently add
weight or alter those supplied forces.

Every entry shares the original supported-geometry/work checks, nonlinear
solver callback, complete physical residual, deformation limits and tangent
stability acceptance. Invalid density coverage/values refuse before the solver
callback. The native moderate-rotation/reference-frame assumptions remain;
no evolving inventory, deformed inertia or additional solver is inferred.
Focused regressions retain graded mass under nonzero preload, distinguish
kinetic changes from fixed dead-load equilibrium, check mass-scaled loaded
Rayleigh quotients with tethers and independent sections, and reject invalid
inventory or an unsupported solver convergence claim.
