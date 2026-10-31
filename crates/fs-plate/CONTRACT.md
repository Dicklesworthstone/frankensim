# CONTRACT: fs-plate

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
Consistent translational and enriched rotary mass is explicit; `ShellMass`
also exposes a separate nodal lumped choice. Drilling inertia is a coordinate
regularization, not a measured physical loss or rotational material law.

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
  `C = mass_damping_per_s*M + stiffness_damping_s*K`. The four parameters
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
  identify physical damping. It retains the supplied DKT/inertia spatial model,
  linear small-deflection kinematics, and the supplied initial-state chain.
  `examples/moving_plate.rs` advances a force across a clamped DKT mesh and
  differentiates sampled deflection through all four mechanical parameters
  and force amplitude, including consistent initial acceleration. The focused
  transient tests compare complete q/v/a and gradients with independent dense
  stepping, check common-scale invariance, undamped energy conservation and
  damped decay, and exercise admission/cancellation/refusal behavior.

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
  consistent DKT rotary-field and torsional beam rotary inertia remain follow-ups. A second
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
evaluates each 18-by-18 element action once per retained mode and then
projects every cross term. Its conservative per-facet work expression is
`1_000_000 + 324*n + 18*n*n`, plus input coverage; this is an operation-count
improvement, not a measured runtime claim. Native affine membrane-energy
checks verify active directions, null modes and complete modal congruence.

## Shared DKT curvature quadrature

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
