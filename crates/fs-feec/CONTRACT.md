# fs-feec CONTRACT

## Signed circumferential Fourier acoustics on the original meridian mesh

`AxisymmetricPressureFe::new_fourier` resolves one complete signed component
`psi_m(r,z) exp(i*m*phi)` under the original peak `exp(+i*omega*t)` convention,
with `p = i*omega*psi`. All original cells, material fields and non-axis pressure
coordinates remain. Native sparse assembly adds the physical angular form
`K_phi,ij = 2*pi*m^2 integral(phi_i*phi_j/(rho*r) dr dz)` to the original
meridian stiffness. Nonzero orders constrain axis pressure exactly to zero
before the existing harmonic factor and expand results to the original full
coordinates. The original constructor and `m=0` operation path remain unchanged.

The angular integral analytically collapses one triangle coordinate and uses
native positive Gauss--Legendre quadrature in the other. Caller-declared
`q`-to-`2q` comparison, `2 <= q <= 32`, `|m| <= 1024`, and original assembly
limits bound admitted work. No artificial axis radius or logarithmic cutoff is
used. The quadrature difference and algebraic residual do not certify spatial,
geometry or physical-model accuracy; physical mesh refinement remains necessary.

Total kinetic energy includes a separately reported angular contribution.
`probe_fourier` reconstructs radial, azimuthal and axial Euler velocity,
Cartesian position and velocity, and active intensity at an authored azimuth.
Every signed component retains the original full-turn `2*pi*r` mass and work
normalization. Real cosine/sine shapes require explicit coherent positive and
negative components; no degeneracy factor is inferred. Geometry, materials and
boundary coefficients must be independent of azimuth, and separate orders do
not couple. P1 fields enforce finite angular energy rather than exact
pointwise `O(r^|m|)` regularity. Nonzero-order vector probes exactly on the axis
refuse because these fields do not establish a unique Cartesian limit there.

## Purpose and layer

Layer: **L3 FLUX** (deps: fs-rep-mesh L2; fs-la, fs-sparse L1;
fs-math, fs-qty, fs-soa L0). The exterior-calculus core (plan §8.1,
Bet 3): cochains on tet complexes, Whitney P₁Λᵏ forms, discrete Hodge
stars, and exact-sequence operator assembly. The exterior derivative
is fs-rep-mesh's INTEGER incidence operator — dd = 0 is exact by
construction (their contract invariant, re-verified here on the zoo);
every metric statement and every approximation lives in the mass
matrices and stars. Exact incidence eliminates algebraic-complex defects;
it does not alone eliminate spurious pressure/EM modes or checkerboarding.
Those formulation-level claims additionally require a conforming subcomplex,
a bounded commuting projection, correct boundary and gauge treatment,
admissible quadrature and coefficient assumptions, and coercivity or inf-sup
evidence. The deferred mixed and curl-curl solve batteries remain the gate.

## Public types and semantics

- `Cochain` — one f64 per k-cell in the complex's canonical
  (BTreeMap-deterministic) cell order, backed by an fs-soa aligned
  buffer, carrying ONE `fs_qty::Dims` tag for the whole field (the
  workspace Qty precedent) and a zero-copy `RawView` descriptor.
  `Cochain::d(&complex)` applies the exact incidence operator (±1
  signed sums) and preserves the dimension tag.
- `fixtures::{single_tet, two_tets, kuhn_cube}` — deterministic tet
  fixtures; `kuhn_cube(n)` is the Freudenthal 6-tet-per-cube unit-cube
  subdivision: conforming, all stored-order volumes positive, pure
  combinatorics (no RNG, no float decisions), the G1 refinement
  ladder.
- `whitney::element_geometry` — per-tet signed volumes, barycentric
  gradients, gradient Gram matrices; Jacobian determinants and
  inverses run through `fs_la::batched` (`batch_det`/`batch_inv` —
  the batched-small-dense consumer this layout was built for).
  Degenerate tets are structured panics (mesh bugs, not runtime
  conditions).
- `whitney::mass_matrix(complex, geo, k)` — P₁Λᵏ Whitney mass in
  deterministic CSR: k = 0 vertex hats (V/20·(1+δ)), k = 1 edge
  Whitney forms (closed form from Gram + scalar integrals), k = 2
  face forms (2(λ_a ∇λ_b×∇λ_c + cyc)), k = 3 per-cell 1/|V|
  (diagonal).
- `whitney::{deram0, deram1, deram2, deram3}` — field → cochain:
  vertex values; Simpson edge integrals (exact through quadratic
  fields); vector-area face fluxes with edge-midpoint quadrature
  (exact through quadratic); centroid cell integrals (exact affine)
  signed by the stored→sorted parity (`sort_parity`) so cell
  orientation matches d2's sorted convention.
- `assembly::incidence_to_csr` — ±1.0 CSR materialization (per-row
  column sort only, no arithmetic); `assembly::stiffness(d, M)` =
  dᵀ·M·d via fs-sparse transpose + Gustavson spgemm (the P₁ Poisson
  stiffness is `stiffness(d0, M1)` — classical FEM, derived from the
  complex).
- `hodge::galerkin_star` (= the mass matrix, SPD, accuracy-first) and
  `hodge::hodge_diagonal_barycentric` (uniform volume shares over
  primal measures: positive on EVERY valid mesh, low-order —
  monotonicity-first). Tradeoffs stated at the definition.
- `weighted::{CellWeight, WeightedMass, WeightedStar, weighted_mass_matrix}`
  adds cellwise material weights on tetrahedra: positive scalar
  weights for degrees 0–3 and symmetric positive-definite Cartesian tensors
  for the vector proxies of degrees 1 and 2. The Galerkin-star entry point
  `hodge::weighted_galerkin_star` uses the same implementation. Coefficient
  jumps are integrated per element before shared degrees of freedom are
  summed. The ordinary and weighted paths share one closed-form Whitney
  element rule; the ordinary path retains its previous arithmetic order.
  `CellWeight::validate` exposes this same scalar positivity, degree and
  interval-SPD coefficient admission without assembling a global Hodge.
  Matrix-free material tangents share the exact gate; geometry admission
  remains separate, and the assembled owner delegates to it unchanged.
  `CellWeight::AffineScalar` retains four positive scalar samples in stored
  tetrahedron order, with independent values on either side of an interface.
  Exact cubic barycentric moments integrate this spatial density through all
  four Whitney degrees; centroid averaging would change kinetic energy and
  the center of mass. Constant affine data uses the original scalar path.
  The same field is integrated by hierarchical H1 mass and diffusion using
  the existing exact-degree Duffy rule through order six. Reference density
  must be paired with reference geometry by the physical owner; no implicit
  current/reference conversion or deformation law is supplied here.
- `weighted2::{CellWeight2, weighted_mass_matrix_2d}` supplies scalar degrees
  0–2 and SPD 2×2 tensor degree-1 pairings on `TriComplex2`. Canonical edge
  orientation follows the complex. Planar thickness and the affine radial
  factor of the axisymmetric measure enter exact barycentric moments;
  centroid-scaled planar mass is not used for an axisymmetric cell. Top forms
  are normalized by the oriented embedded area, not the thickness/radial
  volume. The topology stays two-dimensional.
  `CellWeight2::AffineScalar` supplies three positive cell-local samples in
  stored face order. Planar and axisymmetric mass integrate density with the
  physical measure using exact cubic and quartic barycentric moments,
  respectively. Axisymmetric top-form normalization retains the same joint
  density/radius integral. Uniform samples preserve the original scalar path;
  Cartesian polynomial energy oracles cover all three form degrees and both
  metrics, with face-order and invalid-density controls.
- `WeightedMass` retains coefficient dimensions and matrix-entry dimensions
  `coefficient_dims * metre^(3-2k)` for both tetrahedra and the physically
  weighted `Metric2` pairings. `output_dims` checks propagation to a dual
  result. Rows are test functionals dual to the primal trial basis, not a new
  primal cochain degree. Tensor frames are the stored Cartesian/meridian
  coordinates; transforming material axes into that frame is the caller's
  responsibility.
- `WeightedAssemblyLimits` bounds cells, degrees of freedom, and staged
  triplets before allocation, including orphan-vertex row-storage costs.
  `WeightedError` reports invalid coefficients/geometry, SI exponent or size
  overflow, exhausted bounds, allocation refusal, nonrepresentable entries,
  and cancellation. Coefficient tensors must be exactly symmetric; positive
  definiteness is admitted only when normalized outward-interval Schur
  pivots are strictly positive. Nonpositive and inconclusive cases refuse;
  there is no tolerance symmetrization or negative-eigenvalue clipping.
  Tensor magnitudes are factored outside element integration to avoid an
  overflowing Gram product when the final entry is representable.

### Weighted-operator execution and scope

Each weighted element is a bounded numerical tile. The caller-provided
cancellation callback is also polled during coefficient/geometry validation,
heap-sort sifts, duplicate accumulation, and canonical CSR publication. A
scoped caller supplies `&mut || cx.checkpoint().is_err()`. Allocations are
fallible and all result storage stays private until a final checkpoint
succeeds. There are no worker tasks to drain. The elementary upper triangle
is mirrored so the returned matrix is bitwise symmetric. Contributions are
reduced in deterministic input-cell order; arbitrary cell permutations can
change floating-point rounding and do not carry a bitwise equivalence claim.

The tests in `tests/weighted.rs` and `tests/weighted2.rs` cover independent
constant-field volume/area energies, anisotropy, coefficient jumps, units,
orientation and coordinate scaling, positivity/refusal boundaries, resource
limits, and cancellation with deterministic retry (G0/G3/G4). The 3-D battery
also solves an anisotropic manufactured diffusion problem through the actual
`d0^T M1 d0` composition and checks second-order refinement (G1).
Affine-density controls check first mass moments, nonuniform kinetic energy,
independent volume quadrature for every Whitney degree, interface jumps,
orientation, the exact constant limit and positivity/cancellation refusals.
`tests/material_highorder.rs` compares the H1 restriction with closed-form
Whitney integration and checks higher-mode energy using an overintegrated rule.

`MaterialOperator::prepare_action` prepares the same constant/affine scalar or
tensor H1 pairing without a global sparse matrix. Oriented local descriptors,
one shared quadrature rule and the consistent diagonal are retained under an
explicit storage-byte bound. `MaterialAction::apply_into` contracts and projects
each quadrature field in work linear in local basis size, with caller-owned
scratch and no warm allocation. Output is published only after complete finite
evaluation and a final checkpoint. The diagonal supports caller-owned Jacobi
preconditioning; boundary/gauge restrictions remain the caller's responsibility.
Focused regressions compare assembled action and diagonal, exact graded mass
and centroid without triplet storage, repeated use, and atomic refusal.

`prepare_action_with_scalar_field` freezes a positive scalar multiplier at
caller-selected Duffy quadrature points. It supports general spatial density,
inverse density and density/heat-capacity products without resampling a source
during warm action. Original source errors propagate unchanged, coefficient
dimensions compose explicitly, and sampled storage is admitted before queries.
Independent coefficient scales remain factored through integration so their
separate maxima cannot spuriously overflow a finite combined operator. Tests
use exact quadratic-density moments and the logarithmic inverse-density energy.
The declared quadrature order is a caller-controlled approximation for a general
field; sampled positivity is not an unsampled material-domain certificate.

This includes the constant/affine scalar bulk-operator slice of
`frankensim-ext-feec-weighted-operators-vxth`. It does not close that bead's
high-order, nonlinear/quadrature-point, constitutive-history, or full
formulation proof lanes. Incidence is unchanged, so its exact algebraic
sequence remains unchanged; no new claim about PDE stability, eigenmodes,
boundary/gauge treatment, commuting projection bounds, or continuum error
follows. No cylindrical connection terms or complete axisymmetric vector PDE
are supplied. Geometry arrays on the legacy tetrahedral path remain
caller-supplied and are checked for finite metrics/topology correspondence,
not independently certified against an embedding. Isolated degrees of
freedom can have zero rows; material positivity does not certify positive
definiteness of an arbitrary global constrained/unconstrained problem.

## Additional public types and semantics

- `betti::{integer_rank, betti_numbers}` — exact i128 fraction-free
  (Bareiss) rank of the incidence operators; rank–nullity Betti
  bookkeeping. Fixture-scale certifier: overflow is a checked panic,
  never a wrong answer.
- `highorder::quad1d` — Gauss–Legendre nodes/weights (Newton on the
  Legendre recurrence, FIXED iteration count for bit-stability),
  Lobatto hierarchical shapes (vertex pair + integrated-Legendre
  bubbles), 1D element mass/stiffness by exact-degree quadrature.
- `highorder::hex::TensorSpace` — Q_r tensor-product H¹ on structured
  m³ hex grids (tfz.6 slice 1): 1D dof lattice n₁ = m·r + 1 per axis
  (bubbles have zero endpoint trace — the Dirichlet logic),
  SUM-FACTORIZED matrix-free Poisson apply (per-element axis
  contractions, O(r⁴) vs naive O(r⁶); fixed element order), assembled
  1D reference operators, exact Kronecker Jacobi diagonal,
  tensor-quadrature load and L2 error, and `pcg_matfree` (P6: never
  assemble what we can apply).
- `TensorSpace::prepare_weighted_mass` freezes a positive, potentially
  nonseparable spatial density/material field into `highorder::hex::TensorMass`.
  Original cell coordinates and physical quadrature points reach the source
  exactly once per sample; source errors are retained. The native hierarchical
  Q_r layout and full consistent mass are preserved. Six rectangular axis
  contractions apply the mass in O(r⁴) work when q = O(r), with O(q³) retained
  samples per cell and no element/global matrix. Squared-basis contractions
  prepare the diagonal with the same complexity. Warm `apply_into` has no
  allocation or source queries, uses `ndof + 2*q³` caller workspace, and publishes
  output only on success. SI matrix dimensions are coefficient dimensions plus
  metre³; original reference density belongs on reference geometry for closed
  inventory. Cell-local discontinuities are not averaged across interfaces.
  Positive Gauss sampling uses `r+1 <= q <= 16`; coordinate-wise polynomial
  coefficient degree k is exact when `2q-1 >= 2r+k`. General constitutive fields
  require quadrature-error assessment. Storage is admitted before source work;
  independent per-cell scaling avoids premature overflow of large densities.
- `highorder::simplex::SimplexSpace` — hierarchical high-order H¹ on
  tet complexes (slice 3): Szabó–Babuška entity hierarchy (vertex λ,
  edge λλ·P_k, face λλλ·P_iP_j, interior bubble products) with the
  DETERMINISTIC GLOBAL ORIENTATION convention — every entity's kernel
  arguments use sorted-global-index vertex order, so shared-entity
  traces agree without sign tables; collapsed Duffy tensor quadrature;
  deterministic COO assembly of stiffness/mass/load; entity-based
  boundary masks. Matrix-free/perf paths are slice-4 scope.
- `highorder::derham::TensorDeRham` — the full tensor de Rham complex
  (slice 2): C_r/D_{r−1} 1D factor pair (Lobatto/Legendre; derivative
  operator G in closed form via the integrated-Legendre identity,
  Legendre mass DIAGONAL), Kronecker-assembled grad/curl/div between
  the component spaces E = ((D,C,C),…), F = ((C,D,D),…), W = (D,D,D);
  curl∘grad and div∘curl vanish to machine cancellation (tested);
  canonical commuting projections π_C (endpoint values + derivative
  Legendre moments) and π_D (Legendre coefficients) with
  d∘π_C = π_D∘d by construction.
- `highorder::vecfam` (bead dcng) — the simplicial VECTOR families:
  first-kind Nédélec H(curl) and Raviart–Thomas H(div) on tets at
  r = 1..4 plus discontinuous L² P_{r−1}, completing the simplicial
  P_rΛᵏ tree. NO hand-derived shape functions: per element, Koszul
  bases in centered/scaled Cartesian monomials (RT_r = (P_{r−1})³ ⊕
  ξ·H̃_{r−1}; N_r with the independent cross-product subset α₀ = 0
  when i = 0 — both exact direct sums), square dof-Vandermonde by
  pivoted LU (normal equations were measured at √ε trace error and
  rejected). Dofs are classical moments in sorted-GLOBAL frames:
  edge Legendre moments signed toward the larger global index
  (deram1's direction), face tangential/normal moments in the
  sorted-triple frame and circulation normal (deram2's normal) —
  conformity by construction. `VecSpace` (mass, canonical
  interpolation, L2 error, per-entity global dof layout), `DgSpace`,
  and the chain maps `grad_matrix`/`curl_matrix`/`div_matrix` as
  global interpolation matrices (derivatives exact in monomial
  coefficients).

- `VecSpace::prepare_material_mass` prepares `VectorMaterialMass` on the
  original Nedelec/RT moment coordinates at r=1..4. Positive scalar, positive
  cell-local P1 scalar and Cartesian SPD tensor weights use the existing
  coefficient admission and r+3 Duffy quadrature, exact for their polynomial
  products. Reference density on original reference geometry retains inertia;
  the original unit mass, sorted-global traces and incidence maps are unchanged.
  CSR assembly is symmetric and triplet-bounded. Warm `apply_into` contracts
  trial DOFs into three scalar polynomials, accumulates weighted polynomial
  moments and projects them back: O(n*s + Q*s) for n local DOFs, s monomials
  and Q points, without a dense element matrix, allocation or source query.
  Shared quadrature, flat DOF maps, volumes and the consistent diagonal are
  storage-bounded; caller scratch has ndof entries and output publishes only
  on success. The existing checked FMA capsule supplies supported x86 codegen.
  Units follow the ORIGINAL moments: RT face/cell moments have field units;
  Nedelec edge/face moments have field-times-length units, while its interior
  moments have field units. `coordinate_dims`, `entry_dims(row,col)`,
  `output_dims(row,field_dims)` and `energy_dims(field_dims)` retain these
  differences. Geometry/chart/layout admission is not an independent proof of
  a caller-modified Vandermonde basis. `tests/material_vecfam.rs` covers exact
  density inventory/centroid/energy, actual force-to-acceleration mass solves,
  consistent CSR/action parity, all block scaling under geometry dilation,
  relabeling, anisotropic energy, bounded work and atomic refusals.

- `VecSpace::prepare_material_mass_with_scalar_field` freezes an authored
  positive spatial scalar multiplier into that same vector mass owner. Cell
  weights retain scalar/P1/tensor meaning; field and coefficient dimensions
  multiply. Original cell/barycentric source queries occur once per quadrature
  point after geometry, units and total retained storage admission. Samples
  and their per-cell scales are retained; warm action, diagonal and symmetric
  CSR assembly share them without re-querying evolving material state.
  Independent field/material scales are restored after integration using the
  existing H1 rescaling arithmetic. Caller-selected r+3..=16 Duffy-axis Gauss
  counts retain the existing polynomial contraction cost; storage adds Q+1
  scalars per cell, with no warm allocations. Arbitrary-field quadrature error
  and unsampled source validity remain caller responsibilities. Reference
  density on reference geometry is closed-body inertia, not mass transport.
  `tests/material_vecfam_field.rs` covers exact quadratic density inventory,
  centroid and energy, source freezing, anisotropic CSR/action/diagonal parity,
  inverse-density quadrature convergence, independent extreme scales, and
  source/cancellation refusals without publishing a partial output.

- `cohomology` module (plan §8.1, bead tfz.7): harmonic cochains and
  the discrete Hodge decomposition — the correct treatment of
  multiply-connected domains. `hodge_decompose` splits any k-cochain
  into exact ⊕ coexact ⊕ harmonic in the diagonal-star inner products
  (matrix-free CG on the projection normal equations; components
  re-sum exactly, M-orthogonality residuals reported).
  `harmonic_basis` computes an M-orthonormal kernel basis whose
  dimension equals `b_k` — cross-checked against the integer-rank
  Betti computation (geometry and physics agreeing). `circulation` is
  the cycle pairing extracting Γ from harmonic components
  (Kutta–Joukowski lift); `deflate_harmonics` supplies the
  orthogonality constraints that make saddle systems on handled
  domains well-posed. `fixtures::masked_cube_grid` builds the
  multiply-connected zoo (rings, multi-hole slabs, hollow shells) with
  compacted vertices.

- `differential_characters` (RA.2a) fixes the finite relative-character
  vocabulary consumed by later constructive work: versioned oriented
  primal/dual complexes and mapping-cone relative pairs; integral lattices,
  finite torsion, real cochains, `R/Lambda` characters, and real curvatures with
  `Lambda`-period constraints as distinct nominal sectors. A nonempty relative
  subcomplex `A` requires one degree-`k - 1` trivialization over all of `A`, in
  addition to the separately named terminal-component trivializations.
  Hopkins--Singer character degrees may reach `dim(X) + 1`; raw cellular
  cocycles/cochains stop at `dim(X)`. The vocabulary covers curvature,
  characteristic class, gauge, cup product, holonomy, boundary restriction,
  connecting maps, and both standard short exact sequences.
  Every unary map carries an explicit typed domain, codomain, support,
  coefficient system, and degree shift. Bilinear maps carry both typed inputs,
  their typed output, and an explicit coefficient rule. Cup products require a
  named, versioned coefficient bilinear map plus its immutable artifact
  identity; equality of coefficient type names never invents a multiplication.
  Holonomy pairs with canonical dimensionless integral cycles and returns
  `R/Lambda`; it does not reuse a physical flux or charge lattice as cycle
  coefficients. `delta^2 = 0` is reserved for cellular cochains and `d^2 = 0`
  for de Rham representatives. Exact-sequence schemas expose their image/kernel
  obligations as `RequiresConstructiveWitness`; schema creation never presents
  an unexecuted exactness check as proof. Budgeted canonical schema bytes have a
  deterministic 256-bit replay address, explicitly non-cryptographic until
  wrapped by admitted ledger authority.

- `terminal_relative` (I13.2a, feature `terminal-relative`, [F]) is the first
  constructive physical winding schema. `FiniteCellComplex` admits bounded
  oriented signed-unit chain-complex incidence with coefficients exactly `-1`
  or `+1` and proves `boundary^2 = 0` in `i128`. `TerminalRelativePair` binds a conductor,
  an explicit contained relative subcomplex, insulation, full-dimensional
  conductor components, codimension-one boundary terminal patches, phase and
  orientation roles, complete electrical `PortSchema` values, current/flow
  coordinate selection, and presented Machine-IR graph/port/effort/flow keys.
  Machine references remain presented data: the L6 adapter, not this L3 crate,
  owns graph-membership authority. Pair identity is a strong fs-blake3 semantic
  identity over the complete canonical payload and is declaration-order
  independent. `IntegralWindingRepresentative` is an exact integral relative
  one-cycle representative, deliberately not a homology-class witness.
  `RealCurrentAmplitude`, `DistributedCurrent`, and `GeometricCoil` are
  separate nominal sectors; the only admitted cross-sector relationships are
  named winding-realization and current-realization map declarations carrying
  immutable artifact references. `TerminalRelativeComplexReceipt` retains the
  typed identity receipt, coefficient and incidence domains, current units,
  terminal/PortSchema and presented Machine-IR bindings, admitted conversion
  families, schema version, authority status, and explicit no-claim states.
  Integral chains and `IntegralRelativeCochain` values use the phase-owned
  component quotient basis; exact boundary, integral coboundary, and evaluation
  pairing expose the integer Stokes identity without a metric/Hodge coercion.
  `TerminalRelativeSignedRelabel` admits a representation witness only when a
  complete signed cell bijection commutes with every exact incidence and maps
  conductor, relative, insulation, component, and terminal supports exactly.
  `TerminalRelativePhysicalRelabel` additionally carries complete explicit
  component, phase, and terminal bijections plus a per-phase current-coordinate
  sign; no semantic permutation or orientation compensation is inferred.

- `integral_topology` (I13.2b, feature `moonshot-integral-topology`, [M]) begins the
  exact integral topology checker behind a second default-off gate. Its first
  tranche admits bounded dense row-major `i128` matrices and verifies a complete
  caller-supplied Smith witness only when explicit integer `U`, `U^-1`, `V`,
  and `V^-1` matrices prove both inverse orders and `U A V = D` exactly.
  Canonical `D` is diagonal, nonnegative, has a zero suffix, and each positive
  invariant factor divides its successor. The opaque verified value retains
  the exact source and all witness matrices and is unconditionally tagged
  `AbstractAlgebraOnly`; it is not yet an `IntegralRelativeTopologyReceipt`, a
  homology computation, or physical R3 winding authority.
  Its second tranche extracts an opaque exact boundary matrix bound to the
  admitted pair identity, phase, component, degree, and canonical source/target
  quotient bases. Rows are degree `k-1`, columns are degree `k`, and the
  unaugmented edge maps retain distinct `0 x dim(C_0)` and `dim(C_top) x 0`
  shapes. Terminal orientation, port trivialization, and phase-current signs
  do not alter intrinsic cellular incidence.
  Its third tranche consumes adjacent pair-bound maps `A_k`, `A_(k+1)` plus a
  complete verified `U A_k V = D` witness, computes `V^-1 A_(k+1)`, requires
  the first `rank(A_k)` rows to vanish exactly, and retains only the lower
  kernel-coordinate image. Those rows are coordinates in the retained Smith
  witness's kernel basis, not canonical cells or canonical homology generators.
  Its fourth tranche constructs `D`, `U`, `U^-1`, `V`, and `V^-1` by a
  deterministic minimum-magnitude, row-major Euclidean reduction. It repairs
  divisibility against the entire trailing block rather than accepting a merely
  diagonal matrix, caps live storage, elementary operations, entry work, and
  every source/generated coefficient, and publishes only after the independent
  first-tranche verifier rechecks the complete witness. Construction work stays
  attached to the returned construction receipt alongside the opaque verified
  authority.
  Its fifth tranche binds a complete Smith authority byte-for-byte to the
  retained lower kernel image and owns both complete inputs. For lower image
  rank `s` in a verified kernel of rank `q`, it publishes free rank `q-s`, all
  positive presentation factors, and only factors greater than one as
  nontrivial torsion summands. This is exact phase-local cellular chain
  homology, not generator, period, naturality, embedding, or physical-R3
  authority.
  Its sixth tranche lifts the nontrivial torsion and free presentation columns
  into the original pair-bound `C_k` basis as
  `V_A[:, rank(A_k)..] U_L^-1[:, unit_factors..]`. It retains the corresponding
  torsion filling columns from `V_L`, then checks `A_k G = 0` and
  `A_(k+1) H_j = d_j G_j` exactly for every torsion column before publication.
  Unit-factor columns are omitted. The complete prior Smith authorities remain
  attached, making finite orders exact while keeping representative bytes
  explicitly witness-relative.
  Its seventh tranche constructs free integer cocycle columns as
  `V_A^-T[:, rank(A_k)..] U_L^T[:, rank(L)..]` in the induced dual of the
  retained pair-bound `C_k` basis. It checks `A_(k+1)^T Omega = 0` and retains
  the exact generator-by-cocycle period matrix `G^T Omega = [0 | I]^T`. This is
  only the witness-relative free `Hom(H_k,Z)` dual; integer periods kill the
  torsion generators and do not stand in for torsion cohomology or linking.
  Its eighth tranche accepts one strongly pair/phase/degree-bound integral
  `C_k` chain, proves `A_k x = 0`, and classifies the cycle in the retained
  generator basis. Nontrivial torsion coordinates are canonical nonnegative
  residues modulo their positive invariant factors; free coordinates are
  independently replayed as `Omega^T x`. It retains an exact original-basis
  `C_(k+1)` adjustment and verifies the complete decomposition
  `x = G gamma + A_(k+1) h`. Consequently a zero coordinate column carries an
  explicit filling witness, while a nonzero coordinate column is a successful
  homology classification rather than a refusal.

## Invariants

- dd = 0 EXACTLY: integer path (i64 `apply`) and f64 CSR path (sums
  of ±1 are exact) — both tested per fixture.
- Orientation coherence: edges u→v (u<v), faces by sorted a→b→c
  circulation, cells by sorted-order volume sign; the de-Rham maps
  COMMUTE with d (R∘grad = d₀∘R₀, R∘curl = d₁∘R₁, R∘div = d₂∘R₂),
  which pins every sign convention against the analytic operators.
- Whitney spaces contain constants: interpolating a constant field
  reproduces its L2 energy exactly (⟨Rc, M_k Rc⟩ = |c|²·V) for all k.
- Mass matrices are symmetric positive definite; the stiffness
  composition is symmetric and annihilates affine fields on interior
  rows (patch test).
- Assembly is deterministic: canonical cell order + fs-sparse's
  push-order-independent COO accumulation.
- Relative-character construction rejects degree overflow, coefficient-sector
  coercion, primal/dual mixing, missing or mistyped whole-`A` and terminal
  trivializations, non-subcomplex cell counts, and malformed exact-sequence
  composition. An ambient dimension of 255 is inadmissible because the `u8`
  degree representation could not encode the required `dim(X) + 1` character
  degree. Curvature codomains retain the character's lattice-period constraint.
  Boundary and terminal declarations are canonicalized by identity before
  hashing, so insertion order cannot change semantics.
- Terminal-relative admission requires the relative subcomplex to be contained
  in the conductor and insulation to be explicitly contained in that relative
  subcomplex. Components cover the conductor while partitioning its top cells
  (lower-dimensional closures may be shared), and every terminal patch lies in
  its component and the relative subcomplex, avoids insulation, selects
  electrical flow/current, and meets exactly one component top cell. Every
  phase carries exactly one driven plus one return/reference terminal, both
  current directions, and a common version/shape/basis/frame/clock/power/
  conservation/voltage-reference convention. Its two terminals name exactly one
  conductor component, every component has exactly one such phase binding, and
  two phases cannot silently claim the same component. Geometric-coil rows must
  name the component owned by their phase.
  Relative boundary and real/integral coboundary use the canonical phase-local
  quotient basis; integral accumulation and chain/cochain evaluation are exact
  and refuse publication overflow. If the admitted relative subcomplex equals
  the conductor, every corresponding quotient basis is empty: zero chains and
  cochains remain valid typed values, their maps remain empty, and their exact
  pairing is zero. This is a zero quotient-chain-complex statement, not a
  computed claim about homology.
- A signed relabel is bound to exact source and target pair identities. Its
  canonical identity is independent of input row order; inverse/composition and
  integral chain/cochain/representative transport retain exact cell signs,
  basis reindexing, boundary/coboundary commutation, and evaluation pairing.
  This first relabel lane preserves component, phase, terminal, port, and
  Machine-reference semantic identities rather than silently permuting them.
- A physical relabel validates the component-support, phase/component, and
  terminal/phase/component squares against its signed cell action. Preserve
  keeps terminal role, physical orientation, and trivialization sign; Reverse
  flips all three. Port physics, presented Machine authority, and the explicit
  voltage/current reference artifacts remain exact; only mapped terminal-local
  nominal bundle identities may differ. Generic integral chains/cochains
  receive only the cellular sign, while a winding representative receives the
  combined cellular/current sign and a real current amplitude receives only
  the declared current sign. An already-physical `DistributedCurrent` receives
  only the cellular sign and must be re-admitted with two target constraint
  receipt IDs that reuse neither source receipt; those presented IDs are not
  themselves verified authority, and a target current-realization map remains
  a separate explicit declaration. A `GeometricCoil` is only re-declared on
  the explicitly mapped phase/component with caller-supplied connectivity and
  manufacturing artifact IDs that reuse neither source artifact; no geometry
  is transported, and a target winding-realization map remains separate. A
  physical-map redeclaration accepts only caller-supplied target endpoints on
  the exact mapped pair and phases, preserves each endpoint's nominal sector
  and the source map family, and requires a fresh map ID and artifact. It does
  not synthesize target objects; those come from the sector-specific transport
  and redeclaration APIs above.
  Checked combination occurs before coefficient application so two reversals
  cancel exactly without a spurious intermediate `i64` overflow.
- Vector families (vecfam battery): dimension counts r = 1..4 match
  the closed forms (N: 6/20/45/84, RT: 4/15/36/70 per tet) and the
  exact-sequence Euler identity Σ(−1)ᵏ dim Vᵏ = 1 holds on
  contractible fixtures; dof-Kronecker unisolvence to 5.8e-13 and
  mass SPD at r = 2, 3; tangential (2.3e-12) and normal (1.4e-11)
  trace continuity across shared entities at r = 2..4; curl∘grad ≤
  1.7e-14 and div∘curl ≤ 1.4e-13 through the discrete chain at
  r = 2, 3; canonical-interpolation ladders at order r (slopes
  1.90/1.88 at r = 2, 2.92/2.89 at r = 3); the two-tier G3 battery —
  edge dofs transform with definite parity (−1)^{k+1} under
  relabeling (1.1e-16), interpolated fields are label-invariant
  (6.6e-13); r = 1 members reproduce the Whitney forms to 4.4e-16.

## Error model

Structured panics on shape mismatches, out-of-range degrees, and
degenerate tets (programmer/mesh errors). `TensorSpace` seals its extent-bearing
fields; tensor-product lattice, local scratch cube, element-matrix, and global
DOF extents are checked by its constructor before allocation and panic rather
than wrapping. Public lattice-index helpers reject indices outside that sealed
space.
`betti` panics on i128 overflow rather than degrade exactness. The legacy mesh
and operator paths remain panic-on-invalid-programmer-input because upstream
certificates own their validity. The relative-character schema is an
agent-facing admission boundary and therefore returns structured
`CharacterError` refusals instead of panicking.
The feature-gated terminal-relative boundary returns structured
`TerminalRelativeError` refusals for malformed incidence, support, topology,
port, phase, coefficient, identity, budget, and declared-map semantics. It
never rounds an integral cycle test or silently coerces nominal sectors.
The feature-gated exact topology checker returns `IntegralTopologyError` for
matrix shape, retained/workspace entry, scalar-work, allocation, cancellation,
checked-arithmetic, inverse, transform, and canonical-diagonal refusals. No
overflow or incomplete witness is reinterpreted as rank, torsion, or success.
Pair-boundary extraction additionally distinguishes unknown phase/excess degree
counterexamples from typed component/incidence budget, invariant-loss,
allocation, and cancellation Unknowns.
Kernel-coordinate transport refuses pair/phase/component/degree/basis or Smith-
source mismatches and any exact nonzero nonkernel coordinate. Extent, binding,
retained-storage, scalar-work, arithmetic, allocation, invariant-loss, and
cancellation failures remain typed Unknowns.
Constructive Smith reduction treats coefficient growth, `i128` arithmetic,
live-entry, operation, entry-step, allocation, cancellation, and runtime-
invariant failure as typed Unknown. In particular, `i128::MIN` is refused
because its positive invariant magnitude is not representable. A constructor-
generated witness that fails independent structural verification is an internal
Unknown, not a Refuted judgment against the caller's valid source matrix.
Quotient-homology binding refuses lower-image Smith shape or source-byte drift
as Refuted. Binding/retention limits, cancellation, and opaque-input invariant
loss remain Unknown; factors equal to one reduce boundary rank but never appear
as nontrivial torsion.
Generator lifting consumes only opaque verified homology. Aggregate output,
retention, scalar-work, allocation, arithmetic, cancellation, and internal
shape/equation failures are therefore typed Unknowns, not counterexamples to
caller-supplied topology. No partially checked representative is published.
Free-cocycle construction likewise consumes only opaque verified generators;
output/retention/work, allocation, arithmetic, cancellation, invariant, and
coboundary/pairing failures are typed Unknowns. It never narrows `i128`
cocycles into the existing `i64` relative-cochain type implicitly.
Integral-cycle classification treats a strong pair, phase, or degree mismatch
and an exact nonzero `A_k x` row as Refuted. Opaque-authority shape drift,
coefficient-arity loss after successful binding, output/workspace/retention or
multiply-accumulate budget refusal, allocation, checked arithmetic,
cancellation, period disagreement, and decomposition disagreement remain typed
Unknowns. No partial coordinates or adjustment are published.

## Determinism class

Bit-deterministic cross-ISA by construction: fixture generation is
pure combinatorics; element geometry is batched fixed-order mul_add
(inheriting fs-la::batched's batch-membership invariance); assembly
uses canonical orders; the only transcendental in the crate is
`fs_math::det::sqrt` in the diagonal star's primal measures. Golden
FNV-64 over all four mass matrices, the Poisson stiffness, Betti
numbers, and element volumes on kuhn(2): `0xa973_ca6b_07c3_9639`,
recorded on Apple M4 Pro (aarch64), verified identical on
Threadripper (x86_64). No libm anywhere in the hashed pipeline.
Const-P and runtime-P contractions execute the same ascending-order fused
operations, including zero coefficients. Exceptional field values therefore
propagate consistently instead of being hidden only by one degree path.
Relative-character canonicalization is integer/byte-only: version tags,
length-framed identifiers, exact dimensions/counts, normalized coefficient
bits, sorted boundary identities, and explicit budgets feed a stable four-lane
digest. That digest is deterministic replay evidence, not collision-resistant
or authenticated authority. The encoder accounts for the logical byte length
while refusing to extend its materialized buffer past the configured canonical
byte limit.
Terminal-relative pair and representative identities likewise use sorted
integer/byte canonical frames and strong fs-blake3 schema domains. Floating
current/cochain values are intentionally outside those identities in this
slice; finite-value checks are admission evidence, not cross-ISA numerical
equivalence evidence.
Integral-topology witness verification is integer-only, row-major, and
fixed-order. Exact product traversal, error precedence, invariant-factor order,
and scalar-work counts are deterministic across ISA and thread counts. The
first tranche deliberately emits no persistent identity or promotion receipt.
Pair-boundary extraction traverses canonical `CellRef` and admitted-incidence
order and binds the complete pair identity, so declaration order cannot change
its bases, signs, matrix bytes, or work count.
Kernel transport compares bindings in canonical order and evaluates every
`V^-1 A_(k+1)` coordinate row-major with checked `i128` dot products. Different
valid Smith witnesses may intentionally produce different lower matrices.
Constructive Smith reduction selects the least unsigned-magnitude active entry
with row-major tie-breaking, uses checked Euclidean quotients, and records exact
operation/entry counts. Identical source bytes and budgets replay the complete
constructed witness and counters byte-for-byte; relabeled inputs are promised
only the same canonical invariant factors, not the same transforms.
Homology binding compares lower-image source bytes and verified invariant
factors in fixed row-major/factor order. It retains no separately normalized
generator basis, so witness-dependent lower coordinates remain replay authority
rather than being mislabeled canonical representatives.
Generator lifting visits torsion then free presentation columns in their
verified Smith order and evaluates every checked `i128` product row-major.
Identical retained authorities and budgets replay generator, filling-chain,
work-count, and scalar-count bytes exactly. Alternate valid Smith witnesses may
and generally do select different representatives of the same quotient.
Free-cocycle lifting visits the original dual basis and free presentation rows
in fixed order, then checks incoming columns and generator columns row-major.
The cocycle and generator representatives transform together under alternate
valid Smith witnesses, while the retained zero-then-identity period matrix and
all counters replay exactly for identical authorities.
Integral-cycle classification visits outgoing rows, outgoing-Smith
coordinates, lower-presentation coordinates, free periods, incoming adjustment
coefficients, and decomposition rows in fixed order. Torsion normalization uses
`div_euclid` and `rem_euclid` with positive retained factors. Identical
authorities, typed chain bytes, and budgets therefore replay class coordinates,
adjustment, work counts, and checked multiply-accumulate counts byte-for-byte.

## Cancellation behavior

Bounded synchronous assembly loops; chunking a large mesh to tile
quanta with Cx poll points between chunks is the fs-exec driver's
job (the fs-la/fs-simd discipline). No internal threading.
Relative-character schema construction is bounded by explicit cell,
boundary-component, canonical-byte, and coefficient-product budgets and
performs no unbounded search or solver work; the downstream RA.2b/RA.2c
constructive checkers own Cx polling, drain, checkpoint, and witness budgets.
Terminal-relative construction and exact chain maps are synchronous and bounded
by explicit cell, incidence, component, terminal, and canonical-byte limits.
This schema slice performs no solver, search, or unbounded iterative work.
Integral-topology verification preflights all six retained matrices, scratch
entries, and exact scalar work before allocation. It polls before every output
scalar and again before final publication, so at most `max(rows, cols)` checked
inner-product terms occur between polls. Cancellation publishes only a typed
refusal; no partially verified value is constructible.
Pair-boundary extraction separately caps component visits, incidence visits,
matrix extents/entries, and retained basis-plus-matrix entries. It polls each
deterministic visit and once after the final incidence before publication.
Kernel transport caps extents, output and total retained entries, binding
comparisons, and exact dot terms. It polls every comparison/output scalar,
around allocation, and immediately before allocation-free publication.
Constructive Smith reduction polls every scan item and before each elementary
operation; a swap updates at most `2 * max(n + 2m, m + 2n)` destination
entries, while shears and sign normalization update at most half that bound.
It polls at the verification handoff and then reuses every poll, including final
publication, from the exact witness verifier. Cancellation discards the private
partial reduction and publishes only a typed Unknown.
Homology decomposition preflights source/factor comparisons and retained
authority entries, polls every comparison, and polls again immediately before
allocation-free publication.
Generator lifting preflights both output matrices, total retained entries, and
all exact scalar work before allocation. It polls before allocation, before
each output coefficient or chain-equation row, and immediately before final
publication. Between polls it performs at most
`max(cycle_rank, dim(C_k), dim(C_(k+1)) + 1)` checked scalar terms; cancellation
discards both private matrices and returns exact completed/planned counters.
Free-cocycle construction preflights its cocycle and period matrices, retained
authority, and scalar work before allocation. It polls each cocycle coefficient,
incoming-boundary equation, and generator pairing plus final publication, with
at most `max(cycle_rank, dim(C_k))` checked terms between polls.
Integral-cycle classification preflights output, retained authority, scratch
coordinates, and every checked multiply-accumulate term before allocation. It
polls each cycle equation, coordinate, period, adjustment, and decomposition
row plus preflight, allocation, and final publication. At most
`max(dim(C_k), cycle_rank, boundary_rank, generator_count + dim(C_(k+1)))`
checked terms occur between polls; cancellation discards all private
coordinates and returns exact completed/planned counters.

## Unsafe boundary

One registered leaf capsule:
`src/highorder/fma/mod.rs` uses an x86-64
`#[target_feature(enable = "avx2,fma")]` function so explicit
`f64::mul_add` operations in the const-P apply compile to native fused
instructions instead of baseline-x86 libm calls. Its safe dispatcher checks
both CPU features immediately before the only unsafe call; the numerical body
is otherwise safe slice code. The boundary and precondition are documented in
`src/highorder/fma/SAFETY.md` and registered in `unsafe-capsules.json`. All
other crate code inherits the workspace `unsafe_code = "deny"` policy.

## Feature flags

- `terminal-relative` — default-off [F] I13.2a schema and exact-incidence lane;
  enables the optional `fs-blake3` and domain-neutral `fs-ports` dependencies.
- `moonshot-integral-topology` — default-off [M] I13.2b exact-integer checker; implies
  `terminal-relative` and cannot promote abstract algebra to physical evidence.

## Conformance tests

`tests/cohomology.rs` emits schema-validated `fs-obs` JSONL under suite
`fs-feec/cohomology`. The four fixed-input aggregate identities are `ch-001`
through `ch-004`, each with input seed zero. `ch-003/measurement` is a distinct
Custom companion retaining the coarse/refined circulation and lift values; it
is emitted at the original diagnostic boundary before the final lift gate, and
the `ch-003` aggregate is emitted only after that gate succeeds. The canonical
measurement is replay evidence for the existing G1 tolerances, not an upgrade
of their authority or discretization claim.

`tests/feec_battery.rs` (12 tests, 28 canonical JSONL rows): dd = 0 on the
fixture zoo (single/two-tet, kuhn 1–3) over random integer cochains
AND via CSR spgemm; Kuhn fixture positively oriented, conforming,
unit total volume; Betti (1, 0, 0, 0) on all ball fixtures via exact
integer rank; de-Rham commutation for quadratic-gradient, affine-curl
and affine-divergence fields to 1e−13 (two fields for div, one
divergence-free); constant-field energy exactness for k = 0..3 +
SPD spot checks; stiffness affine patch test (interior residual
< 1e−12) + symmetry < 1e−14; G1 MMS primal Poisson (sin πx·sin πy·
sin πz oracle, PCG solve at 1e−12, M₀-weighted L2) with measured
orders ≈ 2.0 on the n = 4/8/16 Kuhn ladder; Hodge star positivity +
dual-volume partition (Σ = |Ω|); Cochain container semantics (dims
tag through d, container dd = 0, 128-byte alignment); cross-ISA
golden hash.

The primary battery emits schema-validated `fs-obs` JSONL under suite
`fs-feec/battery`. Its 28 retained rows use collision-free identities:
`dd-zero/{fixture}` for the five fixture-zoo rows; `kuhn-fixture/n-{1,2,3}`;
`betti/{fixture}`; `hodge/k-{0,1,2,3}`; and the singleton cases
`dd-zero-propcheck`, `commutation`, `whitney-patch`, `stiffness-patch`,
`mms-order`, `cochain`, and `unit-cube-boundary`.
The DD rows record input root 5, kernel `0xFEEC`, tile 1, fixture ordinal, and
three trials while preserving the single continuous stream across the zoo.
The Whitney patch records input root 9, kernel `0xFEEC`, and tile 2. The
generated DD case records property seed `0xFEEC_4A48_0001`; every other fixed
input aggregate uses seed zero.

The three original MMS diagnostics remain info-only `Custom` events named
`mms-poisson` at `mms-order/measurement/n-{4,8,16}`, after each solve and before
the later ladder and slope gates. `feec-golden/measurement` remains an info-only
`Custom` event named `feec-golden`, recording the formatted legacy value plus
actual and expected hashes before the unchanged equality assertion; no passing
golden aggregate is synthesized. Measurement payloads are object-shaped,
encode non-finite raw floating values as `null`, and explicitly record the lack
of an execution seed. All rows retain their original emission boundaries and
legacy detail prefixes. In particular, each `hodge/k-{k}` row covers only that
degree's positivity and shape gates and remains before the terminal dual-volume
partition assertion. These events add replay evidence without upgrading the
existing G0/G1, cross-ISA, golden, topology, convergence, or performance claims.

Bead 4nh8 adds a separate 600-case G0 property suite rooted at seed
`0xFEEC_4A48_0001`: it selects the same deterministic fixture zoo and cycles
shrinkable integer seeds into exact vertex/edge cochain extents (empty seeds
mean the zero cochain), then requires integer `d₁d₀ = 0` and `d₂d₁ = 0`
exactly. The historical three Philox trials per fixture and golden
`0xa973_ca6b_07c3_9639` remain unchanged; this generated suite is scoped to
the fixture zoo rather than claiming exhaustive topology coverage.
`tests/highorder_battery.rs` (slice 1) emits schema-validated `fs-obs` JSONL
under suite `fs-feec/highorder`. Its aggregate identities are
`gauss-legendre`, `lobatto`, `sumfact-vs-assembled/m-2/r-1`,
`sumfact-vs-assembled/m-2/r-3`, `sumfact-vs-assembled/m-1/r-5`,
`sumfact-vs-assembled/m-3/r-2`, `jacobi-diag`, `mms-ho-order/r-1` through
`mms-ho-order/r-6`. Gauss–Legendre exactness to degree 2n−1 (n = 1..10) plus
node symmetry, Lobatto endpoint structure, the Jacobi diagonal check, and the
G1 MMS Poisson ladders are fixed-input cases with seed zero. The sum-factorized
reference fixtures and golden-hash measurement record their actual `StreamKey`
base seed 6; their replay metadata retains kernel `0x460` and tiles 51, 53, 45,
62, and 90 respectively.

Each sum-factorized aggregate has a same-prefix `/measurement` `Custom`
companion retaining `(m, r)`, deviation, scale, and stream provenance. The
twelve MMS diagnostics use distinct
`mms-ho-order/r-{r}/measurement/m-{m}` identities: `m = 4, 8` for `r = 1, 2`
and `m = 2, 4` for `r = 3..6`. They remain at the original per-ladder point
after PCG convergence and before the order gate; the corresponding aggregate
is emitted only after that gate passes. `ho-golden/measurement` retains the
observed and expected hashes at the original info-only boundary before the
unchanged equality assertion; no post-assert aggregate is synthesized. Every
companion and aggregate is failure-linted, serialized, wire-validated, then
printed; floating measurements serialize as finite JSON numbers or `null`.
These events preserve the existing 1e−12 sum-factorized roundoff gate, slope
gates ≥ r + 0.6 for r = 1..6 (measured ≈ r + 1), and frozen-bit golden
semantics; they do not upgrade the G1, reference-fixture, performance,
cross-ISA, or golden authority.

`tests/derham_battery.rs` emits schema-validated `fs-obs` JSONL under suite
`fs-feec/derham`. Post-assert identities are `dd-highorder/m{m}-r{r}`, `dims`,
`commute-1d/m{m}-r{r}`, `commute-3d`, `proj-orders/r{r}`, and
`legendre-mass`. Fixed-input rows carry seed zero; `dd-highorder` carries input
root 8 and records stream kernel `0xDE4A` plus both tile coordinates. Loop rows
remain at their original per-iteration gates. The frozen hash remains an Info
Custom event named `derham-golden` under `derham-golden/measurement`, before
the terminal equality gate, recording root 8, kernel `0xDE4A`, tiles 300–302,
and actual/expected hashes. It is diagnostic replay evidence, not a pre-gate
passing verdict or an upgrade of the existing G1 claims.

`tests/ho_probe.rs`: per-mode convergence regression — the diagnosis
that single-cell symmetric fixtures superconverge at even r (a metric
trap, so MMS ladders start at m ≥ 2).
`tests/perf_lane.rs` (bead cwjn, release + `--ignored` only):
sum-factorized apply throughput vs the fs-roofline MEASURED machine
peak — the ≥30%-of-peak target at p = 4 plus the throughput-vs-p sweep
(r = 1..6), emitted as retained JSON. A positive gate additionally requires
one frozen authority-admitted baseline snapshot; plain or refused inputs are
reported as candidate measurements only.
`tests/vecfam_battery.rs` (bead dcng): vec-001 dims + Euler
alternating sum; vec-002 dof-Kronecker + mass SPD; vec-003
tangential/normal conformity; vec-004 dd = 0 through grad/curl/div;
vec-005 interpolation ladders; vec-006 two-tier G3 relabeling;
vec-007 Whitney cross-checks + frozen golden (cross-ISA row LEDGERED
PENDING, same policy as the simplex golden).

The vector-family battery emits schema-validated `fs-obs` JSONL under suite
`fs-feec/vecfam`. Fixed-input conformance identities are `vec-001`, `vec-002`,
`vec-003`, `vec-004/r{2,3}`, `vec-005/{nedelec,rt}-r{2,3}`, `vec-006`, and
`vec-007`, all with seed zero. The legacy golden diagnostic remains an
info-only `Custom` event named `vec-007-golden` at
`vec-007-golden/measurement`, recording actual and expected hashes before the
unchanged terminal equality gate. All rows retain their original emission
boundaries; this migration adds replay evidence without upgrading the existing
G1/G3, golden, or cross-ISA claims.

`tests/simplex_battery.rs` (slice 3) emits schema-validated `fs-obs` JSONL
under suite `fs-feec/simplex`. Post-gate aggregate identities are `counts`,
`unisolvence`, `conformity`, `mms-simplex-order/r-1` through
`mms-simplex-order/r-4`, `orientation-r2`, and `orientation-r4`. The retained
legacy `orientation-r2` identity covers the existing r = 3 signed-permutation
fixture; its payload records `r = 3` explicitly. Counts,
unisolvence, and the MMS ladders are fixed-input cases with seed zero.
Conformity records input root 11, kernel `0x51E3`, and tile 1. The r = 3
orientation measurement records both the vector stream (root 11, tile 20) and
the permutation stream (root 12, tile 7); its aggregate seed is the vector
root. The r = 4 orientation case records permutation root 12 and tile 9. These
are logical input streams, not execution randomness; this suite has no Cx or
execution seed.

The eight original per-ladder MMS diagnostics remain info-only `Custom` events
at `mms-simplex-order/r-{r}/measurement/m-{m}` after each solve and before its
slope gate: `m = 4, 8` for `r = 1` and `m = 2, 4` for `r = 2..4`.
`conformity/measurement`, `orientation-r2/measurement`, and
`orientation-r4/measurement` retain the numerical diagnostics beside their
post-assert aggregates. `simplex-golden/measurement` retains actual and
expected hashes plus input root 11, kernel `0x51E3`, and tile 40 at the original
info-only boundary before the unchanged equality assertion; no passing golden
aggregate is synthesized. All measurements use object-shaped JSON, encode
non-finite floating values as `null`, and are failure-linted, serialized,
wire-validated, then printed. Aggregate detail retains each legacy text prefix
and appends the applicable replay coordinates.

The unchanged gates cover Duffy weight/volume exactness and P_r dimension
counts (local dofs = C(r+3,3), r = 1..6); unisolvence via mass SPD (r = 1..5);
shared-face conformity to 1e−13; G1 MMS slopes greater than r + 0.6;
signed-permutation operator equivariance at r = 3 to 1e−12 relative; and G3
physics invariance at r = 4 with vertex deviation below 1e−9 and L2 relative
deviation below 1e−6. The canonical events are replay evidence for those
existing fixtures, tolerances, and frozen bits; they do not upgrade G1/G3,
cross-ISA, golden, completeness, or performance authority.

`tests/derham_battery.rs` (slice 2): curl∘grad and div∘curl ≤ 1e−13
relative on four (m, r) fixtures; exact-sequence dimensions
(χ = 1 for m = 1..3 × r = 1..6); commuting diagram 1D
(G·π_C f = π_D f′ ≤ 1e−11) and 3D on product fields; projection
G1 ladders for both 1D families at r = 1..6 (C: order ≥ r + 0.6
gate, measured ≈ r + 1; D: ≥ r − 0.4 gate, measured ≈ r) — these
drive all four 3D tensor space types' rates; Legendre mass closed
form; its own golden hash.
`tests/differential_characters.rs` (RA.2a, G0/G3 schema battery): circle,
torus, relative interval/solid-torus, and Moore-like torsion metadata;
curvature, characteristic, relative-boundary, product, holonomy, and gauge
domain/codomain checks; explicit pending kernel/image obligations; coefficient,
degree, primal/dual, whole-subcomplex/terminal, and budget falsifiers;
dimension-`dim(X) + 1` character coverage; lattice-preserving curvature and
dimensionless-cycle holonomy checks; and deterministic identity replay plus
semantic mutations. This is object-semantics evidence only; the constructive
exactness, torsion algebra, naturality, and refinement checkers are RA.2b/RA.2c
scope.
`tests/terminal_relative.rs` (I13.2a, G0/G3, feature
`terminal-relative`): exact triangle `boundary^2`, malformed incidence and
subcomplex refusals, explicit relative containment, terminal/insulation and
phase/orientation failures, complete PortSchema and presented Machine-IR
identity mutations, declaration-order replay, nominal separation of integral
representative/current/distributed-current/geometry sectors, named conversion
families, real coboundary typing, nonfinite refusal, and flow-coordinate plus
trivialization enforcement; a terminal-cut loop graph checks exact integral
  Stokes pairing and three relative cycles, while disconnected and shared-
  closure two-phase graphs check owned top-cell restriction, explicit phase
  tagging, and the refusal of ambiguous component or conversion bindings.
  Parallel-edge and orientation-reflection relabel fixtures pin signed basis
  transport, exact naturality, inverse/composition, canonical row ordering, and
  fail-closed non-chain/support maps.
  Multiphase/component and driven/return permutation fixtures pin complete
  semantic squares, explicit current-sign compensation, generic-versus-physical
  transport separation, current-times-winding invariance, cell-natural
  distributed-current transport with fresh nominal constraint receipts,
  geometric-coil redeclaration with fresh nominal realization artifacts, and
  declaration-only target physical-map rebuilding for both realization
  families. Sequential-versus-direct map redeclarations use separately
  transported/redeclared target objects and pin pair/phase/sector checks before
  fresh map ID/artifact checks, with invalid twins for each boundary. A fully
  relative interval fixture pins empty quotient-basis algebra and zero winding
  replay, while schema-version, canonical-preimage, relative-support, port-time,
  and presented MachineGraph-version mutations pin the identity boundary. A
  two-dimensional square-CW surrogate exercises full-dimensional component
  closure, opposite codimension-one terminal patches, multi-degree exact
  boundary/coboundary pairing, and terminal/insulation overlap refusal. It is
  cellular schema evidence, not a torus, linking, manifold, or embedding claim.
`tests/integral_topology.rs` (I13.2b tranche 1, G0/G4, feature
`moonshot-integral-topology`) verifies a nontrivial rank-one Smith witness and refuses
off-diagonal, negative, zero-order, divisibility, fake-inverse, transformed-
source, shape, storage, exact-work limit+1, and checked-overflow twins. Every
cancellation poll is injected transactionally; the empty exact matrix remains
representable and every success is tagged abstract-algebra-only.
The tranche-2 tests pin the terminal-cut-loop matrix, pair/phase/component and
basis binding, declaration-order replay, independent boundary/coboundary action,
both unaugmented rectangular edge maps, every budget limit-minus-one,
unknown-phase/excess-degree refusals, and cancellation through final
publication.
Tranche 3 adds a centered cellular-surface chain with a nontrivial verified
kernel transform, a bottom-edge shear proving `V^-1` rather than cell labels,
Smith-source mutation refusal, every resource limit-minus-one, and exhaustive
cancellation through final binding/scalar counts.
Tranche 4 constructs and independently verifies rank-deficient, rectangular,
negative, empty, and divisibility-repair Smith forms, including the repair
cascade `diag(6,10,15) -> diag(1,30,30)`. It pins byte-identical replay, every
data-dependent limit-minus-one, unrepresentable coefficient refusal,
cancellation across construction and final verification, and an exhaustive
625-matrix `2 x 2` oracle over entries `-2..=2` using independent gcd/determinant
invariant-factor identities.
Tranche 5 binds the centered-surface lower image to `[1,1,1,1]` and free rank
one, plus an explicitly abstract finite-CW pair whose lower image
`[[1,1],[-1,1]]` yields presentation factors `[1,2]` and torsion `[2]`.
Same-shape mutation, empty-shape substitution, both resource limit-minus-one
cases, complete authority retention counts, and every publication poll are
pinned.
Tranche 6 lifts the centered free representative and abstract `Z/2` torsion
representative into their original pair-bound chain bases, retains an exact
torsion filling chain, and pins cycle/order equations, output/retention/work
counts, every resource limit-minus-one, and every cancellation point. A second
valid sheared Smith witness distinguishes `U_L^-1` from the incorrect `U_L`
lift and demonstrates explicitly noncanonical but byte-replayable output. A
connected relative `H_0 = 0` fixture preserves the distinct `dim(C_0) x 0`
generator, `dim(C_1) x 0` filling, dual `dim(C_0) x 0` cocycle, and `0 x 0`
period shapes with zero arithmetic work.
Tranche 7 pins the centered free cocycle `[1,-1,0,0,1,1]`, its unit period,
torsion-only `2 x 0` cocycles, a mixed `Z/2 + Z` zero-then-identity period
matrix, and a sign-changed complete Smith witness whose generator and cocycle
both change while their period remains one. Every resource limit-minus-one and
all 14 cancellation polls through final publication are covered.
Tranche 8 pins the centered free generator, an explicit face boundary, and a
boundary-shifted free representative; the abstract `Z/2` generator, its exact
order-two filling, and negative-residue normalization; and a zero-homology
bottom chain with an explicit filling. Strong binding/noncycle refusals, every
resource limit-minus-one, full decomposition replay, and all 26 centered
cancellation polls through final publication are covered.

## Perf-lane observations (bead cwjn: authority-admitted both-ISA gate open)

- BOTH reference ISAs have historical operator-baseline candidate observations
  on the committed register-accumulator/index-loop contraction kernels
  (2026-07-11, plain baseline stores in `perf-baselines/`):
  macos-aarch64 (M4 Pro, fingerprint 80cb534fbaf60b50): p = 4 attainment
  0.404 (20.8 GFLOP/s); linux-x86_64 (5975WX ts1, fingerprint
  614b09f101a1e33b): p = 4 attainment 0.439 (15.6 GFLOP/s), quiet host. These
  numbers remain useful measurements, but their plain stores do not carry the
  authority admission required for a citable claim.
- The x86 journey, for the record: 0.026 (baseline libm-fma calls) →
  0.046 (a55x fma capsule, scalar vfmadds) → 0.091 (register-array
  accumulators) → 0.44–0.59 (const-bound index loops — the shape x86 SLP
  finally packs). Every step bit-identical: the ascending-l per-element
  order is pinned and the sf-kron golden never moved.
- Sweep from the historical ts1 candidate run (18p⁴+3p³ flop model): r = 1: 9.9,
  r = 2: 8.7, r = 3: 20.9, r = 4: 23.9, r = 5: 21.4, r = 6: 22.0 GFLOP/s
  (low orders are gather/scatter-bound — recorded follow-up, not gated).
- The current lane distinguishes comparison from authority. A plain
  `FRANKENSIM_BASELINE_STORE`, absent/partial authority configuration, or any
  denied, revoked, tampered, missing-source, or cross-machine attestation emits
  a report-only receipt and cannot produce a positive gate. A fully configured
  attested store uses `FRANKENSIM_PROMOTION_AUTHORITY_POLICY` plus
  `FRANKENSIM_RETAINED_SOURCE_RECEIPTS`, captures one atomic authority decision,
  and embeds the full frozen pre/post snapshot in the final gate JSON so the
  measured claim cannot be detached from its admission decision. A positive
  gate additionally requires `FRANKENSIM_ROOFLINE_LEDGER`: the lane atomically
  records the exact admission receipt and exact final-gate JSON through the
  shared `fs-roofline` external-gate protocol before it emits
  `citation_eligible:true`. A missing or empty ledger path keeps the completed
  measurement report-only; a ledger write or exact re-read failure fails closed
  and cannot emit a positive gate.
  The citable policy owns its clock: an unavailable clock or an epoch-day
  rollover between mint and post-probe invalidates attested evidence and ends
  the lane as `environment_invalid`; configuration refusals alone remain
  measured report-only observations.
  The retained-source file is a protected hash-inventory declaration; this lane
  does not fetch or independently prove availability of the named bytes.
  Its conformance matrix removes one named receipt to prove that missing source
  evidence stays report-only, then re-endorses the identical baseline under a
  rotated key and proves that the baseline hash stays fixed while the key and
  authority-policy receipt move.
  Contaminated axes
  remain `environment_invalid`, neither pass nor fail (observed live: the first
  ts1 attempt was refused during its own compile burst).

## No-claim boundaries

- CIRCUMCENTRIC diagonal star deliberately absent: on
  non-well-centered meshes (including the Kuhn fixtures) it produces
  negative dual measures, so shipping it without well-centeredness
  machinery would be a certificate without evidence. Follow-up scope
  together with mesh quality certificates.
- Simplicial H¹ now reaches high order (slice 3, MMS-gated at
  r = 1..4; the basis is valid to r = 6 by the unisolvence/count
  gates). Simplicial H(curl)/H(div)/L² families now ship at r = 1..4
  (bead dcng, vecfam battery); r ≥ 5 needs conditioning work on the
  monomial Vandermondes (recorded successor). Full 3D VECTOR MMS
  (curl-curl / mixed Darcy solves) need tfz.10 — canonical
  interpolation ladders stand in, labeled. The tensor-product side
  covers all four space types (slices 1–2). Unstructured-hex
  orientation is later-slice scope; the cross-ISA ≥30%-peak perf gate remains
  open (bead cwjn, see Perf-lane observations). `HexComplex` incidence is still
  not consumed.
- MMS covers the PRIMAL Poisson form; the mixed-form MMS (flux
  variable through M₂/d₂) joins the solver-stack lane (tfz.10) where
  saddle-point solvers live.
- No FrankenNumpy live views (descriptor only — wf9.5 precedent);
  no DWR/goal-oriented estimates (their own bead); no adjoint hooks
  yet (fs-adjoint, tfz.24); no boundary-condition abstraction beyond
  test-side Dirichlet pinning (the BC/trace machinery is CutFEM/
  solver-stack scope).
- `betti` is a fixture-scale certifier, not persistent homology at
  scale (fs-topo).

## No-claim boundaries (cohomology)

- Stars are the DIAGONAL barycentric variant; Galerkin-star (full mass)
  decompositions and LOBPCG spectral harmonic solvers land when the
  eigensolver bead's machinery is consumed here — the projection route
  needs no eigen machinery and is deterministic.
- Γ recovery is G1-graded (midpoint-sampled line integrals; the
  coexact component of a sampled irrotational field is discretization
  noise, not exactly zero) — tolerances are measured and ledgered, not
  asserted tight.
- The thin-airfoil WING benchmark (BEM + Kutta condition) belongs to
  fs-bem-fmm (tfz.20), which consumes this module's circulation
  functional; the fixture here is the cylinder-vortex ring.

## No-claim boundaries (relative differential characters)

- RA.2a fixes finite object and map semantics only. It does not yet construct
  Hopkins--Singer representatives, compute kernels/images, prove either exact
  sequence, implement Smith normal form or torsion linking, verify gauge-orbit
  equality, or establish subdivision/refinement naturality; those are the
  explicit RA.2b/RA.2c gates.
- Cell counts and named relative components are schemas, not a proof that one
  cell set is an actual subcomplex of another. A constructive checker must
  validate incidence closure and exact arithmetic before any exactness status
  can be promoted.
- A `CoefficientProductSchema` records an immutable coefficient-map artifact
  identity and its typed signature, but does not execute that artifact or prove
  associativity, graded commutativity, unit compatibility, or physical
  normalization. Product-law witnesses remain constructive-algebra evidence.
- No smooth/continuum approximation, electromagnetic force, flux
  quantization, material, topology-change, or performance claim follows from
  these types. The deterministic algebra ID detects replay drift but is neither
  collision-resistant nor an authority receipt.

## No-claim boundaries (terminal-relative winding schema)

- I13.2a admits exact relative-cycle representatives; it does not compute or
  identify relative homology classes, boundaries modulo chains, torsion,
  periods, Smith/Hermite normal forms, or subdivision/refinement equivalence.
  General integer cellular incidence beyond the oriented signed-unit `-1/+1`
  lane is not admitted. Signed incidence plus `boundary^2 = 0` is not a
  certificate of simplicial arity, regular-CW attaching maps, manifoldness, or
  geometric embedding.
- Component declarations are full-dimensional physical conductor partitions;
  no connectivity search is executed. Presented Machine-IR digests and keys are
  not proof of graph membership, port ownership, causalization, or executable
  coupling. Declared conversion artifacts are not executed or verified here.
- `Driven` and `ReturnReference` are DC-oriented first-slice roles. No neutral,
  protective-ground, phasor/RMS, phase-sequence, multi-winding permutation,
  polarity-inference, or topology-event semantics are claimed. Terminal-free
  and single-terminal relative pairs are not admitted by this first slice.
- Phase-local support currently admits exactly one conductor component per
  phase. Parallel-path/branch/netlist aggregation requires a later explicit
  phase-to-component-set schema; an integral cochain is a dual algebraic object,
  not a real field, metric dual, cocycle class, or geometric linking witness.
  Components may share lower-dimensional closure cells, and those cell
  generators intentionally occur in each owning phase's separately tagged
  basis; this slice does not construct a tagged-copy/direct-sum cell complex.
- One admitted signed cell relabel does not establish general representation
  naturality. Phase or terminal permutations, current-orientation compensation,
  tree/cotree or cut changes, non-unimodular/partial maps, refinement, remesh,
  topology events, and homology/cohomology induced-map authority require their
  own explicit schemas and evidence.
- An admitted physical relabel covers only its enumerated bijections. It does
  not authenticate a MachineGraph or netlist equivalence, authenticate fresh
  distributed-current constraint receipts or geometric-coil artifacts,
  transport arbitrary real fields or geometric data, infer phase sequence or
  polarity, or establish
  refinement/remesh/cut/topology-event naturality.
- A re-declared conversion map is a nominal target-side declaration, not map
  execution or a naturality certificate. The relabel does not prove that its
  caller-supplied endpoints are mathematical images of the source endpoints,
  authenticate either map artifact, use the cell sign or current-coordinate
  sign to execute a physical law, or prove a realization-square commutes.
- No field transfer, current-density solve, electromagnetic force, material,
  thermal, manufacturability, geometric embedding, cancellation-latency,
  performance, or authority-receipt claim follows from these types.

## No-claim boundaries (integral topology)

- I13.2b tranche 1 verifies a supplied Smith witness. Tranche 4 now computes a
  deterministic Smith form, emits the same complete witness, and passes it back
  through that verifier. It still does not construct Hermite form, quotient
  homology/cohomology, free or torsion generators, periods, linking pairings,
  long exact sequences, or induced maps.
- I13.2b tranche 2 proves faithful extraction of one admitted phase-local
  quotient incidence matrix. Independently Smith-reducing adjacent boundary
  matrices is not a homology proof: the incoming image must first be
  transported into the outgoing map's verified kernel coordinates.
- I13.2b tranche 3 proves that exact incoming incidence lies in one verified
  outgoing kernel and records its witness-relative coordinates. It does not
  itself claim homology, torsion, free generators, periods, linking, long exact
  sequences, naturality, or physical winding. The generic tranche-4 constructor
  can reduce that lower matrix, and tranche 5 now owns and byte-compares both
  values before publishing only the quotient invariant decomposition.
- I13.2b tranche 5 computes free rank and torsion invariant factors of the
  admitted phase-local cellular quotient. Tranche 6 now publishes checked free
  and torsion cycle representatives and torsion filling chains in the original
  pair-bound quotient-chain bases. They are Smith-witness-relative, not
  canonical across valid witnesses. No integer periods, cocycles, linking
  pairings, long exact sequences, induced maps, relabel/remesh naturality,
  embedding, physical realizability, or physical winding follows.
- I13.2b tranche 7 supplies witness-relative integer cocycles dual to the free
  generator sector and an exact integer period matrix. It is not complete
  integral cohomology: the universal-coefficient `Ext` sector, torsion
  characters, `Q/Z` or `R/Z` linking, cup/Poincare pairings, harmonic or de
  Rham representatives, continuum units, naturality, and physical circulation,
  flux, current, winding, or embedding authority remain absent.
- I13.2b tranche 8 classifies one admitted pair-bound integral cycle and retains
  one exact representative-plus-boundary decomposition in the selected Smith
  witness basis. Its torsion residues are canonical only relative to the
  retained ordered invariant factors; its free coordinates and adjustment are
  witness-relative, not minimal, harmonic, geometric, or canonical across
  alternative valid witnesses. It adds no relabel/remesh naturality, induced
  map, linking, continuum-period, embedding, flux, current, or physical-winding
  authority.
- `AbstractAlgebraOnly` is load-bearing. Synthetic CW/Moore/lens-space matrices
  may test the algebra kernel but cannot establish a conductor, terminal,
  material, embedding, winding, flux, force, or machine claim. Physical R3
  applicability requires a separate receipt binding an admitted
  `TerminalRelativePair`, embedding/boundary assumptions, terminal-subspace
  premise, checker identity, budgets, and independent replay.
- Integer inverse matrices prove unimodularity only for the retained exact
  transform. They do not prove that a caller's complex, basis, cut, relabel,
  refinement, remesh, or topology event has the claimed physical semantics.

## Interior Helmholtz acoustics (`acoustics::helmholtz`, bead pexfp.93 M2.4-W1a, first slice)

Home: this crate, as M2.4-W0 prefers (the dependency `fs-feec -> fs-modal` is L3 -> L1; fs-modal's own dev-dependency chain reaches
fs-feec, a dev-only cycle that cargo allows and the layer check exempts).

- `HelmholtzFe::new(complex, positions, order, Medium)`: the hierarchical H1 space of `highorder::SimplexSpace` at order
  `1..=MAX_ORDER` (4), straight-sided tets, stiffness `K` and mass `M`, with typed refusals (invalid medium, order cap and cap + 1,
  bad position table, degenerate tet). `boundary_faces`, `boundary_mass` (`C = int phi_i phi_j dS`) and `boundary_load`
  (`b = int phi_i dS`, area) take boundary faces as sorted or unsorted vertex triples and refuse a face that is not a boundary face of
  the complex. Face quadrature is the collapsed Gauss-Legendre map, exact to degree `2r + 1`.
- `modes(window_hz, opts)`: rigid-wall (Neumann) eigenpairs of `(K, M)` by `fs_modal::slice_window`, so the in-window count is
  inertia-certified and each mode carries its certified residual; a negative low end includes the constant mode; a window holding more
  than `MAX_MODES` (512) refuses instead of truncating.
- `ReducedHelmholtz`: a Galerkin reduction of the driven problem `(K - k^2 M + j k sum_b y_b C_b) p = j k rho c v_in b_src`
  (time convention `e^{+j w t}`; `y_b = rho c / Z_b` the dimensionless specific admittance of a locally reacting boundary; a source port is a
  prescribed inward particle velocity; every other face is rigid) onto the modes plus one M-orthonormalized static vector `(K + M)^-1 b_port`
  per port. Without the static vectors a boundary-forced modal expansion converges only like `1/N`; the test logs both. A frequency is one
  dense complex solve of the reduced dimension (at most `MAX_REDUCED_DIM`, 640). `response` returns the outward volume flow and mean
  pressure of every port; `transfer_to` is the volume-velocity transfer function (outward flow over inflow).
- Tests (in-module, `cargo test -p fs-feec --lib acoustics`): `helmholtz_rigid_box_modes` (a 1.0 x 0.7 x 0.4 m box, ladders r = 2: n = 3, 4,
  6, 8; r = 3: n = 2, 3, 4, 6; r = 4: n = 2, 3, 4, 6; certified count 21 with the constant mode at the finest rung; pre-registered gates: the
  slope of the lowest nonzero frequency error over the finest pair of rungs is at least 2r - 0.2 for every order, and the 20 modes at the
  finest rung of r = 4 are within 1e-4; every rung prints its dofs and times and every pair of rungs its slope before a shortfall fails),
  `helmholtz_duct_driven_transfer_matches_plane_wave_theory` (exact 1-D duct with an impedance end, a ladder over mesh and over
  static vectors), `helmholtz_refusals`, `helmholtz_modes_repeat_bitwise`.
- Not claimed: isoparametric (curved) walls, the closed-cylinder and resonator tests, the p-ladder, viscothermal boundary-layer admittance
  (M2.4-W1b), a radiation solve (the admittance is an input law), multiport extraction, a direct complex solve for cross-checking the
  reduction at large size. The bead stays open.

### Structural boundary traces on the original pressure basis

HelmholtzFe::boundary_linear_load integrates a caller-supplied linear normal
displacement/velocity trace against every original hierarchical pressure
basis function on selected boundary faces. It uses the existing face rule
at order + 1; the product is integrated before reduction. Original vertex
order is mapped explicitly, repeated faces and nonfinite traces refuse,
and boundary_geometry determines outward normals from the owning volume
cell independently of the caller's triangle winding. The original K/M
operators and modal solver remain unchanged.

The public geometry/order/count accessors allow higher-layer structural
couplers to verify exact boundary ownership and budget the projection before
construction. Tests compare linear traces against analytic triangle moments,
the original uniform boundary load, vertex permutations and volume-oriented
normals through every admitted pressure order.


### Isothermal gas diffusion on the original acoustic FE

`HelmholtzFe::isothermal_diffusion_modes` restricts the original consistent
stiffness/mass pencil to the exact hierarchical trace constraints of explicitly
selected wall faces. Other faces retain natural zero heat flux. It returns all
positive modes through the declared unit-diffusivity eigenvalue ceiling, with
original inertia counts/residuals and full-coordinate shapes. Up to `DENSE_FULL_PENCIL_DOFS` (512) free DOFs the whole pencil is solved once, densely, and the ceiling only filters
(`DiffusionModes::dense`): the eigenvalues, shapes and residuals of the modes kept are then bit-identical for every ceiling, which
the material caller needs because its ceiling is `max_decay_rate/diffusivity` and moves with gas pressure. Above that size the
sparse window owner is used and its last bits follow the window. The material
caller owns kappa/(rho*cv), source coupling, memory/work budgets and physical
validation. Duplicate/nonboundary faces, empty thermal spaces, incomplete or
oversized windows and cancellation refuse. The sparse modal owner is reused;
no heat stepper, penalty boundary or copied geometry is introduced. The original
single eigensolver call has its existing cancellation granularity. Box/mixed-wall
modal and normalization regressions supplement the higher-level acoustic tests.

### Bounded complete coexact edge projection

`cohomology::coexact::PreparedCoexact1` selects the original discrete curl's
free-edge columns D and the original positive barycentric edge metric M. It
prepares P = D M^-1 D^T after checking the sparse-product work budget and
removing zero face rows. P is positive semidefinite; D x is a compatible right
hand side. The caller owns the bounded checked solve, true residual and
physical acceptance. Recovery M^-1 D^T w is M-orthogonal to the entire kernel
of D, including harmonic components, without inferring topology from a random
basis. No uniqueness of w or full-space SPD claim is made.

The fs-em vector-potential owner uses this same preparation for complete
current compatibility and flux-preserving gauge recovery. Its six focused
physics tests include a multiply connected source refusal, field refinement,
and bounded solver/cancellation behavior. Existing cohomology interfaces and
algorithms retain their owners.

### Thermoviscous wall layers on the original interior pressure space

`acoustics::helmholtz::boundary_layer` adds the frequency-domain Wentzell
wall form of Berggren, Bernland and Noreland (JCP 371, 2018, 633–650;
https://arxiv.org/abs/1801.04177). With the existing `exp(+i omega t)`
convention, the operator contribution is
`(-1+i)/2 [delta_v S + k^2 (gamma-1) delta_t C]`. `S` integrates the
tangential gradient of the original hierarchical pressure trace using exact
barycentric derivatives and the original face quadrature; `C` is the original
boundary mass. Normal gradients and zero-trace volume bubbles do not create
wall loss. No scalar impedance substitutes for the surface derivative.

`ReducedHelmholtz::new_with_boundary_layers` assembles and projects those
two matrices once into the existing modes/static-enrichment basis. Its
ordinary `response` then includes both wall mechanisms. The lossless
constructor retains its arithmetic, and zero transport gives bitwise-identical
responses in the focused test. `boundary_layer_power` returns separate
viscous and thermal power from the Hermitian quadratic forms, nonnegative up
to roundoff, without clipping or repairing the calculated energy, and
the physical factor `1/(2 omega rho)`. This closes input work together with
the pointwise acoustic-port power, including nonuniform port pressures.

Transport is supplied by the source material owner; density and sound speed
are the original FE medium. Geometry must supply the actual minimum smooth
curvature radius and half-gap, independently of the mesh. The computed
frequency interval enforces `delta/min(radius,half_gap)` and `k delta` below
the caller's declared ratio, itself at most 0.05. This is an asymptotic
admission condition, not a certified bound on impedance, resonances or Q.
Invalid or overlapping/duplicate wall and acoustic-port faces refuse.
The standalone surface assembly also offers cancellation between faces;
reduction and eigensolver cancellation retain their existing granularity.

Tests exercise exact affine triangle energies for every admitted polynomial
order, normal-gradient exclusion and face orientation, layer scaling and band
refusal, full-FE input/work closure, the zero-loss limit, and an axial box-mode
ladder against analytic viscous/thermal surface integrals. The M2.4-W1b bead
remains open: no curved isoparametric map, narrow-gap LNS, molecular wall/gas
coupling, corner-layer resolution, Bessel-cylinder acceptance or final
instrument-accuracy claim is made. Finite wall thermal response is supplied
by the extension below. Patch edges use natural
zero tangential conormal flux. This frequency-domain wall approximation is not
advertised as a passive time-domain realization.

### Finite-conductivity and finite-capacity acoustic walls

`boundary_layer::thermal_wall` extends the original Wentzell wall through
temperature and normal heat-flux continuity. `ThermalWall` selects the original
isothermal condition, an adiabatic thermal wall, or a homogeneous conductive
half-space/finite slab. A slab explicitly selects zero temperature or zero
normal heat flux at its backside. With `q=(1+i)/delta`, its exact normal heat
admittance is `k*q*coth(q*d)` or `k*q*tanh(q*d)`. The gas temperature deficit
`Y_wall/(Y_gas+Y_wall)` multiplies the original thermal pressure form. Density,
heat capacity and directional conductivity come from the material owner.

`ReducedHelmholtz::with_thermal_wall` retains the original pressure basis,
projected surface matrices and viscous tangential form. Ordinary `response`
and `boundary_layer_power` include the new thermal response. The default
isothermal operation order remains unchanged. The thermal power is total
gas-plus-solid diffusion loss; `ThermalWallResponse` exposes their separate
nonnegative entropy-production factors and complex interface temperature.
The imaginary operator coefficient uses these positive forms directly,
avoiding cancellation of leading reactive terms for thin insulated walls.
No loss is clipped. Scaled admittance division and small-depth analytic
series avoid overflow/cancellation without introducing a fitted damping law.

The original gas band remains enforced. Solid penetration depth, capped at
finite wall thickness, must also be small against the real curvature radius,
explicit minimum tangential variation length and acoustic wavelength. The
declared tangential length includes pressure, material, thickness and backing
variation, independently of face size. This is asymptotic applicability, not
a certified observable-error bound. Properties are frozen at one common
gas/solid equilibrium temperature; the selected heat capacity constraint and
normal conductivity direction remain the caller's responsibility.

Focused tests independently solve the solid normal heat equation with P1
elements on 32/64/128 cells for both backside conditions, requiring order
1.8--2.2 and finest relative error below 2e-5. They also cover half-space and
thin/thick limits, gas/solid entropy partition, very small positive insulated
wall loss, full original FE input/work closure, unchanged isothermal bits and
invalid/unresolved-profile refusal. No molecular relaxation, thermoelastic
heating, imperfect interface contact, lateral conduction, multilayer solid,
mean thermal evolution or time-domain realization is supplied. Existing
curved-map, narrow-gap, corner-layer and Bessel-cylinder limitations remain;
M2.4-W1b is not closed by this extension.


## Spatial work-conjugate interior ports

`acoustics::helmholtz::multiport::InteriorMultiport` projects real piecewise-linear
normal-velocity shapes onto the original hierarchical pressure FE. A shape psi
on physical area A defines inward velocity psi*Q/A and conjugate pressure
P=integral(psi*p)/A; mean acoustic input work is Re(sum(conj(P)*Q))/2. The
uniform shape recovers ordinary volume flow and mean pressure. A zero-mean
shape retains its reactive/transverse response; its Q coordinate is not net
volume flow. Shapes may share faces and need not be orthogonal. Their ordering,
normalization and physical interpretation remain explicit caller data.

Original boundary trace quadrature, volume mass/stiffness, modal basis and
residual-flexibility enrichment are reused; each enrichment uses the actual
spatial load. `prepare` retains one original cancellable complex LU for repeated
simultaneous flow excitations. `impedance` extracts every matrix column from
that factorization without symmetrizing. `scattering` adds the work-conjugate
Robin load directly before factorization, avoiding an intermediate impedance
pole. With positive real reference R, a=(P+R*Q)/(2*sqrt(R)) and
b=(P-R*Q)/(2*sqrt(R)); input work is (norm(a)^2-norm(b)^2)/2.
The frequency, port and reduced-dimension caps refuse invalid work. Cancellation
is checked around projection, between modes, and in the original LU/substitution;
no partially prepared factor or extracted matrix is returned. The existing
volume assembly/static-enrichment cancellation granularity is unchanged.

`new_with_boundary_layers` and `with_thermal_wall` attach the same viscous and
gas/solid thermal operators to these spatial ports. Port and wall faces cannot
overlap. The original `boundary_layer_power` closes their work using the same
solved pressure coordinates. No new damping law or duplicate heat exchange is
introduced. Cached projection and one factorization per frequency are algorithmic
properties, not a measured runtime guarantee.

Four focused native tests cover an independent triangle work moment, parity
with the original uniform-flow FE, multimode reciprocity and input/wall work,
power scattering and replay, invalid-input/late-cancellation refusal, and an
independent straight rectangular guide with a retained below-cutoff cosine mode.
For P2 meshes with 27/125/343 DOFs, plane-impedance relative errors decrease
1.755559e-4 -> 1.832614e-5 -> 4.180726e-6, while evanescent errors decrease
0.2690877 -> 0.03035866 -> 0.005153092. These are fixture-specific FE/trace
convergence results, not general accuracy bounds.

The original exp(+i omega t) convention is retained. Real-reference power
coordinates are not propagating/evanescent traveling-wave amplitudes, whose
normalization and reciprocity metric differ. No evanescent channel is discarded
or incorrectly assigned a Euclidean traveling-wave passivity test. Analytic
cross-section eigenfunctions, curved geometry, resolved thin-gap LNS, transverse
truncation bounds, port de-embedding, qualified tone-hole junctions and complete
M2.4-W2 acceptance remain outside this extraction capability. Spatial traces,
volume mesh/order and reduced basis require their own refinement.

### Joint molecular calorics and resolved spatial heat conduction

`acoustics::helmholtz::molecular::CaloricInterior` reuses the original pressure
FE, consistent volume matrices, spatial port loads and cancellable multiport
factor/loaded-scattering owners. Temperature occupies every original FE
coordinate except exact isothermal wall trace constraints. The supplied
pressure basis remains explicit; automatic static enrichment is not applied.
With exp(+i omega t), the material supplies passive complex cp with Im(cp)<=0:

`[kappa K_T + i omega rho cp M_T + Y_wall C_T] T = i omega H p`;
`A_p = K - omega^2 rho/p0 M + omega^2 rho/T0 H^T T`.

One heat factorization eliminates temperature and one original pressure
factorization serves all simultaneous modal flow excitations per frequency.
Direct loaded power scattering avoids forming a closed-cavity impedance pole.
The pressure/temperature coupling has the same exact discrete work pairing;
the thermal loss is
`[kappa T* K_T T - omega rho Im(cp) T* M_T T + Re(Y_wall) T* C_T T]/(2 T0)`.
Molecular, gas-conduction and solid-conduction terms are separately positive
forms, with no clipped loss and no second complex bulk-modulus correction.
The response exposes actual original-coordinate temperatures, pressure/flow
coordinates and the separate loss terms. Optional viscosity uses only the
original tangential Wentzell form, with its original thin-viscous-layer band;
the original approximate gas thermal wall term is disabled in this lane.

Thermal walls are isothermal, adiabatic, or a source-defined conductive slab or
half-space. `ConductiveThermalWall::normal_admittance` reuses the existing slab
series and solid normal-profile admission; it does not impose a thin gas
thermal layer. Other faces, including acoustic mouths, have zero heat flux.
Gas and solid share a fixed uniform equilibrium temperature. Thermal dimension
is explicitly capped at 640; pressure/port caps and original LU cancellation
remain in force. Assembly polls between wall faces and pressure modes; dense
quadratic forms are bounded by the thermal cap. No partial prepared result is
published after cancellation or a material/numerical refusal.

Four focused tests cover independent uniform molecular compression, exact
input/work closure, actual spatial modes and reciprocal real-reference power
scattering, source/cp/cancellation refusal, and a P2 1/2/3 heat-mesh ladder
against an independent analytic slab with complex molecular cp. Isothermal
relative reaction errors decrease 0.00464670 -> 0.00197543 -> 0.000519832;
finite insulated solid errors decrease 0.000254508 -> 0.000105061 -> 0.0000275360.
The slab uses a declared uniform-pressure Galerkin coordinate; these numbers
qualify that thermal fixture, not arbitrary three-dimensional wave accuracy.

Overlapping thermal layers are represented when the temperature mesh resolves
them. Viscous overlapping layers still require the resolved LNS owner. Volume
mesh/order, thermal penetration depth, pressure basis and port traces require
independent refinement. No curved isoparametric geometry, corner-layer
acceptance, mean flow/heating, evolving composition, lateral solid conduction,
time-domain realization, real-instrument validation or runtime performance
benchmark follows from this frequency-domain extension.

### Resolved three-dimensional thermoviscous acoustics

`acoustics::thermoviscous::ThermoviscousFe` supplies the mixed harmonic
linearized Navier–Stokes owner for a uniform stationary ideal gas on straight
tetrahedra. Velocity uses three original P2 H1 components, pressure uses P1,
and temperature uses P2. The element pair is fixed; an unstabilized equal-order
pair is not offered. Axial diffusion and full deviatoric/dilatational stresses
remain in the volume equations. No thin gas-layer approximation is imposed.

Using exp(+i omega t) and theta=T'/T0, the system has the form
`R + G + i omega H`, with exactly antisymmetric pressure/divergence coupling,
consistent inertial/caloric mass, and physical viscous/thermal dissipation.
Source-evaluated complex cp enters the temperature equation itself; Im(cp)<=0
adds molecular dissipation. Independent bulk viscosity remains explicit and
must not duplicate the same molecular reservoir loss.

Every boundary face belongs to a stationary no-slip wall or a pressure port.
An explicit wall subset is isothermal; other walls and all pressure ports have
zero conductive heat flux. Existing piecewise-linear `BoundaryMode` traces
define total traction `-P psi n` and inward flow `Q=-integral(psi v.n)`.
Port work is `Re(sum(conj(P) Q))/2`. P includes viscous normal stress and is
not a point sample of thermodynamic pressure. Uniform traces have net-volume
flow semantics; nonuniform traction modes require a matching dual map before
connection to a prescribed-velocity modal basis.

Static assembly is reused across frequencies. One original cancellable complex
LU serves repeated drives and every admittance column. Direct real-reference
power scattering adds `sum R_m B_m B_m^T` before factorization and avoids
inverting an intermediate impedance/admittance. At most two correction solves
reuse that LU, with the original `fs_math::Dd` owner evaluating residuals of
the unscaled physical equations. This corrects loss of small flow/work digits
near the static limit without another dense matrix. Final componentwise
backward error above 5e-13 refuses; this is a numerical solve criterion, not
a physical accuracy bound. Cancellation is polled during element assembly,
scaling, residual evaluation, original factor/substitution, and field work.

Limits include explicit element and tetrahedron-aspect bounds, a default
1024/hard 1536 mixed-unknown cap, and `48*ndof^2` bytes for peak dense arrays
(geometry and vectors are extra). The aspect proxy is longest edge divided
by smallest altitude; admission does not prove mesh robustness or convergence.
The response retains original-coordinate pressure, velocity and temperature,
all port flows, and independently integrated shear, bulk, conduction and
caloric work. Losses and work residuals are not clipped.

`PreparedThermoviscous::element_dissipation` and the same method on the loaded
`PreparedThermoviscousScattering` independently integrate the declared full P2
velocity and actual temperature-perturbation coefficients in original tetrahedron
order. They retain the immutable prepared frequency, complex Cp, equilibrium,
transport and affine/quadratic geometry. Each `ThermoviscousElementPower` reports
actual quadrature volume and original `rho*V` mass, shear and independent bulk
viscous heating, conductive entropy production, conductive availability loss,
and complex-Cp caloric availability loss. The same physical quadrature kernel
serves existing response power; its original global point accumulation order is
retained independently of the element totals. Ordinary responses do not allocate
an additional element result vector or refactor the mixed system.

Only `viscous_heat_w()` (shear plus bulk) is a direct mean caloric heat source.
Fourier conduction contributes `T0*entropy production` to the acoustic
disturbance-energy balance; complex-Cp caloric loss is likewise reported as
availability destruction. Neither is silently deposited into a mean heat
equation without its own supported mean-energy closure. Pressure work, boundary
conduction, structural losses and balance residuals are excluded from these
volume records. Full native coefficient lengths and finite values are admitted;
original element/quadrature limits and cancellation apply. An arbitrary admitted
field observation does not certify the field satisfies the mixed equations or
boundary conditions. No loss clipping, fitted attenuation, or vertex-only
gradient reconstruction is used.

`ThermoviscousFe::prepare_stress_at` retains the original P2 physical velocity
gradients and P1 pressure interpolation at an admitted point of an affine or
quadratic tetrahedron. Warm harmonic and instantaneous observations read ten
velocity and four pressure coefficients without allocation or a solve. The
result separates `2 mu dev(sym(grad v))`, independent `zeta div(v)`, and
thermodynamic pressure. Physical tensor shear uses xx/yy/zz/xy/yz/zx ordering.
Harmonic heating includes the cycle-average 1/2; instantaneous heating does not.
Traction normalizes a finite nonzero normal direction and computes tangential
stress from the deviator alone, so a large pressure cannot erase small wall
shear. Rigid spin and isotropic expansion produce zero shear. Selected-element
stress traces are unsmoothed and inherit spatial discretization error; these
point observations do not replace the original weak moving-wall reactions or
constitute an equilibrated interface stress or mesh-accuracy certificate.

`ThermoviscousFe::prepare_field_at` additionally retains the original P2 basis
values and thermal gradients. Harmonic and instantaneous point observations
return pressure, velocity, temperature perturbation, stress and Fourier heat
current without allocation, another solve, or a plane-wave impedance assumption.
The material-bound `ThermoviscousInterior` exposes the same preparation. Signed
disturbance-energy flux separates pressure work, viscous work and conductive
available-energy transport: `J=p*v-tau*v+(T'/T0)*q`, `q=-kappa*grad(T')`, with
`Re(a*conj(b))/2` for peak harmonic fields. This is a uniform-background linear
fluid energy corollary, not first-order heat flux or a streaming calculation.
These selected-element traces retain the original spatial error and are not
substitutes for the weak terminal or moving-wall work coordinates.

Six focused point-field tests and six retained stress tests pass in a native
exact-source leaf: hierarchical field/gradient interpolation, harmonic versus
instantaneous period averaging, Couette and conduction dissipation identities,
standing-wave quadrature, signed flux projection, and input/cancellation refusal.
Parent FE assembly carriers in that leaf are compile-only nominal types; no FE
factory or solve is claimed executed there. Existing full-FE harmonic loss and
manufactured curved-geometry fixtures also exercise the new factory, but those
two augmented regressions await the full dependency build.

Four focused native tests cover traction-work closure and all live loss
channels, reciprocity, simultaneous drives, direct power scattering, replay,
boundary/material/budget/cancellation refusals, and independent continuum
checks. Low-frequency compliance relative errors are 3.33714e-6 isothermal,
1.87563e-7 adiabatic and 1.80514e-7 relaxing, within the unchanged 2e-5 gate.
A long finite rectangular duct on 126/332/638 mixed unknowns has Poiseuille
conductance errors 0.0584342 -> 0.0161579 -> 0.00611767. Finite-end effects
remain separate from the fully developed reference. The finest fixture error
is about 0.612%, not the M2.4-W1c infinite-slit/tube 1e-4 qualification.

M2.4-W1c remains open: dynamic slit propagation and characteristic impedance,
circular-tube Bessel and asymptotic overlap acceptance, automatic conforming
Helmholtz/LNS interface assembly, curved geometry, finite conductive solid walls
in this owner, moving walls, mean flow, nonlinear acoustics, a sequential LNS
variant, time stepping and large-domain iterative solves are not supplied.


## Resolved thermoviscous gas with finite conductive solid patches

The original `acoustics::thermoviscous` P2 velocity/P1 pressure/P2 temperature owner now accepts `wall::ConductiveWallRegion` through `ThermoviscousFe::with_conductive_walls`. Each region retains original no-slip faces, the original `ConductiveThermalWall`, solid profile geometry and common mean temperature. Multiple solids may have different conductivity, capacity, depth and backing. Patches must be disjoint from ports, exact isothermal faces and each other. Unassigned walls remain insulated; replacing the numeric patch set by an empty set restores insulation. Source material binding belongs to fs-couple.

For the original exp(+i omega t) convention and theta=T'/T0, the thermal weak term is `T0 Y_w(omega) C_w theta`. The original normal solid profile evaluates Y_w once per region/frequency. Sparse original P2 trace masses retain region identity; no extra dense matrix or gas boundary-layer approximation is added. The same equilibrated LU and double-double unscaled residual include the solid term. Independent pointwise integration reports `solid_conductive_w = Re(Y_w)/(2 T0) integral |T'|^2 dS`, separately from gas conduction, shear/bulk viscosity and molecular calorics; no loss clipping.

`PreparedThermoviscousScattering` and `PreparedInteriorScattering` retain original real-reference loaded factors for arbitrary simultaneous incident waves and actual original fields. They avoid inversion of an unloaded impedance/admittance at its poles. Existing scattering entry points reuse these owners. Original physical face coordinates, outward normals and spatial shape values are retained internally for subsequent verified trace matching.

Face partition/common-temperature admission occurs at construction. Original solid conductivity/capacity/depth/profile admission occurs at the requested frequency. Cancellation remains bounded at face, quadrature, row and original solver boundaries. Dense limits are unchanged; trace storage is sparse and bounded by the original admitted mesh.

Focused native fallback: three new `thermoviscous::tests::conductive` cases passed, with four retained full-LNS and four retained spatial multiport cases. The full-LNS thermal slab ladder with 2/4/8 elements gives relative compliance errors 2.1337053211184807e-4 / 1.2150913899507234e-5 / 1.0970749045440163e-6 for an isothermal-backed finite solid, and 2.2647154985204132e-4 / 1.4463910901272932e-5 / 1.9246113513980495e-6 for an insulated-backed solid. Finest gate 2e-4; work closure 2e-8; reciprocity/loaded wave work 1e-9. These are low-kL thermal slab results with finite-length/traction-stress defect retained, not the infinite-slit propagation acceptance. Rust 1.100.0-nightly (2026-08-19) original-source focused seam; pinned toolchain and full constellation gate unavailable.

Scope: stationary ideal gas, fixed common mean temperature, perfect thermal contact and original local normal homogeneous solid conduction. No lateral/interpatch conduction, multilayer solid mesh, contact resistance, evolving mean temperature, moving wall or time-domain realization is inferred. IRP M2.4-W1c remains open for its separate propagation and interface qualifications.


## Hybrid pressure / resolved thermoviscous spatial regions

`acoustics::hybrid::HybridInterior` connects original `InteriorMultiport` and `ThermoviscousFe` owners at explicitly selected spatial modes. The retained map is `T[l,h] = integral(psi_L[l] psi_H[h])/A_H[h]`, so `P_H=T^T P_L` and `Q_L=-T Q_H`. It is computed from original immutable physical trace metadata, including actual coordinates, shape values, normals and areas. Coincident triangulation, opposite normals, unique selected indices, equal retained rank and nonsingular pairing are required. Vertex numbering/winding may differ. Same-area substituted shapes cannot supply a forged map. All unselected coordinates remain external in H-then-L original order.

A prepared hybrid retains the two original loaded volume factors and one small interface factor, bounded by the original 64-mode-per-owner caps. Full off-diagonal scattering and simultaneous original pressure/velocity/temperature fields reuse these factors; no combined volume matrix or intermediate impedance/admittance inversion is introduced. `HybridPower` reports independent physical wall/volume losses plus both signed interface work terms, and responses retain pressure/flow continuity residuals. No loss clipping or passivity repair. Original cancellation, resource and typed solver refusals remain intact.

Three focused native tests pass: original nonorthogonal trace map and fields with finite-solid losses, independent mixed Z/Y network closure, interface pressure/flow/work, reciprocal full scattering, internal-reference invariance and replay; a three-rung H propagation ladder under the actual retained LNS load; and geometry/rank/input/cancellation refusals. The H axial 2/4/8-cell ladder (45/81/153 pressure DOFs) gives maximum pressure/flow relative errors 2.601847292052690e-3 / 1.439797900435083e-4 / 7.393154502284279e-6 against independent analytic duct propagation loaded by original LNS admittance. Finest gate 1e-4; interface coordinate/work gates 1e-10; independent total work gate 2e-8. The ladder isolates H propagation and does not certify the LNS continuum solution. Native fallback uses the same original-source seam and Rust 1.100.0-nightly described above, not a full workspace or pinned-toolchain gate.

This is modal Galerkin total-normal-traction/flow continuity, not full pointwise conforming P2 trace continuity. Nonmatching mortar meshes refuse. The LNS thermal interface is insulated; cuts through important thermal/viscous layers retain a separate physical model error. Interface placement, retained-mode convergence and geometric refinement remain explicit obligations. The default preparation's cp parameter affects only LNS; the molecular preparation below supplies the source's complex K to the pressure region as well. These conditions agree with the normal-stress/normal-acceleration and adiabatic assumptions described in COMSOL 6.3, Acoustic–Thermoviscous Acoustic Boundary (https://doc.comsol.com/6.3/doc/com.comsol.help.aco/aco_ug_multiphysics_couplings.14.12.html). IRP M2.4-W1c/.95 remains open for its full conforming and slit/tube/overlap qualifications.


### Original distributed momentum sources across the hybrid

`PreparedHybrid::response_with_volume_sources` carries the existing
`ThermoviscousBodyForceSource` through both original loaded regions. Force
amplitudes [N] multiply the source's original signed Cartesian density [m^-3],
including zero-resultant distributions and the complete moving-wall lifting.
External incident waves, structural forces and volume/surface heat sources are
coherent drives of the same reciprocal interface closure.

The existing native combined-source response supplies one affine outgoing offset;
the retained interface factor closes its return waves, then each original owner
recovers its full field. Source shapes, material coefficients and passive
scattering are unchanged during warm response. Original body-force work enters
`HybridPower::source_input_w` once, alongside structural and thermal source work.
The native fluid reaction retains the physical body load. Signed source work,
external port work, internal transfer and every physical loss remain distinct;
no scattering norm defect supplies heat or fitted damping.

The existing wall/heat-only method delegates with an empty force slice. Exact
source identity, force-count limits, finite amplitudes and cancellation retain
their original admission. Two new regression cases and the extended original
moving-wall case specify independent physical Z/Y closure on nonorthogonal
traces, complete fields, signed source/emission work, coherent forcing, internal
reference invariance, sealed excitation, binding and replay. These definitions
have not been compiled or executed in the disconnected environment.

This extension retains the existing modal interface's spatial and thermal
limitations. It does not provide full pointwise P2 Helmholtz/LNS trace continuity,
nonmatching mortar, mean-flow sound generation or an inferred body-force law.

### Original pressure-region sources across the hybrid

Complete-coordinate `InteriorMultiport` models now prepare the existing
`VolumeVelocitySource` and `BoundaryVelocitySource` owners on their original
`HelmholtzFe`. The same volume/surface quadrature, signed six-node motion and
work-conjugate observations are reused. A private retained allocation identity
admits that immutable geometry even after movement and refuses a geometrically
identical substitute. Source identity additionally binds the exact multiport;
reduced-basis models cannot reinterpret these full-coordinate loads.

Volume amplitudes are peak m^3/s multiplying signed densities [m^-3]; quadratic
boundary amplitudes are peak inward m/s. The unscaled pressure equation receives
`i omega rho b`. All source columns add coherently before the retained loaded or
unloaded LU solve. Source preparation is not repeated by warm responses. Boundary
support must be disjoint from all port/interface and boundary-layer faces, even
for zero shapes/amplitudes. Original source-count, quadrature, finite-value and
cancellation admission remains in force.

`PreparedHybrid::response_with_all_sources` closes both original affine outgoing
source offsets with the existing interface factor, while keeping LNS body force,
structural and volume/surface heat drives live. Native pressure responses retain
separate signed boundary/volume source-power subtotals; `input_power_w` remains
port work. The hybrid counts those source subtotals once alongside original LNS
sources. Independent viscous, conductive, caloric, wall and bulk losses retain
their original meanings; source work is never deposited as heat.

Three focused regression definitions cover exact full-P2 curved-volume
compression/work and surface moments, signed source fields against the original
heterogeneous equations, and bilateral nonorthogonal hybrid closure assembled
independently in physical pressure/flow coordinates. They include source emission,
coherent interference, internal-reference invariance, original identity/support
refusals, count limits, cancellation and replay. Existing array-valued pressure
checks in the source fixtures were corrected to use their actual [real, imaginary]
representation. These changes have not been formatted, compiled or executed in
the disconnected environment.

The existing selected-modal interface, insulated thermal cut, pressure-coordinate
budgets and physical placement limits remain. This extension supplies prescribed
first-order acoustic injection; it does not infer mass/heat-release production,
mean-flow sound or complete pointwise P2 Helmholtz/LNS continuity.

## Molecular compression across both hybrid acoustic regions

The original pressure multiport now offers `prepare_scattering_with_bulk_modulus(frequency_hz, K, references, cancelled)` with finite Re(K)>0 and Im(K)>=0 in exp(+i omega t). It assembles `K_FE - rho omega^2/K M_FE` on the original pressure basis and geometric mass, reusing the same spatial port loads and loaded LU. The original real-medium preparation keeps its existing arithmetic. `PreparedInteriorScattering::bulk_dissipation_w` independently evaluates `-omega Im(1/K) p* M_FE p / 2`, checks cancellation per row, and adds neither a dense matrix nor a fitted attenuation/loss clip.

`HybridInterior::prepare_with_bulk_modulus` supplies complex K to that pressure region and complex cp to the original resolved LNS heat equation. `HybridPower::helmholtz_bulk_w` remains separate from LNS caloric, conductive, viscous and finite-solid work, and enters the global physical balance. No second K-based loss is added inside LNS. The source adapter evaluates the original molecular material once for both coefficients. The standalone K API admits a frequency sample, not a broadband causal realization.

An existing H boundary-layer configuration refuses this complex-K path: its equilibrium thermal correction has not been derived as the joint molecular/thermal response. The original spatial pressure/temperature `CaloricInterior` remains the separate joint owner. Resolved LNS finite solids remain supported. This extension otherwise retains the hybrid trace, insulated-interface, mode-truncation and placement qualifications above.

Three new pressure tests pass: exact complex closed-volume compression/work in an explicitly declared constant-pressure Galerkin coordinate; full P2 duct propagation versus independently evaluated complex-k analytic pressure/flow; and passive-phase/material/frequency/cancellation/unsupported-wall refusals. The P2 mesh 1/2/3 ladder (27/125/343 pressure DOFs) gives maximum field errors 7.483059143923117e-4 / 6.720372378422063e-5 / 1.496891498044034e-5, finest gate 3e-5. Original real-medium parity, reciprocity, loaded-wave work and positive volume dissipation pass. The nonorthogonal hybrid case additionally keeps H bulk, LNS caloric and finite-solid heat losses live and jointly balanced.

Final focused native verification on the completed extensions: 25 distinct cases PASS across selected invocations, 14 new and 11 retained. Pressure multiport 7; resolved LNS/wall 7; FE hybrid 3; source hybrid 3; source finite walls 2; retained resolved source/duct 3. All newly added cases executed. Rust 1.100.0-nightly (2026-08-19); FE tests use the original-source Cargo seam. A redundant complete material rebuild was killed by memory pressure; coupled tests then compiled current FE/source files against the already verified complete original dependency artifacts from the same compiler, with coherent split-metadata/rlib pairs. No numerical substitutes or production dependency changes. DSR/RCH, the pinned toolchain and the canonical sibling constellation remain unavailable; this is not a full-workspace or instrument-render qualification.

## External surface heat forcing on original wall traces

`ThermoviscousFe::surface_heat_source` projects a signed original P2 surface
density [m^-2] through the same physical boundary mass used by conductive walls.
A drive [W] supplies inward heat flux [W/m^2]. All global hierarchical
coefficients, including edge bubbles and zero-net spatial shapes, survive;
midpoint samples are not substituted for hierarchical coefficients. Original
curved geometry and declared surface quadrature supply the physical measure.
The dimensionless shape integral uses the full vertex partition of unity,
independently of eliminated temperature coordinates; no area normalization occurs.

The source is the existing `ThermoviscousHeatSource` bound to the exact FE owner.
Harmonic, loaded-scattering, conforming-interface and causal source paths reuse
its cached weak load, the existing complete mixed operator and original signed
`theta*q` work, with `theta=T'/T0`. Harmonic peak work has the factor one-half;
causal work integrates the actual trajectory. Neither acoustic available-energy
input nor first-order supplied heat is inferred from a residual or deposited
as mean caloric heat. Original volume, solid and structural losses stay separate.

Faces must be unique original no-slip walls and may carry the existing full
structural motion. Exact isothermal faces, completely constrained temperature
traces, nonfinite coefficients and nonzero support outside the selection refuse.
On a conductive wall the source supplies additional heat at the common fluid/solid
interface; the original passive solid and its temperatures remain. This does not
prescribe a remote backing temperature, infer a heat-production law, evolve the
mean equilibrium or move the reference geometry. Per-face assembly bounds
cancellation latency and temporary storage; no source sampling enters warm solves.

Five focused regressions specify independent P2 and curved surface moments,
signed Fourier flux and absolute work on three original mesh rungs, thermally
generated sound under independently closed acoustic loads, full molecular/solid/
moving-wall trajectory subdivision, source binding and transactional refusal.

## Finite solid heat fields in resolved thermoviscous time evolution

`ThermoviscousFe::prepare_time_with_thermal_fields` retains the actual finite
slabs in the original full-field quadratic time owner. Each source region has
an explicitly requested P1 mesh through its physical thickness, tensorized with
every original P2 surface coefficient. The interface temperature is the fluid
coordinate itself. Consistent wall mass includes interface/interior cross terms;
the wall conductivity and volumetric capacity supply the physical Fourier
operator and storage. The original isothermal or insulated backing is retained.
The same solve advances fluid, molecular temperatures, finite-wall temperatures
and structural coordinates, with a separate integrated solid-conduction channel.

The complete state budget includes every added coordinate; no surface modes or
material reservoirs are discarded. The full slab must satisfy the original
normal-profile approximation down to DC, and the clock must sample the declared
upper band. These admissions do not certify depth-mesh or waveform resolution.
Refining the depth mesh converges to the original harmonic slab admittance;
finite elements are not an exact finite-state replacement for that continuum.
Unmeshed finite walls and semi-infinite walls refuse. The source adapter
`prepare_time_with_wall_diffusion` retains the original wall samples and fluid
calorics. Mean-temperature evolution, lateral wall conduction, changing geometry
and wall thermoelastic expansion remain outside this fixed-equilibrium model.


## Causal connection of complete resolved thermoviscous regions

`PreparedThermoviscousTime::connect` joins two original resolved time owners
through their retained spatial P1 traction traces. Construction checks coincident
physical triangulations, opposite normals, both full trace ranks and equality of
the actual retained spaces at every triangle vertex. The original trace Gram
and cross-Gram determine C with `psi_left C=psi_right`; equal surface areas
alone do not establish a connection. Original inward-flow signs give
`p_left=C lambda`, `p_right=lambda` and
`C^T Q_left+Q_right=0`.

The shared `CommonFlowProjection` constrains the unchanged block Hamiltonian,
and the existing `PreparedQuadraticFlow` evolves all original velocity,
pressure, temperature, molecular, finite-solid and structural coordinates
together. Unselected original pressure modes and every original structural
force remain external drives. The existing source/cavity preparation methods
already return this time owner and inherit the connection directly.

Each region retains its original five physical heat forms. Signed interface
work is independently integrated from the same instantaneous pressure reaction
and original flow, including the reaction's held-input term. It is not inferred
from an energy remainder. Original per-region work balances, internal transfer
cancellation and the enclosing pair balance are checked separately; initial or
endpoint constraint inconsistency refuses without changing physical energy.
The whole coupled state and all observations count against explicit preparation
limits. Rejected or cancelled steps publish neither region.

This is a fixed-equilibrium Galerkin normal-traction/normal-flow interface with
insulated thermal cuts and free tangential traction. It does not impose full
P2 shear/temperature continuity, advect composition, mix gases or certify
interface-placement/trace-refinement convergence. Complete production and
unit-test code type-checks with genuine dependencies. Four full-field tests
cover the original constrained equations, transmission with complete caloric
fields and independent work, time subdivision, and transactional refusals;
their native FE execution remains unperformed in this environment.

### Complete conforming resolved-fluid harmonic interface

`PreparedThermoviscous::connect_conforming` couples the complete original P2
Cartesian velocity and translational-temperature traces at one frequency.
It shares geometric matching, original surface mass, active-coordinate selection
and physical weak loads with the resolved time interface below. The original
P1 pressure fields remain independent; balancing total vector traction carries
both tangential shear components. Inward conductive flux is dual to `T'/T0`.

Preparation builds and equilibrates one bounded interface Schur factor from
the retained, physically corrected region LUs. It holds only one volume response
column at a time. Subsequent independently phased pressure, structural-force,
body-force and heat-source drives reuse those factors, including complex heat
capacity and disjoint moving/conductive walls. Corrections check the original
physical trace equations, with separate velocity and temperature scales that
retain the free-field reference when a stationary constraint arrests a trace.

Each returned region retains its original volume and wall loss quadrature.
Signed mechanical and conductive interface work are reported separately and
cancel across the pair; neither is dissipative heat. The external power balance
excludes internal transfers. A failed continuity/work residual, cancellation,
rank loss, overlapping moving trace, mismatched triangulation/equilibrium or
exceeded trace/dense budget refuses the candidate. This is a fixed-geometry,
conforming resolved-fluid connection; it supplies no mesh, mean-flow or material
accuracy certificate.

Four executed tests of the actual Schur module and original complex LU cover
analytical impedances/work, independent primal vector/thermal solves over
sixteen orders of trace-basis scaling, one-sided stationary constraints and
cancellation/reuse. Native tetrahedral regressions additionally specify
permuted geometry, different fluids, complex heat capacity, phased body/heat
sources and original primal-coefficient identification; their full FE execution
remains unperformed in this environment.

### Complete conforming resolved-fluid time interface

`PreparedThermoviscousTime::connect_conforming` joins matching original boundary
triangulations through complete hierarchical P2 Cartesian velocity and
translational-temperature traces. Geometric vertex/edge identities match across
different node numbering; opposite original outward normals and the complete
original P2 surface mass must agree. Every original pressure mode touching the
cut is removed from the external inputs and must lie entirely on that cut.
Unselected pressures and all original structural forces remain explicit inputs.

The original checked normalized mixed-mass inverse maps actual P2 Cartesian
traction and conductive-flux loads to the existing common-flow constraint. The
same original pair owner, `CommonFlowProjection` and `PreparedQuadraticFlow`
advance both complete velocity/pressure/temperature fields, all molecular
reservoirs, finite solid fields and structural coordinates. Pressure remains in
each original P1 volume space: total viscous traction is balanced; pressure
continuity is not an additional constraint. No second integrator, artificial
interface resistance, fitted trace or model reduction enters.

Traction [Pa] is work-conjugate to area-integrated Cartesian velocity. Inward
conductive flux [W/m2] is conjugate to area-integrated `theta=T'/T0`; the original
thermal storage already carries `T0`, so the flux input has no additional `T0`.
Independent same-trajectory quadratics integrate mechanical `t.v` and thermal
`theta*q` interface work separately. Both signed transfers cancel across the
pair and each region retains its original five physical volume-loss channels.
`ThermoviscousTimeWork::conductive_boundary_work_j` is acoustic available-energy
transfer, not caloric heat or dissipation. An actual signed `integral(q dA dt)`
is exposed separately only for a complete active thermal trace and never added
to the quadratic acoustic ledger.

Stationary constraints present on one side also constrain the other side.
Coordinates fixed on both sides are removed from the active trace test space;
reported tractions/fluxes there are only Galerkin representatives modulo
work-null loads, not claims of zero physical stress or heat flux. When the
thermal trace is incomplete, `caloric_heat_transfer_j` is `None`: the actual
caloric integral cannot be inferred from the available-energy transfer.
Moving structural traces overlapping interface DOFs, mismatched geometry,
rank loss, mixed internal/external pressure modes, inconsistent initial traces,
different clocks/equilibria and exceeded dense budgets refuse without publishing
either state. The original normal-flow `connect` retains its insulated thermal
and free-tangential interpretation through the same generalized pair owner.

Scope remains fixed geometry and one common reference pressure/temperature,
bounded dense exact held-input time evolution, and source-owned fixed material
coefficients. There is no mean-flow/species transport, mean thermal-reservoir
update, moving mesh or unqualified acoustic-resolution claim.

Focused regressions compare the projected dynamics with independent primal
identification of the original physical trace coefficients and original weak
mass/force equations, check the actual thermal-flux units, permute geometric
node numbering, retain all molecular fields, test shear and thermal transfer,
exercise one-sided/shared stationary constraints, and check refusal/retry.

## Spatial conduction in original finite solid-wall fields

`ThermoviscousFe::prepare_time_with_spatial_wall_fields` extends the original
finite-wall time owner with explicitly declared `TangentialWallConductivity`.
Each homogeneous patch is a planar extrusion of its actual boundary triangles:
construction checks coplanarity and the common outward normal. The constitutive
declaration is transverse isotropy, with independent positive normal and
in-plane conductivities and zero normal/tangent cross term. The existing normal
scalar is never assumed to supply the tangential coefficient.

The original hierarchical P2 surface gradient owner and P1 depth mass assemble
`T0*k_t*K_surface tensor M_depth`; normal conduction remains
`T0*k_n*M_surface tensor K_depth`. All original surface/depth coefficients,
consistent front cross mass, molecular fields and structural coordinates are
retained. This changes no storage coordinate or mass matrix. Both conduction
forms enter the same original resistance and trajectory-integrated solid heat
channel, with the same held-input step and energy admission.

Outer slab sides are naturally insulated. Separate regions have no inferred
tangential bond. Backing and exact front-temperature constraints retain their
original meanings. Actual finite thickness is required; curvature, multilayers,
general anisotropic cross terms and evolving mean/material state are outside
this path. The normal-only path retains its original asymptotic admission;
spatial conduction uses the planar geometry admission instead. Complete mesh,
clock, memory/work and cancellation bounds still apply. A stated source band
alone establishes no spatial or waveform accuracy.

## Shared planar coefficient admission for matrix-free fields

`weighted2::CellWeight2::validate` exposes the existing scalar positivity and
outward-interval tensor SPD admission to original matrix-free consumers. It
rejects unsupported form degrees and preserves the same coefficient checks as
assembled Hodge operators. Admission proves neither mesh validity nor a global
boundary/gauge condition; those remain with each original field owner.

## Resolved passive equivalent-fluid volumes

`acoustics::helmholtz::heterogeneous::PreparedEquivalentFluid` extends the
original full pressure FE space to complex effective density and bulk modulus
at one explicit frequency. In `exp(+i omega t)` it retains
`div((1/rho_eff) grad p) + omega^2 p/K_eff = 0` on original affine or quadratic
tetrahedra. Material callbacks own their constitutive law, source domain and
spatial variation; jumps must conform to original element faces. The original
boundary forms, exact pressure-release constraints, equilibrated complex LU
and physical residual corrections are reused. No modal reduction is introduced.

Positive-real passive rho/K admission fixes the loss signs. Each original
element retains separate positive density and compression loss forms, evaluated
on the coherent pressure field as `Im(1/rho_eff)|grad p|^2/(2 omega)` and
`-omega Im(1/K_eff)|p|^2/2`. Their quadrature volumes and powers remain in
original tetrahedron order. Boundary absorption stays separate; no net-power
subtraction supplies heat. Explicit total volume-sample, retained local-matrix,
free-DOF and dense-factor limits precede material callbacks. Cancellation
publishes no partial preparation, pressure field or heat table.

`tests/equivalent_fluid.rs` specifies three-rung lossy plane-wave convergence,
independent density/compression heat integrals, complex-interface reflection,
the original lossless pressure-release limit, material signs, limits and
cancellation. A single-frequency constitutive sample supplies neither a
dispersive stored-energy law nor a causal transient realization. Mean flow,
elastic porous-frame waves, negative-parameter media and experimental material
qualification remain outside this numerical formulation.

### Original distributed sources in lossy fluid volumes

`PreparedEquivalentFluid` prepares the existing `VolumeVelocitySource` and
`BoundaryVelocitySource` owners through the same original volume/face projection
loops as the real heterogeneous solver. Volume shapes retain physical positions,
original tetrahedron indices and full pressure-basis loads; quadratic boundary
motion retains all six nodal values and the original curved metric. Signed
shapes are not normalized or reduced to their integral. Each source borrows its
exact prepared frequency/material/wall owner; shared geometry alone does not
permit cross-wiring, even at zero amplitude.

`solve_with_sources` coherently combines raw inward velocity, prepared quadratic
boundary motion, distributed volume injection and original wall tractions in the
retained full pressure/motion/memory factor. Boundary amplitudes are peak m/s,
volume amplitudes are peak m^3/s and traction amplitudes are peak Pa, all in
`exp(+i omega t)`. Independently evaluated boundary/volume source work remains
signed and is already included once in total fluid input work. It is never
material heat. Native bulk, Robin and mechanical dissipation remain separate.
Empty prepared-source slices preserve the original raw-drive arithmetic; warm
solves invoke no source or material callbacks. Original support, source-count,
quadrature, finite-value and cancellation admission remains in force.

`tests/equivalent_fluid_sources.rs` specifies uniform complex compliance and
signed coherent work, zero-net spatial injection, midside boundary excitation
with all wall memories, curved volume/face moments, the original lossless source
limit, strict owner binding and cancellation without partial publication.

### Reciprocal mechanical walls in complex equivalent-fluid volumes

`PreparedEquivalentFluid::new_with_coupled_mechanical_walls` retains original
wall displacement, normal velocity and every Maxwell-memory coordinate in the
same complete complex-material pressure factor. Real and complex bulk fluids
share the original wall trace/material preparation and `CoupledPencil`. The same
consistent face matrix transfers pressure force and wall volume velocity in
opposite work accounts. No wall-impedance inverse, modal reduction or artificial
damping is introduced at a local wall-impedance zero.

Original Robin, mechanical and pressure-release faces must be disjoint;
pressure release on adjacent faces removes no mechanical coordinates. The
equivalent-fluid harmonic DOF and dense-byte limits count all active pressure,
wall-motion and memory coordinates before any volume-material callback. Cached
traction profiles remain bound to their exact preparation and act as physical
pressure forces, independently of prescribed fluid velocity sources.

The response exposes separate signed fluid source work, fluid-to-wall transfer,
external traction work, boundary heat, density/compression volume heat and wall
resistance/relaxation heat. Fluid, wall and complete-system balances retain these
independently evaluated terms. Wall storage follows the original mechanical law;
no dispersive bulk stored energy or frequency derivative is inferred from a
single-frequency rho/K sample. Original rigid porous frames remain rigid; the
added boundary motion is not an elastic Biot interior or a moving pore skeleton.

Both prepared harmonic owners expose `wall_dissipation` in original patch order.
`MechanicalWallPhasors` retains each solved energy-normalized extension
`e_j = sqrt(k_j)*w - z_j` directly, as does its analytic frequency tangent.
Independent wall heat reuses the original positive forms `R*v^H*W*v/2` and
`sum_j e_j^H*W*e_j/(2*tau_j)` used by the native aggregate accounting; it never
reconstructs a tiny extension by subtracting rounded displacement and memory.
Complete state/branch dimensions and finite values are checked, and cancellation
returns no partial heat table. The caller owns coherent state/frequency binding.

Each `MechanicalWallPower` retains the original complete patch area and
`surface_density*area` mass. Area is the same physical trace quadratic form on
the hierarchical constant: one on vertex DOFs, zero on all edge/face bubbles.
This includes original curved-face weights and avoids the incorrect sum over
every hierarchical matrix entry. Pressure exchange, applied traction, Robin
loss, bulk heat and aggregate work residuals do not supply wall heat. Thermal
receivers still require an explicit physical mapping and original caloric law.

`tests/equivalent_fluid_mechanical.rs` specifies G1 lossy-layer convergence to
independent continuum pressure, wall motion/memory and heat integrals; G3
coherent fluid/traction work and lossless agreement with the original wall owner;
and G0/G4 local zero-impedance, adjacent pressure-release, source binding,
complete-state admission and cancellation behavior. Its original-patch heat
case specifies high-order area/mass normalization and local coherent scaling;
native `frequency/coupled/tests.rs` specifies heat that remains positive when
rounded displacement/memory subtraction erases the physical extension.

## Heterogeneous resolved pressure acoustics

`acoustics::helmholtz::heterogeneous::HeterogeneousHelmholtzFe` borrows the
original hierarchical H1 tetrahedral geometry/trace owner and leaves its
homogeneous medium and geometric stiffness/mass unchanged. It evaluates one
coherent density/sound-speed pair per physical volume quadrature point and
assembles both physical operators together:
`K_rho = integral((1/rho) grad(phi_i).grad(phi_j))` and
`M_beta = integral(phi_i phi_j/(rho c^2))`. Density stays inside the pressure
divergence, preserving the weak pressure/normal-velocity interface conditions.
The lower FE layer imports no material model: the caller retains original
source identity, interpolates primitive fields through that source, and
supplies checked coherent samples. Original source errors propagate unchanged.
True coefficient jumps require conforming tetrahedron faces; continuous fields
retain separate geometry, h/p and quadrature error axes.

Rigid modes use the original sparse modal owner with explicitly named
`omega_squared_rad2_s2` eigenvalues and the physical compressibility mass.
The prepared harmonic owner factors the FULL original pressure space with
`fs_la::eigen_complex::lu_complex_with_cancel`. Physical outward specific
admittance and spatial inward velocity enter
`(K_rho - omega^2 M_beta + i omega C_Y) p = i omega b_v` in `exp(+i omega t)`;
there is no reference-medium or port-area factor. Velocity traces retain their
caller-specified face/vertex order. Passive admittance groups require
nonnegative real part and disjoint boundary faces. Sources may share an
admittance face, with `v_out = Y p - v_imposed`.

One factor serves subsequent drives without new source queries. Symmetric
diagonal equilibration and at most two corrections through that same factor
use residuals of the original separate sparse physical terms in double-double
arithmetic. Preparation has explicit DOF/dense-storage limits (hard cap 1536
pressure DOFs; defaults 1024 DOFs and 64 MiB for the two scalar dense matrices).
The reported backward error is algebraic and is not a field-accuracy bound.
Input work, passive boundary dissipation, kinetic energy and compression
energy retain their physical forms; the signed power-balance residual is
reported rather than set to zero. Only negative quadratic-form residues
bounded by assembly roundoff are zeroed. Evanescent propagation introduces
no bulk heat.

### Analytic frequency response of the full pressure space

`PreparedHeterogeneous::{solve_with_frequency_derivative,
solve_with_sources_frequency_derivative}` returns the original pressure and
its complex derivative per **rad/s**, with the same equilibrated LU. The
original sparse pressure equation is differentiated analytically:
`A p' = i b + i omega b' + 2 omega M_beta p + Im(C_Y + omega C_Y') p - i Re(C_Y + omega C_Y') p`.
Linear boundary traces, prepared quadratic boundary motion and volume sources
retain their coherent amplitudes and independently supplied amplitude/phase
slopes. A zero primary amplitude does not erase a nonzero source derivative.
Exact pressure-release coordinates remain zero in both vectors. Derivative
coefficients are pressure slopes, not a second physical field or energy state.

Geometry, density and compressibility are held fixed. Existing preparations
hold their complete complex boundary admittances fixed (`C_Y' = 0`). Explicit
same-frequency samples supplied through `prepare_with_dispersive_admittances`
retain their signed analytic boundary slopes; an ordinary sampled uniform or
spatial law never implies an unspecified dispersion derivative. The tangent
introduces no finite difference, modal reduction, artificial loss, fitted pole
or extra factorization. Additional solve storage is linear in the original DOF
count; explicit boundary derivatives retain sparse face matrices. Each solve
retains at most two residual corrections.

The sparse solve/tangent leaf retains double-double physical forcing and checks
the separate original operator terms. The tangent residual denominator includes
the separate source and operator-derivative terms, preserving coherent
cancellation. Original source admission, residual tolerance and cancellation
apply; a refusal leaves the retained factor reusable and publishes no partial
response. Algebraic residuals do not establish spatial, geometric or resonance
accuracy, a modal census, an isolated linewidth, or a decay quality factor.

Seven native G1/G3 tests exercise the complete production numerical leaf using
the original `fs-sparse`, `fs-math` and complex-LU implementations: a high-Q
resonator at its peak, sealed low-frequency compliance, convergent pressure-soft
duct standing waves, coherent delayed sources with complex boundary terms,
zero-amplitude source slopes, physical-equation rescaling, and cancellation
with factor reuse. The full FE crate and its thin geometry/source adapter were
not built in that validation environment because the constellation runtime's
registry dependencies were unavailable; these native tests do not claim that
broader integration gate or instrument qualification.

`HeterogeneousHelmholtzFe::survey_pressure_resonances` binds the existing
shared driven-resonance survey to an original tetrahedral pressure observation
under constant coherent velocity traces and fixed uniform admittances. Each
frequency uses the original full-pressure solve and analytic tangent. The
result identifies grid-resolved driven peaks; it supplies neither a free-mode
census nor linewidth/decay-Q evidence. Two original-FE integration regressions
specify a loaded duct's spatial resonance and silent constrained observations,
input refusals and cancellation. These full-adapter tests were authored but not
executed without the complete constellation dependency environment.

`HeterogeneousHelmholtzFe::probe` prepares a one-sided observation using the
original hierarchical basis and physical barycentric gradients retained from
assembly. It queries the caller's original frozen material once at the exact
element/position, applies the same coherent density/sound-speed checks as
assembly, and retains that local `Medium`; it never uses the geometry owner's
uniform medium or averages across an interface. Caller-owned source identity
remains the caller's responsibility. Subsequent fields reuse the local stencil
without allocations, material queries or numerical differencing. Harmonic
samples recover pressure, its gradient, `v = i grad(p)/(omega rho)`, and active
and reactive parts of `p conj(v)/2`. Coherent amplitudes are combined before
forming intensity. Transient samples require a state from the exact original
operator and recover `v = -grad(psi)/rho` with instantaneous flux `p v`.
Reactive evanescent flow is retained and is not interpreted as heat or normal
transported power. Pointwise gradient/flux error has its own h/p convergence
obligation; probes do not impose continuity or correct the solved field.
Spatial interpolation and derivative products now retain double-double
residuals through pressure, velocity and their quadratic intensity. The same
rule applies to harmonic and transient surface-power quadrature. This preserves
weak signals at destructive nodes and finite transported flux between large
nearly cancelling components; recomputing intensity from the separately rounded
display phasors need not reproduce that retained flux. No allocation or extra
material query is introduced. The original binary64 coefficients, geometry and
material accuracy still bound the observation.
Focused checks cover analytic P1–P4 gradients on skew/reordered tetrahedra,
one-sided density-jump velocity, propagating/standing/evanescent flux, a solved
Snell intensity mesh ladder, and transient gauge-invariant velocity/flux.

`HeterogeneousHelmholtzFe::surface_probe` integrates signed outward power on
original element faces. The selected tetrahedron fixes both material side and
normal; duplicate physical faces refuse. Harmonic results retain active and
reactive `integral(p conj(v_out)/2)` separately; transient results use
`integral(p v_out)` without the half factor. Prepared normal stencils require
no further source queries and accumulate with fixed-order double-double sums.
Caller-visible face/point budgets and quadrature order bound preparation.
High-order field gradients and continuous-material surface quadrature require
their own convergence; coarse interface-flux discrepancies are not corrected
into continuity. Reactive power is not dissipated heat.

Focused regressions cover homogeneous mode recovery, a matched traveling
wave and its energy/work, same-speed density-interface reflection, analytic
continuous sound-speed travel time, subcritical and evanescent Snell fields
on conforming three-dimensional meshes, an independent continuous-coefficient
quadrature axis, reciprocal spatial traces, and source/limit/cancellation
refusals. These fixtures use analytic coefficients rather than claiming a
new material dataset. The model remains stationary, linear, isotropic,
lossless in the bulk and without mean flow or thermal/species evolution.
Homogeneous thermoviscous and molecular wall/volume models are not silently
reused for these spatially varying fluids. Assembly and harmonic operations
poll cancellation; modal cancellation brackets the existing single-window
modal owner, which has no internal callback.

### Distributed sources in heterogeneous fluids

`HeterogeneousHelmholtzFe::volume_velocity_source` projects a prescribed signed
spatial density [m^-3] against every original hierarchical pressure test.
It retains the original operator identity and explicit source quadrature; a
volume-flow amplitude [m^3/s] gives the continuity source `s` [1/s]. The shape
is neither normalized nor replaced by a singular point. Source sampling occurs
once, independently of material assembly, and original callback errors survive.

`PreparedHeterogeneous::solve_with_volume_sources` and the sparse time owner's
`step_with_volume_sources` compose these sources with the original boundary
velocities before solving. They reuse the same factors, physical residuals,
pressure observations and work checks. For `exp(+i omega t)` the pressure RHS
is `i omega integral(phi s)`; the potential-time RHS is `integral(phi s)`.
No extra density, sound speed, volume or face-area factor enters either path.
The signed volume-source power/work subtotal uses the same solved pressure and
is already included in total input work. It is external work, not bulk heat.
Coherent sources interfere before field and intensity evaluation; reversing
injection can extract compression energy.

Focused G1/G3 regressions cover an independently solved distributed source in
an anechoic duct with mesh refinement and power balance, sealed compression
across a density interface, simultaneous boundary forcing, coherent cancellation,
weak and reversed drives, source reuse, and transactional refusals. These are
prescribed linear acoustic sources, not a resolved bubble, thermal expansion,
mass-transfer or combustion model. Callers must derive source amplitudes from
their physical owners and resolve source shape, h/p and observation bandwidth.
At most 128 simultaneous volume-source columns are admitted per query; the
bounded exact owner's original augmented-state limits apply additionally.

### Broadband heterogeneous pressure waves

`heterogeneous::time` advances the same full spatial operator in pressure
potential `psi` [Pa s] and actual pressure `p = psi_dot` [Pa]. With the frozen
original density, velocity is `-grad(psi)/rho`, and the physical equations are
`M_beta p_dot + C_Y p + K_rho psi = b_in`. Storage is therefore
`(p^T M_beta p + psi^T K_rho psi)/2` joules. A pressure/pressure-rate mechanical
energy surrogate is not used. Every original hierarchical coordinate remains
live, with no modal truncation, dense mass inverse or subsequent material query.

Fixed-step implicit midpoint prepares one symmetrically equilibrated sparse
`M_beta + h C_Y/2 + h^2 K_rho/4` factor through the original LDLT owner.
Admission defaults to 4096 pressure DOFs and 64 MiB of native factor/front
payload; sparse assembly, symbolic analysis and runtime overhead are additional.
Real nonnegative boundary admittance is instantaneous and memoryless. Reactive
or frequency-dependent impedance requires a causal state and is not admitted.
Inward velocity traces use the original face projection at the step midpoint.

Returned candidates preserve the caller's state on refusal. Independent
double-double residuals of the original physical equations, at most two
corrections through the same factor, and a signed physical work balance admit
each step. Input work and boundary dissipation use the same midpoint pressure
with their original forms. Energy change is evaluated as a symmetric bilinear
difference; its tolerance scales with actual work and evaluated roundoff,
without a fixed SI floor or an invented heat channel. Cancellation polls bracket
the existing trace/factor/solve owners and occur between sparse rows; those
owners have no internal callback.

Focused regressions specify exact polynomial compression and pressure work at
two amplitude scales, rigid/zero-state behavior, actual Cayley phase and its
second-order time convergence, a broadband density-interface pulse with
independent reflection/transmission amplitudes and delay, and resource/source
binding refusals. These new tests have not been executed in the current
environment; pinned rustfmt parsing and source/API review do not replace native
validation. Midpoint stability does not certify phase accuracy: time-step,
geometry, h/p and material quadrature refinement remain separate obligations.

### Exact complete-field heterogeneous acoustic evolution

`heterogeneous::time::exact` prepares the same full `M_beta`, `K_rho` and
real boundary-admittance forms for exponential held-input evolution. It
retains every original pressure/potential coordinate, including the constant
potential gauge. Canonical momentum is `pi = M_beta p`; diagonal coordinate
scales condition the original mass inverse, storage and interconnection.
The inverse is checked against every original mass entry. No eigenmode is
truncated, frequency matched, retuned or pinned to remove the gauge.

`prepare_exact_time_with_volume_sources` admits the same original-volume
source loads alongside the boundary ports. `step_with_volume_sources` takes
separate held boundary velocities [m/s] and volume flows [m^3/s], retaining
their distinct integrated pressure observations [Pa m^2 s] and [Pa s]. Both
source families share the original complete trajectory and independent loss
Gramian; their signed work appears exactly once in total input work. No new
dynamic states or source callbacks enter a warm step. Analytic compression at
air/water SI scales and a joint boundary-drive/volume-extraction polynomial
check the complete pressure field and the integrated signed source work.

`ExactVelocityPort` binds a dimensionless piecewise-linear inward-velocity
shape to the original boundary quadrature. Step amplitudes have units m/s.
Different ports may overlap and interfere; a source may share a real-admittance
face. Port work is the held amplitude times the actual integrated
`integral(shape p dA)` pressure, with boundary area appearing exactly once.
Absorption independently integrates `p^T C_Y p` along the same cached
trajectory through the existing pHS Gramian owner. It is neither an endpoint
quadrature nor an energy remainder. Returned pressure/potential states retain
the original material/geometry identity and work with existing field probes.

Preparation is a bounded dense small-domain path, using the existing full-state
exponential owner. The current cap is
`2*(2*pressure_dofs + boundary_ports + volume_sources + 1) <= 256`,
in addition to the caller's original-DOF, dense-cache byte and multiply-work
limits. With zero or one port this permits at most 63 pressure DOFs. Larger
fields refuse without reducing the basis or changing propagators. Warm steps
perform cached contractions and allocate result/scratch vectors; they perform
no material query, factorization or matrix exponential. Cancellation is polled
between preparation rows and around bounded original numerical operations,
and before any candidate is returned. All refusals preserve the accepted state.

Exact means the complete semidiscrete held-input solution up to numerical
exponential error, without midpoint's temporal frequency transformation.
Spatial h/p, geometry and material quadrature still require convergence;
source holding, recording bandwidth and antialiasing remain separate axes.
This path covers frozen isotropic fluids with no mean flow or bulk loss and
real memoryless boundary admittance. It does not promote the sparse midpoint
route, an under-resolved mesh, or a small-domain fixture to instrument acceptance.

### Causal distributed mechanical boundaries

`heterogeneous::time::mechanical` composes the original `fs_phs::WallPin`
surface mass, stiffness and resistance with that same full pressure solver.
Every active original H1 boundary coordinate retains displacement `w` and
outward normal velocity `v`; no uniform piston or diagonal surface mass replaces
the pressure trace. The original symmetric face quadrature supplies
`W = integral(phi_i phi_j)`. Its area enters exactly once in fluid feedback
`+W v`, wall storage `(sigma v^T W v + K w^T W w)/2`, resistor loss
`r v^T W v`, and pressure work `p^T W v`. Wall coefficients are specific to
unit area, matching the original `WallPin` admission and units.

For constant coefficients on a patch, midpoint eliminates the original law
`sigma v_dot + r v + K w = p` to a positive effective admittance
`1/(2 sigma/h + r + h K/2)` and an explicit history. The existing sparse
factor and residual-correction owner solve the condensed acoustic equation;
all wall histories then reconstruct from that same solution. The effective
admittance is an algebraic coefficient, never a physical heat account.
Independent checks retain the uncondensed fluid and wall equations, acoustic
and wall storage changes, actual resistor loss, and signed fluid-to-wall work.
Lossless walls can return stored energy. Imported and returned wall states
remain bound to their exact prepared configuration, and refusal never advances
an accepted participant. Admission adds a default 4096 summed wall-trace DOF
cap, each with two physical coordinates, to the original acoustic limits.

Mechanical and real Robin faces must be pairwise disjoint; imposed velocity
may share either. Adjacent mechanical patches have independent local states
at shared edges. This is a frozen, linear, locally reacting wall with fixed
fluid geometry, without lateral shell/bending waves or inferred material data.
`MechanicalBoundary::relaxation` admits the same positive surface-specific
`AcousticWallRelaxation` spectrum used by the original acoustic chain owner.
Each original trace coordinate retains every energy-normalized Maxwell memory
`z_j`, with `e_j = sqrt(k_j) w - z_j` and `tau_j z_dot_j = e_j`.
These branches add physical force `sqrt(k_j) e_j`, storage `e_j^T W e_j/2`
and dissipation `e_j^T W e_j/tau_j`. Original consistent face mass is used
for all three quantities, including off-diagonal and high-order trace terms.
The branch equilibrium stiffness is not counted twice in `WallPin`.

With `a_j = tau_j/(tau_j+h/2)`, branch condensation adds
`h sum(a_j k_j)/2` to the same positive wall impedance and retains the old
memory in the projected history. The original `RelaxationBranch` scalar
midpoint owner reconstructs each memory after the joint pressure solve; no
extra global factor, split relaxation stage, fit, or new matrix coordinate is
introduced. Preparation and stepping remain linear in the number of retained
wall memories beyond the existing sparse pressure solve. A separate default
16384-coordinate memory limit admits the sum of trace DOFs times arm counts.
Zero-arm walls retain their original arithmetic and physical storage.

`state_with_wall_memory` imports the complete original arm histories.
`state_with_walls` refuses a hereditary wall rather than silently inventing
its past. `zero_state` and `state` explicitly start all walls and memories at
rest. `wall_relaxation_j` is separate stored energy; the Maxwell dissipation
subtotal is already included in total wall dissipation and must not be added
again. Original uncondensed fluid, mechanical and memory residuals and the
three signed energy balances admit the complete candidate transactionally.
Material source, geometry, valid frequency band and time-resolution limits
remain with the caller or material adapter. Midpoint's passivity does not
establish accuracy for time steps larger than a material relaxation time.

Focused hereditary regressions add continuous constant-velocity solutions
with independent storage/heat and area/amplitude scaling, nonuniform P2
prestress release from zero wall displacement, complete-state/budget/refusal
checks, and complex reflected phase/amplitude against the independent Maxwell
impedance over three time-step rungs. These are numerical-model checks, not
experimental qualification of an arbitrary wall material.

Five focused tests specify exact polynomial wall/fluid motion and pressure work
under independent area and amplitude scaling, a nonuniform P2 wall's analytic
fourth-order energy moment and release, reflection phase across the original
wall resonance with three temporal rungs, exact original evolution with no
mechanical walls, and resource/overlap/source-binding refusals. Native execution
remains pending in the current environment; source
review and pinned rustfmt parsing establish neither numerical test passage nor
full material or shell-model validity.

## Exact pressure-release acoustic modes (M2.4-W1a / pexfp.93)

`HelmholtzFe::modes_with_pressure_release` and
`HeterogeneousHelmholtzFe::modes_with_pressure_release` impose exact homogeneous
Dirichlet pressure on explicitly selected original boundary triangles. Every
nonzero hierarchical face trace is eliminated, including vertex, edge and face
coordinates through P4. All remaining faces are rigid. The common constraint
map restricts the original sparse stiffness and consistent mass; it introduces
no penalty, large admittance, resampled medium or spurious boundary eigenvalue.
The heterogeneous pencil retains `1/rho` inside stiffness and `1/(rho c^2)` in
mass, including conforming material interfaces.

The original modal solver supplies free-equation residuals and inertia counts.
Returned shapes use full original coordinates with exact zeros on constrained
entries; constrained equation rows represent boundary reactions. Distinct
`PressureReleaseModes` and `HeterogeneousPressureReleaseModes` wrappers retain
canonical faces and free-DOF counts. They cannot be passed directly to the
existing rigid-wall reduced constructor, whose static vectors are unconstrained.
Empty face sets preserve the original rigid solver's eigenpairs. Repeated,
internal or missing faces, an empty free space, and cancellation refuse without
mutating the original operator. Cancellation brackets the native eigensolver;
its existing single-call cancellation granularity is unchanged.

Focused tests specify quarter/half-wave duct spectra over three P2 mesh rungs,
mixed 3-D box modes, a closed-form density-interface resonance, original-mass
normalization, P1--P4 trace elimination, winding invariance and refusal paths.
This boundary is an explicit ideal `p=0` condition, not an opening-radiation or
end-correction model. Curved/isoparametric geometry and the remaining M2.4-W1a
qualification requirements remain outside this increment.

`HeterogeneousHelmholtzFe::prepare_with_pressure_release` applies the same
original trace constraints to driven harmonic fields. It factors only the free
pressure coordinates and lifts exact zeros into the complete original field.
The retained density and compressibility operators, passive complex admittances,
spatial velocity traces and distributed volume sources remain physical; neither
material interfaces nor sources are resampled. Dense factor limits count the
free coordinates, while original sparse operators and lifted vectors retain
their separately documented storage. Empty constraints preserve the original
preparation and solution arithmetic.

An ideal pressure-release face cannot simultaneously own an admittance or an
imposed velocity. Adjacent faces may share constrained edge or vertex traces.
Residual correction checks only the free equations; eliminated rows carry
boundary reactions. Stored energy and signed input work use the full original
operators, and the zero-pressure reactions do no work. Cancellation leaves a
prepared factor reusable and does not publish a partial field.

Driven regressions specify mesh-converged standing-wave pressure and both stored
energy components, a conforming density/sound-speed interface with independent
transfer-matrix pressure and passive boundary heat, P4 traces, original volume
sources, dense-memory bounds and conflicting-law refusals. An ideal zero-pressure
termination is still not a finite-opening radiation or end-correction model.

`HeterogeneousHelmholtzFe::prepare_time_with_pressure_release` extends the same
trace elimination to sparse midpoint evolution. The original positive-definite
effective operator is factored only on free coordinates, with its native sparse
memory limits; pressure, potential, sources and physical work retain their full
original coordinates. Free-equation residual corrections lift exact zeros into
the constrained pressure trace. No dense mass inverse or extra modal basis is
introduced, and the existing mechanical-wall path keeps its empty constraint
map and original ordering.

This time owner explicitly admits zero pressure AND zero pressure-potential
trace. Since `p = psi_dot`, zero pressure alone would permit a time-independent
nonzero potential trace; arbitrary stationary flow components are outside this
entry point. Imported states must choose a compatible potential gauge and are
refused rather than projected when a constrained coefficient is nonzero. State
identity, passive boundary work, original volume-source work and transactional
cancellation checks remain in force. Empty release sets preserve the original
warm-step path without extra source-face sorting or repeated finite-state scans.

Transient regressions separate spatial mode accuracy from midpoint temporal
phase refinement, verify conserved acoustic energy, and compare a forced P2
polynomial pressure/potential solution with independent kinetic, compression,
volume-input and passive-wall work. They also exercise strict high-order state
admission, conflicting laws, free-DOF limits and reusable factors after refusal.

`prepare_time_with_walls_and_pressure_release` composes these same constraints
with the original causal mechanical and Maxwell-memory wall owner. Mechanical,
real-admittance and pressure-release faces are disjoint, while adjacent patches
may share original trace coordinates. Every wall displacement, velocity and
memory coordinate remains present even on an edge shared with a constrained
fluid trace. The original consistent face mass, wall storage, actual resistor
and relaxation heat, and signed fluid-to-wall work remain unchanged.

The condensed acoustic solve uses the free fluid coordinates and its existing
factor budget; wall and memory limits still admit their complete trace spaces.
The mechanical step explicitly refuses velocity forcing on release faces before
calling the common projected solve. Empty release maps use the original wall
path. Tests specify an exact positive-stiffness spring-balanced wall with an
accelerating fluid at independent area/amplitude scales, a nonuniform Maxwell
memory pulse adjoining a release edge with irreversible heat and acoustic
energy transfer, and full-state/overlap/resource/cancellation refusals.

### Distributed volume excitation of causal mechanical walls

`PreparedMechanicalTime::step_with_volume_sources` composes the existing
original-volume source loads with boundary velocities and all distributed wall
and Maxwell-memory coordinates in the same sparse midpoint solve. Held spatial
source shapes are never collapsed to their net injection, reprojected or
requeried during a warm step. Amplitudes are midpoint volume flows [m^3/s];
negative values are admitted. The same full pressure field supplies source
work, the uncondensed fluid residual and the reciprocal fluid-to-wall exchange.
`MechanicalTimeWork::volume_source_work_j` is a signed subtotal already included
once in `input_work_j`, not another energy term. Existing source count, model
identity, finite-value and transactional cancellation checks apply. Empty source
sets retain the original wall step and arithmetic; no new factor or dynamic
state is introduced. Pressure-release constraints also remain active.

Regressions specify an exact continuum polynomial with simultaneous boundary
forcing, volume injection/extraction and steady Maxwell wall motion, with
independent area/amplitude scaling and analytical storage/heat/work. A signed
source with zero net volume excites a resolved pressure field and a wall beside
a release face; source cancellation and rejected cross-model/nonfinite/cancelled
inputs exercise reuse of accepted states. These native regressions still need
execution on a provisioned worker; parser checks and analytical work identities
do not replace native numerical validation.

### Exact trajectories with pressure-release boundaries

`HeterogeneousHelmholtzFe::prepare_exact_time_with_pressure_release` applies the
shared original hierarchical trace constraints to exact held-input evolution.
Canonical momentum is `pi_f = M_ff p_f`: the consistent mass, stiffness, real
boundary conductance and source columns are restricted BEFORE mass inversion.
Every free coordinate remains, with the existing quadratic-flow, source
observations and loss Gramian. Complete original pressure/potential fields are
recovered with exact zeros on the release trace. Physical energy and signed
boundary/volume-source work use the original forms.

Both pressure and pressure-potential traces must be zero on imported states,
including states from another preparation on the same model. Conflicting
admittance or velocity-port faces refuse; adjacent faces may share constrained
coordinates. Dense pressure and augmented-state limits count free coordinates.
Empty release sets preserve the original arithmetic and reject oversized full
fields before allocating the identity map. Cancellation publishes no candidate.
This remains bounded, frozen-medium propagation with lossless bulk and real
memoryless boundary admittance. Ideal pressure release supplies no exterior
radiation or opening correction.

Four native regressions specify heterogeneous P2 mode phase without midpoint
frequency warping; independent analytical tetrahedron dynamics with boundary
injection, volume extraction and Robin heat; a P3 field with 64 original/eight
free coordinates plus strict trace/resource/cancellation refusals; and bitwise
empty-map compatibility. Pinned Rust parsing and the independent analytical
heat/work oracle were checked; native execution remains pending because the
current environment kills the upstream dependency build before fs-feec.

### Quadratic boundary motion across the original acoustic solvers

`HeterogeneousHelmholtzFe::quadratic_velocity_source` prepares a signed six-node
inward-normal velocity shape on original straight boundary triangles. Values
are dimensionless samples at corners 0/1/2 and edge midpoints 01/12/20, in the
supplied vertex order. Original H1 trace quadrature retains midside motion and
higher-order pressure tests, including shapes with zero net volume flow. No
area normalization, face-average replacement or corner reduction is applied.
Each prepared `BoundaryVelocitySource` retains the original model, canonical
support and load [m^2], with a signed area integral available for observation.

Harmonic `solve_with_boundary_sources`, sparse-time and mechanical-wall
`step_with_boundary_sources` combine original linear traces, these prepared
boundary columns and original volume columns before the same physical solve.
Prepared boundary amplitudes are velocities [m/s]; volume amplitudes remain
volume flows [m^3/s]. Their signed conjugate-pressure work subtotals are already
included once in total input work. Every distributed wall and Maxwell memory
coordinate remains present. `MechanicalTimeState::acoustic_state` borrows the
original typed field for the existing point and surface probes without a copy
or loss of spatial/material identity.

`prepare_exact_time_with_boundary_sources` appends the same cached columns
after linear velocity ports. Combined boundary `port_count`, held velocities
[m/s] and integrated pressure observations [Pa m^2 s] follow that order; volume
flows and observations [Pa s] stay separate. Every column enters the existing
augmented-state budget, forcing matrix and same-trajectory work observation
exactly once. The original exponential, independent loss Gramian and restricted
consistent mass are unchanged. No extra dynamic state or propagator is added.

Preparation projects each shape once. Repeated drives perform no source/material
query or boundary quadrature. At most 128 prepared boundary sources are admitted;
exact evolution also applies its stricter augmented-state limits. Model identity,
finite data, physical support, pressure-release conflicts and cancellation are
checked, including for zero amplitudes/shapes. Empty prepared families retain
existing paths, and empty release sets avoid repeated support-face traversal.
This extends prescribed excitation on straight geometry; structural laws,
curved geometry and radiation remain separate physical owners.

Ten native regressions specify independent triangle moments and all vertex
permutations, P3/P4 pressure-face coupling, zero-integral motion, coherent signed
harmonic work, simultaneous volume input, midside-only transient/wall excitation,
full probe access and exact tetrahedron dynamics with integrated Robin heat.
Parser/format checks and independent exact-rational moment/work oracles passed.
Native execution remains pending under the upstream dependency-memory blocker;
these source checks do not establish numerical test passage or an instrument
acceptance claim.

### Complete pressure-coordinate multiports and exterior loading

`InteriorMultiport::new_full` and `new_full_with_boundary_layers` retain every
original `HelmholtzFe` pressure coordinate. Original CSR stiffness/mass entries
are copied directly into the existing dense loaded-factor storage, and signed
port columns keep their native geometric-area normalization. No eigenmodes,
identity-basis projection, static enrichment or discarded pressure DOFs are
introduced. A caller-specified DOF cap and hard ceiling of 1024 refuse oversized
fields before dense allocation. Unused legacy per-port dense face matrices are
not allocated. The existing tangential viscous/thermal wall forms can likewise
retain their complete original coordinates.

The native loaded scattering and field recovery serve these full fields as
well as the existing caller-supplied Galerkin path. Responses from `new_full`
contain all original FE coefficients in original order; positive real power
references are numerical coordinates. The source-bound orchestrator in
`fs-couple::cavity_modal::interior::exterior` uses this complete constructor to
couple finite openings to the original BEM self/mutual impedance, with physical
pressure/flow work and original port residual checks. Dependency direction is
unchanged: fs-feec does not depend on fs-couple or fs-bem.

Four focused coupling regressions specify an independently assembled original
CSR FE/BEM mixed block, complete P2 fields, signed zero-net opening modes,
reference invariance, finite reactive loading, native wall heat and binding
refusals. Pinned Rust parsing and independent interface algebra were checked;
native numerical execution remains pending under the environment's upstream
dependency-memory blocker. Pressure, trace, geometry and exterior discretization
still require convergence, and no causal exterior realization is implied.

## Original quadratic geometry throughout pressure acoustics (M2.4-W1a)

`HelmholtzFe::new_quadratic` retains one physical midside point per canonical
mesh edge and the original hierarchical H1 pressure coordinates. Shared edges
give adjacent elements the same quadratic face map. A nonsingular corner
tetrahedron is required; the existing `fs-query` interval/Bernstein TET10
orientation owner must additionally admit the complete closed curved element.
Sampling only integration points never admits geometry. Stored negative corner
orientation is reordered for the proof alone. Query work/storage exhaustion,
unresolved geometry, nonfinite metrics and cancellation refuse explicitly.
This is local orientation evidence, not global injectivity or mesh-overlap
certification.

One shared map supplies physical position, pointwise volume determinant and
inverse-transpose barycentric gradients. The same cofactor supplies pointwise
outward face normals and surface measures. Original complete stiffness/mass,
heterogeneous material sampling, linear and prepared quadratic velocity loads,
Robin conductance, mechanical-wall trace mass and viscous tangential surface
stiffness all retain this geometry. Tangential forms discard the omitted
coordinate's normal derivative exactly. Existing modal, harmonic, sparse-time,
mechanical and complete-field exact solvers consume the resulting original
forms without a new pressure reduction or propagation engine. Diffusion and
joint caloric owners likewise reuse these original geometric forms.

Geometric quadrature is explicit in `order+2..=16`; curved default face forms
honor the same minimum. Geometry, pressure order, material/source quadrature
and time remain independent convergence axes. Curved stiffness and surface
measures are generally nonpolynomial, so a quadrature count is not an accuracy
certificate. The original affine constructor keeps its established assembly
arithmetic; its original metrics are cached once. The curved constructor does
not assemble an affine replacement pencil first. Point geometry requires no
allocation or repeated orientation proof.

Point pressure/velocity/intensity probes and distributed volume sources sample
the original mapped physical position. Surface observations use the pointwise
normal and measure, keeping harmonic active/reactive power and transient
signed flux distinct. `probe_at_position` uses bounded projected Newton trials
and backtracking with explicit SI residual tolerance in a caller-selected
original element; failed location publishes no probe or material query.
Prepared stencils and source columns retain geometry-free repeated evaluation.
`quadratic_velocity_source_with_quadrature` permits independent source-metric
refinement. Immutable `support_faces` and `has_nonzero_projection` preserve
physical support, including zero-valued faces, and distinguish a live signed
zero-net shape from a zero projected source.

Legacy APIs returning one triangle and one normal admit only genuinely affine
selected faces, including planar openings in otherwise curved interiors.
Nonaffine selected faces refuse instead of exporting their corner triangle.
Resolved LNS consumes this same map through the complete mixed-fluid extension
specified below. This geometry owner does not imply curved BEM panels, a
curved-mesh importer or closed-cylinder/instrument qualification. The broader
bead stays open.

Twelve native regressions specify analytic volume/gradient moments, all 24
tetrahedron vertex orders, shared curved-face reciprocity, rigorous fold and
subdivision-budget refusal, affine recovery through P4, physical material
sampling, exact-field evolution, bounded physical inverse locations,
variable-Jacobian source moments and curved surface/source/viscous quadrature
ladders. Pinned Rust parsing, targeted formatting, dependency-direction checks
and independent analytical moment/metric/inverse oracles passed. Native Rust
compilation and tests remain unrun under the existing upstream dependency
memory failure and full workspace; these source checks do not establish
numerical test passage.

## Quadratic geometry in the complete resolved viscous/thermal fluid

`ThermoviscousFe::from_geometry`, `from_geometry_with_thermodynamics` and
`from_geometry_with_moving_wall` consume the original admitted P2 Helmholtz
geometry. They retain its full-element orientation proof, original topology,
midside map and declared quadrature, and require the same equilibrium density
and adiabatic sound speed. The complete P2/P1/P2 velocity, pressure and
temperature coordinates remain. No replacement mesh, extra geometry proof,
pressure reduction or second fluid operator is introduced.

Variable volume determinants and physical inverse-transpose gradients enter
every original inertial, viscous, conductive, pressure/divergence and
thermodynamic storage block. Independent loss integration uses that same map.
The existing full-vector moving-wall congruence therefore retains pressure,
shear, inertia and tangential no-slip work on curved reference faces. The
original affine constructors retain their established arithmetic. Element
scratch buffers are reused across quadrature points instead of allocating
value/gradient vectors at each point.

Original heat-source callbacks receive mapped physical positions and variable
volume weights once during preparation. Their independently requested order
must meet the geometry's declared minimum. Signed source moments, cached full
temperature loads, exact model identity and existing harmonic/time work pairing
remain unchanged. Original molecular reservoirs use the same full mapped P2
mass and quadrature as the resolved thermal field; they are not collapsed to
bulk viscosity or a fitted heat capacity.

Finite conductive wall assembly reuses the original mapped boundary mass.
Independent pointwise wall heat uses cached physical surface Jacobians at the
same original local-face quadrature points, with explicit permutation into
canonical temperature traces. This preserves work even when stored tetrahedron
and face orders differ. The solid's existing local normal-profile model still
requires its declared curvature and tangential-scale admission; this does not
supply a curved solid-volume mesh.

Pressure traction ports, conforming planar interfaces and spatial slab owners
continue to refuse nonaffine selected faces where their original physical
representation requires a plane. Curved wall geometry does not silently become
a flat corner triangle. Original mixed-space, dense-memory, mesh-aspect,
constitutive, cancellation and source-identity refusals remain active. The
aspect limit measures corner tetrahedra; curved determinant quality belongs
to the consumed geometry's whole-element proof. Geometry, quadrature and the
mesh resolving viscous/thermal layers remain separate convergence obligations.

Eight focused native regressions specify affine matrix/load identity;
manufactured curved stress, storage and divergence; driven reciprocity and
power; full-vector moving-wall reaction; mapped heat-release moments; analytic
coupled frozen-Cv/molecular exchange; and conductive surface quadrature/work
under both tetrahedron orientations, plus admission refusals. Pinned Rust
parsing/format checks and independent analytical moment/surface calculations
were performed. Native compilation and test execution remain unverified under
the existing upstream dependency-memory failure and full workspace. No
closed-cylinder, instrument, spatial-resolution or broader bead acceptance is
inferred.

## Spatial resistance in broadband pressure and wall evolution

`HeterogeneousHelmholtzFe::spatial_resistive_boundary` prepares the local,
frequency-independent law `v_out = Y(x) p`, with nonnegative specific
conductance in m/(Pa s). It uses the original affine or quadratic surface map,
normal and complete pressure trace at explicit quadrature points. The retained
form is `C_ij = integral(Y(x) phi_i phi_j)`, preserving spatial pressure/material
correlation. No face-average conductance or replacement surface is introduced.
This typed law is distinct from a frequency-local `SpatialBoundaryAdmittance`:
real admittance at one frequency alone is not a causal broadband law.
Its `at_frequency` method reuses the same prepared form for harmonic fields
and their peak-phasor power without another material sampling pass.

Sparse midpoint, complete-field exponential, and both causal mechanical-wall
evolutions accept this same prepared form. It enters the dynamics and the
independent boundary-heat calculation: `h p_mid^T C p_mid` for midpoint and
`integral(p^T C p dt)` for exact held-input evolution. Mechanical resistance,
Maxwell relaxation and reversible fluid/wall work remain separate. Every
original pressure and wall coordinate remains present. Warm steps retain their
existing cost and do not resample material or geometry.

Uniform, spatial, mechanical and pressure-release patches must be disjoint;
inward sources may share non-release faces. Foreign-model fields, nonfinite or
active samples, invalid quadrature and cancellation refuse without publishing
partial state. Existing uniform paths preserve their arithmetic when the
spatial patch list is empty. G1/G3 regressions use exact triangle moments for
nonuniform wall heat, uniform-limit recovery, independent midpoint/exponential
refinement, and separate coupled fluid/wall work balances. The law is prescribed,
linear and memoryless; dispersive materials require their own causal states.
Geometry, material quadrature, spatial and temporal convergence remain distinct.

## Causal spatial admittance with original boundary memories

`HeterogeneousHelmholtzFe::spatial_relaxing_boundary` composes explicitly
authored nonnegative spatial fields into the positive-real law
`Y(x,s) = Y0(x) + sum a_j(x)/(1+s tau_j)`, with one to 64 finite positive
constant relaxation times on each physical patch. Every coefficient retains
the same original geometry/material owner, canonical faces and quadrature.
The fields need not be proportional. A negative residue, arbitrary frequency
fit, spatially averaged pole or complex resonance is not inferred. The owner
represents parallel resistive/inertive branches, not every possible wall law.
Harmonic reuse retains the full `1/(1+i omega tau_j)` response in the original
interior convention, including its negative susceptance, without resampling.

`time::relaxing` evolves every original pressure-potential/pressure coordinate
and every geometric trace memory `tau_j q_dot + q = p`. Zero or semidefinite
residue forms are admitted without inversion or coordinate removal. Adjacent
pressure-release faces constrain fluid fields only: independent memories at
shared trace DOFs remain live and retain their original free-row coupling.

Midpoint eliminates these states into ONE original sparse factor with
`Ceff = C0 + sum h/(2 tau_j+h) Cj`; its history enters the original projected
load with a negative sign. No factorization, material callback or surface
assembly recurs during warm steps. Original linear/quadratic boundary motion
and volumetric sources retain their physical amplitudes and work pairings.
The condensate's algebraic loss is never published as physical heat.
Independent accounts retain memory storage `sum tau_j q_j^T Cj q_j/2`,
instantaneous heat, relaxation heat and signed fluid-to-memory exchange.
All three fluid, memory and total work balances and the original uncondensed
equations are checked before a candidate state is returned.

The seven focused regressions include independent triangle moments for
spatial heat/storage, harmonic/trapezoidal consistency, a separately specified
three-state ODE refinement oracle, energy return from release-shared memories,
zero-residue geometric support, source work, limits, ownership and cancellation.
Midpoint is passive but requires a time-resolution axis for stiff poles and
wave phase. The complete-state held-source exponential path is available in
`time::relaxing::exact`, as specified below. Combined mechanical and spatial
admittance memories are available through the mechanical owners' explicit
`prepare_*_with_walls_and_relaxing_boundaries` methods.

### Analytic dispersive mechanical walls in full-pressure resonance queries

`BoundaryAdmittanceFrequencySample` binds each uniform complex specific
admittance and its analytic derivative per rad/s to the same declared frequency
in Hz. `prepare_with_dispersive_admittances` admits the primary passive law on
original, disjoint faces and retains its signed sparse derivative forms. The
existing coherent boundary/volume-source tangent now includes
`(B + omega B') p - i (G + omega G') p`, reusing the original LU, pressure-release
elimination, double-double forcing, residual correction and cancellation. Only
the primary conductance is a heat form. Existing preparations still freeze
their admittances; bulk density, compressibility and geometry remain fixed.

`BoundaryAdmittanceFrequencySample::mechanical` evaluates the existing
`MechanicalBoundary`'s original `WallPin` and every positive Maxwell arm:
`Z = r + i(omega sigma - K/omega) + sum k/(1/tau + i omega)`,
`Y = 1/Z`, `Y' = -Z'/Z^2` under `exp(+i omega t)`. Double-double spring/mass
cancellation and rate-form complex divisions avoid unnecessary high-Q loss of
precision and `omega*tau`, `tau^2` or denominator-square overflow. This exact
local harmonic condensation requires finite nonzero impedance. In particular,
a lossless wall-impedance zero refuses even if the uncondensed coupled problem
is regular; it is never clipped or supplied artificial loss. Reported field
energies remain fluid energies, without a claim of complete wall storage.

`survey_pressure_resonances_with_mechanical_walls` resamples these original
wall laws and their analytic slopes at every shared-survey query. Fixed
admittances and coherent velocity phasors remain constant. Original point
observations, limits and residual gates are retained. The result is a
grid-resolved driven-peak survey, not a free-mode census, complete coupled-field
qualification, half-power width or decay Q. Local condensation failures and
cancellation publish no partial survey.

Native validation executed the complete production pressure/tangent module
against real `fs_math`, `fs_sparse` and complex-LU owners: 10 tests passed. New
tests compare pressure and wall-velocity slopes to independent uncondensed
cavity/displacement/Maxwell state equations, cover high-Q behavior and lossless
local-pole refusal, and retain signed wall slopes, coherent source derivatives
and exact pressure-release zeros. Two focused FE adapter/survey regressions are
authored; they are not claimed executed because the full `fs-feec` dependency
build is unavailable in this environment. No broader IRP bead is closed.

### Complete spatial causal-admittance exponential trajectories

`time::relaxing::exact` prepares a held-input exponential trajectory for every
original free pressure-potential/pressure coordinate and every geometric trace
memory of the admitted `SpatialRelaxingBoundary` laws. Preparation retains
the original consistent fluid mass, stiffness, instantaneous conductance and
all spatial residue matrices. Each raw memory obeys `tau_j q_dot + q = p`;
the coupled fluid equation contains the original `Cj q_j` reactions. No
residue matrix is inverted, diagonalized, fitted or assumed positive definite.
Zero residues and pressure-release-shared memory coordinates remain live.

The complete physical storage includes `sum tau_j q_j^T Cj q_j/2`. Admission
independently checks the full generator identity
`H A + A^T H + 2(C0_heat + Cj_heat) = 0`. Additional coordinate scales serve
only numerical conditioning; they do not invent physical storage in a residue
null space. Original sparse matrices independently evaluate endpoint energy
changes and the resulting acoustic, memory and total balances.

One cached `PreparedQuadraticFlow` trajectory supplies endpoint states,
instantaneous heat, relaxation heat, signed fluid-to-memory exchange, and
integrated pressure observations for every held physical source profile.
Linear velocity ports, prepared boundary sources and volume sources retain
their original load vectors, units and distinct source-work subtotals. Heat
is an independently integrated trajectory form, never an assigned energy
deficit. Candidate states publish only after the physical balances pass.

Preparation is a bounded dense operation over the complete state and held
inputs. Every memory and augmented coordinate counts against the declared
limits. Warm steps reuse cached endpoint and observation maps without material
callbacks, surface assembly or factorization. Exponential accuracy controls
remain separate from spatial/quadrature convergence and source hold accuracy.
This path handles the stated linear, time-invariant laws; it does not imply
nonlinear boundary physics, arbitrary frequency-law fitting or unbounded
large-mesh dense scalability.

Five focused tests cover an independently specified physical ODE and its
integrated heat/work, interval-partition parity, spatial midpoint refinement,
stiff zero-residue release-shared memories, coherent boundary/volume-source
work and no resampling, and bounded admission, ownership and cancellation.
These regressions are authored, not claimed executed: the native dependency
build exhausted the shared execution environment, which subsequently became
unavailable before the final formatting command could be confirmed.

### Simultaneous mechanical and spatial-admittance wall dynamics

`prepare_time_with_walls_and_relaxing_boundaries` and
`prepare_exact_time_with_walls_and_relaxing_boundaries` compose the original
mechanical/Maxwell wall states and spatial causal-admittance memories against
ONE complete acoustic evolution. Every original displacement, velocity,
Maxwell memory and admittance pressure memory remains explicit. Physical
patches have disjoint faces; independent histories at shared edges, including
pressure-release-shared coordinates and zero residues, remain live.

Sparse midpoint eliminates both boundary laws into the same original acoustic
factor and reconstructs both histories from the same midpoint pressure. The
combined memory cap counts the complete Maxwell and admittance trace extents.
The exact path retains the existing fluid/mechanical pHS coordinates and adds
raw admittance memories to one complete-generator exponential cache. No
residue matrix is inverted or assigned fictitious positive storage. Dense
admission counts every original state and source; its physical identity audit
includes all four heat channels and excludes both reversible exchange forms.

The public accounts distinguish fluid storage, mechanical kinetic/elastic and
Maxwell storage, and `admittance_memory_j`. Work separates instantaneous Robin
heat, Kelvin resistance heat, Maxwell heat, admittance heat, mechanical
exchange, admittance exchange, and the original signed source work. Acoustic,
mechanical, admittance and total balances are independently checked before a
candidate is returned. Admittance q [Pa] is never confused with the mechanical
energy-normalized Maxwell coordinate z. `state_with_all_memories` imports both
histories explicitly; legacy wall-history imports refuse omitted admittance
histories. The explicit zero/rest constructors initialize every history to zero.

Existing preparation entry points retain their empty-admittance behavior.
Warm steps require no new factorization or material sampling. The exact owner
keeps the original four observation forms when no admittance arms are present,
and six forms plus the same source observations when they are present. Source
work remains paired to original linear/quadratic boundary and volume profiles.

Eight focused regressions cover independent surface moments for all physical
accounts, an independent six-state physical ODE with seven integrated accounts,
full interval-partition parity, coupled spatial midpoint refinement, energy
return from pressure-release-shared histories into fluid and mechanical walls,
zero residues, complete memory/state limits, mixed signed sources, ownership
and cancellation. Pinned selective formatting and whitespace checks passed.
Independent numerical reference calculations support the test fixtures; the
Rust regressions are not claimed executed while shared native-build resources
remain exhausted. These are the stated linear local wall laws; they do not
add shell bending, a moving fluid mesh or nonlinear wall constitutive physics.

### Complete harmonic pressure and causal mechanical wall states

`prepare_with_coupled_mechanical_walls` retains the original pressure space,
every original wall displacement and outward velocity, and every Maxwell memory.
The original consistent face mass `W=integral(phi_i phi_j)` enters both
`+i omega W v` in the fluid pressure equation and `-W p` in wall momentum.
No face-average piston, diagonal surface mass, modal reduction or reciprocal
wall impedance replaces these equations. Exact pressure-release elimination
applies only to fluid pressure; all adjacent wall coordinates remain present.
Uniform memoryless `ResistiveBoundary` and original `MechanicalBoundary`
inputs preserve the time owner's material admission and physical units.

The complete sparse blocks form `A0+i omega A1-omega^2 A2`. One bounded
original complex LU, with the existing equilibration, double-double residual
correction and analytic frequency-tangent owner, solves the complete system.
Limits count every free pressure, displacement, velocity and memory coordinate,
with the existing harmonic hard cap. Original coherent linear, prepared
boundary and volume-source projection, source derivatives and signed source
work are reused. Derivatives are per fundamental angular frequency in rad/s.

The invertible internal coordinate `e=sqrt(k)w-z` preserves small physical
relaxation extensions without removing a memory. Original `z` and `dz/domega`
are returned; at large `omega*tau` the original memory ODE recovers them without
subtracting nearly equal phasors. Independent physical accounts report wall
kinetic, equilibrium elastic and Maxwell storage; Kelvin and Maxwell heat;
and signed `Re(p^H W v)/2` exchange. Stored energies use the peak-phasor factor
1/4 and powers 1/2. Total and separate fluid/wall power residuals are observable,
without defining heat from an energy or input-power deficit.

A local lossless `Z=0` no longer causes refusal when the complete coupled
operator is regular. Genuine coupled singularities still fail the native
factor/residual gate; no artificial loss or discarded dark wall state hides
them. `survey_pressure_resonances_with_coupled_mechanical_walls` uses the existing
shared peak search with fixed coherent sources and original point observations.
It can cross a regular local `Z=0`; its grid-resolved peaks do not certify a
complete mode census, linewidth or decay Q. The older condensed API retains its
explicit finite-nonzero-impedance scope. Bulk medium/geometry and the declared
local wall coefficients remain frozen; lateral shell bending, nonmatching mortar
and exterior radiation are not inferred.

Native validation executed the actual complete numerical modules against real
`fs_math`, `fs_sparse` and complex-LU owners: 14 tests passed. Four new tests cover
regular local `Z=0` motion/storage/tangent, genuine coupled singularity refusal,
non-diagonal surface mass and an independent original-z Maxwell formulation
with every storage/heat/work channel (including nearly frozen memory recovery),
and complete-state limits, release constraints and cancellation replay. Two
public FE regressions cover local-pole survey traversal and P2 equivalence to the
admitted condensed equations; they are authored, not claimed executed while
the full `fs-feec` dependency build is unavailable. MR56 and broader IRP
qualification beads remain open.

### Distributed applied traction on coupled acoustic walls

`time::mechanical::traction::WallTractionProfile` supplies a dimensionless
spatial force shape in the ORIGINAL GLOBAL hierarchical wall trace coordinates.
Its amplitude is pressure [Pa], positive outward like the fluid pressure on
the wall. Omitted coefficients are zero; repeated, nonfinite or off-patch
coordinates refuse. Pressure-release constraints remove only fluid pressure
coordinates, never adjacent mechanical traction or memory coordinates.

`PreparedMechanicalTime::prepare_wall_traction` binds the shape to the original
wall/geometry preparation. `step_with_wall_tractions` combines coherent applied
forces in DD before the same sparse midpoint solve and reconstructs every
original displacement, velocity and Maxwell memory. The physical equation is
`sigma*v_dot + r*v + K*w + sum sqrt(k_j)*e_j = p + traction`. The original
consistent face mass projects both force and conjugate velocity. There is no
prescribed wall motion, patch-average piston, additional fluid volume source
or alternate mechanical solve. Prepared shapes occupy O(wall trace DOFs);
warm forcing requires no surface assembly or factorization.

`prepare_exact_time_with_wall_tractions` appends held Pa traction columns after
the original boundary and volume inputs in the complete exponential owner.
The native wall momentum receives W*shape; independently assembled observations
integrate shape^T W*v over the same trajectory. Its `step_with_wall_tractions`
returns `integrated_wall_traction_displacement_m3`: amplitude times this volume
displacement is external work. All traction columns count against the native
augmented-state limits; missing prepared amplitudes refuse. Spatial-admittance
memories remain part of the same complete evolution.

Both owners retain `acoustic_input_work_j` and `wall_traction_work_j` separately,
including negative mechanical work. Total input includes each once;
`wall_work_j` remains solely internal fluid-to-wall exchange. Weak acoustic work
is not recovered by subtracting two large total/mechanical values. Fluid, wall,
admittance and total physical balances are checked before candidate publication.
Original locally reacting wall scope remains: no lateral elastic propagation,
shell bending, moving fluid mesh or nonlinear traction law is inferred.

Eight focused regressions specify a continuum translating-wall solution,
zero-net-force spatial excitation with both memory families, original P2
surface work and weak coherent force cancellation, an independent physical ODE
with nine work/heat/displacement accounts, interval partition, full-state time
refinement, source units and weak acoustic work, pressure-release-shared forcing,
ownership, extent and cancellation. Selective pinned formatting and whitespace
checks passed, and reciprocal source/physics review found no concrete blocker.
Native Rust compilation and test execution remain unavailable in the exhausted
shared build environment; these regressions are authored, not claimed executed.

### Harmonic applied traction and analytic frequency response

Complete harmonic mechanical walls accept reusable owner-bound
`WallTractionProfile` projections in the original global hierarchical trace,
including pressure-release-shared mechanical coordinates. Complex peak traction
amplitudes are Pa, positive outward; their actual weak momentum RHS is
`W * shape * amplitude`. `PreparedMechanicalHarmonicTraction` retains this
signed consistent projection, so repeated responses reuse the original complete
LU without face reassembly, impedance inversion, modal condensation or added
damping. The same complete pressure, displacement, velocity and Maxwell memory
coordinates remain present at a local wall-impedance zero.

`solve_with_wall_tractions` combines the original fluid sources and mechanical
forces coherently. `solve_with_wall_tractions_frequency_derivative` accepts
explicit Pa/(rad/s) traction-amplitude slopes. External forces and their slopes
remain in double-double through the direct-RHS solve, analytic operator tangent
and original-equation residual correction. Mechanical force is never divided
by frequency to disguise it as a fluid velocity source.

`wall_traction_power_w` is the independent signed
`Re((W * traction)^H velocity)/2` account. `acoustic.input_power_w` remains
solely fluid-source work, and `fluid_to_wall_power_w` remains solely internal
pressure exchange. Total and wall balances include external mechanical work
once; the fluid balance excludes it. Existing fluid-only APIs delegate with
no traction and retain their behavior. Prepared columns occupy O(wall trace
DOFs); repeated forcing performs no surface assembly or factorization.

Profiles and spectra refuse invalid coordinates, duplicate or nonfinite data,
foreign owners and cancellation. Active traction counts retain the existing
128-source bound. The original locally reacting wall scope is unchanged; this
does not infer shell bending, lateral elastic waves or nonlinear forcing.

Four focused regressions specify independent P1 acoustic/Maxwell equations and
analytic source slopes, P2 surface work and strong coherent cancellation,
regular coupled response at a local lossless wall-impedance zero, and source
admission, binding and cancellation. Bounded independent source/physics review,
pinned formatting and whitespace checks passed. Native Rust compilation and
these four regression groups remain unexecuted in this environment.

### Complete material damping on resolved thermoviscous moving walls

`ThermoviscousFe::with_moving_wall_damping_channels` retains complete passive
material/joint channels in the original mass-normalized structural coordinates.
The admitted operator is `C = sum l_k l_k^T`, including signed cross terms and
rank-deficient loss null spaces. Harmonic impedance, causal dissipation, and
structural heat use this same operator; harmonic heat evaluates the original
channel amplitudes directly with double-double dot products. Fluid reaction
continues to use only the fluid operator, and resolved viscous, thermal, and
caloric losses remain separate from structural heat.

The caller must provide channels in the wall's original physical modal basis;
a matching mode count alone does not establish that identity. At most 64
complete channels are admitted. Empty channels explicitly select zero material
dissipation. Invalid extents, nonfinite entries, nonrepresentable Gram entries,
nonzero existing modal damping, duplicate attachments, and cancellation refuse.
Omitting the attachment preserves the existing diagonal modal-damping behavior.
The Gram matrix is prepared once; no material query or reassembly enters a
harmonic solve or causal step.

Four focused regressions specify common-motion loss null spaces, differential
motion with an independent damped solution, transient integrated heat and
energy balance, signed coordinate changes, and admission/cancellation. Source
and reciprocal-physics review plus pinned formatting and whitespace checks
passed. Native Rust execution remains unverified because the shared dependency
build exhausted the bounded validation environment.

## Original structural relaxation in resolved moving walls

`with_moving_wall_relaxation_channels` attaches every original signed coupling
row c [1/s] and relaxation time tau [s] in the same mass-normalized wall basis.
The structural modes contain RELAXED stiffness. Each retained energy coordinate
x [sqrt(J)] obeys x_dot = c dot v - x/tau, applies resisting force c*x, stores
x*x/2, and dissipates x*x/tau. No channel, off-diagonal term, or supplied history
is fitted, discarded, or converted into frequency-local diagonal damping.
Separate instantaneous structural damping may coexist without counting these
same material losses twice. One complete memory owner is admitted, with explicit
channel, projection-work, and backing-storage limits.

Harmonic preparation caches the full complex matrix c*c^T/(1/tau+i*omega)
once per frequency. The original solve and unscaled residual use that same cache;
reported fluid traction still excludes structural forces. MovingWallField reports
the complex material coordinates, mean stored energy |x|^2/4 and mean heat
|x|^2/(2*tau); ThermoviscousPower includes structural_relaxation_w once.

Causal preparation appends every original x to the same full mixed-field state
and original inverse-mass congruence. Explicit initial preload and dark-channel
history survive; zero_state declares relaxed rest. The existing exact held-input
trajectory integrates structural_relaxation_heat_j, including mechanical,
heat-source, volume-source and two-fluid interface paths. Complete state/port
admission refuses oversized models instead of truncating source material history.
The source owner remains responsible for its original validity band and for
excluding these branches from its instantaneous loss model.

## Free damped cavity and mechanical-wall resonances

`PreparedMechanicalHarmonic::resonance_pole` and `resonance_spectrum` query the
same full physical pressure/wall/Maxwell polynomial used by native forced
response: `T(s) = A0 + s*A1 + s^2*A2` under `exp(s*t)`. Original heterogeneous
fluid mass/stiffness, real-Robin dissipation, consistent wall trace inertia,
reciprocal pressure work, and every original Maxwell extension remain present.
No wall impedance or material-memory denominator is inverted; an isolated
material relaxation rate remains available when its original full state is a
mode. The reference real frequency required to prepare the existing harmonic
owner changes only fixed numerical equilibration, not these coefficients.

Local refinement and contour discovery use only the existing shared `fs-la`
Newton/Beyn owners. Physical seeds follow free pressure and then each wall's
original displacement, velocity and Maxwell-extension blocks. Results lift
pressure-release coordinates to exact zeros and retain displacement, velocity,
original Maxwell memory and extension at one arbitrary common modal amplitude.
Each original free equation passes an independently evaluated, coefficient-row
backward residual using double-double sums. The returned signed rate is not
clipped; pure real relaxation and modes invisible to a pressure observer are
retained without assigning them a forced peak or a half-power Q.

Shared solver limits cover numerical dimension, work and payload. Separate
`MechanicalPoleLimits` bound sparse callback/audit evaluations and callback plus
result payload, and set the original-row residual tolerance. Cancellation and
any refused mode publish no partial spectrum or accepted-state mutation.

These are numerical modes of the admitted finite spatial/material realization.
Nested contour quadrature, count and rank agreement do not prove continuum
completeness. Exact zero-rate pressure/gauge classification is excluded, and
the requested contour must exclude zero. Repeated or defective roots may
refuse through the common solver. This route neither extrapolates a measured
material validity domain nor certifies experimental decay or receiver audibility.

Focused G1/G3/G4 tests specify a manufactured coupled cavity/Maxwell decay,
reference-frequency invariance, an exact unobservable equal-time memory pole,
source-independent contour discovery, budgets and cancellation. The public FE
adapter is checked against independent P1 tetrahedron volume and surface
integrals with pressure-release constraints and material-dependent Robin decay.

### Complete cavity source participation

`PreparedMechanicalHarmonic::prepare_resonance_response` binds a reusable
simple-pole factor to the original pressure, mechanical-wall and Maxwell
realization. It rechecks the supplied physical coordinates, original row
residuals and pole correction; an unrelated pole report cannot bypass those
checks. The shared `fs-la::pole_response` owner supplies the nonsymmetric
bordered factor and cancellable triangular solves.

`PreparedMechanicalResonance::forcing_residue` projects the original coherent
inward velocity, boundary injection, volume injection and applied wall traction
sources. Amplitudes are their analytic values at the complex Laplace pole:
fluid forcing is `s*b`, mechanical forcing is the consistent `W*traction`.
The original source/traction binding and pressure-release rules remain active.
Numerical congruence scales only the solve; returned pressure, wall motion and
every Maxwell memory/extension residue use original coordinates and units.
Repeated queries reuse the factor at quadratic cost.

The returned Laurent residue is one pole's contribution, not a driven field,
source power or linewidth. With the native `exp(+i*omega*t)` convention its
angular-frequency residue is `-i` times the Laplace residue. Coherent source
cancellation and dark material modes can have zero participation. Repeated or
defective poles can refuse; no local wall impedance is inverted, no modal
damping is inferred, and no real-frequency-only source is analytically
continued without caller-supplied data.

### Original wall-material resonance sensitivities

`PreparedMechanicalResonance::wall_parameter_derivative` differentiates the
complete cavity/wall/Maxwell polynomial along explicitly supplied original
surface-density, stiffness, resistance, branch-stiffness and relaxation-time
directions. One existing bordered factor supplies the signed pole and all
physical mode-coordinate derivatives. Original Maxwell-memory recovery
includes the square-root stiffness coordinate's chain rule, including dark
material modes; subnormal stiffness directions retain representable changes.

Numerical congruence and the reported physical normalization pivot stay fixed.
The full original differentiated equations are independently audited with
double-double products. Additional physical payload is O(N), recovery/assembly
is O(N+nnz), and the reused solve is O(N^2); `MechanicalPoleLimits` bounds the
extra payload and one physical audit. Cancellation or refusal returns no
partial tangent.

The caller supplies coefficient derivatives per its own parameter. Geometry,
trace operators, bulk fluid and source laws stay fixed. The API admits a local
simple-pole derivative and does not infer a temperature law, material validity
extension, finite design step, repeated-mode sensitivity or driven linewidth.


## Complete physical-drive cavity and wall resonance searches

`HeterogeneousHelmholtzFe::driven_resonance_with_coupled_mechanical_walls`
resolves an observed maximum and optional two-sided half-height width for the
complete original pressure, distributed wall and Maxwell-memory equations.
`MechanicalResonanceDrive` combines fixed coherent velocity traces, owner-bound
boundary and volume velocity sources, and original `WallTractionProfile`s.
Each actual-frequency preparation binds the traction to its own full-system LU;
that factor is reused for the complete analytic frequency derivative. Fixed
fluid injections still differentiate their `i*omega*b` RHS. Fixed applied
tractions remain direct wall momentum forces, with no reciprocal impedance.

`MechanicalResonanceObservation` selects pressure at a physical point, a signed
coherent sum of up to 128 physical pressure samples, or dimensionless weights on
one wall's original displacement or normal-velocity coefficients. Source and
receiver cancellation occurs before amplitude powers. The matching
`survey_resonances_with_coupled_mechanical_walls` discovers grid-resolved peaks
under the same experiment; it does not certify a complete mode census.

The original geometry/material snapshots, source ownership, pressure-release
constraints, work accounting, cancellation and resource limits remain active.
Each caller-supplied local interval must isolate a maximum and, if requested,
both half-height crossings. Silent modes, notches, singular samples and missing
crossings do not become invented peaks or widths. These are observation-specific
amplitude resonances, not pole-decay Q or emitted acoustic power.

Two native coupled/search tests cover independent pressure, displacement and
velocity peak/width limits, hereditary loss, and coherent cancellation. Two
public FE integration tests additionally exercise applied traction, integrated
volume injection and foreign-owner refusal.

## Coherent cavity and wall antiresonance bands

`HeterogeneousHelmholtzFe::driven_antiresonance_with_coupled_mechanical_walls`
uses the same immutable physical drive and receiver as the peak adapter, with
the shared zero-safe minimum/recovery search. Fixed integrated fluid injections,
applied wall tractions and signed pressure/wall observations retain their units,
original source ownership, coherent phases and spatial projection. Every actual
frequency reuses one complete pressure/wall/Maxwell factor for the response and
its analytic tangent; recovery crossings skip the tangent. No pressure receiver
or wall impedance is inverted and no internal motion or memory is discarded at
a cancellation. Explicit recovery levels use the observation's squared units
(Pa², m² or m²/s²), require both crossings, and imply no material or free-mode Q.
Original residuals, release constraints, cancellation and query limits remain.

The added public P1 tetrahedron regression independently derives a coherent
volume/traction cancellation at 1.7 rad/s and its two recovery edges. It checks
that a zero pressure observation retains all three wall coordinates, positive
stored wall energy and 0.003625 W of actual traction work and wall heat. Its
rational oracle and derivative were independently evaluated. The public FE test
is authored but unexecuted while the full dependency build remains unavailable;
shared search and original rotor/wave numerical tests are separately verified.
