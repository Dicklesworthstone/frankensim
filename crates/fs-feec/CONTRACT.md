# fs-feec CONTRACT

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
  enables the optional `fs-blake3` and `fs-couple` dependencies.
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

`tests/feec_battery.rs` (10 cases, JSON logging): dd = 0 on the
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
