# fs-feec CONTRACT

## Purpose and layer

Layer: **L3 FLUX** (deps: fs-rep-mesh L2; fs-la, fs-sparse L1;
fs-math, fs-qty, fs-soa L0). The exterior-calculus core (plan §8.1,
Bet 3): cochains on tet complexes, Whitney P₁Λᵏ forms, discrete Hodge
stars, and exact-sequence operator assembly. The exterior derivative
is fs-rep-mesh's INTEGER incidence operator — dd = 0 is exact by
construction (their contract invariant, re-verified here on the zoo);
every metric statement and every approximation lives in the mass
matrices and stars. This is what kills spurious pressure/EM modes and
checkerboarding structurally instead of by stabilization folklore.

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

## Error model

Structured panics on shape mismatches, out-of-range degrees, and
degenerate tets (programmer/mesh errors). `betti` panics on i128
overflow rather than degrade exactness. No Result-based paths in v1 —
the inputs are meshes that upstream certificates (fs-rep-mesh,
fs-mesh audits) have already validated.

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

## Cancellation behavior

Bounded synchronous assembly loops; chunking a large mesh to tile
quanta with Cx poll points between chunks is the fs-exec driver's
job (the fs-la/fs-simd discipline). No internal threading.

## Unsafe boundary

None. `unsafe_code = "deny"` via workspace lints.

## Feature flags

None.

## Conformance tests

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
`tests/highorder_battery.rs` (slice 1): Gauss–Legendre exactness to
degree 2n−1 (n = 1..10) + node symmetry; Lobatto endpoint structure;
sum-factorized apply vs the dense assembled Kronecker reference to
1e−12 relative on four (m, r) fixtures — the acceptance roundoff
gate; Jacobi diagonal vs operator columns; G1 MMS Poisson through the
matrix-free Jacobi-PCG path with slope gates ≥ r + 0.6 for r = 1..6
(measured ≈ r + 1); its own golden hash.
`tests/ho_probe.rs`: per-mode convergence regression — the diagnosis
that single-cell symmetric fixtures superconverge at even r (a metric
trap, so MMS ladders start at m ≥ 2).
`tests/simplex_battery.rs` (slice 3): Duffy weight/volume exactness
and P_r dimension counts (local dofs = C(r+3,3), r = 1..6);
unisolvence via mass SPD (r = 1..5); conformity across the shared
two-tet face (traces from both elements agree ≤ 1e−13 at sampled
face barycentrics); G1 MMS r = 1..4 with slope gates ≥ r + 0.6
(ladders [4,8] for r = 1, [2,4] above); ORIENTATION battery (G3):
signed-permutation operator equivariance at r = 3 (edge P₁ kernels
flip sign under sort-order reversal, K′(Su) = S(Ku) ≤ 1e−12) and
physics invariance at r = 4 under random relabeling (vertex point
values ≤ 1e−9, L2 error relative deviation ≤ 1e−6 — face bases mix
under relabeling, so operator-level identity is the WRONG gate
there); its own golden hash.
`tests/derham_battery.rs` (slice 2): curl∘grad and div∘curl ≤ 1e−13
relative on four (m, r) fixtures; exact-sequence dimensions
(χ = 1 for m = 1..3 × r = 1..6); commuting diagram 1D
(G·π_C f = π_D f′ ≤ 1e−11) and 3D on product fields; projection
G1 ladders for both 1D families at r = 1..6 (C: order ≥ r + 0.6
gate, measured ≈ r + 1; D: ≥ r − 0.4 gate, measured ≈ r) — these
drive all four 3D tensor space types' rates; Legendre mass closed
form; its own golden hash.

## No-claim boundaries

- CIRCUMCENTRIC diagonal star deliberately absent: on
  non-well-centered meshes (including the Kuhn fixtures) it produces
  negative dual measures, so shipping it without well-centeredness
  machinery would be a certificate without evidence. Follow-up scope
  together with mesh quality certificates.
- Simplicial H¹ now reaches high order (slice 3, MMS-gated at
  r = 1..4; the basis is valid to r = 6 by the unisolvence/count
  gates). Simplicial H(curl)/H(div)/L² high-order families (the
  Nédélec/RT tree) remain OPEN scope — split assessment pending. The
  tensor-product side covers all four space types (slices 1–2); full
  3D VECTOR MMS solves need tfz.10 — projection ladders stand in,
  labeled. Unstructured-hex orientation and the ≥30%-peak perf gate
  are later-slice scope. `HexComplex` incidence is still not
  consumed.
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
