# fs-adjoint CONTRACT

## Purpose and layer

Layer: **L3 FLUX** (deps: fs-solver/fs-feec L3, fs-ad/fs-la/fs-sparse
L1, fs-math L0, fs-rep-mesh L2). GRADIENT TRUTH (plan §8.7, Bet 4's
FLUX realization): discrete adjoints through solvers via the implicit
function theorem — ONE transposed solve sharing the primal's
infrastructure, never a differentiated Krylov iteration; time-dependent
adjoints under revolve checkpointing; Hadamard boundary AND
density/SIMP volumetric shape derivatives; Sobolev (H¹) gradient
smoothing; and the gradient-verification gate ci-gauntlet wires so a
solver without a passing gradient check cannot merge.

## Public types and semantics

- `ift_gradient_matfree(jacobian, ∂J/∂u, ∂J/∂p, (∂R/∂p)ᵀ·, tol, …)` —
  dJ/dp = ∂J/∂p − (∂R/∂p)ᵀλ with Jᵀλ = ∂J/∂u solved by fs-solver
  transposed GMRES; returns (gradient, `AdjointReport` with the
  adjoint residual — the honesty check).
- `DensityPoisson`/`DensityOp` — the density-parameterized Poisson
  family K(ρ) = Σ_t ρ_t·|V_t|·G_t (matrix-free per-cell apply) with
  `density_pullback(λ, u)[t] = λᵀK_t u` — the EXACT volumetric (SIMP)
  chain rule.
- `sobolev_smooth(M, K, α, g_raw, tol)` — the H¹ Riesz step
  (M + αK)·g = M·g_raw through fs-solver CG. The metric α is a real
  tradeoff (measured: α = h² beats stronger smoothing, which clamps
  the signal toward the interior-reduced operator's zero boundary).
- `hadamard::{volume_shape_gradient, compliance_shape_gradient,
  boundary_faces}` — boundary-integral shape derivatives on the
  complex's boundary faces (outward normals derived from the owning
  tet). Compliance sign is PLUS: dJ[V] = +∫(∂u/∂n)²(V·n) dA, pinned
  by the 1D closed form and by the FD-consistency gate that CAUGHT
  the minus-sign first draft.
- `HeatAdjoint`/`heat_initial_gradient` — backward-Euler heat with
  terminal-misfit gradient w.r.t. the initial condition by a
  revolve-checkpointed reverse sweep (O(log N) memory; every reverse
  step is a transposed solve; the forward-recompute count is returned,
  not hidden).
- `verify::verify_gradient(j, p, gradient, directions, eps, tol)` —
  central-FD directional checks returning a `GradientVerdict` (worst
  relative error + per-direction pairs). The gate itself is tested to
  REJECT a corrupted gradient — a gate that cannot fail is not a gate.

- `transpose` module (addendum Proposal 1, bead bk0o.1; [F], behind
  the `ledger-transpose` feature until its Gauntlet tier + kill metric
  are green): TRANSPOSE THE LEDGER. `VjpRegistry` makes the op-spec
  amendment executable — every op registers a VJP or declares itself
  non-differentiable WITH color consequences; `Tape::transpose` chains
  the VJPs across seams in deterministic reverse order (bit-equal
  re-runs). A missing VJP or a declared op inside a differentiation
  path is a STRUCTURED, LOUD block — never a silent zero (the Goodhart
  trap the review named). `check_transpose` is the ⟨Av,w⟩=⟨v,Aᵀw⟩ G0
  battery; `fd_falsifier` is conditioning-aware (the FD self-error at
  two step sizes widens the band so ill-conditioned seams don't fire
  false hits). `CheckpointStore` is the content-addressed spill
  contract shared with Proposal 2's cache discipline: `spilled_adjoint`
  reproduces BIT-EQUAL gradients with and without spill (f64↔bytes
  round-trips are exact), verified against the real fs-ledger CAS via
  dev-deps.

- `mitigate` module (addendum Proposal 1, bead bk0o.2; [F], behind
  `diff-mitigations` → `ledger-transpose`): the three meshing
  mitigations IN ORDER. (1) DIFFERENTIABILITY AS ROUTING: the Rep
  Router's fitness gains a differentiability term at the cost-oracle
  seam (`DiffAwareOracle` penalizes non-differentiable edges under a
  gradient request — no fs-geom changes; the oracle IS the fitness),
  so SDF/spline paths win when gradients are requested while plain
  queries keep the cheap mesh path. (2) HADAMARD boundary forms as the
  mesh-free path (base `hadamard` module, verified against
  perturbation-resolve). (3) UNAVOIDABLE remesh in the path →
  `GradientGrade::EstimatedWithDiscontinuity`: Proposal-3 Estimated
  color with INFINITE dispersion plus the crossing edges named — never
  a silently-verified gradient across a topology event. `grade_ops` is
  the tape-level twin. Deterministic tie-breaks inherited from the
  router.

- `certs` module (addendum Proposal 1, bead bk0o.3; [S], behind
  `gradient-certs` → `diff-mitigations`): GRADIENTS ARE CLAIMS and get
  colors. `adjoint_residual_bound` computes a VERIFIED enclosure of the
  transpose-consistency residual in outward-rounded fs-ivl arithmetic;
  `fd_spot_checks` runs the mandatory falsifier pairing
  (`adjoint-gradient` → `finite-difference-spot-check`, declared in
  fs-evidence's standard registry) along seeded directions with
  conditioning-aware tolerances; `certify` assigns colors —
  smooth+bounded = Verified(residual), flagged remesh = Estimated
  (inherited, never upgraded: that would be laundering), anchored =
  Validated(regime, dataset), evidence-free = Estimated (a gradient
  without a certificate is folklore); `merge_gate` is the CI
  gradient-gate discipline extended across seams: missing or failing FD
  checks refuse with teaching text.

## Invariants

- No differentiation through Krylov iterations anywhere: adjoints are
  IFT transposed solves at converged primals.
- Every gradient this crate produces is FD-verified in its battery
  (the adjoint-vs-FD legs of the acceptance triangle; the dual-number
  leg lives where duals can reach — fs-opdsl/fs-ad).
- Reverse time sweeps share the forward sweep's operators and
  tolerances.

## Error model

Structured panics on dimension mismatches and failed inner solves
(fixture-scale SPD systems — a failure is a modeling bug). Adjoint
solve quality is REPORTED (`AdjointReport.adjoint_residual`), never
assumed.

## Determinism class

Bit-deterministic through the fs-solver/fs-tilelang reduction
discipline; golden FNV-64 over IFT, Sobolev, and heat-adjoint
gradients: `0x0896_7e37_81b3_c044`, recorded on Apple M4 Pro,
verified on Threadripper (x86_64).

## Cancellation behavior

Bounded synchronous computation over fs-solver's iteration-granular
solvers; revolve sweeps checkpoint by construction (pause = keep the
snapshot set). Cx wiring is driver scope (workspace discipline).

## Unsafe boundary

None. `unsafe_code = "deny"`.

## Feature flags

None.

## Conformance tests

`tests/adjoint_battery.rs` (8 cases): IFT source-parameter gradient
vs central FD along random directions (rel ≤ 1e−6, adjoint iterations
reported); the SIMP density chain rule FD-verified (rel ~1e−6);
Sobolev smoothing DEMONSTRABLY rescuing a grid-noise-corrupted
gradient (cosine alignment 0.49 → 0.93, numbers logged; the α
tradeoff measured); Hadamard volume gradient vs perturb-and-resolve
FD to 1e−7 relative (nonzero-divergence velocity ON PURPOSE — a
divergence-free first draft produced two zeros and a meaningless
comparison, kept as a comment); Hadamard compliance CONSISTENCY gate
(gap to FD shrinks under refinement, 0.84 → 0.66 — the P1 one-sided
normal trace is low-order, so the volumetric form is the production
path at lowest order); revolve-checkpointed heat adjoint vs FD
(recompute counts logged); the verification gate accepting a correct
gradient and REJECTING a corrupted one; cross-ISA golden hash.

## No-claim boundaries

- Ledger-backed checkpoint SPILL (degrade-to-storage for
  long-horizon adjoints) is recorded follow-up — revolve's O(log N)
  recomputation tier is what ships here.
- The Hadamard boundary form at P1 is consistency-verified, not
  tight: production shape gradients at lowest order should use the
  volumetric/density form (exactly verified); the boundary form
  earns tight tolerances with high-order traces (fs-feec tfz.6
  spaces — follow-up integration).
- Free-surface LBM adjoints deferred [M] (the bead's honesty note);
  gradient-based fluid shape steps use continuous-adjoint NS when
  that lane lands.
- No second-order adjoints (Hessian-vector products via
  adjoint-of-adjoint — 7tv.3's TR-Newton-Krylov consumes them; the
  fs-ad `second_directional` dual path covers small parameter
  counts).
- The ci-gauntlet WIRING of `verify_gradient` (merge-blocking) is
  that pipeline bead's scope; the gate helper and its self-test live
  here.

## No-claim boundaries (transpose)

- The registry chains VJPs the ops SUPPLY; per-solver adjoint
  correctness (IFT, revolve schedules, Hadamard forms) is the base
  modules' contract, not re-verified here.
- Solver VJPs must be TRANSPOSED solves; nothing in this module
  differentiates through Krylov iterations, and non-symmetric transposed
  preconditioner plumbing (BDDC transposes) lands with the solver-dd
  bead.
- The kill criterion (adjoint-driven optimization beating derivative-free
  baselines on ≥70% of wedge tasks) is a QUARTERLY measurement owned by
  governance; this module ships the machinery and the falsifier.
- Tape recording is caller-driven (the forward code applies ops and
  records); automatic capture from live ledger op rows lands with the
  Proposal-2 integration.

## No-claim boundaries (mitigations)

- The discontinuity FLAG declares the event; at the event itself the
  conditioning-aware FD falsifier deliberately widens (Richardson
  self-error) rather than false-firing — the exploded self-error IS the
  measurable discontinuity signature, and sub-gradient/smoothing
  treatments of topology events are out of scope.
- The non-differentiable edge set is declared by converter owners (the
  registry of record is the VJP registry's declarations); automatic
  detection of topology events inside meshing kernels is fs-mesh's
  future contract.
- Hadamard applicability (volume, boundary compliance) follows the base
  module's fixtures; general goal functionals need their own boundary
  forms.

## No-claim boundaries (certs)

- The interval residual certifies the TRANSPOSE PAIR's consistency,
  not the objective's differentiability — path smoothness is the
  routing grade's claim (mitigate module).
- Anchoring evidence (Proposal 11 assimilation) is caller-supplied;
  this module records it, it does not validate the dataset.
- The CI wiring of merge_gate into the Gauntlet runner is the
  ci-gauntlet bead's contract; this module ships the gate function and
  its refusals.
