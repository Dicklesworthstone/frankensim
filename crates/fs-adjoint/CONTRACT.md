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
  adjoint residual). Every component of the returned vector is the
  TOTAL derivative or the call refuses: `∂J/∂p` and the pullback are
  length-asserted against each other, because combining them with a
  `zip` silently truncated to the shorter and returned a vector of the
  expected length whose tail carried only the explicit partial.
  NO-CLAIM on the adjoint residual: it is measured against the operator
  the caller supplied through `LinearOp::apply_transpose`, so it is an
  honesty check on the SOLVE, not on the transpose. fs-solver does not
  verify symmetry, and a nonsymmetric Jacobian that inherits the
  default `apply_transpose` yields a small residual for the wrong
  system (see the fs-solver CONTRACT's `LinearOp` no-claim).
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
  relative error + per-direction pairs). The finite-difference experiment is
  exactly `p ± eps * direction`; `eps` is a literal scalar step and no hidden
  point/direction-norm rescaling occurs. Relative error uses
  `max(|analytic|, |fd|, 1e-12 * ‖direction‖∞)`, making the comparison floor
  homogeneous under paired rescalings of a direction and the inverse step.
  The complete experiment still refuses rescalings whose perturbations or
  arithmetic are not representable. Non-finite points, gradients,
  directions, objective values, or arithmetic intermediates, and non-finite or
  non-positive steps/tolerances, fail closed with `pass=false`, a deterministic
  positive-infinity error sentinel, and only the finite directional-pair prefix.
  Empty/all-zero directions and finite perturbations that round back to the
  unperturbed coordinate are likewise refused as vacuous evidence, matching the
  certificate path's perturbation preflight. A direction is only counted as
  evidence when it carries SIGNAL — when `max(|analytic|, |fd|)` clears the
  same `1e-12·‖direction‖∞` floor. A probe whose analytic and finite-difference
  derivatives are BOTH at that floor (both exactly zero is the reachable case:
  a forward path that drops the parameter produces a zero cotangent AND two
  bit-identical objective values) scores a relative error of `0.0` — the
  PERFECT score — whether the gradient is right or silently zero, so it is no
  evidence at all. `GradientVerdict::informative_directions` reports the count,
  and a verification with none of them fails closed with the `+∞` sentinel
  instead of passing. The gate itself is tested to REJECT corrupted,
  non-finite, and vacuous evidence — a gate that cannot fail is not a gate.

- `transpose` module (addendum Proposal 1, bead bk0o.1; [F], behind
  the `ledger-transpose` feature until its Gauntlet tier + kill metric
  are green): TRANSPOSE THE LEDGER. `VjpRegistry` makes the op-spec
  amendment executable — every op registers a VJP or declares itself
  non-differentiable WITH color consequences; `Tape::transpose` chains
  the VJPs across seams in deterministic reverse order (bit-equal
  re-runs). A missing VJP or a declared op inside a differentiation
  path is a STRUCTURED, LOUD block — never a silent zero (the Goodhart
  trap the review named). VJP SHAPES ARE CHECKED, not assumed: arity
  AND per-input cotangent length are asserted against the recorded
  primals on every sweep (as is the seed against the output node), and
  a violation panics — accumulation used to `zip`, which truncates to
  the shorter and made a masked cotangent's dropped tail
  indistinguishable from a genuine zero. `"leaf"` is a RESERVED op
  name that terminates the sweep without a registry lookup, so
  `Tape::apply` rejects it rather than letting an unregistered op end
  the chain silently. `check_transpose` is the ⟨Av,w⟩=⟨v,Aᵀw⟩ G0
  battery: it REFUSES zero probes and empty dimensions (both of which
  returned the perfect score `0.0` from an experiment that never ran)
  and returns the probe count and pairing scale alongside the absolute
  residual so a caller's threshold can be relative. `fd_falsifier` is
  conditioning-aware (the FD self-error at two step sizes widens the
  band so ill-conditioned seams don't fire false hits) and its band is
  HOMOGENEOUS in the derivative scale: `base_tol·max(|adjoint_dd|,
  |fd_fine|, noise_floor) + 3·self_error`, where `noise_floor` is the
  difference quotient's own rounding floor. The old absolute
  `max(…, 1.0)` term made the band never tighter than `base_tol` in
  absolute units, so every derivative below `base_tol` accepted `0` and
  the exact sign flip. When both derivatives sit at or below the noise
  floor the verdict is `falsifiable = false` — "no signal at this
  step", NOT agreement — and `consistent` is false, which the merge
  gate refuses with that reason. `CheckpointStore` is the content-addressed spill
  contract shared with Proposal 2's cache discipline: `spilled_adjoint`
  reproduces BIT-EQUAL gradients with and without spill (f64↔bytes
  round-trips are exact), verified against the real fs-ledger CAS via
  dev-deps.

- `mitigate` module (addendum Proposal 1, bead bk0o.2; [F], behind
  `diff-mitigations` → `ledger-transpose`): the three meshing
  mitigations IN ORDER. (1) DIFFERENTIABILITY AS ROUTING: the Rep
  Router first plans over a smooth-only subgraph under the request's ORIGINAL
  cost/error budgets. Only when no differentiable route is admissible does it
  plan over the complete graph and downgrade a winning remesh path. This is a
  deterministic lexicographic policy, not a synthetic wall-cost penalty, so an
  expensive but admissible smooth path still wins and an unavoidable affordable
  remesh is not made infeasible by policy arithmetic. Plain queries keep the
  ordinary cost-optimal mesh path. (2) HADAMARD boundary forms as the
  mesh-free path (base `hadamard` module, verified against
  perturbation-resolve). (3) UNAVOIDABLE remesh in the path →
  `GradientGrade::EstimatedWithDiscontinuity`: Proposal-3 Estimated
  color with INFINITE dispersion plus the crossing edges named — never
  a silently-verified gradient across a topology event. `grade_ops` is
  the tape-level twin. Deterministic tie-breaks inherited from the
  router.

- `certs` module (addendum Proposal 1, bead bk0o.3; [S], behind
  `gradient-certs` → `diff-mitigations`): GRADIENTS ARE CLAIMS and get
  colors. `adjoint_residual_bound` computes an outward-rounded diagnostic for
  the exact seeded transpose-consistency probes evaluated; it validates sparse
  dimensions, storage, indices, finite weights, generated vector components,
  and the aggregate forward-plus-transpose sparse-entry visits before probe
  allocation or indexed execution, and returns structured errors on refusal.
  `fd_spot_checks` runs the mandatory falsifier pairing
  (`adjoint-gradient` → `finite-difference-spot-check`, declared in
  fs-evidence's standard registry) along seeded directions with
  conditioning-aware tolerances and bounded vector/work budgets, returning a
  sealed `FdCheckBatch`. Its domain-separated context digest binds a
  structurally valid leaf objective identity, exact point/gradient bits, seed,
  direction count, and every verdict bit. Before evaluating the objective, the
  exact coarse/fine plus/minus perturbations are checked: each coordinate must
  remain finite and every nonzero direction component must change the stored
  coordinate at both signs and both step sizes. Fixed steps that round away at
  large coordinates and an all-zero sampled direction therefore refuse with a
  location-bearing error rather than yielding a meaningless zero-difference
  pass. `GradientCertificate` fields are private and read-only through
  accessors; `certify` consumes an optional sealed batch and
  retains valid residual/caller-anchor metadata as diagnostics but emits only
  `Estimated`: sampled transpose consistency does not prove the gradient
  formula, and a raw public dataset/regime cannot authenticate itself. Flagged,
  contradicted, malformed, and evidence-free gradients also remain Estimated;
  `merge_gate(cert, expected_context_digest)` is the CI
  gradient-gate discipline extended across seams: missing or failing FD
  checks refuse with teaching text.

- `dwr_accept` module (addendum Proposal 9, bead lmp4.4; [F], behind
  `dwr-accept` → optional fs-verify dep): the GOAL-ORIENTED accept
  test. `dwr_integral_qoi` is the 1-D reference DWR (enriched dual on
  the once-refined mesh, per-element fine-remainder indicators) for the
  explicit functional
  `J(u) = ∫_{[w_lo, w_hi] ∩ [0, 1]} u(x) dx` and returns
  `Result<DwrOutput, DwrError>`. The problem arrives as fs-verify's immutable,
  fallibly admitted `MmsProblem`: its canonical class identity binds the exact
  solution and derived forcing, its canonical domain is exactly `[0,1]`, and
  DWR shares fs-verify's 1..=6 coefficient envelope rather than advertising a
  divergent class. Finite strictly increasing raw windows may extend outside
  that domain: partial overlap clips numerically to the intersection and a
  disjoint window denotes the exact zero functional. The raw `w_lo`/`w_hi`
  binary64 bits are nevertheless retained, identity-bound, and compared
  bitwise, so equal clipped intersections (including two empty intersections
  or `-0.0` versus `+0.0`) are not interchangeable replay claims.
  `DwrQuery` schema 2 is the first unambiguous schema to name this meaning as
  `DomainClippedWindowedIntegral`; schema 1's ambiguous windowed-integral
  meaning is superseded and is not accepted or reinterpreted. Schema 2 binds
  homogeneous canonical `+0.0` Dirichlet endpoints and carries six dimension
  exponents in fs-qty order.
  The reference scalar API seals `[0; 6]`; callers cannot relabel its output as
  dimensional. DWR admission validates
  2..=1,000,000 coarse mesh nodes through fs-verify's shared node cap,
  candidate shape/finite values and bit-canonical homogeneous `+0.0`
  endpoints, strict finite cell geometry, representable midpoints, a finite
  increasing raw QoI window, and a conservative mesh×polynomial aggregate
  budget of 100,000,000 work units before refined allocation. The
  discontinuous QoI indicator is never sampled over a whole cell: primal P1
  values and dual P1 loads use the closed-form formula on each
  cell/window/domain intersection, so narrow non-mesh-aligned windows cannot
  alias to zero. The formula is evaluated in binary64, not exact-expansion or
  interval arithmetic. A computed-zero nonempty overlap is accepted only when
  the original P1 data prove an identically-zero field, or prove exact
  cancellation from additive-inverse cell endpoints over a clipped interval
  exactly symmetric about that cell's midpoint; every other such zero is
  refused until
  exact-expansion integration can distinguish true cancellation from binary64
  underflow/rounding. Empty domain intersections are exactly zero by the
  functional definition. Assembly, residual quadrature, forcing,
  elimination pivots, slopes, residuals, and outputs refuse on any non-finite
  derived value; zero-interior-DOF dual systems are handled explicitly rather
  than indexed through an empty solver. `DwrOutput` construction is sealed;
  it retains its complete typed functional and v7 evidence identity, and
  `accept` cross-checks the query schema, raw window bits, endpoint convention,
  semantics, and dimensions against that retained binding before considering
  the tolerance. `accept` encodes the color logic MECHANICALLY: all current
  accepts carry ESTIMATED. `AcceptOutcome::refused` is true only when malformed
  public inputs prevented an accept/reject decision; it is false for both a
  valid acceptance and a valid over-tolerance rejection. `Bracket` fields
  are sealed and `Bracket::cauchy_schwarz` reruns the equilibrated verifier on
  bounded exact problem/candidate inputs, but the resulting energy-product is
  diagnostic only: the typed functional query does not yet carry a
  re-verifiable proof that the independently supplied dual problem is the
  exact dual of that functional. Consequently a bracket can neither promote
  nor veto acceptance and a public/forged `VerifierReport` is never consumed
  as authority. Non-finite/non-positive tolerances and non-finite/negative DWR
  estimates refuse with structurally valid infinite-dispersion colors. Machine
  estimator ids remain separate from human QoI/audit labels. In the solved
  G1 fixture, refinement indicators concentrate where the QoI error lives;
  the broader algebraic-indicator limitation is stated below.

- `explain` module (addendum Proposal B, bead knh1.5; [F], behind
  `explanation-objects`): explanation OBJECTS — a tree of
  `(channel, contribution, bound, color, evidence, derivation_digest,
  payload_fingerprint, fingerprint_version)` nodes plus a versioned top-level
  receipt. Node claim fields are private and read-only through accessors; a
  private origin tag distinguishes built-in derivations from unretained caller
  payloads and is fingerprint-bound. Built-in derivation digests bind exact
  floating-point input bits, operator/domain versions, masks or edit position,
  and full history/problem roots. Three engines:
  `adjoint_attribution` (the exact discrete bilinear identity
  `J₁−J₀ = −∫Δa·u₀′·u₁′` over channel masks, with its unproved floating-point
  solve/accumulation allowance honestly `Estimated`),
  `provenance_attribution` (input-bound caller tuples with one-ulp subtraction
  envelopes, always unretained `Estimated` until a ledger replay receipt is
  authenticated),
  and the far-field drag flagship (private-construction elliptic `LiftingLine`
  Trefftz wake integral + viscous strip + an explicit wave-model declaration;
  all three are `Estimated`. The measured O(1/n) midpoint trend is a
  conformance diagnostic, not an outward-rounded discretization proof).
  Public construction/execution is fallible: `ExplanationNode::new`,
  `finalize`, `Elliptic1d::{new,solve,compliance}`, `adjoint_attribution`,
  `provenance_attribution`, `LiftingLine::{elliptic,cl,
  induced_drag_coefficient,aspect_ratio}`, and `drag_decomposition` return
  `Result<_, ExplanationError>`. `Elliptic1d` dimensions and all `LiftingLine`
  fields are sealed; caps, shapes, indices, identities, colors, and finite
  inputs are checked before allocation/indexing, while assembly, pivots,
  accumulation, subtraction envelopes, lift/drag normalization, and receipt
  arithmetic fail closed on non-finite derived values.
  `finalize` is THE HONESTY GATE: an unattributed residual above its
  effective limit produces `Refused` (the partial tree is forensics, not a
  claim). Only built-in-origin `Verified` node bounds discharge certified
  residual coverage. Unretained caller nodes never do and self-declared
  `Verified` caller colors demote to unbounded `Estimated` in the aggregate;
  `Estimated` bounds never certify coverage, and `Validated` nodes fail closed
  until the module has a retained regime-membership witness. Duplicate node
  fingerprints or derivation digests are rejected, preventing clone-based
  multiplication of trusted coverage. Built-in nodes also carry a private,
  fingerprint-bound batch digest/index/size; a receipt must contain exactly one
  complete built-in attribution batch, so trusted nodes from separate engine
  calls cannot be recombined to inflate coverage. Contribution summation uses
  an explicit outward-rounded enclosure (zero aggregation roundoff for one
  node), and the effective limit is the stricter of the requested threshold and
  certified coverage. Node colors compose under Add semantics, so the receipt's
  aggregate color cannot outrank its weakest term. The receipt root binds the
  outcome variant, ordered node/derivation digests, observation, residual,
  requested threshold, certified coverage, effective limit, aggregation
  roundoff, aggregate color, and color-algebra version. `reconciles` is true
  only for a structurally valid `Explained`; `is_structurally_valid` separately
  validates honest refusals. `render_narrative` remains explicitly
  NON-AUTHORITATIVE. Builders reject non-finite arithmetic, negative bounds or
  thresholds, blank/control-bearing or oversized identities, malformed color
  envelopes, digest/version mismatch, and malformed fixture dimensions.
  Version-2 node fingerprints use exact length-framed binary fields,
  `Color::canonical_bytes`, and the in-tree domain-separated BLAKE3 owner.
  The owner-local `fs-adjoint:explanation-node` declaration classifies every
  `ExplanationNode` and `ExplanationNodeAuthority` field, names the exact
  domain/version/encoder/consumers, and links an independent mutation for each
  semantic input. Evidence-list coverage distinguishes count, order, and item
  mutations rather than treating the vector as one opaque case; the declared
  artifact domain is independently mutated and a digest recomputed under a
  foreign domain is refused by the live verifier. The companion
  `fs-adjoint:explanation-receipt` declaration exhaustively classifies the
  receipt fields and binds the outcome variant, ordered node count/order and
  each retained node-identity component, observed value, residual, receipt
  policy arithmetic, color-algebra version, and exact aggregate-color bytes.
  Its guards refuse unknown receipt/color-algebra versions and foreign-domain
  roots even when the supplied root is internally recomputed. `xtask
  check-identities` compares both declarations
  to the generated registry and their golden-coupling surfaces. Fingerprint and
  receipt-root outputs are derived and recomputed; natural-language rendering
  and `Color::payload_json()` are presentation-only and forbidden as identity
  inputs. The unretained, adjoint, provenance, lifting-line, and drag
  derivation grammars remain explicit follow-up registry rows; shared v1
  derivation domains that currently host more than one grammar are not claimed
  collision-separated by a role tag until a deliberate versioned rotation
  lands.
  Finalized trees contain 1..=1024 unique nodes. Adjoint channel names are
  unique and masks are nonempty, duplicate-free, and mutually disjoint (partial
  coverage is allowed so the honesty gate can expose omitted channels).
  Provenance histories contain 1..=1024 uniquely named finite edits and must
  telescope exactly between adjacent states (bitwise, with signed zero
  canonicalized). Elliptic explanation fixtures contain 1..=65536 interior
  nodes; elliptic lifting-line fixtures contain 1..=4096 stations, bounding the
  fixture solver allocations and quadratic wake summation.

## Invariants

- No differentiation through Krylov iterations anywhere: adjoints are
  IFT transposed solves at converged primals.
- Every gradient this crate produces is FD-verified in its battery
  (the adjoint-vs-FD legs of the acceptance triangle; the dual-number
  leg lives where duals can reach — fs-opdsl/fs-ad).
- Reverse time sweeps share the forward sweep's operators and
  tolerances.

## Error model

DWR-accept budget contract (bead sj31i.6): `Bracket::cauchy_schwarz`,
`accept`, and `dwr_integral_qoi` take a caller-supplied deterministic
`fs_exec::TimeSource` and admit the ambient `cx.budget()` plus the
preflighted work plan through `fs_exec::AdmittedBudget` before any work
(one logical work unit = one cost unit): expired deadlines
(`Budget::ZERO` included) and over-quota cost plans refuse at admission,
every checkpoint enforces cancellation → deadline → poll quota, and each
advanced work unit charges cost. Poll/deadline/cost exhaustion is a
typed `DwrError::BudgetRefused`/`BracketError::BudgetRefused` retaining
the accountant's refusal verbatim; real cancellation keeps the
structured `Cancelled` shape with exact completed/planned work.

The hardened public certificate, DWR, and explanation paths return structured
`GradientCertError`, `DwrError`, and `ExplanationError` refusals for malformed
inputs, resource caps, invalid indices, unrepresentable perturbations,
non-finite derived arithmetic, and failed fixture pivots. They do not use panic
as caller-input validation and do not return partially authoritative objects on
error. Older base adjoint fixture modules may still treat impossible internal
dimension/solver failures as assertions; their solve quality is REPORTED
(`AdjointReport.adjoint_residual`), never assumed.

## Determinism class

Bit-deterministic through the fs-solver/fs-tilelang reduction
discipline; golden FNV-64 over IFT, Sobolev, and heat-adjoint
gradients: `0x0896_7e37_81b3_c044`, recorded on Apple M4 Pro,
verified on Threadripper (x86_64).

## Cancellation behavior

Bounded synchronous computation over fs-solver's iteration-granular
solvers; revolve sweeps checkpoint by construction (pause = keep the
snapshot set). The public DWR, bracket, and accept paths own explicit `Cx`
wiring, checked complete work plans, caller poll quotas, and phase/publication
checkpoints. Their variable scans poll at an identity-bound 256-item stride;
independent variable-length identity streams are separated by a checkpoint so
the stride cannot be reissued across a seam. Cancellation returns typed phase
and exact completed/planned work and publishes no partial authoritative
object. The retained DWR execution identity binds work-plan version 2,
poll-policy version 3, evidence identity version 7, typed-functional schema
version 2, the canonical domain-clipped-windowed-integral semantics tag, and
the unmodified window endpoint bits.

## Unsafe boundary

None. `unsafe_code = "deny"`.

## Feature flags

All OFF by default per the Ambition-Tag rule (huq.18 reconciliation —
this section previously said "None" while the manifest declared five):
- `ledger-transpose` [F] — the ledger-DAG transposition layer (VJP
  registry); disabled until its Gauntlet tier + kill metric are green.
- `explanation-objects` [F] — explanation objects (knh1.5, Proposal B).
- `diff-mitigations` [F] — non-differentiable meshing mitigations
  (smooth-first routing, Hadamard path, estimated+flag downgrade); implies
  `ledger-transpose`.
- `gradient-certs` [S] — gradient certificates (colors + interval
  residual bounds + the FD-falsifier merge gate); implies
  `diff-mitigations`.
- `dwr-accept` [F] — the DWR goal-oriented accept test.
Each gates its own integration target (required-features declared).

## Conformance tests

`tests/adjoint_battery.rs`: IFT source-parameter gradient
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

`tests/mitigate.rs` (8 cases, feature `diff-mitigations`): ordinary routing
takes the cheap mesh path while gradient routing selects an admissible smooth
path; unavoidable remeshing remains affordable under the ORIGINAL budget and
is downgraded to Estimated with an explicit discontinuity; Hadamard and direct
path checks retain their falsifiers; tie-breaking is deterministic; an
admissible smooth route wins even above the former fixed-penalty scale; an
over-budget smooth route falls back honestly; and an empty graph preserves the
router's structured no-path refusal. Smooth-first and fallback planning consume
one spec-scoped oracle snapshot, so a live cost backend cannot change authority
between the two policy passes.

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

- The interval arithmetic encloses each sampled transpose residual, not the
  unsampled operator domain and not the objective's gradient error. It cannot
  mint `Verified`; path smoothness is separately only the routing grade's
  claim (mitigate module).
- Anchoring metadata (Proposal 11 assimilation) is caller-supplied and retained
  only after bounded structural validation. It cannot mint `Validated` until a
  sealed, independently checkable anchored-source certificate exists.
- Probe count, direction count, dimensions, sparse entries, total generated
  scalar components, and aggregate sparse-entry visits are explicitly capped.
  One transpose probe accounts for two visits per stored entry (forward and
  transpose apply), with at most 16,777,216 visits per call, so individually
  admissible entry/probe counts cannot amplify into an unbounded scan.
  Malformed indices/non-finite values return location-bearing
  `GradientCertError`s before allocation or indexing.
- The sparse-entry-visit limit is a conservative deterministic work proxy, not
  a wall-time, energy, memory-bandwidth, or cancellation-latency certificate.
  It bounds this synchronous v0 diagnostic while Cx-aware tiled probing remains
  driver/integration work.
- FD probes use fixed absolute coarse/fine steps. When either signed
  perturbation is non-finite or rounds back to the original coordinate, the
  batch refuses; adaptive representable step selection is future work, not a
  silent fallback.
- The batch prevents raw-verdict fabrication and makes cross-context replay
  observable/refusable when the orchestrator supplies its expected digest. The
  objective identity itself is caller-declared, not authenticated: a malicious
  caller can evaluate a different closure under the same name. External
  authority still requires a typed objective receipt whose owner retains and
  checks the expected context digest.
- The CI wiring of merge_gate into the Gauntlet runner is the
  ci-gauntlet bead's contract; this module ships the gate function and
  its refusals.

## No-claim boundaries (dwr-accept)

- The reference estimator is the current 1-D manufactured elliptic class;
  the 2-D cutfem DWR lives in fs-dwr (tfz.23) and plugs into
  the same accept/color logic.
- The Cauchy–Schwarz bracket is sharp only up to the product's
  pessimism; sharper goal-oriented equilibrated bounds are the
  verifier's growth path, not this module's claim.
- No current DWR path emits `Verified`. Schema 1 is superseded and refused;
  schema 2 is the first query schema that binds the exact domain-clipped
  functional requested. Promotion additionally requires a
  re-verifiable relation proving that the supplied dual problem is that
  functional's mathematical dual—not two reports (even two genuine reports
  for unrelated problems) and not a caller-authored guarantee bit.
- The retained v7 identity provides exact local payload/execution identity; it
  is not external authority. `accept` does not yet take separately trusted
  expected problem, candidate, or execution-policy identities and reexecute
  `dwr_integral_qoi` from those inputs before deciding.
- The primal-equilibrium acceptance gate currently compares the maximum
  interior assembled residual to an absolute `1e-9` threshold. This is not yet
  a scale-safe componentwise/backward-error theorem. The global `eta` includes
  the coarse algebraic residual contribution, but the per-cell refinement
  indicators currently expose only the fine-remainder contributions and do
  not localize that algebraic term.
- `DwrError` is an execution refusal, not an Estimated answer. The aggregate
  work cap bounds the current scalar-work model; it is not a wall-time or
  memory certificate, and the reference estimator remains limited to the 1-D
  manufactured elliptic class.
- Falsifier cadence (how often the high-fidelity spot check runs) is
  the budget allocator's decision (Proposal 6); this module ships the
  check and the pairing.

## No-claim boundaries (explain)

- The adjoint engine's exactness is the compliance/self-adjoint case;
  general QoIs get first-order attribution with remainder bounds — the
  growth path.
- The lifting-line flagship is the incompressible far-field fixture
  (wake integral vs the analytic elliptic envelope). Its state is private and
  its O(1/n) allowance is still only a measured/heuristic diagnostic; it cannot
  create a built-in `Verified` channel or certified residual coverage.
  Full-CFD far-field decomposition rides
  fs-bem/fs-vpm when their wake machinery lands.
- The viscous strip and wave-model channels are `Estimated` by construction, so
  the aggregate drag explanation is also `Estimated`. A Mach value plus an
  evidence-label string never certifies exact zero physical wave drag.
- Raw provenance tuples authenticate neither ledger membership nor replay and
  never create built-in/Verified authority. The v0 adjoint rounding allowance
  is likewise a heuristic Estimated diagnostic; Verified attribution requires
  outward-rounded solve and accumulation evidence.
- Provenance one-ulp diagnostics require both adjacent representable values and
  both gaps to be finite and positive. Finite extrema such as `f64::MAX` refuse
  instead of silently dropping the infinite side and reporting a false
  symmetric dispersion.
- `ExplanationError` reports invalid public inputs or non-finite derived
  arithmetic; it carries no partial authority object. Successful `Result`
  construction still provides payload integrity rather than external evidence
  authentication.
- BLAKE3 roots provide exact payload integrity and deterministic replay
  identity, not external authority. This crate does not retain expected roots,
  resolve evidence links, validate caller datasets, or sign receipts. An
  external checker must replay the content-addressed inputs and compare against
  a separately retained root before treating a digest as authoritative.
- Caller-created `ExplanationNode::new` values carry an explicitly unretained
  payload-derived derivation digest. They are locally tamper-evident after
  finalization but do not prove that the named evidence exists, cannot
  contribute certified residual coverage, and cannot produce a `Verified`
  aggregate claim merely by self-labeling their color.
- Downwash sign convention (Katz & Plotkin) is load-bearing and was
  caught by the analytic envelope during development — the conformance
  test is the regression guard.
