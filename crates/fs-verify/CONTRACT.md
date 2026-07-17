# fs-verify — CONTRACT

## Purpose and layer

L3 (FLUX-adjacent). The certified-speculation VERIFIER (Proposal 9,
the addendum's SINGLE research bet): a certified accept test that lets
untrusted proposers be maximally aggressive, because correctness never
depends on the learned component. For the elliptic class,
equilibrated-flux a-posteriori estimators (Prager–Synge) give
GUARANTEED, CONSTANT-FREE upper bounds on the energy-norm error; the
remaining risk — floating-point rounding — is retired by
OUTWARD-ROUNDED interval evaluation over mathematically exact
quadrature.

## Public types and semantics

- `interval::Iv` — outward-rounded interval arithmetic: every op
  widens one ulp each direction (bit-nudge directed rounding), so
  enclosures are RIGOROUS, not to-nearest-plus-slack. Small on
  purpose (`add/sub/mul/div_pos/sq/sqrt/scale_pos`). Invalid or
  zero-containing positive divisors, wholly negative square-root domains,
  reversed intervals, and nonpositive/non-finite positive scales produce the
  entire real line and therefore force verifier refusal.
- `fem1d` — the v0 class testbed: P1 elements for `−u″ = f` on (0,1),
  polynomial manufactured solutions of degree ≤ 5 (so 5-point Gauss
  quadrature is MATHEMATICALLY exact for every integrand the verifier
  meets — rigor rests on rounding alone), fallible Thomas solves, a
  fallible high-resolution oracle, and the toy nonlinear class
  (`−u″ + u³ = f`, Newton). `Fem1dError` is the shared structured
  validity boundary. `Poly::new` normalizes signed and trailing zeroes;
  `MmsClass::new` admits the homogeneous exact solution and exclusively
  derives `f = -u''` plus the rounded zero-constant antiderivative;
  `MmsProblem::from_class` admits and owns the mesh. All semantic fields are
  private and immutable. `MmsProblem::new` fully admits/canonicalizes its mesh
  before constructing the class, and `with_mesh` does so before fallibly
  copying the retained class data and rebuilding its identity.
  Construction, followed by defensive operation-level
  validation, checks canonical domain/BCs and derived polynomials, finite
  strictly increasing mesh, exact nodal shape,
  finite values, and bounded resources. v0 caps are 1,000,000 mesh
  nodes, six polynomial coefficients (exact-solution degree at most
  five) after canonical zero trimming, 4,096 raw coefficients before
  normalization, 10,000 Newton updates, 4,096 class-name bytes, 4,441
  class canonical-identity bytes, 8,004,626 meshed-problem
  canonical-identity bytes, and 50,000,000 conservative scalar-work units per
  synchronous call. Solver, refinement, ordering, candidate work vectors, and
  canonical identity construction use fallible reserve. The degree-five envelope is
  load-bearing: the equilibrated integrand squares an antiderivative
  of degree `deg(u)-1`, so five-point Gauss exactness requires
  `deg(u) <= 5` even though P1 load assembly alone admits higher degree.
- `fem1d::solve_p1` and `fem1d::true_energy_error` return `Result`;
  numerical failure cannot hide in a NaN/empty ordinary value. The
  correctly rounded ordinary-f64 Gauss constants are independently
  locked between adjacent-f64 truth brackets. `solve_nonlinear`
  returns `NonlinearSolveReport { solution, iterations,
  residual_norm, converged }`, checks finite assembly/update/pivots,
  re-evaluates the residual after the last admitted update, and keeps
  finite nonconvergence distinct from success. `max_iter=0` is an
  initial-residual probe, not a zero-cost success.
- `estimator::verify(problem, candidate, tolerance)` — THE VERIFIER:
  the equilibrated flux is `σ = c − F` (any `c` is sound because
  `σ′ = −f` exactly; the FREE CONSTANT is optimized in plain f64 for
  TIGHTNESS — rigor from structure, tightness from optimization). The
  certified path intervalizes the candidate difference and division,
  mesh subtraction and affine quadrature map, quadrature weights, and
  every coefficient division in the exact antiderivative of the
  authoritative forcing. Gauss nodes and weights use correctly rounded
  full-precision literals, then widen one ulp; independent adjacent-f64
  truth brackets lock that the irrational constants are enclosed. The
  rounded `MmsClass::rounded_forcing_antiderivative` value is replay metadata
  and a tightness aid, never silently treated as a point enclosure of
  exact antiderivative coefficients. The bound `‖σ − u_h′‖` is then
  interval-evaluated. Accept ⟺ `bound.hi ≤ tolerance`; an accept
  carries `Color::Verified`; a reject carries no color.
- `estimator::verify_with_receipt` emits an immutable private-field production
  receipt only after successful verifier publication. Canonical retained bytes
  decode to the distinct `PresentedVerifierReceipt`, whose private production
  record exposes routing identity only: it has no bounds, acceptance bit,
  statement, color, or raw-`VerifierReceipt` accessor. Only
  `admit_verifier_receipt`'s exact independent replay returns
  `AdmittedVerifierReceipt`, whose replay-authenticated result getters and
  `color()` yield `Verified` authority for an accepted receipt and `None` for
  an ordinary above-tolerance reject. A direct in-memory consumer of the
  freshly returned production result may pass its finite enclosure through
  fs-evidence's guarded constructor without replay; that producer trust does
  not survive retention.
- `VerifierWorkPlan::for_inputs` is the constant-time public shape preflight.
  For mesh nodes `n`, cells `c=n-1`, and canonical coefficient counts `u`, `f`,
  and `a`, it charges validation `3+2n+c+3u+3f+2a`, tightness `c`,
  equilibrated construction `c`, flux hashing `7+2n+f+a`, finalization `1`,
  and their checked `u128` sum. Hashing binds the hash domain, flux-identity
  version, optimized `c_star`, both mesh and candidate lengths and values, and
  both polynomial lengths and coefficients; this is distinct work from validation.
  Public accessors and `identity_fields()` expose those six counts in stable
  phase order. Plan version is `VERIFIER_WORK_PLAN_VERSION = 2`.
- `estimator::verify_with_checkpoint` is the generic cancellation/composition
  seam. Its callback receives immutable invocation-global `VerifierProgress`
  at every phase entry, each nonzero global multiple of 256 logical work units,
  a structured-refusal flush, and the final publication gate. Poll schedule
  version is `VERIFIER_POLL_POLICY_VERSION = 1` and stride is
  `VERIFIER_POLL_STRIDE_WORK_UNITS = 256`. A callback error is returned
  unchanged, wins over a pending scientific refusal, and prevents report
  publication. Shape refusal happens before the driver exists and invokes no
  callback. `verify` is the bitwise-equivalent infallible no-op wrapper.
- Public admission precedes indexing or arithmetic. Meshes contain
  `2..=1_000_000` finite, strictly increasing nodes with bit-canonical
  `+0.0` and `1.0` endpoints. Candidates have exactly the mesh length,
  finite values, and bit-canonical homogeneous `+0.0` endpoints.
  Tolerances are finite and positive. Each polynomial input is nonempty,
  finite, bounded to 4,096 raw coefficients, and has at most six canonical
  coefficients (degree at most five) after signed/trailing-zero normalization.
  Manufactured `u(0)` is canonical `+0.0`; `u(1)` is checked by a bounded
  exact binary superaccumulator over the stored finite f64 coefficients.
  Point Horner can hide a nonzero residue, while interval containment only
  proves that zero is possible, so neither is boundary authority. Exact
  binary-rational cancellations remain admissible even when ordinary Horner
  rounds to a nonzero value. Canonical `f = -u''` and its rounded
  zero-constant antiderivative are computed from `u` at the sole class
  construction boundary and exposed only by borrowed accessors. A versioned
  `fs-obs::ReplayIdentity` binds the class name, exact solution, both derived
  polynomials, and their explicit schema; the problem identity additionally
  binds the generic child frame/version/root, the complete class canonical
  bytes, and the exact canonical mesh bytes. Producer schema v2 deliberately
  rotates the formerly unversioned artifact kinds to
  `fs-verify/fem1d-mms-class.v2` and
  `fs-verify/fem1d-mms-problem.v2`; this is distinct from fs-obs's enclosing
  replay-frame v1. `MmsClassIdentityReceipt` and
  `MmsProblemIdentityReceipt` retain producer version, exact bytes, and root.
  Receipt capture reserves fallibly. Retained receipt parts are untrusted
  transport inputs: construction checks the exact producer version before
  inspecting or hashing bytes, then enforces the byte cap and self-consistent
  FNV root. Admission against a live object defensively repeats those checks
  before exact byte/root equality. Stale/future v1/v3 receipts fail closed
  rather than being guessed. Canonically equivalent
  signed/trailing-zero inputs therefore share identity, while changing any
  admitted semantic field changes the corresponding identity. The bounded
  identity builder checks u64 framing and the complete producer-specific byte
  cap, reserves each field before mutation, and maps refusal into the existing
  `Fem1dError::ResourceLimit` or `AllocationFailed` boundary; no partial replay
  identity escapes.
- `VerifierReport::refusal` carries a closed structured reason. A
  refusal preserves the requested tolerance and estimator family but
  returns `[-∞,+∞]`, `accept=false`, no color, and flux hash zero. A
  valid finite reject instead has `refusal=None`. The flux hash binds
  the selected constant and both recomputed polynomials. Ledger rows
  JSON-escape problem names, encode non-finite numbers as `null`, and
  distinguish `refused` from a valid `reject`.
- `estimator::hierarchical_estimate` — the INDEPENDENT second family
  (refined-mesh comparison; not guaranteed; the falsifier's
  cross-check, never a color source). It is fallible and bounds the
  doubled-mesh allocation before refinement.
- `estimator::warm_start` — the honest nonlinear fallback: candidates
  are WARM STARTS with measured iteration savings and an ESTIMATED
  color, never certificates (the R1 boundary). Both cold and warm
  runs must explicitly converge; otherwise the call returns an error
  and no savings record exists. `effectivity` likewise propagates
  oracle/report failure instead of mapping NaN to the ideal value 1.

- `zoo` (bead lmp4.2) — the PROPOSER ZOO behind one `propose()`
  trait: hot-swap `Registry` (register/deregister without touching
  consumers), `speculate()` ordering candidates by ADVISORY confidence
  (descending, NaN last, deterministic name tie-break — confidence
  never enters any accept decision or certificate), and the
  TYPE-LEVEL safety invariant: `CertifiedAnswer` has no public
  constructor — it exists only when the verifier says yes. v0
  proposers: NEIGHBOR EXTRAPOLATION (nearest certified run, Taylor
  correction with a cached sensitivity, zeroth-order degradation,
  smaller-θ equidistant tie-break) and COARSE-RUNG PROLONGATION
  (halved-mesh solve, linear prolongation, honest decline on tiny
  meshes). `quantize_f16` demonstrates the precision discipline:
  speculate LOW, verify HIGH — the certificate inherits the
  VERIFIER's precision. Proposal construction and `speculate()` are
  fallible: malformed queries, poisoned built-in cache coordinates or
  vectors, allocation failure, and coarse-solve failure propagate
  before an ordinary miss can be reported. The registry, neighbor
  cache, and telemetry each admit at most 4,096 entries/keys; counter
  overflow is atomic and structured. `ZooTelemetry` tracks
  per-proposer per-regime accept rates with the AUTO-DEMOTION hook
  (collapse ⇒ disabled in that regime) and ledger rows. First-pass
  rejected candidates are retained behind a sealed outcome: neither
  all-rejected warm starts nor accepted-run drift accounting invokes a
  stateful proposer twice or forgets an earlier rejection.

- `economics` (bead lmp4.3; COMPLEMENTARY to the standalone
  fs-spececo policy crate, which owns the abstract decide/telemetry/
  drift logic — this module is the INTEGRATED loop driving real
  solves through the zoo, and fs-ledger v3 supplies the persistence
  fs-spececo's no-claim defers; drift-logic consolidation between the
  two is a follow-up) — the accept/reject CONTROL LOOP:
  `run_speculative` accepts OUTRIGHT on a certified pass (no solve at
  all); otherwise the best rejected candidate (smallest certified
  bound) WARM-STARTS the true solve with savings MEASURED and recorded
  CLAMPED at ≥ 0 (a worse-than-cold start is never a win; the raw
  negative delta stays in the ledger). `DriftGuard` is the drift
  detector that falls out of the telemetry: an accept-rate collapse in
  a regime demotes the proposer THERE (localized), only after a
  minimum sample count, with probation hysteresis whose evidence bar
  DOUBLES per failed probation (no flapping). Zero-telemetry regimes
  report the conservative 0.0 prior. `solve_node_record` emits the
  four-field schema amendment `(proposer_id, accepted, bound,
  iterations_saved)` stored in fs-ledger's v3 `speculation` extension
  table. Dashboards are kernel × regime × proposer with median
  savings. Drift state admits at most 4,096 proposer/regime keys and
  retains only the most recent 1,024 savings samples per key; recovery
  clears stale savings. Policy scalars and evidence counts are
  validated. The control loop itself returns `Result`; iteration caps
  are checked before proposal work, and solver nonconvergence or
  numerical refusal produces no cold/warm decision, savings telemetry,
  or drift observation for the failed solve.

## Invariants

1. THE UPPER-BOUND PROPERTY (G1 MMS): the bound dominates the oracle
   truth on every battery case INCLUDING adversarially perturbed
   candidates — Prager–Synge holds for ANY conforming candidate,
   which is exactly what makes untrusted proposers safe (ver-001,
   120/120).
2. EFFECTIVITY: median bound/truth = 1.000 on the Galerkin battery
   (band ≤ 3; the ~30% accept-rate kill criterion is unreachable with
   loose-but-sound bounds — soundness alone does not close the
   economy), zero tightness failures (ver-002).
3. Interval soundness: near-ulp enclosure widths; NaN/∞ candidates
   FAIL CLOSED; wild-but-finite candidates stay finite and rejected
   (ver-003).
4. G5: verdicts, bound endpoints, and flux hashes are BITWISE
   reproducible; equality accepts are sound by domination; single-
   and zero-interior-DOF meshes bound truthfully (ver-004).
5. Certify-the-certifiers: an injected UNSOUND estimator (bound/10)
   undershoots truth and is CAUGHT by the harness (a fooled bound is
   a Sev-0 wrong answer wearing a badge); the hierarchical family
   stays within its stated band (ver-005).
6. The warm start saves ≥ 1.5× Newton iterations with an ESTIMATED
   color and complete ledger rows (ver-006).
7. Zoo: answers exist only through the verifier; confidence is
   advisory in BOTH directions (a NaN-confidence good proposer still
   wins; a confidence-1.0 adversary never does) (zoo-001).
8. Warm adjoints beat zeroth-order >2×; both degrade into verified
   accepts; equidistant ties are deterministic (zoo-002).
9. Coarse-rung candidates accept at honest tolerances, reject at
   tight ones, and fp16-quantized candidates still verify (zoo-003).
10. THE FALSIFIER: an adversarial surrogate lands ZERO incorrect
    accepts despite satisfying the nodal boundary conditions, its rate
    collapse auto-demotes it, and demoted proposers stop being consulted
    (zoo-004).
11. The economics loop stays sound end-to-end with per-proposer
    per-regime rows shipped to the ledger (zoo-005).
12. Outright accepts ship without a solve; warm starts measure real
    savings on the strongly nonlinear class (econ-001); antithetical
    warm starts clamp to zero recorded savings with the raw delta
    preserved (econ-002).
13. Drift demotion is LOCALIZED to the collapsed regime (econ-003);
    hysteresis prevents flapping and genuine recovery re-promotes
    (econ-004); priors are conservative, decisions deterministic, and
    the dashboard ships via fs-obs (econ-005).
14. Hostile public inputs — short and oversized meshes, invalid
    domains/order, length mismatches, non-finite values, noncanonical
    boundary data, invalid tolerances, excessive polynomial degree,
    tampered derived polynomials, and invalid MMS endpoints — refuse
    before compute with no evidence authority. Cancelling valid MMS
    endpoints remain admissible, and hostile problem names cannot
    inject ledger JSON (ver-007).
15. The interval arithmetic battery independently gates positive
    division using the sign of an FMA residual. A double-double oracle
    separately proves that element width, midpoint, half-width,
    candidate slope, mapped nodes, and mapped weights enclose the
    exact values represented by their f64 inputs; removing input
    widening fails the regression.
16. fem1d public entry points reject empty/short/oversized meshes,
    nodal shape/BC/non-finite errors, invalid canonical polynomial
    state, reciprocal-width failure, work overflow, allocation
    failure, and unusable pivots before returning an ordinary value.
    A zero-update nonlinear run reports finite nonconvergence unless
    its initial residual already satisfies the tolerance.
17. Zoo and economics propagate malformed queries and finite
    nonconvergence as errors. No failed cold/warm run can become a
    solve count, iteration saving, drift observation, or dashboard row
    (zoo-006, econ-006).
18. Coarse-rung prolongation uses one monotone segment cursor, not a
    full mesh search per fine node; a 4,097-node strongly nonuniform
    regression preserves every injected coarse value bitwise
    (zoo-008).
19. Stateful proposers run once per speculation. All-rejected
    economics consumes the retained first-pass candidate (econ-007),
    while an outright winner preserves every earlier rejection in
    drift telemetry (econ-009). Drift/zoo key counts, counters, rolling
    savings windows, and policy inputs fail closed at fixed bounds.
20. G4 verifier traces lock exact phase entries, invocation-global 256-unit
    boundaries, structured-refusal flushes, callback-error precedence, and the
    distinct final publication callback. Sparse-callback execution is bitwise
    equivalent to the legacy convenience wrapper.

## Error model

Fail closed is the error model: no exceptions cross the boundary; a
malformed input or unevaluable intermediate is a structured REFUSAL
with an unbounded sentinel, zero hash, and no color, never a panic and
never an accept. A well-formed finite bound above tolerance is a
distinct, ordinary REJECT. Non-authority numerical APIs return
`Fem1dError`; finite Newton nonconvergence remains explicit in the
low-level report and is promoted to an error by every savings caller.
`verify_with_checkpoint` adds a separate execution-control result layer:
callback errors are returned unchanged and no report is published. Shape and
scientific failures remain fail-closed `VerifierReport` refusals when no
callback error intervenes; work-plan overflow or an executed/preflighted-count
mismatch is a structured verifier refusal.

## Determinism class

Fully deterministic and sequential in v0: fixed quadrature, fixed
Thomas elimination order, bit-nudge rounding. Bit-identical reports
across runs (ver-004). Thread-count independence rides the
deterministic-reduction contract when the tile-kernel form lands
(no-claim below).

## Cancellation behavior

The equilibrated verifier is synchronous and bounded by mesh, coefficient, and
checked work caps. `verify_with_checkpoint` exposes bounded logical progress to
a caller-provided callback without depending directly on `fs-exec`; a parent
may map callback invocations to `Cx::checkpoint`, another cancellation token,
or retained trace collection. Callback failure prevents publication. Other
legacy fem1d/zoo/economics loops do not thereby gain `Cx` cancellation; their
existing synchronous resource caps remain their boundary.

## Unsafe boundary

None. `#![forbid(unsafe_code)]` via workspace lints; no capsules.

## Feature flags

`certified-speculation` — the whole verifier, OFF by default per the
Ambition-Tag gating rule until its Gauntlet tier and the six-month
kill-metric (accept rates > ~30% at customer tolerances, warm-start
savings ≥ 1.5×) stay green.

## Conformance tests

`tests/conformance.rs`, cases ver-001..ver-007; `tests/zoo.rs`, cases
zoo-001..zoo-008; `tests/economics.rs`, cases econ-001..econ-009; and
the fem1d unit refusal/convergence/Gauss batteries. The identity battery
independently mutates producer domain/version, every class transport field,
both nested problem bindings, and mesh bytes; a one-ULP exact-solution change
with identical six-digit display still moves the root. It locks the v2 class
root `0x959a77719f308c27` at 281 bytes and problem root
`0x7148ea04d6605664` at 490 bytes, and refuses stale/future receipts. Each
completed aggregate case emits a canonical fs-obs `ConformanceCase`: the
randomized ver-001 battery records its literal LCG input seed
`0x1001_2026_0707_0091`, while fixed ver-002..ver-007 cases record zero.
Existing fs-obs companion events retain the effectivity table, ledger row, and
structured refusal row. Assertions and expectations that abort during setup or
per-input checking precede aggregate verdict emission and remain ordinary Rust
test diagnostics; the suite does not claim a canonical aggregate record for
those early exits. Any reimplementation must pass the suite unchanged.
The G4 additions specifically cover exact work-plan formulas, sparse callback
traces, legacy equivalence, invocation-global boundaries, refusal progress,
callback errors, and publication distinct from a boundary at the same count.
The retained-receipt unit regressions lock exact replay before positive or
rejected result access. The `PresentedVerifierReceipt` rustdoc compile-fail
cases lock the absence of bounds, acceptance, and color APIs on presented
transport data.

## No-claim boundaries

- v0 is the 1D elliptic class with polynomial data (quadrature
  exactness is the rigor backbone). The 2D/3D FEEC H(div) patchwise
  equilibration (Braess–Schöberl / Ern–Vohralík) rides fs-feec's
  Whitney machinery as the successor — the architecture (accept test,
  colors, falsifier, fail-closed) is class-independent and lands here
  unchanged.
- The `u` endpoint check is manufactured-solution metadata
  consistency, not a premise used to issue the PDE certificate. The
  certificate is for the homogeneous-boundary problem defined by the
  recomputed canonical forcing and the conforming candidate; it does
  not claim that rounded Horner evaluation is symbolic algebra.
- Canonical bytes and their owner-local receipts are replay identities for the
  producer-v2 MMS schemas. They do not authenticate provenance or replace a
  ledger signature; consumers must still use the evidence and ledger trust
  boundaries.
- `ReplayIdentity::clone` and its formatted display remain allocation-infallible
  lower-layer conveniences; this contract's fallible construction guarantee
  covers MMS admission and identity minting, not arbitrary downstream cloning.
- Variable diffusion coefficients, non-polynomial data (with data-
  oscillation terms and explicit Poincaré constants), and quadrature
  ERROR bounds for transcendental integrands are the same successor.
- The roofline tile-kernel form of patchwise equilibration (with the
  stated ISA acceptance bands) belongs to the perf lane.
- Interval arithmetic here is local; unification with fs-ivl is the
  workspace-wide interval consolidation.
- Accept-rate telemetry at customer-realistic tolerances (the kill
  measurement) accumulates once the first physics vertical is live.
- The NEURAL surrogate proposer (FrankenTorch fp16/fp8 over
  FrankenNetworkx graphs) is fs-surrogate's (bead 7tv.8); the zoo's
  trait and adversarial gating are ready for it.
- The asupersync speculative RACE (proposer vs target, loser drained
  request→drain→finalize at a tile boundary, zero leaks) is the
  concurrent form — v0 speculates sequentially before the solve.
- Logical work and the 256-unit callback stride are stable accounting
  semantics, not instruction, wall-clock, memory, deadline, drain-latency,
  pause/resume, or 200-microsecond responsiveness claims. The verifier report
  itself does not bind work-plan/poll-policy versions or stride; downstream
  evidence that claims this execution policy must bind the plan fields and all
  three public policy constants explicitly.
- The generic callback seam is not itself `Cx` integration and does not prove
  request→drain→finalize behavior, leak freedom, or cross-ISA/thread-count bit
  stability for a surrounding async or future parallel scope.
