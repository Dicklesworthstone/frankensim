# fs-constraint — CONTRACT

## Purpose and layer

L4 (ASCENT). The constraint CALCULUS (plan §9.1, patch Rev F):
constraints with SEMANTICS, not anonymous `g(x) ≤ 0` — a typed kind
taxonomy with per-kind optimizer treatments, evidence-typed
evaluations, and infeasibility DIAGNOSIS (minimal unsat cores, ranked
calibrated repairs) that turns optimizer failures into design
conversations. fs-opt hosts the expression graphs; this crate owns
what the constraints MEAN.

## Public types and semantics

- `ConstraintKind`: `Hard` (never traded → feasibility restoration),
  `Soft(PenaltyLaw::{Quadratic, Hinge})` (→ penalty term), `Chance
  {level, ChanceEstimator::MonteCarlo{samples, delta}}` (→ estimate
  then act on the BOUND), `Robust{half_widths}` and `Certification
  {ProofKind::{Interval, Sos}}` (→ prove or escalate), `Fabrication
  {process}` and `Code{standard}` (→ domain check; semantics named
  for the ledger). `treatment()` is the routing table.
- `evaluate(problem, spec, x, noise) -> ConstraintEvidence`: status
  (`Satisfied`/`Active`/`Violated`/`NeedsProof`/`Proven`/
  `BoundNotCleared`), EXACT violation certificates for algebraic
  graphs, active-set role, penalty per law. Chance kinds compute a
  Hoeffding lower confidence bound and report satisfied ONLY when the
  bound clears the level — `BoundNotCleared` exists precisely for the
  case where the raw empirical rate clears but the bound does not
  (the validity machinery is the feature). Certification kinds refuse
  "satisfied" pointwise REGARDLESS of how good `g(x)` looks.
  The v1 host gate first requires exactly one `Rn` variable and an `x`
  whose length equals that declared point dimension. Zero/multiple-variable
  problems and Sphere, SO(3), or Stiefel hosts return typed `InvalidDomain`
  before graph evaluation or a chance-noise callback can observe work.
  Chance draws are additive Euclidean point-coordinate offsets only; this
  boundary neither projects nor retracts draws for manifold variables.
  Public policy fields are admitted before graph evaluation: active-set
  tolerances, robust half-widths, and soft weights must be finite and
  non-negative; robust width count must equal point dimension; chance
  levels/delta/sample counts must satisfy their declared domains, and
  `1 - delta` must remain strictly inside `(0,1)` in binary64. A
  chance noise draw must match the point dimension exactly; extra
  components cannot be silently ignored. A
  computed soft penalty must also remain finite and non-negative. No
  malformed value may turn a violation into satisfaction or a negative
  objective reward. The currently declared fixed-horizon Monte-Carlo
  estimator emits a `HalfWidth` certificate; substituting an unbound
  e-value or no statistical certificate downgrades the ledger row rather
  than changing the estimator after the fact.
- `interval_eval` (the in-house prover): dependency-aware interval
  inclusion formulas over fs-opt graphs per node; refuses (teaching
  reasons) on division through zero, domain violations, negative powers,
  and PDE/stochastic nodes. The formulas are conservative over exact
  endpoint operations, while theorem-strength binary64 endpoint authority
  remains subject to the explicit outward-rounding no-claim below. Before
  memo allocation it checks the sealed root depth and
  aggregate admission-work receipt against fs-opt's default cap
  schedule, and the walk itself is EXPLICIT-STACK (reachability
  worklist + bottom-up arena-order sweep; bead frankensim-xf8v7), so a
  graph built under looser caps refuses typed and no admitted graph can
  overflow the call stack; the exact max-depth boundary is a G4 fixture.
  The public `interval_eval` boundary admits a zero-variable constant
  problem with an empty binding, or exactly one `Rn` host with exactly
  one interval per component. It refuses non-finite, reversed,
  unrepresentable-span, missing, or extra boxes before graph evaluation
  and never aliases several host variables onto one binding vector.
  `prove_interval`
  accepts only `Certification { proof: Interval }`; it refuses SOS and
  non-certification requests rather than relabeling an interval result.
  Proof subjects accept the same zero-variable constant or one-`Rn`
  host shapes; multi-variable/manifold or dimension-mismatched subjects
  refuse rather than sharing bindings.
  A cleared structural bound returns the current `Proven` schema status
  plus a sealed `ProofArtifact` carrying both interval ends and a
  `ProofSubject` bound to the
  full-width admitted fs-opt problem semantic identity, exact node, and
  ordered endpoint bit patterns. Every `Proven` evidence value also
  retains those exact proof-bound bits; the artifact/evidence verifier
  requires the subject and both bound endpoints to match. Robust kinds
  use the same admitted uncertainty-box path and carry ENCLOSURE
  certificates; default theorem authority awaits the named outward
  successor rather than treating fp-slack sampling as a proof.
- `diagnose_infeasibility(problem, specs, domain, cx) -> Diagnosis`:
  elastic-relaxation solve (multi-start projected subgradient descent
  on total hinge violation, deterministic LCG starts) classifies
  feasibility and yields a witness or an unsat core seeded from the
  elastic support and refined by the DELETION FILTER — the core is
  MINIMAL (dropping any member restores feasibility). Repairs
  (relax-bound at graded slacks; drop-soft for soft members) come
  RANKED by Monte-Carlo feasible-volume estimates.
  Domain admission precedes solver allocation and evaluation: exactly
  one `Rn` variable, one range per point coordinate, finite ordered
  endpoints, and a finite span. Equal endpoints are valid fixed
  coordinates.
- `serialize_specs`/`parse_specs`: canonical `fscon v2` line form
  (floats as bit patterns), identical round-trips for admitted specs,
  and line-numbered refusals. The infallible writer preserves untrusted
  public field bits; the parser is the policy-admission boundary, so
  writer output from an invalid in-memory spec is intentionally refused.
  Dynamic tokens use canonical uppercase UTF-8 byte-percent
  encoding: literal percent and space are distinct (`%25` versus
  `%20`), unsafe bytes cannot split records, and `%` is the unambiguous
  empty-token sentinel. The parser rejects malformed/noncanonical
  encodings, legacy v1 ambiguity, trailing fields, and non-finite or
  negative raw-bit policy values.
- `ConstraintEvidence::to_ledger_row` and `Diagnosis::to_json`: the
  Rev S ledger row and the agent-facing diagnosis payload. All string
  fields receive complete JSON escaping. Publicly forged non-finite,
  negative, structurally inconsistent, or unbound required evidence
  emits a canonical `valid:false`/`status:no-claim` representation with
  a deterministic first reason and no positive satisfied/feasible
  headline; malformed numbers become `null` only after that downgrade.
  Diagnosis admission also requires its finite component violations to
  sum to the reported elastic total, and `elastic_solve` recomputes that
  published total from the final published component vector in canonical
  constraint order. A stale optimizer-carried aggregate therefore cannot
  authorize a feasible or infeasible headline.

## Invariants

1. All seven kinds map to their declared treatments; spec sets
   round-trip identically; ledger rows validate through fs-obs
   (fscon-001).
2. Statuses, roles, exact violation certificates, and both penalty
   laws evaluate as declared (fscon-002).
3. The chance BOUND decides: on an analytic uniform-noise fixture the
   raw rate clearing the level while the Hoeffding bound does not
   yields `BoundNotCleared`, never `Satisfied`; the half-width
   travels as a `StatisticalCertificate` (fscon-003).
4. Certification kinds refuse without artifacts; interval proofs
   succeed exactly on provable domains and refuse honestly otherwise;
   interval containment holds over random nonlinear boxes (G0);
   robust kinds carry enclosures (fscon-004).
5. Unsat cores are MINIMAL against enumeration: the core is jointly
   infeasible, every single deletion restores feasibility, bystanders
   are excluded, feasible systems return witnesses, and elastic
   verdicts match grid enumeration on the seeded fixture family
   (fscon-005).
6. Repairs are ranked by feasibility estimate, soft members offer
   drop actions, and estimates are CALIBRATED against exact
   enumeration (worst gap < 0.05 on the worked mass/strength
   example); the diagnosis payload ships through fs-obs (fscon-006).
7. Malformed elastic domains refuse before allocation/evaluation,
   while fixed axes remain admissible; hostile constraint/repair text
   cannot escape its JSON string and every emitted payload remains
   valid JSON.
8. Policy admission is identical for direct evaluation and raw-bit
   parsing; interval proof engines cannot satisfy a different declared
   proof kind; proof subjects change when problem, node, endpoint order,
   or even signed-zero endpoint bits change, and artifact/evidence
   verification includes both proof-bound endpoint bits.
9. Percent encoding is injective and writer/parser bytes are a fixed
   point. Invalid public evidence always loses positive claim authority
   with one stable reason rather than preserving a claim beside `null`.
10. Direct evaluation admits its single-`Rn` host and exact point length
    before graph/noise work. Unsupported product or manifold hosts cannot
    turn ambient-coordinate chance draws into an accidental manifold claim.

## Error model

`ConError` teaching errors: `NotScalar`, `Eval` (fs-opt errors carried
through), `NotProvable{why}` (an honest gap with escalation advice,
not a failure), `BadParam`, `ProofKindMismatch`,
`InvalidDomain(DomainError)`,
`Parse{line, what}`. `DomainError` distinguishes host variable count and
manifold, point-dimension representation, point/domain component-count
mismatch, and an axis-specific `InvalidRange` reason. The interval engine's
`IvalError` names each refusal reason, including the aggregate cap name,
observed count, and enforced limit.

## Determinism class

Fully deterministic: LCG-seeded multi-starts and Monte-Carlo streams,
canonical constraint ordering, bitwise float serialization. Identical
inputs give identical diagnoses, estimates, and bytes.

The conformance aggregate records distinguish input generation from
execution provenance. Randomized fscon-003, fscon-004, and fscon-005
carry their literal base input seeds `0x1001_2026_0707_0003`,
`0x1001_2026_0707_0004`, and `0x1001_2026_0707_0005`; fscon-003 derives
sample stream `s` as `base ^ s.wrapping_mul(0x9E37_79B9_7F4A_7C15)`.
Fixed-input cases use aggregate seed zero. Scoped-runtime cases
separately use the fixed Cx execution seed `0xC0C0`; that execution
identity is not misreported as randomized input provenance.

## Cancellation behavior

`elastic_solve` (and therefore `diagnose_infeasibility`) polls
`cx.checkpoint()` before solver allocation/evaluation and between scalar
graph evaluations at fixed strides through active-set construction,
constraint totals, restart/step/dimension loops, final evidence,
deletion filtering, and repair sampling. It returns the carried
`Cancelled` teaching error and does not publish a partial success. A
single synchronous fs-opt scalar evaluation is not yet internally
cancellation-tiled. Full evaluator tiling, cost/poll/deadline admission,
memory leases, and retained work-consumption receipts are explicitly
tracked by P0 successor
`frankensim-constraint-restoration-budget-receipts-x5sev`.

## Unsafe boundary

None. `#![forbid(unsafe_code)]` via workspace lints; no capsules.

## Feature flags

None.

## Conformance tests

`tests/conformance.rs`, cases fscon-001..fscon-006 — completed aggregate
cases emit `fs_obs::EventKind::ConformanceCase` with Info/Error severity,
failure-record linting, JSONL wire validation, and printing before the
aggregate assertion. Seeded LCG randomness follows the mapping above.
The object-shaped Custom companions for ledger rows and the full
diagnosis payload remain wire-validated and printed, alongside G4
shared-DAG/work and exact/max+1 depth-boundary fixtures for interval
evaluation. Adversarial fixed cases cover malformed domain
bounds/dimensions, fixed axes, negative/non-finite policy values and
their exact boundary neighbors, interval proof-kind confusion, exact
artifact subject binding, hostile/injective wire and JSON text,
positive-claim downgrade, deterministic round trips, and pre-cancelled
evaluation ordering. Any reimplementation must pass the suite
unchanged.

## No-claim boundaries

- The interval prover rounds to-nearest; outward-rounded arithmetic
  joins with fs-ivl (containment carries an fp-slack caveat until
  P0 theorem/kernel successor `frankensim-zup19` is green). SOS
  certificates are REPRESENTED (`ProofKind::Sos`), not executable —
  fs-sos is a later bead. Evaluation retains the explicit
  `NeedsProof { Sos }` intent, while the interval-only sealed artifact API
  refuses to mint or accept an opaque caller-authored SOS reference; this
  removes false authority without removing the planned SOS capability.
- `ProofArtifact` subject binding prevents accidental cross-problem,
  cross-node, or cross-domain reuse and its sealed fields prevent
  caller-side struct-literal minting. It is not external authenticity,
  independent proof checking, or a substitute for the future fs-sos
  verifier/trust channel.
- The elastic solver is small-fixture machinery (multi-start FD
  subgradient); the production feasibility-restoration solver is a
  later ASCENT bead. Nonconvex fixtures can defeat it — verdicts are
  cross-checked against enumeration only at conformance scale.
- Cancellation is now polled at bounded logical strides, but the
  current `Cx::checkpoint` path alone does not interrupt the interior of
  one scalar graph evaluation, enforce ambient deadline/poll/cost quotas,
  charge allocations to a lease, or retain a resource receipt. Those
  claims remain unavailable until the named P0 resource successor is
  green.
- Chance estimation is Monte-Carlo/Hoeffding v1; e-process anytime
  validity and richer estimators join with the UQ beads. Its current
  synchronous closure API has no `Cx`, work admission, cancellation, or
  resource receipt and must not be used for untrusted large sample
  counts; P0 successor `frankensim-oxyjg` owns the budgeted, replayable,
  anytime-valid redesign rather than hiding this gap behind an arbitrary
  cap.
- Fabrication/Code kinds carry semantics to the ledger; process
  models (fs-fab) and code-check rule packs bind in their beads.
- Repair generation covers bound relaxations and soft drops; material
  /topology switches (the patch's richer vocabulary) need fs-xform
  and fs-fab integration.
- Direct-evaluation and elastic-solver hosts are single-`Rn`-variable v1.
  Their public boundaries enforce that no-claim explicitly; multi-variable
  and manifold-variable domains generalize with the restoration solver.
- Assertions and expectations reached before an aggregate verdict are
  ordinary Rust test diagnostics; an early abort cannot claim that a
  canonical aggregate conformance record was emitted.
