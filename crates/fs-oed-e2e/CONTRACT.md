# CONTRACT: fs-oed-e2e

SensorForge — optimal experimental design that knows when to stop. Layer L4
(ASCENT).

## Purpose and layer

Composes `fs-assimilate` and `fs-exec` (one explicit cancellation context and
private invocation-global poll ledger), `fs-voi` (EVPI + recommend),
`fs-toleralloc` (first-order budget), `fs-evidence` (Estimated lineage), and
`fs-blake3` (bounded canonical identities). Dependencies point downward.

## Public types and semantics

- `Candidate::new(...) -> Result<Candidate, CandidateError>` checks and then
  seals the name, finite truth/prior, non-negative prior variance, and positive
  finite sensor noise/cost. Read-only accessors expose the declaration.
- `run_campaign(&[Candidate], threshold, max_sensors, &Cx) ->
  Result<OedReport, OedError>` validates the threshold, unique candidate
  identities, and explicit synchronous work caps, then CANONICALIZES the
  candidate menu (name-ascending) EXACTLY ONCE at admission (bead
  sj31i.62) before greedily placing VoI-chosen sensors under
  deterministic cancellation checkpoints. Every derived belief/estimate/
  action/report sequence follows canonical order, so campaign outcomes
  and retained artifacts are INVARIANT under caller menu permutation.
- Canonical EVPI representation (bead sj31i.62): estimates live in an
  opaque `CanonicalDesignMenu` whose identity order is verified once
  (one O(n) window scan per posterior refresh, `CanonicalOrderViolated`
  on breakage) and never re-sorted; EVPI evaluates through
  `fs_voi::evpi_by` with NO per-call clone or sort, and the predictive
  quadrature substitutes the target's mean/statistical uncertainty
  through a non-owning, lifetime-bound `MeanOverrideView` (validated
  index and finite payload, `OverrideInvalid` otherwise) — stale
  restoration or cancellation-corrupted scratch state is
  unrepresentable because nothing is mutated. Target lookup is
  O(log n) binary search on the canonical order. The old per-node
  clone-and-sort work was never charged by the work plan; eliminating
  it makes the retained record-visit accounting HONEST without
  changing any admitted/realized work formula, pinned work constant,
  or poll boundary.
- `demo_candidates() -> Result<Vec<Candidate>, CandidateError>` builds the four
  checked candidates whose close top two drive the decision.
- `OedReport` includes the native EVPI trace, posterior summaries, one
  instrument-bound `assimilation_color` per placement, and input-bound final
  variance/EVPI colors so consumers do not need to transcribe the loop.
- Both retained report estimators use identity version 6 (v6 binds the
  CANONICAL candidate sequence — bead sj31i.62). Their canonical
  preimages bind every candidate declaration, threshold and placement cap,
  execution mode, stream/budget context, admitted and realized work shapes,
  poll/planning policy versions and strides, quadrature rule constants, every
  realized output sequence, and both color dispersions.
  Prefixes are `sensorforge-posterior-variance:v6:` and
  `sensorforge-evpi:v6:`; the manifest also locks planning policy v3, poll
  policy v2, `fs-voi` EVPI semantics, record stride 256, and action stride 1.

## Invariants

- Sensor planning uses the same scalar Kalman variance model as execution:
  `P' = PR/(P+R)`, evaluated with an overflow-safe form. The action effect is
  rebuilt after every placement from the current `P` and the candidate's
  declared noise `R`; no fixed sensor reduction is permitted.
- Action value integrates posterior-mean outcomes with a retained deterministic
  nine-point normal Gauss-Hermite rule. Sensors therefore land on candidates
  that are both decision-relevant and informative at their declared noise and
  cost; each completed placement shrinks total posterior variance.
- The campaign evaluates STOP before its placement cap, including when
  `max_sensors == 0`. `decision_robust` is true only when final modeled EVPI is
  at or below the checked threshold; a no-useful-action stop above threshold is
  not mislabeled robust.
- The posterior variance is `Estimated`: the scalar Kalman formula is exact for
  its declared linear-Gaussian model, but neither floating-point roundoff nor
  model-form assumptions carry an interval certificate. The EVPI stop is
  `Estimated`. Both bounded estimator identities commit to the complete ordered
  candidate declarations, threshold, placement cap, realized placement and
  posterior sequences, and canonical assimilation colors.
- Zero total prior variance has fractional reduction `0.0`, not NaN. A
  zero-sensitivity candidate receives `+infinity` allocation, the exact
  unconstrained first-order optimum; positive-sensitivity allocations must be
  positive and finite.
- Deterministic (the worked-campaign readings hit each truth; Kalman variance is
  observation-free; planning quadrature is fixed-order). Equal-mean EVPI inputs
  are canonicalized by candidate identity before the current top-two
  approximation, so caller menu order cannot alter the sensor policy. This does
  not upgrade that approximation into full multi-alternative EVPI.
- With `n` candidates and admitted placement cap `s`, preflight charges setup
  `5n`, each complete placement `11n^2+5n+2`, maximum finalization
  `18n+6s+5`, and their checked sum. For `p` completed placements, `a` action
  rounds, `m` positive-prior candidates, and `d=a-p` (necessarily zero or one),
  realized accounting replaces the maximum tail with
  `13n+5m+6p+5` and charges an incomplete action round, when present, as
  `d(11n^2+2n)`. Publication requires exact equality with that realized shape
  and never exceeds the admitted shape.
  Admission also requires `11n^2*s <= MAX_CAMPAIGN_EVALUATIONS`.

## Error model

No documented input panic. Candidate and campaign rejection is structured as
`CandidateError` / `OedError`; lower-layer assimilation failures retain their
`AssimError` source. `Cancelled`, `AssimilationCancelled`, ordinary nested
`Assimilation`, and `WorkPlanMismatch` remain distinct, so nested quota or
cancellation cannot be laundered as a scientific refusal. Resource admission
caps candidates, placements, and the quadratic action/design work multiplied
by the retained expectation-rule cost before campaign allocation or iteration.
Derived posterior variances, posterior means, expected EVPI, and value-per-cost
must remain finite.

`BudgetRefused` (bead sj31i.6) retains the ambient accountant's typed refusal
verbatim: `run_campaign` admits `cx.budget()` plus the admitted work plan
through `fs_exec::AdmittedBudget` before any work (an expired deadline —
`Budget::ZERO` included — a deadline without an ambient time source, or a
plan exceeding the cost quota refuses at admission), observes
cancellation/deadline at every campaign checkpoint, and accrues completed
work as retained cost consumption. The poll ledger remains the raw shared
counter because nested assimilation consumes the same slice; its exhaustion
keeps the structured `Cancelled` shape.

## Determinism class

Same-process deterministic for a fixed implementation/toolchain manifest (G5).
Equal value-per-cost actions use their canonical action identity as the
order-independent tie-break. This is not a cross-ISA identity promise.

## Cancellation behavior

Synchronous and cancellation-aware. A private invocation-global ledger derives
its poll quota once from `Cx`; only campaign checkpoints and the nested
assimilation transaction may decrease it. `Cx` is polled at deterministic
admission, action (every action), 256-record, assimilation-commit, refresh,
finalization, identity, and publication boundaries. Campaign cancellation,
nested assimilation cancellation, or quota exhaustion never publishes a
partial report. Deadline and cost-quota enforcement remain cross-workflow
follow-up scope: those budget fields are identity-bound provenance but only the
poll ledger is consumed locally. A scratch posterior is committed only after
nested assimilation and the placement-commit checkpoint both succeed; final
publication additionally requires both identity manifests and exact realized
work reconciliation.

## Unsafe boundary

None; `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/oed.rs` plus the crate unit tests (G0/G3/G4/G5): decision-relevant
placement and stopping; low-noise versus high-noise ordering under adversarial
menu permutations; predicted/realized Kalman variance agreement and extreme
finite noise limits; initial STOP at a
zero placement cap; full-report determinism; adversarial candidate/campaign
input rejection; zero-variance behavior; cancellation/poll bounds;
unmeasured-input evidence binding; instrument-bound assimilation lineage;
admitted/realized work locks; hostile-maximum admission; nested and finalization
quota sweeps; execution/budget/work identity binding; exact quadrature-bit
semantics; and sealed-output identity movement.

## No-claim boundaries

- Scalar (diagonal) beliefs and a Gaussian EVPI model; the greedy VoI policy is
  not claimed globally optimal. Nine-point Gauss-Hermite outcome integration is
  an Estimated deterministic approximation and currently has no certified
  quadrature-remainder bound. The precision budget uses a prior-std sensitivity
  proxy.
- `decision_robust` is a modeled EVPI criterion, not a physical decision
  certificate or independent validation. The worked campaign currently injects
  declared truth as its deterministic reading; a production measurement
  provider and stochastic outcome stream are separate required work.
- Logical-work and poll counts do not claim instruction counts, wall-clock or
  memory bounds, pause/resume support, deadline/cost enforcement, drain latency,
  or a 200-microsecond cancellation guarantee. Identity roots are
  replay/integrity bindings, not authenticated provenance or proof of global
  policy optimality.
- The G5 battery is same-process and fixed-manifest evidence; it does not claim
  cross-version or cross-ISA floating-point identity. As of identity v6 the
  report identity IS caller-menu-permutation-invariant (the preimage binds the
  canonical candidate sequence); v5 artifacts that bound a non-canonical
  declaration order remain valid only under their own version prefix. As of
  identity v7 the preimage additionally binds the admitted byte plan and the
  byte-accounting policy version; v6 artifacts remain valid only under their
  own version prefix.
- Full-menu decision algebra (bead sj31i.5, planning policy v4, byte policy
  v2): the STOP gate, initial/final EVPI, and the trace are the FULL
  multi-alternative `fs_voi::expected_opportunity_loss_by`; action valuation
  uses the SAME evaluator through the override view (one algebra — an
  action's value is exactly the full loss it removes). One logical work unit
  is one record's participation in one bounded evaluation (the quadrature
  depth is the fixed `fs_voi::EOL_QUADRATURE_PANELS`); byte charges scale
  robustness-bearing and action-tile scans by `FULL_EOL_SCAN_SWEEPS`. The
  demo campaign now EARNS global robustness: it spends late single
  placements on the non-contender alternatives whose residual optimality
  probability blocks the certificate, instead of declaring robustness while
  ignoring them.
- Byte accounting (identity v7, byte policy v2): every bounded seam —
  admission scan and one-time canonical sort, belief/estimate/menu builds and
  their window verifications, EVPI scans, action construction, each
  quadrature-view action evaluation (view construction reads plus per-node menu
  scans), chosen-action lookup, observation, placement commit and retention,
  finalization scans, allocation and posterior retention, and both report
  identity preimages/digests — charges a deterministic formula-based BYTE
  bound evaluated on the actual shape, with the worst-case envelope preflighted
  at admission under checked arithmetic. Charges refuse typed
  (`ByteBudgetExceeded`) rather than exceed the plan; a refused charge
  publishes nothing. The consumed ledger is an accounting UPPER BOUND on bytes
  visited/compared/hashed/retained, not measured allocator traffic, and the
  retained subset is exactly the bytes live in the published report. Retained
  campaign outputs (placements, colors, EVPI trace, posterior summaries, the
  canonical candidate menu) reserve their admitted capacity fallibly
  (`OutputAllocationRefused`) before scientific work; transient scratch keeps
  ordinary infallible allocation and is covered by the charge ledger instead.
