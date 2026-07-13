# CONTRACT: fs-surrogate

Learned accelerators with guarantees: surrogates permitted only inside certified
validity bands — ML proposes, certified numerics disposes.

## Purpose and layer

Layer L4 (surrogate / ROM). The default core is dependency-free and pure Rust
(an in-house symmetric eigensolver for the method of snapshots); the optional
ladder feature depends downward on `fs-evidence`, `fs-exec`, `fs-alloc`, and
asupersync for evidence, bounded tile execution, memory admission, and live
task cancellation.

## Public types and semantics

- `pod(&[Vec<f64>], energy_threshold) -> Result<Pod, SurrogateError>` — a POD
  reduced-order model via the method of snapshots (correlation matrix `SᵀS`,
  symmetric eigendecomposition, modes `φₖ = Svₖ/σₖ`), retaining the fewest modes
  capturing `energy_threshold` of the mean-centered energy.
- `Pod` — `rank`, `energy_captured`, `project`, `reconstruct`,
  `reconstruction_error` (the reduced-vs-full error).
- `conformal_band(residuals, alpha) -> ConformalBand` — the distribution-free
  split-conformal band (the `⌈(1−α)(n+1)⌉`-th smallest residual); `covers`,
  `half_width`. `empirical_coverage(&ConformalBand, &[(pred, truth)])`.
- `certify_or_escalate(&ConformalBand, in_validity_domain, decision_tolerance)
  -> Decision` — `UseSurrogate` iff inside the domain AND the band is at least as
  tight as the decision tolerance, else `Escalate`.
- Root `SurrogateError` — `NoSnapshots` / `DimMismatch` / `BadThreshold`.

- `ladder` module (addendum Proposal A, bead knh1.4; [F], behind
  `abstraction-ladder`): a bounded abstraction ladder whose present
  authority is ESTIMATED, not certified. `TruthModel` defines the P1
  full-order elliptic family's DECLARED level-0 semantics; its f64 solve
  is not an enclosure. `RbLevel` uses offline snapshots, an
  energy-orthonormal basis, and online Galerkin evaluation of the
  textbook residual/coercivity a-posteriori estimator. `ConceptLevel`
  uses interpolation with total dispersion calibrated at finite probes
  as `|concept − lower RB| + lower RB QoI estimator`; admission also
  evaluates that quantity at the actual query and takes the larger value.
  `Ladder::at_level(k)?.query(μ, tol)` performs AUTOMATIC BOUNDED
  DESCENT: an RB/concept rung answers only when its estimator is within
  tolerance; otherwise the leak is recorded and the query descends.
  `ladder::SurrogateError` names all ladder refusals. `RbCoveragePlan`
  retains the exact ordered battery and ladder shape. `rb_coverage_scoped`
  executes one logical parameter tile per μ and returns a complete Estimated
  fraction or an incomplete no-claim plus deterministic prefix/progress
  receipt. `rb_coverage` remains the smaller synchronous compatibility oracle.

## Invariants

- POD reproduces an exactly-representable (low-rank) snapshot set to roundoff;
  its modes are orthonormal; the retained rank captures `>= energy_threshold`.
- The conformal band achieves at least its nominal `(1−α)` empirical coverage on
  exchangeable held-out data.
- `certify_or_escalate` uses the surrogate ONLY when trustworthy (in-domain +
  band tight enough), so a fleet of queries costs strictly less than
  all-high-fidelity whenever any query is served by the surrogate.
- Every ladder-emitted color is `Estimated` and passes the shared
  `fs-evidence` payload validator. RB answers carry the f64-evaluated
  QoI estimator as dispersion; concept answers carry the larger of the probe
  maximum and query-local cross-rung discrepancy PLUS lower-rung QoI
  dispersion, so agreement with an inaccurate RB cannot erase its
  uncertainty. Level 0 carries
  infinite dispersion because an unproved floating-point solve makes no
  spread claim.
- Ladder state is sealed. Truth dimension, training range, basis,
  calibrated dispersion, rung collection, family identity, and answer
  evidence cannot be mutated or forged through public fields. Every rung
  is bound to one identity containing the truth dimension and exact
  floating-point range endpoints.
- Public ladder arithmetic and lookup operations are fallible. Queries
  reject non-finite/non-coercive/out-of-range inputs before lookup, and
  generated training/probe grids must be strictly increasing in f64.
- Ladder construction preflights nonempty, capped, strictly decreasing
  requested RB dimensions plus checked aggregate memory/work budgets
  before the first snapshot. After orthogonalization, actual retained
  dimensions must also strictly decrease before a rung is stored.
  Coverage batteries are nonempty and capped on both axes, their Cartesian
  product, and conservative aggregate work. Each parameter performs at most
  one descent (including at most one truth fallback), and the resulting RB
  estimators classify every requested tolerance without repeating solves.
- Production coverage validates and retains exact IEEE-754 input bits, rung
  dimensions, family, Cartesian count, and conservative work before execution.
  A worst-case live-scratch charge is reserved on the caller's bounded
  `OperationMemoryLease` before slots or workers are created. Results commit to
  unique parameter slots only after their final checkpoint; final aggregation
  is ascending by parameter index and ignores out-of-order completions beyond
  the first gap.
- Incomplete coverage cannot expose a fraction through its type. Its receipt
  contains only the longest fully committed parameter prefix and the first
  unfinished parameter's operational progress, together with an absorbing
  `NumericalCertificate::no_claim()`.

## Error model

Structured `SurrogateError`. Ladder construction, energy, compliance,
training, lookup, level selection, querying, and coverage return named
errors for invalid shapes/values/ranges/grids, singular or non-finite
derived arithmetic, and resource excess. The non-ladder conformal helper
still panics on nonsensical inputs (empty residuals, `α ∉ (0,1)`).

## Determinism class

Fully deterministic for synchronous work. A completed scoped coverage result is
bit-identical across worker counts and steal schedules because each parameter
has a unique slot and final integer aggregation is in parameter order. The
semantic receipt excludes timing/steal data. A cancellation triggered at the
same logical checkpoint replays the same retained prefix; wall-time-triggered
cancellation may naturally observe a different prefix and remains no-claim.

## Cancellation behavior

The POD/conformal API and `rb_coverage` compatibility wrapper remain
synchronous; `rb_coverage` has the explicit smaller
`MAX_SYNCHRONOUS_COVERAGE_WORK_UNITS` cap and no interruption claim.

`rb_coverage_scoped` requires a live asupersync `Cx`, `TilePool`, external
`CancelGate`, declared `RunId`/per-tile `Budget`, and bounded operation memory
lease. It polls both the tile gate and ambient task at every parameter/rung
phase boundary and after at most 256 logical scalar updates inside allocation
initialization, validation, assembly, Thomas solves, dense elimination,
reconstruction, residual/Riesz work, and tolerance classification. A request
raises the shared gate, all workers drain and join, scratch is released, and
only then is an incomplete no-claim returned. A final ambient checkpoint occurs
after drain and before a coverage fraction can be published.

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

- `abstraction-ladder` [F] (default OFF) — the estimated, leak-alarmed
  abstraction ladder (knh1.4/y6yv, Proposal A; `dep:fs-evidence`,
  `dep:fs-exec`, `dep:fs-alloc`, `dep:asupersync`); gates the `ladder`
  integration target.

## Conformance tests

`tests/surrogate.rs` (8 cases): POD reproduces a low-rank set exactly;
orthonormal modes; energy-based rank + reduced error; bad-input rejection; the
conformal band achieves nominal coverage; certify-or-escalate uses the surrogate
only when trustworthy; the policy reduces cost vs all-high-fidelity;
determinism.

`tests/ladder.rs` (feature-gated): f64 RB estimator containment on the
elliptic fixture, bounded descent, Estimated-only payload authority,
deterministic replay, structured hostile-input refusals, representable
grid and family binding, requested/retained fidelity descent,
lower-rung uncertainty inheritance, pre-training memory/work limits,
bounded coverage batteries, complete scoped replay across worker counts,
pre-cancel and budget-exhaustion storms, deterministic retained prefixes,
bounded observation latency, memory/arena quiescence, and successful pool reuse.

## No-claim boundaries

- v0 is the CLASSICAL ROM core (POD via method of snapshots) + the conformal /
  certify-or-escalate guardrail. NEURAL OPERATORS (Fourier neural operators,
  DeepONets via FrankenTorch), DEIM nonlinear-term interpolation, BALANCED
  TRUNCATION for LTI subsystems, and KOOPMAN/DMD are the fuller deliverable,
  staged.
- The eigensolver is a small dense Jacobi for the snapshot correlation matrix;
  the production path is fs-la randomized/TSQR SVD over large snapshot matrices.
- The conformal band is SPLIT-conformal (exchangeable data); the anytime-valid
  e-value formulation with online recalibration under drift is the
  conformal-hardening follow-on.
- Continuous training from the ledger, versioned/model-carded surrogate
  artifacts, and design-family-respecting splits are downstream integrations.

## No-claim boundaries (ladder)

- The beachhead covers the AFFINE-PARAMETRIC ELLIPTIC regime (1-D
  fixture family here); nonlinear/transient coarse levels are the
  research frontier and enter only as estimated-color concept rungs.
- Level 0 is the declared FE semantics, but neither floating-point solve
  error nor FE discretization error is enclosed here. Its Estimated
  color therefore has infinite dispersion; there is no zero-error claim.
- The RB residual/Riesz/solve path is evaluated in round-to-nearest f64
  without outward rounding or independent linear-solve certificates.
  Its textbook estimator is useful and tested for containment on this
  fixture, but it does not authorize `Color::Verified` and descent is not
  called certified.
- Compliance dispersion includes both the squared energy estimator and the
  floating reduced solve's computable Galerkin defect
  `|f(u_rb) - a(u_rb,u_rb)|`; exact orthogonality is never assumed.
- The concept rung's dispersion is a finite-probe MAXIMUM of
  `|concept − lower RB| + lower RB QoI estimator`, augmented by the same
  query-local quantity. Neither is an enclosure over the continuous range.
  The Estimated color is load-bearing.
- The synchronous `rb_coverage` helper is only a small compatibility oracle and
  makes no interruption claim. Production-scale authority requires
  `rb_coverage_scoped` and its explicit execution/memory inputs.
- Scoped coverage's static scratch reservation covers its bounded vector,
  matrix, slot, and receipt payload plan plus executor-charged arena/root
  metadata. Thread stacks, allocator bookkeeping, the immutable pre-existing
  ladder, and OS scheduling latency remain outside the numerical claim. A
  callback or dependency that ignores its supplied cancellation protocol would
  likewise be outside claim; the current sealed arithmetic path does not do so.
- A completed coverage fraction still classifies round-to-nearest f64 RB
  estimators. It is `Estimated` with infinite dispersion, never a certificate
  that the continuous model error lies below the requested tolerances.
- The eventual certificate destination is an outward-rounded residual,
  Riesz solve, reduced solve, coercivity floor, and QoI enclosure whose
  complete arithmetic path is independently checkable. Only that path,
  once admitted by the Gauntlet, may upgrade a rung to `Verified`.
- Per-REGION (spatial) RB decomposition and the fs-ir at_level query
  integration are the named growth seams.
