# CONTRACT: fs-assimilate

Data assimilation as a living belief (plan addendum, Proposal 11): checked
sensor readings produce a content-addressed posterior candidate for a proposed
regime. Assimilation itself does not claim experimental validation.

## Purpose and layer

Layer L4 (inference/UQ). Depends on `fs-evidence` (`Color` and
`ValidityDomain`), `fs-exec` (explicit `Cx`, mode, stream identity, and
budgets), `fs-ivl` (outward-rounded PSD admission certificates), and the
in-tree `fs-blake3` identity primitive. This is the checked, deterministic
linear-Gaussian core of weak-constraint assimilation.

## Public types and semantics

- `Belief` — a Gaussian state belief with private state. `new(..., &Cx)`,
  `scalar`, `diagonal(..., &Cx)`, and `validate(&Cx)` return `Result`;
  admission requires a non-empty finite mean and a finite, square, symmetric,
  positive-semidefinite covariance within the `MAX_DENSE_STATE_DIM` v0
  envelope. `mean` and
  `covariance` expose read-only slices; `component_mean` and `variance` are
  bounds-checked.
- `diagonal_belief_invocation_resources(dimension)` returns the checked typed
  work, poll, cost, evaluation, conservative live-memory, and retained-output
  grant for diagonal construction. `Belief::diagonal_budgeted(..., &Cx,
  &mut ChildBudget)` consumes only that affine parent-issued grant.
- `Observation` — a private, checked scalar reading
  `value = operator·state + noise`. `Observation::new` rejects empty, all-zero,
  and non-finite operators; non-finite readings; non-positive noise; and
  malformed instrument leaf identities.
- `point_sensor(component, dim, value, instrument_noise, instrument) -> Result`
  — the REGISTRATION-FREE observation (`operator = e_component`), rejecting
  zero dimensions and out-of-range components.
- `scan_observation(..., registration_var, ...) -> Result` — adds a finite,
  non-negative registration variance to strictly positive instrument noise and
  rejects overflow of the total (R8).
- `misfit(&Belief, &[Observation], &Cx) -> Result<f64, AssimError>` —
  `Σ (h·mean − y)² / r`; empty sets and non-finite arithmetic are refused.
- `assimilate(&Belief, &Observation, &Cx)` / `assimilate_all(&Belief,
  &[Observation], &Cx)` — scalar Kalman fusion. Aggregate fusion uses canonical
  observation-content order, not caller order. Aggregate admission enforces
  `MAX_DENSE_OBSERVATIONS` and the multiplicative
  `MAX_DENSE_UPDATE_CUBIC_WORK` envelope before sorting or updating.
- `RobustObservation` — a private-state record constructed as available,
  missing, censored, saturated, or delayed. Pathology constructors validate
  instrument identity through the checked `Observation`, finite thresholds and
  ranges, saturation-at-endpoint, and positive sensor time constants.
- `ObservationBatch::new(records, covariance, &Cx)` — a checked robust batch.
  Covariance rows belong, in declaration order, only to available records;
  instrument identities must be unique; each diagonal must exactly equal the
  corresponding checked scalar observation variance; and the matrix must be
  finite, exactly symmetric, interval-certified PSD, and then strictly
  Cholesky-positive. The dense correlated dimension is capped at
  `MAX_CORRELATED_OBSERVATIONS = 256`.
- `assimilate_observation_batch(&Belief, &ObservationBatch, &Cx) ->
  RobustAssimilation` — excludes missing readings, returns a no-posterior audit
  for censored/saturated/delayed records, delegates an exactly diagonal usable
  batch to `assimilate_all` without changing its observations, and
  Cholesky-whitens a correlated batch before using the same scalar Joseph
  updater. `ObservationAudit` retains every instrument disposition, the
  effective observation degrees of freedom, terminal `BatchOutcome`, and a
  domain-separated receipt. A refused result cannot expose a partial posterior.
- `FiniteDifferenceSettings::new(relative_step,
  max_standardized_innovation)` — checked positive finite policy for
  scale-aware central differencing and the local innovation gate.
  `NonlinearObservationSpec::new(...)` binds a checked scalar reading, noise,
  instrument, caller-declared nonlinear model identity, and those settings.
- `linearize_nonlinear_fd(&Belief, &NonlinearObservationSpec, predict, &Cx)` —
  on success evaluates the supplied scalar model exactly
  `1 + 2 * state_dimension` times, retains the actual representable plus/minus
  states and model values in `FiniteDifferenceProbe` records, and constructs a
  local affine `Observation` only when
  `|reading - h(mean)| / sqrt(H P H^T + R)` is within the declared gate.
  `NonlinearLinearization` retains the prediction, innovation, local variance,
  standardized innovation, Jacobian row, probes, disposition, and
  domain-separated receipt.
- `assimilate_nonlinear_fd(...) -> NonlinearAssimilation` — routes an admitted
  affine observation through the established scalar `assimilate` updater.
  `DemotedLargeInnovation` exposes neither an affine observation nor a
  posterior; demotion is an audited no-publication result, not a failed
  Gaussian update.
- `assimilate_colored(&Belief, &[Observation], regime_param, lo, hi, &Cx) ->
  AssimilatedPosterior { belief, color, misfit_before, misfit_after }` — a
  read-only, regime-tagged **estimated candidate**. Access is through getters.
- `colored_assimilation_invocation_resources(...)` returns the checked typed
  grant and conservative payload envelope for a colored update.
  `colored_assimilation_invocation_resources_for_shape(...)` derives the same
  envelope from state dimension, observation shapes, and `ExecMode`, allowing a
  parent to seal admission before constructing a validated belief.
  `assimilate_colored_budgeted(..., &Cx, &mut ChildBudget)` reserves that
  envelope before temporary allocation, accounts the whole update through the
  same non-cloneable child, and publishes retained output only on success.
- `assimilate_colored_with_shared_poll_quota(..., &Cx, &mut u32)` is the
  compositional seam for a parent workflow that owns one monotonically
  decreasing poll slice. It rejects a supplied slice above the ambient quota;
  a raw counter does not authenticate provenance or prevent caller reissue.
- `AssimError` — structured invalid-state, dimension, bounds, identity,
  empty-input, noise, covariance, innovation, and non-finite-computation
  refusals. It implements `Display` and `Error`.

## Invariants

- Checked `Belief` and `Observation` values cannot be mutated into malformed
  states through the public API.
- At the external belief-admission boundary, covariance validation is an
  `O(n^3)` diagonal-pivoted Schur-complement certificate after dimensionless
  correlation scaling, with exact zero-variance-row handling. Every 2-by-2
  principal minor is first checked by an exact comparison of the input `f64`
  values as binary-rational products, so a sub-ULP negative determinant cannot
  be rounded onto the correlation boundary. One- and two-dimensional PSD
  inputs, including exact rank-one singular matrices, are thereby decided
  exactly. For active dimension three or greater, every scaled entry and Schur
  update also carries an outward-rounded `fs-ivl::Interval`; admission requires
  every pivot enclosure to be finite and wholly positive. A rounded positive
  point pivot is never authority. Singular PSD boundaries, sufficiently
  ill-conditioned strictly SPD matrices, and unresolved indefinite cases are
  refused with
  `CovarianceCertificationUnresolved`, distinct from a certified indefinite
  principal minor, rather than tolerance-clamped into a PSD claim or falsely
  described as indefinite.
- The scalar Kalman covariance update uses Joseph form
  `(I-KH)P(I-KH)^T + KRK^T`, computes one triangle and mirrors it for bit-exact
  symmetry, then passes the result through the same full validator as an
  externally supplied `Belief`. It does not silently clip eigenvalues or
  correlations. The checked dense implementation is `O(n^3)`.
- In exact arithmetic, fusing valid observations cannot increase component
  variances and the batch posterior cannot increase the weighted measurement
  misfit. Floating results are checked for finiteness, not interval-certified.
- `assimilate_all` and `misfit` canonicalize the complete observation records,
  making results bit-stable across input permutations while retaining duplicate
  readings and therefore their multiplicity.
- Scan noise is instrument noise plus non-negative registration variance; point
  sensors carry no registration term.
- Aggregate APIs reject empty or oversized observation sets. Dense update APIs
  additionally reject `observation_count * state_dimension^3` above the public
  work envelope before canonical sorting, identity materialization, or Joseph
  updates. Instrument and regime-axis identities use `fs-evidence`'s bounded
  leaf grammar; regime bounds must be finite and ordered.
- State and observation dimensions are capped at `MAX_DENSE_STATE_DIM` before a
  linear-size input can trigger quadratic dense allocation. The cap is 256 for
  this synchronous `O(n^3)` v0 core.
- A robust batch has one observation-noise authority. Its covariance diagonal
  and the available scalar observations' `noise_var` values must be bit-equal
  after signed-zero canonicalization. Unique instrument identities make
  covariance-row ownership unambiguous. The PSD certificate is the same
  interval-backed admission used by `Belief`; the additional strict Cholesky
  gate refuses singular PSD and unresolved positive-definite cases because
  whitening requires an invertible factor.
- Missing records contribute no synthetic value and reduce effective degrees
  of freedom. A censored, saturated, or delayed record refuses the complete
  publication boundary: otherwise-available records are reported as withheld,
  and no partially assimilated belief escapes.
- Robust audit receipts use
  `robust-observation-audit:v1:<64 lowercase hex>` and bind the prior, ordered
  records and pathology metadata, complete observation covariance,
  dispositions, terminal outcome, and published posterior if present.
- A nonlinear central-difference probe for component `i` requests
  `relative_step * max(1, |mean_i|)` and refuses unless binary64 addition and
  subtraction produce two distinct finite states, each distinct from the
  linearization point, whose difference is also finite and positive. The
  derivative denominator is the difference between those actual retained
  states, not twice the nominal step. Non-finite base or probe predictions and
  non-finite derivatives are typed refusals.
- Nonlinear admission is fail-closed. An admitted local model publishes
  `H = finite_difference_jacobian` and affine value
  `reading - h(mean) + H·mean`, so applying that observation at the
  linearization point preserves the original innovation. A disposition outside
  the declared standardized-innovation gate publishes no `Observation`; the
  wrapper consequently cannot publish a posterior. A linear scalar model with
  exactly representable probes reduces exactly to the existing linear updater.
- Nonlinear receipts use
  `nonlinear-fd-linearization:v1:<64 lowercase hex>` and bind the complete
  prior, reading, noise, instrument, caller-declared model identity, finite-
  difference and gate settings, prediction, innovation and variance, Jacobian,
  every actual probe value, and disposition. The arbitrary callback itself
  cannot be authenticated by this crate.
- The candidate's `Color::Estimated.estimator` is
  `assimilation-candidate:v4:<64 lowercase hex>`. Its BLAKE3 input uses typed
  length framing and binds the full prior, canonical observation multiset
  (including operator, reading, noise, instrument, and multiplicity), and
  proposed regime. The identity is deterministic, collision-resistant, bounded,
  delimiter-unambiguous, and changes with every bound semantic field.
- Candidate identity v4 additionally binds PSD admission policy
  `exact-2x2-interval-schur:v1` and its numeric policy version, execution mode, logical stream
  identity (all four `StreamKey` fields), every ambient budget field (deadline
  presence/value, poll quota, cost-quota presence/value, and priority), the
  effective poll-quota slice, every component of the checked work plan, and
  poll policy `fixed-stride:v3` with scalar stride 256, record stride 16,
  canonical-comparison byte stride 1,024,
  and identity-hash byte stride 1,024. The work plan itself has no separately
  claimed public “v4” version.
- For state dimension `n`, observation count `m`, canonical record sizes `b_i`,
  `b_max = max(b_i)`, and `L = ceil(log2(m))`, preflight accounts exactly for
  observation validation `m(n+4)`, record materialization `sum(b_i)`, canonical
  merge ordering `2m + mL + mL*b_max`, and each misfit pass
  `m(2n+4)`. An updating call also accounts for the prior clone `n+n^2`, each
  Joseph update `n^3+6n^2+8n+2+(n(n+1)/2)(n+1)`, and each PSD revalidation
  `n^3+6n^2+2`; colored calls add the exact bounded canonical-identity byte
  count. All sums/products are checked before work begins.
- Dense belief validation separately declares `6n+8n^2+n^3` scalar units.
  `misfit` enables one misfit pass without update/hash; plain assimilation
  enables update without misfit/hash; colored assimilation enables two misfit
  passes, update, and hashing. Successful data-dependent paths may consume less
  than the declared fail-closed bound.
- Typed planners retain logical work and abstract cost in distinct dimensions
  with equal numeric values, declare one scientific evaluation, and compute
  poll and byte envelopes before an affine child can spend. Diagonal memory
  covers the retained belief plus active-index, row-descriptor, scaled-cell,
  and pivot-column PSD workspace. Colored assimilation memory covers retained
  output and regime-row payload; canonical records, record descriptors, order
  and merge-scratch indices; old/new belief overlap; Joseph vectors/matrices;
  and PSD revalidation workspace. Output is a separate retained capacity.
- Budgeted entry points charge work, cost, evaluation, every poll, memory, and
  output to one borrowed `ChildBudget`. Scientific failures are latched into
  the parent invocation receipt. Unused capacity returns only when the child
  is consumed by `finish`; this crate cannot clone, increase, or reissue the
  parent authority.
- Candidate dispersion is `+infinity`, the shared explicit no-spread-claim
  sentinel. No API in this crate directly constructs `Color::Validated`.

## Error model

Structured `AssimError`; valid public calls return refusals rather than indexing
panics. Empty/ragged/indefinite/unresolved/non-finite beliefs, malformed observations,
invalid identities and regimes, non-positive or non-finite noise, dimension
mismatches, oversized aggregate count/work requests, degenerate innovations,
and finite-input arithmetic overflow are refused. `WorkPlanOverflow`,
`WorkPlanExceeded`, `PollQuotaExceedsAmbient`, `InvocationBudget`,
`BudgetRefused`, and `Cancelled` distinguish planning, accounting,
compositional-quota, typed invocation, ambient-budget, and
observed-cancellation failures; no partial belief/candidate is returned.
Ambient-path entry points (bead sj31i.6) admit `cx.budget()` plus the
preflighted work plan through `fs_exec::AdmittedBudget` before any work:
an already-expired deadline (`Budget::ZERO` included), a deadline on a
context without an ambient time source, or a plan exceeding the cost quota
refuses at admission; deadline and poll quota are enforced at every
checkpoint (cancellation takes precedence) and each advanced work unit
charges one cost unit. Poll exhaustion is a typed
`BudgetRefused(PollsExhausted)`, no longer conflated with `Cancelled`. Scientific preflight and domain refusals are latched fail-closed into
the invocation. Allocation failure remains
Rust's process-level behavior inside the admitted public resource envelope.
Preflight checks both retained covariance and temporary interval-certificate
matrix byte shapes before the initial checkpoint.
Robust-batch errors additionally distinguish covariance row/shape/finite/
symmetry/noise-authority failures, certified non-PSD input, unresolved PSD
certification, strict Cholesky failure, duplicate instrument identity, and
invalid pathology metadata. Robust covariance certification and the delegated
posterior update use the existing admitted core paths. The bounded Cholesky,
whitening, and audit-framing wrapper has cancellation checkpoints but does not
yet expose a separate typed invocation resource planner.
Nonlinear errors additionally distinguish invalid finite-difference/gate
settings, an unrepresentable perturbation, a non-finite model evaluation, and a
non-finite Jacobian coefficient. A large but finite standardized innovation is
a successful typed demotion rather than an error. The bounded nonlinear wrapper
has no separate typed invocation resource planner.

## Determinism class

Deterministic for a fixed execution manifest: fusion, misfit, candidate
identity, and posterior are pure functions of semantic inputs. Observation
permutation is canonicalized, so it does not change result bits or identity.
Duplicate observations are retained.
Changing the canonical schema or numerical algorithm requires a candidate
identity version/domain bump.
For robust batches, declaration order is semantic because it owns covariance
rows. Reordering records requires the matching row/column permutation; v1 does
not claim bit-identical Cholesky rounding or receipt identity across such
equivalent permutations. Replay of the same ordered batch is deterministic.
For nonlinear finite differences, state-component order fixes probe evaluation
and Jacobian order. Replay is deterministic only when the supplied callback is
a pure deterministic implementation of the declared `model_id` under the same
binary64 implementation/toolchain manifest. The crate cannot enforce callback
purity or prove that two callbacks sharing an identity are equivalent.

## Cancellation behavior

Synchronous and cancellation-aware. Long-running public APIs perform bounded
shape/work-plan admission before taking an explicit `Cx` checkpoint, then poll
before numerical work and before publication and within bounded scalar, record,
comparison-byte, and hash-byte tiles. Poll-quota exhaustion and
`Cx::checkpoint` cancellation return exact completed/planned logical work and
publish no partial result. The shared-quota seam supports a parent-owned
invocation ledger but cannot authenticate the raw counter it receives.
The budgeted forms instead poll their child authority, which checks an absolute
clock and the originating cancellation gate before consuming each poll. They
do not publish typed output after a deadline, cancellation, resource, or
scientific refusal.
Robust construction polls before admission, at each Cholesky row, and through
the reused PSD validator. Correlated whitening polls at each observation row;
audit hashing polls at each record and matrix row and once before finalizing.
Assimilation polls again before returning a result, so no posterior or audit is
published after an observed cancellation.
Nonlinear linearization polls before the base evaluation, before each
component's two probes, and before receipt publication. The wrapper polls again
after the delegated scalar update and before publishing its result. A model
callback is one bounded synchronous evaluation and cannot be interrupted from
inside this crate; callers needing finer cancellation must implement it inside
their callback.

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/assimilate.rs` (Proposal 11, G0/G3/G4/G5): mean shift and variance
reduction;
misfit reduction; point-versus-scan noise; exact permutation stability; honest
estimated color and bounded identity; empty/ragged/non-finite/indefinite belief
refusals; checked indexing and point-sensor bounds; malformed observation and
instrument refusals; scan-noise failures; empty aggregates; regime failures;
delimiter-collision and multiplicity identity probes; finite-input overflow;
strict rejection of near-boundary negative curvature; Joseph-form replay of a
cancellation-induced invalid-posterior counterexample; seeded randomized
post-update validation; dense-state, observation-count, and multiplicative-work
resource admission; deterministic replay; pre-cancel/final-publication
cancellation; exact quota sweeps through validation, ordering, update, PSD,
hash, and commit; hostile maximum-work cancellation; shared-quota depletion;
typed planner shapes, equality of validated-belief and pure shape preflights,
affine budgeted diagonal/colored execution, and receipt integrity; and
execution-mode, every budget field, every stream field, work, and poll identity
binding.

`tests/robust.rs` (bead `f85xj.12.4`, G0/G3/G4): exact replay of the existing
diagonal update; missing-record exclusion and effective-DOF accounting;
audited no-posterior refusals for censored, saturated, delayed, and all-missing
batches; correlated-noise posterior difference from a naive diagonal update;
posterior invariant revalidation; covariance PSD/noise-authority/identity
refusals; pathology metadata validation; and pre-cancellation.

`tests/nonlinear.rs` (bead `f85xj.12.4`, G0/G3/G4): exact reduction of a
dyadically probed linear model to the established scalar updater; actual
finite-difference probe retention; deterministic replay and model-identity
receipt sensitivity; nonlinear polynomial Jacobian recovery; large-innovation
demotion with no affine observation or posterior; malformed policy and identity
refusals; non-finite model/Jacobian and unrepresentable-step refusals; exact
callback count; and pre-cancellation.

## No-claim boundaries

- v1 is the LINEAR-GAUSSIAN assimilation with linear observation operators;
  full weak-constraint 4D-Var over a nonlinear forward model uses Proposal 1's
  CERTIFIED adjoints (adjoint-certs) and is the fuller deliverable — this crate
  is the Kalman/Gauss-Newton core.
- Robust v1 does not implement censored or clipped likelihoods, temporal state
  augmentation for delayed sensors, robust M-estimation, or
  sparse/shared-factor covariance beyond the dense cap.
  Censored, saturated, and delayed records are typed audited refusals, not
  silently Gaussianized data. A time constant is retained as metadata but is
  not a delay correction.
- The nonlinear v1 seam is one local central-finite-difference linearization,
  not Gauss-Newton iteration, `fs-ad` automatic differentiation, a derivative
  certificate, or a nonlinear-posterior accuracy guarantee. It does not select
  or adapt the step, compare multiple steps for truncation/roundoff stability,
  prove model smoothness, retain an AD authority, or authenticate callback code.
  Its standardized-innovation threshold is caller policy for local demotion,
  not a validity certificate, calibrated hypothesis test, outlier likelihood,
  convergence criterion, or `Color::Validated` promotion. Multi-iteration
  solvers, automatic Jacobians, typed child-resource planning, temporal
  augmentation, and calibrated full-rig replay remain separate work.
- The observation covariance gate certifies PSD using the crate's existing
  bounded interval policy, then requires a successful binary64 Cholesky factor
  for whitening. This is not a condition-number certificate, an exact
  correlated-posterior proof, or authority for covariance/calibration
  provenance. Singular PSD covariance is refused rather than pseudoinverted.
- `effective_dof` is the count of scalar readings actually passed to the update;
  it is not a statistical rank, residual degrees-of-freedom estimate, or
  likelihood-complexity certificate.
- The robust audit receipt is content-addressed replay evidence, not an
  authenticated calibration, ledger admission, or experimental-validation
  receipt. The wrapper's Cholesky/whitening/hash work has no typed child-budget
  planner yet.
- `Color::Estimated` is deliberate. Consuming calibrated measurements does not
  by itself prove model-form validity in a regime. Promotion to
  `Color::Validated` requires an external admission boundary that verifies
  calibrated dataset provenance, retained validation evidence, and an
  authenticated authority. No generic authenticated validation capability is
  currently owned by this crate, so it exposes no promotion API.
- The semidefinite admission check uses outward-rounded intervals to certify
  the signs of the particular scaled-elimination pivots it admits. This is not
  a complete PSD decision procedure, an eigenvalue enclosure, or a theorem
  about conditioning, and it can return `CovarianceCertificationUnresolved`
  for a mathematically singular PSD matrix, a sufficiently ill-conditioned
  strictly SPD matrix, or an indefinite matrix whose sign it cannot resolve.
  Public admission requires a successful certificate, not merely mathematical
  positive semidefiniteness. The
  Joseph update and other Kalman arithmetic remain floating-point and are not
  interval-certified. No ambiguous or negative pivot is tolerance-clamped into
  an admitted zero pivot.
- The dense v0 state cap is a resource contract, not a scientific limit. Larger
  states require a future sparse or matrix-free implementation with explicit
  memory budgets and tile-boundary cancellation rather than raising this cap.
- The observation operators are SUPPLIED restriction-map rows; deriving them
  from the sheaf's trace maps is fs-geom's. The registration variance accepted
  by scan observations is likewise a caller-supplied, separately calibrated
  variance. The current fs-asbuilt residual RMS is a global fit diagnostic,
  not transform covariance, and MUST NOT be passed as that variance. A typed
  fs-asbuilt covariance or spatial-uncertainty handoff remains future work.
- Ledgering the update onto the per-regime model-form posterior (Proposal 3's
  maps), authenticating calibration receipts, and admitting a validated claim
  are fs-ledger/fs-package integration work; this crate produces the estimated
  candidate and proposed validity domain.
- Logical-work totals and fixed poll strides are accounting/cancellation
  semantics, not claims about instructions, wall-clock time, allocation peaks,
  pause/resume state, deadline/cost enforcement, drain latency, or a
  200-microsecond cancellation bound. G5 is same-process replay for a fixed
  implementation/toolchain manifest, not cross-ISA bit stability.
- Typed planner byte counts are conservative semantic payload envelopes, not
  allocator-overhead or process-RSS measurements. `CostUnits` is abstract and
  is not a wall-time, currency, or energy certificate. A planner describes a
  grant but does not itself admit an invocation; the parent `fs-exec` issuer
  owns admission, deadline/capability/accuracy identities, and the terminal
  receipt.
