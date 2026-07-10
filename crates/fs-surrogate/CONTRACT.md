# CONTRACT: fs-surrogate

Learned accelerators with guarantees: surrogates permitted only inside certified
validity bands — ML proposes, certified numerics disposes.

## Purpose and layer

Layer L4 (surrogate / ROM). No dependencies — pure Rust (an in-house symmetric
eigensolver for the method of snapshots).

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
- `SurrogateError` — `NoSnapshots` / `DimMismatch` / `BadThreshold`.

- `ladder` module (addendum Proposal A, bead knh1.4; [F], behind
  `abstraction-ladder`): the certified abstraction ladder.
  `TruthModel` (the P1 full-order elliptic family — level 0's DECLARED
  semantics: "truth" means the FE model, discretization honesty stated
  here rather than hidden in a bound), `RbLevel` (offline snapshots +
  energy-orthonormal basis; online k×k Galerkin with the textbook
  a-posteriori bound ‖u−u_rb‖_a ≤ ‖r‖_{V′}/√α_LB via the exact Riesz
  representer and the exact affine coercivity floor; compliance QoI
  bound = energy bound SQUARED by Galerkin symmetry), `ConceptLevel`
  (interpolation lookup, ESTIMATED color, dispersion calibrated by
  cross-rung discrepancy probes), `Ladder`/`at_level(k).query(μ, tol)`
  (AUTOMATIC CERTIFIED DESCENT: a rung answers only when its
  certificate meets the tolerance; leaks are recorded and descended
  past — invisible until it leaks), `rb_coverage` (the kill
  measurement).

## Invariants

- POD reproduces an exactly-representable (low-rank) snapshot set to roundoff;
  its modes are orthonormal; the retained rank captures `>= energy_threshold`.
- The conformal band achieves at least its nominal `(1−α)` empirical coverage on
  exchangeable held-out data.
- `certify_or_escalate` uses the surrogate ONLY when trustworthy (in-domain +
  band tight enough), so a fleet of queries costs strictly less than
  all-high-fidelity whenever any query is served by the surrogate.

## Error model

Structured `SurrogateError`; the only panics are nonsensical conformal inputs
(empty residuals, `α ∉ (0,1)`).

## Determinism class

Fully deterministic: the eigensolver, POD, band, and policy are pure functions
of their inputs.

## Cancellation behavior

None (synchronous pure functions).

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

- `abstraction-ladder` [F] (default OFF) — the certified abstraction
  ladder (knh1.4, Proposal A; `dep:fs-evidence`); gates the `ladder`
  integration target.

## Conformance tests

`tests/surrogate.rs` (8 cases): POD reproduces a low-rank set exactly;
orthonormal modes; energy-based rank + reduced error; bad-input rejection; the
conformal band achieves nominal coverage; certify-or-escalate uses the surrogate
only when trustworthy; the policy reduces cost vs all-high-fidelity;
determinism.

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
- Level 0's bound is zero BY DECLARATION (the FE model is the truth
  semantics); the FE discretization error is a separate ledger entry,
  not this module's claim.
- The concept rung's dispersion is a probe MAXIMUM, not a bound — the
  Estimated color is load-bearing.
- Per-REGION (spatial) RB decomposition and the fs-ir at_level query
  integration are the named growth seams.
