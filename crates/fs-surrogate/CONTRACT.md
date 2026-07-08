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

None.

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
