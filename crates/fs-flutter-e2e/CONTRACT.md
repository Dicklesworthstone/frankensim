# CONTRACT: fs-flutter-e2e

FlutterCert — a PROVEN fluid-structure stability boundary. Layer L4 (ASCENT).

## Purpose and layer

Composes `fs-sos` (Lyapunov certificate), `fs-spectral` (independent eigenvalue
cross-check), `fs-couple` (partitioned FSI solver), `fs-ivl` (outward-rounded
enclosure of the published decay rate), `fs-evidence` (Verified).
Deps point downward.

## Public types and semantics

- `operator(mu) -> [[f64;2];2]` — the 2-DOF coupled operator `A(μ)`.
- `numerical_abscissa(mu) -> f64` — largest eigenvalue of the symmetric part
  `(Aᵀ+A)/2` via `fs-spectral`; this IS the Lyapunov `P=I` condition, so it is an
  implementation cross-check, not an independent criterion.
- `spectral_abscissa(mu) -> f64` — the largest real part of `A(μ)`'s ACTUAL
  eigenvalues, `−1 + √(max(μ−1, 0))`, in round-to-nearest `f64`. A DIAGNOSTIC:
  the value can land up to about one ulp on either side of the exact abscissa,
  so it carries no bound authority.
- `spectral_abscissa_interval(mu) -> fs_ivl::Interval` — the same quantity as a
  CERTIFIED outward-rounded enclosure (the `μ−1` subtraction included, so the
  enclosure covers both the ideal operator and the `f64` matrix `operator`
  builds). Non-finite `μ` yields `Interval::WHOLE`, the no-claim answer.
- `run_campaign(lo, hi, steps) -> FlutterReport` — sweeps `μ`, records the
  Lyapunov certificate, the independent spectral abscissa, and naive/Aitken
  partitioned-solver convergence at each point. The two `*_boundary` fields are
  largest sampled stable points, not standalone boundary locations.
  `boundary_bracket` exists only when an increasing sweep witnesses the same
  stable-to-unstable adjacent sample pair under both criteria and their
  classifications agree at every sample; `boundaries_agree` is true exactly
  when that bracket exists.

## Invariants

- The system is asymptotically stable iff `μ < 2`; `lyapunov_certifies_stability
  (A(μ), I)` recovers exactly that boundary → the certificate is `Verified`.
- INDEPENDENT CROSS-CHECK: the Lyapunov `P=I` proof is only SUFFICIENT (it equals
  `−1+μ/2 < 0`). The necessary-and-sufficient eigenvalue criterion `−1+√(μ−1) < 0`
  is a genuinely different function of `μ` that reaches the same boundary, so the
  certificate is TIGHT. The sampled campaign reports `boundaries_agree` only
  when both criteria actually cross within the same retained bracket; equal
  maxima from a sweep truncated to one side cannot mint that claim. Separately,
  `fs-spectral`'s numerical abscissa recomputes the Lyapunov condition, so its
  per-sample agreement with `fs-sos` (`impl_consistent`) is an implementation
  cross-check.
- The naive partitioned solver quits early; Aitken relaxation converges strictly
  past it (`aitken_beats_naive`), giving a witness `μ` that is certified stable,
  beyond the naive solver, and Aitken-computable.
- NAMED, OUTWARD-ROUNDED WITNESS CLAIM: `witness_decay_rate_color` is
  `Verified{lo, hi}` over exactly one quantity — the LARGEST eigenvalue real part
  of `A(witness_mu)`, i.e. the asymptotic decay rate — and both endpoints are
  `spectral_abscissa_interval`'s, never the round-to-nearest `spectral_abscissa`.
  It is NOT an enclosure of the operator's spectrum: for `μ > 1` the second
  eigenvalue's real part `−1 − √(μ−1)` lies strictly below `lo` and is outside
  the claim. No color is minted when the enclosure is not finite.
- Deterministic (a fixed `μ` grid; no RNG).

## Error model

Panics only when `steps < 2`.

## Determinism class

Fully deterministic (G5).

## Cancellation behavior

None (a synchronous batch).

## Unsafe boundary

None; `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/flutter.rs` (6): the boundary is proven and cross-checked (μ*≈2, naive
quits ~0.85, Aitken reaches 2.45, witness ≈1.05 Verified); beyond μ*=2 nothing is
certified; equal co-truncated maxima and a descending sweep do not claim a
boundary; the witness decay-rate enclosure is checked against an `fs_math::dd`
(~106-bit) oracle at a reachable `μ = 1.3` where the round-to-nearest evaluation
overshoots the exact abscissa, and the second eigenvalue is confirmed outside the
named claim; the enclosure is sound on both branches and mints no bound for a
non-finite `μ`; determinism.

## No-claim boundaries

The model is a 2-DOF linear operator and a scalar added-mass interface map — a
minimal flutter/added-mass surrogate, not a full aeroelastic model; the Lyapunov
certificate uses `P=I` (a sufficient witness), not a synthesized `P`.
The exact `μ*=2` result comes from the analytic operator identities. A sampled
campaign only localizes that boundary to its retained adjacent-sample bracket;
it makes no location claim when the sweep stays on one side or runs in reverse.
The decay-rate enclosure is a claim about the LARGEST eigenvalue real part at one
sampled `μ` only. It bounds neither the whole spectrum, nor the transient
response, nor any decay rate away from that sample, and the `lyapunov_stable`
flag it accompanies remains `fs-sos`'s numerically thresholded verdict rather
than an outward-rounded one.
