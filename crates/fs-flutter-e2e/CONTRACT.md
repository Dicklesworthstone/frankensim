# CONTRACT: fs-flutter-e2e

FlutterCert — a PROVEN fluid-structure stability boundary. Layer L4 (ASCENT).

## Purpose and layer

Composes `fs-sos` (Lyapunov certificate), `fs-spectral` (independent eigenvalue
cross-check), `fs-couple` (partitioned FSI solver), `fs-evidence` (Verified).
Deps point downward.

## Public types and semantics

- `operator(mu) -> [[f64;2];2]` — the 2-DOF coupled operator `A(μ)`.
- `spectral_abscissa(mu) -> f64` — largest eigenvalue of the symmetric part.
- `run_campaign(lo, hi, steps) -> FlutterReport` — sweeps `μ`, records the
  Lyapunov certificate, the independent spectral abscissa, and naive/Aitken
  partitioned-solver convergence at each point; reports the proven boundary, the
  cross-check agreement, the solver boundaries, and a witness.

## Invariants

- The system is asymptotically stable iff `μ < 2`; `lyapunov_certifies_stability
  (A(μ), I)` recovers exactly that boundary → the certificate is `Verified`.
- CERTIFYING THE CERTIFIER: the independent spectral abscissa `−1 + μ/2 < 0`
  gives the SAME boundary (`boundaries_agree`).
- The naive partitioned solver quits early; Aitken relaxation converges strictly
  past it (`aitken_beats_naive`), giving a witness `μ` that is certified stable,
  beyond the naive solver, and Aitken-computable.
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

`tests/flutter.rs` (3): the boundary is proven and cross-checked (μ*≈2, naive
quits ~0.85, Aitken reaches 2.45, witness ≈0.95 Verified); beyond μ*=2 nothing is
certified; determinism.

## No-claim boundaries

The model is a 2-DOF linear operator and a scalar added-mass interface map — a
minimal flutter/added-mass surrogate, not a full aeroelastic model; the Lyapunov
certificate uses `P=I` (a sufficient witness), not a synthesized `P`.
