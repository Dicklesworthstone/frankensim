# CONTRACT: fs-spectral

Spectral health monitoring (plan addendum, Proposal 5): the sheaf-Laplacian
λ-gap as a runtime health metric, with low-confidence propagation and a router
conditioning term. Owns risk R5.

## Purpose and layer

Layer L1. Depends only on `fs-evidence` (UTIL, the `Color` lattice). Pure
deterministic numerics + epistemics.

## Public types and semantics

- `symmetric_eigenvalues(&[Vec<f64>]) -> Result<Vec<f64>, SpectralError>` —
  ascending eigenvalues of a small dense symmetric matrix (in-house cyclic
  Jacobi); rejects empty/non-square/non-symmetric.
- `spectral_gap(&[f64]) -> Option<SpectralGap>` — the Fiedler gap above the
  smallest eigenvalue and its dimensionless `ratio = gap / spread` in `[0, 1]`.
- `GapHealthMonitor::new(degrade_below, restore_above)` + `update(ratio) ->
  Health` / `health()` — a HYSTERESIS classifier (`Healthy` / `Degraded`); a
  region degrades below the lower threshold and recovers only above the upper
  one, so it cannot flap. `new` panics on an inverted band.
- `propagate(color, health) -> Color` — a `Degraded` gap DEMOTES a
  verified/validated color to estimated and NEVER promotes; `Healthy` passes
  the color through.
- `compose_conditioning(&[f64]) -> Result<f64, SpectralError>` — per-op
  amplification factors multiply into an end-to-end conditioning estimate
  (empty → 1.0); rejects a negative/non-finite factor.
- `route(&[RouterPath], conditioning_weight) -> Option<&RouterPath>` — picks the
  path minimizing `base_cost + weight · ln(max(conditioning, 1))`.

## Invariants

- `propagate` is monotone-down under `Degraded`: confidence can only fall, so
  merge/triage in a degraded region always surfaces low confidence.
- The health monitor never flaps (`degrade_below < restore_above`).
- Conditioning composes multiplicatively; a well-conditioned pipeline is `1.0`.
- The eigensolver is deterministic and validates symmetry.

## Error model

Structured `SpectralError`; the only panic is a programmer error (inverted
hysteresis band).

## Determinism class

Fully deterministic: eigenvalues, gap, health transitions, and routing are pure
functions of their inputs.

## Cancellation behavior

None (synchronous pure functions).

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/spectral.rs` (Proposal 5, 9 cases): Jacobi recovers known spectra (2×2,
antisymmetric-off-diagonal, diagonal, tridiagonal Toeplitz) and rejects
malformed matrices; the gap ratio reflects separation and collapse; hysteresis
health + inverted-band panic; low-confidence propagation (demote, never
promote); multiplicative conditioning + bad-factor rejection; the router's
cheapness-vs-conditioning trade; determinism.

## No-claim boundaries

- The eigensolver is a small DENSE Jacobi for the monitoring path; production
  tracks the gap on SPARSE Laplacians with a few Lanczos/LOBPCG vectors,
  WARM-STARTED across edits (fs-la) — this crate provides the health/gap/
  conditioning logic, not the production sparse solver.
- `propagate` implements the color DEMOTION; wiring it into Proposal 10's merge
  adjudication outputs is the merge crate's integration.
- The router conditioning TERM is provided here; the full Rep Router fitness
  (cheapness + conditioning + other terms) lives in the router.
- Per-op amplification factors are supplied by each op's error model; this
  crate composes them, it does not estimate them.
