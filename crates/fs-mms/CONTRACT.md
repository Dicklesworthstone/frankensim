# CONTRACT: fs-mms

The G1 harness (bead frankensim-epic-gauntlet-6nb.2): manufactured-
solution refinement ladders, deterministic convergence-order fitting,
the 0.2 slope gate for primal and adjoint orders, and the declared
lintable MMS battery matrix.

## Purpose and layer

Layer L1 (depends on fs-math for `det::ln`, keeping order fits
bit-identical across ISAs). Discretizations at any layer feed it
(h, error) ladders; it fits, gates, and records.

## Public types and semantics

- `RefinementLadder::new(hs, errors)` — ≥ 3 rungs, strictly decreasing
  positive `h`, positive finite errors; accessors `hs()`/`errors()`.
- `fit_order(&ladder) -> OrderFit { observed, intercept, rms_residual }`
  — least-squares slope of log(error) vs log(h) via `fs_math::det::ln`.
- `OrderGate { theoretical }` and `ORDER_GATE_TOLERANCE = 0.2` —
  `check(case, side, &ladder)` returns the passing `OrderVerdict`
  (JSON-line renderable) or the `mms-order-gate` refusal carrying the
  failing record; deviation ABOVE theoretical fails too
  (superconvergence claims need their own theory).
- `LadderSide::{Primal, Adjoint}` — identical gate, visible distinction
  (dual consistency is verified, not assumed).
- `MmsMatrix`/`MmsMatrixRow`/`Coverage::{Covered, Gap}` — the battery
  matrix declared in data; `gaps()` is the lint; `json_lines()` renders
  coverage.

## Invariants

- Order fits are pure and deterministic: `det::ln` logs, fixed
  summation order, no RNG — bit-identical across ISAs and runs.
- A ladder needs at least three rungs; two-point slopes are refused
  (`mms-ladder-shape`), as are non-monotone h and non-finite or
  non-positive values.
- The gate is two-sided: |observed − theoretical| > 0.2 fails in either
  direction, and every refusal carries the full failing record.
- Matrix rows are either covered by a NAMED test or an explicit
  reasoned gap; there is no third, silent state.

## Error model

Typed `MmsError { rule, detail }` with stable rule slugs
(`mms-ladder-shape`, `mms-ladder-order`, `mms-ladder-domain`,
`mms-gate-domain`, `mms-order-gate`). No panics on the public surface.

## Determinism class

Deterministic (cross-ISA bitwise): the only transcendental is
`fs_math::det::ln`; everything else is fixed-order IEEE arithmetic and
`sqrt` (correctly rounded by IEEE-754).

## Cancellation behavior

None: fitting and gating are synchronous over caller-materialized
ladders; running discretizations to PRODUCE ladders is the caller's Cx
scope.

## Unsafe boundary

No `unsafe` (workspace forbids it; nothing here needs it).

## Feature flags

None.

## Conformance tests

`tests/mms.rs` — exact slope recovery on pure power laws (1e-10),
two-sided gate behavior at 1.85/1.75/2.31 vs 2.0, adjoint-side records,
every named ladder refusal, a REAL 1-D finite-difference Poisson MMS
(u* = sin(πx), Thomas solve, four-rung ladder) gating green at second
order, and the declared battery matrix with lintable gaps.

## No-claim boundaries

- The harness fits and gates orders it is HANDED: it does not
  discretize, mesh, or solve.
- No symbolic forcing-term generation: manufactured-solution synthesis
  from operator definitions is the fs-opdsl integration (tfz.4),
  declared as an explicit matrix gap until it lands.
- No result ledgering or machine-fingerprint binding (fs-obs/fs-ledger
  scope); verdicts are JSON lines for the caller's log.
- The 0.2 tolerance gates asymptotic-range ladders; it cannot detect a
  ladder taken entirely in the pre-asymptotic regime (the rms_residual
  envelope is reported so callers can see bent ladders).
- Migrating existing per-crate slope tests (fs-feec, fs-cutfem, fs-iga)
  onto this harness is adoptive per-crate work, tracked by the matrix.
