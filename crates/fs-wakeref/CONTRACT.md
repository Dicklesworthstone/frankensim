# fs-wakeref — CONTRACT

Layer: L3. Bead: frankensim-wf-root-guzez.5.22 (E4.9b, Wright Flyer program).
Spec: COMPREHENSIVE_PLAN_FOR_REAL_TIME_WRIGHT_FLYER_SIM_WITH_FRANKENSIM.md
(ROUND 6 steady state), V-08b1.

## Purpose and layer

Layer L3 (independent referee leaf). What this crate IS:

The INDEPENDENT unsteady prescribed-wake referee — the dense reference
leaf the E4.3b3 campaign compares the A1 FOM/ROM lane against. Its whole
value is structural independence: it depends on `fs-math` + `fs-blake3`
ONLY, shares no solver code with fs-airfoil (Duhamel indicial machinery)
or fs-wing (strip solver / prescribed-wake operator), and even the
Biot–Savart segment kernel is written here from the formula. The battery
pins the dependency closure.

## Invariants

Model tier, DECLARED (these are the referee's structural laws):

Single-chordwise-ring unsteady vortex lattice over canard + wing:
- PRESCRIBED wake — rings shed each step convect rigidly downstream at
  `convection · V`; never a free wake.
- Exact flat-ground vortex images (mirror + reversed traversal) when a
  ground case is on.
- Thin-airfoil closure; unsteady lift = Kutta–Joukowsky + half-chord
  apparent-mass term. The impulsive-start step therefore carries a
  non-circulatory spike at step 0 (physical, declared), and the
  single-ring lattice's circulatory starting deficiency is SHALLOW
  (~0.92 at s≈1 vs 0.5 for resolved 2-D Wagner) — recorded as the
  referee's own character in the receipt; A1 is judged against THIS
  receipt, never against 2-D Wagner.
- Wake memory capped (`MAX_WAKE_ROWS`, oldest dropped) — the cap is in
  the receipt; truncation is never silent.

## Public types and semantics

| Entry | Contract |
|---|---|
| `wright_geometry_v1` | registered 1903 reference geometry (flyer-reference lineage) + content digest |
| `v08b1_cases_v1` | the registered case set: impulse / step / chirp / reversal × free-air / flat-ground |
| `run_case` | one dense march → per-step wing lift, canard lift, hinge moment (40 %-chord axis, declared) + bitwise series digest |
| `emit_v08b1_receipt` | the V-08b1 receipt: per-case digests, steady values, Wagner-convention ratio (step rows sample step 3 — after the non-circulatory spike), receipt digest (golden-pinned) |

## Error model

Typed refusals (`Refusal { code, message, ranked_repairs }`), caps
tested at cap AND cap+1. Vocabulary:

`referee-case-invalid` (caps at cap AND cap+1: v [5,40], |alpha0| ≤ 0.3,
rho (0.5,2.0), dt (0,0.05], steps [1,2400], convection [0.5,1.5], ground
strictly below both surfaces); `referee-system-singular`.

## No-claim boundaries

- No real-time claim — dense and slow on purpose; never a sim-loop path.
- No historical-fidelity claim: the geometry is the registered reference,
  the closure is thin-airfoil; the receipt exists for A1-vs-referee
  DELTAS (E4.3b3), not absolute truth.
- Hinge moment is the quarter-chord lift about the 40 %-chord axis — a
  declared tier, not the E4.2b mechanism model.

## Determinism class

Bit-exact: fixed-shape marches, `fs-math det::` transcendentals only,
no RNG, no wall-clock. Every series carries a bitwise `fs-blake3`
digest; the V-08b1 receipt digest is golden-pinned (measure-then-pin,
golden-bump protocol on any change).

## Cancellation behavior

None: `run_case` is a bounded synchronous march (steps capped at 2400)
with no I/O, no threads, and no checkpoint surface. Callers budget by
choosing `n_steps`; the cap refuses, never truncates.

## Unsafe boundary

No `unsafe`. No FFI. The dependency closure (`fs-math` + `fs-blake3`
only) is pinned by the battery — structural independence IS the crate's
value and any new dependency is a contract change.

## Feature flags

None. No cargo features alter semantics.

## Conformance tests

`tests/` battery: dependency-closure pin, geometry digest stability,
case-cap refusals both sides, per-fixture series goldens, the receipt
golden (c439dbce class), and the Wagner-convention sampling rule.
Repro: `cargo test -p fs-wakeref`.
