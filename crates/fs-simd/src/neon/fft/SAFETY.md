# SAFETY: fs-simd/src/neon/fft/mod.rs

Registered in unsafe-capsules.json (split from the elementwise capsule
under the 300-line cap, bead 8nfp); enforced by `cargo run -p xtask --
check-unsafe`.

## Invariants
`r4qrun_f64`'s `vld2q`/`vst2q` accesses touch exactly 4 f64 at offset
4·q2 with 4·q2 + 4 ≤ s2 (asserted, with s2 % 4 == 0 on the vector
path); the four output rows live at disjoint offsets j·s2 within `out`
(len 4·s2, asserted). Runs whose length is not a multiple of 4 f64
delegate WHOLE to the scalar twin in safe code.

## Aliasing assumptions
Four shared input runs and one exclusive output block; the borrow
checker rules out overlap at the façade.

## Alignment assumptions
None (`vld2q_f64`/`vst2q_f64` tolerate unaligned addresses).

## Lifetime assumptions
No pointer escapes; lifetimes are the borrowed slices'.

## Panic behavior
The run-length assert fires BEFORE any unsafe block.

## Cancellation behavior
No poll points: one bounded pass over the runs.

## Concurrency behavior
No shared state, no atomics.

## Miri coverage
Miri routes to the scalar twin (dispatch layer). Compensating checks:
the tier battery runs r4qrun bitwise against the twin over even/odd
run lengths, both directions, special values; fs-fft's golden hash is
tier-invariant (it did NOT move when this capsule replaced the inline
loop).

## Model-checking coverage
N/A.

## Fuzz/property coverage
`tier_equivalence_battery` (see Miri coverage).

## Proof obligations discharged by callers
None; the façade is safe and total.
