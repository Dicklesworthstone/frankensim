# SAFETY: fs-simd/src/neon/gemmf32/mod.rs

Registered in unsafe-capsules.json; enforced by `cargo run -p xtask --
check-unsafe` (registration, <300 lines, this file).

## Invariants
`btile4x4pf32` walks eight stream pointers of the form base(t) + l·mb
+ 4·q over PACKED l-contiguous operands (A i-major at ((i0+t)·k)·mb, B
j-major at ((j0+t)·k)·mb): the leading assert bounds the maximal
dereferenced offset (t ≤ 3, l ≤ k−1, 4q ≤ mb−4) inside both packed
buffers, every access is exactly 4 f32, and the per-quad rewind
(−k·mb + 4) never leaves the borrowed allocations (provenance
preserved). Lane counts not divisible by 4 take the scalar twin whole
in safe code.

## Aliasing assumptions
Two shared input slices, one exclusive output; the borrow checker
rules out overlap at the façade.

## Alignment assumptions
None (`vld1q_f32`/`vst1q_f32` tolerate unaligned addresses).

## Lifetime assumptions
No pointer escapes; lifetimes are the borrowed slices'.

## Panic behavior
The bounds assert fires BEFORE any unsafe block.

## Cancellation behavior
No poll points: one bounded pass per lane quad.

## Concurrency behavior
No shared state, no atomics.

## Miri coverage
Miri routes to the scalar twin (dispatch layer). Compensating checks:
the tier-equivalence battery runs this kernel bitwise against its twin
over lane counts covering the quad path and the twin-delegation tail,
with offset tiles and special values.

## Model-checking coverage
N/A.

## Fuzz/property coverage
`tier_equivalence_battery` (see Miri coverage).

## Proof obligations discharged by callers
None; the façade is safe and total.
