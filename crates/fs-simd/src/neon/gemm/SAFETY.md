# SAFETY: fs-simd/src/neon/gemm/mod.rs

Registered in unsafe-capsules.json (split from the elementwise capsule
under the 300-line cap, bead 8nfp); enforced by `cargo run -p xtask --
check-unsafe`.

## Invariants
`mk8x4_f64` derives raw pointers from the packed panel slices with
computed offsets; the leading `assert!(a.len() >= kc*8 && b.len() >=
kc*4)` bounds every offset used (max a: kk·8+6 with a 2-lane read ⇒
kc·8; max b: kk·4+2 with a 2-lane read ⇒ kc·4), and the accumulator
loads/stores are 2-lane accesses at offsets 0 and 2 of `[f64; 4]` rows
whose extent the type system proves. Padded panel tails are the
CALLER's packing invariant (fs-la zero-pads); soundness never depends
on it — only the asserted lengths matter.

`btile4x4_f64` walks eight stream pointers of the form base + l·stride
+ 2·p (a) / base + l·k·stride + 2·p (b): the leading assert bounds the
maximal dereferenced offset (ti ≤ 3, l ≤ k−1, 2p ≤ mb−2) inside both
plane buffers, every access is 2 lanes, and the per-pair rewind
(−k·stride +2 / −k²·stride +2) never leaves the borrowed allocations,
so pointer provenance is preserved. The odd batch-lane tail routes
through the scalar twin in safe code.

## Aliasing assumptions
Inputs are `&[f64]`, outputs `&mut` slices/arrays; the borrow checker
guarantees no mutable aliasing at the façade.

## Alignment assumptions
None. `vld1q_f64`/`vst1q_f64` support unaligned access architecturally.

## Lifetime assumptions
No pointer escapes its function; all lifetimes are the borrowed slices'.

## Panic behavior
Bounds asserts fire BEFORE any unsafe block. No unwinding between an
intrinsic load and its paired store.

## Cancellation behavior
No poll points: bounded, allocation-free kernels (microseconds).

## Concurrency behavior
No shared state, no atomics; `Send`/`Sync` are the slices' properties.

## Miri coverage
Miri cannot interpret NEON intrinsics; the dispatch layer routes to
scalar under Miri. Compensating checks: the tier-equivalence battery
runs both kernels bitwise against their scalar twins on every native
test run (kc ∈ 0..17 ∪ {256} with special values and nonzero starting
accumulators for mk8x4; k/mb/offset sweeps for btile4x4), and fs-la's
GEMM/batched goldens are tier-invariant.

## Model-checking coverage
N/A (no concurrency).

## Fuzz/property coverage
`tier_equivalence_battery` in fs-simd's tests (see Miri coverage).

## Proof obligations discharged by callers
None; the façades are safe and total (length violations panic before
unsafe code).
