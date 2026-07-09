# SAFETY: fs-simd/src/neon/mod.rs

> THE exemplar capsule (obligation transferred from the unsafe-safety-cases
> bead). Registered in unsafe-capsules.json; enforced by
> `cargo run -p xtask -- check-unsafe` (registration, <300 lines, this file).

## Invariants
All `unsafe` is confined to NEON load/store/arithmetic intrinsics operating
on pointers derived from `as_chunks::<2>()` fixed-size arrays (and
`.add(2)` within `as_chunks::<4>()` blocks) over live `&[f64]`/`&mut [f64]`
slices. Every access is inside the borrow-checked allocation, correctly
typed, and exactly 2 lanes wide. Tails are handled by the scalar twin in
safe code; no partial-lane intrinsic access exists.

`mk8x4_f64` (the GEMM register microkernel, bead xdgf) additionally derives
raw pointers from the panel slices with computed offsets; the leading
`assert!(a.len() >= kc*8 && b.len() >= kc*4)` bounds every offset used
(max a: kk·8+6 with a 2-lane read ⇒ kc·8; max b: kk·4+2 with a 2-lane read
⇒ kc·4), and the accumulator loads/stores are 2-lane accesses at offsets 0
and 2 of `[f64; 4]` rows whose extent the type system proves. Padded panel
tails are the CALLER's packing invariant (fs-la zero-pads); soundness here
never depends on it — only the asserted lengths matter.

`btile4x4_f64` (the batched-GEMM tile kernel, bead 9ekv) walks eight
stream pointers of the form base + l·stride + 2·p (a) / base + l·k·stride
+ 2·p (b): the leading assert bounds the maximal dereferenced offset
(ti ≤ 3, l ≤ k−1, 2p ≤ mb−2) inside both plane buffers, every access is
2 lanes, and the per-pair rewind (−k·stride +2 / −k²·stride +2) never
leaves the borrowed allocations, so pointer provenance is preserved. The
odd batch-lane tail routes through the scalar twin in safe code.

## Aliasing assumptions
Input slices are `&[f64]`, outputs `&mut [f64]`; Rust's borrow rules already
guarantee no mutable aliasing at the façade. No function both reads and
writes overlapping memory except through a single `&mut` chunk (`axpy`,
`scale`), where load-then-store to the same chunk is a plain read-modify-write
of exclusively-borrowed memory.

## Alignment assumptions
None. `vld1q_f64`/`vst1q_f64` support unaligned access architecturally; the
capsule imposes NO alignment precondition (fs-alloc's 128-byte policy is a
performance choice upstream, not a soundness requirement here).

## Lifetime assumptions
No pointers escape the loop iteration that derives them; all lifetimes are
those of the borrowed slices.

## Panic behavior
Length-mismatch `assert_eq!` fires BEFORE any unsafe block (programmer-error
contract, documented in CONTRACT.md). No unwinding can occur between an
intrinsic load and its paired store.

## Cancellation behavior
No poll points: every function is a bounded, allocation-free loop over its
input (microseconds at kernel sizes). Tile-level cancellation happens in the
callers that chunk work (fs-exec discipline).

## Concurrency behavior
No shared state, no interior mutability, no atomics. `Send`/`Sync` are the
slices' own properties. Nothing here can data-race.

## Miri coverage
Miri cannot interpret NEON intrinsics; under `cfg(miri)` the dispatch layer
routes to the scalar twin and this module is not exercised (documented
limitation). Compensating checks: the tier-equivalence property battery runs
the capsule against the scalar twin on every native test run, including
subnormals, NaN, ±0, and lengths 0..67 covering all tail shapes.

## Model-checking coverage
N/A (no concurrency — see Concurrency behavior).

## Fuzz/property coverage
`tier_equivalence_battery` in fs-simd's tests: seeded LCG inputs across
special values and every tail length; elementwise ops bitwise-equal to the
scalar twin (both fused per the FMA policy); reductions within the
documented cross-shape envelope and bit-stable per tier. `mk8x4_f64` is
battery-covered bitwise vs its twin for kc ∈ 0..17 ∪ {256} with special
values and nonzero starting accumulators (the KC-chunk fold path), and
transitively by fs-la's GEMM golden hash, which is tier-invariant.

## Proof obligations discharged by callers
None. The façade functions are safe and total for all slice inputs
(mismatched lengths panic by documented contract before any unsafe code).
