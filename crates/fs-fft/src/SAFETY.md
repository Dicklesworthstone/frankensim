# SAFETY: fs-fft/src/simd_view.rs

Registered in unsafe-capsules.json; enforced by `cargo run -p xtask --
check-unsafe` (registration, <300 lines, this file).

## Invariants
Two functions, each a single slice reinterpretation from `&[C64]` /
`&mut [C64]` to `&[f64]` / `&mut [f64]` of twice the length. `C64` is
`#[repr(C)] { re: f64, im: f64 }`: size 16, align 8, no padding, so the
byte range of the source slice is exactly the byte range of the result
and every byte pattern is a valid f64.

## Aliasing assumptions
The borrow (shared or exclusive) and its lifetime transfer to the view
unchanged; no new aliasing is possible.

## Alignment assumptions
f64 alignment (8) is implied by C64 alignment (8).

## Lifetime assumptions
The views borrow the input slices; nothing escapes.

## Panic behavior
None (no arithmetic, no indexing).

## Cancellation behavior
No poll points: O(1) pointer casts.

## Concurrency behavior
No shared state; the views inherit the slices' Send/Sync.

## Miri coverage
Fully Miri-checkable (no intrinsics); exercised by every fs-fft test
through the transform path.

## Fuzz/property coverage
Transitively by the fs-fft oracle battery and golden hash (the views
feed fs-simd's r4qrun kernel, whose tier battery is bitwise), and the
golden did NOT move when the capsule path replaced the inline loop.

## Proof obligations discharged by callers
None; total for all slices.
