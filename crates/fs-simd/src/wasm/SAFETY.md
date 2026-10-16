# SAFETY: fs-simd/src/wasm/mod.rs

> WASM32 SIMD128 Tier-1w capsule (bead frankensim-wf-root-guzez.1.5, E0.5).
> Registered unsafe capsule; enforced by `cargo run -p xtask -- check-unsafe`.

## Invariants
All `unsafe` is confined to WebAssembly 128-bit vector (`v128`) intrinsics in `core::arch::wasm32`
operating on pointers derived from `as_chunks::<2>()` fixed-size arrays over live `&[f64]`/`&mut [f64]`
slices. Every access is within bounds, correctly typed, and exactly 2 lanes wide. Tails are handled
by the scalar twin in safe code.

## Aliasing assumptions
Input slices are `&[f64]`, outputs `&mut [f64]`; Rust borrow rules guarantee exclusive mutable access.

## Alignment assumptions
WASM SIMD `v128.load` and `v128.store` support unaligned memory access in the WebAssembly specification.

## Lifetime assumptions
No raw pointers escape the local loop scope.

## Panic behavior
Pre-condition length checks fire before unsafe blocks.

## Miri coverage
Under Miri, dispatch routes to scalar twins.

## Equivalence coverage
Tested against scalar reference outputs on Biot-Savart, BEMT, AXPY, scale, mul_elem, FMA, dot, and sum.
