# SAFETY: fs-fmm/complex/avx2

## What is unsafe here

One `unsafe fn` (`dot_avx2`), unsafe only because `#[target_feature]`
requires it, and one `unsafe { }` call site in the dispatcher (`dot`). All of
it is compiled only for x86-64.

## Invariants

- `dot_avx2` is called only after `is_x86_feature_detected!("avx2")` returned
  true on the running CPU, immediately before the call. That is the only
  precondition `target_feature` imposes.
- The body (`dot_body` in `complex.rs`) is ordinary safe Rust: `chunks_exact`,
  slice zips, and `f64` multiplies, adds and subtracts. No intrinsics, raw
  pointers, transmutes or inline assembly.
- Bit identity: every product is a separate IEEE multiply, and every lane
  accumulates in the same order. AVX2 does not enable FMA, and Rust emits no
  contraction of a separate multiply and add, so each lane performs the same
  correctly rounded operations under both codegens. Only the vector width
  (two or four lanes per instruction) differs.

## Aliasing assumptions

N/A: safe shared slices only; the body writes nothing but local accumulators.

## Alignment assumptions

N/A: no pointer arithmetic or explicit vector loads; ordinary slices of `C64`.

## Lifetime assumptions

N/A: no lifetimes are erased; arguments are borrowed for the call only.

## Panic behavior

Identical to the portable path. The body cannot panic for any inputs: zips
stop at the shorter slice, and the remainders come from `chunks_exact`.

## Cancellation behavior

N/A: a pure function of two slices; cancellation is polled by the callers
between passes, exactly as with the portable body.

## Concurrency behavior

The function is pure and touches no shared state. The FMM's declared workers
may call it concurrently on read-only slices. `Send`/`Sync` of the facade
types are unchanged.

## Miri coverage

Miri does not model x86 target features; the capsule's only effect is codegen,
and the same body runs under Miri through the portable dispatch branch.

## Model-checking coverage

N/A: no atomics, locks or shared mutable state.

## Fuzz/property coverage

`complex/tests.rs` (`avx2_dot_capsule_is_bitwise_the_portable_body`) compares
the dispatched `dot` with the portable `dot_body` bitwise. It covers lengths
0 to 1000, including every remainder class, with mixed-magnitude random
complex inputs. Every FMM potentials test then runs through the dispatched
path on an AVX2 host.

## Proof obligations discharged by callers

None: `dot` is safe to call on any CPU.

## Blast radius

If feature detection were wrong, the CPU would fault on the first AVX2
instruction; there would be no silent corruption. Results cannot differ,
because the arithmetic is operation-for-operation identical.
