# SAFETY: fs-lbm/d3q19/fma

## What is unsafe here

One `unsafe fn`, `collide_cell3_fma`. It is unsafe only because
`#[target_feature]` requires it. There is one `unsafe { }` call site, in the
private safe thunk `fma_checked`. All of it is compiled only for x86-64.

## Invariants

- **Precondition:** `collide_cell3_fma` is called only from `fma_checked`, and
  the selector `collision_kernel3` is the only code that names `fma_checked`.
  It does so only after `is_x86_feature_detected!("fma")` returned true on the
  running CPU. The result is cached for the process, because a CPU's features
  don't change while a process runs. That is the only precondition
  `target_feature` imposes.
- **Dispatch runs once, not per cell.** Hot loops call `collision_kernel3`
  before the loop and hold the returned kernel. `fma_checked` is private, so a
  caller can't reach it except through that selector.
- **The body is ordinary safe Rust.** `collide_cell3_with_projection`, and the
  moment helpers it inlines in `d3q19.rs`, use fixed-size array indexing,
  `f64::mul_add` and plain arithmetic. There are no intrinsics, raw pointers,
  transmutes or inline assembly.
- **Bit identity:**
  - `f64::mul_add` is one correctly rounded IEEE-754 fused operation under both
    codegens: a library call on baseline x86-64, and `vfmadd` here.
  - Rust never contracts a separate multiply and add into a fused one, and it
    never reassociates floating-point reductions.
  - So every value is identical, operation by operation. The moment-space
    helpers are `#[inline(always)]` so that they generate inside this function.
    A helper that failed to inline would only lose speed, never change a value.

## Aliasing assumptions

N/A: arguments and results are passed by value (`[f64; 19]`, `[f64; 3]`, a
`Copy` model). There are no references.

## Alignment assumptions

N/A: there is no pointer arithmetic and no explicit vector load. Stack arrays
have their natural alignment.

## Lifetime assumptions

N/A: no references and no erased lifetimes.

## Panic behavior

Identical to the portable path.
- Admissibility failures return `CollisionError3`, the same variant on both
  paths (pinned).
- The only panic is the internal `unreachable!` for a BGK model routed to the
  moment helpers, which the model match makes impossible on both paths.
- No state outside the returned value exists, so an unwind cannot leave a
  partial result.

## Cancellation behavior

N/A: one cell collision is a pure function with no poll points. Callers poll
at their existing tile and step boundaries, unchanged.

## Concurrency behavior

The function is pure and reentrant. Pooled grid steps call it from many worker
threads at once, each through the kernel it resolved before its loop. The
selector keeps that kernel in a `std::sync::OnceLock`, so concurrent first calls
run one initializer and every reader sees its result.
- Test-only: a `PORTABLE_ONLY` atomic flag, read with SeqCst ordering by the
  selector before the `OnceLock`, returns the portable kernel. It is compiled
  only under `cfg(test)`, and the pins that set it are serialized by a mutex. It
  takes effect at the next selector call, which a grid makes each step.
- `Send`/`Sync` of the public types are unchanged.

## Miri coverage

Miri does not model x86 target features. The capsule's only effect is codegen,
and the same body runs under Miri through the portable branch.

## Model-checking coverage

N/A (single-threaded per call): the only shared state in production is the
write-once `OnceLock` holding the selected kernel, and its ordering is std's.

## Fuzz/property coverage

`fma/tests.rs` pins two things bitwise on an x86-64 FMA host:
- **Single cells:** the dispatcher against the portable-codegen body, over
  4000 pseudo-random states per model. The models are BGK (tau 0.503 and 0.8),
  central moments, and reduced cumulants (two rate sets, one of them the gate-1
  jet's), forced and unforced. Refusals (negative density, NaN population) are
  pinned too.
- **A grid run:** every population bit after 60 `BoundaryGrid3` steps, with
  the run forced down each path in turn. The grid is 16 x 16 x 4 with a solid
  cylinder, in a forced periodic case and an open velocity/pressure duct case.

The crate's existing collision tests (moment inverse, forcing, viscous rate)
run through the dispatcher unchanged.

## Proof obligations discharged by callers

None: `collision_kernel3`, the kernel it returns, and `collide_cell3` are safe
to call on any CPU.

## Blast radius

If feature detection were wrong, the CPU would fault on the first FMA
instruction; there would be no silent corruption. Results cannot differ,
because the arithmetic is identical operation for operation.
