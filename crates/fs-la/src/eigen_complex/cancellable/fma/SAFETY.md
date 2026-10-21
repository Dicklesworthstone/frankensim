# SAFETY: fs-la/eigen_complex/cancellable/fma

## What is unsafe here

Two `unsafe fn` (`factor_fma`, `solve_fma`), unsafe only because
`#[target_feature]` requires it, and one `unsafe { }` call site in each
dispatcher. All of it is compiled only for x86-64.

## Invariants

- A `*_fma` function is called only after `is_x86_feature_detected!("fma")`
  returned true on the running CPU, immediately before the call. That is the
  only precondition `target_feature` imposes.
- The bodies (`factor_body`, `solve_body` in `cancellable.rs`) are ordinary
  safe Rust: slice indexing, `split_at_mut`, `f64::mul_add`, `C64` division and
  magnitude. No intrinsics, raw pointers, transmutes or inline assembly.
- Bit identity: every complex product is `a.re.mul_add(b.re, -(a.im * b.im))`,
  `a.re.mul_add(b.im, a.im * b.re)`, exactly `C64`'s `Mul`. `f64::mul_add` is one
  correctly rounded IEEE-754 fused operation under both codegens (a library
  call on baseline x86-64, `vfmadd` here), so each value is identical.

## Aliasing assumptions

N/A: safe references only. The elimination borrows the pivot row and the
target row through `split_at_mut`, so the borrow checker proves exclusivity.

## Alignment assumptions

N/A: no pointer arithmetic or vector loads; ordinary slices of `C64`.

## Lifetime assumptions

N/A: no lifetimes are erased; arguments are borrowed for the call only.

## Panic behavior

Identical to the portable path: indexing panics or a panicking cancellation
callback unwind out of the factor/solve before any factor or solution is
published (`lu_complex_with_cancel` returns a fresh factor only on success;
`solve_with_cancel` copies into `b` only after the body returns).

## Cancellation behavior

The bodies contain the original poll points (each pivot, each elimination
row, each substitution row); the callback runs inside the FMA-compiled code
exactly as in the portable code. The crate's cancellation tests count polls
and check atomic publication at every checkpoint on both paths.

## Concurrency behavior

N/A (single-threaded): no shared state; `Send`/`Sync` of the facade types are
unchanged.

## Miri coverage

Miri does not model x86 target features; the capsule's only effect is codegen,
and the same bodies run under Miri through the portable dispatch branch.

## Model-checking coverage

N/A (single-threaded, no atomics).

## Fuzz/property coverage

`eigen_complex/cancellable/tests.rs` pins factor data, permutation and
solutions BITWISE against the frozen original operation graph on structured,
random (n = 64, 200) and BEM-like (e^{ikr}/r, n = 150) complex matrices. On an
x86-64 FMA host that compares the capsule with baseline codegen directly.

## Proof obligations discharged by callers

None: `factor` and `solve` are safe to call on any CPU.

## Blast radius

If feature detection were wrong the CPU would fault on the first FMA
instruction (no silent corruption); results cannot differ because the
arithmetic is operation-for-operation identical.
