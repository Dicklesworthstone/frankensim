# SAFETY: fs-feec highorder fma capsule

## What is unsafe here

One `unsafe fn apply_mono_x86` (only because `#[target_feature]`
requires it) and one `unsafe { }` call site in the dispatcher,
compiled only for x86-64.

## Why it is sound

- The dispatcher verifies `avx2` and `fma` via
  `is_x86_feature_detected!` immediately before the call — the only
  precondition `target_feature` imposes.
- The function body is `TensorSpace::apply_mono_body`, ordinary safe
  Rust (slice indexing, `mul_add`), marked `#[inline(always)]` so its
  code is GENERATED inside this function under the enabled features;
  no intrinsics, no pointers, no transmutes anywhere.
- Bit-identity: `f64::mul_add` is a single IEEE-754 fused operation
  under both codegens (libm call on baseline, `vfmadd` here) — same
  value per element, gated by the crate's sf-kron golden hash.

## Blast radius

None beyond the capsule contract: if feature detection were wrong the
CPU would fault on the first FMA instruction (no silent corruption);
results cannot differ because the arithmetic is operation-for-
operation identical.
