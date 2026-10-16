# SAFETY — fs-cheb x86 FMA-codegen capsule

Scope: `crates/fs-cheb/src/fma/mod.rs` (`cheb_eval_x86`,
`dsq_into_x86`, `matvec_into_x86`, `os_matmul_x86` and their
dispatchers), compiled only on `target_arch = "x86_64"`.

## Invariants the `unsafe` calls rely on

1. **Feature availability.** Every `target_feature(enable = "avx2,fma")`
   function is called exactly once per dispatch, immediately after
   `is_x86_feature_detected!("avx2") && is_x86_feature_detected!("fma")`
   verified the CPU supports it. No other call sites exist (the fns are
   private to this module).
2. **The bodies are safe code.** `Cheb1::eval_body`, `dsq_into_body`,
   `matvec_into_body`, and `orr_sommerfeld::matmul_body` are
   `#[inline(always)]` safe slice arithmetic; the `unsafe` here is ONLY
   the `target_feature` calling contract, never memory access.
3. **Bit identity.** `mul_add` compiles to the native fused instruction
   under the enabled feature and to libm `fma()` without it — both are
   correctly-rounded single-rounding IEEE fused ops, so results are
   bit-identical; the crate's batteries and the workspace goldens gate
   this. Every chain shape (Clenshaw backward recurrence, j-inner dot
   chains, k-outer skip-zero saxpy) is untouched — LLVM does not
   re-associate strict FP reductions, so enabling the feature cannot
   reorder them.

## Twin

The portable path IS the twin: the same `inline(always)` body compiled
without the feature. The dispatcher falls back to it on every non-FMA
x86 host and on every non-x86 target, so the capsule can never be the
only implementation of anything.
