# SAFETY — fs-sparse x86 FMA-codegen capsule

Scope: `crates/fs-sparse/src/fma/mod.rs` (`spmv_x86`, `spmm_x86` and
the format/product dispatchers, and the direct LDLT solve dispatcher),
compiled only on `target_arch = "x86_64"`.

## Invariants the `unsafe` calls rely on

1. **Feature availability.** Every `target_feature(enable = "avx2,fma")`
   function is called exactly once per dispatch, immediately after
   `is_x86_feature_detected!("avx2") && is_x86_feature_detected!("fma")`
   verified the CPU supports it. No other call sites exist (the fns are
   private to this module).
2. **The bodies are safe code.** the format/product bodies and `LdltFactor::solve_body` are
   `#[inline(always)]` safe slice arithmetic; the `unsafe` here is ONLY
   the `target_feature` calling contract, never memory access.
3. **Bit identity.** `mul_add` compiles to the native fused instruction
   under the enabled feature and to libm `fma()` without it — both are
   correctly-rounded single-rounding IEEE fused ops, so results are
   bit-identical; the crate's cross-format bitwise suites gate this.
   The reduction shape (ascending-column chain from +0.0) is untouched.

## Twin

The portable path IS the twin: the same `inline(always)` body compiled
without the feature. The dispatcher falls back to it on every non-FMA
host and on every non-x86 architecture.

## Direct LDLT solve

The direct solve uses the same block/panel order and fused reductions in both
arms. Pivot choices, permutations, diagonal 1x1/2x2 solves, and gather/scatter
operations are unchanged. The existing direct-solve residual fixtures compare
the dispatched result bitwise with a non-inlined portable caller, including
indefinite and delayed-pivot cases. A separate RHS-scale case covers signed
zeros and small/large finite values. This extends code generation only; it
does not assert a new cross-ISA qualification or make the direct solver the
instrument engine's default.
