# SAFETY — fs-sparse x86 FMA-codegen capsule

Scope: `crates/fs-sparse/src/fma/mod.rs` (`spmv_x86`, `spmm_x86` and
the format/product dispatchers, the direct LDLT solve dispatchers,
and the positive-factor whitening dispatchers),
compiled only on `target_arch = "x86_64"`.

## Invariants the `unsafe` calls rely on

1. **Feature availability.** Every `target_feature(enable = "avx2,fma")`
   function is called exactly once per dispatch, immediately after
   `is_x86_feature_detected!("avx2") && is_x86_feature_detected!("fma")`
   verified the CPU supports it. No other call sites exist (the fns are
   private to this module).
2. **The bodies are safe code.** the format/product bodies,
   `LdltFactor::solve_body`, `LdltFactor::solve_in_place_body`, and
   `PositiveWhitening::{apply_body,transpose_body}` are
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

## Cooperative sparse products

The `checkpointed_dispatch` entries for CSR, BSR, SELL and compact CSR use
the same AVX2/FMA feature check and safe, inlined arithmetic bodies. Their
fallible callbacks are ordinary safe Rust calls; unwinding crosses no foreign
ABI and the staged output is dropped normally. Shape validation precedes
dispatch. Checkpoints do not reorder any output's fused accumulation chain.
The portable body and the dispatched body are compared directly in the local
product tests, in addition to comparison with the original kernels.

## Direct LDLT solve

The allocating direct solve uses the same block/panel order and fused reductions in both
arms. Pivot choices, permutations, diagonal 1x1/2x2 solves, and gather/scatter
operations are unchanged. The existing direct-solve residual fixtures compare
the dispatched result bitwise with a non-inlined portable caller, including
indefinite and delayed-pivot cases. A separate RHS-scale case covers signed
zeros and small/large finite values. This extends code generation only; it
does not assert a new cross-ISA qualification or make the direct solver the
instrument engine's default.

The in-place entry uses that same dispatch boundary. Its precomputed cycles
move permutation values without floating-point arithmetic, and direct global
updates follow the original block-then-panel fused reduction order. The
existing residual fixtures compare both solve entries bitwise, including
indefinite and delayed-pivot factors. The factor stays immutable; each caller
owns its output slice, and no mutable shared workspace is introduced.

## Positive-factor whitening

The two whitening actions reuse the same admitted feature check and the same
safe slice arithmetic in both dispatch arms. Scalar positive pivots are
admitted before construction; 2×2 blocks are explicitly refused. Permutation,
triangular-column order, inverse square roots, callback order and failure
semantics do not change at dispatch. No raw memory operation is introduced.
The portable and dispatched actions are compared bitwise for both transpose
directions and small, ordinary and large finite right-hand sides. Independent
checks verify Sᵀ A S = I and S Sᵀ b = A⁻¹ b on a nontrivial AMD permutation.
