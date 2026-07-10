# SAFETY: fs-simd/sme2 (bead wf9.3, feature `frontier-sme2`)

## What the unsafe does
Two `asm!` blocks:
1. `rdsvl` — reads the streaming vector length. Architecturally defined
   whenever FEAT_SME is present; the OS probe (sysctl subprocess /
   /proc/cpuinfo — no FFI) runs FIRST and gates the call. Read-only,
   `nomem/nostack/pure`.
2. The 16×16 streaming-mode GEMM microkernel: `smstart` … `fmopa` on
   za0.s … horizontal-slice stores … `smstop`, all inside ONE block.

## Why it is sound
- STREAMING-MODE CONTAINMENT: `smstart`/`smstop` both live inside the
  same asm block, so the compiler never schedules Rust FP/SIMD code
  while streaming mode is active, and ZA state never escapes.
- REGISTER DISCIPLINE: all v0–v31 are declared clobbered (z registers
  alias them). rustc allocates neither SVE predicate nor ZA registers
  under any enabled codegen feature, and the block touches only p0/p1
  and za0.s beyond declared operands. w12 is bound because ZA slice
  indices architecturally require w12–w15.
- BOUNDS: the safe facade computes `k * TILE` with checked arithmetic
  and asserts panel/tile lengths before capability dispatch or the block;
  the loop performs exactly k 64-byte loads per panel and 16 64-byte
  stores to c — all within the asserted allocations. The kernel
  overwrites the output tile; it does not read or accumulate the prior
  contents of `c`.
- CAPABILITY GATING: `sme2_available()` (OS probe AND SVL == 64 bytes)
  guards every kernel entry with a loud assert; on non-SME hardware the
  flag is inert and nothing here executes.

## Compensating checks
- G0 equivalence battery vs the scalar mul_add twin (identical k-order)
  on the same panels; capability-absent fallback test runs everywhere.
- Miri: the module is compiled out under Miri (asm is outside its
  model), matching the neon/x86 capsule convention.

## Residual risk
Assembler/OS support for SME2 is young (the bead's named risk). The
mitigation IS the tier: exploratory flag, runtime gate, NEON committed.
