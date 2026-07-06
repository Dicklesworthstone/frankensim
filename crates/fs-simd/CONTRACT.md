# CONTRACT: fs-simd

## Purpose and layer
SIMD tiers behind safe façades (plan §5.1, patch Rev Q): scalar reference,
NEON capsule (aarch64), AVX2/AVX-512 capsule (x86-64), one-shot dispatch.
Layer: L0.

## Public types and semantics
- `ops() -> &'static Ops` — function table resolved EXACTLY once from
  fs-substrate's tier; fields: axpy, scale, mul_elem, fma3 (all fused),
  dot, sum (fixed per-tier reduction shapes); `tier` for ledger keys.
- `scalar::*` — the semantic definition of every primitive (Tier 0).
- `neon::*` / `x86::*` — registered unsafe capsules (SAFETY.md beside each);
  all public capsule functions are SAFE (NEON is architecturally guaranteed;
  x86 façades re-verify CPU features and fall back to scalar).
- `is_cache_line_aligned`, `TernaryOp`.

## Invariants
- Elementwise ops match the scalar twin BITWISE on every tier (FMA policy:
  fused everywhere via mul_add — coordinated with fs-math's contraction
  policy).
- Reductions: fixed shape per tier (same tier + same input → same bits);
  cross-tier differences bounded by the documented envelope (machine
  identity, G5's domain — never run-to-run jitter).
- Tails handled by the scalar twin inside each function; no partial-lane
  intrinsic access exists.
- Length mismatches panic BEFORE any unsafe code (programmer-error contract).

## Error model
No fallible APIs; length mismatch = loud assert (documented programmer error).

## Determinism class
Deterministic per tier. Cross-tier: elementwise bitwise; reductions
envelope-bounded (feeds the G5 cross-ISA report).

## Cancellation behavior
Bounded allocation-free loops, no poll points; callers chunk work at tile
granularity (fs-exec discipline).

## Unsafe boundary
Two registered capsules: src/neon/mod.rs (THE exemplar — obligation from the
unsafe-safety-cases bead) and src/x86/mod.rs; both <300 lines with full
SAFETY.md files; enforced by `xtask check-unsafe`. Under Miri, dispatch
routes to scalar (intrinsics outside Miri's model; compensating equivalence
battery documented in the SAFETY files).

## Feature flags
None yet; `experimental-portable-simd` (Tier 2, nightly std::simd) arrives
when a consumer wants it — never load-bearing. `frontier-sme2` is the
separate fs-simd-sme2 bead.

## Conformance tests
tier_equivalence_battery (lens 0..67 × seeds, subnormal/NaN/±0/1e18 values,
bitwise + envelope), dispatch singleton + tier match, known answers
(bit-exact), alignment helper, loud length mismatch. VERIFIED EXECUTION:
aarch64-apple NEON (M4 Pro, local) and x86-64 AVX2 (Threadripper PRO 5995WX,
trj) — both green. Miri lane green (scalar dispatch).

## No-claim boundaries
- AVX-512 EXECUTION unverified (Zen 3 lacks it; compile-checked for both
  x86 targets; runs when a Zen 4/Sapphire-Rapids-class runner exists —
  ci-gauntlet-pipeline bead).
- x86 tier v1 covers axpy/dot/sum; scale/mul_elem/fma3 fall back to scalar
  there until fs-la's packing kernels demand them (<300-line capsule cap).
- No f32 variants yet (arrive with their consumers).
- No performance claims (roofline harness owns those).
