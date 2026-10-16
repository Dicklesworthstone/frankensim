# CONTRACT: fs-simd

## Purpose and layer
SIMD tiers behind safe façades (plan §5.1, patch Rev Q): scalar reference,
NEON capsule (aarch64), AVX2/AVX-512 capsule (x86-64), one-shot dispatch.
The table also carries operation-specific tier receipts: f64 GEMM reports
AVX2 on an AVX-512-capable host because its current audited microkernel is
AVX2/FMA, and reports scalar when either AVX2 or FMA is unavailable, preventing
tune evidence from claiming a tier that did not execute. This exact predicate
is resolved once while constructing the process-wide table; hot GEMM tiles call
the selected microkernel pointer without repeated feature detection. The x86
GEMM submodule owns selection and keeps its unchecked target-feature thunk
private, so safe sibling modules cannot bypass that one-shot proof.
Layer: L0.

## Public types and semantics
- `ops() -> &'static Ops` — function table resolved EXACTLY once from
  fs-substrate's tier; fields: axpy, scale, mul_elem, fma3 (all fused),
  dot, sum (fixed per-tier reduction shapes), mk8x4_f64 (the 8×4 f64
  GEMM register microkernel over packed k-fastest panels — BITWISE
  across tiers by the k-ascending fused-order contract; NEON capsule
  on aarch64, AVX2+FMA capsule on x86 with an exact-feature scalar fallback —
  bead xdgf), btile4x4_f64 (the batched-GEMM 4×4 entry-tile
  microkernel over plane-SoA batches, bead 9ekv — same BITWISE
  contract, lanes are independent matrices), trn1c64 (8×8-tiled
  out-of-place transpose of an n₁×n₁ interleaved-complex matrix —
  fs-fft's six-step tile pass, bead 27d3: pure exact 16-byte moves,
  one complex per q-register on NEON, BITWISE trivially; scalar on
  x86 until a measured need); `tier` for ledger keys.
- `scalar::*` — the semantic definition of every primitive (Tier 0).
- `neon::*` / `x86::*` — registered unsafe capsules (SAFETY.md beside each);
  all public capsule functions are SAFE (NEON is architecturally guaranteed;
  x86 façades re-verify CPU features and fall back to scalar).
- `sme2::*` (bead wf9.3, feature `frontier-sme2`, [F]) — the EXPLORATORY
  streaming-mode GEMM prototype: `sme2_available()` (runtime probe —
  sysctl subprocess on macOS / cpuinfo on Linux, NEVER compile-time
  assumed; also requires SVL = 512 bits, the fixed prototype shape),
  `streaming_vl_bytes()`, `gemm_tile_f32` (16×16 fmopa outer-product
  microkernel on za0.s, one self-contained smstart…smstop asm region),
  and the scalar mul_add twin. NEVER in the `ops()` table — NEON stays
  the committed path; promotion requires beating NEON across the
  autotuner shape sweep (the xdgf perf lanes' call, unclaimed here).
  MEASURED on Apple M4 Pro (ledgered): G0 equivalence vs the scalar
  twin is BITWISE (worst 0 ULP over k ∈ {1,3,17,64,257}); 263 GFLOP/s
  vs the release-build scalar twin's 8.9 (29.7×) at 16×16×1024 —
  evidence, not a perf gate; the honest NEON comparison belongs to the
  autotuner sweep. NO cross-ISA determinism-mode claim (the bead's
  explicit non-goal) until the G5 report characterizes streaming-mode
  accumulation across SVL classes.
- `is_cache_line_aligned`, `TernaryOp`, `Mk8x4`, `Btile4x4`.

## Invariants
- Elementwise ops match the scalar twin BITWISE on every tier (FMA policy:
  fused everywhere via mul_add — coordinated with fs-math's contraction
  policy).
- Reductions: fixed shape per tier (same tier + same input → same bits);
  cross-tier differences bounded by the documented envelope (machine
  identity, G5's domain — never run-to-run jitter).
- Tails handled by the scalar twin inside each function; no partial-lane
  intrinsic access exists.
- Length and checked-extent mismatches panic BEFORE any unsafe code
  (programmer-error contract), including `4·run_len` for radix-4 FFT output.

## Error model
No fallible APIs; length mismatch = loud assert (documented programmer error).

## Determinism class
Deterministic per tier. Cross-tier: elementwise bitwise; reductions
envelope-bounded (feeds the G5 cross-ISA report).

## Cancellation behavior
Bounded allocation-free loops, no poll points; callers chunk work at tile
granularity (fs-exec discipline).

## Unsafe boundary
Three registered capsules: src/neon/mod.rs (THE exemplar — obligation from
the unsafe-safety-cases bead), src/x86/mod.rs, and src/sme2/mod.rs
(feature-gated [F]; streaming-mode containment + full register discipline
documented in its SAFETY.md); all <300 lines with full SAFETY.md files;
enforced by `xtask check-unsafe`. Under Miri, dispatch
routes to scalar (intrinsics outside Miri's model; compensating equivalence
battery documented in the SAFETY files).

## Feature flags
None yet; `experimental-portable-simd` (Tier 2, nightly std::simd) arrives
when a consumer wants it — never load-bearing. `frontier-sme2` is the
separate fs-simd-sme2 bead.

## Conformance tests
The shared `fs-casebook` runner records three bounded, replay-complete cases:
exact primitive known answers, the exhaustive 16-row pure
`mk8x4_f64_tier_for` admission truth table, and resolved dispatch-table versus
scalar-twin bit equivalence over 30 vector fixtures plus eight packed 8×4 tile
fixtures. Records name the global dispatch tier and the effective 8×4 tier
separately, preserving operation-specific fallback truth. A disclosed
`0xF51D_0001` seeded reference corruption proves both the structured red record
and `assert_green` refusal path. These records are portable G0 evidence suitable
for separately retained cross-host execution; they do not replace the broader
source-level battery or themselves constitute a G5 cross-ISA comparison.

That battery covers lengths 0..67 × finite seeded values for reduction
envelopes, explicit subnormal/qNaN/sNaN/infinity/invalid-product/±0
microkernel cases, dispatch singleton + tier match, alignment, and loud length
mismatch refusal. VERIFIED EXECUTION: aarch64-apple NEON (M4 Pro, local) and
x86-64 AVX2 (Threadripper PRO 5995WX, trj) — both green. Miri lane green
(scalar dispatch). The structured Casebook tranche remains proof-pending until
the central batch records its requested package run.

## No-claim boundaries
- AVX-512 EXECUTION unverified (Zen 3 lacks it; compile-checked for both
  x86 targets; runs when a Zen 4/Sapphire-Rapids-class runner exists —
  ci-gauntlet-pipeline bead).
- x86 tier v1 covers axpy/fma3/fmacc/dot/sum; scale and mul_elem fall back to
  scalar there (<300-line capsule cap).
- No f32 variants yet (arrive with their consumers).
- No performance claims (roofline harness owns those).
