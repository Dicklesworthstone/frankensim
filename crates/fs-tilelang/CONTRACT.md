# fs-tilelang CONTRACT

## Purpose and layer

Layer: **L0 SUBSTRATE** (deps: fs-substrate, fs-tilelang-macros). The
safe tile-kernel DSL runtime (plan patch Rev C): one restricted
`kernel!` body lowers to a scalar reference variant, a lane-shaped
variant, per-kernel arithmetic-intensity metadata (P6: every kernel
ships its intensity analysis), and auto-generated G0/G5 twin tests —
so hot loops stop being hand-maintained N times and "the primal
changed but the SIMD variant didn't" rot becomes structurally
impossible. fs-opdsl (tfz.4) generates INTO this layer.

## Public types and semantics

- `kernel!` (re-exported macro — see fs-tilelang-macros for grammar):
  generates a module with `META: KernelMeta`, `run_scalar(…)`,
  `run(…)` (dispatches ONCE on the resolved tier), and
  `run_lanes::<LANES>(…)` (monomorphized lane-grouped loops for
  LANES ∈ {2, 4, 8}).
- `KernelMeta { name, flops_per_elem, bytes_per_elem, halo,
  reduction, determinism }` with `intensity()` (FLOP/byte, the
  roofline x-axis) and address-free JSON `descr()`. Dynamic (uparam)
  halos record 0 in META — the metadata is static by design.
- `DeterminismClass::{BitwiseAllTiers, PerTier}`;
  `ReductionKind::{None, DeterministicSum, FastSum}`.
- `lane_width()` — 1/2/4/8 from fs-substrate's cached tier dispatch
  (Scalar/NEON/AVX2/AVX-512); resolved once, never in hot loops.
- `deterministic_sum` / `fast_sum` — the reference reduction
  combiners; `REDUCTION_CHUNK = 64` fixes the deterministic tree
  shape as a function of LENGTH ONLY (never tier, never thread).

## Invariants

- Scalar and lane variants are BITWISE-EQUAL by construction: the
  lane variant runs token-identical per-element arithmetic in
  ascending index order (lane grouping shapes the loop for the
  autovectorizer without touching arithmetic), and reductions fold
  identical fixed chunks in identical order.
- Generated kernels write only declared buffers over the halo-clipped
  range; reads/writes never alias (macro-enforced).
- `DeterministicSum` results are a pure function of the input values
  and length. `FastSum` currently coincides but contractually MAY
  reassociate — its bit pattern must never enter a golden.

## Error model

Buffer-length mismatches and out-of-range shift/gather indices panic
(safe slice indexing — loud, never UB). All shape/grammar errors are
compile-time (see fs-tilelang-macros).

## Determinism class

Bit-deterministic across tiers, runs, and ISAs for everything except
`FastSum` (envelope-bounded vs the deterministic variant). Generated
twin tests enforce G0 tier equivalence at ALL lane widths (not just
the resolved one) and G5 repeat determinism.

## Cancellation behavior

Generated kernels are bounded synchronous loops over one tile's
worth of data; drivers own chunking to tile quanta and Cx poll
points between tiles (the fs-exec discipline; same policy as
fs-la/fs-simd).

## Unsafe boundary

None. `unsafe_code = "deny"`; the macro REJECTS `unsafe` in kernel
bodies.

## Feature flags

None.

## Conformance tests

`tests/tilelang_battery.rs`: the three acceptance-criteria reference
kernels written once — batched axpy (map), 3-point and 3D 7-point
stencils (literal and stride-uparam shifts, halo discipline: NaN
canaries untouched), SDF-style trilinear grid sample (gather form) —
each checked against hand-written oracles bitwise; a deterministic
dot-product reduction equal at every lane width and equal to the
fixed-shape reference combiner; META flop/byte counts asserted and
the per-kernel intensity table logged (roofline food). AUTO-GENERATED
twin tests (`__twin_tests` modules) run alongside for every kernel
the macro can drive (no gather, no uparams — those drive their own
twins in the battery, the adoption-policy boundary).
`tests/compile_fail.rs`: 10 rejection fixtures through the in-house
offline harness (alias, unsafe, loops, allocation, unknown reduction,
acc mismatches, unassigned writes, undeclared gather targets,
reserved names).

## No-claim boundaries

- v1 lowering is shaped-safe-Rust + LLVM autovectorization on
  monomorphized lane loops — NO explicit NEON/AVX-512 intrinsics and
  NO performance claims. If measured parity vs hand-written kernels
  is insufficient, intrinsic microkernels join the fs-simd capsule
  lane (xdgf/9ekv precedent). Performance parity is DOCUMENTED when
  measured, never assumed.
- No accelerator targets (capability-gated future scope per the
  bead), no multi-buffer typed layouts beyond f64/u32, no Qty flow
  through kernel bodies yet (unit checking joins fs-opdsl's typed
  IR), no Morton-tile geometry binding (kernels see flat index
  spaces; tile identity is the driver's).
- `FastSum` has no reassociated implementation yet — it is the
  deterministic fold plus a contractual right to change.
- Hand-written kernels remain allowed; they must hand-supply META
  and twin tests (the adoption-policy lint is future scope).
