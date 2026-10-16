# SAFETY — fs-simd NEON transpose capsule

Scope: `crates/fs-simd/src/neon/transpose/mod.rs` (`trn1c64`, `gath8c64`,
and `scat8c64`), compiled only on `target_arch = "aarch64"` where NEON is
baseline. No runtime feature gate is needed for `vld1q_f64`/`vst1q_f64`.

## Invariants the `unsafe` block relies on

1. **Matrix extent.** `checked_trn1c64_len` constructs `2·n1²` without
   wrapping, and matrix slice lengths are asserted equal to it before any
   pointer math. In `trn1c64`, loop structure guarantees `i < n1` and
   `j < n1`, so every element index `j·n1 + i` and `i·n1 + j` is `< n1²`.
   The dereferenced f64 offsets `2·idx` and `2·idx + 1` (one q-register =
   two f64) are therefore strictly inside both slices.
2. **Strip extent.** `gath8c64` and `scat8c64` validate
   `g.checked_add(8)` and require the resulting exclusive end to be at most
   `n1` before creating a pointer. Thus `g + c < n1` for every `c in 0..8`
   without relying on wrapping arithmetic. Their dense buffer length is
   exactly `2·8·n1` f64 values. For each matrix row `i`, `i·n1 + g + c` is
   `< n1²`; for each buffer column `c`, `c·n1 + i` is `< 8·n1`. Doubling
   either complex index and loading or storing two f64 values stays within
   the asserted slice. The near-`usize::MAX` regression test verifies
   refusal happens at this checked boundary before the NEON pointer path.
3. **Aliasing.** Every source `&[f64]` and destination `&mut [f64]` pair is
   held as simultaneous shared and exclusive borrows, so the slices cannot
   overlap; loads and stores never alias.
4. **Alignment.** `vld1q_f64`/`vst1q_f64` are unaligned-tolerant on
   AArch64; no alignment precondition is assumed.
5. **No floating-point arithmetic.** The capsule performs exact 16-byte
   moves only; bitwise equality with the scalar twin is move equality, gated
   in `tier_equivalence_battery`. Transpose coverage includes square,
   non-multiple-of-8, and `n1 = 1` edges. Strip coverage includes first,
   middle, and final eight-column groups, plus scatter-of-gather identity.

## Twin

`scalar::{trn1c64,gath8c64,scat8c64}` implement the identical tiling or
strip iteration order with bounds-checked indexing. The battery asserts
bitwise equality of every full destination for each tested shape. Both
tiers also share the same checked near-`usize::MAX` strip refusal contract.
