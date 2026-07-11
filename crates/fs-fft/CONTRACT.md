# CONTRACT: fs-fft

## Purpose and layer
Fast Fourier transforms for FrankenSim: 1D complex Stockham autosort FFT,
real-input transform (r2c), and DCT-II/III via FFT folding (the Chebyshev
transform path fs-cheb builds on). Layer: **L1 BEDROCK**. Depends only on
fs-math (strict-mode twiddles). Plan §6.3.

v1 was correctness-first radix-2; the 27d3 perf lane has since extended it
with the r2c **inverse** (c2r), **N-dimensional** (2D/3D) separable pencil
transforms, **mixed radix-8/4/2 Stockham** stages (three log₂ bits per
full-array pass; one radix-4-or-2 residue), and NEON/AVX2 q-run capsules for
the large-stride radix-4 stage. Remaining perf scope — cache-blocked pass
ordering, copy-back fusion, executor-tiled pencils, the still-unmet roofline
gate — see No-claim boundaries.

## Public types and semantics
- `C64 { re: f64, im: f64 }` — minimal complex scalar. `norm_sq` uses a fused
  multiply-add. Local to this crate until a shared complex home exists.
- `Fft::new(n)` — plan for power-of-two `n ≥ 1`; precomputes the half-length
  twiddle table `w[k] = exp(−2πik/n)` from `fs_math::det::{sin, cos}`.
  Plans are immutable and reusable across calls and threads.
- `Fft::forward(&mut data, &mut scratch)` — **unnormalized** DFT,
  `X[k] = Σ_j x[j]·exp(−2πijk/n)`. Both slices must have length `n`.
- `Fft::inverse(&mut data, &mut scratch)` — inverse DFT **scaled by 1/n**, so
  `inverse(forward(x)) = x` (round-trip tested).
- `RealFft::new(n)` / `RealFft::forward(&[f64]) -> Vec<C64>` — real input of
  power-of-two length `n ≥ 2`; returns the `n/2 + 1` non-redundant bins via
  half-size complex packing + untangling (half the complex work; oracle-tested
  against the embed-into-complex path).
- `RealFft::inverse(&[C64]) -> Vec<f64>` — c2r: reconstructs the `n` real
  samples from the `n/2 + 1` bins, the **exact algebraic inverse** of `forward`
  (solves the 2×2 untangle system, then the half-size complex inverse; Hermitian
  symmetry assumed per the standard c2r contract). Verified by r2c→c2r
  round-trip and against the full-size complex IFFT of the Hermitian completion.
- `FftNd::new(&[usize])` — plan an N-dimensional complex FFT over a **row-major**
  buffer (last axis contiguous); every axis a power of two ≥ 1. `shape()`,
  `total()` (product of dims = required buffer length), `forward(&mut [C64])`,
  `inverse(&mut [C64])` (`1/total`-normalized so `inverse(forward)=id`). The
  transform is **separable**: the planned 1D `Fft` applied along each axis in
  turn (row–column / pencil algorithm), deterministic by construction.
- `dct2(&[f64])` — unnormalized DCT-II, `X[k] = Σ_j x[j]·cos(πk(2j+1)/(2n))`,
  via even/odd folding + one complex FFT.
- `dct3(&[f64])` — DCT-III with the k=0 halving convention such that
  `dct3(dct2(x)) · 2/n = x` (round-trip tested).

## Invariants
1. Forward/inverse are exact round-trip inverses up to floating-point error
   (tested < 1e-12 relative at n ≤ 512).
2. Output agrees with the naive O(n²) DFT oracle to < 1e-12 relative error for
   all power-of-two n ≤ 512 (exhaustive over sizes; random inputs).
3. Parseval, linearity, impulse/constant, and circular-shift identities hold
   (tested).
4. Butterfly execution order is a pure function of (n, stage structure): no
   data-dependent branching, no threading in v1.
5. c2r inverse is an exact round-trip inverse of r2c forward (tested < 1e-11)
   and matches the full-size complex IFFT of the Hermitian-completed spectrum.
6. N-D transforms match a fully independent naive N-D DFT oracle to < 1e-12
   relative (2D and 3D), round-trip to identity, satisfy separability (row-then-
   column) and N-D Parseval, obey the 2D circular convolution theorem, and are
   bitwise deterministic across runs.

## Error model
Size violations (non-power-of-two, mismatched scratch length) panic with
structured messages naming the size and the remedy. These are programmer
errors, not runtime conditions — silently computing a wrong-size transform
would be worse. No other fallible paths.

## Determinism class
**Bit-deterministic cross-ISA by construction**: twiddles come from fs-math's
strict det functions; the transform body is fixed-order +, −, ×, mul_add.
Evidence: FNV-64 golden hash over 16 batches of n=128 forward outputs is
`0xbd55_68d2_33f4_b4bc`, recorded on aarch64-apple (M4 Pro) and required to
match on x86-64 (Threadripper) in the crate's own test suite. Golden-evidence
policy: the hash may only change with a semantically justified bump.

## Cancellation behavior
v1 transforms are single-tile, O(n log n), and short; no cancellation poll
points inside a single transform. Executor-tiled multi-dimensional transforms
(follow-up bead) will poll at pencil boundaries per Decalogue P7.

## Unsafe boundary
None. `unsafe_code` denied; no capsules.

## Feature flags
- `frontier-sixstep` [F] (default OFF, bead 27d3): the cache-blocked
  six-step path for n ≥ 2¹⁶ with even log₂ (square out-of-place tiled
  transposes + cache-resident √n row sweeps with fused twiddle, six
  full-array passes). CORRECT (cross-path agreement, transform laws, own
  golden 0x79aa_108f_a517_012f — new bit territory, the small-n stage
  golden is untouched by the pure-function-of-n dispatch, asserted both
  ways) but still MEASURED SLOWER than the radix-8/4/2 stage walk on M4.
  The fs-simd `trn1c64` NEON transpose capsule (one interleaved complex
  = one q-register; bounds checks hoisted) roughly HALVED the deficit
  (2026-07-11 interleaved relative lane at n = 2²⁰: forward ratio
  0.73 vs the pre-capsule ~0.44) without moving any bits (pure exact
  moves, twin-gated bitwise). Ships OFF until a quiet-host run shows a
  WIN (`tests/sixstep.rs` `--ignored` relative lane is the instrument);
  the remaining lever is fusing the ping-pong copy-back. Enabling
  changes large-n output bits by design.

## Conformance tests
In-crate suite (`cargo test -p fs-fft`): naive-DFT oracle sweep (n = 1..512),
impulse/constant/linearity, Parseval + shift theorem, r2c vs embedded-complex
oracle, c2r round-trip + full-IFFT oracle, DCT-II/III vs naive definitions +
round-trip, N-D (2D/3D) vs independent naive N-D DFT + round-trip + separability
+ N-D Parseval + 2D convolution theorem + determinism, determinism + golden
hash, structured rejection of bad sizes. Any reimplementation must pass this
suite bit-for-bit on the golden-hash case.

## No-claim boundaries
- **The ≥40% roofline gate is NOT met** (measured 2026-07-10, radix-8/4/2,
  corrected traffic model): aarch64 M4 gated sizes 0.225–0.248 attainment
  (raw throughput +11–18% over the radix-4/2 formulation — elems/s is the
  truth; the tighter model raised the roof), x86-64 ts2 0.165–0.184 (AVX2
  capsule only serves the residual radix-4 stage now). Remaining levers per
  the bead: make the gated six-step WIN (its transposes need the SIMD-tile
  capsule treatment — the scalar version measured slower and ships off,
  see Feature flags), fusing the ping-pong copy-back into the final stage
  pass, executor-tiled pencils with cancellation polls. Higher-radix golden bumps are pre-authorized with justification —
  bumped twice so far (radix-2→4/2→8/4/2), each recorded at the golden and
  the current value four-quadrant verified (M4 + ts2 × debug + release).
- The N-D transform is CORRECT and separable but not yet cache/execution
  optimized: it gathers each pencil into a temporary line (allocated per axis,
  reused across pencils) rather than blocking transposes or tiling on fs-exec.
- N-D real transforms (r2c/c2r with only the first axis half-length) are not
  shipped; N-D is complex-in/complex-out. Callers pack real fields as `C64`.
- No general-n (mixed-radix / Bluestein) support; power-of-two only.
- `C64` is not a public complex-arithmetic library; only what the FFT needs.
