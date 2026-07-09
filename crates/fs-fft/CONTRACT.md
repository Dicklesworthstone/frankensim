# CONTRACT: fs-fft

## Purpose and layer
Fast Fourier transforms for FrankenSim: 1D complex Stockham autosort FFT,
real-input transform (r2c), and DCT-II/III via FFT folding (the Chebyshev
transform path fs-cheb builds on). Layer: **L1 BEDROCK**. Depends only on
fs-math (strict-mode twiddles). Plan §6.3.

v1 is correctness-first radix-2, now extended (bead fs-fft-perf-multidim) with
the r2c **inverse** (c2r) and **N-dimensional** (2D/3D) transforms via separable
pencil decomposition. The remaining perf scope — radix-4/8 kernels, SIMD lanes,
cache-blocked transposes, executor-tiled pencils, and the roofline gate — stays
recorded follow-up; see No-claim boundaries.

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
None.

## Conformance tests
In-crate suite (`cargo test -p fs-fft`): naive-DFT oracle sweep (n = 1..512),
impulse/constant/linearity, Parseval + shift theorem, r2c vs embedded-complex
oracle, c2r round-trip + full-IFFT oracle, DCT-II/III vs naive definitions +
round-trip, N-D (2D/3D) vs independent naive N-D DFT + round-trip + separability
+ N-D Parseval + 2D convolution theorem + determinism, determinism + golden
hash, structured rejection of bad sizes. Any reimplementation must pass this
suite bit-for-bit on the golden-hash case.

## No-claim boundaries
- **No performance claims yet**: the kernels are scalar radix-2 (1D and the N-D
  pencils). Radix-4/8 higher-radix stages, SIMD lanes (fs-simd Ops), cache-
  blocked transposes, executor-tiled pencils with cancellation polls at pencil
  boundaries, and the roofline gate (≥40% of memory-bound peak on both reference
  ISAs, denominator = fs-substrate STREAM triad) all remain the perf follow-up.
  A higher-radix kernel that changes twiddle-application order may legitimately
  bump the golden hash — bump once, justify here, re-verify cross-ISA on trj.
- The N-D transform is CORRECT and separable but not yet cache/execution
  optimized: it gathers each pencil into a temporary line (allocated per axis,
  reused across pencils) rather than blocking transposes or tiling on fs-exec.
- N-D real transforms (r2c/c2r with only the first axis half-length) are not
  shipped; N-D is complex-in/complex-out. Callers pack real fields as `C64`.
- No general-n (mixed-radix / Bluestein) support; power-of-two only.
- `C64` is not a public complex-arithmetic library; only what the FFT needs.
