# CONTRACT: fs-fft

## Purpose and layer
Fast Fourier transforms for FrankenSim: 1D complex Stockham autosort FFT,
real-input transform (r2c), and DCT-II/III via FFT folding (the Chebyshev
transform path fs-cheb builds on). Layer: **L1 BEDROCK**. Depends on fs-math
(strict-mode twiddles) and fs-simd (resolved exact-move/stage capsules).
Plan §6.3.

v1 was correctness-first radix-2; the 27d3 perf lane has since extended it
with the r2c **inverse** (c2r), **N-dimensional** (2D/3D) separable pencil
transforms, **mixed radix-8/4/2 Stockham** stages (three log₂ bits per
full-array pass; one radix-4-or-2 residue), NEON/AVX2 q-run capsules, and a
frontier six-step path whose logical transpose/copy structure is fused into
two full-array gather/scatter passes. Remaining perf scope — executor-tiled
pencils, the current x86 six-step verdict, and the still-unmet roofline gate —
is recorded under No-claim boundaries.

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
- `FftNd::{forward,inverse}_pooled(&mut [C64], &impl fs_exec::KernelRunner,
  &CancelGate)` — executor-tiled N-D passes (bead 27d3): per-axis pencil work
  tiled over outer blocks (or axis-0 column groups behind row mutexes), with
  per-pencil arithmetic and order EXACTLY the serial path's, so output is
  bitwise identical to `forward`/`inverse` at every worker count (the P2 law,
  gated by conformance). Tiling granularity (bead 3f6c, v2 kernels): both
  kernels enforce a ~4096-element per-tile work floor — outer blocks are
  grouped consecutively per tile (~8 tiles/worker target) and the column
  group cannot shrink below the floor at high worker counts. Grouping is
  timing-only; it never reorders pencils or changes bits.
- `FftNd::{forward,inverse}_pooled_observed(..., &mut dyn FnMut(NdPassReport))`
  (bead 3f6c) — the same transforms, additionally reporting each axis pass's
  geometry (`axis`, `kernel`, `n`, `stride`, `outer`, `tiles`, `completed`,
  `workers`) and `wall_ns` to the observer in execution order. `NdPassReport`
  is MEASUREMENT ONLY — envelope-class like `fs_exec::RunReport`; results
  never depend on it, and on cancellation the interrupted pass is still
  observed with `completed < tiles` (0 for a pre-requested gate).
- `dct2(&[f64])` — unnormalized DCT-II, `X[k] = Σ_j x[j]·cos(πk(2j+1)/(2n))`,
  via even/odd folding + one complex FFT.
- `dct3(&[f64])` — DCT-III with the k=0 halving convention such that
  `dct3(dct2(x)) · 2/n = x` (round-trip tested). DCT folding and
  post-rotation are part of `fs-fft:transform-bits=1`.

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
7. DCT-II/III match their direct definitions, round-trip, and replay bitwise
   within one build. There is no standalone cross-ISA DCT golden in fs-fft;
   registered fs-cheb and vessel goldens provide downstream change detection.

## Error model
Size violations (non-power-of-two, mismatched scratch length) panic with
structured messages naming the size and the remedy. These are programmer
errors, not runtime conditions — silently computing a wrong-size transform
would be worse. No other fallible paths.

## Determinism class
The complex stage path is **bit-deterministic cross-ISA**: twiddles come from
fs-math's strict functions and the transform body uses fixed-order `+`, `−`,
`×`, and `mul_add`. The current radix-8/4/2 FNV-64 golden over 16 batches of
n=128 outputs is `0x22dd_b617_266e_a792`, reproduced on aarch64 M4 Pro and
x86-64 ts2 in debug and release. History: radix-2
`0xbd55_68d2_33f4_b4bc`; radix-4/2 `0x0506_a4a0_955d_cf8e`.

The frontier six-step arithmetic has its own frozen hash
`0x79aa_108f_a517_012f`; two-pass fusion and vectorized exact moves preserved
it. DCT-II/III are fixed-order and same-build bit-replay tested, but fs-fft has
no standalone DCT bit golden. Current downstream sentinels are the fs-cheb
aggregate `0xeea0_4b0a_01de_46cd` and vessel smoke
`0x4541_d7f3_2926_1082`. The fs-cheb aggregate has current aarch64 and x86-64
debug/release reproduction; vessel smoke has aarch64 debug/release evidence
while its current x86-64 row remains pending. Golden hashes may change only
with a semantically justified bump.

## Cancellation behavior
v1 transforms are single-tile, O(n log n), and short; no cancellation poll
points inside a single transform. The executor-tiled N-D path polls at pencil
boundaries per Decalogue P7 and is a SCRATCH-TRANSFORM API: on `Err` the
buffer contents are unspecified (some axes/pencils may be transformed);
callers needing transactional output stage a copy first. The structured
`RunError` and the drained pool are the guarantees.

## Unsafe boundary
No local unsafe code. `unsafe_code` is denied in fs-fft. The fused six-step
uses fs-simd's safe dispatch facade; its AArch64 gather/scatter leaf is audited
in `crates/fs-simd/src/neon/transpose/SAFETY.md` and bitwise twin-gated.

## Feature flags
- `frontier-sixstep` [F] (default OFF, bead 27d3): the cache-blocked
  fused six-step path for power-of-two n ≥ 2¹⁶ with even log₂. Stage A
  gathers columns, runs cache-resident √n transforms plus fused twiddles,
  and scatters into scratch; stage B transforms scratch rows and scatters
  them into final output columns. This is **two** full-array passes, with no
  materialized transpose or copy-back. CORRECT (cross-path agreement,
  transform laws, exact-move twin battery, own golden
  `0x79aa_108f_a517_012f`) but still MEASURED SLOWER than the stage walk on
  M4 after vectorized gather/scatter. At n = 2²² on 2026-07-11, six-step was
  0.0822–0.0852 s versus 0.053–0.055 s for stages (ratio 0.64–0.67).
  x86 verdict (2026-07-11, ts2 5975WX, n = 2²², 3 runs): ratio 0.80–0.87 —
  weaker prefetch narrows the gap but does not invert it; golden verified on
  x86 debug + release. PERMANENTLY FRONTIER on both reference ISAs — every
  recorded lever (radix-8, transpose capsule, 6→2 fusion, vectorized strips)
  is exhausted; the relative lane stays armed as the instrument of record.
  Enabling it changes large-n output bits by design. Its roofline
  evidence version is `27d3-6s-fused2`, separate from the transform-bit
  version because exact-move optimization changed traffic without moving bits.

## Conformance tests
In-crate suite (`cargo test -p fs-fft`): naive-DFT oracle sweep (n = 1..512),
impulse/constant/linearity, Parseval + shift theorem, r2c vs embedded-complex
oracle, c2r round-trip + full-IFFT oracle, DCT-II/III vs naive definitions +
round-trip + same-build bit replay, N-D (2D/3D) vs independent naive N-D DFT +
round-trip + separability + N-D Parseval + 2D convolution theorem +
determinism, determinism + golden hash, and structured rejection of bad sizes.
`tests/conformance.rs` additionally registers three fixed-order exact Casebook
records: the n=2 complex KAT and exact round-trip; the current 16-batch n=128
stage fingerprint `0x22ddb617266ea792`; and forward/inverse serial-vs-pooled
bit equality over five shapes and worker counts `1, 2, 3, 7`, using an explicit
fixed two-CCD/eight-core topology rather than a host probe. Every failure
retains the complete canonical input frame, including that topology, plus
exact mismatch bits or the structured run error. A disclosed `0xF5FF0001`
seeded corruption flips bit 0 of the stage oracle and proves both the typed red
report and `assert_green` refusal paths. These records make central failures
replay-complete; one local record is not, by itself, a fresh dual-ISA G5 proof.
The declared G3 adopter `forward-signal-scale-equivariance`
(`tests/metamorphic.rs`, seed `0x2ACE_0002`, 384 cases) applies non-identity
power-of-two rescalings to exactly eight generated complex samples and checks
every output component from `Fft::forward` at `1e-13` absolute-relative
tolerance. The fixed oracle, theorem, round-trip, and golden pins above remain
authoritative and independent.
The frontier dispatch battery also rejects non-power-of-two lookalikes. The
performance lane binds its two-pass traffic count to evidence version
`27d3-6s-fused2`. Its historical-axis input is report-only unless the baseline
is an attested envelope accepted under the configured promotion-authority
policy from `FRANKENSIM_PROMOTION_AUTHORITY_POLICY`, with all named source
receipts supplied in `FRANKENSIM_RETAINED_SOURCE_RECEIPTS`. The lane captures
one atomic authority decision and embeds the full frozen pre/post snapshot in
the final gate JSON so the measured claim cannot be detached from its decision.
A positive gate additionally requires `FRANKENSIM_ROOFLINE_LEDGER`: before
emitting `citation_eligible:true`, the lane atomically records the exact
admission receipt and exact final-gate JSON through the shared `fs-roofline`
external-gate protocol and verifies the stored bytes. A missing or empty ledger
path leaves the completed measurement report-only; a write or re-read failure
fails closed and cannot emit a positive gate.
The citable policy owns its clock; clock failure or an epoch-day rollover
between mint and post-probe invalidates attested evidence and ends the lane as
`environment_invalid`. Configuration refusals alone remain measured
report-only observations.
Missing, partial, malformed, denied, revoked,
tampered, or cross-machine authority inputs never produce a positive gate. Any
reimplementation must pass this suite bit-for-bit on the golden-hash cases.
The retained-source file is a protected hash-inventory declaration; the lane
does not fetch or independently prove availability of the named source bytes.
Its conformance matrix removes one named receipt to prove that missing source
evidence stays report-only, then re-endorses the identical baseline under a
rotated key and proves that the baseline hash stays fixed while the key and
authority-policy receipt move.

## No-claim boundaries
- **The ≥40% roofline gate is NOT met** (measured 2026-07-10, radix-8/4/2,
  corrected traffic model): aarch64 M4 gated sizes 0.225–0.248 attainment
  (raw throughput +11–18% over the radix-4/2 formulation — elems/s is the
  truth; the tighter model raised the roof), x86-64 ts2 0.165–0.184 (AVX2
  capsule only serves the residual radix-4 stage now). Remaining levers per
  the bead: executor-tiled pencils with cancellation polls and further
  stage-path locality work. The alternative six-step path has already landed
  two-pass fusion and vectorized strip moves; it remains a measured negative
  on M4 rather than an unimplemented lever, and its current x86 verdict is
  not claimed. Higher-radix golden bumps are pre-authorized with justification
  — bumped twice so far (radix-2→4/2→8/4/2), each recorded at the golden and
  the current value four-quadrant verified (M4 + ts2 × debug + release).
  The explicit perf lane reports that target and, on an authority-admitted
  snapshot, enforces only a 15% anti-collapse floor. The committed plain
  baseline files make historical rows candidate observations, not
  authority-admitted citations; an
  environment-invalid measurement fails closed rather than returning a green
  test with no admissible evidence.
- The N-D transform is CORRECT and separable but not yet cache/execution
  optimized: it gathers each pencil into a temporary line (allocated per axis,
  reused across pencils) rather than blocking transposes or tiling on fs-exec.
- N-D real transforms (r2c/c2r with only the first axis half-length) are not
  shipped; N-D is complex-in/complex-out. Callers pack real fields as `C64`.
- No general-n (mixed-radix / Bluestein) support; power-of-two only.
- `C64` is not a public complex-arithmetic library; only what the FFT needs.
