# CONTRACT: fs-fft

## Purpose and layer
Fast Fourier transforms for FrankenSim: 1D complex Stockham autosort FFT,
real-input transform (r2c), and DCT-II/III via FFT folding (the Chebyshev
transform path fs-cheb builds on). Layer: **L1 BEDROCK**. Depends on fs-math
(strict-mode twiddles) and fs-simd (resolved exact-move/stage capsules).
Plan §6.3.

`SmallFft` is the bounded, cold-planned radix-2 owner for real-time receiver
partitions (complex lengths 16 through 256). It allocates only at preparation,
uses deterministic twiddles and ordinary multiply/add arithmetic, and scales
the inverse by 1/N. Its rounding is separate from `Fft`'s transform-bit contract;
adding it does not alter that contract or regenerate existing goldens. Its
focused check compares every supported length to `Fft` and checks round trips;
receiver consumers additionally compare the resulting convolution and timing
against their direct FIR. This is not a general large-transform replacement.

## Exact-length physical-record transforms

`GeneralFft` admits positive arbitrary 1D complex record lengths, including
odd and prime N. Power-of-two records call the original `Fft` directly and
preserve its transform bits. Other sizes use Bluestein's chirp convolution
with M = next_power_of_two(2*N-1), retaining the original N output bins and
their spacing 1/(N*dt). Internal convolution padding does not extend the
physical observation window or resample its clock. Forward transforms remain
unnormalized and inverse transforms include 1/N.

The original Stockham numerical owner supplies both convolution transforms.
The chirp and transformed convolution kernel are cached at preparation;
subsequent calls reuse caller scratch. Integer j*j modulo 2*N is formed in
u128 before conversion to floating point, then reduced to an exact quadrant
and small residual angle before strict trigonometry. This preserves weak
quadrature components near quarter turns at large indices. Bluestein inputs
are normalized by their largest component before convolution and restored to
their physical scale afterward, avoiding unnecessary intermediate range loss.
This does not extend f64's exponent range or certify arbitrary cancellations.

`GeneralFftBudget` bounds N and the combined numeric plan/minimum-scratch
payload before allocation. `required_workspace_bytes` performs allocation-free
admission planning; `scratch_len` states the reusable complex buffer length.
The count includes preparation temporaries and any recursive six-step row
buffers, but excludes transformed data, extra caller scratch, stack and
allocator metadata. Default evaluation allocates nothing; enabling the existing
`frontier-sixstep` route retains its explicitly counted transient allocations.
Work and storage are O(M log M) and O(M). This synchronous owner inherits the
existing transform's lack of mid-kernel cancellation, not a deadline guarantee.

Zero, inexact-integer or unaddressable sizes, sample/storage limits, malformed
buffers and nonfinite inputs refuse. Buffer/input refusals preserve both
buffers. Numerical overflow is reported; the direct in-place route may have
modified data, while the Bluestein route publishes only complete finite output.
The new G0/G1/G3 tests compare direct DFT sums, odd-grid highest positive tones,
prime-length Parseval/round trips, physical amplitude scales, deterministic
replay, original power-of-two bits, modular chirps and admission boundaries.
These are finite-record DFT checks; no leakage removal, alias prevention,
automatic tone discovery, PSD calibration, continuum-frequency error bound,
new cross-ISA golden or throughput/roofline claim is added.

`harmonic_window` supplies analytic normalized-window means and real-product
moments for peak phasors `Re(a exp(-i omega t))`. Positive-time consumers
conjugate their amplitudes. `HarmonicWindow` retains start and duration
separately; `HarmonicProduct` caches each frequency pair's midpoint rotations
and centered cosine/sine Gram weights for allocation-free reuse at every
spatial point and component. Both beat and sum-frequency terms are retained,
including DC, signed frequencies and subperiod exposures. Short-window sine
energy is evaluated with a `sinc-1` series rather than subtracting rounded
numbers near one. Strict trigonometry and double-double phase/products retain
small physical terms without changing any FFT transform bits.

`product_kernel_dd` and `proper_product_kernel_dd` additionally accept retained
double-double angular frequencies. They normalize finite high/low components
before the original Gram/rotation implementation, so a resolved offset smaller
than a binary64 step of its common carrier is not erased at the API boundary.
The proper kernel subtracts the retained frequencies before applying the clock.
The real DD kernel likewise forms its beat before interval multiplication and
rotates the high and low frequency components separately; scaling a large
carrier first can otherwise erase a representable slow beat.
Spatial consumers can use the same algebra for wavenumbers and lengths. The
original scalar APIs keep their arithmetic and normalization; focused controls
check sub-ulp real beats, exact proper common-carrier cancellation, invalid
components and bitwise compatibility when the low input components are zero.

Carrier rotations and the sine/cosine inputs of centered Gram weights use
`fs_math::det::sin_cos_dd` without an intermediate f64 rounding. Independently
rounded sinc values can otherwise make the mean-square pressure of two nearly
equal opposite-phase tones negative, even when the final sum uses DD. The
focused beat regression integrates the real difference waveform through its
stable sine-product identity, with uniform and triangular shutters. No negative
result is clipped and no relative guarantee at arbitrary nulls is added.

G0/G1/G3 in-module tests compare real-signal quadrature, complete-cycle peak
normalization, subperiod energy/flux, near-DC products, clock offsets and typed
refusals. Input clocks, amplitude units, geometry, and source coherence remain
the domain caller's responsibility. Invalid/nonfinite and unrepresentable
inputs return `WindowError`; no sampling-rate assumption, inferred coherence,
interval certificate, or relative-accuracy guarantee near arbitrary nulls is
claimed. Each cached scalar contraction is bounded constant work; callers own
the line-count/pair-work budgets and cancellation of larger field batches.

`HarmonicWindow::proper_product_kernel` additionally integrates explicitly
authored jointly proper complex covariance factors: `E[z_r]=0`,
`E[z_r conj(z_s)]=delta_rs`, and `E[z_r z_s]=0`. Its cached
`ProperHarmonicProduct` returns `1/2 Re[a conj(b) K(omega_a-omega_b)]`, where
`K` is the selected profile's normalized average of the negative-time complex phase.
The factor 1/2 is included. Properness eliminates sum-frequency terms even
for a subperiod shutter; cross-frequency correlations can remain nonstationary.
Deterministic means still use `HarmonicProduct`, with their sum terms intact.
The input factorization supplies square-root statistical weights and positive
semidefinite covariance; this numerical owner does not estimate either.

The compensated frequency difference is formed before clock multiplication.
Both cached kernel components and the product accumulation remain double-double;
half-angle `cos-1` and the `sinc-1` series preserve weak destructive beats that
round to zero if midpoint cosines are first rounded to f64. Equal-frequency
normalization is exact even when individual carrier phases are unrepresentable.
Focused tests use an independent four-phase real-field ensemble and time
quadrature, a dyadic tiny-beat series, common-carrier/clock invariance and typed
refusals. These are exposure-averaged instantaneous second moments, not the
variance of the window-averaged field, random realizations, detector uncertainty,
Gaussianity, source qualification, or arbitrary-null relative certificates.

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
- `FftNd::{forward,inverse}_pooled_scoped(&mut [C64], &TilePool,
  &asupersync::Cx, &CancelGate)` — the same N-D pencil transforms under the
  caller's live asupersync task scope (bead lx0e). Each parallel pass uses
  `TilePool::run_scoped`, so task cancellation or budget exhaustion reaches
  the pool's ordinary drain boundary; serial-pencil boundaries checkpoint the
  same task before transforming. Arithmetic and output bits remain those of
  the local pooled API. This is the task-owned path. The generic
  `KernelRunner` entry points remain explicit local-worker APIs and do not
  imply task ancestry.
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
would be worse. Separate harmonic-analysis APIs have the typed refusals
described in their respective sections.

### Known-frequency harmonic estimation

`harmonic::HarmonicFitPlan` jointly fits a constant and the caller's declared
positive harmonic orders to a uniformly sampled record. Frequencies are not
rounded to FFT bins; partial-period and non-power-of-two records are admitted.
The shared `fs-la::factor::qr` Householder owner solves the complete design,
without normal equations or independent nonorthogonal tone projections.
`HarmonicFit::phasors` are peak amplitudes under `exp(-i omega t)` with time
zero at the first sample. The jointly fitted DC and unexplained sampled RMS
retain the input's units. Reusing `fit_into` needs only caller-provided sample
scratch and phasor output; it performs no allocation.

The sampled basis and frequency/pole derivatives retain `frequency_hz * dt * j`
in DD until whole cycles are removed, then evaluate the reduced phase with the
shared deterministic scalar trigonometry. The low radian component is applied
before rounding the design columns. Harmonic frequencies are authored once as
`fundamental_hz * order`, as in the physical propagation consumers: a decimal
clock does not distinguish that front from identical explicit tones or
zero-decay poles. This prevents record-length-dependent phase drift from
polluting weak partials and the residual; it adds no storage or fit-time work.
Focused G1/G3 regressions compare long records with the independent DD
observation clock and check exact fit parity between all three fronts.

Admission checks the explicit storage budget, sample count, harmonic ordering,
strict Nyquist limit and a caller-set reciprocal-condition floor. Conditioning
is measured using the complete triangular `R` inverse's 1-norm, not a diagonal
ratio. Input centering and scaling preserve small changes on a large preload
and avoid norm overflow across physical unit scales. Invalid data, insufficient
storage and unresolved columns return `HarmonicFitError`; failed fits do not
partially overwrite caller phasors.

Preparation is O(N K^2 + K^3) and O(N K) storage for N samples and K tones;
reuse is O(N K + K^2). This is a known-frequency, finite-dimensional estimate,
not frequency identification, an anti-alias filter, a statistical uncertainty
bound or evidence of physical model validity. Omitted tones, noise, transients
and frequency error can bias fitted components; the sampled residual discloses
unexplained content without certifying its cause. Existing FFT bit semantics
are unchanged. This synchronous path has no mid-factorization cancellation;
same-build replay is tested, with no new cross-ISA certification claim.

G1/G3 tests in `tests/harmonic_fit.rs` recover weak even harmonics on fractional
period records, verify phase-origin shifts and unit scaling, preserve small
harmonics on a large preload, measure an analytically known omitted harmonic's
RMS, and refuse underresolved/invalid fits. The sinusoidal model and the reason
to fit its constant jointly follow the [NIST sinusoidal-model discussion](https://www.itl.nist.gov/div898/handbook/eda/section3/eda33r3.htm);
the record-length dependence of FFT leakage is documented in
[Sandia SAND1993-0518C](https://www.sandia.gov/research/publications/details/the-use-of-dft-windows-in-signal-to-noise-ratio-and-harmonic-distortion-com-1993-03-01/).

### Inharmonic and independent tones

`harmonic::ToneFitPlan` accepts explicit positive, strictly increasing physical
frequencies and uses the same QR, scaling, residual, admission and buffer owner
as `HarmonicFitPlan`. Stiff-string dispersion and noncommensurate source tones
retain their actual frequencies; no integer harmonic labels are invented.
Close tones retain their relative phase and finite-record cross terms. A short
record that cannot separate them at the requested conditioning floor refuses.
Both interfaces retain the same preparation/reuse complexity and storage
budget. G1/G3 cases cover manufactured stiff-string partials, beating tones
closer than one FFT bin, phase origin, shared numerical results, and unresolved
tone refusal. This is known-frequency estimation, not blind source separation
or a claim that arbitrary close tones can be resolved in noisy records.

`refine_tone_frequencies` and `refine_tone_frequencies_channels` now identify
locally seeded independent frequencies without imposing harmonic spacing.
Every tone has an explicit ordered, separated, sub-Nyquist search interval.
The existing joint DC/sinusoid QR eliminates amplitudes, while a second QR
solves all projected frequency derivatives together, retaining close-tone
cross terms. For vector records, each channel retains independent phasors and
its original common-unit amplitude weight. Streaming augmented derivative and
residual blocks through compact QR factors costs O(E C (N K² + K³)) work and
O(N K + K² + C K) storage; no normal equations or full channel-stack matrix
is formed. The derivation is the projected/Kaufman Gauss--Newton form of the
[O'Leary--Rust variable-projection algorithm](https://www.cs.umd.edu/users/oleary/software/varpro.pdf).
Its omitted residual-dependent Jacobian term does not change the stationary
gradient, but no global convergence or full-Newton claim is made.

`ToneFrequencyBudget` separately gates the sinusoid design and the common-scale
projected frequency Jacobian, and bounds the complete storage and candidate
count. Unresolved/absent tones, active bounds, exhausted work and stagnation
refuse. The convergence criterion is every **untruncated** frequency correction,
not clipped motion or a small residual. Tests in `tests/tone_frequency.rs`
recover beating unisons with a weak inharmonic partial, preserve phase origin
and physical unit scaling, rotate vector components, retain omitted energy,
and exercise numerical-rank, storage, work and domain refusals. This remains
a local stationary-spectrum estimate: it supplies neither global discovery,
tone count, damping, time-varying tracking nor measurement uncertainty.

### Damped tones and double decay

`harmonic::DampedToneFitPlan` fits a constant and the caller's supplied
oscillatory decay poles through the same `FitDesign` Householder QR, full
triangular conditioning, amplitude scaling and allocation-free `fit_into`
owner. Each `DampedTonePole` supplies a positive frequency in hertz and a
nonnegative amplitude decay rate in reciprocal seconds. Its two sampled
columns are `exp(-decay*t) cos(2 pi f t)` and `exp(-decay*t) sin(2 pi f t)`.
Frequency means the damped oscillation frequency, not the undamped natural
frequency. Phasors use `a + i b` under `exp(-i omega t)` and refer to the first
sample, before exponential attenuation; DC and residual RMS retain input units.

Caller pole order is preserved. Equal frequencies with distinct decay rates
can retain double decay, and distinct frequencies retain beating; neither is
replaced by added stationary harmonics. Duplicate poles refuse. The existing
sample-count, strict Nyquist, finite-value, storage and conditioning admissions
also apply. Very close poles or envelopes that vanish before the oscillation
is sampled refuse when the complete design is unresolved. Storage includes
both retained physical parameters for each pole before allocation. Preparation
and reuse have the same asymptotic costs as the stationary fit; zero decay
reproduces the stationary design and fit bits exactly.

The poles are inputs, not inferred parameters. This does not discover modes,
identify damping, certify noisy-record separability, or remove aliased energy.
Omitted modes or incorrect envelopes can bias fitted components and remain in
the sampled residual. `tests/damped_tone_fit.rs` uses an independent damped
oscillator initial-value solution, same-frequency double decay, close beating
tones, time-origin attenuation/phase covariance, stationary-limit bit parity,
physical-unit scaling and preloads, direct residual reconstruction, and typed
rank/domain/storage refusals.

### Joint ringdown pole refinement

`harmonic::refine_damped_tone_poles` and
`refine_damped_tone_poles_channels` locally refine the supplied frequencies
AND amplitude decay rates of `DampedToneFitPlan`. The original Householder
amplitude/DC fit, bounded backtracking and streaming channel-QR owner are
shared with stationary tone refinement. The damped path retains the complete
variable-projection Jacobian: for `r = y - A*c`, the negative residual
derivative is `P_perp*D*c + (A^+)^T*D^T*r`. With `A = Q1*R` it is evaluated
in Q coordinates as `[R^-T D^T r; Q2^T D*c]`, without forming a projector,
pseudoinverse or normal equations. The residual-dependent first block is
retained even away from an exact fit. See
[O'Leary and Rust, §2.3](https://www.cs.umd.edu/users/oleary/software/varpro.pdf).

`DampedToneBounds` supplies a finite closed rectangle for each pole, with
positive sub-Nyquist frequency and nonnegative decay intervals of nonzero
width. Rectangles must be strictly disjoint in at least one coordinate;
separated decay intervals admit same-frequency double decay. Seed, pole,
phasor and correction order is preserved. Frequencies remain damped
oscillation frequencies, not undamped natural frequencies; decay acts on
amplitude, and all peak phasors refer to the first sample. Common-unit
channel weights remain physical under orthogonal changes of vector
components. Constant channels retain their DC without inventing pole data.

`DampedToneRefinementBudget` bounds all candidate evaluations and conservative
numeric storage before allocation. Every untruncated frequency AND decay
correction must satisfy its own tolerance. Trials are limited to a quarter
cycle and half an amplitude e-fold across the record, with monotone residual
backtracking. Active bounds, stagnation and exhausted work refuse. Both the
full Jacobian and its projected physical sensitivity are independently gated
by the declared reciprocal-condition floor, using coordinates `f*duration`
and `decay*duration` without independent column normalization. The reported
condition is their minimum: residual curvature cannot make an unexcited
declared pole identifiable. Projected sensitivity must also exceed a roundoff
floor in globally normalized signal units: a uniformly roundoff-sized
Jacobian can have a good condition number. This floor is invariant to common
physical unit changes and retains physically weak resolved records.
Stationary APIs retain their original Kaufman
derivative path; fixed-pole zero-decay fits retain stationary design/fit bits.

G1/G3 checks in `tests/damped_tone_poles.rs` use an independent damped
oscillator initial-value solution, beating unisons with a weak inharmonic
partial, equal-frequency double decay, physical-unit and first-sample clock
changes, vector-channel rotations, omitted-motion residuals, absent poles,
and work/storage/domain refusals. An in-module finite-difference check
compares the complete profiled Jacobian with independently refitted records
at nonzero residual; a separate omitted-tone case prevents residual-only
curvature from passing the observability gate.

This supplies a ringdown observation building block for
`frankensim-irp-v5-instrument-realism-pexfp.279` (two-polarization decay) and
`frankensim-irp-v5-instrument-realism-pexfp.280` (string loss and late decay).
It does not complete either acceptance bead, certify material parameters,
discover mode count, estimate uncertainty, fit growth/chirps, remove aliasing
or guarantee global/noisy-record identifiability. Unmodeled motion stays in
the residual and can bias the fitted poles. Work remains
O(E C (N K² + K³)), storage O(N K + K² + C K); no new throughput or
mid-factorization cancellation claim is made.

### Linear frequency drift and crossing chirps

`harmonic::ChirpedToneFitPlan` jointly fits DC and constant peak phasors at
supplied linear instantaneous-frequency trajectories. Each `ChirpedTone`
declares the frequency at the first and last samples in hertz. For duration
`T=(N-1)*dt`, its phase in cycles is `f_start*t+(f_end-f_start)*t²/(2*T)`.
Phasors use `exp(-i phase)` and refer to the first sample. Tone order is
retained, including crossing trajectories and shared endpoint frequencies;
duplicate trajectories and numerically unresolved designs refuse.

Both endpoints must be positive and strictly sub-Nyquist. Linearity then
admits the full instantaneous-frequency path, without claiming a bandlimit
or anti-alias filtering. Endpoint coordinates avoid an underflowed Hz/s drift
rate or an overflowed time squared. The integrated phase is evaluated as a
DD convex combination of endpoints times the original DD sample clock, then
reduced in cycles before deterministic trigonometry. Zero drift calls the
original stationary phase owner and preserves its design and fit bits.
All plans share the existing `FitDesign` QR, conditioning, physical scaling
and allocation-free `fit_into` owner.

`refine_chirped_tones` and `refine_chirped_tones_channels` infer both endpoint
frequencies locally inside explicitly supplied disjoint search rectangles.
They extend the original bounded refinement owner, retaining the complete
variable-projection Jacobian and the independent projected-sensitivity
identifiability gate. Coordinates are endpoint frequency times duration;
relative channel amplitudes remain in their common physical units. Constant
channels retain DC without inventing frequency information. No normal
equations, FFT-bin rounding, resampling or piecewise stationary fits are used.
Both untruncated endpoint corrections must meet the frequency tolerance.
Quarter-cycle steps, backtracking, rejected-attempt accounting, active-bound
refusal, rank gates and complete storage admission stay shared with the
stationary/damped solver. Work is O(E C (N K² + K³)); numeric storage is
O(N K + K² + C K).

`tests/chirped_tones.rs` checks independent quadratic-phase records, crossing
trajectories, multi-sensor nodal components, frame/unit changes, record-origin
phase, stationary bit parity, representable endpoint laws with underflowed
Hz/s rates, omitted energy and typed rank/domain/work refusals. In-module
checks exercise a billion-sample phase clock without allocating its record,
and compare the full profiled Jacobian with independent refits at nonzero
residual. These are numerical checks, not measurement-uncertainty coverage.
This local linear-law owner does not discover tones, infer arbitrary tracks
or envelopes, extrapolate a chirp, or reinterpret changing frequencies as
stationary modal poles. It supplies a primitive for chirp analysis; it does
not establish an Euler-disc mechanism or close its terminal-asymptotic bead.

### Joint amplitude decay and frequency drift

`harmonic::DampedChirpedToneFitPlan` fits DC and initial peak phasors at
supplied linear frequency trajectories with exponential amplitude decay.
`DampedChirpedTone` retains the first- and last-sampled frequencies in hertz
and the nonnegative amplitude decay rate in reciprocal seconds. The basis is
`exp(-decay*t) * (cos(phase), sin(phase))`, using the original DD endpoint
phase law over `T=(N-1)*dt` and the original deterministic decay exponential.
Phasors use `exp(-i phase)` and refer to the first sample before attenuation.
This shares the original `FitDesign` QR, scaling, conditioning and
allocation-free transactional `fit_into` owner. Zero decay preserves the
supplied chirp plan's fit bits; zero drift preserves the supplied damped-tone
plan's fit bits. Both frequency endpoints remain strictly positive and
sub-Nyquist. Equal frequency laws with distinct decays remain separate when
the complete sampled design resolves them; duplicate triples refuse.

`refine_damped_chirped_tones` and `_channels` identify both endpoint
frequencies AND amplitude decay in one variable-projection solve. The shared
bounded optimizer admits three nonlinear coordinates per tone, expressed as
`[f_start*T, f_end*T, decay*T]`. It retains the complete residual-dependent
Jacobian, independent projected-sensitivity rank and amplitude floors,
common-unit channel weights, monotone backtracking and rejected-attempt
accounting. Every untruncated frequency and decay correction must satisfy its
physical tolerance; active bounds are not convergence. Endpoint corrections
remain quarter-cycle limited and decay corrections half-e-fold limited.
All derivative centering removes only fitted basis columns; the actual
envelope remains `exp(-decay*t) <= 1` and is not centered into a growing
exponential. Disjoint three-coordinate search boxes retain original tone
identity, including components distinguished only by decay.

Storage admission covers all three retained parameters, candidate conversions,
all three derivative columns per tone and both full/projected channel QR
folds. Joint residual dimensions must support all `3*K` nonlinear coordinates
after the `1+2*K` linear amplitudes/DC columns. Work remains
O(E C (N K² + K³)), numeric storage O(N K + K² + C K). No normal equations,
frequency/decay alternating approximation or new cancellation claim is added.

`tests/damped_chirped_tones.rs` exercises decaying crossing trajectories,
same-frequency-law double decay, nodal channel rotations and unit scaling,
record reslicing with phase AND attenuation, exact supplied-law limits,
endpoint laws at clocks where Hz/s underflows, and typed work/storage/rank/
domain refusals. An in-module finite-difference test checks all three full
profiled derivatives against independent refits with substantial omitted
motion. These are local numerical observations supporting instrument-plan
§7.1/7.2 and the decay observations needed by `pexfp.279`/`pexfp.280`;
they do not close those specimen/physics acceptance beads. Nonlinear drift,
nonexponential envelopes and omitted tones remain model mismatch that can
bias fitted parameters. This does not discover modes, identify material
loss, provide uncertainty coverage, remove aliasing or turn changing
frequencies into stationary modal poles.

### Local fundamental-frequency refinement

`harmonic::refine_harmonic_frequency` refines an independently supplied pitch
seed inside a positive frequency interval. Every candidate jointly refits DC
and harmonic amplitudes through the same Householder QR owner. Its analytic
sinusoid derivative is projected outside the fitted columns; a scalar projected
Gauss--Newton correction uses monotone backtracking and a quarter-cycle step
cap across the record. This follows the variable-projection separation of
linear and nonlinear parameters described by
[O'Leary and Rust](https://www.cs.umd.edu/users/oleary/software/varpro.pdf).

The explicit budget bounds all candidate QR evaluations and conservative peak
numeric storage, including refinement scratch. Ill-conditioned or unresolved
trial candidates consume an evaluation and backtrack to a shorter step; an
inadmissible seed refuses immediately. Exhaustion reports the last accepted
iterate. Nyquist admission covers the entire search interval. Normalized fitting
and search arithmetic preserve physical unit scaling. Missing fundamentals can
be inferred from the declared upper harmonic pattern; no synthetic fundamental
amplitude is introduced.
An active bound, unidentifiable frequency, exhausted work, or floating-point
stagnation produces a typed refusal rather than a nominally converged spectrum.

This offline local estimator costs O(E (N K^2 + K^3)) for at most E fits and
O(N K) storage. It does not discover pitch octaves, certify uncertainty, track
time-varying frequency, or remove inharmonic energy. Frequency tolerance is a
numerical correction criterion. Tests recover a missing-fundamental record
whose coarse pitch hid weak harmonics, preserve phase/DC across unit scales,
retain an omitted inharmonic residual, and enforce storage, search and work
limits. Same-build deterministic arithmetic is retained; no new cross-ISA or
mid-factorization cancellation claim is made.


### Common-frequency vector records

`harmonic::refine_harmonic_frequency_channels` fits simultaneous channels on
one physical sample clock with a shared candidate frequency and QR. DC and
phasors remain independent per channel; the frequency minimizes their total
squared residual in common physical units. Original amplitude weights are
retained, so rotating Cartesian force components preserves the Euclidean
objective. Constant channels retain their DC without supplying frequency
information. This neither averages independently estimated pitches nor gives
a weak channel an invented unit-amplitude weight.

Channels stream through two shared sample buffers; all retained candidate
fits and channel scales are included in storage admission. For C channels,
work is O(E (N K^2 + C N K + K^3)) and storage O(N K + C K + N).
The scalar API delegates to this same owner. Tests check scalar identity,
orthogonal frame changes and phasor rotation, amplitude weighting, common
unit changes, constants, channel shape and complete storage refusal. Inputs
with unrelated units or noise metrics require explicit caller conversion;
the local stationary and uncertainty boundaries above remain unchanged.

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
`tests/frankenscipy_fft_oracle_casebook.rs` adds the test-only external
FrankenScipy transform seam. It binds fs-fft's forward complex and real
power-of-two transforms to fsci-fft's independent naive DFT under explicit
strict, backward-normalized, single-worker options. Analytic impulse evidence,
fixed seeded finite fixtures, same-run replay, canonical input frames,
observed output-bit digests, and a disclosed exact-reference corruption are
separate structured records. `fsci-fft` is a development dependency only;
this tranche is central package-proof pending.
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
- The serial N-D transform is separable and gathers each pencil into reusable
  line buffers rather than blocked transposes. Legacy calls allocate the two
  buffers once per call; prepared `FftNdScratch` calls retain them across calls.
- N-D real transforms (r2c/c2r with only the first axis half-length) are not
  shipped; N-D is complex-in/complex-out. Callers pack real fields as `C64`.
- `Fft`, packed real transforms and N-D plans retain their power-of-two
  admission. `GeneralFft` supplies arbitrary-length 1D complex transforms;
  arbitrary-length packed-real, DCT and N-D paths are not supplied.
- The FrankenScipy Casebook is agreement evidence for its disclosed bounded
  fixtures, not a forward-error certificate. It does not establish Python
  SciPy execution, arbitrary-size, inverse/c2r, N-D, or DCT equivalence,
  performance parity, or fresh cross-ISA bit identity. FrankenScipy's naive
  oracle uses platform trigonometry, so its observed output digests are
  same-run diagnostics rather than cross-ISA goldens.
- `C64` is not a public complex-arithmetic library; only what the FFT needs.

## Extended harmonic field contractions

Deterministic rectangular observations now retain double-double carrier
rotations, including their small cosine decrements. The Hz-authored
`mean_phase_hz`, `product_kernel_hz`, `proper_product_kernel_hz` and
`phase_hz_from_clock_dd` methods reduce cycle products before multiplication by
2*pi. Scalar high/low products and source-clock components are reduced
separately, so neither the irrational conversion error nor adding a subperiod
offset to a huge carrier phase creates a false beat. Sum/difference frequencies
are formed in Hz before the centered Gram calculation; sinc numerators use the
same cycle reduction. The small-width sinc-1 and disparate-frequency divided
difference retain their original cancellation-stable formulas. Signed
frequencies, DC and uniform/triangular sensitivity remain supported; APIs and
radian-authored clocks are unchanged. This remains bounded, allocation-free DD
arithmetic, not an arbitrary-null relative-error or interval certificate.

G1/G3 regressions compare an optical-carrier quadrature pair against its
independent 1 Hz real beat at large epochs, preserve subperiod exposures under
complete-cycle epoch translations, and retain a DD quarter turn on a huge
integer-cycle source clock. Complete-cycle exposure zeros are also checked.

Both prepared product types expose `mean_real_product_parts_dd`, accepting
[real, imaginary] double-double field components. Coalesced low components
therefore survive later products without expanding every contraction into four
C64 calls. Existing C64 entry points delegate to the same arithmetic.
Contractions remain allocation-free and reuse cached trigonometry. Transcendental
leaves remain deterministic f64 evaluations, not certified enclosures or a
relative accuracy guarantee at every null.

Four focused regressions cover real-signal quadrature, a resolved high-carrier
beat, a fractional carrier on its original clock, retained low amplitudes and
cosine decrements, and invalid inputs. The complete 16-test module passes
native Rust type analysis (`rustc --test -Z no-codegen`); these checks do not
execute tests or constitute a full workspace build.


## Normalized shaped detector exposures

`HarmonicWindow::with_profile` admits explicit `HarmonicWindowProfile::Uniform`
or `Triangular` sensitivity on the original source clock. `new` retains the
uniform profile and its existing arithmetic. The profile is retained by the
Copy window and returned by `profile()`. Start and duration remain separate.

The triangular sensitivity is `2*(1-|2*u-1|)/T` for
`u=(t-start)/T` in [0,1], with zero outside the interval. It is nonnegative,
has unit integral, vanishes at both endpoints and peaks at 2/T at the midpoint.
The same single sensitivity weights the instantaneous observable: tapering
both input fields would square the weight and describe a different detector.
DC normalization and proper equal-frequency half-power normalization are
unchanged.

All linear, deterministic and proper kernels, including Hz-authored paths,
honor the profile. The centered triangular phase transform is
`sinc(omega*T/4)^2`. Deterministic cosine/cosine and sine/sine weights use the
convolution of two uniform clocks with half the support: from their Gram
entries C,S, the triangular entries are C^2+S^2 and 2*C*S. This preserves
subperiod and disparate-frequency sine energy without differencing nearly
equal transforms. Original DD phase rotations, small sinc decrements and
retained-amplitude contractions remain in use. Preparation and contraction
have bounded constant work, no quadrature, no sampling grid and no allocation.
Domain owners retain their original line/work limits and cancellation.

These are detector-weighted averages of fields, products, heat or storage
rates. Multiplying by T gives a sensitivity-weighted exposure. It gives the
unweighted energy/impulse only for Uniform. For Triangular, the average of a
local storage derivative is generally not the endpoint storage change divided
by T. A pointwise physical balance remains valid under the common weight;
an unweighted endpoint balance must use Uniform. No sampled electronic filter,
detector bandwidth calibration, new source statistics or arbitrary-null
relative certificate is supplied.

Six focused G0/G1/G3 regressions cover normalized DC and unchanged defaults,
weighted real-field and exact four-phase proper-ensemble quadrature, tiny sine
energy, vastly disparate frequencies, sub-ulp destructive beats, analytic
transform zeros, clock translation, high-carrier Hz beats and observable
weighting versus endpoint storage. Native formatting and type analysis
(`rustc --test -Z no-codegen`) passed for the complete 22-test module; native
test execution remains pending in the storage-exhausted authoring environment.
An independent standard-library Python reconstruction of four deterministic
fixtures and their four-phase ensembles agreed with weighted real-time Simpson
integration to maximum absolute errors 4.45e-15 and 2.23e-16 respectively.
Those are reference calculations, not executed Rust results.

## Time-dependent quadratic product spectra

`harmonic_window::PreparedHarmonicProductSpectrum` prepares a finite analytic
line spectrum from explicit `HarmonicProductInput` bilinear contractions.
Peak amplitudes use exp(-i omega t); signed input frequencies are supported.
For `HarmonicProductLaw::Real`, each product retains one half of A*conj(B)
at the difference frequency and one half of A*B at the sum frequency. Real
also describes the expected product sharing a zero-mean real unit-variance
coefficient. ProperComplex assumes a shared coefficient with E[z]=E[z^2]=0
and E[|z|^2]=1 and retains only the difference term. Cross-factor products
may be omitted only when the caller's authored joint moments make them zero.

Output consists of a DD DC value and positive-frequency peak tones.
Negative-frequency terms are conjugated. Frequencies retain the DD sum or
difference of the supplied f64 Hz values, including low sum-frequency bits.
Exactly equal DD frequencies coalesce in authored order with DD amplitude
sums; no frequency tolerance or amplitude threshold drops physical beats.
Callers must combine coherent source fields before forming products.
Preparation takes O(K log K) work and O(K) storage for K generated terms.

Explicit limits cover supplied products, pre-coalescing generated terms and
the logical vector payload at its admitted capacity; inputs, iterator state,
stack and allocator overhead are excluded. A fixed 131072 generated-term
ceiling bounds the synchronous sort, with cancellation before and after it.
Iteration and coalescing poll cancellation. Domain, budget, allocation or
cancellation refusal publishes no partial spectrum.

Scalar-time evaluation and Uniform/Triangular detector-window means allocate
nothing and cost O(distinct output tones). They preserve the source clock,
DD frequency and amplitude parts and the existing shared detector kernel.
The result is an actual deterministic or expected quadratic waveform,
not PSD/Hz, a sampled FFT, detector bandwidth, random realizations or the
fourth-moment noise spectrum of the observable. No null-relative certificate
or stronger source/propagation accuracy is claimed.

Five focused tests cover direct real-time products, the exact four-phase
proper ensemble, both detector profiles, signed frequencies, retained
high-carrier sum bits, weak amplitudes, admission and cancellation.
The exact shared window and spectrum modules passed all 27 of their tests
in a focused native Rust harness using the real fs_math library. This scoped
execution does not substitute for the full repository build and quality gate.

### Original angular clocks in quadratic spectra

AngularHarmonicProductInput and from_angular_frequencies retain source-authored radians/second without rounding each carrier through f64 Hz. The original Hz constructors remain available. One shared preparation forms DD sum/difference frequencies and coalesces exactly in the authored unit. Each tone retains DD Hz and angular frequency; angular_frequency_dd exposes the value used directly by temporal kernels. The logical tone payload therefore includes both frequency coordinates. This preserves narrow beats and sum-frequency low parts on large angular carriers. Both frequency coordinates must remain representable. Two regressions cover an original 0.25 rad/s beat at a 2^50 rad/s carrier, a half-ULP angular sum, signed frequencies, both detector profiles and cancellation/domain refusal.


## Causal integration of interference loads

The prepared quadratic spectrum exposes `causal_integral` and a cancellable
variant. For start t0, duration h and nonnegative decay lambda they return
both integral_0^h f(t0+s) ds and
integral_0^h exp[-lambda*(h-s)] f(t0+s) ds in observable-units * seconds.
The result also retains exp(-lambda*h) and its one-minus decrement in DD.
The domain adapter supplies physical units, initial state and relaxation;
the generic kernel does not infer heat from a signed acoustic or EM flux.

Source start and duration remain separate. Native angular frequencies and
their DD low components survive phase construction, including an interval
too short to change the rounded absolute endpoint. All DC, sum and beat terms
participate. Zero decay makes the two integrals identical; zero duration
returns zero input with unit decay. A bounded local divided-difference series
retains small time/frequency/decay corrections; scaled complex quotients avoid
squaring large frequencies. No sampled time step, aliasing cutoff or per-tone
allocation is introduced. Cost is O(distinct tones), with cancellation at
bounded blocks and explicit refusal for invalid or unrepresentable arithmetic.
The result inherits fs_math deterministic transcendental accuracy and finite
DD dynamic range; it is not an interval or relative-error certificate at nulls.

Six new causal regressions cover constant loads and stiff/short decay,
independent real-time product quadrature, accepted-interval composition,
weak sine energy beside cancelling DC/carrier terms, native angular beats
and clocks with unrepresentable rounded endpoints, and transactional refusal.
Together with the two new angular-admission regressions and existing
temporal/spectrum tests, all 35 definitions compiled, linked and ran in a
focused native Rust executable using the real fs_math library, with no
nominal physics boundaries. This is not the full workspace quality gate.


## Spatially contracted spectra and resolved nonlinear forcing

`AngularHarmonicSpectrumTerm` admits an already-contracted REAL observable
tone in the existing negative-time convention, with a DD angular frequency
and DD peak amplitude. Its amplitude already includes physical quadrature
weights and product normalization. `from_angular_terms[_cancellable]` retains
DC, conjugates negative frequencies and uses the same deterministic coalescing
path as the original field-product constructors. Native DD pairs are normalized
before grouping. This lets domain adapters contract physical space first,
without storing a separate tone for every quadrature point. Explicit term
and vector-payload limits, the existing fixed tone ceiling and cancellation
apply before publication. No statistics, units, spatial accuracy or dissipation
are inferred from an already-contracted term.

The causal integral also accepts a DD source-clock start through
`causal_integral_from_clock_dd[_cancellable]`. A source epoch plus accepted
elapsed time therefore remains resolved when a single f64 absolute timestamp
would erase later substeps. The f64 API delegates to the same implementation;
input normalization, validation, decay and physical-integral semantics agree.

`input_ripple_bound[_cancellable]` bounds how far cumulative input can depart
from the endpoint chord when the source is replaced by its exact interval mean.
For F(s)=integral_0^s f(t0+u)du on [0,h], it returns B with
|F(s)-(s/h)F(h)| <= B independently of t0. DC contributes zero. For each tone,
a conservative amplitude envelope A gives the minimum of A*w*h^2/8, A*h/2
and 2*A/w. The first is the linear interpolation remainder of the primitive
(NIST DLMF https://dlmf.nist.gov/3.3.E5); the other two follow from the bounded
input and bounded harmonic primitive. This detects unresolved energy excursions
even when the net input over a complete carrier cycle is zero.

Frequency intervals and every positive arithmetic operation are rounded
outward with IEEE f64 predecessor/successor operations. The four absolute
DD amplitude components are scaled before summation to avoid premature
envelope overflow. Positive underflow is rounded upward, while exact darkness
remains zero. The result encloses the retained analytic spectrum's input
excursion, not source-model, propagation, spatial-quadrature or nonlinear
state-solver errors. A domain integrator must still control its thermal/mechanical
integration error and physical validity separately. All queries are O(distinct
tones), allocation-free and cancellable; an unrepresentable bound refuses.

Eight new tests cover contracted signed/DC/weak and DD-frequency terms,
atomic admission, late DD-clock substeps and partition composition, direct
primitive samples, full-cycle alias witnesses, small-window refinement,
high-carrier suppression, weak cancellation, extreme finite amplitude scaling
and cancellation. All 43 shared temporal/spectral/causal tests compiled,
linked and executed with the real fs_math library and no nominal physics
boundaries. This scoped execution is not the full workspace quality gate.

## Physical-unit scaling and balanced causal integration

`PreparedHarmonicProductSpectrum::causal_integral_from_clock_dd_scaled` and
its cancellable counterpart apply a finite positive output-unit conversion to
DC/tone amplitudes before the existing allocation-free causal kernels.
Returned integrals use converted-output-units times seconds; decay and
one-minus-decay factors are unchanged. Existing DD-clock methods delegate with
scale one. This lets a tiny authored observable expressed in large physical
units retain a representable physical input integral without first forming
an underflowed unscaled quantity.

Harmonic contributions prefer physical time gains: duration multiplies the
normalized kernel before the source amplitude. Their natural scale for long
or stiff intervals is inverse frequency or inverse decay, so a small
normalized kernel does not erase a finite accumulated contribution.
If that ordering is unrepresentable for a very short interval, the complete
tone falls back to the original ordering: combine the normalized real and
imaginary source products in DD, then apply duration. Exact-zero amplitude
components skip unnecessary gain formation.

Both orderings check their multiplicative steps. Nonzero products that become
fully zero, and nonfinite products, refuse when no admitted ordering remains.
Genuine signed cancellation between representable contributions stays valid.
Early physical-amplitude overflow or underflow can still refuse even when a
different multiplication order outside this API could represent a later
integral; this API does not claim an expanded DD exponent range. The existing
phase, exponential, divided-difference and local-series kernels are unchanged.
No stronger cancellation-relative or interval accuracy claim is introduced.

Focused causal regressions cover DC and harmonic unit restoration, signed
phase, ordinary and exponentially weighted input, explicit underflow/overflow
refusals, signed cancellation, identity scaling, domains and cancellation.
Long-window tests retain an amplitude of 1e-200 over 1e200 seconds and a stiff
weighted input beside a finite ordinary integral. The short-window complex
source case exercises the normalized-projection fallback.

The native nonlinear enthalpy consumer takes the already-converted integral
directly as source joules, including input below its original f64 phase-state
resolution; any unapplied part remains visible in the numerical energy
residual instead of being relabeled as a dark source.

All 49 shared regressions compiled, linked and executed with real fs_math in
the same native run as the 13 original-owner nonlinear thermal regressions:
62 passed, 0 failed. Six new causal tests cover the unit and arithmetic-range
cases above. The shared math is exercised directly, without nominal numerical
or constitutive substitutes. Native rustfmt passed; the scoped executable is
not the full repository quality gate.

### DD-clock observable endpoints and storage increments

`PreparedHarmonicProductSpectrum::evaluate_from_clock_dd` retains the original
source epoch and low result components. `increment_from_clock_dd` computes
f(start+duration)-f(start) tone by tone, using separate source and half-duration
phase rotations and -2i*sin(omega*duration/2). DC cancels algebraically before
summation; no rounded endpoint subtraction or f64 start+duration is formed.
The existing f64 evaluator shares the same endpoint implementation.

Duration must be finite and nonnegative; zero duration yields exact zero.
Nonfinite or unrepresentable phases, including a nonzero phase increment lost
to underflow, refuse. Work is O(distinct tones), allocation-free and cancellable.
The original deterministic transcendental and DD finite-range limits remain;
this is not an interval or near-null relative-error certificate. An endpoint
change is neither a normalized detector average nor an input time integral.

Six native tests pass against the complete production harmonic-window/product
modules: independent multitone endpoint values, short intervals beside a huge
DC and unresolved f64 epoch, DD-clock subdivision, retained low storage and
f64 continuation, zero/invalid/underflow cases, and cancellation/replay.
The finite-volume acoustic observer consumes these endpoints and changes in
its separate surface/volume conservation tests. Full workspace validation was
not executed in the delivery environment.

## Recorded multichannel coherence: joint Welch cross spectra

`welch::WelchCrossSpectrum` estimates synchronized real channels using the
original deterministic FFT. Every channel uses the same sample clock,
segment boundaries, detrending and rectangular or periodic Hann window.
Separate channel slices and existing interleaved Cartesian records share the
same owner; the latter requires no copied or transposed input history.

Each positive-frequency bin retains one complex channel vector per segment.
Their outer-product sum is the one-sided cross-spectral density
`S_ab=sum_r B_ar conj(B_br)` in input-unit squared per hertz. The factors use
the physical `exp(-i omega t)` convention, preserve complex cross-channel
phase, and include the actual window energy, segment average, sample interval
and the interior-bin factor of two. DC and Nyquist are not doubled. This is
Welch's averaged modified periodogram, not a coherent-gain tone fit:
<https://doi.org/10.1109/TAU.1967.1161901>.

Keeping factors avoids a dense cross-spectral matrix and any numerical
eigenfactorization. Storage is O(bins * segments * channels); the caller must
admit this retained rank explicitly. Channel/sample, FFT-count and storage
budgets are checked before allocation; storage accounting includes retained
factors, FFT data/scratch, window, recursive plan tables and six-step temporary
buffers. Any exact segment length in 4..=65536 is admitted. Power-of-two
segments retain the original direct FFT arithmetic and hot path; other lengths
reuse the cached `GeneralFft` Bluestein plan and caller scratch. Its internal
convolution padding does not change the physical window, sample clock or bin
spacing `1/(N*dt)`. The complete chirp/convolution workspace is charged before
reading samples. `maximum_transforms` counts physical channel/segment DFTs;
non-power-of-two DFTs use two convolution FFTs of length at most 131072.
Cancellation polls input admission and bounded work between channel DFTs.
Incomplete tails are reported and still checked for nonfinite values; no
physical samples are padded or interpolated.

Only even segment lengths have a Nyquist endpoint. The highest positive bin
of an odd segment retains both complex components and the interior-bin factor
of two, including when borrowed for physical peak-phasor propagation. Focused
regressions compare odd/prime/composite cross spectra and integrated power to
direct DFT and Parseval oracles, preserve highest-bin quadrature/cancellation,
and check batch/stream identity, cancellation retry and convolution budgets.

Per-segment scaling and double-double reductions protect representable
fluctuations, common offsets and spectral products. Nonfinite input, lost
nonzero normalized components and unrepresentable factors refuse explicitly.
These checks do not turn finite precision into an interval or arbitrary
relative-error certificate near a cancellation null.

One segment has rank one and cannot establish a source's physical coherence.
Overlapping segments are algebraic factors, not a declaration of independent
realizations. No degrees of freedom, confidence interval, stationarity,
cross-frequency covariance or unresolved-band power is inferred. The caller
retains responsibility for record duration, leakage and spectral resolution.

Focused native regressions use the original FFT/SIMD/math implementation and
compare all cross-spectral entries with an independent direct DFT. They cover
complex phase, rectangular/Hann windows, detrending, Parseval with endpoint
weights, coherent and partially correlated channels, layout equivalence,
numerical range, resource refusal, cancellation and replay. This scoped
execution is not the full workspace or distributed repository gate.

### Coherent transfer of recorded channels

`WelchCrossSpectrum::transfer_channels` applies an explicitly supplied complex
frequency-dependent map `y=H(f)x` to every retained segment factor before any
power contraction. The resulting joint spectrum has `S_out=H S_in H*`; its
clock, window, segment count, bin width, tail accounting and one-sided
normalization remain the original record's. Existing physical peak-bin
consumers can propagate the result without an additional FFT or another
peak/RMS or bandwidth conversion. Channel order, source identity, physical
units, calibration and transfer validity remain explicit caller obligations.

The transfer uses the existing negative-time phasor convention without an
implicit conjugation: a delayed signal has gain `exp(+i omega tau)`, while a
readout written `w* x` supplies conjugated weights. DC and even-length Nyquist
coefficients must be real; the highest positive odd-length bin retains its
complex transfer. This operates on the estimated spectral bins. It does not
reconstruct a causal finite-record filter, remove taper leakage or edge
transients, or infer cross-frequency covariance and statistical confidence.

Complex products and complete channel contractions remain double-double until
each output factor is rounded once, protecting weak residuals through strong
cancellation. Coefficients are evaluated once per bin/output/input and reused
across all segments. One frequency matrix is cached; no dense CSD or
eigenfactorization is constructed. `channel_transfer_requirements` reports
exact logical coefficient/product counts, output factors and numerical storage.
The complete output plus reusable matrix is admitted before callbacks or
allocation. Input and callback-owned storage and allocator overhead are
excluded. Bounded cancellation, callback errors, nonfinite coefficients,
endpoint violations and unrepresentable factors publish no partial result and
leave the input untouched.

Six focused regression groups compare real calibration against independently
mixed original records, complex maps against dense `H S H*`, compensated
delayed-channel reinforcement/cancellation, odd-bin peak-power normalization,
quiet fields and exact product residuals, preflight work/storage refusal,
coefficient caching and transactional cancellation/replay. These are numerical
regressions, not an arbitrary-null relative-error certificate.

### Reference-assisted interference transfer and held-out residuals

`WelchCrossSpectrum::fit_reference_transfer` learns a complex per-bin map
from explicitly selected reference channels to one or more separate targets.
It minimizes the segment-factor residual of
`target = sum(gain * reference)` in the original negative-time convention,
with no implicit gain conjugation. The immutable `WelchReferenceTransfer`
retains DD gains, original channel order, numerical ranks and training segment
count. It can be applied unchanged to a compatible held-out record through
`residual_spectrum`; that record may contain a different number of segments.

The fit uses pivoted modified Gram-Schmidt with repeated projection, directly
on the original factor columns. Each nonzero reference is normalized before
rank selection, so changing reference units does not discard an independent
weak channel. The caller supplies a relative residual-norm tolerance in (0,1);
it controls numerical rank, not a physical noise floor. Zero and dependent
channels contribute no new direction. Rank-deficient cases return a
deterministic basic solution, without a minimum-norm or source-identifiability
claim. Neither a channel covariance matrix nor its inverse is constructed.

Application computes `target - sum(gain * reference)` on every segment
factor before any power contraction. DD products and sums retain small
residuals which subtraction of rounded predicted/total power would erase.
All output cross spectra therefore retain the original PSD factorization,
Welch normalization, bin widths, clock, endpoint reality and tail accounting.
The same physical peak-bin owners can observe those residuals. Training and
application share the exact channel layout, Welch configuration and sample
interval; channel identities, units and physical interpretation remain caller
responsibilities.

`reference_fit_requirements` reports the full-rank conservative work bound,
retained coefficients and peak reusable QR storage. For M references, T
targets and S segments, scratch is O(M*S + M²), retained gains are
O(bins*T*M), and fitting work is O(bins*(M²*S + T*M*S + T*M²)).
The complete output/model and scratch are admitted before allocation.
Application admits only its new factors and reference products and performs
no new FFT or refit. Cancellation polls at most every 256 vector/algebra
entries. Invalid sets, incompatible records, arithmetic range loss and budget
refusals publish no partial result. Across-segment normalization refuses to
erase a nonzero component outside its representable range; an absent
reference subspace leaves even a wide-range target untouched.

This is a finite-record linear prediction/measurement operation. It does not
establish that correlated content is unwanted, recover physical independent
sources, implement a causal time-domain filter, deconvolve taper leakage,
supply cross-frequency covariance or confer statistical confidence.
Same-record fitting can overfit; a one-segment fit can cancel its training
factor while failing on independent data. Model selection and physical
stationarity must be justified separately.

Seven focused new tests and all 26 existing Welch/transfer tests pass using
the original current FFT core, exact-length owner, streaming owner and
SIMD/math dependencies. They cover real-record complex-gain recovery and
held-out multichannel cross spectra, dependent/zero/rescaled references,
quiet residual power, nearly parallel references and rank choice, endpoint
phase laws, single-segment overfit, work/storage admission, range refusal,
cancellation and exact replay. This native execution does not replace the
unavailable full Cargo/DSR workspace gate.

### Shared physical-bin admission for propagation owners

`peak_bin_hz` and `peak_bin_rad_s` borrow a `WelchPeakBin` without copying its
retained factors. They admit only positive interior bins at exactly the
prepared physical frequency. The angular-frequency entry point compares
`TAU*bin_hz` directly, avoiding a lossy inverse conversion back to hertz.
`peak_factor` applies `sqrt(2)*sqrt(df)` once and refuses nonfinite, overflowing
or silently underflowing physical source components. The resulting peak
factors feed the original half-weighted harmonic statistics. Raw borrowed
factors remain explicitly in source units per square root hertz.

Physical consumers still own source identity/order, units, material, work and
cancellation admission; this view invents none of those. The original 2D
Curle consumer now uses this shared conversion, which also supports the 3D
and BEM propagation owners. Eleven native estimator tests and the four
original physical 2D consumer tests pass after rebuilding/relinking; new
coverage includes exact angular-frequency matching and representable
normal/subnormal conversion alongside explicit range refusals.

`HarmonicProduct::centered_gram_dd` exposes the existing centered cosine/cosine
and sine/sine double-double weights without changing their numerical owner.
Joint spatial/temporal observers can compose independent symmetric windows
and apply one combined center phase, preserving weak translated interference
nulls instead of separately rounding spatial and temporal rotations. Callers
still own the original center phase and source clock; these weights alone
are not the translated product or an additional phase-average assumption.

### Streaming synchronized interference measurements

`welch::WelchStream` acquires a declared record through synchronized channel
slices, interleaved chunks or borrowed dynamic layouts through `push_with`.
Its `(channel, sample)` accessor must present the same original value on every
read, retaining channel identity and the declared clock without a transpose.
Construction admits total channel/sample/FFT counts, the final factor bank,
original FFT/window workspace and two
segment-sized input buffers. Batch and streaming share one segment algorithm
and the final segment count: chunk boundaries preserve factor bits, overlap,
complex cross-channel phase, detrending and one-sided PSD normalization.
Completed input is discarded; chunks are borrowed without a record-sized copy,
and small chunks append without repeatedly copying the entire overlap.

Each push commits accepted samples/segments only after complete validation,
transforms, tail preparation and the final cancellation poll. Invalid input,
numerical refusal and cancellation leave accepted history unchanged. Retry
overwrites unaccepted factor slots, including zero channels. `finish` requires
the complete declared record and moves the existing immutable spectrum out
without a second factor bank, transpose or padding. Retained-factor storage
remains O(bins * segments * channels), with O(segment_samples * channels)
input storage. This is bounded-record acquisition, not open-ended snapshots,
inferred segment independence or cross-frequency coherence.

Five new native regressions and eleven existing Welch regressions pass against
the production FFT/SIMD/math owners. Coverage includes bitwise chunk-partition
equivalence across hop/window/detrend/layout choices, late cancellation and
changed-data retry, exact overlap/tail counts, extreme scaling, input/budget
refusal and incomplete finalization. Full workspace gates remain unrun.

## Direct channel access for joint Welch observations

`WelchCrossSpectrum::estimate_with` publicly exposes the existing callback-based
estimator used by channel-slice and Cartesian observations. A caller can present
a fixed record in its native layout through `(channel, sample)` access without
allocating a channel transpose or histories for identically zero channels. The
callback must present consistent samples from that record. All existing finite-
input, clock, Welch-configuration, work/memory budget, cancellation, window,
detrending, tail and factor-normalization rules remain unchanged; the estimator's
arithmetic is unchanged. The distributed three-axis jet consumer uses this seam
with source-major `[0,Fx,Fy,Fz]` channels and the original physical sample clock.


## Reusable cancellable N-dimensional pencils

`FftNd::scratch` creates two longest-axis buffers retained by the caller.
`forward_with_scratch_cancellable` and `inverse_with_scratch_cancellable`
reuse those buffers and return the caller's unchanged interruption reason.
Checkpoints occur at entry/exit and before/after each 1D pencil. The individual
1D transform is the smallest cancellation unit; bounded-latency consumers
must bound axis length.

Legacy serial forward/inverse use the same pencil loop with an infallible
callback. Axis order, gather/scatter order, butterfly arithmetic and inverse
normalization are unchanged. Scratch from a sufficiently large plan can be
reused. Cancellation can leave transform data partially updated, so consumers
that require atomic publication must stage their data; scratch can be reused
after a refused attempt.

The separable driver allocates nothing during a prepared call. The existing
optional six-step 1D route can still allocate internally for sufficiently large
axes; this interface makes no broader allocation or performance guarantee.
Independent rectangular DFT/round-trip and interruption/retry fixtures cover
the shared path. The hardness-contact consumer uses bounded padded axes below
that six-step threshold and retains its own original physical residual checks.

## Positive short-band uniform-exposure interference

`HarmonicWindow::prepare_uniform_band_hz` provides a cached positive temporal
factor for a coherent analytic field, using the original authored frequencies
in hertz and the existing double-double shutter clock. It removes only a common
carrier before phase reduction. Each retained Legendre mode is projected from
the complex source amplitudes before any power is formed, so a weak multitone
null does not depend on subtracting rounded bright pair powers. The companion
complex cross moment uses the same modes for mutual coherence. Independent
source banks remain separate caller-owned contractions.

Admission is deliberately bounded: a uniform window, 3 through 64 tones, and
`pi * duration * (max_hz - min_hz) / 2 <= 0.5`. Unsupported bands return
`None` so existing pair kernels retain their original role. The allocation-free
`uniform_band_mode_count_hz` preflight reports the retained mode count before
preparation. The factor uses at most 17 modes, O(tones * modes) retained storage,
and fixed observation scratch. Logical work is `(tones + 1) * modes` for power,
`(2 * tones + 1) * modes` for a complex cross moment, and
`(2 * tones + 4) * modes` for a two-axis coherency. Each projection has at most
1088 terms; consumers own work admission and cancellation checks around this
bounded leaf.

`project_parts_dd` returns unweighted coefficients, omitting each mode's common
`(-i)^mode` phase only because it cancels from second moments.
`mode_cross_parts_dd` applies the complete integer weight
`2 * mode + 1` and complex conjugation. Scalar powers, complex crosses and
Jones outer products therefore share the same contraction. The weight is
applied before scalar multiplication can underflow, scaling the smaller
operand to avoid unnecessary overflow. These are analytic-field moments:
there is no implicit real-peak factor of one half. Shape, index, nonfinite-input
and overflowing intermediate/result values return explicit errors. This path
does not add a new underflow-refusal policy.

The exact-arithmetic power in omitted Legendre modes is bounded by
`5e-50 * (sum_i |A_i|)^2` for admitted bands. This is only the modal-tail bound.
The separate 17-term small-argument spherical-Bessel series, floating
arithmetic, phase reduction and cancellation within individual modal
coefficients are not included in it. No relative-accuracy certificate at an
arbitrary null, general shutter approximation or change to the two-tone path
is claimed.

Five focused regressions and all 67 existing harmonic-window/spectrum tests
pass in native execution against the actual math and harmonic owners. They
cover the positive fourth-order three-tone null down to a 1e-12-second shutter,
independent quadrature near the admission boundary, complex Hermitian and
shared-modal identities, large-epoch/common-carrier invariance, work counts,
large amplitudes, weighted subnormal power, boundary fallback and refusals.
This targeted execution does not replace the unavailable full Cargo/DSR gate.

### Native angular clocks and actual real fields without allocation

The native-rad/s `prepare_uniform_band_rad_s` and
`uniform_band_mode_count_rad_s` retain source-authored angular frequencies
without a rounded conversion through hertz. The Hz and rad/s `_in` preparation
methods borrow `UniformHarmonicBandWorkspace`, whose fixed numerical payload
is 19,456 bytes. Its `UniformHarmonicBandRef` shares the cached owner's basis
builder, projection and contractions. Changing windows can therefore prepare
one factor and project several pressure, velocity, electric, magnetic or
reservoir components without heap allocation. Scratch can change on refusal;
a later successful preparation overwrites every admitted entry.

The explicit `prepare_uniform_real_band_hz` and `_rad_s` constructors, with
their corresponding preflight and `_in` forms, evaluate actual real-peak
fields. Admission requires 3 through 32 original tones and
`max(abs(omega)) * duration / 2 <= 0.5`; Hz uses its original cycle clock.
The absolute carrier remains. Even temporal modes retain the rotated real
amplitude and odd modes its imaginary amplitude. This is the conjugate-pair
frequency expansion evaluated by parity, preserving both sum- and
difference-frequency terms without duplicating sources or halving tiny
amplitudes. At least `min(2 * tones, 17)` modes are retained. The mode-zero
real coefficient is the actual window mean. One/two-tone, wider-band and
nonuniform-window queries retain the existing pair-kernel fallback.

Real factors return full real products. Analytic factors with explicitly
proper complex peak statistics use positive scale 0.5 through
`mode_cross_parts_dd_scaled`, `mean_square_parts_dd_scaled` or
`mean_cross_parts_dd_scaled`. The scale enters each modal weight before
multiplication, preserving finite half-power when the unscaled squared norm
would overflow. Weights below one scale the larger operand to preserve tiny
cross operands; weights at least one scale the smaller operand. Positive
finite scale admission uses existing error variants. Signed cross moments
remain signed; this representation adds no clipping, inferred independence,
carrier-period average or new underflow-refusal policy.

The existing bounds on retained modes, modal-tail power and floating accuracy
remain. Work methods count original source columns, including the optimized
real path. Source-pair caps in consumers remain structural admission bounds;
the additional modal work is separately bounded by 64 source columns and
17 modes, with cancellation polled around preparation and projection leaves.

Five additional regressions bring the actual harmonic-window/spectrum native
suite to 77 passing tests. They exercise native-angular proper and actual-real
fourth-order nulls, real field means, independent real-time quadrature with
sum tones, source-clock preservation at a large epoch, fixed/cached equivalence,
workspace reuse after refusal, finite proper power for amplitude 1.5e154,
and subnormal real/cross components. The full workspace Cargo/DSR gate remains
unavailable.

## Causal original-state harmonic decay moments

`PreparedHarmonicDecayInterval` advances the original scalar equation
`q_dot + q/tau = Re(sum drive_n exp(-i omega_n t))`, retaining the supplied
initial state, original positive relaxation time and signed angular source
frequencies. The normalized DD source clock and physical duration remain
separate. Endpoint state, direct endpoint increment, mean total-state square
and independently contracted mixed work are returned. Optional checked
integrals defer the duration multiplication when its result would underflow.

Up to 64 authored tones are admitted before allocation. Equal and conjugate
real carriers combine before quadratic products. Short intervals use the
original ODE polynomial with positive shifted-Legendre squared modes and an
explicit analytic remainder admission. Longer intervals use closed DD
pair/decay moments and refuse unresolved cancellation. Large retained states
are decayed without first rounding a tiny exponential factor. Preparation and
queries poll cancellation and admit explicit payload/work bounds; queries do
not allocate. Existing settled causal-product arithmetic remains unchanged.

There is no settled-memory reset, negative-energy clamp, sampled time grid or
interval/null-relative floating-point certificate. Exponent-range refusals and
the accuracy of the original deterministic transcendental functions remain.
Nine G1/G3/G4 regression tests specify independent real ODE integration,
cosine-peak startup, tiny direct increments, subdivision and large DD clocks,
coherent duplicate/conjugate carriers, stiff decay, budgets and cancellation.
They were authored and reviewed but not compiled or run: the execution
environment became unavailable. Independent equation checks are not native
Rust test evidence.

## Physical peak normalization across numerical range

Product-spectrum preparation and proper-field pair contraction preserve the physical factor of one half when an unnormalized intermediate would overflow. Finite ordinary DD component sums retain their existing contraction order, including representable subnormal sums. Only overflowing intermediates move the half into individual products, scaling the larger operand to preserve a small partner. The exact half-weight modal cross path also combines subnormal components before rounding away each half. Native clocks, real sum tones, signed cross moments, source statistics and admission limits are unchanged; truly unrepresentable physical results remain refusals. `harmonic_window/spectrum/normalization_tests.rs` covers finite proper power, real difference/sum signs, signed pair fallback, tiny products and genuine overflow.

Instantaneous native angular phases are exposed by `phase_rad_s_from_clock_dd` and `phase_rad_s_dd_from_clock_dd`. Both retain the DD source epoch and negative-time convention without a finite-shutter surrogate or rounded Hz conversion. The DD-frequency entry also preserves authored difference-frequency low parts. Signed frequencies, DC, invalid inputs, small cosine residuals and sub-ulp epochs have focused regressions; the existing deterministic radian reduction and finite-range refusals remain authoritative.

## Retained source epochs for complete observation windows

`HarmonicWindow::from_clock_dd` and `with_profile_from_clock_dd` retain a
normalized DD source epoch independently of the positive finite duration.
`start_clock_dd` returns that epoch; `start_s` remains its rounded binary64
accessor. Existing constructors retain their source units and window profiles.
All original radian/Hz phase, deterministic/proper pair and positive temporal
band paths now consume the retained epoch. Start and interval rotations remain
separate; an interval need not advance the rounded timestamp. Invalid source
parts, nonpositive duration and unrepresentable physical endpoints refuse.

The coherent optical shutter exposes `try_from_clock_dd` and `start_clock_dd`.
Original EM aperture center/first-moment and integrated heat handoffs retain
both epoch parts. Single/multitone tensor-Prony energy and causal observation
handoffs also preserve that same source epoch. No new source coherence,
material history, temporal sampling or shutter profile is inferred.

Ninety-eight actual shared harmonic tests pass with warnings denied, including
six new groups for source normalization/admission, independent fractional-cycle
Hz phases, real/proper products, native-radian field-phase translation, owned
and borrowed positive bands, and spectrum interval partition/continuation.
These exercise uniform and triangular profiles and clocks whose low component
is smaller than one f64 step. Downstream handoff types and source paths were
inspected; complete optical, EM and material integration execution remains
unrun. DD arithmetic is not an arbitrary-null relative-error certificate.

## Positive coherent factors for triangular observations

The profile-aware `band_mode_count_hz/rad_s`,
`real_band_mode_count_hz/rad_s`, `prepare_band_hz/rad_s`,
`prepare_real_band_hz/rad_s` and their borrowed `_in` constructors admit
Uniform and Triangular normalized detector profiles. `PreparedHarmonicBand`,
`HarmonicBandRef` and `HarmonicBandWorkspace` share the original storage and
contractions. Existing uniform-only constructors retain their original domain.

Triangular sensitivity is integrated once against the observable. Its positive
factor uses orthogonal polynomials for `(1-|x|) dx` on [-1,1], derived from exact
rational moments `mu_(2k)=2/((2k+1)(2k+2))`. A positive three-term recurrence
forms Taylor projection coefficients before original fields are combined.
No bright frequency-pair matrix, Cholesky repair, time raster or field taper is
introduced. The modal weight remains `2m+1`; the zeroth actual-real coefficient
is the normalized weighted field mean.

Original Hz/rad-s clocks, DD source epochs, real-carrier sum terms and proper
half-normalization remain explicit. The existing 3..64 analytic / 3..32 real
source and half-phase-width <=0.5 admission supplies at most 17 modes.
Orthogonal triangular projection has exact-arithmetic omitted power below
`1e-49*(sum |A_i|)^2`: its best weighted residual is at most twice the uniform
Legendre residual. Twenty Taylor terms per mode suppress the separate analytic
projection remainder. Floating phase/coefficient error, source uncertainty and
arbitrary-null relative accuracy are not certified.

Retained workspace stays 19,456 numerical bytes. Triangular preparation adds
6,864 bounded numerical scratch bytes for the recurrence and Taylor rows.
It allocates no heap in the borrowed path; observations retain O(F*M) field
projection and the same work counts. A normalized shaped mean does not by itself
give unweighted physical interval energy or an endpoint storage change.

All 102 current shared harmonic tests pass with warnings denied, including four
new groups: proper fourth-order and real even/odd nulls through 1 ps; independent
weighted time integration of real and complex cross moments; retained DD epochs
and proper carrier shifts; and scale/admission/uniform-compatibility behavior.
Owned and borrowed factors agree in both original frequency units.
