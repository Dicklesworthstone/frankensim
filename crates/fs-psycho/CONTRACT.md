# fs-psycho — Contract

## Purpose and layer

Layer L3. Standards-anchored psychoacoustic metrics (bead
frankensim-fsim-psychoacoustic-4y0v8, IN PROGRESS — this crate is the
first slice): perceptual objective axes for Pareto fits and
listening-evidence tables. THESE METRICS ARE NEVER A SUBSTITUTE FOR
HUMAN LISTENING (the program's listening law, pinned as data and by
test).

## Public types and semantics

- `loudness_stationary(levels, field)` — stationary Zwicker loudness
  (ISO 532-1 clause-5 third-octave-level method), ported
  statement-for-statement from the standard's own free reference
  implementation; data tables in [`tables`] transcribed MECHANICALLY
  by script from that C source (standards.iso.org free electronic
  insert, Annex A.4). Input: 28 third-octave levels 25 Hz..12.5 kHz
  in dB SPL (already absolute); free or diffuse field. Output: total
  sones + specific loudness on the 0.1-Bark grid. Sones here; the
  phon conversion lives in `signal::phon_from_sone` (a from-memory
  formula was once written and REMOVED; the conversion landed only
  after the reference's own `f_sone_to_phon` was found).
- `sharpness_din(specific)` — DIN 45692 sharpness [acum] over the
  specific-loudness pattern: `S = 0.11 sum(N' g(z) z dz)/N`, `g = 1`
  below 15.8 Bark, `0.15 exp(0.42 (z-15.8)) + 0.85` above (the
  standard's weighting, cross-checked against the Apache-2.0 MoSQITo
  reference). Level-RELATIVE (no calibration). Silence REFUSES
  (sharpness of nothing is undefined, not zero).
- `log_attack_time(pcm, rate, window)` — timbre-toolbox
  `log10(t90 - t10)` on a moving-average envelope. Level-relative.
- `spl_from_pcm_rms(pcm, calibration)` — the ONLY PCM-to-absolute
  bridge in v1: REFUSES without a [`Calibration`]
  (`UncalibratedAbsolute`, the honest-scope law in a type).
- `Calibration::from_pressure_scale(pa_per_full_scale, reference_pressure_pa)`
  converts an explicit linear pressure-per-sample scale to the existing
  **peak-one sine** calibration, including its RMS factor exactly once.
  `spl_from_pressure_rms` observes a physical acoustic pressure record directly,
  with an explicit RMS reference (20 µPa is supplied as
  `AIR_REFERENCE_PRESSURE_PA`). Both level paths use the same scaled
  sum-of-squares kernel and logarithmic reconstruction: finite weak signals
  are not silently called zero, and intermediate squaring cannot overflow.
  No source independence, peak normalization, DC removal, weighting filter,
  microphone calibration authority, or absolute thermodynamic pressure meaning
  is inferred. Coherent sources must be combined before observation.
- `roughness::roughness_dw_block` — Daniel-Weber roughness [asper]
  per 8192-sample analysis block (Pa input), ported from the
  Apache-2.0 MoSQITo reference with tables transcribed mechanically
  ([`dw_tables`]); checked against the reference re-run standalone at
  the same block length with its FFT frequency axis corrected to k*fs/N,
  threshold interpolation aligned to that axis, and DC excluded from audible
  tones. The seven AM-modulation pins and eight sideband cases at 44.1/48/96
  kHz agree within 1e-12 relative. Correct physical frequencies feed ear
  weighting, hearing thresholds, Terhardt slopes and Bark-channel placement;
  the former k+1 labels changed a 1 kHz / 40 Hz AM result by 6.6%.
  Executed port lessons recorded in the module: the
  reference's empirical 1.42 one-sided factor (a 2.0 doubling read
  2.59 asper at the anchor), its zero negative-frequency H half
  (mirroring + doubling read ~4x), and — review-caught — its
  floor()-truncated H support with the 502 Hz cutoff REUSED for
  H16/H21/H42 (interpolating to each table's own endpoint read only
  7-9 matching digits behind a loose 1e-3 pin). Disclosed residual
  deviations (edge-channel numpy wraparound, corrcoef NaN) in the
  module doc. Modulation weights remain zero at and above Nyquist even
  when a low acquisition rate intersects the tabulated band, so the
  factor-two reconstruction never double-counts negative-frequency bins.
  No resampling, DC detrending, alias-removal or perceptual acceptance claim
  is added.
- `signal::loudness_stationary_from_pcm` — the reference's
  stationary-from-signal path: 28-band 48 kHz filterbank (biquad
  difference tables in [`filter_tables`], transcribed mechanically),
  mean-square levels after `time_skip`, then the clause-5 chain.
  REPRODUCES the Annex B.3 published tone values (4.019 / 14.655
  sone) that the level-vector path missed by 13% — the filterbank
  leakage is real sones.
- `signal::loudness_time_varying` — the clause-6 time-varying method:
  filterbank, squaring + three frequency-dependent smoothing
  low-passes, 2 kHz decimation, per-frame core loudness, the
  nonlinear temporal decay (24x virtually upsampled two-capacitor
  element), slopes, 0.47/0.53 dual-low-pass temporal weighting;
  outputs the loudness series plus Nmax and the standard's N5
  percentile (the reference's own estimator, ported exactly).
- `signal::phon_from_sone` — the reference's `f_sone_to_phon`
  (verified source; previously deliberately absent): 10 log2(N) + 40
  above 1 sone, `40 (N + 0.0005)^0.35` floored at 3 phon below, with
  the reference's small branch step at 1 sone kept and documented.
- `signal::pareto_metrics` — the batch set for Pareto/listening
  consumers (stationary sones + phon, N5/Nmax, sharpness, mean
  DW-block roughness with block count, log-attack-time, ECMA T-TNR
  and T-PR on the complete native-length recording, with the original
  sample count disclosed) from one 48 kHz Pa block. AGGREGATION-EXACT by contract: every field is
  bitwise the standalone call's value (pinned by test); any
  component's typed refusal propagates — no partial results.
  TNR and PR share one windowed FFT and one candidate-screening pass;
  each method retains its own critical-band grouping and prominence rule.
- `tonality::{tone_to_noise_ecma, prominence_ratio_ecma}` — ECMA-74
  annex D tonality (TNR and PR with prominence classification, totals
  over prominent tones per ECMA TR/108), ported from the Apache-2.0
  MoSQITo reference with a corrected physical FFT frequency axis
  (Bray smoothed-spectrum screening, 1/24-octave criteria, annex
  D.8/D.10 critical bands, D.7.1 hearing threshold). Detection band
  (89.1, 11200) Hz; exact native-length block, including odd and
  non-power-of-two counts through `fs-fft::GeneralFft`. Bluestein's internal
  convolution padding does not change the record, sample clock or frequency
  spacing. Power-of-two records retain the original Stockham numerical route.
  Coefficient k is located at k*fs/N Hz: the reference's former one-bin
  offset is corrected before screening, critical-band placement and reported
  tone frequencies. These are spectral-bin centers, not sub-bin pitch
  estimates. Even records exclude Nyquist; odd records include their last
  positive-frequency coefficient. Independent reference pins apply the same
  axis correction to the upstream Python routines. Focused G1/G3 tests
  exercise 48000/48001-sample records, tones arriving after the former
  truncated prefix, physical frequencies across 44.1/48/96 kHz clocks,
  DC and the upper positive-frequency bins, and bitwise replay.
  Disclosed deviations in the module doc, all
  on paths the reference cannot execute deterministically or
  correctly: uninitialized smoothed-spectrum head/tail
  (deterministic completion), zero-magnitude floor at the
  reference's own -140 dB amp2db substitution, out-of-band LTH
  clamp, EMPTY 1/24-octave band REFUSES (the reference's bookkeeping
  corrupts — review-executed silent divergence at 96 kHz/8192 before
  this gate), and the screening left walk clamps at the band edge
  where numpy index-wraps (review-executed usize-underflow panic
  before the fix).
- `harmonic::HarmonicAnalysisPlan` — physical, known-frequency harmonic
  measurements over the existing `fs_fft::harmonic::HarmonicFitPlan`. Joint
  DC/sine/cosine fitting retains noninteger-period records, physical Hz, peak
  phasors referenced to the first sample, and the original sampled residual.
  The prepared front adds no heap allocation: callers reuse N scratch scalars
  and K phasors, and the returned view keeps those coefficients paired with
  their plan. It inherits the numerical owner's clock, workspace, conditioning,
  and transactional output admission; fixed-block analysis is synchronous.
  Fundamental-relative levels, THD, even/odd contributions and adjacent-odd
  even-harmonic deficits use period-averaged fitted harmonic amplitudes.
  Missing orders/categories are absent, exact zero levels have explicit
  infinite/undefined states, and an explicit minimum fundamental peak in input
  units withholds normalization when the reference is unresolved. Scaled norms
  preserve weak/large signals without squaring physical amplitudes directly;
  nonzero linear ratios outside finite f64 range refuse. Sparse-order THD is
  labelled incomplete through its highest retained order. DC and residual are
  never counted as harmonic distortion, and no finite-record energy partition,
  noise/alias certificate, frequency estimate, SPL calibration or listening
  verdict is inferred. `tests/harmonic.rs` supplies G1/G3 analytic amplitude,
  phase, gain, even-deficit and unmodeled-partial checks. This supplies only the
  measurement primitive for X-Pp1 (`frankensim-irp-v5-instrument-realism-pexfp.342`),
  not its resolved-jet, specimen, Fletcher-comparison or listening obligations.
- `partials::track_partials` — full-record physical partial trajectories over
  the shared independent-tone, damped-pole, linear-chirp and damped-chirp
  refiners. Fixed caller-authored search domains preserve partial identities.
  Stationary models estimate independent frequencies; ringdown models also
  estimate amplitude decay and can separate equal-frequency, different-decay
  components. Chirped models identify each first- and last-sample frequency
  over the physical `(frame_samples - 1) * timestep_s` span, retaining linear
  within-frame drift. Damped-chirp models jointly identify both endpoint
  frequencies and exponential amplitude decay; disjoint decay intervals may
  separate components with matching frequency trajectories. Each frame retains
  the original DC, residual and first-sample peak phasors; its exact integer
  source sample index distinguishes
  changing phase origins without rounding a large epoch to floating seconds.
  End frequency is absent for nonchirped fits and decay is absent outside
  ringdown and damped-chirp fits, rather than silently measured as zero.
  Successful chirp endpoints predict the next frame using the actual relative sample gap,
  including shortened terminal hops and unresolved frames. Only those initial
  guesses are clamped to the immutable domains; fitted endpoints must meet the
  original untruncated-correction gate. Damped-chirp predictions retain the
  last fitted decay; each frame independently refits its initial amplitude
  and phase at the new source origin, retaining exponential envelope loss
  instead of reporting a window-averaged amplitude. A terminal overlapping full
  frame covers the last sample without padding, attack skips or tail truncation.
  Unidentifiable, ill-conditioned and unconverged frames retain typed refusals
  and leave the last successful seeds and their source offset unchanged.
  Results stream one frame at a time; frame count and candidate work are
  admitted explicitly, seed storage (both anchor and prediction for chirps)
  is charged to the fit's numeric budget, and callbacks may stop between
  frames. No mode discovery, arbitrary track reassignment, noise-floor rule,
  uncertainty, calibrated level or instrument verdict is inferred. Nonlinear
  drift, attacks and growth require a local-model/frame-length check; the
  selected local law does not identify dynamics outside its model.
  `tests/partials.rs` covers inharmonic motion, coherent double decay,
  unresolved recovery, exact clock
  indices, bounded work/storage, deterministic replay, a three-rung stationary
  pitch-glide residual ladder, linear-chirp endpoints and phase, exact
  warm-start gaps, bounded prediction, dropout-anchor retention, simultaneous
  drift and decay, local initial amplitudes, and reidentification after a
  change in the decay law. This supplies native trajectory observations for
  instrument-plan §7.1/7.2; it
  does not replace the Python grader or establish the struck/plucked Q2 decay
  gate.
- `fluctuation::fluctuation_strength_frame` — fluctuation strength
  [vacil] per 131072-sample frame, implemented from the PUBLISHED
  equations of Osses/Garcia/Kohlrausch (ICA2016-113; NO code ported —
  the only open implementations are CC BY-NC or unlicensed):
  Eq. 1 over 47 critical-band channels with Terhardt excitation
  patterns, generalised modulation depth with the paper's 3:1
  compression above 0.7, `H(fmod)` band-pass, cross covariances
  i±2, `g(z)` high-frequency taper, exponents 1.7. Disclosed
  realisation choices (module doc): power-of-two frame (2.73 s vs
  the paper's 2 s), zero-phase Butterworth-style `H` magnitudes with
  corners FITTED to the published AM curve per the paper's own
  methodology (1st-order HP 1.7 Hz, 3rd-order LP 13 Hz), `C_FS`
  re-derived on the paper's 1-vacil reference, DW threshold table
  reuse, edge-channel single-factor rule. Component thresholds and
  excitation slopes use RMS SPL (`sqrt(2) |X[k]|/N` for interior
  real-signal bins), not peak pressure. Ear transmission and spreading
  scale the original complex FFT coefficients directly, preserving pressure
  and phase without repeating magnitude calculations for all 47 channels.
  The existing reference fixture pins carrier RMS SPL; full modulation
  raises its long-term overall SPL by about 1.76 dB. That disclosed fixture
  convention is retained, rather than claiming the paper supplies its exact
  carrier-to-total normalization recipe.
- `LISTENING_LAW` — the not-a-substitute statement as data, so its
  removal breaks a test.
- `receipt::ListeningReceipt` (music bead 3ez8g.1.2) — the human half of
  the listening law as a durable artifact: who listened, to exactly which
  content-addressed rendered artifact, judging exactly which question,
  with what verdict, with metrics attached as EVIDENCE FIELDS (context,
  never the verdict). `ListeningVerdict::Unadjudicated` is first-class —
  a rendered fixture awaiting an ear is a real, recordable state — but
  structurally useless as pass evidence: only an adjudicated Pass answers
  `supports_pass()`, so a claims-registry row can never cite an
  unlistened receipt as approval. Calibration honesty is structural:
  `AttachedMetrics` fields are `Option`s that stay `None` when
  uncomputed — without a `Calibration`, only level-relative
  log-attack-time is attachable; loudness/sharpness/SPL stay `absent` in
  the canonical bytes, never fabricated. Canonical line encoding
  (`frankensim-listening-receipt-v1`, strict round-trip decoder) carries
  the LISTENING_LAW verbatim as its terminal line — a receipt that
  rewrites the law refuses to decode. A FAIL verdict without
  observations refuses (a human who failed a fixture owes the next agent
  the reason); an adjudicated verdict behind a placeholder listener name
  refuses. No wall-clock in the bytes: `session` is a caller-declared
  label. Digesting the bytes happens where fs-blake3 is in scope (this
  crate's dependency set stays fs-fft + fs-math only).

## Invariants

1. EXACT-path fidelity: ISO Annex B.2 test signal 1 reproduces the
   compiled reference binary's 83.29566042436214 sone within 1e-9
   RELATIVE. Two executed lessons live in this gate: the `-.25f`
   transcription mangling (caught at 0.001 sone), and the
   float32-literal table rounding (C's `-0.6f` widened into a double
   is not 0.6) which a 0.001-sone gate NEVER saw — surfaced only by
   the signal-path exactness work and fixed by round-tripping every
   f-suffixed constant through f32 in the generator (tables module
   doc).
2. Cross-path tones (Annex B.3 wav references 14.655/4.019/1.549
   sone): within 15% + 0.1 sone — the single-band level vector
   under-represents filterbank leakage that loudness compression
   turns into real sones (measured 13% at 1 kHz; the exact-path pin
   is invariant 1).
3. The sone scale's anchors behave: 1 kHz 40 dB within 0.15 of
   1 sone; monotone in level with ~doubling per 10 dB (1.7..2.4).
4. Sharpness rises with frequency (4 kHz > 2x 250 Hz at equal
   level); the 1 kHz 60 dB pattern lands within 0.7..1.3 acum of the
   DIN anchor.
5. Mutation: silencing the 11 low bands (the DLL/LCB machinery)
   moves test signal 1 outside 3x the standard tolerance — the
   low-band tables are load-bearing.
6. Calibration: uncalibrated absolute SPL REFUSES by name;
   calibrated half-scale sine reads -6.02 dB re full scale.
   G1/G3 amplitude tests also pin 20 µPa RMS = 0 dB, 1 Pa RMS =
   93.9794 dB, the +3.0103 dB square-versus-sine crest-factor difference,
   coherent doubling versus quadrature energy addition, explicit-reference
   scaling, finite subnormal/large levels, and typed non-finite refusals.
   Log-attack-time normalizes its envelope before accumulation, preserving
   its level-relative meaning across changes of signal units and gain.
7. Diffuse field differs measurably from free field (DDF live).
8. Determinism bitwise; typed refusals for shape/NaN/degenerate.
9. Signal-path fidelity: on bit-identical fixture PCM (1 kHz 60 dB
   tone, 250 Hz 80 dB tone, 500 ms 1 kHz 70 dB pulse) the stationary
   value, the time-varying Nmax/N5, and rise/plateau probe frames
   match the compiled reference binary within 1e-9 relative; the
   pulse's decay tail carries a one-time branch-flip offset (the
   nonlinear element's 1e-5 equality band under ulp libm differences;
   full-series review measurement: peak -7.3e-9 near the flip, then a
   constant -2.9e-9, zero frames beyond 1e-8) and gates at 1e-7. The
   stationary signal path lands on the Annex B.3 published values
   within 0.05 sone. Temporal behavior: silence before onset,
   loudness persisting >15% of Nmax 90 ms after offset.
10. Phon conversion: 1/2/4 sone map to 40/50/60 phon exactly; the
    sub-sone branch and its 3-phon floor match the reference formula;
    negative and non-finite inputs refuse by name.
11. Signal-path refusals: non-48 kHz rates refuse by name
    (`UnsupportedRate` — the reference ships 48 kHz tables only,
    resampling would be a silent claim), non-finite samples, too-short
    signals, out-of-range or signal-consuming `time_skip`, and
    finite samples that OVERFLOW the squaring stage (1e200 Pa once
    read Ok(inf/NaN) — review-caught, now a typed refusal) all
    refuse typed.
12. Tonality fidelity: on bit-identical fixture PCM the extracted
    reference with its FFT axis corrected reproduces — two-tone fixture:
    TNR 12.44/51.33 dB and
    PR 61.49/50.70 dB at the exact bin frequencies, all prominent,
    ratios pinned at 1e-9 relative; noise-only fixture: TNR finds NO
    tones and PR finds EXACTLY the reference's 10 weak candidates,
    every frequency bit-pinned (none prominent, totals exactly 0 —
    the candidate-set equality is the sharp gate: it caught the
    even-b band-edge bug the strong tones sailed through); low-SNR
    fixture: the 500 Hz tone detected, prominent, pinned; close-pair
    fixture exercises the TNR proximity-MERGE branch (990+1010 Hz in
    one band -> one merged tone, pinned); far-pair fixture exercises
    the far-apart branch's ACCEPT path (bin-exact tones; a weak
    non-prominent tone with total exactly 0). A boundary pair at
    500.244140625/535.400390625 Hz is separated by 35.15625 Hz,
    outside the reference's 33.5293 Hz proximity limit. The factor 1.2
    multiplies the logarithm's 1.8 power; putting it inside that power
    incorrectly widens the limit to 36.0856 Hz and merges this pair.
    Tonality refusals
    (short, NaN, bad rate, unresolvable detection
    band, too-coarse resolution for 1/24-octave smoothing incl. the
    review's panic configuration) typed and executed.
13. Fluctuation strength: the 1-vacil reference stimulus reads
    EXACTLY 1.0 (calibration regression, bitwise-stable); the
    published AM-tone curve (Fastl-Zwicker values via the paper's
    Table 1) is reproduced within +-(0.05 + 25%) per point —
    about 12-13% low at 4-8 Hz, a LEVEL-GROWTH deficit (60->70 dB
    growth weaker than published) disclosed with direction asserted;
    shape gates (bandpass peak at 4-8 Hz, monotone rise, 32 Hz
    collapse) discriminate flat models; FM values are REPORTED with
    the overestimation disclosed (the paper's own model
    overestimates FM > 4 Hz; this realisation reads ~1.6x at the
    4 Hz FM point) and gated qualitatively; the Bark formula
    cross-checks the independent Daniel-Weber BARK_FREQS table
    within 0.35 Bark; AM BROADBAND NOISE (the decorrelated-channel
    class that exercises the Eq 9 cross covariance — review-executed:
    a k := 1 mutant passed every tone gate) reads 1.845 vacil at the
    published 1.80 point (envelope 0.9..2.7 + bandpass shape); the
    paper's 3:1 compression worked example (0.85 -> 0.75) is pinned
    directly (the tone battery never crosses the knee —
    review-executed clamp mutant); unmodulated/silent inputs are
    smooth/zero; refusals typed; determinism bitwise.
14. Roughness: the 100% AM 1 kHz 60 dB sweep PEAKS in 55..90 Hz (the
   published ~70 Hz signature), ALL SEVEN sweep values match the
   standalone reference with the physical frequency-axis correction within
   1e-12 relative (R(70 Hz) = 1.0693 asper at the published ~1-asper anchor),
   falls off both sides, and an unmodulated tone is far smoother. Independent
   sideband pins cover low carriers and 44.1/48/96 kHz clocks; modulation
   support tests keep DC and the negative-frequency half zero. Roughness
   refusals (short block, NaN sample, NaN/zero/negative/INFINITE
   sample rate) are typed and executed. (This item was numbered 9
   before the signal path landed.)

## Error model

Typed `PsychoError` (shape, non-finite, uncalibrated-absolute,
degenerate signal, unsupported-rate). No silent degradation, no
fabricated values.

## Determinism class

Deterministic: table-driven, fs_math::det transcendentals, no RNG.

## Cancellation behavior

Synchronous, milliseconds-class. No `Cx` integration (workspace
`frankensim-ccmn`).

## Unsafe boundary

`#![deny(unsafe_code)]` via workspace lints; no capsules.

## Feature flags

None.

## Conformance tests

`tests/fluctuation.rs` (9 + 2 ignored probes): reference
calibration; published AM curve + shape; FM behavior; AM broadband
noise (cross-covariance discriminator); compression-knee pin;
smoothness; Bark cross-check; refusals + determinism; subthreshold versus
audible two-tone beats at known component RMS levels. Two private stage tests
also check known-RMS FFT coefficients across phase and physical-pressure
reconstruction through ear/spreading gains and an inverse FFT.

`tests/psycho.rs` (plus an ignored provenance tool): ISO signal-1
exactness; cross-path tone
references; anchor + monotonicity; sharpness behavior + silence
refusal; low-band mutation; calibration refusal + calibrated value;
log-attack-time ordering + linear-ramp closed form; bitwise
determinism + typed refusals; diffuse-vs-free; listening-law pin;
AM-roughness sweep with seven-point exactness pins; roughness typed
refusals; signal-path stationary exactness (+ Annex B.3 published
values); time-varying steady-tone and pulse exactness pins with
temporal-asymmetry behavior; phon conversion + signal refusals;
`dump_reference_signals` (ignored) regenerates the bit-identical PCM
the reference binary consumed to mint every signal-path pin; batch
aggregation-exactness (bitwise) + refusal propagation; tonality
exactness on six fixtures (including the proximity-boundary pair and both
multi-tone branches) + typed tonality refusals; native-length tail-tone
coverage, physical-frequency/sample-rate covariance and bitwise replay.
The private spectrum test also checks DC and the final odd/even positive bin.

Tonality oracle provenance: Apache-2.0 MoSQITo 1.2.1,
commit `d990c33f94f1b2050b2d811dea492219ebe97b30`.
The unchanged `comp_spectrum`, `_tnr_main_calc`, `_pr_main_calc` and their
upstream helpers run on the raw little-endian f64 fixture PCM generated by
`dump_reference_signals`. The uncorrected run first reproduces the former
pins; the corrected run replaces only the frequency axis with
`numpy.arange(len(spectrum)) * fs / len(pcm)` before the TNR/PR routines run.
Updated pins retain the existing 1e-9 relative ratio tolerances (and the
existing 1e-6 weak-candidate/total tolerances). Python is a development oracle,
not a production dependency.

`examples/psycho_descriptors.rs` (4 unit tests + IRP M7.5-W1 Python
bridge `test_psycho_bridge.py`): bitwise equality between parsed
`fs-psycho-descriptors/1` JSON output and direct Rust calls; refusal
without `--pa-per-full-scale`; partial per-metric refusal on short
blocks while preserving valid metrics; RIFF WAV parser rejection of
multichannel and oversize inputs.
The report also emits unweighted full-record `spl_db` with an explicit
`spl_reference_pa` of 20 µPa, retaining the existing calibration-source and
`level_unverified` fields. This observation includes attack and release and is
independent of the stationary loudness filterbank's skip interval.

## No-claim boundaries (the bead's remaining scope — OPEN)

- Fluctuation strength: single-frame API v1 (no 90%-overlap frame
  train); level growth with SPL is weaker than published (measured
  about 12-13% at 70 dB mid-band — the AM envelope gates carry it,
  disclosed); FM absolute values overestimate (reported, gated
  qualitatively); NOT in `pareto_metrics` (its 2.73-s frame exceeds
  the batch block contract — recorded decision).
- Roughness time-series over long signals (block averaging /
  overlap): single-block v1.
- Tonality: TNR/PR are the ECMA-74 discrete-tone metrics; Aures/
  Sottek psychoacoustic tonality (sharpness-of-hearing model) is a
  different metric and is NOT implemented.
- The ISO 532-1 PCM filterbank and Pareto batch require 48 kHz, without
  resampling. Standalone ECMA tonality honors the declared sample rate,
  subject to its detection-band and resolution admission. WAV parsing is
  the caller's job in library APIs (inputs are
  Pa slices; the `psycho_descriptors` CLI example parses mono 48 kHz
  PCM16/float32 WAVs for grader integration).
- Per-frame specific loudness is computed internally but not
  returned by `loudness_time_varying` (memory scale); time-varying
  sharpness is a recorded follow-up.
- Mono, free/diffuse-field assumption per ISO 532-1; no binaural
  claims.

## Synchronized partial trajectories

`partials::track_partial_channels` fits a shared physical frequency, decay and/or
linear-chirp trajectory to synchronized, equal-length, common-unit channels.
Each channel retains its original DC, complex peak phasors and sampled residual;
no channel is averaged, individually normalized, reordered or replaced by zero.
A partial invisible at one sensor can remain identifiable through the array.
The objective is the sum of original-channel squared residuals; different units,
sensor-noise weights or delay alignment require explicit caller conversion.

Both scalar and array entry points use the same schedule, immutable identity
domains, successful-anchor warm starts and actual-gap chirp predictor. An
unresolved shared fit returns one typed failure for the whole array. Channel
shapes and samples are admitted before callbacks. The original fit workspace
budget includes warm seeds and reusable channel-slice views; a single channel
uses a stack view and retains the scalar API's original storage budget.
Callbacks borrow one frame at a time and can stop subsequent work between fits.

`tests/partial_channels.rs` exercises sensor nodes and opposite-sign responses,
rotation of common-unit vector observations with nonzero residual, coherent
double decay, chirp and damped-chirp dropouts, shortened terminal hops,
deterministic replay and channel-shape/storage/callback refusals. This is a
shared-law observation, not mode discovery, automatic channel calibration or
instrument acceptance.
