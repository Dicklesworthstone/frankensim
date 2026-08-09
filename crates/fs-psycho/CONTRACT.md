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
- `roughness::roughness_dw_block` — Daniel-Weber roughness [asper]
  per 8192-sample analysis block (Pa input), ported from the
  Apache-2.0 MoSQITo reference with tables transcribed mechanically
  ([`dw_tables`]); numerically EXACT against the reference re-run
  standalone at the same block length (better than 1e-12 relative at
  every AM modulation frequency, all seven values pinned in the
  test). Executed port lessons recorded in the module: the
  reference's empirical 1.42 one-sided factor (a 2.0 doubling read
  2.59 asper at the anchor), its zero negative-frequency H half
  (mirroring + doubling read ~4x), and — review-caught — its
  floor()-truncated H support with the 502 Hz cutoff REUSED for
  H16/H21/H42 (interpolating to each table's own endpoint read only
  7-9 matching digits behind a loose 1e-3 pin). Disclosed residual
  deviations (edge-channel numpy wraparound, corrcoef NaN) in the
  module doc.
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
- `LISTENING_LAW` — the not-a-substitute statement as data, so its
  removal breaks a test.

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
12. Roughness: the 100% AM 1 kHz 60 dB sweep PEAKS in 55..90 Hz (the
   published ~70 Hz signature), ALL SEVEN sweep values match the
   standalone reference run within 1e-12 relative (exactness pins;
   R(70 Hz) = 1.0448 asper at the published ~1-asper anchor), falls
   off both sides, and an unmodulated tone is far smoother. Roughness
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

`tests/psycho.rs` (16 + 1 ignored provenance tool): ISO signal-1
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
the reference binary consumed to mint every signal-path pin.

## No-claim boundaries (the bead's remaining scope — OPEN)

- Fluctuation strength: not implemented (roughness IS — see above).
- Roughness time-series over long signals (block averaging /
  overlap): single-block v1.
- Tonality/harmonicity: not implemented.
- Batch API for Pareto consumers: not implemented.
- Sampling rates other than 48 kHz: refused, not resampled (no
  claim); wav parsing is the caller's job (inputs are Pa slices).
- Per-frame specific loudness is computed internally but not
  returned by `loudness_time_varying` (memory scale); time-varying
  sharpness is a recorded follow-up.
- Mono, free/diffuse-field assumption per ISO 532-1; no binaural
  claims.
