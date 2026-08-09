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
  sones + specific loudness on the 0.1-Bark grid. Sones ONLY: a phon
  conversion is deliberately absent (the reference outputs sones; a
  from-memory formula was written and REMOVED — no unverified
  claims).
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
  standalone at the same block length (10+ digits at every AM
  modulation frequency). Executed port lessons recorded in the
  module: the reference's empirical 1.42 one-sided factor (a 2.0
  doubling read 2.59 asper at the anchor) and its zero
  negative-frequency H half (mirroring + doubling read ~4x).
- `LISTENING_LAW` — the not-a-substitute statement as data, so its
  removal breaks a test.

## Invariants

1. EXACT-path fidelity: ISO Annex B.2 test signal 1 reproduces the
   reference implementation's 83.29566 sone within 0.001 sone (the
   compiled reference binary was run to establish it; hiding a port
   bug inside the standard's 5% compliance tolerance is exactly what
   caught the executed `-.25f` transcription mangling — see the
   tables module doc).
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
9. Roughness: the 100% AM 1 kHz 60 dB sweep PEAKS in 55..90 Hz (the
   published ~70 Hz signature), R(70 Hz) matches the standalone
   reference's 1.0448 within 1e-3 (exactness pin), falls off both
   sides, and an unmodulated tone is far smoother.

## Error model

Typed `PsychoError` (shape, non-finite, uncalibrated-absolute,
degenerate signal). No silent degradation, no fabricated values.

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

`tests/psycho.rs` (11): ISO signal-1 exactness; cross-path tone
references; anchor + monotonicity; sharpness behavior + silence
refusal; low-band mutation; calibration refusal + calibrated value;
log-attack-time ordering + linear-ramp closed form; bitwise
determinism + typed refusals; diffuse-vs-free; listening-law pin;
AM-roughness sweep with the exactness pin.

## No-claim boundaries (the bead's remaining scope — OPEN)

- PCM third-octave filterbank (the wav path) and TIME-VARYING
  loudness: not implemented; loudness input is band levels.
- Fluctuation strength: not implemented (roughness IS — see above).
- Roughness time-series over long signals (block averaging /
  overlap): single-block v1; the wrapper is trivial once the
  time-varying loudness lands.
- Tonality/harmonicity: not implemented.
- Phon conversion: deliberately absent pending a verified source.
- Batch API for Pareto consumers: not implemented.
- Mono, free/diffuse-field assumption per ISO 532-1; no binaural
  claims.
