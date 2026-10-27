# Calibrated native material comparisons

The native wind/string panel retains two rendered results, not two browser
synthesizers. Material substitution, atmosphere evaluation and physical time
integration remain in the selected native executable.

## Compare a change

Render the baseline with the existing local native service. Give it a label and
choose **Keep current result as A**. Use **Replace…** on an attached source to
swap its original bytes while retaining its declared path. The performance and
any gas/profile associations therefore keep referring to the same path; no
physics coefficients are translated into browser controls. The file hash changes
and the existing WAV becomes visibly stale until another successful render.
Atmosphere records can also be edited directly in the performance.

Render the modified performance. It becomes **B**; A keeps its WAV, source text,
attachment identities, options and calibration. An unsuccessful or cancelled
render does not replace either result. **Play comparison** starts both accepted
WAVs at the same simulation time and sample offset. **Listen A/B** switches on
that running clock with a complementary 5 ms crossfade; it does not restart the
note. Seeking pauses the pair; resume explicitly with Play. Listening volume is
one common control, initially 25%, and never a per-result normalization.

Each capture exposes its submitted source and calibration. Both complete native
WAVs can be downloaded unchanged. The browser holds only one reference and one
current result; replacing/clearing them releases their old object URLs. Leaving
the page releases playback and captures. They are not a persistent library.

## What is measured

The service reads `music_render`'s native provenance sidecar or
`material_strings`' final writer summary. The returned pressure receipt carries
the actual output sample rate/count, full-scale pressure and clipped-sample
count, bound to the WAV's SHA-256. The client verifies the receipt against the
returned bytes. Neither layer infers calibration from a material name, a
performance string, a peak, or a nominal input clock.

The comparison displays RMS and peak pressure for A and B, B/A RMS level in dB,
and RMS sample difference. These are **quantized PCM pressure measurements**,
not unquantized solver diagnostics, measured instrument accuracy, or sound
pressure calibrated at the listener's physical headphones. Pressure reconstruction
uses the native encoder's `i16::MAX` scale (32767). Web Audio playback uses its
normal signed PCM conversion (32768) equally on both sides.

A and B must have exactly the same native pressure full scale and no reported
encoder clips. Different scales or clipping disable calibrated comparison;
changing the playback slider cannot repair lost pressure. Missing legacy
calibration still permits ordinary WAV playback but not calibrated A/B.
Malformed emitted calibration is a visible refusal, never a silent default.

There is no onset alignment, retuning, time stretching, loudness matching or
source-gain compensation. Native propagation delays remain in the samples.
For unequal durations, playback and statistics explicitly use only the common
prefix from simulation time zero; complete downloads are unaffected. A silent
reference produces an undefined level ratio rather than a made-up finite dB
value. The audio device may resample playback; measurements always use the
original 48 kHz mono PCM16 samples.

Each WAV retains the native service's 16 MiB limit. Two captured WAVs and their
bounded decoded playback copies are retained in browser memory. The service's
existing source, attachment, diagnostics, subprocess and cancellation limits
remain unchanged.

## Focused checks

```sh
node --test apps/instrument-lab/tests/native-pressure.test.mjs \
  apps/instrument-lab/tests/native-comparison-playback.test.mjs
```

The HTTP tests use explicitly synthetic executable output, not native physics.
The playback tests exercise the actual player with an observable AudioContext
boundary: common starts, offsets and windows; complementary switching; retained
amplitude; cancellation during resume; and refusal before pair replacement.
The panel also received an offline Chromium DOM/Web Audio smoke check with
synthetic source/response fixtures. Neither this workflow nor those tests
qualifies any material model or instrument recording as physically accurate.
