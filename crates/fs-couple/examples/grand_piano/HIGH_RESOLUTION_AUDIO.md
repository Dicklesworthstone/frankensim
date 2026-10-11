# Preserve the physical performance in high-resolution audio

`piano_audio` encodes the existing unquantized receiver-pressure CSV as **24-bit
PCM or 32-bit floating-point WAV**, without rerunning the expensive mechanics.
It does not up-convert an already quantized WAV. Quiet pressure that rounded to
zero in the original PCM16 preview can now remain present in the exported file.

```sh
cargo run --release -p fs-couple --example grand_piano -- \
  --preset steinway-d --sample-rate 48000 --pcm-full-scale-pa 64 \
  --receiver-pressure-csv performance.csv --render preview.wav
cargo run --release -p fs-couple --example piano_audio -- \
  --input performance.csv --sample-rate 48000 --full-scale-pa 64 \
  --format pcm24 --output performance-24.wav
```

The first command still writes the original PCM16 preview. Its existing
headroom check still applies: choose a larger explicit `--pcm-full-scale-pa`
if the modeled peak reaches the declared scale. This only changes preview
encoding gain; receiver CSV samples remain the actual modeled pascals. The
number 64 above is an encoding choice, not acoustic calibration or an automatic
peak fit. Existing MIDI, performance, microphone and physical-model options
remain available on the first command.

For floating-point output, change `--format pcm24` to `--format float32` and
choose a fresh output path. The exporter uses one explicit Pa-per-full-scale
value for both channels. It preserves frame order and sample clock, checking
every CSV timestamp against `--sample-rate`. Missing frames, changed clocks,
invalid/nonfinite samples, unknown formats and existing destinations refuse.

## Precision and headroom

PCM24 uses packed signed little-endian samples and symmetric scaling to
+/-8388607. It has approximately 256 times finer quantization than the existing
PCM16 path at the same pressure scale. It refuses over-range input rather than
silently clipping. Float32 retains finite representable normalized values even
outside +/-1, as well as much smaller tails. **A player may still clip float
values beyond +/-1**; explicitly attenuate before fixed-point playback. Float
export rounds f64 to f32 and is not a bit-exact archive of the original trace.

Both files use WAVEFORMATEXTENSIBLE with explicit sample precision and mono or
left/right channel masks. Floating output includes a frame-count `fact` chunk;
odd PCM payloads include RIFF padding. There is no peak normalization, per-channel
level adjustment, dither, resampling, EQ, synthetic tail or harmonic generation.
Retain the supplied pressure scale with the file when comparing renders.

The batch converter reads bounded rows and retains at most 33554432 channel
values (about 349 seconds of stereo at 48 kHz). The shared library API
`fs_couple::pcm_wav::high_resolution::encode_pressure_wav` also accepts an
existing in-memory pressure history directly; the old PCM16 APIs are unchanged.

## Focused checks

```sh
cargo test -p fs-couple --lib pcm_wav::high_resolution
cargo test -p fs-couple --example piano_audio
```

Six encoder regressions cover exact headers, signed 24-bit packing, padding,
floating headroom, quiet values, stereo order and refusals. Four converter
regressions cover the real receiver CSV layout, clock/row admission, explicit
options and fresh-output behavior. Independent Python reference files in all
four format/channel combinations were decoded successfully by SciPy; that is
format-reference evidence, not execution of the new Rust encoder. Native tests
and paired listening are still needed. Export precision does not correct
missing or inaccurate mechanical/acoustic physics.
