# Preserve upper partials through acoustic propagation

The geometric piano's mono and stereo Rayleigh microphones now select the
shared `BaffledPressure::new_multirate_bandlimited` receiver. Both channels still
observe one mechanical trace, with the same loaded board projection, physical
areas, distances, signed spatial cancellation and pressure units. Diagnostic
volume output and finite-body exterior BEM are not changed.

The original `BaffledPressure::new` and `new_multirate` constructors remain
available as the two-tap linear-delay image. Their numerical path is unchanged.

## What changes

The former propagation kernel interpolated fractional travel times using two
adjacent **audio-rate** samples, even when mechanics ran at 4x. At a half-sample
fraction, its magnitude is `abs(cos(pi*f/fs))`: approximately -2.01 dB at 10 kHz
and -5.11 dB at 15 kHz for 48 kHz output. This is a numerical loss, separate
from physical string/wood damping, spatial interference and anti-aliasing.

The new image builds a 65-tap Blackman-windowed sinc for each patch's fractional
travel time. It normalizes DC gain, retains the original SI Rayleigh coefficient,
and combines the patch contributions into the existing modal FIR bins. No new
mechanical integrator, invented overtone, post-render EQ or pressure normalization
is involved. Every sample uses the original allocation-free observer loop.
Kernel construction is more expensive and its merged support is wider;
`multiply_adds_per_sample` includes the actual merged kernel width.

The independent Python/libm reference sweep covers 129 fractional positions and
2001 frequencies from DC through 0.45 of the output rate. It found maximum
complex error 0.000345 and magnitude error below 0.003 dB. Those are sampled
**interpolation-only** results, not a continuum certificate, compiled-Rust
measurement or an error bound for the full spatial pressure sum. Cancellation
can amplify a small relative component error. The finite kernel is not exact
at Nyquist and is not a replacement for the existing anti-alias decimator.

## Timing and compatibility

The sinc image adds **32 output samples** of fixed interpolation latency
(0.667 ms at 48 kHz). At 4x oversampling the existing anti-alias filter contributes
44 samples, so `filter_delay_samples()` now reports **76** for the piano,
excluding distance-dependent travel. `interpolation_delay_samples()` reports
32 separately. Do not add that 32 twice. The help text's 44-sample anti-alias
latency still refers to that filter alone.

`delay_samples` describes the complete propagation FIR support, not a group
delay. No tap occurs before the floor of the physical travel delay, and the
full wider support must fit the original bounded history. Modal-pressure traces
still sum to emitted pressure on the WAV clock. Both stereo channels receive
the same fixed extra latency; their geometric time differences are retained.
Allow for the extra latency when comparing output with older renders.

This removes one avoidable receiver-side loss. It does not supply missing
soundboard modes, improve the structural mesh or bridge impedance, calibrate
felt, add lid/room scattering, or establish perceptual similarity to a real
piano. Paired physical renders and listening remain necessary.

## Focused native checks

```sh
cargo test -p fs-couple --lib pcm_wav::baffled
cargo test -p fs-couple --example grand_piano
```

The new regressions cover complex fractional-delay response, independent
sinusoidal acceleration/pressure, integer-delay equivalence, causal support,
DC units, modal reconstruction, rejected-step rollback, buffer retention,
history bounds and activation in the played piano. Existing block/stereo
regressions exercise the same changed receiver through the audio stream.
