# Preserve upper partials through acoustic propagation

The geometric piano's mono and stereo Rayleigh microphones select the shared
`BaffledPressure::new_multirate_differentiating` receiver, combining bandlimited
acceleration with fractional-flight interpolation. Both channels still
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
44 samples and bandlimited differentiation another 8, so
`filter_delay_samples()` reports **84** for the piano,
excluding distance-dependent travel. `interpolation_delay_samples()` reports
32 separately. Do not add that 32 twice. The acceleration-input
`new_multirate_bandlimited` comparison image retains its 76-sample delay.

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

## Refine the cubic soundboard's radiated field

Flat edge-cubic boards now accept the same explicit spatial refinement as P1
boards. For a geometric pressure render, combine:

```text
--edge-cubic-board-mass --acoustic-refinement-levels 2
```

Levels 0 through 3 subdivide each structural triangle into `4^level` equal-area
cells, using three positive integration points per cell. Every new point
evaluates the original cubic field, including its physical slope coordinates.
The structural eigenproblem, mass, damping, bridge projections and modal
coordinates are unchanged. This allows finer integration of propagation phase
and spatial cancellation without rebuilding the mechanical mesh or replacing
the cubic field with nodal P1 interpolation. Level 0 preserves the original
quadrature path; the total 120,000-point receiver limit still applies. Crowned
shells and diagnostic volume output remain outside this flat-panel control.

An independent spatial integration regression prescribes a quadratic normal
acceleration on an oblique triangle, observes two off-axis receivers, and checks
1.5, 4, 8 and 12 kHz. It compares the actual cubic shape and production
quadrature with a tensor-Simpson integral that evaluates the prescribed
polynomial and the finite-distance Rayleigh kernel directly.

| Refinement level | Complex relative RMS error across those cases |
| --- | ---: |
| 0 | 189.97% |
| 1 | 79.66% |
| 2 | 10.85% |
| 3 | 0.194% |

The independent reference changed by 8.30e-8 relative RMS when its resolution
doubled. This is a manufactured spatial-integration result, not a piano
recording, universal error bound or selected production tolerance. A signed
field can cancel strongly at a receiver; relative errors can exceed 100%.
Separate tests compare all nine signed cubic displacement/slope area moments
with the existing owner's analytic integral and check positive weights and
panel area. All six tests in the actual-source native closure passed with the
pinned compiler, including the unchanged cubic-owner tests. The added full
board-preparation equivalence test and a complete piano render remain unrun.

## Focused native checks

```sh
cargo test -p fs-couple --lib pcm_wav::baffled
cargo test -p fs-couple --example grand_piano
cargo test -p fs-couple --example grand_piano acoustic_quadrature
cargo test -p fs-couple --example grand_piano cubic_acoustic_refinement
```

The new regressions cover complex fractional-delay response, independent
sinusoidal acceleration/pressure, integer-delay equivalence, causal support,
DC units, modal reconstruction, rejected-step rollback, buffer retention,
history bounds and activation in the played piano. Existing block/stereo
regressions exercise the same changed receiver through the audio stream.
