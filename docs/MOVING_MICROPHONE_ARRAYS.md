# Independently moving microphones on one physical frame clock

`pcm_wav::observation::moving::frames::MovingPressureFrameRenderer` composes
already prepared `MovingCompactRenderer` observers. Each microphone keeps its
actual receiver trajectory and its own retarded source times, Doppler admission,
clearance and direct/reflected propagation paths. It is not a panner or a stereo
copy of a single pressure signal. No velocity-dependent amplitude gain is added.

Construct the individual observers through the existing moving-source APIs and
borrow the same immutable accepted source histories when observing the same
physical scene. All channels must have the same receiver-time origin, sample
rate, finite horizon and next sample index. Different bandwidth declarations,
clearance limits and root tolerances remain with their individual admitted paths.

The frame constructor takes a fixed-size array of existing observers and bounded,
unique channel IDs. It reuses their interpolants and callback buffers with no
additional allocation. `MovingPressureFrameBudget` checks aggregate scalar output
values, source-path/receiver samples, root-iteration work and logical retained
payload for the full original horizon. This aggregate check complements, rather
than replaces, each observer's original admission. Channel IDs are the mapping to
retain beside exported WAV bytes; no physical speaker layout is inferred.

```rust,ignore
use fs_couple::pcm_wav::observation::{DecimatedFrameRenderer, moving::frames::*};
use fs_couple::pcm_wav::stream::{Pcm16WavStream, render_pressure_frames_pcm16};

// left and right already describe their own actual accepted microphone paths,
// observed at 96 kHz. Both start at sample zero with complete source coverage.
let array = MovingPressureFrameRenderer::new(
    ["left microphone", "right microphone"],
    [left, right],
    array_budget,
    &mut cancelled,
)?;
let mut audio = DecimatedFrameRenderer::new(array, 96_000, 48_000, 32)?;
let frames = audio.output_frames_for(source_frames)?; // refuses fractional groups
let mut wav = Pcm16WavStream::new_interleaved(sink, 48_000, 2, full_scale_pa, 32)?;
render_pressure_frames_pcm16(&mut audio, &mut wav, &gate, &mut [0.0; 64], frames)?;
```

## Pressure, clocks and continuation

Outputs are interleaved **frames**, each containing one Pa value from every
microphone. Channel order and signs are preserved; the writer applies one common
Pa-to-PCM scale. Sample times remain the original integer-indexed receiver clock,
independent of callback partitioning. An array can also join already observed
mono streams at an exactly aligned prefix; `into_channels()` returns the original
observers at their common accepted position without cloning or resetting them.

Array callbacks stage every channel before publishing any output or advancing
any cursor. A late physical failure or cancellation leaves the complete callback
buffer and every channel cursor untouched. A retry overwrites staging scratch.
Read-only channel access prevents changing a child cursor behind the array.

## Causal multichannel decimation

The reusable `observation::DecimatedFrameRenderer` accepts any native
`PressureFrameRenderer`, including these moving arrays and harmonic arrays. It
uses the existing `Decimator`, its original coefficients and independent channel
delay lines. Only pressure is filtered, not mechanical motion, source strengths,
events or forces. One rate group is buffered, independent of callback size and
duration; sources admitting one-frame callbacks are supported. Rate one delegates
the original callback without extra arithmetic.

The input rate must already resolve each actual receiver's Doppler and amplitude
modulation. Filtering cannot repair aliasing at that input rate. The existing
filter's passband is 0 through 0.45 of output rate; the transition and rejection
properties are those of `pcm_wav::decimate`, not an ideal brick-wall filter.
`info()` retains the common causal group delay and first output group's source
index. Delay is not compensated, and no missing startup history, padding or tail
is synthesized. A new filter refuses a source that has already advanced.

Cancellation of the existing PCM stream occurs between complete callbacks. Keep
the **same** array, decimator and writer to resume, preserving every channel's
filter history. Unlike the array's stateless observation transaction, a failed
source/filter callback inside a generic decimator may have consumed a prefix.
That adapter is then poisoned and refuses continuation; discard the failed
callback's output. It never pretends it can roll back arbitrary source mechanics.
Input-shape and horizon admission failures occur before advancement and remain
retryable.

The moving-source, material, acoustic-linearity, finite-source compactness and
sampling-convergence obligations are unchanged. Ideal microphones are not a head
scattering model, and these adapters add no HRTF or measured-SPL calibration.
