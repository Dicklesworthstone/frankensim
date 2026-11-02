# Continuous physical forces in a modal performance

The existing broadband modal renderer accepts prescribed sinusoidal loads at
its physical force ports. Run the supplied example with a built `music_render`:

```sh
music_render modal examples/modal_continuous_tones.performance continuous.wav --block 37
```

The example is an authored numerical fixture: one damped 500 Hz structural
mode, two signed actuator columns, three continuous force tones, a held load
released halfway through the record, and an instantaneous velocity receiver.
It does not reconstruct an instrument or claim a measured radiation operator.
The renderer retains physical pressure scaling; it does not normalize the peak.

Place the following source records after the ordinary `force` events and before
the required `broadband` receiver block:

```text
tone_limits MAX_TONES STORAGE_BYTES PREPARATION_WORK STEP_WORK
tones COUNT
tone VOICE_INDEX PORT_INDEX OMEGA_RAD_S COSINE_N SINE_N
```

Each row prescribes

`F(t) = COSINE_N*cos(OMEGA_RAD_S*t) + SINE_N*sin(OMEGA_RAD_S*t)`

in newtons at the original signed, work-conjugate port. Time starts at the
performance's first interval. Tones remain continuous throughout the declared
window and across callback boundaries. Frequencies can be harmonic, inharmonic,
repeated, or close enough to beat; amplitudes are physical loads, not output
audio partials. Existing force assignments continue as held loads alongside
them. Static preload uses only the pre-window held loads, preserving the
authored initial-state meaning.

`MAX_TONES` includes the held-load DC column and is at most 32 per affected
voice. Storage and work limits cover the affected voices together, with driver
buffers reserved and the remaining mechanics allowance divided equally.
Projection work is charged against the existing `compile_limits`. The original
mechanical clock, Nyquist guard, displacement, velocity, energy, pressure, and
receiver limits remain active; no tone is silently dropped.

The prepared continuous-force transition advances mechanics, physical force
work, viscous heat, and the attached broadband receiver together. Its source
phase advances only after an accepted step. The sidecar records the number of
physical tones, the phase origin, and the continuous-forcing integrator. A
complex constant legacy quadrature transfer is refused for these drives;
provide a causal broadband operator or a real instantaneous velocity gain.

## Continuous force attacks, releases and signed amplitude changes

`examples/modal_enveloped_tones.performance` adds sample-clock force envelopes
to the same model and receiver. Render it with the same command and a new output
path. The optional records sit after all `tone` rows and before `broadband`:

```text
envelopes COUNT
envelope TONE_INDEX POINT_COUNT
gain SAMPLE GAIN
```

Each envelope refers to a zero-based tone row in its original file order,
including when tones for different voices are interleaved. A tone can have one
envelope. Its first knot must be at sample zero, subsequent sample indices must
strictly increase, and every knot must lie in `0..=samples`. Gains are finite,
signed physical force multipliers. The gain varies linearly between knots and
holds its final value afterward. Tones without an envelope retain unit gain;
ordinary held-force events and static preloads keep their original meanings.

For example, knots `(0,0), (100,0), (400,1), (900,1), (1200,0)` prescribe silence,
a continuous attack, a plateau and a release. The force is the declared tone
multiplied by that continuous gain. It is integrated within each mechanical
interval together with motion, physical work, heat and the causal receiver.
The interval-end gain uses the same interval-start carrier basis; it is not a
phase reset. A zero-gain interval keeps the carrier running, so a later attack
has the phase implied by the original performance clock. Reversing gain sign
reverses the physical force continuously through zero. Turning off the force
does not erase existing mechanical vibration or receiver history.

The shared parser admits at most 65536 knots in total. Knot payload, curve
cursors and additional endpoint coefficient arrays consume the existing
`tone_limits` storage budget; endpoint scaling also consumes its step-work
allowance. The mechanics/filter preparation keeps its own bounded continuous
affine operator. Malformed knots and insufficient storage/work refuse; gains
are never replaced with held samples to fit the budget. Only accepted samples
advance carrier phase and envelope cursors. Files with no envelopes keep the
existing constant-amplitude propagation path. The sidecar records
`physical_force_envelopes` and identifies continuous affine-envelope integration.
