# Continuous force tones on a coupled modal assembly

`music_render modal examples/coupled_continuous_tones.performance OUT.wav --block 512`
renders the supplied version-2 mechanical assembly under physical sinusoidal
actuator forces, alongside the existing held-force events. This example is an
authored numerical assembly, not a validated instrument or measured radiation
model. Its explicit `full_scale_pa` controls PCM scaling; output is not normalized.

Append this section after the ordinary v2 `force` records:

```text
coupled_tone_limits MAX_TONES MAX_DIMENSION STORAGE_BYTES MULTIPLY_TERMS
tones COUNT
tone VOICE_INDEX PORT_INDEX OMEGA_RAD_S COSINE_N SINE_N
```

Without gain or frequency curves, each tone applies `COSINE_N*cos(OMEGA_RAD_S*t) + SINE_N*sin(OMEGA_RAD_S*t)`
newtons through the voice's original signed force-port shapes. Time starts at
the beginning of the performance, and the tones continue over its entire
declared window. Frequencies can be inharmonic or arbitrarily off an FFT bin.
Repeated frequencies interfere coherently across ports and bodies; nearby
frequencies remain distinct. They are prescribed external physical loads, not
authored output harmonics or acoustic partials.

The complete spring/dashpot network and every original retained modal state
advance together under continuous forcing. The same interval supplies external
work, each physical heat channel, and the accepted pressure sample. Held-force
assignments still apply before their numbered sample and form a DC input during
that interval. Callback boundaries never reset phase. A refused physical step
publishes no new phase, event cursor, mechanical state or work/heat record.

`static-preload` keeps its original v2 meaning: all components must request it,
and the complete static spring equilibrium uses only the pre-window held forces.
The tones begin at time zero; their steady state is not silently assumed.

The tone count limit includes a DC column and is at most 32. Dense dimension is
at most 256, storage at most 134217728 bytes, and preparation work at most
100000000000 scalar multiply-add terms. Source records, frequencies and all phase
buffers are reserved from `STORAGE_BYTES` before the remaining allowance is
passed to the complete dense physical owner. Existing mechanical, callback,
event, projection and guard-frequency limits remain binding. Both quadratures
charge the original `compile_limits` projection budget. Larger assemblies or
tone sets can refuse; no coordinate is dropped to meet a budget.

Without an explicit receiver, the file requires real instantaneous
pressure-per-modal-velocity maps (`H_IM = 0` in every mode row). Their signed
pressure gains observe the actual forced endpoint velocity. A nonzero imaginary
constant map refuses because natural-mode quadrature cannot describe arbitrary
continuous forcing. Ordinary files without the new tone section keep their
existing path. Nonlinear/contact v3-v5 schemas remain outside this linear owner.
`--decimate` uses the existing explicit pressure-observation filter when
mechanics run above 48 kHz.

## Physical force attacks, releases and re-attacks

`music_render modal examples/coupled_enveloped_tones.performance OUT.wav --block 512`
adds continuous gain envelopes to the same coupled broadband example. Gain
multiplies the physical actuator force before mechanical evolution; it never
scales the emitted audio. A released source leaves mechanical, material,
receiver and propagation histories intact, so physical ringdown remains.

Place an optional knot section after the tone records and before the receiver:

```text
envelopes COUNT
envelope TONE_INDEX POINT_COUNT
gain SAMPLE GAIN
```

Repeat `gain` for each point and `envelope` for each declared curve. Tone indices
refer to the original zero-based `tone` row order, including separate rows at
the same frequency. Each curve starts at sample zero; subsequent samples must
strictly increase and lie within the declared performance horizon. Gains are
finite signed dimensionless force multipliers. They interpolate linearly
between knots and hold their final value. Undeclared tones keep unit gain;
`envelopes 0` keeps the original constant-tone path. Held DC events are separate.
The common parser caps the aggregate at 65536 knots and reserves their payload
from the existing source byte allowance before allocation.

The affine preparation integrates both endpoint force coefficients on each
sample interval using four rotating source coordinates per distinct modal
tone. Both endpoints share the interval-start carrier basis. Original work,
each heat channel, and the causal pressure filters use that same trajectory.
Independent unscaled carrier phases continue through zero gain and negative
gain; a re-attack never divides by the old gain or restarts its phase. Curves,
events and carrier phase advance only after the complete observation accepts.
The original force guard bounds the entire interval using endpoint coefficient
norms. Dense dimensions and work grow explicitly; no modes or tones are dropped.

Piecewise-affine interpolation is exact for the declared force history, within
floating-point error. Approximating a smooth physical control requires refining
its knots independently of the mechanics clock. Carrier-frequency admission
does not certify that envelope sidebands or sharp knot transitions are
bandlimited; an adequate observation clock and the existing decimation filter
remain necessary. An amplitude-only section keeps the source frequency fixed.
Nonlinear contact, jet, bow and hammer excitations retain their physical owners.

## Phase-continuous physical force-frequency sweeps

`music_render modal examples/coupled_swept_tones.performance OUT.wav --block 512`
combines gain envelopes, three independently swept physical forces, the same
coupled mechanical assembly, and its causal signed receiver. Two sources share
one prepared carrier and still retain distinct frequency histories. This is a
prescribed shaker/actuator load: structural frequencies, tension and material
properties retain their original values and physical update rules.

After the optional gain envelopes and before the receiver, supply:

```text
frequency_sweeps COUNT
sweep TONE_INDEX POINT_COUNT MAX_FORCE_ERROR_N
omega SAMPLE OMEGA_RAD_S
```

Repeat `omega` for each knot and `sweep` for each declared source. Indices refer
to the original tone rows. Each curve begins at zero, increases strictly within
the performance horizon, interpolates angular frequency linearly, and holds
its final frequency. Frequencies are finite, nonnegative, and below the
**original coupled model's guarded Nyquist limit**. Zero frequency holds the
current source phase. An undeclared row retains its original carrier frequency;
`frequency_sweeps 0` preserves the existing source path exactly. Gain and
frequency knots may have different sample indices.

Phase advances by the exact integral of each linear frequency interval,
`h * (omega_start + omega_end) / 2`, including while gain is zero. The force is
reconstructed with affine cosine/sine coefficients around the original fixed
carrier. The existing exact propagator integrates this reconstructed force,
its work, each dissipation channel and the causal receiver together. No matrix
exponential is rebuilt during a sweep, and there is no basis or resonant-state
transfer. A modulation proposal is computed once per source interval, then
reused for all signed modal projections. Phase, gain/frequency cursors and held
events publish only after the complete physical and acoustic step accepts.

`MAX_FORCE_ERROR_N` is a positive absolute allowance for **each original physical
tone row**, before modal projection. Admission bounds the interpolation defect
over the complete declared history. For original quadrature magnitude `R`, gain
bound `G`, maximum gain change per sample `dg`, maximum frequency change per
sample `dw`, and maximum detuning angle `D = h * max|omega - carrier|`, it uses:

```text
force_error <= R * (2 * dg * D + G * (h * dw + D * D)) / 8
```

The bound follows from the second derivative of the complex gain envelope.
Using independent extrema across the history is conservative when the largest
gain and frequency changes occur in different segments. The parser refuses an
insufficient allowance before returning a renderer. A finer mechanics clock or
a fixed carrier closer to the sweep reduces the reconstruction error; changing
the prepared carrier preserves the intended physical frequency history.

The bound covers force reconstruction in real arithmetic. Floating-point phase
accumulation, acoustic fit error and mechanical response sensitivity have their
own scopes. Multiple coherent force rows can add their errors; the reported
largest per-row bound is not a total pressure error. The CLI records
`physical_force_frequency_sweeps`,
`maximum_tone_force_reconstruction_error_n` and `source_reconstruction`.

Sweep knots and retained source metadata use the original tone byte allowance;
aggregate sweep knots are capped at 65536. The existing affine-source dimension
and work limits still apply. Source-frequency admission does not certify
sideband bandwidth or sharp knot transitions. Choose the original mechanics
and observation clocks accordingly, using the existing anti-aliasing decimator
when exporting a higher-rate mechanical run.

## Causal broadband pressure from the complete coupled trajectory

`music_render modal examples/coupled_broadband_tones.performance OUT.wav --block 512`
adds a causal receiver with signed frequency-dependent gains and a three-sample
propagation delay. The mechanical assembly, physical tones, preload and held
events are identical to the first example. Both examples are authored numerical
fixtures; their transfer samples do not claim measured or BEM-validated radiation.

Append one complete receiver after the last tone:

```text
coupled_broadband MAX_DIMENSION STORAGE_BYTES MULTIPLY_TERMS
transfer PRODUCER_ID DELAY_S PRODUCER_REL_ERROR ORDER ITERATIONS DELAY_HALF MAX_REL_ERROR TRAINING_COUNT HELDOUT_COUNT
frequency OMEGA_RAD_S H0_RE H0_IM H1_RE H1_IM ...
withheld OMEGA_RAD_S H0_RE H0_IM H1_RE H1_IM ...
```

Supply exactly `TRAINING_COUNT` frequency rows followed by `HELDOUT_COUNT`
withheld rows. Each row contains every mode, first all modes of voice 0, then
all modes of voice 1, and so on. The complex coefficients are receiver pressure
per generalized modal velocity, in Pa s/(m sqrt(kg)), under `exp(-i*omega*t)`.
They include the declared propagation delay. Preserve the same signed modal
basis and normalization as the mechanics; transfer rows are never magnitudes
or force-to-pressure coefficients. Training frequencies must increase, and the
independent withheld grid must be disjoint. All frequencies must be positive
and below the mechanical Nyquist frequency.

The receiver fits a real stable proper continuous rational transfer for each
modal velocity. Signed direct and residue contributions are retained. The
explicit delay is removed before fitting and restored using the existing
sampled-delay kernel. The held-out gate checks the actual fitted operator with
that kernel. `MAX_REL_ERROR` bounds these sampled relative errors; it is not
a continuum enclosure. The producer's separate discretization estimate and
the fit/delay controls remain in the CLI provenance as `EstimateOnly`.

Filter states follow the original within-sample source trajectory through the
shared exact exponential. Mechanical modes, Maxwell histories and reciprocal
acoustic states stay coupled; the filter is not driven by a held endpoint
velocity. The complete physical owner still computes its original work and
each heat channel. This receiver replaces legacy pressure maps and adds no
reaction, physical storage or duplicated radiator. Its pressure is subject to
the original network pressure ceiling. A rejected observation leaves source,
filter, delay, phase and event cursor at their last accepted values.

The receiver's filter and delay histories start at zero. Static mechanical
preload does not imply a pre-existing acoustic history. Integer delays use
exact sample lags. Fractional delay is an explicit bandlimited approximation;
its declared centered kernel must fit within a causal nonnegative lag history
and pass the same held-out error gate. This does not certify that transients
are bandlimited. Output rate changes still use the existing decimation filter.

The separate `coupled_broadband` allowance covers conservative fitter workspace,
retained rational models, continuous filter matrices/states, delay history and
exact observer preparation. Its dimension cap includes the complete source,
filter and forcing coordinates and is at most 256; storage is at most 128 MiB
and exponential preparation work at most 100 billion multiply-add terms. Fit
order and iterations are independently bounded by the transfer controls. The
original physical source keeps its own `coupled_tone_limits` budget. No mode,
tone or fitted pole is dropped to meet either allowance.

The second example's delay-free analytic transfers are
`H0 = 0.25 + 1.5*900/(900 - i*omega)` and
`H1 = -0.4 - 0.6*900/(900 - i*omega)`, multiplied by
`exp(i*omega*3/48000)` for the declared delay. This deliberately exercises
negative direct pressure gain, coherent signed summation, independent withheld
frequencies and the same physical tone/event path as the first example.
