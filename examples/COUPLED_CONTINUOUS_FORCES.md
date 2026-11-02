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

Every tone applies `COSINE_N*cos(OMEGA_RAD_S*t) + SINE_N*sin(OMEGA_RAD_S*t)`
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
