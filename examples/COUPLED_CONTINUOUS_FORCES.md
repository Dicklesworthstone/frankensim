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

This initial file path requires real instantaneous pressure-per-modal-velocity
maps (`H_IM = 0` in every mode row). Its signed pressure gains observe the actual
forced endpoint velocity. It does not interpret natural-frequency quadrature
as a causal broadband observer, add acoustic delays, or derive radiation from
air or geometry. Nonzero imaginary maps, nonlinear/contact v3-v5 schemas, and
trailing receiver blocks refuse explicitly. Ordinary files without the new
section keep their existing path. `--decimate` uses the existing explicit
pressure-observation filter when mechanics run above 48 kHz.
