# Higher-rate mechanics, explicit 48 kHz observation

```bash
cargo run -p fs-couple --bin music_render -- \
  modal crates/fs-couple/examples/free-striker-192k.performance \
  /tmp/free-striker-192k.wav --decimate --block 37
```

This runs the supplied free striker and two-mode receiver at **192 kHz**, then
filters the observer pressure before producing **48 kHz** PCM. It uses the same
contact, modal integrator, force scheduler and WAV encoder as the ordinary path.
The causal decimator is the former grand-piano implementation, now shared as
`fs_couple::pcm_wav::decimate::Decimator`; the piano example re-exports it rather
than retaining a duplicate filter. Its coefficient arithmetic is unchanged.

The example has 19204 mechanics samples and therefore exactly 4801 output
samples. A force starts before mechanics sample 1201, **inside** an output
interval, not rounded to the next audio sample. Every force and contact reaction
acts on mechanics before the pressure observation. Changing `--block` changes
only the output callback partition, not the event clock or the filter history.

## Clocks are declared, not inferred

Without `--decimate`, the command continues to require a 48000 Hz model. With
`--decimate`, the input must declare a strictly higher integer multiple of
48000 Hz. The current file-reader rate cap admits 96000, 144000 and 192000 Hz;
the native adapter supports integer ratios from 1 through 16. No arbitrary-rate
resampling, upsampling or automatic rebuilding of a lower-rate model occurs.

The file's `sample_rate_hz`, `samples` and all `force SAMPLE` indices are
**mechanical** quantities. To compare the same physical performance at four
times the rate, multiply its horizon and intended event indices by four. The
new example intentionally offsets the first event by one additional mechanics
sample to exercise sub-audio event timing. Model parameters, port shapes,
forces, pressure scaling and numerical tolerances retain their original units.

The mechanical horizon must be exactly divisible by the conversion ratio.
Otherwise the command refuses **before creating either output artifact**. It
never drops a final fraction of an audio interval or supplies synthetic motion
past the declared horizon. A short final callback containing complete output
samples is retained, as in the original command.

## Filter delay and boundaries

The filter starts with zero pressure history. It is causal: there is no lookahead
or reflected edge. The first emitted frame ends at source index `ratio - 1`;
mechanical samples themselves are observations after the corresponding step.
The group delay is retained, not compensated: 40 output samples at ratio 2,
80 at ratio 3, and 44 at ratio 4. No extra filter tail is flushed after the
declared mechanical window. Extend the physical simulation horizon explicitly
when a longer observed decay is needed; do not treat omitted tail energy as a
mechanical dissipation term.

Only pressure is filtered. The mechanical state and its storage/work/loss audit
are untouched. An observer filter is not a material loss, actuator smoothing,
normalization, or a replacement acoustic radiation law. Peak/RMS and clip counts
in the WAV provenance describe the **filtered output pressure**. The nested
`modal_input.observation` record includes source/output clocks, source horizon,
ratio, filter profile, delay, initial-history and tail policies. The exact input
hash still binds the authored mechanical performance. Ordinary unfiltered
48 kHz renders do not acquire an observation record or a changed PCM algorithm.

The inherited filter design targets 0..0.45 of output rate as its retained band,
with foldover rejection evaluated from 0.55 of output rate. The transition near
Nyquist is not a certified flat audio band. Its finite-grid response tests are
not a proof between frequency samples. Nonlinear contact can generate energy
above the **mechanics** Nyquist limit; no later filter can undo aliasing already
created there. Refine the physical step and compare the relevant observables.
This path is not an alias-free, physically validated or hard-real-time claim.

## Native streaming and continuation

`pcm_wav::observation::DecimatedRenderer` owns the scheduled mechanics and complete
filter history. `pcm_wav::stream::render_pressure_pcm16` feeds it into the same
incremental PCM writer. Cancellation occurs between output callbacks; the
in-flight mechanics/filter callback drains to the sink. Resume with the **same
adapter and stream** and a fresh cancellation gate. Constructing a new adapter
from a partially advanced source is rejected because the missing pressure
history cannot be reconstructed from its current modal state alone.

An invalid output shape or expanded-clock request is rejected before motion.
A physics or filter failure may leave a callback prefix advanced and therefore
poisons the adapter. A partially failed I/O write retains the stream's existing
poisoning rule. Previously completed output callbacks can be finalized explicitly;
they are not labelled a complete requested performance.

Focused native checks:

```bash
cargo test -p fs-couple --lib pcm_wav::decimate
cargo test -p fs-couple --test decimated_render --test music_render_decimated
cargo test -p fs-couple --bin music_render
cargo test -p fs-couple --test wav_stream --test music_render_free_mass
```

These native checks were not executed in the authoring environment, which lacks
Cargo/rustc. The separate numerical references do not execute the Rust code.
