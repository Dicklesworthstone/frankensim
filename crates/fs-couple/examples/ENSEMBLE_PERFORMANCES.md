# Render independent physical models together

The ensemble command combines existing file-driven models without combining
or retuning their mechanical solvers. For example, render the 192 kHz free
striker/contact model alongside the 48 kHz geometry/material-derived plate:

```sh
cargo run --release -p fs-couple --bin music_render -- \
  ensemble /tmp/contact-and-plate.wav --full-scale-pa 1 --decimate \
  --modal crates/fs-couple/examples/free-striker-192k.performance \
  --plate crates/fs-couple/examples/plate-mesh.performance --block 37
```

The existing examples both cover exactly 4801 output intervals. The first
advances all 19204 mechanics samples; the second advances all 4801 plate
samples. Their force assignments stay on their original sample clocks.
The contact object, its contact law, and the plate's mesh, sections, supports,
force footprint and modal reduction all come from the existing source loaders.
No new oscillator, collision law, integrator or material-name preset is added.

Repeat `--modal FILE` and `--plate FILE` in the desired summation order. Every
supported modal performance version retains its existing bilateral, normal
contact, friction, force-port and preload semantics. Sources are independent:
this is acoustic superposition, **not mechanical coupling between files**.
Put components that exchange forces in one existing coupled/contact input.
There is no new bow/reed assembly-file parser; their programmatic pressure
producers can use the same library composition API.

## Clocks, latency and duration

The output is 48000 Hz. Every source must already use this rate, or explicitly
request `--decimate` for an integer multiple supported by the existing
observation owner (ratios 2 through 16). Each source format retains its own
stricter rate limits. Upsampling, noninteger conversion and conversion flags
on an entirely 48 kHz ensemble refuse. No source clock is changed.

Only pressure is filtered. The original causal Blackman–Harris decimator
keeps its filter history and disclosed passband. A pure integer output delay
aligns each path to the largest decimator group delay in the ensemble. For
192 kHz plus 48 kHz, the common delay is **44 output samples**: the decimated
path already has that delay; the 48 kHz path receives 44 samples of pure delay.
For exclusively 48 kHz inputs the common delay is zero and no filter is used.

The latency remains in the WAV. Pre-window pressure history is explicitly zero;
there is no lookahead, latency removal, tail flush or fictitious force history.
Choose source windows long enough to include the desired ringdown AND output
latency. The complete declared mechanical windows are simulated; the final
latency's worth of pressure remains in observation history, not in the WAV.
Retain a library renderer to continue its history rather than rebuilding it.

All source files must describe exactly the same physical duration and contain
complete output intervals. A longer source is not silently truncated and a
shorter one is not silently padded or frozen. Admission failures happen before
creating either output artifact. Callback size partitions work only.

## Observer and PCM scale

The caller must supply pressures for the **same observer and physical time
origin**, with compatible sign conventions. The loaders do not infer a room,
relative source placement, propagation delays or observer compatibility.
Filter-delay alignment does not replace those physical acoustic transfers.
All existing compactness, bandwidth, material and model-validity qualifications
still apply. In particular, the plate input retains its approximate compact
baffled observer and does not gain radiation reaction or room propagation.

`--full-scale-pa` is required and applies once, after summing pressures. Each
source file's PCM full-scale value is retained in metadata but is NOT a gain.
There is no peak normalization, per-part limiting, mixing attenuation or pitch
correction. The shared PCM16 writer counts actual clipped samples. Existing WAV
and provenance files are never overwritten.

The command admits at most 16 files of 4 MiB each, 128 source components,
8192 retained modes, 262144 compiled outer controls and 600 output seconds.
Every individual loader's own stricter budget remains enforced. Reduction and
source-state memory are unchanged; pressure/PCM staging is callback-sized.
These are workload bounds, not a real-time throughput claim.

## Library and checks

`pcm_wav::observation::ensemble::PressureEnsemble` consumes explicitly
constructed `DecimatedRenderer` parts. Use `Box<dyn PressureRenderer>` to
combine different producer types, including nonlinear impact observers and
finite `EnsembleRenderer` compositions. It validates the entire finite window,
exposes per-part clocks and delays, and implements the existing streaming
`PressureRenderer` interface. `render_pressure_pcm16` supplies callback-boundary
cancellation with an exact emitted prefix; resume the same mixer and stream.
A physical refusal poisons the mix; discard the failed callback, not just its
last source. Source and filter state may have partially advanced.

Focused native checks:

```sh
cargo test --release -p fs-couple --test pressure_ensemble --test bowed_ensemble
cargo test --release -p fs-couple --bin music_render --test music_render_ensemble
```

These tests compare real physical trajectories, source-clock controls, separate
filtering and indexed latency alignment against the mixed path, including actual
contact-plus-mesh WAV output and a material-change comparison. Native execution
is still required; source review is not compilation or physical validation.
