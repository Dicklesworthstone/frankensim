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

Repeat `--modal FILE`, `--plate FILE`, `--bow FILE`, `--reed FILE`,
`--valve FILE`, `--flue FILE`, `--lip FILE`, `--voice FILE` and `--resolved FILE`
in the desired summation order. Every
supported modal performance version retains its existing bilateral, normal
contact, friction, force-port and preload semantics. Sources are independent:
this is acoustic superposition, **not mechanical coupling between files**.
Put components that exchange forces in one existing coupled/contact input.
`--bow FILE` loads an authored string, friction/body parameters and one canonical
bow trajectory; see `BOWED_PERFORMANCES.md`. `--lip FILE` and `--voice FILE`
use the same source files and physical pressure producers as the standalone
`music_render lip` and `music_render voice` commands. The voice's internal
mechanical substeps and both sources' mouth-flow and microphone histories stay
inside their original owners. They contribute exterior pressure, not bore or
subglottal pressure. Programmatic pressure producers use the same library
composition API.

All these source loaders receive the shared `--gas-component SPECIES.pack
NASA9.pack TRANSPORT.pack` registry, with explicit `--gas-water none|ID` and
`--gas-phase frozen|liquid-water` declarations. Each performance keeps its own
`ambient` record. Changing its pressure or composition, or changing only the
source pack data, therefore reaches its flow, duct losses and emitted pressure
before mixing. Every part with an exterior medium must agree on the listener's
complete gas state; incompatible pressure or mixture declarations refuse before
output creation. No ensemble option substitutes air for an unresolved source.

## Flue phrases alongside the other native sources

`--flue FILE` accepts the existing `frankensim-flue-performance-v1` network
source and `frankensim-flue-dg-performance-v1` DG source used by
`music_render wind`. Each complete source keeps its original jet card,
geometry, wall/key material packs, regional gas records, blowing/angle/key
control history, physical loads and selected microphone. Material and supplied
MIMO-radiation files resolve relative to that performance file, as they do for
the other native ensemble sources. The shared gas-component registry also
reaches both flue loaders. No note-triggered reconstruction interrupts their
jet, bore, material-memory or gesture histories.

For example, combine two existing DG performances covering the same 2048
mechanical intervals (1024 output samples). Each contributes the pressure from
its own explicitly described microphone:

```sh
cargo run --release -p fs-couple --bin music_render -- \
  ensemble /tmp/two-flues.wav --full-scale-pa 20 --decimate \
  --flue crates/fs-couple/examples/flue-dg-baffled.performance \
  --flue crates/fs-couple/examples/flue-dg-side-openings.performance --block 37
```

Replace either path with another complete, equal-duration source and repeat
any supported source flag to include reeds, bowed strings or other physical
parts. Matching exterior gas states remain mandatory. The graph source owns
its internal mechanical substeps and multirate microphone; the DG source keeps
its same-clock pressure producer and requires explicit output decimation when
its clock is above 48 kHz. Both report their original numerical microphone
delay separately from physical travel time.

A DG file explicitly selecting `observation inlet-pressure` still contributes
its internal diagnostic trace and is identified as having no exterior
microphone. Ensemble admission does not turn it into radiated sound.
The source jet's exploration-only status and every existing physical-model
limitation remain in the shared source receipt; composing a phrase does not
qualify a measured flute.

## Resolved fluid sources alongside structural and acoustic sources

`--resolved FILE` accepts the same `frankensim-resolved-dg-performance-v1`
description as `music_render resolved`. It advances the original resolved
fluid, physical boundary controls, live DG duct and material wall history.
The original microphone, exterior gas propagation and acoustic radiation load
remain attached throughout playback. For example, the two existing resolved
fixtures have matching observation windows and exterior gas:

```sh
cargo run --release -p fs-couple --bin music_render -- \
  ensemble /tmp/two-resolved-sources.wav --full-scale-pa 20 --decimate \
  --resolved crates/fs-couple/examples/resolved-dg.performance \
  --resolved crates/fs-couple/examples/resolved-dg-pressure.performance --block 37
```

Other source kinds can share the mix when their declared observation durations
and listener gas agree. Material source files resolve relative to each resolved
performance; thermochemical gas packs use the same explicit registry as the
other ensemble sources. Missing packs and incompatible listener media refuse
before either output is created.

The resolved fluid's physical time step comes from its mesh and material gas.
It is not rounded to the ensemble's clock. Its original conservative temporal
projection supplies pressure on the file's observation clock, followed by the
existing explicit decimator if requested. Metadata retains the exact physical
step, actual physical-step count, observation clock and receiver/filter delays.
The reported physical sample rate can therefore be nonintegral. The physical
prefix can lead its observation window by less than one coupled interval,
as in standalone playback; temporal projection still requires refinement.

A file selecting `observation inlet-pressure` contributes its explicitly
internal diagnostic pressure and has no declared exterior microphone. Resolving
the fluid and mixing its pressure does not add measured-source qualification,
room acoustics or mechanical coupling between independent performance files.

## Clocks, latency and duration

The output is 48000 Hz. Every source must already use this rate, or explicitly
request `--decimate` for an integer multiple supported by the existing
observation owner (ratios 2 through 16). Each source format retains its own
stricter rate limits. Upsampling, noninteger conversion and conversion flags
on an entirely 48 kHz ensemble refuse. No source clock is changed.

Only pressure is filtered. The original causal Blackman–Harris decimator
keeps its filter history and disclosed passband. A pure integer output delay
aligns each path to the largest combined source-microphone and decimator group
delay in the ensemble. Existing differentiating-microphone latency is read
from the original receiver, not recomputed or counted twice. For
192 kHz plus 48 kHz, the common delay is **44 output samples**: the decimated
path already has that delay; the 48 kHz path receives 44 samples of pure delay.
For exclusively 48 kHz inputs there is no additional sample-rate filter;
their original microphone filters may still require alignment delays.

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
Every individual loader's own stricter budget remains enforced (bow, lip,
voice, flue and resolved files are limited to 1 MiB). Bow controls are counted inside their hosted scheduler, not
misreported as an empty outer schedule. Reduction and
source-state memory are unchanged; pressure/PCM staging is callback-sized.
These are workload bounds, not a real-time throughput claim.
Lip and voice pressure histories execute inside their source owners; the
ensemble counts their retained initial values/events or pressure breakpoints
against the aggregate control budget without creating a new outer schedule.
Their mechanical and pressure clocks and sample counts are reported separately.

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
contact-plus-mesh WAV output, continuous graph/DG flue phrases and wall-material
substitution, flue/lip/voice microphone-delay alignment, arbitrary
source-backed gas changes, resolved-fluid clock and pressure preservation,
and incompatible listener-medium refusals. Native execution
is still required; source review is not compilation or physical validation.
