# Physical scale and stretch tuning in live playback

Live playback can now consume the same SI scale CSV and octave-tension tuner
as the offline piano. Preparation happens before `LivePiano` or `PianoAudio`
builds its reference Bank, partial tables, static bridge loading, hammer and
felt geometry. The sound comes from those physical courses, not a pitch-shifted
recording or independently assigned oscillator frequencies.

The original embedded equal-tempered scale remains the default. Neither this
selection nor imported geometry is promoted to a measured instrument model.

## Browser selection

After rebuilding the piano WASM through the existing build path, pass this
object as the existing `options` in either piano worker's init message:

```js
const options = {
  boardCoordinates: 32,
  substeps: 4,
  stretchTuning: {
    concertAHz: 440,
    bassOctaveWeights: [0, 1, 1, 0],
    trebleOctaveWeights: [1, 0, 0, 0]
  }
};
```

`concertAHz` is required and must be finite in 430..450 Hz. Each optional weight
array has four nonnegative entries for lower:upper partial matches **2:1, 4:2,
6:3, 8:4**, with at least one positive entry. An omitted array means `[1,0,0,0]`.
These example weights are an explicit local octave objective, not a calibrated
recommendation for a particular piano. C4..B4 stay anchored; the bass and treble
are solved outwards through supplied octave neighbours.

An optional positive `maxPartialHz` limits the neighbour target partials used
by that fit. It is not the audio bandwidth, a mode limit, or a soundboard
accuracy certificate. An octave with no retained positive-weight pair refuses.
The original tuning owner remains in
`crates/fs-couple/examples/grand_piano/stretch_tuning.rs` and `partial_tuning.rs`;
the live library compiles these exact files through `refmods`, without copying
their equations or introducing another frequency law.

Both workers already share `createPianoEngine`. Their board-selection probes,
separate benchmarks, and chosen engine copies therefore retain the same scale,
tuning and substep selection. No new score or note event type is required.
Clock/band limits and 1x versus higher-rate failure policies are unchanged.
Older WASM builds lacking `withPreparation` reject a requested setup instead
of silently playing the original equal-tempered scale. `engine-piano/build-wasm.sh`
builds and installs the WASM.

## Play an already prepared or authored scale

Produce a CSV with the existing `piano_tune` CLI or author the canonical SI
format. Load its text into `options.scaleCsv`. For a file served by the app:

```js
const response = await fetch('/voiced-scale.csv');
if (!response.ok) throw new Error(`scale fetch failed: ${response.status}`);
const options = {
  boardCoordinates: 64,
  substeps: 2,
  scaleCsv: await response.text()
};
```

Without `stretchTuning`, supplied tensions are retained bit for bit. This
avoids an extra retuning pass discarding a previously prepared scale. Supplying
both controls explicitly means **import first, then retune**. The other Course
fields remain authored: lengths, linear densities, bending rigidities, unison
spread, strike geometry, hammer masses, felt dimensions and duplex lengths.
The existing CSV parser owns physical admission, key uniqueness and ordering.
Sparse scales are legal; sparse tuning still requires connected octave chains.
The native and browser ingresses bound scale text at 256 KiB before parsing.

The board file's modal coordinates, full damping and acoustic surface are not
changed by scale preparation. The live engine then rebuilds its string-loaded
mechanical basis from the resulting courses; it does not reuse a previous
scale's static string loading. This is cold construction, not a hot retuning
of an already running engine or a continuation of its hammer histories.

`info().scalePreparation` identifies imported text, course count, selected
reference/weights/band, octave-pair count and largest weighted RMS mismatch.
The RMS is the shared local frequency objective, not loudness or perceived
roughness. Existing key geometry and partial-frequency observations read the
actual prepared courses. Copies retain the setup report and physical state.
Original constructors with no preparation retain their previous metadata.

## Native composition

```rust
use fs_live_piano::{Config, Description, PianoAudio, preparation::{prepare, Tuning}};

let mut description = Description::embedded()?;
let policy = Tuning::new(442.0, [0.0, 1.0, 1.0, 0.0], [1.0, 0.0, 0.0, 0.0], None)?;
let report = prepare(&mut description, None, Some(policy))?;
let piano = PianoAudio::new(&description, Config::live(48_000), 4)?;
```

`prepare(&mut description, Some(csv), None)` imports without tuning. Neither
option preserves even the original course allocation. Both operations are
transactional together: invalid text, missing neighbours, inadmissible tension
or an empty target band leaves the original Description unchanged.

## Atmosphere loading of source-material strings

`PianoAudio::new_with_materials` and the browser's canonical material-pack path
use the supplied atmosphere for string drag as well as microphone propagation.
The outermost radius comes from each actual concentric section, including all
of its declared layers. Its solid line mass remains the material integral.
The shared `fs-couple::air_path::oscillating_cylinder_air_resistance_per_length`
owner evaluates the small-amplitude transverse-cylinder resistance at each
actual partial frequency. Resistance divided by line mass supplies the gas-only
velocity-loss rate. Pressure, temperature or composition can therefore alter
the solved string/bridge motion and decay as well as received pressure.

The sampled rate is frozen during construction and enters the existing passive
loss flow. Sounding strings and the silent-string bank use the same prepared
coefficients at every supported mechanical clock. Physical gas loss is included
once in the original `string_loss_j` account. Source Kelvin–Voigt bending loss
keeps its separate constitutive origin; imported solid coefficients are not
modified to impersonate the atmosphere. The same admitted gas sets Rayleigh
density and propagation speed. Copies retain both that gas and the drag rates.
`string_gas_drag_rates_per_s(key)` exposes the granted speaking partials' gas
coefficients; `info().stringGasLoading` names the selected reduction.

The silent bank keeps its existing explicit linear-bridge trajectory. Gas drag
acts on absolute modal velocity there too: in the bank's relative-to-chord
coordinate it supplies the constant force `-gas_rate * beta * bridge_velocity`.
The original damped transition advances that forced equation by an equilibrium
shift. Gas work and pad work use their distinct absolute and relative velocity
integrals, and the same transfer enters the existing coupling-work account.

This is the existing sampled, linear, unconfined, circular-cylinder approximation
on retained absolute-velocity modes. It does not add fluid inertia, viscous
memory, inter-string or end interactions, or the unresolved chord's fluid load.
It does not retune elastic frequencies. Finite-amplitude, rarefied and confined
flows need their own admitted models. The generic causal `CylinderFlow` remains
a separate model whose fluid history needs explicit transfer through the live
engine's changing voice basis before it can be used here.

Effective scale CSVs do not declare outer radius or separate their aggregate
loss into material and gas contributions. `new_with_gas` retains its explicit
receiver-only meaning for that path; radius is never guessed from density and
an additional air loss is never stacked onto the aggregate preset. Source-pack
construction supplies the complete geometry and loss separation needed above.

Focused additions to `--test scale_preparation` check independent SI resistance
values, pressure/composition substitutions, actual bridge and decay response
after dividing out receiver density gain, silent-bank use and cloned playback.
Native execution of these additions remains pending.

## Scope and verification

The tuner fits nominal pinned-string partials, not the coupled board/string
poles, time-discretization shifts, a global temperament, or a measured Railsback
curve. Unstruck strings still have the live engine's existing static image;
this does not add full sympathetic resonance, calibrated felt/damping data,
spatial soundboard convergence, or a listening-quality guarantee.

```sh
cargo test --manifest-path apps/instrument-lab/engine-piano/Cargo.toml --test scale_preparation
node --test apps/instrument-lab/tests/piano-preparation.test.mjs
```

Six native tests cover unchanged defaults, exact shared-solver parity,
physical-field preservation, import/retune order, atomic failure and actual live
partial/pressure composition at 1x/2x/4x; all six pass (macOS, release). The
Node constructor-selection tests use a recording mock, not a WASM piano. With the
built WASM in Node, `stretchTuning: { concertAHz: 440 }` prepares 88 courses
through 76 octave pairs (largest weighted RMS mismatch 4.5e-13 Hz) and plays
finite output, as does 442 Hz with 2x. No new numerical solver, filter,
dependency, lockfile or workflow is introduced. Listening remains outstanding.
