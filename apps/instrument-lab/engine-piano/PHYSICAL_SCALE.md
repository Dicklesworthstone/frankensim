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

## Nonlinear string motion in browser playback

`options.stringStretching` accepts the complete reference `--string-stretching`
file as text, bounded to 64 KiB. It selects the existing native `linear`,
`stretch` and `waves` models per key. This is separate from `stretchTuning`,
which prepares static string tensions. Rust validates the specification against
the final imported/tuned courses and constructs its strain storage and forces
before any note plays. The same selection reaches standalone and scheduled
workers, board timing probes, playback copies, and the 1×/2×/4× physical clocks.
When using `stringMaterials`, the specification's axial-rigidity field must be
`section`, so the existing parser takes EA from the same source-bound section
used by the other string physics and refuses an independent numeric value.

```js
const response = await fetch('/physical-string-extension.fspx');
if (!response.ok) throw new Error(`string specification fetch failed: ${response.status}`);
const options = { substeps: 4, stringStretching: await response.text() };
```

Use `stretch` for the existing averaged geometric-extension model. `waves`
additionally needs board data containing geometry-derived axial participation
rows for its keys; a board without that data refuses the request. No axial
projection is inferred from the vertical bridge row. The model remains the
declared live reduction, with its numerical work in the physical ledger.

The axial coordinates consume the shared exact linear-force transition from
`fs-modal`. Their viscous heat is integrated along the same continuous motion,
including during resonant buildup, rather than estimated from endpoint powers.
The live path accepts the file's complete damping-ratio range `[0, 1]`, including
critical damping. Transitions and heat coefficients are prepared before playback;
the axial update allocates nothing. A refused axial trial preserves that segment's
accepted coordinates and force history and is reported by `take_error()`.
This segment-level refusal does not roll back the complete piano timestep.

One extension/longitudinal state now belongs to each physical retained string,
including an unplayed sympathetic string. Its geometric strain uses the union
of its active prefix and remaining bank tail. Both the transverse and axial
bridge reactions return to the same board, and axial viscous heat is counted
once. The state continues while a key is released or its active voice retires.
When ownership changes without changing that union, the exact same axial
displacement, velocity, previous force and split-energy baseline are retained.

A changed partial budget projects only the transverse prefix. The unchanged
axial sine coordinates keep their physical displacement and velocity; the new
current strain starts the next linearly interpolated force. The actual change
of linear, chord and nonlinear storage is accounted at that boundary. A net
removal goes to `retired_j`; a net increase goes to signed representation work
in `coupling_j`. Nonlinear cross terms can increase extension energy when a
partial is removed, so this projection is never called physical heat and no
amplitude normalization is applied. Failed rebuilds restore the nonlinear
state, drive history and ledger along with the existing linear state.

`LivePiano::string_extension_energy_j(midi)` observes the course's extension
storage and axial kinetic energy across both owners. Total stored energy
includes that physical state once. The conservative nonlinear coupling still
uses the declared explicit kick; its time-discretization work remains in
`coupling_j`. Shared ownership does not turn it into an implicit or exactly
energy-conserving nonlinear integrator. Mode-band and timestep convergence
remain separate requirements.

The optional final argument on the WASM preparation constructors carries the
unchanged text. The factory checks `supportsStringStretching()` first, so older
WASM cannot silently ignore the selection. Rebuild using
`engine-piano/build-wasm.sh`; omission retains the original constructors and
physics. `info().stringStretchingSelected` reports the explicit selection,
including a supplied all-linear file.

## Per-key damper preparation

`options.damperSpecification` accepts the complete reference damper file as
text, bounded to 1 MiB. Every course in the final imported or tuned scale must
have exactly one `pad` or `free` row. The existing Rust specification owner
checks coverage, finite values, and that each pad lies strictly inside its
speaking string before the engine is created.

```text
frankensim-piano-dampers-v1
pad,69,0.20,0.24,0.7
free,100
```

This example covers a two-key scale. Supply all keys when using the built-in
88-key scale. Pad endpoints are distances from the string's fixed end in metres;
the last field is total viscous drag in N s/m **per speaking string**. A `free`
row declares no damper, independently of MIDI pitch. Pads never touch duplex
segments.

The same admitted geometry and drag govern sounding-string pad projections,
sympathetic-bank decay, damper eligibility, and sustain/sostenuto behavior.
Sustain retains the current squared engagement law. The selected text survives
timing probes and playback copies and composes with `stringStretching`, source
string materials, source gas, and all three mechanical clocks. Omitting it uses
the reference owner's existing estimated pads.

Each station now uses the complete retained string velocity: active partials
and a remaining sympathetic-bank tail share one contact, including the chord's
bridge motion. Its exact rank-one viscous flow mixes those partials and applies
the reciprocal bridge impulse. Modal combinations that leave the physical
station still dissipate no heat. Pre/post station sweeps reverse order around
the mechanical step; shape, ownership and coefficient tables are prepared
before those sweeps, with no audio-step allocation. Exact nonnegative station
heat is booked in `damper_loss_j`, separately from conservative coupling work.

Finite nonnegative pad drag is no longer restricted by the silent oscillator's
underdamped transition domain. Strong pads use the same exact dissipative
flow without clipping their drag.

The model remains finite-footprint viscous drag with the declared quadrature.
Falling damper action and felt impact are separate physical models.
Independent analytic station tests cover nodal
cancellation, reciprocal recoil, integrated viscous heat, and modal ownership;
native engine execution and WASM cost/quality measurements are still required
before claiming new performance or audio-comparison results.
Rebuild WASM with `engine-piano/build-wasm.sh`;
the factory requires `supportsDamperSpecification()` before forwarding this
selection. `info().damperSpecificationSelected` reports an explicitly supplied
file, including an all-free file.

## Reciprocal sympathetic strings

When `Config::silent_bank_hz` is enabled, all retained bank partials enter the
same conservative bridge Schur solve as sounding strings. Each mode uses the
shared modal owner's exact undamped zero-order-hold map; its average string
displacement and the accepted average bridge displacement determine reciprocal
forces within the same step. This resolves coupled normal modes and sympathetic
beating without a delayed bank reaction or a numerical energy source in
`coupling_j`. That account still records the separate nonlinear stretching
split and any extension energy introduced by a changed reduction when selected.

Active prefixes and bank tails share absolute modal displacement and velocity.
Activation and retirement transfer ownership without copying or resetting the
retained motion. Each retained mode appears once in the conservative and
intrinsic-loss solves, and each string's unresolved chord drag is counted once
over the retained union. Existing physical board-basis transfer and failed
rebuild rollback remain in force; finite pads use the same union of velocities.
With string stretching selected, that union also drives the shared nonlinear
state described above. With no stretching, the existing arithmetic is unchanged.

Intrinsic loss uses the same reciprocal midpoint half-flow as sounding strings:
material strain drag retains its bridge cross terms, while gas drag acts on
absolute modal velocity. Critical and overdamped finite nonnegative source
rates are admitted without clipping; nonrepresentable prepared flows refuse.
Their time accuracy, like the existing active loss scheme, depends on the
mechanical timestep. The implementation prepares a board-sized inverse and
per-course projections, then performs allocation-free modal passes and a
board-sized solve. It does not form a dense operator over all partials.

The small production modal adapter has native analytic tests for reciprocal
excitation, energy-conserving resonant beating with second-order frequency
convergence, and the full material/gas/chord resistance and physical heat,
including strong resistance. Public engine regressions cover the same partials
under either owner, ownership changes during ringing, retirement/restrike,
and a strongly damped source-material receiver. Those full engine regressions
and native/WASM cost measurements remain unexecuted in this change; no new
throughput or measured audio-comparison result is claimed.

Additional nonlinear regressions cover a ringing receiver represented entirely
in the bank or with a borrowed active prefix, continued axial motion through
retirement and restrike, and restoration after a refused prefix expansion.
An analytic projection case retains nonzero axial displacement and velocity
while removing a cancelling transverse strain component; it checks the resulting
increase of extension storage and the subsequent free axial decay. These added
tests have syntax validation here; full native integration execution remains
outstanding.

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
The radius is each string's finished radius: its outermost helical wrap, or
else the outer layer of its bonded core. Its line mass remains the material
integral, core plus wraps.
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

The silent bank and sounding strings share the same reciprocal intrinsic-loss
solve. Gas drag acts on absolute modal velocity, while material strain loss
includes the appropriate bridge cross terms. The common bridge recoil and
modal kinetic changes determine physical `string_loss_j`; gas is never
extrapolated into an invented loss on the residual chord. Pad heat comes from
the separate shared station flow, including its reciprocal bridge reaction.

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
curve. With the optional silent bank disabled, unstruck strings retain the
live engine's existing static image;
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
