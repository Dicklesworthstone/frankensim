# Full physical preparation in finite-body piano playback

Both `piano_exterior render` and `render-loaded` accept the physical controls
already used by `grand_piano`. The former remains one-way; the latter retains
passive radiation feedback. No new piano engine or acoustic solver is selected.
Played output requires one or two receivers. The harmonic `response` and
`admittance` commands support arrays of up to 64; larger played-output requests
refuse before loading the structural inputs or preparing BEM/receiver fits.

```sh
cargo run --release -p fs-couple --example piano_exterior -- \
  render-loaded settled.fss strings.csv body.obj acoustics.fspe loaded.wav 6 \
  performance.mid --modes 128 --substeps 8 \
  --hammers materials.fsh --hammer-footprints faces.fshp --dampers pads.fspd
```

Each control is optional. Old positional MIDI invocations and the default
4-substep, 24-partial source-shank image remain available unchanged. Without a
score the existing A4 study is used, or the nearest admitted key for a partial
scale lacking A4. `--note` and `--velocity` explicitly select a supplied-key
study and post-escapement hammer velocity in m/s. Outputs remain 48 kHz physical-Pa
PCM, without automatic normalization. The render duration remains 0.05–60 s.
If the selected full-scale pressure would clip any PCM sample, the render
refuses before publishing a WAV and reports the clip count and physical peak.
Choose explicit headroom above that peak and rerender; the signal is never
silently attenuated or published with saturated samples.

## Persistent finite-body pressure blocks

The prepared finite-body playback path now exposes
`exterior_audio::ExteriorStream`. Both `render` and `render-loaded` use this
same stream for their existing WAV output. A host can keep it alive across
arbitrary blocks and apply physical key/pedal gestures between calls:

```rust,ignore
// piano, score and baked are fully prepared before constructing the stream.
let mut stream = exterior_audio::ExteriorStream::new(&mut piano, score, &baked)?;
let mut block_pa = vec![0.0; 256 * stream.channels()]; // allocate before playback
stream.render_interleaved_block(&mut block_pa)?;
stream.instrument_mut().set_sustain(0.5)?;
stream.render_interleaved_block(&mut block_pa)?;
```

Blocks contain physical pressure in Pa, in the supplied receiver order. Mono
also accepts `render_block`; stereo uses frame-interleaved L,R output and
refuses implicit downmixing. The mutable piano reference remains with the
stream, and receiver runtimes borrow the immutable baked coefficients. Its
instrument accessors support physical gestures and accounting inspection;
changing the prepared modal basis or clocks requires fresh preparation.

Each output frame dispatches controls on one 48 kHz sample clock, advances the
piano once, differentiates the complete mechanical-substep board trace, and
feeds one shared causal decimator into the fitted receivers. Each receiver
retains its own state-space and propagation-delay history. Block boundaries
reset none of these states. The arithmetic and receiver order match the
previous offline path, including passive radiation feedback when attached to
the instrument. Empty calls consume no events or time, and a live gesture
applies before the next output sample. Scheduled events retain their exact
sample positions even inside a block.

All stream buffers are allocated during construction. Successful block calls
perform no allocation, fitting, BEM solves, file access or logging in this
composition layer. They retain the existing decimator and flight delays;
`sample_position` counts emitted output frames and `decimator_delay_frames`
reports observation latency in addition to the receivers' flight delays.
This is a reusable output API; no audio-device integration or measured
real-time deadline is claimed.

An incomplete stereo frame refuses before advancing controls or mechanics,
zeros the supplied buffer, and may be retried with a complete buffer. A
mechanical, control or receiver execution error preserves only the successful
frame prefix and zeros the entire failed frame and remaining suffix.
`BlockError.completed_frames` counts that prefix in the current block;
`BlockError.sample` gives the absolute failed frame. Such an error latches
the stream: later nonempty calls fail without advancing it. In particular,
an observer can fail after mechanics or an earlier receiver advanced, so
recovery requires a fresh instrument/stream instead of resuming mismatched
histories. The WAV wrapper returns no output from a failed stream or from
clipped PCM conversion.

## Physical preparation options

`--modes 1..512` changes the per-string retention ceiling, including unplayed
sympathetic strings and duplex segments. Bass strings are no longer restricted
to the 24-partial demonstration through this front door. Existing partials retain
their geometry/material-derived frequencies. The original 21.6 kHz retention
ceiling still stops inadmissible modes; raising a count is not a full-band or
convergence certificate. It changes the retained endpoint mass and therefore
the loaded board basis: the acoustic solve is prepared from THAT same new bank.

`--substeps 1..16` selects the mechanical rate, 48 kHz times this count. Both
receivers observe every substep on one performance clock; the existing decimator
is prepared with the same factor. Raising the rate does not silently add string
partials or change the frequency band in the acoustic specification. Every
acoustic pole and coupling-strength guard still applies. In particular, the
8-times-band storage pole requires the top fit frequency below
`2700 * substeps` Hz, in addition to the independent acoustic/output guards.
The 10.8 kHz figure in the original feedback description is the DEFAULT 4x case.

`--hammers materials.fsh` uses the complete per-key WoolFelt/Prony input described
in `HAMMERS.md`. `--hammer-footprints faces.fshp` uses `HAMMER_FOOTPRINTS.md` to
select point or finite longitudinal contact faces. Each site keeps independent
felt/relaxation history, sharing its physical hammer, rather than filtering a
point-contact waveform. `--dampers pads.fspd` uses the supplied spatial viscous
pads from `DAMPERS.md`; `--dampers estimated` explicitly selects the existing
approximate spans and drag values. These are not new measured Steinway data.

With the `steinway-d` scale and source hammer cards, opt-in
`--rt0425-hammer-stiffness` selects the already implemented published per-string
felt stiffness. Adding `--rt0425-hammer-dissipation` replaces the estimated
Prony/crush loss with the published per-note R_H relaxation term. The second
flag requires the first; supplied `--hammers` cards and other scales refuse.
Both `render` and `render-loaded` use the same contact mechanics and retain
their original acoustic distinction. The default hammer law is unchanged.
These source values are not calibration of this particular piano or recording;
compare bridge motion and pressure at held-out keys before accepting them.

`--rt0425-string-damping` independently selects the report's per-key `R_u`
and `eta_u` intrinsic string losses for the `steinway-d` scale. The existing
scalar stiff-string modes use a reduced damping projection; the flag does not
import the report's complete higher-order string model. It reaches both
one-way and radiation-loaded playback, including finite hammer footprints,
and the harmonic `admittance` model. The estimated common loss remains the
default. The CLI requires the source scale for this selection in every command.

All supplied cards must cover EVERY admitted scale key, not just the notes in
the score. Missing files, incomplete cards, duplicate options and invalid
budgets refuse before structural/BEM preparation; no source-default fallback
occurs. Constitutive and spatial admission remain with their existing owners.
The output report identifies the chosen controls, actual mechanical rate,
retained string-coordinate count and contact-site count.

`response` and `admittance` accept `--modes`, `--substeps`, the flat-board inertia
options, `--string-polarization` and `--rt0425-string-damping` after their output
path. Use the same scale, geometry, frame card and retention choices as playback.
Both transverse string directions then contribute their actual endpoint inertia
and reciprocal bridge forces in one loaded board basis. The admittance model
also retains the selected intrinsic loss law for all unison and duplex segments.
Its applied bridge force and reported bridge velocity use the primary hammer
direction; the secondary strings respond through the coupled board.

The complete frame card is admitted before board preparation and projected from
the same retained P1, cubic or crowned motion used by playback. Missing frames
or inconsistent primary geometry refuse. The supplied lateral damper ratio is
still part of the card, but the harmonic model has no key-damper contacts.
`--lossless-structure` remains an admittance-only comparison and refuses an
explicit simultaneous `--rt0425-string-damping` selection.

`response` gives pressure per prescribed modal acceleration. Its input basis
includes the selected directional string mass loading; intrinsic damping does
not change this acoustic motion-to-pressure transfer into a structural force
response. Use `admittance` for damped bridge mobility, receiver pressure per
bridge force, and the wood/string/radiation power balance. Hammer, damper,
score and nonlinear-extension controls still require played simulation.

More retained partials, contact sites and substeps increase work and memory.
Nothing here certifies real-time performance, spatial convergence, calibrated
materials, full keyboard action, or pressure accuracy outside the sampled band.

## Force-driven performances and explicit MIDI mappings

`--performance events.csv` selects the existing `sample,event,key,value` format.
This reaches the source shank's actual jack-force port, not a velocity alias:

```text
sample,event,key,value
0,sustain,0,0.5
24,jack_staccato,69,70
240,sostenuto,0,1
1200,note_off,69,0
1440,sostenuto,0,0
1680,sustain,0,0
```

Sample indices refer to the **48 kHz output clock**, independently of substeps.
The selected mechanical clock resolves the physical force pulse and let-off.
`jack_staccato` is a 7 ms pulse and `jack_legato` a 100 ms pulse, with the CSV
value giving peak newtons at the published jack station. The engine determines
hammer acceleration, shank bending, felt contact and escapement. `note_on`
instead retains the post-escapement velocity interface in m/s. `note_off`,
sustain travel, sostenuto and una-corda use their existing physical owners.
Equal-sample rows retain file order. This is the existing force-driven action
fragment, not a complete keyboard/repetition-action reconstruction.

```sh
cargo run --release -p fs-couple --example piano_exterior -- \
  render-loaded settled.fss strings.csv body.obj acoustics.fspe force.wav 2 \
  --performance events.csv --substeps 8 --modes 128
```

`--midi score.mid` is equivalent to the legacy positional score path. Add
`--midi-channel 1..16`, `--midi-velocity-max-m-s V` and `--midi-half-pedal`
to select the existing importer's mappings. Defaults remain channel 1,
velocity 127 -> 4.5 m/s and switched CC64. The velocity maximum must be in
(0,8] m/s; it is an explicit uncalibrated hammer-speed mapping, not output gain.
Continuous CC64 uses travel `value/127`; CC66/67 retain sostenuto/una-corda.
The report exposes selected note counts, ignored channels/messages, skipped
SysEx, end sample and synthesized end releases. No pitch-wheel string retuning
or unsupported controller emulation is inferred.

CSV and MIDI are alternatives; demonstration overrides cannot accompany either.
MIDI mapping flags without a MIDI score refuse rather than being ignored.
The entire score is admitted before structural/BEM work: missing files, invalid
keys/units/times and events outside the render window never produce a partial
WAV or select a demonstration instead. Gesture admission does not guarantee
that every trial will satisfy physical contact/energy/rate limits; those remain
owned by the transactional mechanical runtime. A very slow hammer launched
below the strings can miss under gravity, and is not forced into contact merely
to produce audible output.
