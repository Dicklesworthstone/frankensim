# Rhythmic physical-player scores

Both `--stick-force-file` and `--second-stick-force-file` accept a bounded
`frankensim-stick-score-v1` file as an alternative to the unchanged two-column
SI force CSV. A score adds reusable hand-force gestures, individual accented
strokes, rolls and tempo changes. It compiles into the **same force program**
before the first mechanical step. It is not a drum sample sequencer, MIDI
importer, prescribed-impact-velocity controller or replacement contact solver.

## One file per physical hand

```text
frankensim-stick-score-v1
# Explicit rhythmic clock: beat index, beats per minute.
tempo,0,240
# Authored push/lift profile: name, offset SECONDS, force NEWTONS.
# These illustrative forces are NOT measured hand or contact data.
shape,tap,0,0
shape,tap,0.003,2
shape,tap,0.006,0
shape,tap,0.009,-1
shape,tap,0.012,0
# start beat, spacing in beats, count, shape, force multiplier
roll,0,1,4,tap,1
# An accent adds physical force; it is not an output-level adjustment.
stroke,2,tap,0.25
```

The example's four gestures begin at 0, 0.25, 0.5 and 0.75 seconds. The extra
stroke adds a quarter-strength gesture at 0.5 seconds. The last tail ends at
0.762 seconds. A left-hand file can use `roll,0.5,1,4,tap,0.8` to interleave
its onsets; its tail ends at 0.887 seconds. Supplied examples use exactly these
right/left schedules and can be passed directly to the existing drum host:

```sh
cargo run --release -p fs-couple --example percussion -- \
  drum-modal-mic 48000 20 --strike-speed-m-s 0 \
  --second-stick-position-m -0.05 0.02 --second-stick-speed-m-s 0 \
  --stick-force-file crates/fs-couple/examples/percussion/alternating-right.score \
  --second-stick-force-file crates/fs-couple/examples/percussion/alternating-left.score \
  > alternating.wav
```

Choose a new output path; shell redirection has its normal overwrite behavior.
The one-second requested pressure render leaves time after both force tails,
not a guarantee that all physical resonance has decayed. The original mechanical
and BEM work limits apply. A successful native render is not claimed by this
example's presence.

Initial stick speeds and strike stations remain independent inputs. Without an
explicit zero first speed, the original 0.8 m/s launch still occurs. The second
hand must have its own admitted physical station. Neither score invents an
additional stick, changes its mass, moves its contact point, clears a head,
reinitializes felt/snare/air history or bypasses a rejected mechanical step.

## Records and timing

`tempo,beat,bpm` records are mandatory, start at beat zero and have strictly
increasing nonnegative finite beat positions and positive finite tempos.
Tempo is piecewise constant. For example, adding `tempo,2,120` to the score
moves beat three to 1.0 seconds; the beat-two onset remains at 0.5 seconds.
Gesture offsets remain **seconds**, even across a tempo change: the musical
clock changes onset spacing, not a supplied physical force profile's duration.

`shape,name,offset_s,force_n` records define piecewise-linear signed forces.
Names use 1..64 ASCII letters, digits, underscores or hyphens. Each shape needs
at least two strictly increasing finite offsets, starts at zero seconds with
zero force, and ends at zero force. Positive pushes toward the head; negative
lifts the existing stick. Unused shapes are validated too. Blank lines and `#`
comments are permitted. All fields and record counts are exact; unknown records
and unknown shape references refuse rather than falling back to a demo.

`stroke,beat,name,scale` begins one gesture at a nonnegative finite beat with a
finite nonnegative force multiplier. Zero is an explicit silent gesture. This
multiplier is not MIDI velocity, a hammer-speed mapping or acoustic gain.
`roll,start_beat,spacing_beats,count,name,scale` repeats that gesture at positive
finite beat spacing; count is a positive integer. Shapes may be declared before
or after strokes. Strokes need not be ordered; equal-beat strokes add rather
than being deduplicated. Tempo records and each shape's own knots must be ordered.

## Overlap, work and refusal

Overlapping gestures on one hand are superposed on the union of their physical
time breakpoints. Quiet gaps stay zero; the last force stays exactly zero.
They are not concatenated, clipped or reduced to the loudest stroke. The original
interval-integrated force path retains knots between mechanical ticks and
accounts for the resulting work through the actual mass-normalized player port.
Both hands share the existing retryable clock and one mechanical acceptance.
No score parsing, expansion or additional allocation occurs in the audio loop.

The complete last gesture tail must fit the requested duration. Source admission
retains the existing 4 MiB limit and caps shape definitions at 128, tempo records
at 1024, source/expanded force knots at 65,536, and expanded strokes at 32,768.
Overlap compilation admits at most 2,000,000 local interpolation terms. These
are setup work limits, not changes to material or contact laws. Nonfinite sums,
representationally collapsed onsets/knots, excess work and the original combined
generalized-force limits refuse instead of silently truncating or clamping.

A gesture onset is **not a promise of a collision at that beat**. The player's
force acts on the current moving stick; contact, rebound, possible missed
strikes and sound follow from mechanics. This fixed-axis open-loop input is not
a complete human arm/grip model or motion-capture controller. It does not establish
force identification, acoustic fidelity, full-band convergence or real-time cost.

Focused regressions use the existing executable target:

```sh
cargo test -p fs-couple --example percussion mechanics::drive::score_tests
```
