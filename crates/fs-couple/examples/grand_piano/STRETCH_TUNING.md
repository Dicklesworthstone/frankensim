# Physical octave stretch

Select octave stretch directly when playing the piano or preparing exterior
acoustics. The existing tuner changes string tensions, not oscillator frequencies
or output processing. Length, linear density, bending rigidity, unison spread,
and hammer and duplex geometry remain unchanged.

```sh
cargo run --release -p fs-couple --example grand_piano -- \
  --preset steinway-d --stretch-tuning 440 --render stretched.wav
cargo run --release -p fs-couple --example piano_exterior -- \
  render-loaded board.fsb steinway-d body.obj acoustics.fspe stretched.wav 6 \
  --stretch-tuning 440
```

`--stretch-tuning` requires an explicit A4 reference from 430 to 450 Hz. The
same option is available for exterior `render`, `render-loaded`, `response`
and `admittance`; place it after the command's positional arguments. The board,
body and acoustic files in the example must be supplied physical inputs.
Tuning finishes before control cards, string modes, bridge coupling, soundboard
reduction or exterior loading are prepared. MIDI and CSV performances then play
that same physical scale. Response/admittance CSV comments and playback reports
record the selected policy. A prescribed-motion `response` still reports
pressure per modal acceleration, not force-driven structural response.

Defaults are unchanged: the `grand_piano` preset compensates the rounded source
tensions to first-partial equal temperament at A4=440 Hz, and its explicit
`--scale` preserves supplied tensions. Exterior commands preserve supplied/source
tensions unless stretch is explicitly selected. In `grand_piano`, stretch
excludes both `--concert-pitch` and `--raw-tensions` so an equal-temperament pass
cannot silently discard it. The similarly named `--string-stretching` is a
separate nonlinear axial-extension control, which can compose with this static
tension preparation.

## Exporting a prepared scale

The standalone tuner remains useful for inspection and reuse:

```sh
cargo run --release -p fs-couple --example piano_tune -- \
  --preset steinway-d --concert-pitch 440 --output stretched.csv
cargo run --release -p fs-couple --example grand_piano -- \
  --preset steinway-d --scale stretched.csv --render stretched.wav
```

For an authored scale, replace `--preset steinway-d` on the **tuner** with
`--scale input.csv`. Omitting both uses the explicitly estimated demonstration
scale. `--output` is required and must be fresh. Tuning and validation finish
before output creation; an existing destination is never overwritten.

Do **not** pass `--concert-pitch` to the renderer for this prepared scale: that
existing option would retune its first partials to equal temperament, discarding
the stretch. An explicit scale otherwise preserves its tensions, including with
`--preset steinway-d`. The renderer also accepts `--stretch-tuning` directly;
the exported CSV route remains optional. MIDI, pedals and both acoustic render
paths consume the same physical scale format.

## Declared tuning policy

C4 through B4 have equal-tempered first partials at the requested A4 (430 to
450 Hz, default 440). Working outward, a lower string's second partial matches
its octave neighbour's first partial. For each solved partial n, the tuner uses

```text
k_n = n*pi/L
T = mu*(2*pi*f_target/k_n)^2 - EI*k_n^2
```

The bending term is not scaled with tension. A nonpositive/unrepresentable
tension is an error, not a clamp. The complete solve is transactional: no input
course is changed on failure. Input order is preserved. Sparse scales need an
unbroken chain of supplied octave neighbours to the temperament octave; missing
geometry is never invented. Reapplying the same preparation does not accumulate
stretch. With zero bending rigidity this policy recovers equal temperament.

This is one explicit **2:1 partial-matching policy**, not a measured Railsback
curve, perceptual calibration, or a professional multi-interval tuning fit.
It fits nominal pinned-string partials, not loaded-board poles or finite-rate
dispersion. The Model D table remains an approximate source-derived model;
exporting its tuned tensions does not make it a measured instrument. This does
not address the separate high-frequency soundboard/receiver-attack gap.

## Register-specific multi-partial tuning

`piano_tune` also accepts independent `--bass-octave-weights` and
`--treble-octave-weights`. Each takes four comma-separated, finite, nonnegative
weights for lower:upper partial pairs **2:1, 4:2, 6:3 and 8:4**, respectively.
At least one weight must be positive. Scaling all four weights by the same
positive factor does not change the policy. Each register defaults to
`1,0,0,0`, preserving the original 2:1 arithmetic. Here “bass” means keys below
C4 and “treble” means keys above B4; the temperament octave stays fixed.

For example, the following explicitly blends 4:2 and 6:3 below C4 and keeps
2:1 above B4. These weights are a user-selected objective, **not** a calibrated
Model D recommendation or evidence of a universally preferable tuning.

```sh
cargo run --release -p fs-couple --example piano_tune -- \
  --preset steinway-d --concert-pitch 440 \
  --bass-octave-weights 0,1,1,0 --treble-octave-weights 1,0,0,0 \
  --tuning-band 20000 --output voiced-scale.csv
cargo run --release -p fs-couple --example grand_piano -- \
  --preset steinway-d --scale voiced-scale.csv --render voiced-piano.wav
```

Do **not** add `--concert-pitch` **or** `--stretch-tuning` to the renderer when
using this exported scale: either would replace the prepared tensions. The
new weight/band flags belong to `piano_tune`; the renderer continues to consume
its ordinary scale CSV. The same import route applies to exterior commands.
All key, pedal, hammer, board and acoustic behavior remains in the existing
renderer rather than a separate instrument or output filter.

### Objective and physical solve

For each octave the already-tuned neighbour is held fixed. If `n` is the upper
partial number, define `beat_n = f_lower(2*n) - f_upper(n)`. Mixed weights
minimize the fundamental-equivalent mean squared mismatch

```text
J(T) = sum_n w_n * (beat_n(T)/n)^2 / sum_n w_n
```

Only one course's positive **tension** is solved at each outward step. The
frequency owner still computes every partial from that tension, length, mass
and bending rigidity. Neither an assigned oscillator frequency nor a fixed
key-number correction is introduced. A one-pair policy uses the exact existing
partial solve; mixed policies use the monotone derivative of this convex
objective, with bounded bisection and a final stationarity check through the
physical frequency owner. They preserve all non-tension fields, input order,
the anchored octave, transactional failure and repeatable preparation.

The optional `--tuning-band` is a finite positive ceiling in Hz on the
**neighbour's target partial**. Pairs above it are excluded and the remaining
weights renormalized. An octave with no retained positive-weight target is an
error, not a silent fallback. This is not a claim that the final response is
bandlimited, a renderer sample-rate setting or a soundboard validity bound.
Omitting the flag admits all explicitly weighted targets.

The API `stretch_tuning::partials::tune` returns every solved octave's signed
pair mismatches, with excluded/zero-weight pairs marked `None`, and
`sqrt(J)` in Hz. The CLI prints the chosen weights, target band, pair count and
largest reported weighted RMS mismatch. These are frequency comparisons,
**not** perceptual roughness or loudness measurements. The CSV contains the
complete resulting physical scale, not just the report or a set of MIDI notes.

This is a sequence of local adjacent-octave fits, not a simultaneous global
keyboard optimization, fifth/twelfth temperament fit, measured tuning curve,
loaded-board calibration or listening-quality guarantee. Arbitrary weights
can widen some intervals or trade one partial's beating for another's. Compare
explicit policies on the same performance and mechanics before selecting one.

### Multi-partial checks

Nine additional core regressions cover all four individual pair policies,
mixed-objective stationarity and candidate comparisons, exact default parity,
physical invariants, order/idempotence, flexible strings, band exclusions and
transactional refusal. Four CLI regressions cover actual scale preparation,
CSV round trips, option admission and exclusive output after successful tuning.

An independent Python/SciPy reference checked 3,648 octave solves across 48
full-scale configurations: both supplied scales, four A4 references and six
weight policies. Mixed tensions agree with an independent unscaled Brent
solution and improve the declared local objective over the individual-pair
candidates. This is a separate numerical reference, **not execution of the
Rust tests** or a listening result. For one source-scale C1:C2 fixture with
its upper neighbour held fixed, equal 4:2/6:3 weights reduce that specific RMS
metric from approximately 0.04370 to 0.01809 Hz relative to a 2:1-only solution.
The comparison does not establish a corresponding percentage improvement in
sound quality or accuracy against a real piano.

## Focused checks

```sh
cargo test -p fs-couple --example piano_tune
cargo test -p fs-couple --example grand_piano stretch_selection
cargo test -p fs-couple --example piano_exterior stretch_tuning
```

Tests cover the 76 octave pairs, anchored temperament, zero-rigidity limit,
positive-rigidity stretch direction, idempotence, unchanged non-tension fields,
unison tension offsets, unordered/sparse scales, rollback, both built-in scales,
CSV round trips and CLI/output refusal. Frontend tests cover option conflicts,
unchanged default tensions, source/imported scale parity, tuning before control
preparation, and actual prepared Bank partials across bass, anchor and treble
courses including unison tensions. Full frontend execution is still required in
an environment with the workspace toolchain and dependencies. For listening,
use the same performance, mechanics/acoustic settings and pressure full-scale
in both renders; numerical octave matching is not a perceptual quality claim.
