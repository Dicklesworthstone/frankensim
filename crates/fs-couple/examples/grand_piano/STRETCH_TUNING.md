# Physical octave stretch

Prepare an octave-stretched physical scale, then play it through the existing
piano. The tuner changes string tensions, not oscillator frequencies or output
processing. Length, linear density, bending rigidity, unison spread, and hammer
and duplex geometry remain unchanged.

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
`--preset steinway-d`. The main renderer's flags and defaults are unchanged;
`--stretch-tuning` is not a renderer flag. MIDI, pedals and both existing acoustic
render paths continue consuming the same physical scale format.

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

## Focused checks

```sh
cargo test -p fs-couple --example piano_tune
```

Tests cover the 76 octave pairs, anchored temperament, zero-rigidity limit,
positive-rigidity stretch direction, idempotence, unchanged non-tension fields,
unison tension offsets, unordered/sparse scales, rollback, both built-in scales,
CSV round trips and CLI/output refusal. Native execution is still required in
an environment with the workspace toolchain and dependencies. For listening,
use the same performance, mechanics/acoustic settings and pressure full-scale
in both renders; numerical octave matching is not a perceptual quality claim.
