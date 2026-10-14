# Reciprocal exterior pressure and reaction for a supplied kettle

```sh
# Same head, cavity, mallet, material histories and physical microphone as --wav.
# Last argument is an explicit TOTAL acoustic-preparation work allowance.
(set -C; cargo run --release -p fs-couple --example kettle -- \
  crates/fs-couple/examples/kettle/illustrative-parabolic.fsk 4800 \
  --wav-reciprocal 20 40 500 32 12 0.4 0.3 1.2 2000000000000 > kettle-reciprocal.wav)
```

`--wav-reciprocal` takes the original eight numerical audio arguments followed
by `MAXIMUM_WORK`, an integer in 1..100000000000000. It explicitly selects the
existing `fs_bem::reciprocal_radiation` P1/P0 image. `--wav` remains collocation,
with its original 100-billion-work allowance and original CBIE/Burton-Miller
selection. There is no retry on another image when geometry, power or fitting
fails. This is an execution request, not an observed successful recording.

The outer-wall stations, rigid bearing annulus and original head triangles are
unchanged. One source solve supplies BOTH the finite-point pressure and the
complete radiation-reaction matrix. Exact panel-mean pressures are transposed
against the source rows for surface work. Receiver pressure instead evaluates
the original continuous P1 double layer plus P0 auxiliary single layer; it never
manufactures a collocation solution from those means. The shared passive load
fit, independent receiver audit, physical propagation and single PCM path stay.

Tensioning options in `TENSIONER.md` also work: exact mechanical-sample edge
forces still enter the existing nonlinear head/felt/interior-air/exterior-memory
solve. Selecting an exterior operator does not reset motion, change the force
schedule, make the edge actuator a normal radiation source, normalize audio,
linearize head stretching or remove a felt history. The existing reference and
prepared mechanical implementations remain independent of this cold choice.

The work allowance bounds the sum of conservative per-frequency source work,
including the complete mode batch. It is not measured flops or elapsed time.
The reciprocal owner additionally limits pressure unknowns to 4096, retained
geometry to 256 MiB, surface wavelength resolution and angular quadrature.
Its whole-band slab-based interior-resonance exclusion is still enforced: a
large or deep kettle cannot automatically use the same band as a thin cymbal.
No new resonance-safe BEM formulation, moving boundary or continuum accuracy
certificate is supplied. Surface, pressure-space, spectral and time refinement
are still independent. The sample geometry/material coefficients remain estimates.

Three focused Rust tests exercise compatible controls alongside tensioner
assignments, the actual supplied outer-wall pressure representation and signed
work, plus budget/cancellation/resonance refusal. Existing mechanical and WAV
tests remain registered. Native execution was not available in the authoring
environment; no full-render, listening or real-time qualification is claimed.

```sh
cargo test --release -p fs-couple --example kettle reciprocal_audio_tests -- --test-threads=1
```
