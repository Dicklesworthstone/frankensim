# Reciprocal exterior air and close microphones

`--reciprocal-air ANGULAR_REFINEMENT MAX_SWEEP_WORK` explicitly selects the
existing P1/P0 Galerkin acoustic image for kettle `--wav`. Without it, the
original collocation path is unchanged. Neither formulation is a fallback for
an error in the other; no data, source mode, damping or fit gate is repaired.

```sh
(set -C; cargo run --release -p fs-couple --example kettle -- \
  crates/fs-couple/examples/kettle/illustrative-parabolic.fsk 4800 \
  --wav 20 40 500 32 12 0.4 0.3 1.2 \
  --reciprocal-air 1 2000000000000 > reciprocal-kettle.wav)
```

This is an execution request, not an observed successful native render. The
input remains an illustrative, rigid, sealed bowl, not measured manufacturer
geometry. The existing edge-tensioner and sample-indexed force controls work
unchanged; both formulations feed the same mechanics, felt histories, acoustic
memory, decimation, filters and PCM encoder. A selected microphone does not
become a mechanical force input. Full acoustic radiation still loads the head.

## One source representation for pressure and reaction

The shared `fs_bem::reciprocal_radiation::Prepared::solve_velocity_ports` solves
one UNIT GENERALIZED VELOCITY per complete real surface row. It projects
`Z_ij = integral b_i p_j dA` from those same rows and original P1 pressure means.
Signed, dependent and zero rows are retained. No symmetrization, eigenvalue
clipping, invented resistance or unit-amplitude normalization is performed.
The immutable results keep their original source identities for receivers.

Pressure uses the source's actual P1 double layer/P0 auxiliary single layer.
The auxiliary density is NOT velocity. Panel mean pressures are valid for
work against P0 velocity, but are never substituted into a collocation observer.
Only AFTER observing pressure does the kettle convert velocity response to
acceleration response using `1/(-i*omega)=i/omega`. The impedance remains a
velocity-to-resisting-force map, so this conversion cannot add a frequency
factor to the mechanical load. The same gas density/speed feeds both paths.

Microphones are admitted against actual triangle clearance and winding. An
exterior point within the enclosing sphere is legal. Positive clearance is
still required; no microphone is placed inside the bowl. Adaptive quadrature
is prepared once per frequency and shared across all source modes. A lower
travel delay is peeled from the fitted response and restored once. At distances
below two output samples, all propagation remains in the fitted response and
there is no explicit delay line. There is no additional inverse-distance gain.

## Budgets and limits

Angular refinement is explicit in 0..2. Total source-solve/projection work is
bounded by the supplied integer in 1..10^13, divided across every frequency;
these units are NOT flops, milliseconds or a real-time budget. Geometry uses
at most 4096 acoustic unknowns and 256 MiB logical cache. Near-field quadrature
retains the existing default depth, error and work ceilings, separately from
source assembly work. The reciprocal owner rejects bands outside its slab
resonance-exclusion or angular-resolution chart; broader bands do not silently
switch to another operator. All original fit, mode, pole, timestep, energy and
PCM-clipping limits remain in force.

Discrete reciprocity/passivity does not establish continuum accuracy. Acoustic
mesh, angular quadrature, receiver quadrature, modal basis, nonlinear bandwidth
and finite-band fits still need independent convergence and measurements.
No native runtime, full WAV render, listening or deadline result is claimed.
Eight focused Rust regressions cover original work projection, complex-source
power, signed/zero/dependent bases, cancellation and work refusal, command
admission, actual reciprocal receiver phase/units and close microphone geometry.
Cargo was unavailable during authoring, so these are unexecuted regressions.

```sh
cargo test --release -p fs-bem --lib reciprocal_radiation::ports
cargo test --release -p fs-couple --example kettle audio::reciprocal
```
