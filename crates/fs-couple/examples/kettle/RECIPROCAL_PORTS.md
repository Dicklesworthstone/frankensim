# Shared reaction ports and close microphones

The existing `--wav-reciprocal` selection now uses
`fs_bem::reciprocal_radiation::Prepared::solve_velocity_ports`. No second
reciprocal command or automatic acoustic fallback is introduced. The original
`--wav` collocation image and all edge-tensioner controls remain available.

```sh
# A microphone 10 mm above the reference head, inside the enclosing sphere
# but outside the actual closed exterior. This is a requested offline render,
# not an observed successful recording or a measured-instrument specification.
(set -C; cargo run --release -p fs-couple --example kettle -- \
  crates/fs-couple/examples/kettle/illustrative-parabolic.fsk 4800 \
  --wav-reciprocal 20 40 500 32 12 0.12 0 0.29 2000000000000 > close-kettle.wav)
```

One velocity-source batch supplies both the complete resisting impedance and
all original P1/P0 field solutions. The source owner's transpose projection
uses actual surface areas and unchanged signed normal-velocity rows. No mode is
discarded and no pressure trace is repaired. Receiver evaluation uses those
original fields, never substitutes panel means into a different observer, and
prepares its spatial quadrature once per frequency for all modes.

The reaction matrix is force per generalized velocity. Receiver pressure is
converted from p/v to p/a by multiplication by i/omega under exp(-i*omega*t),
AFTER observation. Only that acceleration response loses a declared propagation
delay for later restoration. There is no second spreading gain or injected
pressure source. Both outputs use the same physical medium and reference mesh.

The shared receiver owner admits actual exterior winding and triangle clearance,
not distance outside a bounding sphere. An interior, intersecting or unresolved
receiver refuses. For a clearance travel time below two output samples, the
whole propagation stays in the fitted transfer and no extra delay is invented.
Other receivers peel their clearance-based delay and restore it exactly once
through the existing fractional-delay implementation. Quadrature tolerances and
fixed-reference geometry remain approximations, not exact clearance certificates.

Existing complete-power admission, passive load fitting, independent pressure
fit audits, source addresses, nonlinear mechanics, tuning forces, decimation and
PCM publication are retained. The declared cold work now includes the complete
source/reaction projection; receiver quadrature retains the existing owner's
separate limits. Angular refinement remains the existing bounded 0..2 selection.
The reciprocal source still refuses bands outside its slab-resonance exclusion
or resolution budget. It never falls back to collocation to make a sound file.

Three additional actual-BEM regressions exercise pressure/load unit separation,
signed sources, near microphones and invalid/budget/cancelled input. The existing
reciprocal image-selection and source-field tests are preserved. Four shared
velocity-port tests exercise complete power/work and cancelled projections.

```sh
cargo test --release -p fs-bem --lib reciprocal_radiation::ports
cargo test --release -p fs-couple --example kettle audio::reciprocal
```

Local Cargo attempts could not start: no Rust toolchain. Native compilation,
full BEM-to-WAV, listening and deadline results are not asserted. The independent
108-case synthetic projection/phase check passes, but is not geometric BEM or
Rust execution. These changes do not resolve the previously observed kettle CI
failure by assertion, provide measured materials or certify full-band sound.
