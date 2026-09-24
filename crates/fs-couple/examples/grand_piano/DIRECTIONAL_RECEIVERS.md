# Directional close microphones from the acoustic field

The exterior piano supports ideal first-order pressure/velocity microphones.
Orientation and proximity response come from the solved three-dimensional
acoustic field, not an equalizer, distance gain or presumed point-source ray.
They work in `response`, `admittance`, `render` and `render-loaded`, including
native board skins, supplied acoustic meshes and static rigid assemblies.

Add these rows to the acoustic specification:

```text
receiver-evaluation,near-field
receiver-m,0.2,1,0.15
receiver-m,1.2,1,0.15
receiver-pattern,0,0.5,0,0,-1
receiver-pattern,1,0.5,0,0,-1
```

The existing receiver positions are in metres in the board frame. A pattern
row is `receiver-pattern,index,pressure_fraction,front_x,front_y,front_z`.
The index is zero-based, following `receiver-m` order. The front axis must be
an explicitly supplied unit vector **pointing from the microphone toward its
front source**. Here both microphones point down at the soundboard. Axes are
not normalized, inferred from geometry, or rotated to make the result work.
Rows may appear before their receiver positions; nonexistent indices,
duplicates, nonunit/nonfinite axes and fractions outside [0,1] refuse.

The pressure fraction selects the ideal first-order family: 1 is omni, 0.5
is cardioid, 0 is figure eight, and intermediate values are also admitted.
Omitting a pattern keeps the original omnidirectional pressure observation.
An all-omni selection, including explicit fraction 1 rows, uses the original
scalar quadrature and preserves its numerical output. A directional selection
requires `receiver-evaluation,near-field`; the centroid observer never stands
in for an unavailable vector field. Existing clearance/work limits remain.

## What is calculated

The core `fs_bem::near_field` owner integrates analytic target-position
spatial derivatives of **both** Green layers, on the same adaptive tree as
pressure. Six derivative kernels have separate discrepancy checks; scalar
pressure convergence alone does not admit particle velocity. Boundary pressure
and normal-velocity traces remain the original BEM solution. For the existing
`exp(-i omega t)` convention, linear Euler momentum gives:

```text
particle_velocity = gradient(pressure) / (i * omega * density)
output = alpha * pressure
       - (1-alpha) * density * sound_speed * dot(front, particle_velocity)
```

The minus sign follows the front-axis convention: a plane wave arriving from
the front propagates opposite that axis. The output is **Pa-equivalent**,
normalized to unit on-axis plane-wave pressure sensitivity. It is not literal
scalar pressure for a directional channel, nor a voltage or calibrated model
of a named microphone. CSV and source reports record the selected fractions,
axes and equivalent-pressure interpretation. The existing `full-scale-pa`
controls PCM encoding with that same normalization; there is no peak or
per-channel normalization.

An outgoing spherical field has particle velocity proportional to
`pressure * (1 + i/(k*r))`. Its reactive term therefore creates proximity
response in the gradient part of the microphone. The program does not assume
that a piano is a monopole: actual panel geometry and interference determine
its pressure and velocity. Coincident opposing cardioids sum to the scalar
pressure; opposing figure-eight axes reverse polarity. These identities are
checked on the real BEM path, separately from plane-wave/monopole coupons.

A prepared observer is reused across every modal solution at a frequency.
Source formulation selection, surface pressure, the complete radiation-load
matrix, materials and mechanical coordinates do not depend on microphone
orientation. In feedback playback the same acoustic storage and string/hammer
trajectory feed all receiver filters on one clock. Orientation changes only
the observation. It never changes power injected into the piano.

## Model D study

`model-d-cardioid-mics.fspe` contains a complete **estimated** setup with two
ideal cardioids 15 cm above the reference board plane. These are not measured
recording-session poses or manufacturer microphone parameters.

```sh
cargo run --release -p fs-couple --example grand_piano -- \
  --preset steinway-d --mesh-divisions 4 --dump-geometry model-d.fsb

cargo run --release -p fs-couple --example piano_exterior -- \
  render-loaded model-d.fsb steinway-d board-skin-continuous \
  crates/fs-couple/examples/grand_piano/model-d-cardioid-mics.fspe cardioid.wav 2 \
  --note 69 --velocity 1 --modes 128 --substeps 8
```

Use fresh output paths and the scale corresponding to the supplied equilibrium.
All existing hammer, contact-face, damper, nonlinear-string, MIDI/CSV and rigid-
assembly controls remain. Geometry, wavelength, passive-load and held-out
receiver-fit failures still refuse; a successful geometry preparation is not a
promise that a full Model D acoustic discretization passes those bounds.

## Accuracy boundary

This is a coincident ideal first-order sensor. It does not model diaphragm or
capsule dimensions, electronics, frequency-dependent polar patterns, low-frequency
rolloff, noise, microphone-body scattering or microphone backreaction. Real
pressure-gradient microphones have design-dependent frequency and proximity
responses; this is not a replacement for a measured microphone transfer.

Both field fitting and any proximity prediction are confined to the admitted
sampled band. In particular, the ideal low-frequency velocity response must
not be read as a validated DC microphone response. Gradient cancellation can
amplify relative observation error. The core reports componentwise quadrature
discrepancy estimates, not bounds on spatial discretization or source-model
error. Mixed omni/directional sets may refine scalar pressure more deeply to
resolve the derivative field. No quadrature, fit or energy threshold is relaxed.

References: DPA Microphones, *Microphone technology - the essentials* and
*Proximity effect in microphones explained*; SCHOEPS, *Proximity Effect*.
These explain pressure/gradient sensing and source-dependent proximity. The
implemented ideal pressure/velocity formula is not a model of their products.
