# A single head on a supplied kettle profile

```sh
cargo run --release -p fs-couple --example kettle -- \
  crates/fs-couple/examples/kettle/illustrative-parabolic.fsk 4096 > kettle-motion.csv
cargo test --release -p fs-couple --example kettle
```

This example connects the previously unconsumed `MeridianCavity` to a real
`TensionedDisk` pencil, its statically relaxed nonlinear stretching reduction,
and a four-site `FeltStriker`. It uses the existing analytic prepared Gonzalez
step, not an instrument-local integrator, authored frequency bank or strike pulse.
The constant pressure mode supplies the bulk compliance exactly once; nonuniform
modes feed pressure back onto the same vibrating head. The mallet only touches
the head, never drives gas coordinates directly. Rejected steps do not advance
motion, felt conditioning, recovery or time.

The file supplies the actual piecewise-linear inner wall. Its first station is
on the axis; the last is the circular head rim. Reentrant radial profiles are
allowed when they remain positive radial graphs. The bowl wall and rim are rigid,
sealed and stationary. There is no fictitious second membrane closing its base.
This is a useful kettle/head/air image, not a calibrated complete timpani.

The illustrative file is **not measured manufacturer geometry**. Its 26-inch
clear diameter and 280 mm depth are chosen inputs; its wall samples a parabola.
The PET, installed tension, damping and felt data are editable estimates. A
parabolic/quartic kettledrum wall is supported as a geometry family by US4674390:
https://patents.google.com/patent/US4674390 . No production dimensions, copper
thickness, residual stresses or manufacturing map are inferred from that patent.
Rigid walls intentionally have no invented copper modulus or loss coefficient.

Distributed air loading matters: Davis's *Mathematical modeling of the orchestral
timpani* treats membrane/air coupling rather than obtaining timpani pitches by
retuning an isolated membrane: https://docs.lib.purdue.edu/dissertations/AAI8825522/ .
This example still omits **exterior air reaction**, wall elasticity, the vent,
pedal/rim mechanics, evolving mallet footprint and material identification.
Its internal-pressure columns are not microphone signals or substitute WAV audio.

## Inputs and approximation boundaries

All record values use the SI units annotated in the supplied file. Every record
except `station` and `creep` is mandatory exactly once. `station` repeats in
strictly increasing axial order, with zero radius only at the bottom. Zero
`creep` rows means no Kelvin recovery. `maximum_slope` bounds the existing
moderate-slope membrane model; it is not a pitch-bend parameter. The head is
fixed transversely at the rim, with free bending slope and fixed incremental
in-plane rim motion. Installed tension is included once in its original pencil.

Head mesh/window and cavity meridional/angular mesh/window are independent.
Cold construction refuses more than 32 head modes, 8 complete cavity modes,
512 cavity nodes or the stated shared-owner work limits; nothing is silently
truncated. Both members of retained angular pressure pairs survive. The modal
frequencies printed on stderr are **uncoupled basis frequencies**, not predicted
coupled instrument pitches. Head angular discretization and the axisymmetric
pressure chart approach the same circular geometry under independent refinement.

The mallet is one fixed-axis effective inertia. Its entire flat circular face
must fit inside the actual polygonal head; four positive-area sites each retain
separate material history. Whole-face Kelvin coefficients are partitioned by the
existing pad owner. Four-site quadrature is not an evolving-contact-area model.
The mallet begins at minus its declared gap with its declared velocity; no
parallel Hertz tip, normalizing gain or force envelope is installed.

The mechanical CSV reports physical mallet/head displacement, actual internal
pressure, and the shared energy ledger. The owner enforces a 20 J energy limit,
1e6 generalized-force limit, 1e-9 J absolute / 1e-6 relative balance tolerance,
and the requested accepted-step budget. Clock and spectral caps are admission
limits, not alias-free or real-time certificates. Native compilation, full
renders, listening tests and deadline benchmarks were unavailable during this
implementation. The Rust regressions are authored tests, not observed passes.
