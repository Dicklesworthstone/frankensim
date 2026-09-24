# Physical two-head drum inputs

`--drum-spec instrument.fsd` supplies the actual film materials, installed
head tensions, geometry and numerical window. It is not an oscillator-frequency
table, a pitch shifter or an output envelope. Both head pencils are assembled by
`fs-plate::shell::head`, then reduced by the existing `fs-modal` window solver.
The same dimensions feed the enclosed-air volume/modal cavity and the closed
exterior acoustic mesh. The same head modes feed contact, snare coupling, air
work and radiation. No new solver or material law is introduced.

Without the option the original estimated 14 x 6.5 inch reference remains the
input, including its original damping and nearest-node default strike. The
published dimensional anchors and estimates are described in `README.md`.

## Use

The tracked `estimated_drum.fsd` is an editable starting point, not a measured
instrument. For a nonlinear head/cavity attack study:

```sh
cargo run --release -p fs-couple --example percussion -- \
  drum-stretch 4096 \
  --drum-spec crates/fs-couple/examples/percussion/estimated_drum.fsd \
  --cavity-modes --prepared-nonlinear --strike-speed-m-s 4 \
  --strike-position-m 0.06 0.01 > drum.csv
```

The same file works with `drum`, `drum-modal`, `snare`, `snare-off`, and their
existing `-wav` / `-mic` variants. For example, a fixed-receiver pressure render:

```sh
cargo run --release -p fs-couple --example percussion -- \
  drum-modal-mic 4800 2 0.08 0.05 0.35 \
  --drum-spec crates/fs-couple/examples/percussion/estimated_drum.fsd \
  --strike-position-m 0.06 0.01 > drum.wav
```

Existing restrictions remain: distributed cavity and stretching need their
nonlinear image, not the modal/snare image. Vented exterior radiation is still
unimplemented, so a cavity neck remains CSV-only. A supplied head must physically
contain the unchanged snare span; this option never silently shortens wires,
retensions them or substitutes a different instrument. `--shell-profile` remains
the separate cymbal geometry input and cannot be mixed into a drum command.

With a supplied drum and no explicit strike position, the contact is interpolated
at `(0.35 * clear_radius, 0)` metres, rather than snapping the old fixed 6 cm point
onto an unrelated rim. An explicit `--strike-position-m X Y` uses the existing
geometric point locator and refuses points outside the mesh. Use the same explicit
position when comparing imported and original reference trajectories.

## SI file format

```text
frankensim-drum-spec-v1
geometry,0.1703,0.1651,0.1778
head,batter,0.000254,4000000000,0.38,1390,3000,0.001
head,resonant,0.0000762,4000000000,0.38,1390,1500,0.001
mesh,5,32
band_hz,80,500
```

`geometry` gives vibrating radius [m], cavity depth [m], and rigid outer shell
radius [m], in that order. The clear radius is not the nominal drum diameter.
Each `head` record gives its name, thickness [m], Young modulus [Pa], Poisson
ratio, density [kg/m^3], installed isotropic tension [N/m], and nonnegative modal
damping ratio. The drag coefficient is `2 * ratio * angular_frequency`; zero
means zero declared head drag, not fallback damping. The batter and resonant
head can have different thicknesses, materials, tensions and losses.

`mesh` gives shared radial intervals and azimuthal divisions. `band_hz` gives
the explicit retained eigenfrequency window. The five base records are mandatory;
each head appears exactly once. Record order is free, and `#` starts a comment.
Unknown records, duplicates, missing values, nonfinite parameters and invalid
material domains refuse. Input reads are bounded to 64 KiB. Constitutive
validation remains with the existing section owner.

The bounded example admits 1..32 radial intervals and 8..128 azimuths. Existing
mesh, nonlinear reduction, acoustic panel and total-state budgets still apply;
not every combination within these scalar bounds fits every physical image.
An empty window or more than 63 combined head modes refuses rather than silently
changing the requested window. Distributed air needs additional state capacity.
The window must also satisfy the mechanical Nyquist guard. Audio admission
requires it inside the existing 40..1640 Hz acoustic bake band; broader windows
are CSV-only. These checks do not establish modal convergence or full-band sound.

## Scope and verification

This remains a homogeneous isotropic two-film material model with uniform
installed tension by default, optional equilibrated spatial tensor variation,
and a rigid cylindrical shell/rim. Elastic shell motion, individual tuning-lug
mechanics, layered/coated-head identification and measured damping are not
created by importing a file. The gas properties, distributed-air basis limits,
striker, tip contact and snare properties remain the existing declared inputs.
Exterior radiation is the existing one-way, undeformed-boundary calculation;
it does not add acoustic backreaction or room scattering.

Nine focused tests cover strict import and physical-field identity, command
composition/refusals, original-reference trajectory equivalence with an explicit
strike, actual FEM stiffness/resonance scaling, depth-dependent air feedback,
matching larger-head exterior geometry, and supplied-geometry composition with
stretching, distributed air and prepared nonlinear mechanics. They exercise
production owners, not substitute oscillators. Native execution is required:

```sh
cargo test --release -p fs-couple --example percussion -- --test-threads=1
```

A successful import is neither a measured-specimen certificate nor a real-time
performance claim. It provides the physical parameter path needed for those
subsequent comparisons without editing the solver or the example's Rust source.

## Spatial installed tension

Optional per-head `tension_variation` records now supply equilibrated affine
tensor prestress to the original head pencil before eigenanalysis. See
[TENSION.md](TENSION.md) for the SI law, complete example, tensile-domain checks
and composition with contact, cavity, nonlinear heads and pressure rendering.
Omitted or all-zero variations preserve the existing uniform-tension image.
