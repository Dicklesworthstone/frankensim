# Physical two-head drum inputs

`--drum-spec instrument.fsd` supplies the actual film materials, installed
head tensions, geometry and numerical window. It is not an oscillator-frequency
table, a pitch shifter or an output envelope. Both head pencils are assembled by
`fs-plate::shell::head`, then reduced by the existing `fs-modal` window solver.
The same dimensions feed the enclosed-air volume/modal cavity and the closed
exterior acoustic mesh. The same head modes feed contact, snare coupling, air
work and radiation. No new solver or material law is introduced.

The option is required: drum and snare commands refuse without it, because the
program assumes no built-in head geometry, tension or loss. The tracked
`estimated_drum.fsd` carries the values of the former built-in 14 x 6.5 inch
estimate. A supplied specification strikes at 0.35 of the clear radius unless
`--strike-position-m` is given, where the former no-file input used the node
nearest (0.06, 0). The published dimensional anchors and estimates are
described in `README.md`.

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

Distributed cavity pressure supports both numerical images; geometric
stretching and material history use the joint nonlinear-capable owner. Vented
audio needs the explicit prescribed-vent option described in VENT_RADIATION.md,
and a vent cannot be combined with the closed elastic barrel. A supplied head must physically
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
head,batter,0.000254,4000000000,0.38,1390,3000,stated-unknown
head,resonant,0.0000762,4000000000,0.38,1390,1500,stated-unknown
mesh,5,32
band_hz,80,500
ambient,293.15 101325 0.65
```

`geometry` gives vibrating radius [m], cavity depth [m], and outer shell
radius [m], in that order. The clear radius is not the nominal drum diameter.
The sidewall is rigid by default. `--elastic-barrel wall.fsb` uses these same
dimensions for a finite-thickness shell clamped at both hoops, with supplied
material and damping. See [BARREL.md](BARREL.md) for shared pressure work,
outer-skin radiation, geometry approximation and input limits.
Each `head` record gives its name, thickness [m], Young modulus [Pa], Poisson
ratio, density [kg/m^3], installed isotropic tension [N/m], and where its loss
comes from. The loss field is a token, never a number:
- `stated-unknown`: no citable film or edge loss exists. The head gets zero
  modal drag, and the render prints "head <name> film and edge loss: STATED
  UNKNOWN, entered as 0; its decays are upper bounds".
- `relaxation`: the head's loss is the `--head-relaxation` material law (see
  [HEAD_RELAXATION.md](HEAD_RELAXATION.md)). The head has no modal drag, and
  the file refuses without that option.

An authored damping ratio refuses, because it has no source: a decay fitted to
a recording would otherwise enter as a material property. The batter and
resonant head can have different thicknesses, materials, tensions and losses.

### Material-card heads

A head can obtain its density, Young modulus and Poisson ratio from a compiled
material-card pack instead of the numeric `head` row:

```text
head_material,batter,0.000254,3000,stated-unknown,film.fsmatpk installed-film T=293.15
head_material,resonant,0.0000762,1500,stated-unknown,film.fsmatpk installed-film T=293.15
```

The six comma-separated fields are `head_material`, side, **thickness [m]**,
**installed tension [N/m]**, loss token, and the shared material record
`PACK STATE [AXIS=VALUE ...]`. The last field uses spaces between its tokens.
`STATE` must match the card's explicitly compiled process/state. Query axes use
the pack's declared normalized quantity schema. Pack paths stay relative to
the drum file; the shared native loader admits at most 8 MiB per pack. Each
distinct path is loaded once for the two heads, so both queries use the same
immutable source. A numeric `head` and a `head_material` row cannot both declare
the same side.

The existing material resolver must supply an unambiguous isotropic
`density`, `young_modulus` and `poisson_ratio` at that query. Those exact values
enter the original head mass and bending matrices, geometric stretching,
contact/cavity coupling and sound-rendering path. The complete original pack,
resolved state and source record stay with the executing drum. Its receipt
names the pack identity and state point. Missing properties, wrong state,
out-of-domain queries and later numeric overwrites of resolved coefficients
refuse. No value is inferred from a material name.

This is a frozen material state for each run. Temperature or moisture may
select different source properties, but the declared dimensions and installed
tension remain explicit specimen data; there is no inferred thermal expansion,
retensioning or ongoing heat/moisture evolution. Frequency-, stress-, strain-
or rate-dependent properties require an appropriate constitutive owner and
refuse on this static elastic route. Source-backed elastic coefficients do not
identify damping: `stated-unknown` still supplies zero modal drag. `relaxation`
still requires the existing `--head-relaxation` law with matching equilibrium
Young modulus and a covering frequency band, without adding a second loss term.

`mesh` gives shared radial intervals and azimuthal divisions. `band_hz` gives
the explicit retained eigenfrequency window. `ambient` gives the air inside and
around the drum as one shared gas declaration: the enclosed volume spring (rho
c^2), the distributed cavities, the neck admissions and the exterior radiation
all use it. The six base records are mandatory;
each head appears exactly once, using either head record form. Record order is free, and `#` starts a comment.
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
and rigid hoops. The cylindrical sidewall remains rigid unless an elastic
barrel card is supplied. Individual tuning-lug mechanics, layered/coated-head
identification and measured damping are not created by importing a file.
The gas properties, distributed-air basis limits,
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

## Air

The air is data, never a constant in code. The built-in reference instruments
load the committed record `examples/moist-air-20c-65rh.gas`. An imported description
declares its own `ambient,<gas record>` row, and `--gas-record FILE` overrides
either. The record is the shared declaration (`T_K P_PA RH`, `mixture ...`, or
`thermochemical ...`). A record that does not resolve refuses; there is no air
default. Each run prints which declaration supplied the air.

For `thermochemical` records, supply each normalized species/NASA-9/transport
triple with `--gas-component SPECIES.pack NASA9.pack TRANSPORT.pack` and declare
`--gas-water none|SPECIES_ID --gas-phase frozen|liquid-water` explicitly.
The resolved source record stays with the actual drum; its density and sound
speed supply cavity acoustics and its `rho*c^2` supplies adiabatic bulk
stiffness. Changing a NASA-9 source therefore changes head/air feedback without
retuning either head. Missing or invalid sources refuse. See [AUDIO.md](AUDIO.md#air)
for the shared options and source admission rules.
