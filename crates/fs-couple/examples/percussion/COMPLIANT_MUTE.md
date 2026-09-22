# Moving compliant mutes: contact, inertia, squeeze and retraction

`--compliant-mute file.fsm` adds one felt-lined moving pad on a drumhead or
one/two opposed felt-lined jaws on a cymbal. Physical signed force programs
move real jaw inertias. Each footprint site makes compression-only contact
with the actual vibrating surface. There is no PCM fade, mode reset, prescribed
contact-force pulse, or instantaneous release command. This differs from
`--muffler`, whose fixed bilateral viscous attachment remains available.

```sh
# Opposed jaws on one nonlinear cymbal; retain the original stick and stand felts.
cargo run --release -p fs-couple --example percussion -- \
  splash 4096 --analytic-newton --impact-substeps 8 511 \
  --compliant-mute crates/fs-couple/examples/percussion/estimated-cymbal-mute.fsm \
  > squeezed-cymbal.csv

# The same physical contact before both BEM microphone observers, not stereo panning.
cargo run --release -p fs-couple --example percussion -- \
  splash-mic 4800 20 --analytic-newton --impact-substeps 8 511 \
  --microphone-right -0.08,0.05,0.35 \
  --compliant-mute crates/fs-couple/examples/percussion/estimated-cymbal-mute.fsm \
  > squeezed-cymbal.wav

# Exterior moving pad on a stretching head, retaining distributed cavity loss.
cargo run --release -p fs-couple --example percussion -- \
  drum-stretch 4096 --analytic-newton --impact-substeps 8 511 \
  --cavity-modes --cavity-drag-per-s 20 \
  --compliant-mute crates/fs-couple/examples/percussion/estimated-drum-mute.fsm \
  > played-mute.csv
```

The supplied geometry, masses, gaps, pad laws and playing forces are **editable
estimates**, not measured fingers, skin/tissue parameters, commercial gel pads,
or proof of an acoustic match. These commands demonstrate control syntax;
convergence, a particular choke strength or release time, and real-time speed
are not guaranteed by an input file. The sample force programs end at 0.008 s;
the requested mechanical/audio duration must cover every force knot.

## Explicit physical input

The bounded UTF-8 format starts with `frankensim-compliant-mute-v1`, supports
`#` comments, and contains comma-separated SI records:

```text
surface,shell
site,x_m,y_m,area_m2
jaw,above,mass_kg,drag_Ns_m,gap_m,thickness_m,f_ref_Pa,eps_ref,p,q,crush_fraction,eps_densify,prior_strain
creep,above,stiffness_N_m,viscosity_Ns_m
force,above,time_s,force_N
```

`surface` occurs exactly once: `shell`, `batter`, or `resonant`. One to four
positive-area `site` rows are footprint quadrature points in the reference XY
plane, not a boundary polygon or four copies of the same averaged modal shape.
Their areas partition the actual pad area on **each** jaw. Every site samples
its own unnormalized physical surface participation. Positions outside the
surface or inside the mounting hole refuse; nonmoving head positions refuse.
There is no automatic spatial-resolution or contact-patch convergence claim.

One or two `jaw` rows provide independently translating effective masses and
pads. A shell can have `above`, `below`, or both. A drumhead permits exactly one
exterior jaw: `batter` with `above`, or `resonant` with `below`. Interior jaw/air
displacement is not modelled and cannot be silently omitted. The example's
surface-motion axis is downward; jaw displacement, velocity and player force
are **positive inward on either side**. The pad face conforms to the declared
reference surface, with the same initial gap at every site. A jaw starts at
rest with zero inward travel, not in a manufactured equilibrium.

`drag_Ns_m` is nonnegative grounded drag on the jaw, not damping on the sound
output. `gap_m` is nonnegative. The remaining parameters bind directly to the
existing `fs-material::fiber::WoolFelt` law and conditioning history. Zero to
four `creep` rows per jaw specify series Kelvin elements for that jaw's **whole
footprint**. The compiler partitions both stiffness and viscosity by area;
each site's deformation and irreversible history remain independent.

Every jaw needs at least two `force` rows. Times must be finite, nonnegative and
strictly increasing for that jaw; the first/last forces must be zero. Values
are external inward forces in newtons, with zero outside the supplied interval.
The existing force player integrates linear interpolation over each mechanical
tick. Positive force pushes inward, negative force retracts. **Zero force does
not open a gap or stop an inertia.** Contact ends only when motion and pad
recovery actually remove compression. Neither force knots nor a render horizon
teleport the jaw or reset material history. Absent-jaw force/creep records,
unknown fields, duplicates, invalid cards and files above 64 KiB refuse.

## Coupled playback and observations

All sites on one jaw share exactly one inertia. Equal-and-opposite contact
reactions enter the existing nonlinear impact solve with the shell/heads,
sticks, stand, and enclosed air. Stored pad recovery energy, jaw kinetic energy,
viscous work and felt conditioning loss participate in its energy balance.
The original resonator addresses do not move. Existing contacts, mufflers,
stand pads, swept-volume rows and radiating modes have exact zeros on the
new jaw coordinates. Gas inertia is appended before private Kelvin state.

The same one-clock force staging accepts up to four distinct physical inputs:
two sticks plus two jaws. Existing independent stick force files, supplied
shell/drum cards, spatial mufflers, analytic Newton, internal substeps and
mono/stereo observers compose. A refused mechanical tick consumes none of
these force programs; internal recovery rolls back all jaw/site histories too.
The no-mute path and all original physical/mode/energy ceilings remain.

CSV appends `mute_above_*` and/or `mute_below_*`: inward travel in metres,
inward velocity in m/s and the sum of current site contact forces in newtons.
The latter is an **endpoint constitutive observation**, not the discrete
step-average reaction. `player_work_j` is the current tick's total external
work over all sticks/jaws, not just the stick work or a cumulative sum.

This image is nonlinear even for linearly reduced drumheads. `drum-modal` and
snare commands refuse it instead of converting their different, larger contact
solver or truncating wires. Vented exterior audio remains unsupported. Jaw
radiation and acoustic occlusion/scattering by the mute are omitted: existing
one-way BEM observes the resulting instrument motion on the same boundary.
There is no friction, hand skeleton, rotation, adhesive force, manufacturing
calibration, or full-band cymbal/skin fidelity claim. Those require separate
physical models and measurements, not a more elaborate output envelope.

Focused native checks:

```sh
cargo test --release -p fs-couple --lib render::plate::impact::compliant
cargo test --release -p fs-couple --example percussion compliant_mute -- --test-threads=1
```
