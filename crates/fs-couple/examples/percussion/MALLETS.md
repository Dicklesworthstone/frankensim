# Finite-area felt mallet strikes

`--mallet-spec tip.fsmallet` replaces the first hard stick and its point Hertz
contact with one supplied effective inertia and a circular felt face.
`--second-mallet-spec tip.fsmallet` independently replaces the second stick.
Either or both may be selected; an unselected striker keeps its original law.
These are replacements, never additional felt and Hertz reactions in parallel.

The full face samples actual batter-head displacement at four positive-area
quadrature sites. Each site keeps its own compression, irreversible conditioning
and optional Kelvin recovery. All sites share one mallet inertia and react on
the same head. There is no force pulse inferred from impact velocity, sound
sample, output envelope or microphone gain adjustment.

```sh
# Existing drum mechanics and supplied felt excitation; explicit position.
cargo run --release -p fs-couple --example percussion -- \
  drum 4096 --strike-position-m 0.06 0.01 --strike-speed-m-s 0.5 \
  --mallet-spec crates/fs-couple/examples/percussion/estimated-felt-mallet.fsmallet \
  --analytic-newton --impact-substeps 8 511 > felt-strike.csv

# Two different physical tip cards can be supplied here. This example uses the
# same illustrative card, independently, at two different positions and speeds.
cargo run --release -p fs-couple --example percussion -- \
  snare-mic 4800 20 0.08 0.05 0.35 \
  --strike-position-m 0.06 0.01 --strike-speed-m-s 0.5 \
  --mallet-spec crates/fs-couple/examples/percussion/estimated-felt-mallet.fsmallet \
  --second-stick-position-m -0.05 0.02 --second-stick-speed-m-s 0.3 \
  --second-mallet-spec crates/fs-couple/examples/percussion/estimated-felt-mallet.fsmallet \
  --head-stretching --analytic-newton --impact-substeps 8 511 --cavity-modes \
  --microphone-right -0.12,0.05,0.4 > felt-snare.wav
```

These are supported entry points, not claims of completed, calibrated or
real-time renders. The example card is authored, not identified wool-felt data.
The existing force, energy, densification, slope, step and acoustic-fit limits
can still reject a performance. No guard is relaxed to guarantee a recording.

## Complete physical input

UTF-8 input is bounded to 64 KiB and starts with `frankensim-felt-mallet-v1`.
Blank lines and `#` comments are allowed. Exactly one of each mandatory record
is required; zero to four `creep` records may be present:

```text
frankensim-felt-mallet-v1
geometry,0.02,0.012,0.006,0.00002
felt,100000,0.2,2.2,3,0.15,0.7
conditioning,0
creep,3000,6
```

`geometry,effective_mass_kg,face_radius_m,felt_thickness_m,initial_gap_m` uses
positive mass/radius/thickness and nonnegative gap. The effective mass includes
whatever head/shaft/grip reduction the supplied card represents; it is not added
to the default Z5A mass. No shaft geometry or calibrated mallet label is inferred.
The face is flat, circular, parallel to the undeformed batter, and constrained
to the same vertical fixed axis as the existing stick. Initial tip displacement
is minus the supplied gap, preserving the tip observation convention.

`felt,reference_stress_Pa,reference_strain,loading_exponent,unloading_exponent,
crush_fraction,densification_strain` supplies the existing `WoolFelt` law. The
reference value is a stress, not the total force of the whole face. Every site
receives its physical area times that stress, and its existing material validity
gates apply. `conditioning,prior_maximum_strain` explicitly supplies prior crush
history. There is no implicit pristine reset between strokes.

`creep,stiffness_N_m,viscosity_Ns_m` gives one series Kelvin element for the
WHOLE face. Site values are scaled by their area fractions using the existing
`MovingPads` implementation; splitting the footprint does not multiply total
creep stiffness or viscosity. Omitting creep explicitly selects no recovery
branches. A selected branch must have finite positive coefficients.

The normal playing controls remain separate: `--strike-speed-m-s` and
`--second-stick-speed-m-s` supply actual initial velocities. The first mallet
requires `--strike-position-m X Y`; the second requires an explicit second-stick
position as before. Existing independent stick-force CSV programs push and lift
the corresponding mallet through its supplied inverse-root-mass force port.
A force timestamp does not prescribe a contact or rebound timestamp.

## One shared mechanical and acoustic path

`FeltStriker` reuses `MovingPads` by a cold coordinate permutation and origin
translation. It adds no contact law, time integrator or active force. Four points
at radius `face_radius/sqrt(2)` with equal areas `pi*face_radius^2/4` integrate
constant, linear, quadratic and cubic polynomial fields over the disk exactly.
This does NOT make nonlinear contact stress or a finite-element mode field
exactly integrated. Independent point histories are evaluated before summing
forces; replacing them by an averaged shape would miss local contact when
opposite sides of a mode move in opposite directions.

The entire disk must fit inside every edge of the actual polygonal head, not
just inside its circumcircle or at the four quadrature nodes. Rim crossings
refuse rather than snapping, shrinking or clipping the face. This is a bounded
four-site reference, not a resolved growing Hertz patch or a spatial-convergence
certificate. Each mallet consumes four of the existing sixteen felt-pad slots.
Other pads count against that same total. No pad, wire or mode is discarded.

The tip retains the original first/second striker coordinate. Both head ranges,
all snare wires, moving carrier and acoustic coordinates stay unchanged. The
ordinary full mechanical energy and loss accounting include felt crush and
creep once. `felt_crush_j` and `loss_j` in mechanics CSV expose those totals.
Neither mallet is a direct radiation source; the existing observer receives the
resulting head motion, and the same single trajectory drives mono or stereo.

Available on `drum`, `drum-stretch`, `snare`, `snare-off` and their `-wav`/`-mic`
forms. A snare automatically selects its nonlinear-capable joint solver when
felt is supplied, independently of head/wire stretching. Preparation, analytic
Newton and bounded substeps preserve the same physical model. Supplied head and
wire specifications, head relaxation, carrier forces, fixed mufflers, and the
already-admitted cavity and prescribed-vent options keep their semantics.
The `drum-modal` linear-only image refuses rather than dropping felt history.
Curved cymbal faces, rimshots, mallet rotation, shaft bending, stick-stick
collisions, evolving footprint area, and exterior mallet radiation are outside
this flat fixed-axis chart. Existing snare restrictions on compliant mutes remain.

## Focused tests

Four core tests exercise coordinate/gap/mass scaling, distinct local contact,
actual impact with retained history and exact retry, and malformed admission.
Five example tests exercise the complete file, physical disk moments and rim
clearance, full two-mallet/snare/acoustic addresses, two force programs with
felt/head/air dynamics and retry, and unchanged unselected behavior.

```sh
cargo test --release -p fs-couple --lib render::plate::impact::compliant::striker::tests
cargo test --release -p fs-couple --example percussion mallets::tests -- --test-threads=1
```

Native tests and end-to-end audio require execution. Independent arithmetic and
source checks are not evidence of a Rust test pass, audio fidelity or calibration.
