# Profile-defined kettle: head, felt, interior air and exterior pressure

```sh
# COUNT means accepted mechanical steps without --wav.
cargo run --release -p fs-couple --example kettle -- \
  crates/fs-couple/examples/kettle/illustrative-parabolic.fsk 4096 > kettle-motion.csv

# With --wav COUNT means 48 kHz output frames. Full scale is 20 Pa;
# acoustic band 40..500 Hz, 32 training intervals, fixed fit order 12;
# receiver x/y/z are metres in the supplied wall's coordinate frame.
(set -C; cargo run --release -p fs-couple --example kettle -- \
  crates/fs-couple/examples/kettle/illustrative-parabolic.fsk 4800 \
  --wav 20 40 500 32 12 0.4 0.3 1.2 > kettle.wav)

cargo test --release -p fs-couple --example kettle -- --test-threads=1
```

These are execution requests, not claims of calibrated renders. The audio
command's native run is recorded below; listening tests and deadline benchmarks
have not been performed. The existing percussion CI lane includes this example.

## One physical composition, existing numerical owners

`TensionedDisk` supplies the actual film pencil. Its statically relaxed nonlinear
stretching reduction retains installed tension exactly once. `MeridianCavity`
and its shared `end_disk_interface` supply the actual wall-profile pressure
basis, volume norms and signed overlaps. `CavityCoupling` applies reciprocal
pressure work: the constant mode supplies bulk compliance once; nonuniform modes
retain inertia. No equivalent-volume cylinder or fictitious second head is used.

One `FeltStriker` replaces a hard tip, not an additional Hertz reaction. The
entire flat circular face must fit inside the actual polygonal head. Four
positive-area sites retain independent wool-felt conditioning and Kelvin
recovery while sharing one effective mallet mass. Whole-face Kelvin coefficients
are partitioned by the existing pad owner. The initial gap and velocity are
physical inputs, not an authored force envelope or calibrated tip profile.

The shared analytic prepared Gonzalez step owns all mechanics, material history,
interior pressure and selected exterior memory. Failed steps retain accepted
state/history/time. The mechanical CSV reports actual motion, internal pressure
and energy/loss accounting. Its pressure columns are NOT microphone signals.

## Explicit geometry and material provenance

Every SI record is annotated in `illustrative-parabolic.fsk`. Records other than
repeated `station`, optional repeated `outer_station`, and optional `creep` are
mandatory exactly once. Input is bounded to 64 KiB. Head thickness, E, nu, density,
installed tension, damping, slope cap, air state, meshes, spectral windows,
clock, mallet inertia/face and felt coefficients remain explicit editable inputs.

The example is **not measured manufacturer geometry**. Its 26-inch clear
head diameter and 280 mm depth are chosen inputs; straight wall segments sample
a parabola. Parabolic/quartic kettledrum walls are a documented geometry family:
https://patents.google.com/patent/US4674390 . No production dimensions, hammer
map, copper modulus, residual stress or normal thickness are inferred from that
patent. PET/tension/damping/felt values and the separate outer profile are also
illustrative, not identified material data. The bowl/rim remain rigid and sealed.

Membrane/air coupling is central to orchestral timpani models, rather than
retuning isolated membrane frequencies into a harmonic series. See Davis,
*Mathematical modeling of the orchestral timpani*:
https://docs.lib.purdue.edu/dissertations/AAI8825522/ . The uncoupled frequencies
printed by the mechanics command are basis frequencies, NOT loaded instrument
pitches. This image omits wall elasticity, vent flow, pedal/rim movement, evolving
mallet footprint, skin identification, stand motion and calibration.

## Closed exterior and coupled 48 kHz pressure

Audio requires a separate supplied outer wall: its axis bottom lies below the
inner bottom, its top matches the head plane and its rim lies beyond the film.
Enclosure is checked at the union of inner/outer profile knots; wall crossings
between inner stations cannot be hidden. A rigid annulus joins the actual head
seam to the outer rim. Only the head moves: mallet and interior air coordinates
never become exterior radiation sources.

One closed-boundary triangle BEM sweep supplies BOTH the finite-point receiver
pressure and the complete resisting force/velocity matrix. Transposed source
rows project those same pressure traces into signed cross-mode work.
Reciprocity makes the exact modal impedance symmetric (Z = Z^T for real mode
shapes); collocation breaks it by its discretisation error, about 0.5 % on the
illustrative card, and the skew part of the large reactance then shows up as
+/- eigenvalue pairs of the Hermitian form. The load and the power form
therefore use the reciprocal part Z_r = (Z+Z^T)/2, whose whole power form
Re Z_r is checked. A sweep refuses when the asymmetry exceeds 0.02, the load
fit's own complex RMS limit, or when a frequency's lowest power eigenvalue lies
below -asymmetry*|Z| at that frequency; nothing is clipped. The receipt states
the largest asymmetry |Z-Z^T|/|Z|, the band's lowest power eigenvalue, and the
least headroom above that refusal (the two need not share a frequency).
`--reciprocal-air` is symmetric by construction, but it still fits the fixed
32-term dictionary named below. On the `--wav` path the shared positive-real
matrix fitter retains all accepted residues. Its dictionary is the travel-time
one (`fit_resolved`): resonators every c/(6D) across the band, D the closed
surface's largest extent, because a
kettle's load has detail every c/D in frequency that the fixed 32-term
log-spaced dictionary cannot follow (on the illustrative card it refuses at
complex RMS 0.033 against the unchanged 0.02 limit). `attach_exterior_radiation`
adds the resulting poles to the SAME
nonlinear equation before numerical preparation. Acoustic storage and damping
are already in the total energy ledger, not a post-step force or output EQ.

Receiver transfer fitting uses one explicit order. Even samples alone determine
poles and normalization; odd samples audit complex response (15% peak / 5% RMS,
relative to training scale) without selecting another order. Finite-point BEM
already includes distance spreading and reactive terms. Only a lower-bound travel
delay is peeled and restored once through the existing fractional-delay owner;
no second 1/r gain is applied. The receiver must lie outside the enclosing sphere
for this bounded observation path.

Playback observes accepted head velocity increments, uses the shared anti-alias
`Decimator`, executes the existing causal filters and uses the shared PCM/WAV
encoder. Mechanical rate must be an integer 2..32 multiple of 48 kHz; the
illustrative input uses 768 kHz. FIR and step-average acceleration timing are
retained and disclosed. Any mechanical, fit, power or clipping refusal publishes
no WAV. Discard this offline renderer on failure: a prefix may have advanced,
so it does not promise cross-owner rollback or resumable failed audio exports.

## Limits and validation status

Construction retains complete requested mode sets or refuses: at most 32 head
modes, 256 cavity modes, 512 cavity nodes and 2048 exterior panels here. The
interior air's near-field mass on a head mode comes from cavity modes far above
it, so `air_maximum_hz` must reach well past the head band: on kettle B a
480 Hz cut keeps no cavity mode of azimuthal order 2 or more, and the (2,1),
(3,1) and (4,1) head modes then feel no interior air at all; with a 1500 Hz
cut (68 modes) the coupled (2,1)..(5,1) ratios lie 0.5 to 1.1 % above a
2500 Hz basis. Exterior poles
are admitted while the joint state stays within 1024 scalars
(`radiation::MAX_JOINT_STATES`); acoustic coupling/clock limits still apply.
Each pole adds two states to the dense step, so the exterior dominates its
cost: 0.107 s per step on the illustrative card's 230 poles against about 2 ms
without them (pexfp.376.2 is the structural fix).
The head and axisymmetric pressure discretizations independently approach the
same circular interface under refinement. Four-site contact is not a resolved
moving pressure patch. The 20 J energy ceiling, 1e6 generalized-force ceiling,
1e-9 J absolute / 1e-6 relative balance tolerance and finite solve budgets are
admission guards, not real-time, bandwidth or model-accuracy certificates.

Both exterior fits approximate one stationary linear acoustic scene within the
sampled band; their individual audits do not prove joint observer/load accuracy
or continuous-band convergence. Nonlinear harmonics and loaded modes can extend
beyond it. Wider instruments require explicit head/cavity/surface/grid/clock
refinement. A narrow-band render is NOT a full-band calibrated timpani.

Independent Python equations checked the supplied piecewise-linear chamber's
volume (0.047705031260215386 m3), a P1/Fourier pressure reduction retaining the
uniform mode and first angular pair near 382.11 Hz, and the 768-panel closed
exterior. Polygonal outer volume agreed with the frustum sum to 2.22e-16 relative;
inner/outer head flux agreed to 2.78e-17 in a kinematic fixture. These are NOT
execution of the Rust implementation or measurements. Authored native tests
cover input/rim/enclosure refusal, uniform pressure work, felt impact/retry,
complete source power, transfer phase/audits, zero-load parity and passive load
state/history integration. The audio command above rendered natively at
5cf1fd1fe on ts1 (single-threaded, 1442 s): 4800 frames, peak 0.593 Pa, max
|Z-Z^T|/|Z| = 8.06e-3, lowest power eigenvalue +1.03e-6|Z| at 305.9 Hz,
246 load poles at complex RMS 0.0187 against the 0.02 limit. The margin is
thin since 88b352677 damped the fit's off-band storage poles; before it the
same impedance fitted with 230 poles at 0.0071.
`data/instruments/descriptions/timpani/` holds the measured card.

## Preserved zero-thickness comparison surface

The concurrently added `exterior.rs` and its four regressions remain available
unchanged as a separate zero-wall-thickness comparison. They use the supplied
INNER wall as a closed exterior and translate it about an explicit origin;
they do not infer copper thickness. The `--wav` path deliberately does not use
that approximation: `boundary.rs` uses the separately supplied outer stations,
rigid annulus and the original world coordinates. Its source fields and receiver
coordinates must never be mixed with the translated comparison geometry.

## Physical tensioning while the head rings

[Force-driven edge actuation](TENSIONER.md) adds an explicit effective inertia,
extra hardware stiffness, drag and sample-timed radial forces to the original
nonlinear head. The head reacts back on that coordinate, and its input work is
included in the shared energy ledger. Both CSV and coupled WAV use the same
trajectory. This is a small-strain uniform-edge image, not measured pedal or
hoop geometry; the unselected fixed-rim model remains unchanged.
