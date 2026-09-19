# Geometry to radiated pressure: an offline percussion reference

The example now composes the existing mechanical and exterior-acoustic owners.
**It has not been natively compiled or executed in this implementation session:**
`cargo` and `rustc` were unavailable. The new commands, tests and physical path
are implemented, but no native WAV, listening result, successful BEM bake or
real-time benchmark is claimed. A failed solve/fit aborts export, not a fallback
sound. This is not yet a full-band or calibrated instrument.

## Commands

From the workspace root, on the repository's supported Rust toolchain:

```sh
# One second of 48 kHz mono pressure, 1 Pa mapped to digital full-scale.
cargo run --release -p fs-couple --example percussion -- splash-wav 48000 1.0 > splash.wav
cargo run --release -p fs-couple --example percussion -- drum-wav 48000 1.0 > drum.wav

# Finite-distance microphone: actual point pressure, not a far-field shortcut.
cargo run --release -p fs-couple --example percussion -- splash-mic 48000 1.0 > splash-mic.wav
cargo run --release -p fs-couple --example percussion -- drum-mic 48000 1.0 0.08 0.05 0.35 > drum-mic.wav

# The original mechanical CSV paths remain available at their original 2 us step.
cargo run -p fs-couple --example percussion -- splash 4096 > splash.csv
cargo run -p fs-couple --example percussion -- drum 4096 > drum.csv

# Added skin and acoustic-composition regressions.
cargo test -p fs-plate shell::reduction::radiation
cargo test -p fs-couple --example percussion
```

The WAV count is audio frames, not mechanical steps; maximum 480000 frames.
The full-scale pressure is explicit and strictly positive. No peak normalization
or hidden gain adjustment is applied. Clipped sample count and peak pressure
are printed to stderr. The candidate WAV is written only after every requested
mechanical and acoustic sample succeeds. A redirection shell can still create
an empty destination file if the program refuses. Native runs may be expensive:
41 dense BEM frequency solves precede an allocating reference mechanical solve.

## Physical path and reusable pieces

`ShellReduction::radiation_surface` constructs both physical faces of the same
curved, tapered shell plus its outer rim and mounting-hole walls. Physical axial
rotations contribute through `u_skin = u_mid + theta cross offset`; the side
walls do not silently lose rotary motion. Mode normalization and ordering remain
those of the mechanical reduction. A pressure-force projection is the negative
area-weighted transpose of outward velocity, preserving the sign of interface
power. This projection is available for future loading, but the example does
not yet apply the exterior pressure back to the mechanics.

The drum adapter constructs the two existing film meshes, rigid bearing-edge
annuli and a subdivided outer cylindrical shell at the supplied actual outer
radius and depth. It does not acoustically replace the shell with a cylinder
at the smaller clear film radius. Both mechanical film coordinates point down;
outward top velocity is negative, bottom velocity is positive. The diagnostic
internal pressure is positive for net compression; its former CSV sign was
reversed. The mechanical volume-spring Hamiltonian is unchanged.

The boundary is passed to `fs-bem::helmholtz::solve_radiation_batch` as exact
triangles. One frequency shares its matrix/factorization across all modes.
Normal velocities for a unit generalized acceleration obey `v = i a/omega`
under BEM's `exp(-i omega t)` convention. Materially negative radiated power
beyond the solver's roundoff interval is refused, as are BEM work/resolution
limits. The default splash and drum each have 1024 exterior panels.

The `*-wav` observer at `[1.5, 0.7, 1.5]` metres is an explicit **far-field** receiver,
not a close drum microphone. Its direction is projected before vector fitting,
so it needs one proper filter per retained generalized acceleration rather than
a full spherical-harmonic bank. `fs-vfit` owns fitting, stability and Tustin
realization. Negative-time BEM responses are conjugated and physical frequencies
are warped before fitting; omitting either changes the transfer's phase.

The origin-referenced far field can contain a geometric advance from the near
side of a finite source. It is delayed by the enclosing radius divided by sound
speed before fitting. The existing propagation line then applies only the
remaining `(range-radius)/c` delay and `1/r` once. The complete frequency-domain
phase remains that of the original field propagated by `range/c`. The delay
line uses its existing two-tap fractional interpolation, not an exact all-band
delay. The ten-source-radius admission is a screening rule, not a near-field
error certificate.

The separate `*-mic` commands use `fs-bem::helmholtz::exterior_pressure_at_points`
at an actual finite position. Its Green representation retains near-field
terms as well as outgoing waves; the sound is not approximated by taking a
far-field direction and dividing by distance. The default microphone is
`[0.08, 0.05, 0.35]` metres in the geometry frame, with x/y in the head plane
and z upward. An optional complete x/y/z triple follows frame count and
full-scale. For the drum, the heads are at z = +/- depth/2.

This owner already returns physical pressure INCLUDING spreading and travel.
Only a lower-bound travel delay `(range-enclosing_radius)/c` is peeled before
fitting and restored once by the existing propagation line. **No additional
1/r gain** is applied. A conservative rule requires the microphone to be
outside the source's enclosing sphere with at least two samples of remaining
propagation. Some physically exterior positions closer to a face are therefore
not admitted. This is a point-pressure reference, not a microphone capsule,
proximity-effect/electronics model or moving microphone. The same fit and mesh
accuracy limits still apply.

Mechanical steps are `1/768000 s`; the new renderer observes interval-average
modal acceleration without adding a force pulse. The existing causal integer
`Decimator` filters this observation by 16 before the 48 kHz transfer filters.
It does not filter the contact forces or modify the nonlinear energy state.
Its extra 47-audio-frame latency is retained; interval averages are tagged at
step end. The implicit mechanics remains the original `ImpactSystem`/`fs-phs`
owner. The observer errors are terminal for this example; no cross-owner retry
or rollback is asserted. Export uses the existing pressure PCM encoder.

## What the gates establish, and what they do not

There are 21 training and 20 interleaved held-out BEM frequencies from 40 to
1640 Hz. Per-input maximum and RMS complex errors are normalized by that input's
largest sampled response, with limits 0.15 and 0.05. These are authored numerical
fit tolerances, not measured perceptual tolerances or a continuum certificate.
A refusal requires better resolution/order/sampling, not disabling the gate.
The same mesh supplies the training and held-out data: mesh error does not
vanish because a filter agrees with those data. The reference fit is evaluated
outside its checked band by a transient; **above-band transfer error is unknown**.
The output anti-alias filter does not restrict the input to the 1640 Hz fit band.

The retained mechanical bands are still only 50–1200 Hz for the splash and
80–500 Hz for each head. Missing high-frequency and in-plane modes, nonlinear
transfer convergence, fitted-filter above-band behavior, contact-time-step
convergence and aliasing remain major accuracy gates. A 48 kHz WAV header does
not turn this basis into a full-band cymbal crash.

Exterior radiation is one-way and linear about the undeformed boundary. No
radiation mass/damping back-coupling, deformation-updated normals, microphone
capsule/electronics, air absorption, room, snare wires, vent flow, flexible maple
shell, full stand rocking, hand/grip model or flexible/anisotropic hickory stick
is inferred. The original cavity uses its declared approximate bulk modulus;
exterior dry air comes from `Medium::air()`. A complete timpani needs bowl/head/
air modes and their loading, not just rescaling the example's volume spring.

Geometry, material and contact provenance remain in [README.md](README.md).
The Zildjian manufacturer page confirms Z5A length, diameter, hickory, taper
category and oval tip, not its complete radius profile or tip compliance.
Kaselouris et al. supply the splash's diameter, bell diameter, mounting hole,
edge thickness and literature B20 constants, not a surveyed complete taper or
hammer map. Their paper uses a different, 394 mm stick. These facts were checked
against the primary sources on 19 September 2026:

- https://zildjian.com/products/5a-drumsticks
- Kaselouris et al., *FEM-BEM Vibroacoustic Simulations of Motion Driven
  Cymbal-Drumstick Interactions*, Acoustics 2023, DOI 10.3390/acoustics5010010;
  https://www.mdpi.com/2624-599X/5/1/10

Native regressions have been added for skin closure/rotation, reciprocity, drum
outer dimensions and opposite head motion, acceleration phasors, held-out fit
phase, source-delay decomposition, propagation scale/delay and compression sign.
Independent Python geometry checks reproduce the polygonal drum volume to
2.43e-17 m^3 and head volume derivatives to 1.75e-10 relative error over four
meshes. An independent real rational least-squares oracle reproduces a known
two-state transfer on the withheld warped grid to 1.30e-15; it is **not execution
of fs-vfit**. These checks do not establish compilation, native acoustic quality
or real-time performance.

Three additional Rust regressions cover finite-point admission versus far field,
no double distance gain, and the existing Green evaluator plus delay peeling.
An independent 64-by-128 spherical surface quadrature of a known translating
(dipole) field matches its analytic finite-distance pressure to below 1e-13
relative error at four microphone ranges. At k*r=0.2 its pressure magnitude is
5.099 times the far-field approximation, so the distinction is load-bearing.
That independent calculation is not execution of the Rust BEM or its new tests.
