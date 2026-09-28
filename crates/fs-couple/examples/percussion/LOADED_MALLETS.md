# Flexible shafts with physical felt heads

Use a `frankensim-felt-mallet-v2` card together with the corresponding
`--flexible-stick` or `--second-flexible-stick`. The head's actual mass and
rotary inertia load the original beam pencil BEFORE modal reduction. The
finite felt face then applies reciprocal forces AND moments to that same
loaded basis. There is no second free head mass, parallel Hertz contact,
extra oscillator, recorded attack, or output EQ.

```sh
cargo run --release -p fs-couple --example percussion -- \
  drum-stretch 512 --strike-position-m 0.06 0.01 --strike-speed-m-s 0.5 \
  --flexible-stick crates/fs-couple/examples/percussion/estimated-flexible-stick.fst \
  --mallet-spec crates/fs-couple/examples/percussion/estimated-loaded-mallet.fsmallet \
  --analytic-newton --impact-substeps 8 511
```

Both hands can independently select this composition. Use
`--second-flexible-stick`, `--second-mallet-spec` and an explicit
`--second-stick-position-m` for the second hand. Drum, stretching drum, snare,
snare-off and single-cymbal commands retain their existing CSV, WAV and
microphone forms. Modal-only and paired-hi-hat felt paths still refuse.
These examples are invocation recipes, not claimed completed recordings.

## An explicit physical head, not a reinterpreted effective mass

```text
frankensim-felt-mallet-v2
geometry,0.02,0.012,0.006,0.00002
attachment,0.0000012,0
felt,100000,0.2,2.2,3,0.15,0.7
conditioning,0
creep,3000,6
```

V2 `geometry` means **actual moving head mass [kg]**, face radius [m], felt
thickness [m], and initial clearance [m]. The shaft profile must EXCLUDE this
head. Its centre of mass is at the shaft's declared contact station. The mass
includes moving head/felt material, not an already reduced contribution from
shaft or hand. Material-history coordinates do not add another head inertia.

`attachment` supplies the positive central transverse moment of inertia
[kg m2] and shaft azimuth in the target's XY plane [radians, -pi..pi]. This
moment is about the head's centre, normal to the bending plane; the parallel-axis
term is derived from the actual pivot/contact separation. No homogeneous-head
inertia is inferred from the face radius or felt thickness. Off-axis centres
of mass require a different inertia operator and are not represented here.

All records have the same strict duplicate, finite-value, size and material
validation as v1. `attachment` is required exactly once in v2 and forbidden in
v1. A v2 head without its shaft refuses. A v1 effective-mass card with a shaft
also refuses: its mass may already include shaft/grip reduction and cannot
safely be added again. Existing v1-only performances are unchanged. The supplied
v2 card is a synthetic estimate, not a measured commercial mallet.

## Mass, face rotation, and reciprocal work

The beam owner adds `m N N^T + J N' N'^T` to its source mass pencil. It separates
true rigid rotation using the LOADED mass metric and recomputes the complete
requested elastic slice. Thus initial inertia, frequencies, damping, tip and
hand participation all describe the same shaft/head. Both bare and loaded
preparations must fit the existing input bounds; no mode is discarded to make
an attachment fit. This cold construction is not performed on an audio thread.

Each of the four existing positive-area sites has inward displacement
`u_tip + xi * theta_tip`, where `xi` is its signed offset along the declared
shaft direction. The conjugate generalized force includes the site's moment
`F * xi`. Unequal contact across the face can therefore rotate the head and
bend the shaft; averaging the face into a single force would lose that torque.
The original independent WoolFelt and Kelvin histories remain at every site,
with the same area partition and whole-face creep stiffness/viscosity.

The face starts horizontal at the supplied clearance. Its affine angular
reference is the initial rigid shaft angle: a site uses
`u_tip + xi*(theta_tip-theta_initial)`. This prevents initial rigid positioning
from introducing fictitious unequal preload. It is a fixed reference offset,
not a per-step correction or reset. Both initial displacement and velocity are
placed in the combined rigid mode; no initial bending energy is manufactured.

Drumhead and finite-thickness cymbal-skin rows, including shell force moments,
retain their original signs and geometry. The cymbal face retains its complete
footprint-derived reference plane and individual clearance at each site.
Physical hand forces still use the hand row, not the head centre or footprint
average. Accepted-step rollback, felt conditioning, internal refinement, and
energy/work/loss accounting stay with the existing time owner.

## Model limits and checks

This is a planar, small-rotation, fixed-reference contact image, not unrestricted
3D mallet motion. Only normal motion and one bending-plane head rotation are
represented. In-plane footprint travel, evolving normals, torsion, tangential
friction, a compliant head/shaft joint, shear, off-axis COM, and large deflections
are not inferred. The fixed four-site quadrature and retained bandwidth still
need convergence studies. Loading a head does not calibrate its felt or damping.

No head/shaft coordinate is a new microphone source or a direct cavity-pressure
input. Its audible contribution passes through physical contact with the
existing drumheads or shell. Original head/shell, snare and carrier addresses
remain fixed; cavity placement uses the enlarged structural prefix.

```sh
cargo test --release -p fs-plate --lib shell::stiffened::beam
cargo test --release -p fs-couple --lib render::plate::impact::striker::flexible
cargo test --release -p fs-couple --example percussion mallets::shaft -- --test-threads=1
```

Focused regressions cover loaded mass orthogonality and frequencies, rigid
launch energy, force/moment work, strict v1/v2 admission, reference gaps, actual
felt-contact bending and history rollback, and both heads in complete
snare/cavity and cymbal-source layouts. Native execution is requested in the
existing percussion workflow; it is not established by the independent
numerical mass-pencil check. No measured fidelity or real-time claim is made.
