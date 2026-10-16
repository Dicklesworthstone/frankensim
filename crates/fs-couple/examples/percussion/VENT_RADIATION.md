# Prescribed-flow vent radiation through the existing BEM observer

`--prescribed-vent-radiation` explicitly enables a one-way exterior observer
for a drum/snare with `--cavity-modes` and `--cavity-neck`. It works on existing
`-wav` and `-mic` commands, including stretching heads, two sticks, fixed
mufflers, supplied drum geometry, acoustic drag and stereo receivers. Existing
compliant-mute combinations remain available on the drum paths that already
admit them. Bare vented audio still refuses rather than silently choosing this
approximation. The flag is not accepted on mechanics CSV or cymbals.

```sh
# Illustrative SI inputs, not a measured opening or loss calibration.
(set -C; cargo run --release -p fs-couple --example percussion -- \
  snare-mic 4800 20 0.08 0.05 0.35 --cavity-modes \
  --cavity-neck 0.003 0.008 1000 0.4 0.08 \
  --prescribed-vent-radiation --microphone-right -0.12,0.05,0.4 \
  > vented-snare.wav)
```

Add `--head-stretching --analytic-newton --impact-substeps 8 511` to select the
existing joint nonlinear heads/wires rather than the original linear-head
snare. This does not change which acoustic approximation is selected. The
same original force programs advance on one mechanical clock, independently
of receiver count. A completed render is subject to all original contact,
energy, slope, BEM and held-out fitting gates. These example invocations are
not claims of a completed native render or a real-time deadline.

## What is coupled, and what is deliberately not

The mechanics already solves head/cavity/neck exchange. Its neck coordinate
has acoustic inertance `L = rho * effective_length / area`, physical outward
volume `v = q / sqrt(L)` and flow `Q = p / sqrt(L)`. The observer now obtains
this exact port from the admitted mechanical neck. It never treats a
mass-normalized momentum as a surface velocity or pressure directly.

The exterior's rigid wall strip is retriangulated around a circular aperture.
A prescribed normal-velocity disk replaces the rigid cap there. Its panel
weights integrate to `1/sqrt(L)`, so the discrete outward flux is exactly Q,
including its sign. This is geometric area quadrature, not an output gain or
per-channel normalization. Head triangles and source weights are retained.
The aperture is an additional source in the SAME BEM batch and exterior Green
representation; scattering from the remaining boundary is included. It is not
a free-field monopole waveform summed after head rendering. Both receivers
retain their independent fits, travel times and pressure histories.

**This is prescribed-flow, one-way radiation, not a fully radiation-loaded
vent model.** The mechanical neck still discharges into its declared zero-gauge
reservoir. Computed exterior pressure is NOT fed back to that neck or to either
head. Radiation power is therefore not a new term in the mechanical loss
ledger. The mechanical supplied effective length and resistance remain exactly
as provided; no second end correction or fitted radiation resistance is added.
Use this observer only where that prescribed-motion approximation is justified.
A coupled pressure/flow aperture impedance and radiative backreaction remain
unfinished. The flag does not claim to close that physics gap.

## Geometry and bandwidth admission

The center uses the supplied sidewall azimuth and the existing axial coordinate
measured from the batter into the cavity. In the exterior frame its height is
`depth/2 - axial_position`. The exterior cylinder is faceted. This first
implementation admits a circular opening wholly inside ONE planar sidewall
facet, though it can cross that facet's axial mesh subdivisions. It refuses
angular seams, the head/rim and intersected moving panels instead of snapping
the opening to a different location. The planar center is the specified
azimuth's intersection with that actual faceted wall, not a new smooth-cylinder
surface. Angular mesh refinement and specimen geometry are separate concerns.

The circular polygon has at least 32 segments plus required perimeter rays.
Its area must be within one percent of the supplied throat area. Every original
outer edge vertex is retained exactly; the closed boundary must remain oriented
and watertight. The old 2048-panel and 63-source ceilings apply after insertion.
No head mode is discarded to make room for the aperture.

The compact neck must satisfy `k * max(radius,effective_length) <= 0.3` over the
FULL unchanged 40..1640 Hz BEM fitting window, not only the interior cavity
window. At the existing air sound speed, a 12 mm effective length exceeds that
radiation guard; the example above declares 8 mm rather than bypassing it.
Increasing the allowed neck size needs a distributed acoustic chart, not a
looser fit tolerance. Constant disk velocity, constant-panel BEM, an undeformed
boundary, a narrow fitting band and one-way acoustics remain approximations.
Passing admission does not establish spatial/frequency convergence, full-band
sound quality, measured SPL or calibration.

## Verification

Four focused source tests cover exact head/neck flux, watertight geometry and
area refinement, physical coordinate recovery, actual BEM superposition and a
low-frequency compact limit. Three integration tests cover explicit command
selection, unchanged two-stick/snare mechanics and energy/time, and a complete
BEM-to-fitted-stereo-to-PCM render. They are native regressions to run, not a
claim that an independent arithmetic test executed the Rust program:

```sh
cargo test --release -p fs-couple --example percussion acoustics::aperture -- --test-threads=1
cargo test --release -p fs-couple --example percussion vented::tests -- --test-threads=1
```

Native execution is requested in the existing percussion workflow. A separate
400-case Python geometry check passed oriented-area and exact edge-pairing
checks; it does not validate native compilation, BEM accuracy or audio fitting.
