# Close microphones at their actual positions

Add this explicit row to a `frankensim-piano-exterior-si-v1` acoustic file:

```text
receiver-evaluation,near-field
```

The existing `receiver-m,x,y,z` rows remain physical SI coordinates in the
board frame. The selected evaluator works in `response`, `admittance`, `render`
and `render-loaded`, with supplied BODY meshes, native section skins, and
explicit rigid lid/cabinet assemblies. It permits exterior points **inside the
body's enclosing sphere**, including close positions above or below a board.
No microphone is moved outward or replaced with a far-field direction.

Omission, or `receiver-evaluation,centroid`, keeps the original centroid-panel
observer and enclosing-sphere admission. This remains a deliberately separate
comparison path. An invalid or under-resolved near-field request refuses; it
never falls back to the centroid observer or a different microphone position.

## A source-derived Model D study

A complete estimated setup is included as `model-d-close-mics.fspe`. Its two
receivers are 0.15 m above the reference board plane, not measured factory or
recording-session positions. It retains the earlier study's medium, sampled
band, structural mode window, fit order and physical PCM scale.

```sh
cargo run --release -p fs-couple --example grand_piano -- \
  --preset steinway-d --mesh-divisions 4 --dump-geometry model-d.fsb

cargo run --release -p fs-couple --example piano_exterior -- \
  response model-d.fsb steinway-d board-skin-continuous \
  crates/fs-couple/examples/grand_piano/model-d-close-mics.fspe close-response.csv

cargo run --release -p fs-couple --example piano_exterior -- \
  render-loaded model-d.fsb steinway-d board-skin-continuous \
  crates/fs-couple/examples/grand_piano/model-d-close-mics.fspe close-piano.wav 2 \
  --note 69 --velocity 1 --modes 128 --substeps 8
```

Existing material, hammer-face, nonlinear-string, MIDI/CSV, pedal and rigid-part
options still apply. Reuse the scale that generated any loaded equilibrium.
These commands exercise the physical path, not a promise that a particular
full-instrument BEM discretization or passive fit will pass its existing bounds.
Use fresh output paths. There is no normalization, microphone gain, room model,
or empirical close-mic equalizer added by this option.

## Geometry and propagation

`fs_bem::near_field::Geometry` admits the complete final triangle surface.
It welds only exactly equal coordinates (with signed zero canonicalized),
requires closed, consistently outward components, and checks each receiver's
solid-angle winding against **each component separately**. A point inside a
rigid lid is not exterior merely because that lid has zero prescribed velocity.
Open, inward, coincident and ambiguous cases refuse rather than being repaired.
The input contract still requires disjoint, non-self-intersecting components;
these checks are not a global intersection or cavity-accessibility certificate.

Clearance is the minimum Euclidean distance to an actual triangle face or edge,
not its centroid. The piano wrapper requires at least **two output samples of
flight distance** from every component. At 48 kHz and 340 m/s this is about
14.17 mm. A scale-dependent floating-point guard also rejects unresolved tiny
gaps. The maximum allowed propagation lower bound remains 0.5 seconds.

In near-field mode the retained flight lower bound is the minimum distance to
any panel divided by the sound speed. Including rigid scatterers makes this a
conservative bound on the last propagation leg, not an inferred direct path or
an occlusion ray. The observer removes that delay only for fitting, then restores
it through the existing causal receiver delay. It does not remove the physical
phase or reactive near-field pressure from the solved transfer.

## What is integrated

The boundary solve, its constant panel pressure and velocity traces, and its
radiation-reaction matrix are unchanged. The receiver evaluates the same Green
representation with the fixed `exp(-i omega t)` convention:

```text
p(x) = sum_panels integral [ dG(x,y)/dn_y * p_panel
                            - G(x,y) * i*omega*rho*v_panel ] dS_y
```

Adaptive 4x4/8x8 Gauss rules integrate each retained triangle. A separate
geometric-distance and wavelength condition prevents both rules from missing
a narrow nearby peak or an under-sampled oscillation. Local acceptance uses
relative discrepancy 1e-7 against the integral of absolute kernel magnitudes.
Preparation refuses after 14 subdivision levels or two million kernel evaluations
across all receivers/panels at a frequency. These are explicit cold work bounds,
not unbounded refinement or a new solve. Each prepared receiver row is shared by
all modal fields at that frequency.

The core API exposes propagated quadrature discrepancy estimates for evaluated
pressure. They estimate integration error in the supplied panel traces only;
they do **not** bound boundary-solve, mesh, modal truncation or physical-model
error. Cancellation can make pressure's relative error larger than individual
kernel tolerances. The piano's existing receiver-fit, passive-load, temporal-band
and wavelength checks remain separate and unchanged.

A microphone relocation cannot change the surface solve, modal load, string
state or combined mechanical/acoustic energy. `render-loaded` fits the same
full-matrix impedance irrespective of receiver position. Both stereo channels
observe every accepted mechanical substep on one clock. Tests compare relocated
receivers against an independently advanced copy of the same piano, including
acoustic storage and loss when feedback is selected.

The soundboard is still linearized about its supplied equilibrium; added lid
and cabinet surfaces remain rigid. This adds no microphone diaphragm/directivity,
flexible cabinet, absorbing finish, room, factory mesh, or material measurements.
Closer evaluation makes spatial discretization error more visible, not smaller:
full-band or measured Steinway fidelity is not established by passing this
receiver's geometric and quadrature checks.
