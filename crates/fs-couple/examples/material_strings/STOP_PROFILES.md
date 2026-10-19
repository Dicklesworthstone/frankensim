# Stopped plucks and bows on the same material instrument

The optional `--stop-profile PROFILE` binds fixed compliant stops BEFORE the
existing force preload and bow schedule are constructed. Use it alongside
`--body-profile` in either order. Omitting it keeps the original path.

The current two-plane and rotational-bridge constructions are retained. A
stop uses either the source's transverse/bow direction or an independently
declared normal in the same transverse frame, without adding a normal-contact
torque. A separately selected rotational bridge continues to react through its
original receiver geometry.

After compiling the material packs described in README.md, run from the root:

```bash
cargo run --release -p fs-couple --example material_strings -- \
  /path/to/pluck.material-strings /path/to/new-stopped.wav \
  --stop-profile crates/fs-couple/examples/material_strings/stopped-pluck.stops
```

The output path must not exist. Material/geometry/contact/control/preload
admission finishes before creating it; a later physical or I/O failure may
leave an incomplete prefix, as in the original performer.

## Profile records

```text
frankensim-string-stops-v1
stop support 0 0.43 positive fixed 0 5000 1 0.2 10 0.001
```

The fields after `stop` are: unique ID, zero-based string index, station as a
fraction of its COMPLETE span, side, receiver, unloaded surface clearance [m],
contact stiffness [N/m^alpha], exponent alpha, internal loss [s/m], maximum
reaction [N], and maximum penetration [m]. `positive`/`negative` choose the
closure direction along the source's bow tangent. These are physical inputs
and refusal ceilings, not a selected pitch or clipped contact force.

`fixed` explicitly selects an immovable support. `bridge:0` instead takes the
actual body displacement/force footprint of string zero's sole declared bridge.
This is that receiver's same-basis shape, not an extra force on string zero.
It works with the existing uniform or regional body profiles. Ambiguous multiple
bridge links refuse. Selecting this receiver does not invent a fingerboard or
relocate the original bridge. Arbitrary receiver footprints remain a Rust API.

Contact coefficients are authored declarations, NOT measured finger/fret data
and not derived from bulk elasticity. The supplied profile is synthetic. The
original contact-law producer owns any material correspondence; the Rust API
also accepts the original source-backed obstacle receipts. Strings and bodies
retain their material packs, geometric projections and shared gas.

## Independent pluck and stop directions

A pluck or fingerboard normal need not coincide with the bow tangent. The
existing force record accepts `along DX DY` instead of `centreline` or `surface`:

```text
port 0 0.43 along 0 1 0.001
```

This supplies 0.001 N along transverse plane two, at station 0.43 of string zero.
It is a centreline force: no twisting torque is added. `centreline` still means
the source's bow direction; `surface` still means the original bow-facing
surface, including its torsional velocity and conjugate torque where selected.
An arbitrary surface force would additionally need its lever arm and is not
inferred from a direction alone.

A stop record can replace `positive`/`negative` with `normal NX NY`:

```text
frankensim-string-stops-v1
stop support 0 0.43 normal 0 1 fixed 0 5000 1 0.2 10 0.001
```

All direction components use the common orthonormal transverse frame declared
by the source's `polarizations two` record. They must be finite unit vectors
(squared norm within 1e-12); they are not normalized or treated as gains.
Nonzero second-plane components require a two-plane material source. A
one-plane source refuses them instead of discarding the perpendicular load.
Neither command changes the material modes, bow tangent or bridge direction.
The same material-derived projection applies force and measures displacement.

For `normal NX NY` with `bridge:I`, the selected bridge's scalar axis must be
collinear with that normal in this common frame. Opposite orientation negates
its original receiver row; an oblique or orthogonal axis refuses. A scalar
bridge response cannot supply missing directional body motion. General
receiver maps remain available through the Rust API, with the producer owning
that geometric correspondence. Existing positive/negative records keep their
previous interpretation.

Copy `directional-pluck.material-strings` next to the compiled wire/body packs,
then render the supplied synthetic stopped-pluck fixture:

```bash
cargo run --release -p fs-couple --example material_strings -- \
  /path/to/directional-pluck.material-strings /path/to/new-directional.wav \
  --stop-profile crates/fs-couple/examples/material_strings/directional-fingerboard.stops
```

The fixture uses an X-directed bow, a Y-directed load and stop, and an oblique
bridge. This is a geometry/solver integration example, not measured fingerboard
or instrument calibration. Changing only the stop normal is an actual contact
intervention, not a pitch or output-filter edit.

## Preload and playback

`forces preload` solves the complete stationary network AND all attached normal
contact potentials with the bows released. Contact activity is solved rather
than supplied. All original static force, penetration, energy and sweep limits
apply. The dynamic compliance matrix is not used as a static inverse stiffness.
The returned energy is an initial condition, not measured pre-window work.
Holding that load should retain equilibrium; a sample-zero or later force event
releases it into the same bow/string/body/cavity dynamics. No modal state, contact
law, pressure observation or sample clock is reset to produce the pluck.

The binding also admits the original zero-state geometry before preload.
Excessive initial overlap can therefore refuse even when a hypothetical later
load would relieve it. Existing motion cannot be overwritten by preload.

Profiles are limited to 16 KiB and 32 normal stops; the original combined
normal/friction port and setup limits remain stricter where applicable. Contact
geometry stays fixed. This does not solve moving fingers, friction at the stops,
changing contact area, nonlinear string tension, rigid fretting, or bow normal
load distribution. Truncated modal/contact time refinement and material-domain
checks remain necessary; the existing acoustic observation is inherited.

## Focused native checks (not executed in the authoring environment)

```bash
cargo test -p fs-couple --lib friction::multiple::obstacles -- --test-threads=1
cargo test -p fs-couple --test material_bow_interface stops:: -- --test-threads=1
cargo test -p fs-couple --test material_bow_interface directions:: -- --test-threads=1
cargo test -p fs-couple --example material_strings stopped:: -- --test-threads=1
```
