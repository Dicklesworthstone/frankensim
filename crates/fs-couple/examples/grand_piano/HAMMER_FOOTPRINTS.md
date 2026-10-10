# Finite and profiled hammer footprints

`--hammer-footprints faces.fshp` supplies a finite face on each selected key:
either a uniform span with two or four Gauss sites, or an authored profile with
one to four sites carrying explicit positions, recession, local felt thickness
and area fractions. Sites share the key's original hammer inertia and shank,
but retain independent existing WoolFelt conditioning and Prony recovery.
The actual string and soundboard displacements determine progressive contact.

This selector is available on `grand_piano` renders, `piano_exterior render`
and `render-loaded`, and the engine construction API. It composes with the
preset, supplied scales and hammer materials, supplied flat/crowned boards,
spatial dampers, MIDI/CSV performances, jack-force playing and the existing
acoustic observers. Omitting the selector preserves the original model.

## Complete physical input

The bounded UTF-8 file starts with `frankensim-hammer-footprints-v1`. Every key
in the loaded string scale needs exactly one `point`, `span` or `profile` record, including
unplayed keys. Empty lines and `#` comments are allowed; duplicate, missing,
unknown or malformed records refuse rather than acquiring defaults.

```text
frankensim-hammer-footprints-v1
point,60
span,69,0.012,4
```

This example is complete ONLY for a two-key scale containing 60 and 69. The
full preset needs all 88 records. `point,key` retains that key's existing
strike-point geometry. `span,key,length_m,sites` supplies a finite longitudinal
length in metres and exactly 2 or 4 Gauss sites. The length above is illustrative,
not measured Model D hammer geometry. No width is inferred from a key number,
material label, hammer mass or the source preset.

The face is centred at the scale's existing `strike_fraction * length_m`.
Its entire span must lie strictly inside the speaking string, not merely its
quadrature nodes. No clipping, shifting or automatic shortening is allowed.
For a `span`, the effective face is uniform and flat; changing its length at
fixed total area implies a different effective transverse width.

### Authored crown and local thickness

A `profile,key` declaration is followed by one to four explicit site rows:

```text
frankensim-hammer-footprints-v1
profile,69
site,69,-0.0015,0.0003,0.006,0.15
site,69,-0.0005,0,0.008,0.35
site,69,0.0005,0,0.008,0.35
site,69,0.0015,0.0003,0.006,0.15
```

This illustrative fixture is complete only for a one-key scale containing 69;
it is not measured Model D crown geometry. Each site row is
`site,key,offset_m,recession_m,thickness_m,area_fraction`:

- `offset_m` is the longitudinal offset from the scale's existing strike station.
  Rows must be strictly ordered and every resulting station must lie strictly
  inside the speaking string. Positions are neither moved nor clipped.
- `recession_m` is the nonnegative distance behind the nominal hammer crown.
  Zero marks an unrecessed site. Compression is crown displacement minus local
  string displacement minus this recession. More recessed sites engage only
  when their actual gaps close.
- `thickness_m` is the positive local undeformed felt-column thickness. It is
  independent authored geometry, not inferred from recession or a backing shape.
- `area_fraction` is a positive fraction of one string's allocated felt area.
  Fractions must sum to one within `1e-12`; they are never normalized.

All fields must be finite. A site requires an earlier profile declaration for
its key. Missing sites, a fifth site, nonrepresentable separation or derived
material geometry, duplicate key selections and incomplete scale coverage refuse.
The same profile is applied to each original unison member.

The scale supplies total felt area and hammer mass. It also supplies the
thickness for `point` and `span`; a `profile` supplies each local thickness.
`--hammers` supplies the existing constitutive properties. Area is divided among
the original unison members, then partitioned by the span's positive Gauss
weights or the profile's authored fractions. Force, elastic storage, Prony
stiffness, permanent crush and densification use each site's physical area and
thickness. Adding sites does not multiply total felt area or hammer mass.
Each duplex segment stays outside direct hammer contact.

Profiles represent discrete parallel compression columns on a supplied face.
The engaged set changes through actual unilateral contact, without an imposed
area envelope. They do not resolve a continuous three-dimensional crown,
rotating contact normals, lateral friction or a full hammer action.

```sh
# faces.fshp must cover every key in the supplied strings.csv.
cargo run --release -p fs-couple --example grand_piano -- \
  --scale strings.csv --board-geometry panel.fsb --hammers materials.fsh \
  --hammer-footprints faces.fshp --performance score.csv \
  --dampers estimated --duration 6 --render spatial-hammer.wav

# faces.fshp must instead cover all 88 keys for this full-preset invocation.
cargo run --release -p fs-couple --example grand_piano -- \
  --preset steinway-d --hammer-footprints faces.fshp --midi score.mid \
  --microphone 0.2,0.8,1 --microphone-right 1.2,0.8,1 \
  --duration 30 --render spatial-hammer-stereo.wav
```

These are entry points, not assertions of completed or calibrated renders.
The same original force, densification, convergence and energy limits apply.
Output pressure is not normalized to compensate for a changed physical face.

## Reciprocal mechanics and material memory

For every site the same moving-boundary string representation supplies its
modal shapes and residual bridge lift. The contact compliance retains both
same-string cross-site response and the full reciprocal soundboard contribution.
All site forces enter the existing exact-ZOH/Schur mechanical step together.
The single hammer receives their summed reaction. Averaging string displacement
BEFORE the nonlinear material law would erase local contact near nodes; that
shortcut is not used.

Each site evaluates the existing discrete felt force and its existing series
Prony recovery at its local thickness. Elastic strain is the local overlap less
Prony deformation, divided by that thickness; elastic storage uses local
area times thickness. Static recession belongs only to the overlap, so it does
not change the work-conjugate string displacement or create external work.
These same gaps participate in the simultaneous hammer solve and the optional
nonlinear string-stretching iteration. Histories are independent across
positions and unison members. Released sites keep recovering without erasing
permanent crush.

The published `R_H` selector remains a penetration-based coefficient in
`N s / m^p`, not an inferred strain-rate material. Profiled recession with the
original uniform course thickness can use it, distributing each per-string
coefficient by the existing area fractions. A profile with any changed local
thickness explicitly refuses this source selector. No local-rate profile input
is provided; varying thickness can use the existing supplied WoolFelt/Prony
material cards without guessing a conversion of `R_H`.
Physical work and loss enter the ordinary instrument accounting once; no
separate gain, damping correction or time integrator is introduced. Refusal
restores every contact history, hammer/jack state and mechanical coordinate.

Una corda selects whole unison strings: all quadrature sites on a selected
string stay enabled. The remaining area is NOT renormalized after excluding a
string. Sustain and sostenuto retain their original damper mechanics. Repeated
notes preserve prior contact conditioning and surviving instrument vibration.
Both stereo microphones consume the same substep board trace and score clock.

Two-point Gauss quadrature exactly integrates polynomial fields through degree
three on the uniform span; four-point quadrature through degree seven. This
is not exact integration of sinusoidal modes or nonlinear contact stress.
Authored profile sites carry no automatic quadrature-accuracy claim.
Compare the explicit resolutions, mechanical timesteps and retained mode bands
before claiming convergence. The preparation ceiling is 1,056 contacts (four
sites on three strings for each key), with at most 512 stored shape coefficients
per string site. Dense contact coupling can be costly; no real-time claim is
made. Existing acoustic approximation/bandwidth limits remain unchanged.

## Focused native checks

```sh
cargo test --release -p fs-couple --example grand_piano hammer_footprint -- --test-threads=1
cargo test --release -p fs-couple --example grand_piano footprint_tests -- --test-threads=1
cargo test --release -p fs-couple --example grand_piano footprint_render_tests -- --test-threads=1
cargo test --release -p fs-couple --example grand_piano hammer_profile_render_tests -- --test-threads=1
```

Seven bank tests cover complete input, span quadrature, authored profile
geometry, nodal contact, reciprocal compliance, actual force/work closure and
original point parity. Nine engine tests additionally exercise area and local
thickness scaling, unloaded rest, progressive crown engagement with nonlinear
string stretching, independent material history, una corda, re-strikes,
transactional refusal, published-rate admission and jack drive. Frontend tests
exercise supplied geometry and actual MIDI/CSV-to-stereo-to-PCM playback on one
physical timeline. These checks require native execution; syntax parsing does
not establish Rust compilation, passing dynamics, measured realism or acoustic
convergence.
