# Finite longitudinal hammer footprints

`--hammer-footprints faces.fshp` replaces the point contact on a selected key
with two or four independent felt contact sites along each speaking string.
They share the key's original hammer inertia and shank, but each site retains
its own existing WoolFelt conditioning and Prony recovery history. The actual
string and soundboard displacements determine contact; no attack envelope,
point-response equalizer, direct hammer sound, or extra oscillator is added.

This selector is available on the `grand_piano` render command and its engine
construction API. It composes with the preset, supplied scales and hammer
materials, supplied flat/crowned boards, spatial dampers, MIDI/CSV performances,
jack-force playing and mono/stereo observers. Other example frontends keep their
existing contact selection. Omitting the selector preserves the original model.

## Complete physical input

The bounded UTF-8 file starts with `frankensim-hammer-footprints-v1`. Every key
in the loaded string scale needs exactly one `point` or `span` record, including
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
The effective longitudinal face is uniform and flat; changing its length at
fixed total area implies a different effective transverse width. It is not a
resolved curved crown or a contact area that grows with compression.

The scale STILL supplies the total felt area, thickness and hammer mass;
`--hammers` STILL supplies constitutive properties. Area is first divided among
the original unison members, then partitioned by positive quadrature weights.
Thus the force, elastic energy and Prony branch stiffness at a site all use
that site's physical area. Adding sites does not multiply total felt area or
hammer mass. Each duplex segment stays outside direct hammer contact.

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
Prony recovery. Those histories are independent across positions AND unison
members. Released sites continue recovering without erasing permanent crush.
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
```

Five bank tests cover input/quadrature, nodal contact, reciprocal compliance,
actual force/work closure, and original point parity. Five engine tests cover
area/Prony scaling, contact motion, independent history, una corda, re-strikes,
refusal and jack drive. Two frontend tests cover complete supplied inputs and
actual MIDI/CSV-to-stereo-to-PCM playback with one physical timeline. These
require native execution; independent matrix checks do not establish Rust
compilation, passing tests, measured realism or acoustic convergence.
