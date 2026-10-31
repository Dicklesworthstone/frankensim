# Played string extension and longitudinal waves

`piano_exterior render` and `render-loaded` accept
`--string-stretching strings.fsps`. The same per-key file selects linear
strings, averaged Kirchhoff--Carrier extension (`stretch`), or finite
longitudinal string dynamics (`waves`). Both nonlinear selections act through
the existing string/hammer/bridge mechanics in the same mechanical tick.

```sh
cargo run --release -p fs-couple --example piano_exterior -- \
  render-loaded settled.fss strings.csv body.obj acoustics.fspe piano.wav 6 \
  --string-stretching strings.fsps --modes 128 --substeps 8 \
  --hammer-footprints faces.fshp --dampers pads.fspd --midi score.mid
```

The optional controls compose with the same complete supplied hammer materials,
point/finite hammer faces, jack-force CSV, MIDI mapping, sustain, sostenuto,
una corda, spatial dampers and one-way or passive radiation feedback. Existing
inputs without this option retain the original linear-string image.

## Grand-piano preset and half-space playback

The `grand_piano` renderer accepts the same material format and solver through
`--string-stretching axial.fspx`. For example, with complete supplied scale and
material files:

```sh
cargo run --release -p fs-couple --example grand_piano -- \
  --preset steinway-d --scale strings.csv --string-stretching axial.fspx \
  --midi score.mid --dampers estimated --microphone-right 1,1,1 \
  --duration 6 --render piano.wav
```

The option requires `--render`, a nonempty file path and a record for every
selected scale key. `--note` never narrows the material coverage requirement.
The renderer admits the supplied material against its already-tuned scale and
attaches it before excitation, for every default/preset/imported hammer path.
Unreadable or invalid input refuses instead of silently reverting to linear
strings. The material path cannot also name a render or export output.

This composes with CSV and MIDI performances, finite hammer footprints, supplied
hammer cards, all pedals, spatial dampers and the existing mono/stereo Rayleigh
receivers. The render reports the selected speaking/duplex channel count and
retained longitudinal coordinate count. Zero nonlinear channels means explicit
all-linear selection. All-linear files preserve the original playback
arithmetic, including the Model D preset. Acoustic assumptions remain
those of `grand_piano`, not the finite-body exterior renderer above.

### One material section for transverse and axial mechanics

`grand_piano` also accepts `--string-sections bonded.fsss` together with
`--string-stretching axial.fspx`. Each bonded section supplies its original
layer geometry, density and Young modulus to `fs-material`. The same section
integrals supply linear mass density, bending rigidity EI and axial rigidity EA.
Select `section` in place of the numeric EA in each nonlinear row:

```text
frankensim-piano-string-stretching-v1
waves,69,section,0.2,4000,4,0,0.05
stretch,72,section,0.2
linear,74
```

Changing only a material layer's density or modulus now propagates through
transverse frequencies and mass normalization, nonlinear tension and stored
extension energy, longitudinal wave speed and mode retention, and the existing
reciprocal string/bridge forces. Installed tension remains fixed unless the
caller explicitly requests tuning. The slope and strain limits, longitudinal
band, mode budget and axial loss remain explicit model inputs; the layer data
does not identify those independently.

`section` refuses when no material section was resolved. A numeric EA alongside
a bonded section also refuses, preventing a second source of axial stiffness.
A cached stretching specification must be resolved again after changing its
section EA. Original scale CSVs remain valid for effective linear coefficients;
their twelve columns do not preserve the material binding, so retain and supply
the section file when reloading a material study.

These are solid or perfectly bonded concentric members. A helical winding with
slip requires its own physical model; this option does not treat a piano bass
string's winding as bonded. Material-specific measured damping from RT-0425
remains excluded for substituted sections.

Five focused `grand_piano` regressions cover CLI/source admission, malformed and
incomplete material files, changing actual tension with the original energy
ledger, bitwise all-linear playback, and nonlinear MIDI/CSV stereo equivalence
across different block sizes with finite hammers and spatial dampers. Run them
with the existing target; their existence is not a native passing-test claim:

```sh
cargo test -p fs-couple --example grand_piano string_stretching
```

## Complete per-key SI input

```text
frankensim-piano-string-stretching-v1
# Illustrative three-key study, NOT identified Steinway axial rigidities.
waves,69,150000,0.2,60000,16,0.02,0.05
stretch,72,150000,0.2
linear,74
```

The file must cover **every key in the supplied scale exactly once**, including
unplayed sympathetic courses. A full preset needs 88 rows. `linear,key` explicitly
retains the original linear potential for that key. A `stretch` row supplies
key, positive finite axial rigidity **EA in newtons**, and a finite continuous
slope limit in `(0,0.3]`. The example is valid only with a three-key 69/72/74 scale
whose geometry and selected mechanical rate admit the stated axial band.
Comments after `#` and blank lines are permitted. The file is bounded to 64 KiB.

The complete `waves` row has eight fields:

```text
waves,key,EA_N,maximum_slope,axial_band_hz,mode_cap,damping_ratio,maximum_strain
```

| Field | Meaning and admission |
| --- | --- |
| `key` | One supplied scale key, covered exactly once by `linear`, `stretch` or `waves`. |
| `EA_N` | Positive finite axial rigidity in newtons, or `section` when a bonded material section supplies EA. |
| `maximum_slope` | Bound on the complete transverse slope, in `(0,0.3]`. |
| `axial_band_hz` | Positive finite axial frequency ceiling, at most `0.45` times the mechanical sample rate. |
| `mode_cap` | Integer in `1..32`, applied to each physical speaking or duplex segment. A cap below the modes needed for the supplied band refuses. |
| `damping_ratio` | Supplied axial modal damping ratio in `[0,1]`; zero selects lossless axial oscillators. |
| `maximum_strain` | Positive additional material-strain bound at most `0.1`, applied to the absolute mean plus the sum of retained nonuniform strain amplitudes about the prestressed reference. |

For segment length `L` and supplied linear mass density `mu`, axial frequencies
are `f_n = n sqrt(EA/mu)/(2L)`. The adapter retains consecutive modes from `n=1`
through the selected band. It uses each segment's actual length, including duplex
length; `mode_cap` is an admission limit, not permission to drop modes inside the
band. The mechanical sample rate is `48000 * --substeps` in these renderers.
Selecting a larger band therefore requires both enough modes and enough explicit
time resolution. For the example above, `--substeps 8` admits the frequency
ceiling; the supplied scale still determines whether its mode cap is sufficient.

EA is neither installed tension T, bending rigidity EI, nor spring rate EA/L.
In particular, an effective winding mass and bending rigidity do not determine
the axial rigidity of a wrapped string. The program does not infer or invent
that missing material input. One course's supplied EA applies to all its unison
members and existing duplex segments; the geometry and individual detuning
remain those of the original scale.

On exterior commands, missing files, duplicate/extra/missing keys, nonfinite
values and invalid material limits refuse before structural eigenanalysis or
acoustic BEM preparation. Clock, per-segment retention and dynamic strain checks
remain part of physical preparation and execution.
An explicitly all-linear file is legal and retains the old numerical path.
`response` and `admittance` reject the option: their present linear harmonic
experiments cannot silently stand in for an amplitude-dependent response.

## Shared mechanics

The existing string owner supplies the extension channel. The piano adapter
includes both the fixed-interface string displacement and the moving bridge's
chord slope. Its work-conjugate forces act on the **same string and board
coordinates** as the existing contact solution. The nonlinear extension force
is resolved against all enabled hammer sites, their original felt/Prony
histories and the shared board before accepting one mechanical tick.

The unloaded small-amplitude frequencies, modal mass normalization, unison
layout, material loss and acoustic boundary projection are unchanged. During a
strike, geometric extension stores additional energy and changes physical
string tension. A `waves` selection also carries longitudinal displacement and
velocity and their supplied modal damping. The complete energy balance includes
the retained axial kinetic energy and strain potential once. No stored
energy is turned into an extra dissipative loss or an output pitch envelope.

Hammer/jack clocks, damper flow and acoustic splitting advance once per tick,
not once per nonlinear trial. The engine retains its original convergence,
material, energy and sample-rollback checks. A slope, strain or convergence refusal
does not clip the string or publish a weaker linear substitute. Raising
substeps is an explicit caller choice, never automatic clock repair.

Both exterior receivers observe the same nonlinear mechanical trajectory.
Radiation-loaded playback retains the fitted acoustic storage and reaction;
selecting extension neither changes the fit tolerances nor bypasses a failed
passive fit. The report identifies whether nonlinear channels are active and
the selected material path and number of retained longitudinal coordinates.

## Finite longitudinal dynamics

The `fs-nlmodal::longitudinal::LongitudinalString` owner resolves the
moderate-slope additional material strain about the prestressed reference

```text
epsilon = u_x + (y_x^2 + z_x^2)/2.
```

The supplied `maximum_strain` bounds this increment. It does not add the
installed reference strain `T0/EA`. Total tension is checked separately as
`T0 + EA * epsilon`, including its lower bound over the complete reduced span.

Here `u` is longitudinal displacement coupled at the bridge through the
speaking termination's **axial participation row** `a_n` (below), and `y` and `z`
are the retained transverse directions. A single transverse direction is also
supported. Both polarizations of one physical string share one axial state and
one strain energy; selecting two directions does not double the string. The
moving bridge supplies both the transverse chord and the axial boundary extension
`u(L, t) = sum_n a_n x_n(t)` over the board's loaded coordinates `x_n`, with
exact reciprocal discrete-gradient force pullback (the force on coordinate `n`
is `-a_n * EA * epsilon(L, t)`, plus the retained axial modes' share).

### The axial participation row

`a_n = t . (U_n + Theta_n x r)`: the termination's motion along the string axis
`t` (pointing away from the speaking length, so a positive value extends the
string), for board mode `n` with translation `U_n` and axial rotation `Theta_n`
at the attachment, and the arm `r` from the board's structural site to the string
contact. It replaces the former single authored scalar `BRIDGE_AXIAL_PROJECTION =
0.06`; the row is per mode, changes sign from mode to mode, and is 0 for a
board with a perfectly flat string and no arm. There is no default: a `waves` key
without its row refuses when the bank is configured (`axial_rows` of the
specification), and `grand_piano` derives the rows from the board file when
`--string-stretching` selects `waves` (`BoardGeometry::bridge_attachments`, then
`MotionSurface::project_at`, so it is the same map as the bridge coefficient):

- **The arm** is the height of the outer fibre of the stiffener that contains the
  bridge station and lies on the string side of the panel's midplane
  (eccentricity > 0), `e + h/2`, with `h = sqrt(12 I / A)` read from that row's
  area and bending inertia (a rectangular section is assumed).
- **The string axis** is NOT in the board file. It is taken perpendicular to the
  bridge in plan, with the speaking length toward decreasing `y` (the RT-0425
  frame's front edge), so `t` has `t_y >= 0`. This is an assumption; see
  `u-string-direction` in `data/instruments/descriptions/piano/unknowns.tsv`.
- **Downbearing** is absent from a flat plate, so the axis has no vertical
  component: the `sin(alpha) U_n` term is exactly 0.
- With `--board-reduction` the rotations are those of the reduced basis, built
  from vertical bridge-force snapshots only.

`piano_exterior` and the live engine still use a labelled authored ratio
(`EXTERIOR_AXIAL_RATIO`, `LIVE_BRIDGE_AXIAL_RATIO`, both 0.06 times the vertical
row) until they carry rotation rows of their own.

The retained axial sine modes have inertia, propagate their strain at the
material wave speed, and exchange energy with transverse and bridge motion
during contact and free vibration. Their damping acts on axial velocity. Higher
unretained zero-mean strain channels are statically equilibrated: they have no
resolved inertia or propagation. The mean extension remains present. An empty
retained axial set therefore has the averaged Kirchhoff--Carrier potential; a
small axial band does not establish a complete longitudinal wave solution.

No independent force lag or sample clock is introduced. Trial axial states are
accepted with the same converged string/contact solve, and a refused frame
restores their state together with the original mechanics. The additional-strain
bound is checked across the complete reduced span. Loss of positive total
tension refuses the trial rather than entering an unmodeled slack-string regime.

This is a reduced moderate-slope model with an effective axial-to-transverse
bridge port (`u(L, t) = eta * b(t)`). Both transverse slope forces and axial
tension fluctuations at `x = L` act reciprocally on the bridge coordinates and
close the same mechanical energy ledger. It remains simpler
than the geometrically nonlinear Timoshenko string model in Chabassier,
Chaigne and Joly, [*Time domain simulation of a piano. Part 1* (2014)](https://archive.numdam.org/articles/10.1051/m2an/2013136/),
which motivates resolving longitudinal waves in piano modeling; this feature
does not claim to reproduce that full model.

## Boundaries

`stretch` selects averaged moderate-slope extension; `waves` adds the explicitly
retained longitudinal modes and the axial participation row described above. With `--string-polarization`, both
transverse directions share their total strain and one extension energy per
physical segment; see
`STRING_POLARIZATION.md` for the required geometric input. It does not resolve
independent 3-D rigid-body DOFs of the bridge cap, tension-dependent winding slip, full 3-D hammer/string
contact, or measured piano string materials. The acoustic fit
remains bounded to its admitted sampled band even when a nonlinear attack
contains higher-frequency energy. Retained modes, substeps and a supplied EA
are not a realism, convergence or real-time certificate.

The exterior regression target exercises strict complete-scale admission,
bitwise all-linear jack/pedal behavior, and actual one-way and radiation-loaded
stereo against direct calls to the same mechanical owner. It requires nonzero
extension, changed physical motion, felt/damper loss and the original combined
energy bound; native execution status is reported separately from these tests'
existence.

The longitudinal renderer regressions exercise finite axial state through actual
Rayleigh and finite-body BEM pressure, MIDI/CSV and block-size equivalence,
finite hammer faces, spatial dampers, and the same exterior radiation feedback.
They compare the selected wave model with the same supplied-EA averaged model;
they do not serve as a mesh, axial-band or measured-piano convergence study.

```sh
cargo test -p fs-couple --example grand_piano longitudinal_render_tests
cargo test -p fs-couple --example piano_exterior playback::longitudinal_tests
```
