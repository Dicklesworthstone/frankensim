# Piano damper contact and spatial damping

`--dampers pads.fspd` selects either **finite-gap wool-felt contact** or
**spatial viscous damping**, according to the file header. Both act on the
existing speaking strings and loaded soundboard. `--dampers estimated`
explicitly selects approximate viscous pads. Omitting the option keeps the
original point-drag image.

Felt contact adds physical separation, finite pad travel, compression, crush
and time-dependent relaxation. It requires supplied pad geometry and material
cards; there is no estimated felt selection. Its actuator prescribes pad
height and accounts for the work done against the strings. It does not include
pad mass, gravity, preload or a full keyboard linkage.

```bash
cargo run --release -p fs-couple --example grand_piano -- \
  --preset steinway-d --dampers estimated --note 69 \
  --duration 6 --render release.wav

cargo run --release -p fs-couple --example grand_piano -- \
  --preset steinway-d --dampers pads.fspd --midi performance.mid \
  --duration 30 --render performance.wav
```

The same `--dampers` input also works with `piano_exterior render` and
`render-loaded`. Harmonic response/admittance modes do not admit played
nonlinear damper contacts. Input selection happens before rendering or writing
output, and the existing direct input/output path collision checks apply.

## Finite-gap felt pads

The header is `frankensim-piano-felt-dampers-v1`. Every key in the selected
string scale needs exactly one `pad` or `free` row. Every pad also needs its
own explicit `felt` material card; optional `branch` rows follow that card.

```text
# Illustrative supplied values for a TWO-KEY scale only.
# These values are not identified Steinway damper parameters.
frankensim-piano-felt-dampers-v1
pad,69,0.12,0.14,0.00001,0.004,0,0.001,0.05,1,0,2
felt,69,60000,0.2,2.5,3.2,0.1,0.8,300000
branch,69,200000,0.002
free,90
```

The geometry record is:

```text
pad,key,start_m,end_m,area_m2,thickness_m,closed_gap_m,lift_m,speed_m_s,normal_u,normal_v,stations
```

| Field | Physical meaning |
| --- | --- |
| `start_m`, `end_m` | Pad endpoints measured from the fixed end of the speaking string; require `0 < start < end < length`. |
| `area_m2` | Total positive contact area **per speaking string**, divided across that pad's sites. It is not an area per metre or divided by unison count. |
| `thickness_m` | Positive unloaded felt thickness at every site of this pad. |
| `closed_gap_m` | Nonnegative separation between the rest string and the closed, unloaded pad. Zero permits contact when the string moves toward the pad. |
| `lift_m` | Positive additional separation at fully raised travel. |
| `speed_m_s` | Positive magnitude of prescribed pad motion toward its target height. |
| `normal_u`, `normal_v` | Supplied unit normal in the string's primary/secondary transverse frame, pointing toward the pad. A nonzero secondary component requires two-plane string geometry. |
| `stations` | One to four equal-area midpoint contact sites along the span. At most twelve sites are admitted per three-string unison course. |

Each site has independent wool-felt crush and Prony memory. The `felt` and
`branch` rows have the same units and parameter order as
[HAMMERS.md](HAMMERS.md):

```text
felt,key,reference_stress_Pa,reference_strain,loading_exponent,unloading_exponent,crush_fraction,densification_strain,equilibrium_modulus_Pa
branch,key,modulus_Pa,relaxation_time_s
```

The constitutive owners are the existing `fs-material::WoolFelt` and
`GeneralizedMaxwell` creep realization. No second felt law is introduced.
Zero to eight Prony branches per pad material are allowed. Material validation,
including positive moduli/times and valid felt strain limits, remains with
those owners. A `free` course has no pad and must not be given an unused felt
card. Missing, duplicate, unknown, nonfinite or unsupported inputs refuse
without selecting a viscous fallback. Blank lines and `#` comments, including
comments before or after the header, are accepted. Reads are bounded to 1 MiB.

### Pedal travel, contact and work

Let `s` be sustain travel in `[0,1]`. An unheld, unlatched course targets
`closed_gap + lift*s`. A held key or a key captured by sostenuto targets
`closed_gap + lift`. The prescribed height moves toward that target by at most
`speed*dt` per mechanics substep. Contact depends on the actual local string
displacement relative to that moving height. Partial sustain therefore changes
when the vibrating string meets the felt; it does not scale a drag coefficient.

The contact force is unilateral. A string may separate from a site, return,
and meet the site's retained crushed/relaxed state. The local string port
includes the moving bridge and the supplied transverse normal. Each speaking
member of a unison is loaded independently; duplex segments have no direct pad
contact and interact through the shared bridge.

Felt forces enter the same coupled bank/hammer/string solve at the existing
mechanics rate. Contact iterations use candidate state and do not advance pad
travel or material memory repeatedly. A failed public audio frame restores pad
heights, contact history and work/loss accounting along with strings, hammers
and acoustic state. Configuring the pads is a cold operation before excitation;
an active contact history cannot be replaced by another damper image.

For a positive compressive contact force `F` and height increment `delta_h`,
the actuator work is `-F*delta_h`. Closing against a string can add energy;
lifting a loaded pad can remove it. The signed work is included in
`accounting.input_work_j` and exposed separately as
`accounting.damper_actuator_work_j`. It must not be added to the input a second
time. Pad storage is included in `Instrument::energy_j()` and can be inspected
with `damper_energy_j()`. `damper_relaxation_loss_j` is a subset of
`damper_loss_j`, which is already included in total dissipated energy.

Both renderers report pad/site counts, damper storage, loss, relaxation loss
and signed actuator work. Microphones consume the resulting physical board
motion. There is no extra output envelope, point drag or spatial viscous pad
stacked on the felt contact model. The contact-site cap is a work budget, not
a spatial convergence bound. It does not establish measured damper decay,
stick-slip, acoustic calibration, real-time performance or full piano action.

## Spatial viscous pads

The viscous image remains available with `frankensim-piano-dampers-v1` or the
explicit `estimated` selection. It has no finite gap or moving pad state.

The original point at 35% of string length is blind to partials with a node
there, including the 20th pinned-string partial. A finite pad instead integrates
local velocity **squared** across its span. Opposite-moving portions do not
cancel their dissipation, as they would when damping only the average velocity.
The full modal cross terms and moving-bridge displacement contribution remain.
Each speaking member of a unison is damped independently; duplex segments have
no direct pad force but continue to interact through the shared bridge.

Key hold, sostenuto capture and sustain use the existing engine controls,
including MIDI and SI CSV. A held or latched key lifts its pads. Sustain travel
retains the existing `(1-travel)^2` drag scaling; full sustain removes pad drag.
The spatial model replaces the point model.
Explicit pad/free rows own the upper break and individual drag coefficients.
There is no output envelope, resonance reset, tail cutoff or pressure gain.
The original contact solver, loss accounting and failed-sample rollback remain
in use. Microphones receive the changed physical board velocities.

## Silent key depression and selective sympathetic resonance

An SI performance CSV can hold chosen keys open without striking them. Use
`silent_key_down,key,0`, then the ordinary `note_off,key,0` to release each key.
The same control is available to hosts as `Instrument::silent_key_down(midi)`
and to prepared schedules as `Control::SilentKeyDown { key }`.

For example, save this 48 kHz schedule as `sympathetic.csv` and render a scale
that contains keys 48 and 60:

```csv
sample,event,key,value
0,silent_key_down,60,0
0,sostenuto,0,1
4800,note_off,60,0
9600,note_on,48,2
33600,note_off,48,0
96000,sustain,0,0.5
96000,sostenuto,0,0
144000,sustain,0,0
```

```bash
cargo run --release -p fs-couple --example grand_piano -- \
  --preset steinway-d --dampers estimated --performance sympathetic.csv \
  --duration 4 --render sympathetic.wav
```

Key 60 is held silently before the rising sostenuto edge and remains captured
after its key release. The later strike of key 48 supplies the energy; motion
reaches key 60 through the existing reciprocal string/bridge/soundboard model.
At two seconds the sostenuto release engages the half-pedal drag, and at three
seconds sustain returns to zero. Equal-sample rows execute in file order: a
silent depression after an already-down sostenuto pedal is not captured.

Silent depression changes only the existing key hold. It preserves hammer
position and velocity, felt crush and relaxation memory, and all string/board
motion. Repeated silent depression is harmless, including on a ringing key;
release it before requesting a new strike, under the usual rearming rules.
It works with point, spatial viscous and finite-gap felt dampers. For felt
pads it sets the raised target, which is reached at the supplied travel speed.
An undamped upper key has no pad to lift. There is no change to MIDI:
zero-velocity note-on remains note-off, and CSV `note_on` still requires
positive SI hammer velocity.

The point and spatial viscous images model no work for the key control itself.
The finite-gap felt image accounts for subsequent pad actuator work whenever
its finite-speed motion meets a loaded contact. None includes the complete
mechanical work or inertia of a silent key depression.

## Viscous per-key input in SI units

The header is `frankensim-piano-dampers-v1`. Every key in the supplied string
scale must have exactly one row, even when the performance never strikes it:

```text
frankensim-piano-dampers-v1
# Example ONLY for a two-key scale, with key 69 longer than 0.14 m.
pad,69,0.12,0.14,0.4
free,90
```

`pad,key,start_m,end_m,drag_ns_m` specifies endpoints measured from the fixed
end of the speaking string. Require `0 < start < end < string length` and
finite positive drag up to 1e6 N s/m. Drag is the TOTAL coefficient per speaking
string across the pad, not a coefficient per metre and not divided by unison
count. The local dissipation is `c/(end-start) * integral(v(x)^2 dx)`.
`free,key` explicitly declares no damper on that course. A full preset therefore
needs 88 rows; the two-row example is not a complete full-keyboard specification.
Blank lines and `#` comments are accepted. Duplicate, missing, unknown, invalid
or out-of-span rows refuse without choosing defaults. Reads are bounded to 1 MiB.
For either input format, a file literally named `estimated` can be selected as
`./estimated`; the bare word is reserved for the following estimate.

`estimated` uses a centre at 35% of speaking length, width equal to the smaller
of 40 mm and 8% of length, and 0.4 N s/m per string. Keys above MIDI 88 are free.
These are explicit estimates, NOT Steinway dimensions or measured felt losses.
Supplied geometry, tensions, hammer cards and acoustic observer are unchanged.

## Viscous numerical scope

Positive midpoint quadrature retains at least eight stations per footprint,
with additional stations tied to the highest retained spatial partial. The
limits are 128 stations per string and four million projection terms; they
refuse rather than dropping modes. This guard is not a spatial error bound.
The entire shared-board port sequence is composed forward then backward using
exact dissipative rank-one velocity flows. Point dampers use the existing
half-step before mechanics and its reversed sequence afterward; spatial pads
compose the symmetric sequence inside each damping half-step. Each string
receives the original total drag interval. This gives a passive, second-order split
approximation, not the exact full coupled damping exponential. Stored
energy does not change when installing a viscous specification. Preparation is
cold; stepping adds no allocation in this damping kernel.

The render reports actual pad/station counts and damper energy loss, already
included in total losses. This viscous image makes no falling-pad, preload,
collision or hysteretic contact claim. It also makes no real-time, full-band
or measured-decay claim.

## Regression scope

Existing viscous regressions cover nodal damping, moving-bridge geometry,
energy, pedals, duplex behavior, contact-failure retry and MIDI/CSV render
preparation. The finite-gap renderer regressions are G3 checks of actual
Rayleigh/exterior BEM pressure, contact loss, signed work accounting, direct
owner versus file/inline preparation, and identical playback across block
boundaries. These are executable numerical checks; source presence alone is
not a passing result or an audible/calibrated piano validation.
