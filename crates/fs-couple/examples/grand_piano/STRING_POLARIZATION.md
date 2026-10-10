# Two transverse string directions in played piano

`grand_piano` and `piano_exterior render` / `render-loaded` accept
`--string-polarization bridge-frames.fspp`. This selects the existing vector
string mechanics for every speaking and duplex segment. It supplies the missing
geometric connection from the played board's full motion to the strings' lateral
direction. Without the option, both frontdoors retain their original one-plane
mechanics.

```sh
cargo run --release -p fs-couple --example grand_piano -- \
  --preset steinway-d --board-geometry supplied.fsb --scale strings.csv \
  --string-polarization bridge-frames.fspp --render piano.wav \
  --midi performance.mid --midi-half-pedal --dampers pads.fspd

cargo run --release -p fs-couple --example piano_exterior -- \
  render-loaded supplied.fss strings.csv body.obj acoustics.fspe piano.wav 6 \
  --string-polarization bridge-frames.fspp --midi performance.mid
```

The files must describe the same board and complete scale. The preset may also
supply the board, but its lateral string/bridge frames are still explicit input;
no measured Model D bridge heights, axes or lateral drag are inferred. A modal
CSV or authored demonstration mode list contains no full-vector geometric
motion and cannot supply this capability. Harmonic `response` and `admittance`
commands reject this playback option.

## Physical input

The UTF-8 file is limited to 64 KiB. Blank lines and `#` comments are allowed.
It starts with this header and one attribution row, followed by one `course`
row for **every key in the selected scale**, including keys absent from a score:

```text
frankensim-piano-string-polarization-v1
source,estimated,explicit example frame for a supplied test panel
# course,key,triangle,w0,w1,w2,arm_x_m,arm_y_m,arm_z_m,string_x,string_y,string_z,hammer_x,hammer_y,hammer_z,lateral_damper_ratio
course,69,0,0,0,1,0,0,0.02,0,1,0,0,0,1,0.35
```

This is an **authored example**, not a source of measured piano parameters. It
is complete only for a one-course A4 scale. The example chooses the third vertex
of triangle 0, a 20 mm bridge arm along global +z, a string along +y, a hammer
normal along +z, and a lateral viscous drag 0.35 times the primary drag.

| Field | Meaning |
| --- | --- |
| `source` | `estimated`, `mixed`, `published`, or `measured`, plus a nonempty attribution without commas. This records the caller's claim and does not certify it. |
| `key` | Piano MIDI key 21–108, exactly once for every selected course. |
| `triangle` | Zero-based triangle index in the supplied or generated structural board mesh. |
| `w0,w1,w2` | Barycentric location in that triangle's stored node order. Each lies in [0,1]; their sum is one. |
| `arm_x_m,arm_y_m,arm_z_m` | Supplied vector from that structural interpolation site to the string's bridge attachment, in global metres. Zero explicitly means no arm. |
| `string_x,string_y,string_z` | Global unit vector along the string. |
| `hammer_x,hammer_y,hammer_z` | Global unit vector in the existing hammer contact direction, perpendicular to the string. |
| `lateral_damper_ratio` | Finite value in [0,10] multiplying the primary damper drag for the other transverse direction. Zero explicitly disables lateral pad drag. |

All numeric values must be finite. Axes are checked for unit length and
orthogonality; they are not silently normalized. Triangle indices, sites and
frames must be updated when the structural mesh changes. Missing or duplicate
keys, invalid frames and an unavailable motion source refuse before playback.

## How the geometry enters the mechanics

The board preparation retains translations and **physical axial rotations**
from the same eigensolve that supplies its primary bridge coefficients. For
each bare-board mode, the bridge motion is

\[
u_b = \sum_{i=0}^{2} w_i\bigl(u_i + \theta_i \times a\bigr).
\]

The primary coefficient is the projection onto the supplied hammer axis. It
must reproduce the board's existing bridge coefficient at that key, in every
retained mode, to numerical roundoff. This prevents the new file from silently
moving or reorienting the original hammer-plane coupling. The lateral axis is
`string_axis × hammer_axis`, and its coefficient is its dot product with the
same bridge motion. These rows enter the existing loaded mass coordinates once.
Both string directions exert reciprocal forces on the same board.

A flat DKT board has transverse nodal translation and two physical rotations.
An explicitly supplied bridge arm can therefore provide lateral motion through
bridge rocking. At its midplane, a zero arm can legitimately give a zero
lateral row. The code preserves that zero; it does not add an artificial
coupling to make the result audible. A crowned shell additionally supplies
in-plane nodal translation.

Flat P1 motion supports the original inertia, mass equilibration and consistent
P1 transverse inertia choices. `--edge-cubic-board-mass` instead retains the
existing cubic displacement and its analytic physical rotations
`[dw/dy,-dw/dx,0]`. Bridge frames and acoustic skins evaluate that same field;
they do not interpolate the cubic solution with P1 motion. The primary bridge
compatibility check still applies at every mode, including interior sites.
Both flat choices are available in `piano_exterior` response, admittance,
render and render-loaded. Crowned shells retain their existing full-vector
field and refuse flat-board mass selections. Rayleigh acoustic refinement
does not change the structural interpolation sites or modes and currently
requires P1 motion.

## Played behavior and limits

The original hammer/felt solver drives the primary direction only. Each key
still has one hammer and its original total contact area; finite hammer faces
keep their existing independent contact memories. Source or supplied hammer
materials, shank mechanics, repeated strikes and una corda continue through the
same engine. The lateral response comes from physical board/string coupling.

Key release, sustain travel, half pedal and sostenuto apply the explicitly
supplied lateral drag through the existing point or spatial damper flow.
Duplex segments remain undamped by those pads. When `--string-stretching` is
also selected, both transverse directions share their total geometric strain,
one physical tension and one extension energy per segment.

Mono and stereo receivers observe a single mechanical trajectory. The lateral
coordinates and optional radiation reaction participate in the same mechanical
substeps, energy accounting and rollback as the primary coordinates. No second
audio clock, output effect, detuning law or new oscillator is introduced.

This remains a reduced string/board model. The file supplies one frame per
course, shared by that course's unison and duplex segments. It does not resolve
3-D hammer contact, string torsion, longitudinal waves, bridge material
deformation, falling felt pads or a complete piano action. A supplied frame
and a larger mode budget are not a convergence, real-time or perceptual
similarity certificate.
