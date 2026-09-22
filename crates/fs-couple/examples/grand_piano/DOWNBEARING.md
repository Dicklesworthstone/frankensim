# Static downbearing changes the soundboard that is played

A crowned soundboard is not automatically a loaded soundboard. The existing
crowned-file path now accepts explicit per-course downbearing, solves the full
supported shell equilibrium, and uses its **consistent tangent stiffness** for
modal analysis. String/felt/shank/bridge dynamics and the pressure observer then
consume those modes. No output filter or prescribed frequency shift is added.

## Supply an unloaded reference and the loads

In a `frankensim-crowned-board-si-v1` file, retain the normal geometry, physical
sections, ribs, clamped supports and bridge rows. Add:

```text
preload-reference,unloaded
downbearing-source,estimated,Explicit study loads; not factory measurements
downbearing,69,10
```

The last row is a **10 N total downward load for the entire key-69 course**, not
10 N per string. It is an illustrative study value, not a Steinway measurement.
There must be exactly one `downbearing` row for **every bridge station in the
geometry**, including unplayed keys. Use an explicit zero for an unloaded
course. Missing, duplicate, extra, negative and nonfinite rows refuse. A source
row and the explicit unloaded-reference declaration are mandatory. Omitting all
three row types preserves the original unloaded model.

Alternatively, `bridge-load,key,fx_n,fy_n,fz_n` supplies a full signed Cartesian
course resultant (negative z is downward). It replaces the shorthand for that
key; declaring both is a duplicate error. Horizontal forces are not discarded.

The force acts at the station's actual barycentric position. `bridge_arm`
produces the corresponding moment by `r cross F`; the interpolation is the
transpose of the existing displacement/rotation-to-bridge-motion map. A vertical
arm alone creates no moment under a vertical force. Loads on fixed coordinates
are carried by support reactions. No rib, support or bearing location is guessed.

Render the file through the ordinary instrument path, with source Model D
strings and hammer mechanics retained:

```sh
cargo run --release -p fs-couple --example grand_piano -- \
  --preset steinway-d --board-geometry unloaded-with-loads.fsc \
  --note 69 --duration 2 --render loaded-board.wav
```

MIDI, supplied felt and damper cards, tuning and microphone controls still use
that same path. All bridge loads are static inputs: changing string tensions or
bearing angles requires recomputing their load rows. Importing an already-loaded
measured shape and labelling it unloaded is not made physically correct by the
parser. Manufacturing residual stresses cannot be inferred from surface shape.

## Derive loads from string support geometry

`load-strings` evaluates every string member's tension against the actual two
outgoing string directions. Supply the existing scale CSV and a point file:

```text
frankensim-bearing-points-si-v1
source,estimated,Example schema; replace support coordinates with supplied geometry
reference,unloaded
bearing,69,0,0.0,0.5,0.0,1.0,0.5,0.0
```

The row is `bearing,key,member,front_x,front_y,front_z,rear_x,rear_y,rear_z`.
Coordinates are **metres in the board's Cartesian chart**; member is zero-based.
The two points are the next string supports on either side of the bridge
bearing, not necessarily the remote tuning/hitch pins. The bridge point itself
comes from the supplied barycentric station plus its `bridge_arm`. The numerical
example above is illustrative, not a Steinway support measurement. Every
member of every scale course must have a row, and the scale must cover exactly
the board's bridge keys. There is no default bearing angle or missing-string
fallback.

For each member, the force is `T * (unit(front - bridge) + unit(rear - bridge))`.
Member tension uses the **same** `Course::tension_at_cents` and unison offsets
as the mechanical bank, leaving flexural rigidity unchanged. Full resultants
are summed by course, so asymmetric supports also produce horizontal bridge
force and lever-arm torque. No additional duplex branch is counted: the two
outgoing directions already represent the two sides of that bearing.

```sh
cargo run --release -p fs-couple --example piano_board_import -- \
  load-strings unloaded.fss strings.csv supports.fsbp with-loads.fss
cargo run --release -p fs-couple --example grand_piano -- \
  --preset steinway-d --board-geometry with-loads.fss --scale strings.csv \
  --note 69 --duration 2 --render piano.wav
```

Use the same scale for derivation and playback. An explicit scale without
`--concert-pitch` preserves its supplied tensions. `load-strings` only prepares
force rows; equilibrium is solved later at rendering/modal preparation. It
refuses already-loaded inputs, absent/duplicate members, nonfinite points,
coincident supports, oversized inputs and output overwrites.

Forces are computed at the supplied **unloaded reference**, then held as dead
loads. The utility does not iterate anchor directions or string axial stretch
with the board deflection, infer bridge-pin friction, replace speaking lengths
from the point file, or independently verify that source meshes and scale cards
came from the same instrument. These remain explicit model/data assumptions.

## Physical model and boundaries

The structural owner is `fs_plate::shell::preload`. It uses full-coordinate
Green–Lagrange P1 membrane strains with the existing linear DKT bending,
relative-spin stabilization and bonded eccentric beam operators. In-plane
static degrees of freedom are **not** truncated to the audible mode set. The
linear membrane energy is replaced by its nonlinear counterpart, not counted
twice. The reference mass matrix is retained.

The adapter injects `fs_solver::NewtonKrylovState` with line search. Eight equal
load increments each admit at most 32 Newton iterations; the entire solve admits
20,000 residual/Jacobian evaluations. Every accepted increment independently
checks physical force balance and positive tangent inertia with the existing
sparse LDLT owner. Maximum displacement gradients and nodal rotations are 0.1
and 0.1 rad, respectively. Loss of stability, excessive deformation or exhausted
work refuses instead of inserting artificial stiffness. The final residual,
static stored energy, maximum translation and solve work are reported in the
prepared-board description.

This is **small vibration about a solved dead-load equilibrium**, not nonlinear
soundboard dynamics. Beam laws/arms remain linear in their reference frames;
beam buckling, follower-string forces, glue slip/creep, manufacturing stress and
large-rotation bending remain outside the model. A complete high-load
continuation or buckling-branch calculation is not claimed.

The acoustic observer uses the equilibrium surface's area, normal and XY
position, but still projects its source to an infinite flat baffle. It is not
lid/cabinet/room scattering or acoustic backreaction. Native eigenvalue residuals
do not certify spatial convergence or measured piano realism.

Reference for the distinction between crown, downbearing equilibrium and
vibration about a prestressed configuration: Mamou-Mani, Frelat and Besnainou,
*Numerical simulation of a piano soundboard under downbearing*, JASA 123(4),
2401–2406 (2008), DOI 10.1121/1.2836787. The implementation is not a reproduction
of their measured piano or their full finite-element study.
