# Settle the soundboard at the supplied final tuning

`settle-strings` closes the static geometry/load loop. Unlike `load-strings`,
which evaluates directions on the undeformed board and freezes them, it solves
string directions and the nonlinear shell equilibrium together. Every unison
member uses its own supplied scale tension, including the existing detuning
calculation. Horizontal force and the bearing-arm moment are retained.

```sh
cargo run --release -p fs-couple --example piano_board_import -- \
  settle-strings unloaded.fss strings.csv supports.fsbp settled.fss
cargo run --release -p fs-couple --example grand_piano -- \
  --preset steinway-d --board-geometry settled.fss --scale strings.csv \
  --note 69 --duration 2 --render piano.wav
```

Use the same complete scale for settling and playback, without retuning it via
`--concert-pitch`. Support points use the existing `frankensim-bearing-points-si-v1`
format in `DOWNBEARING.md`: two points per member in the board's Cartesian SI
frame, attributed source, and `reference,unloaded`. Every geometric bridge key
and every member is required. No anchor, angle or material is inferred. Output
must be a new path; inputs are never changed.

## The tension boundary condition is explicit

This command represents a board at **final tuned tensions**: a tuner has brought
each member to the scale's prescribed tension while the board settles. The
point-to-anchor directions change with displacement, but the target tensions
remain fixed. This is not the transient response of fixed-rest-length strings,
nor a claim to calculate how an untended piano detunes as the board relaxes.

The underlying `fs_plate::shell::preload::equilibrate_tethered_shell` also accepts
positive axial stiffness for fixed-rest-length studies. Its tensile law is
`T = max(0, T_ref + k * (length - length_ref))`; slack spans do not push. Zero
stiffness explicitly selects the prescribed-tension boundary used here. Force,
material stiffness and direction-change stiffness derive from one potential.
The existing Newton/Krylov and sparse tangent-stability checks own the solve.
These massless connections do not replace dynamic strings.

## Preserve prestress without counting the strings twice

The coupled static solve returns a tangent that includes its tethers. The piano
must NOT bake that matrix into its board modes and then add the same string
endpoint stiffness again in the existing mechanical bank.

Instead, this command exports **converged force rows on the original unloaded
reference mesh**. It then runs the actual bare-board preparation under those
forces and checks that it recovers the coupled solution's loaded surface. A
mismatch, unstable bare board, excessive deformation or invalid acoustic graph
refuses the export. Changing the reference nodes to their loaded positions and
pretending that shape is stress-free would lose the prestress; this command does
not do that. Its summary reports physical force residual, maximum translation,
and maximum outgoing-span length change.

Playback therefore receives the shell-only tangent about this equilibrium;
the existing moving-boundary string bank supplies its own dynamic stiffness,
mass and reciprocal bridge coupling exactly once. The acoustic observer uses
the loaded surface through the existing projected flat-baffle path.

## Remaining approximation

The dynamic bank still uses its current transverse moving-endpoint string
model, not a full 3-D string/axial system linearized about every supplied anchor.
The support-point lengths do not replace the supplied speaking/duplex lengths.
Changing scale tensions or support geometry requires settling again. Beam
laws and bearing rotation offsets remain linear; no bridge-pin slip/friction,
manufacturing residual stress, full beam buckling, cabinet/lid scattering or
room acoustics is introduced. Numerical equilibrium is not validation against
a measured Steinway specimen.
