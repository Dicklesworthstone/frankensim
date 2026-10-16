# Simultaneous finite-hammer contact

A finite hammer face can load up to twelve separate felt sites: four positions
on each of three unison strings. All those sites act on one hammer inertia,
while nearby sites also share string displacement. Iterating each site against
frozen forces at the others can stall even when the physical contact equation
has a well-resolved admissible solution. A failed solver is not a soft strike,
and increasing felt compliance to make it converge would change the instrument.

Finite faces now use a simultaneous **per-hammer contact block** inside the
existing global contact iteration. No new command-line option is necessary:
`--hammer-footprints` selects the physical geometry as before. Complete point
selections and unselected hammers retain the original scalar arithmetic.
The same engine is used by `grand_piano`, imported-piano playback, and both
one-way and loaded exterior playback wherever finite faces are selected.

## The same physical equation

After the existing string/board/hammer response has been condensed, write the
elastic end overlaps as `e = free - A*f`. Here `A` includes the full signed
mechanical compliance and each site's series Prony compliance on the diagonal;
`free` includes the old Prony strain propagated by the existing material owner.
Each `f_i(e_i)` is the original discrete-average WoolFelt force with that site's
own area, thickness, and committed loading/unloading history.

The block solves `e + A*f(e) - free = 0`. Its Jacobian is
`I + A*diag(df/de)`, using the derivative of the **discrete-average force**, not
an endpoint tangent or a linearized replacement material. The existing
`fs-la::LuWorkspace` solves this at most 12-by-12 system. Overlaps can be negative,
so opening/closing sites need no tensile-force clamp or merged active contact.
Trial line search stays below the original densification bound. A root must
also pass the original force-equation check at the actual mechanical overlap.

Every hammer block updates all other active contacts through the same full
soundboard compliance. Cross-key coupling is not discarded or treated as
independent voices. The original 32 global sweeps and final simultaneous force
check remain; each block has at most 32 Newton updates and 24 line-search trials
per update. Singular, nonfinite, outside-domain or exhausted solves refuse.
This improves one important convergence limitation, not every possible stiff
multi-key performance or an underresolved physical timestep.

## Preserved state, work, and observations

All block iterations modify numerical scratch only. String, board, hammer,
jack, felt-conditioning, Prony, and acoustic states still advance at their
original accepted-sample boundary. Failure, including a later hammer block,
restores the complete previous physical state and accounting. A retry cannot
inherit a failed Newton iterate. LU is allocated at cold construction only;
there is no allocation in the new contact solve or group assembly.

The final forces still enter the original exact-ZOH/Schur mechanical update.
The same felt force performs the same displacement work; loss and recoverable
energy are not recomputed from a different endpoint law. Force, energy,
densification, stiffness, mass, pressure scale and output rate are unchanged.
Whole-string una corda removes all sites of the excluded string; remaining
areas are not renormalized. MIDI/CSV score timing and stereo's single substep
trace are unchanged. No attack envelope, synthetic sound, added damping or
radiation channel is introduced.

## Direct checks and limits

A stiff twelve-site coupon with one shared hammer and local string compliance
leaves an approximately 0.00083 N force-equation defect after the original 32
scalar sweeps. The simultaneous equation has equal forces near 13.04413 N and
converges in six Newton updates in an independent arithmetic transcription.
Those are explicitly authored coupon parameters, not calibrated piano felt.
Native tests also manufacture a root using an actual loaded string-bank matrix,
without changing its source geometry or averaging away the finite face.

Additional native checks compare converged scalar and block trajectories for a
two-note chord, retain distinct loading/unloading memory and signed cross terms,
and exercise whole-string pedal selection and exact retry with acoustic history.
Point-only configurations retain a bitwise comparison against the old image.
The physical finite-face and acoustic limitations in `HAMMER_FOOTPRINTS.md` and
`RADIATION_FEEDBACK.md` still apply; no real-time or specimen-fidelity claim is
made. This solver improvement does not establish modal or temporal accuracy.

```sh
cargo test --release -p fs-couple --example grand_piano felt::block -- --test-threads=1
cargo test --release -p fs-couple --example grand_piano contact_solver -- --test-threads=1
```

Independent Python arithmetic is not native Rust execution. The new authored
tests require native execution before their pass status can be asserted.
