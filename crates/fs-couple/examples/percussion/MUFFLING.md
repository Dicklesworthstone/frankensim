# Spatial drum and cymbal mufflers

`--muffler SURFACE X_m Y_m C_Ns/m` adds an ideal fixed, bilateral viscous
attachment to the actual vibrating surface. Repeat it for several attachment
points. `batter` and `resonant` select the two drumheads; `shell` selects a
`splash` cymbal. Coordinates are in the reference XY plane, in metres, and
resistance is nonnegative N s/m. All attachments act along the same vertical
axis as the existing transverse playing coordinates. The following resistances
are illustrative inputs, **not measured finger, gel or tape material data**.

```sh
# Muffle both heads of the prepared drum, retaining its contact/air mechanics.
cargo run --release -p fs-couple --example percussion -- \
  drum-modal 4096 --strike-speed-m-s 4 --strike-position-m 0.06 0.01 \
  --muffler batter 0.12 0.01 0.8 --muffler resonant -0.10 0.02 0.4 > muffled.csv

# Same controls work before BEM pressure rendering, never as a PCM envelope.
cargo run --release -p fs-couple --example percussion -- \
  drum-modal-mic 4800 20 --muffler batter 0.12 0.01 0.8 > muffled.wav

# Hold a viscous attachment on the nonlinear curved shell, preserving stand felt.
cargo run --release -p fs-couple --example percussion -- \
  splash 4096 --prepared-nonlinear --muffler shell 0.075 0 0.5 > held-shell.csv
```

The controls also compose with supplied `--drum-spec` / `--shell-profile`,
nonlinear `drum-stretch`, snares, two independent sticks and their force files,
and the existing distributed cavity. Existing unsupported combinations remain
unsupported: for example, vented-cavity exterior audio is still refused.

## Mechanical meaning

For the **original, unnormalized** signed shape row `b`, local physical velocity
is `b.v`, force is `-C*(b.v)` and modal reaction is `-C*b*(b.v)`. The complete
resistance is `C*b*b^T`, including cross-mode terms. Its dissipated power is
`C*(b.v)^2`. A mode with a node at that attachment has no direct damping there;
changing attachment position changes which motion is suppressed. Multiple ports
add their resistances. They do not damp modes independently or retune frequencies.

The existing reference/prepared nonlinear owner receives this positive
semidefinite momentum resistance inside its implicit solve. The linear image
compiles the same physical row into the existing simultaneous bilateral port
solve. Dissipation appears in the existing energy reports. Solid mufflers have
exact zero participation on both sticks, every snare wire and appended gas/neck
coordinates. These other bodies can still respond **indirectly** through their
original reciprocal mechanical and pressure coupling.

Head weights use the existing piecewise-linear transverse interpolation; shell
weights use the existing reduction's point port. Outside-mesh, mounting-hole and
nonmoving head stations refuse rather than snapping to a nearby surface. Zero
resistance preserves the original numerical path but still requires a valid
attachment. No option means no new resistance. At most 16 physical ports are
accepted; the prepared image additionally retains its original eight-connection
budget, shared with volume or cavity springs. No limit is silently enlarged.

## Deliberate limits

This is a **fixed point dashpot**, not unilateral hand contact, a finite-area gel
pad, a clamp, added mass, friction, finger preload, or a timed grab/release gesture.
It is active from the first mechanical step. A cymbal's full physical choke needs
those additional contact and hand mechanics; this option does not claim them.
A finite patch cannot in general be replaced by averaging its modal row: its
local squared velocities must be integrated. Sample-rate, modal-bandwidth,
geometry/convergence and calibration restrictions of the original examples
remain. No real-time or full-band cymbal-fidelity claim follows.

Focused native regressions:

```sh
cargo test --release -p fs-couple --lib render::plate::impact::damping
cargo test --release -p fs-couple --example percussion muffling -- --test-threads=1
```
