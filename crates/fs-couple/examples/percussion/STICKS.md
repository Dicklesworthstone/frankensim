# Two sticks, one physical drum or cymbal

`--second-stick-position-m X Y` adds an independent inertial stick contacting
that location on the **same batter head or curved cymbal shell**.
`--second-stick-speed-m-s V` sets its initial tip speed in m/s (0..20, default
**zero**). First-stick controls and the original default 0.8 m/s launch remain
unchanged.

```bash
cargo run --release -p fs-couple --example percussion -- \
  drum-modal 192 --strike-speed-m-s 4 --strike-position-m 0.06 0.01 \
  --second-stick-position-m -0.05 0.02 --second-stick-speed-m-s 2.5 > two.csv

cargo run --release -p fs-couple --example percussion -- \
  snare-mic 48000 20 --cavity-modes --strike-speed-m-s 2 \
  --second-stick-position-m -0.05 0.02 --second-stick-speed-m-s 1.6 > two.wav
```

Both sticks use the existing geometry-derived effective mass, estimated stick
profile and elastic tip law. They have separate displacement and momentum,
separate geometric contact rows, and equal-and-opposite reactions against one
shared resonator. Simultaneous contacts are solved together with head, air and
snare reactions, or with the nonlinear curved shell and its stand felts. No
resonator state, contact impulse, audio waveform or resonator bank is copied,
reset or mixed to manufacture a second hit. Impact timing and rebound
follow from mechanics, not the launch timestamp alone.

The original two head ranges remain fixed. The second stick is placed after
them, before snare wires and acoustic inertia. It has zero direct swept-volume
or radiation participation: the existing BEM microphone observes the resulting
combined head motion. CSV with two sticks appends displacement and velocity
for each tip, in metres and m/s. Unchanged single-stick commands retain their
original columns and physical construction.

For `splash`, the original shell mode range also stays fixed. The second stick
follows the shell modes; the time owner places its private Kelvin coordinates
after both sticks and preserves every felt state. Both contact rows sample the
actual curved-shell reduction at their own reference-XY stations. The BEM
boundary continues to radiate only that same shell, never a direct stick signal.
Stand-felt and muffler rows have exact zeros at the second stick coordinate.

```bash
# Two separate impacts on one nonlinear cymbal, with bounded hard-impact recovery.
cargo run --release -p fs-couple --example percussion -- \
  splash 4096 --analytic-newton --impact-substeps 8 511 \
  --strike-position-m 0.06 0.01 --strike-speed-m-s 0.8 \
  --second-stick-position-m -0.05 0.02 --second-stick-speed-m-s 0.8 \
  --muffler shell 0.075 0 0.2 > two-cymbal.csv
```

The same second-stick controls and force files work on `splash-wav` and
`splash-mic`, including `--microphone-right` stereo and supplied
`--shell-profile` geometry. There is one joint mechanical trajectory regardless
of receiver count. Internal recovery leaves force programs and acoustic sample
timing unchanged; see `SUBSTEPS.md`. Harder strikes may still refuse at physical
or numerical limits rather than silently discarding energy or changing the model.

## Independent player forces

`--stick-force-file left.csv` and `--second-stick-force-file right.csv` apply
separate signed forces to the two actual tip coordinates. Either file may be
used alone. The second file requires an explicit second-stick position; it
never creates a phantom striker or routes a force directly to a head mode.
Both files use the bounded headerless `time_s,force_n` format in `DRIVE.md`.
Each begins and ends at zero force, with zero force outside its own interval.

```bash
cargo run --release -p fs-couple --example percussion -- \
  drum-modal-mic 48000 20 --strike-speed-m-s 0 \
  --strike-position-m 0.06 0.01 --second-stick-position-m -0.05 0.02 \
  --stick-force-file left.csv --second-stick-force-file right.csv > played.wav
```

The example starts both sticks at rest. Positive force pushes toward the head
or shell; negative force lifts. These are player forces, not prescribed impact pulses.
You can supply separated, alternating or overlapping pushes and lifts, but
actual contact times depend on the continuing stick and head motion. Neither
input teleports a stick or resets a ringing head. Initial launch speeds remain
independent inputs; the first stick still defaults to 0.8 m/s unless specified.

All programs are admitted before stepping, then integrated independently over
each mechanical interval and staged into one force vector. One joint owner
solve advances the resonator and both player schedules. On refusal, both input
ticks remain pending alongside the unchanged physical state. The shared `player_work_j`
column reports the step's total external work, not a per-hand or cumulative
quantity. Pressure rendering consumes the same jointly driven mechanical state.

Two sticks work with splash, drum, drum-modal, drum-stretch and snare/snare-off,
including WAV/microphone variants and supplied shell/drum specifications. Both
nonlinear drums (including stretching heads) and prepared linear drum/snare
admit distributed air. `--prepared-nonlinear` retains that same nonlinear two-stick/cavity model.
The existing cavity/neck controls remain drum-only; vented exterior audio is
still refused. The existing modal/energy/slope limits remain limits;
no modes are silently removed to make room for a striker. Positions must be
inside the actual vibrating-surface chart. Missing/outside locations, a cymbal's
mounting hole, and nonmoving head stations refuse rather than snapping to an edge.

This is two fixed-axis effective-mass sticks, not a full hand/arm controller,
finite-footprint collision detector, stick-stick collision, 6-DOF stand model
or rim shot. The stick geometry and contact material remain estimates. There is
no new calibration, full-band adequacy or real-time performance claim. Native example
regressions cover the real two-contact head motion, energy accounting, snare
and cavity coordinate layout, shared-shell backreaction and felt memory,
independent force/work scaling, joint retry and control admission; native Rust
execution is required.
