# Two sticks, one physical drum

`--second-stick-position-m X Y` adds an independent inertial stick contacting
that location on the **same batter head**. `--second-stick-speed-m-s V` sets
its initial tip speed in m/s (0..20, default **zero**). First-stick controls
and the original default 0.8 m/s launch remain unchanged.

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
shared head. Simultaneous contacts are solved together with head, air and snare
reactions. No drum state, contact impulse, audio waveform or resonator bank is
copied, reset or mixed to manufacture a second hit. Impact timing and rebound
follow from mechanics, not the launch timestamp alone.

The original two head ranges remain fixed. The second stick is placed after
them, before snare wires and acoustic inertia. It has zero direct swept-volume
or radiation participation: the existing BEM microphone observes the resulting
combined head motion. CSV with two sticks appends displacement and velocity
for each tip, in metres and m/s. Unchanged single-stick commands retain their
original columns and physical construction.

Two sticks work with drum, drum-modal, drum-stretch and snare/snare-off, including
WAV/microphone variants and supplied drum specifications. Prepared linear drum
and snare also admit distributed air. Two-stick **nonlinear distributed air**
is currently refused. The existing modal/energy/slope limits remain limits;
no modes are silently removed to make room for a striker. Positions must be
inside the physical moving-head chart; missing/outside/rim-only stations refuse.

This is two fixed-axis effective-mass sticks, not a full hand/arm controller,
finite-footprint collision detector, stick-stick collision, cymbal duet or rim
shot. The stick geometry and contact material remain estimates. There is no new
calibration, full-band adequacy or real-time performance claim. Native example
regressions cover the real two-contact head motion, energy accounting, snare
coordinate layout and control admission; native Rust execution is required.
