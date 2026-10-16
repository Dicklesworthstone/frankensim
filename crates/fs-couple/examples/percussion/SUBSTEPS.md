# Hard impacts without changing the playback clock

`--impact-substeps DEPTH ATTEMPTS` enables bounded internal refinement of the
existing nonlinear drum or cymbal step. It implies prepared nonlinear execution.
Add `--analytic-newton` to use analytic tangents; without that flag the original
finite-difference Jacobian is retained. Omitting this option leaves the original
fixed-step realization unchanged.

```sh
# Curved shell, stand felt and a hard physical stick launch.
cargo run --release -p fs-couple --example percussion -- \
  splash 4096 --analytic-newton --impact-substeps 8 511 \
  --strike-speed-m-s 4 > hard-splash.csv

# Both stretching heads, cavity inertia, independent sticks and a spatial pad.
cargo run --release -p fs-couple --example percussion -- \
  drum-stretch 4096 --analytic-newton --impact-substeps 8 511 --cavity-modes \
  --strike-speed-m-s 4 --strike-position-m 0.06 0.01 \
  --second-stick-position-m -0.05 0.02 --second-stick-speed-m-s 2.5 \
  --muffler batter 0.08 0.01 0.4 > hard-drum.csv

# The existing BEM receiver still consumes exactly the original 768 kHz clock.
cargo run --release -p fs-couple --example percussion -- \
  splash-mic 4800 20 --analytic-newton --impact-substeps 8 511 \
  --strike-speed-m-s 4 > hard-splash.wav
```

These are usage examples, not claims that arbitrary hard strikes converge or
that the full bandwidth of a real cymbal is represented. All original geometry,
mode-count, material validity, input-force and energy ceilings still apply.

## What happens inside a sample

The existing prepared Gonzalez solver first attempts the complete mechanical
tick. A Newton-stall or energy-balance refusal may split that interval into two
halves. Accepted leaves use the same contact, nonlinear shell/head, felt and
Kelvin laws and update their actual history. Coarser trials resume only at exact
dyadic boundaries. No interval is omitted and no stick, resonator or felt is reset.

Every leaf sees the original tick's held force, for its own actual duration.
The existing player converts its force program to that held value once per
nominal mechanical tick; refinement does not reinterpret force-file knots.
Both sticks and all head/air coordinates remain in one joint solve. The outer
force schedule advances only when the entire tick succeeds.

Energy allowances are apportioned by leaf duration rather than multiplied by
the number of substeps. Supplied work, viscous loss and irreversible felt loss
are accumulated, and the complete original tick's energy balance is checked.
A failure or cancellation anywhere rolls back the **whole tick**, including
leaves already accepted internally, motion, felt history and sample count.
Only numerical Newton/energy failures are eligible for refinement; physical
invalidity, force/energy ceilings and material densification still refuse.

Depth is 0..10 and bounds the smallest trial to `nominal_dt / 2^DEPTH`.
Attempts is 1..2047 and includes failed parent solves as well as accepted leaves.
For depth 8, at most 256 leaves and 511 total tree nodes are possible. A lower
attempt allowance can stop earlier. These are explicit work ceilings, not a
measured CPU deadline. Exhaustion reports the attempted solves and rolled-back
accepted leaves; it never substitutes silence or partial output.

The output/sample count, CSV timestamps, mechanical-to-audio ratio, acoustic
filter coefficients and PCM rate remain unchanged. The acoustic observer still
uses differences of velocities on that original clock; the new internal samples
are **not** advertised as additional microphone bandwidth or an anti-alias cure.
This is convergence recovery, not local time-error estimation. Temporal, spatial
and modal convergence and instrument calibration require separate validation.

`drum-modal` and snare commands retain their different joint contact solver and
reject this option rather than silently changing models. Supplied geometry,
spatial mufflers, player force files and the existing nonlinear cavity/neck/loss
controls are preserved; vented exterior audio remains unsupported.

Library callers can use `prepared.with_substeps(ImpactSubstepConfig { .. })` at
any accepted sample without resetting the instrument. The returned
`SubsteppedImpactSystem` implements the existing `ImpactSource`; `last_substeps()`
describes the last accepted outer tick, and `into_prepared()` disables refinement
without changing state. Construction preallocates rollback buffers; the stepping
host does not allocate. Storage callbacks and CPU performance remain separate
obligations.
