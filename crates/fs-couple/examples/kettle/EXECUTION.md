# Bounded contact recovery on the unchanged output clock

Kettle playback now prepares the existing analytic impact solver with bounded
Newton backtracking and the existing transactional dyadic substep wrapper. This
addresses a concrete delivery blocker: a fixed Newton attempt can fail during a
real felt strike even when geometry and acoustic preparation succeed. Native
confirmation of the original failing kettle tests is still required.

The default work controls are eight Newton backtracking halvings, substep depth
eight and 511 solve attempts per output tick. Every leaf includes the SAME head,
felt/Kelvin, interior air, edge actuator and selected exterior memory. External
force is held through the leaf duration. No contact is omitted, no impulse is
replaced by a force envelope, and no waveform is normalized to hide a failure.

```sh
cargo run --release -p fs-couple --example kettle -- \
  crates/fs-couple/examples/kettle/illustrative-parabolic.fsk 4096 \
  --impact-substeps 8 511 --newton-backtracking 8 > guarded-motion.csv

# Exact old single-step, undamped iteration selection for comparison.
cargo run --release -p fs-couple --example kettle -- \
  crates/fs-couple/examples/kettle/illustrative-parabolic.fsk 4096 \
  --impact-substeps 0 1 --newton-backtracking 0 > original-motion.csv
```

These options also apply to `--wav`, with either exterior acoustic image, and to
`--edge-tensioner` / `--edge-force`. They do not change the specimen, its force
schedule, the source mode addresses, PCM scale or either sample rate. CSV columns
remain unchanged. The work settings are printed on stderr before preparation.

## What is preserved

Backtracking tries the full Newton direction first, then powers of one half,
using the complete original Gonzalez residual. It reuses the existing scratch;
the accepted-state and energy equations are unchanged. A refused search enters
the same stagnation/best-iterate acceptance as the original solver. Nonfinite or
invalid constitutive evaluations propagate rather than being ignored. The
reference stepper and all unselected `StepWorkspace` callers remain unchanged.

A failed numerical/energy trial may bisect the time interval through
`SubsteppedImpactSystem`. Physical validity failures are not repaired. Material
history advances only for accepted leaves. If any later leaf fails or the work
budget is exhausted, the COMPLETE output tick restores state, felt history,
clock and prior accepted substep report. Energy allowances are divided among
leaves and checked again on the whole tick; subdivision does not multiply them.

Force events are read once from accepted output time, held over internal leaves,
and not consumed on refusal. Audio samples observe only complete mechanical
ticks. No microphone, decimator or radiation filter runs on a partial tick.
The read-only mechanical observations still expose the last accepted substep
count and the last leaf's Newton work counters. These are not execution times.

## Limits and checks

The controls accept depth 0..10, 1..2047 solve attempts and 0..16 backtracking
halvings; zero backtracking preserves the original traversal. Finite workspace,
clock, slope, strain, energy, force and material limits remain owned by the
existing components. This is bounded solve recovery, NOT a local truncation-error
controller, aliasing estimator, unique-root guarantee or real-time certificate.
Large attacks may still refuse. Do not infer native success or a calibrated sound
from the existence of the command.

The original failing strike/history regressions are left intact. Four added
example tests exercise configuration, clock/parity, a felt/creep/head/air/load
strike with retry, and retained edge-force events across refusal. Five shared
solver tests cover nonlinear convergence, legacy parity, work/energy, every
cooperative cancellation boundary, exhausted search and condensed/full solves.

```sh
cargo test --release -p fs-phs --lib prepared::backtrack
cargo test --release -p fs-couple --example kettle -- --test-threads=1
```

Rust/Cargo are unavailable in the authoring environment; these new native tests
have not been executed. Independent Python equations exercise hard Hertz contact
and a driven/damped coupled quartic system, not this Rust kettle or its BEM audio.
