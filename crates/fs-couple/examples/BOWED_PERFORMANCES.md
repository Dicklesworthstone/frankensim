# File-driven bowed strings

Render the supplied string, friction, body and gesture trajectory without
writing Rust fixture code:

```sh
cargo run --release -p fs-couple --bin music_render -- \
  ensemble /tmp/bowed.wav --bow crates/fs-couple/examples/bowed-string.performance \
  --full-scale-pa 1 --block 37
```

Add `--plate crates/fs-couple/examples/plate-mesh.performance` to mix the existing
mesh/material-derived plate at the same observer. Both examples declare 4801
samples at 48000 Hz. Independent mechanical systems do not acquire force
coupling merely by being mixed. `--decimate` explicitly admits higher-rate
sources; see `ENSEMBLE_PERFORMANCES.md` for complete-duration and delay semantics.

## Input

The fixed-order header uses ASCII whitespace-separated fields:

```text
frankensim-bowed-performance-v1
audio RATE_HZ SAMPLES FULL_SCALE_PA
string LENGTH_M TENSION_N LINEAR_DENSITY_KG_M EI_N_M2 ETA_I_N_M2_S MODE_COUNT
damping ZETA_1 ... ZETA_MODE_COUNT
stribeck MU_STATIC MU_KINETIC CHARACTERISTIC_SPEED_M_S
subsamples FRICTION_REFRESH_STEPS_PER_SAMPLE
body AREA_M2 MASS_KG FREQUENCY_HZ DAMPING_RATIO
ambient TEMPERATURE_K PRESSURE_PA RELATIVE_HUMIDITY
listener_m DISTANCE_M
compile_limits MAX_SOURCE_VISITS MAX_COMPILED_CONTROLS
schedule
```

After the LF-delimited `schedule` line, append exact canonical bytes from
`fs_scenario::gesture::GestureSchedule::to_canonical_bytes()`. The example shows
its actual tab-separated representation. Exactly one track targeting bow string
zero is supported; every other target or extra track refuses. Reversals, station
changes, interrupted ramps, lift and re-entry use the existing sampler and
`ScheduledBowedRenderer`, including the original control-tick/audio-clock rule.
Commands after the last observed control tick refuse instead of disappearing.
A ramp that starts in the window may extend beyond it.

The initial gesture retains the solver's positive-load admission. An explicit
step at time zero can lift the bow before any physical sample; later zero loads
release it without erasing string or body ringdown. No note-to-frequency mapping,
amplitude envelope, gain preset or alternate friction integrator is introduced.

The existing pinned Euler–Bernoulli string derives modes from length, tension,
density and flexural stiffness. Damping must include any declared viscous bending
loss. The compact body is still the existing one-way, rigid-bridge approximation,
not a measured soundboard or a reciprocal string/body junction. Friction and
material numbers are caller supplied, not sourced material identity. The shared
moist-air model supplies observer density only; it does not heat the string or
retune it. Its own validity limits and all original solver limits remain.

Admission caps are 1 MiB, 512 string modes, 256 friction substeps, 16384 authored
commands, 262144 compiled controls, 16777216 source-sampling visits and 28800000
mechanical samples, with duration at most 600 seconds and rate at most 192000 Hz.
Smaller authored compilation budgets remain binding. These are finite workload
bounds, not a hard-real-time or continuum-convergence claim.

Library entry point: `bowed_string::runtime::schedule::file::BowedPerformance`.
Its `into_renderer()` returns the existing `ScheduledRenderer`; no source file
is reopened during playback. Source identity and actual inner bow-control count
are retained in ensemble output metadata. Run the direct solver/CLI regressions
with `cargo test --release -p fs-couple --test bowed_file`.
