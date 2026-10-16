# Hear a fitted physical design without rebuilding it by hand

```bash
cargo run -p fs-ascent --features equilibrium-design --bin equilibrium_fit -- \
  crates/fs-couple/examples/equilibrium-playback.model \
  crates/fs-couple/examples/equilibrium-playback.fit \
  --playback-wav /tmp/fitted-release.wav \
  --playback-case tuned-load --playback-samples 4801 \
  --playback-release 37 --playback-full-scale-pa 0.004 --playback-block 37
```

This example fits a physical load to an explicitly supplied static displacement
of an authored damped mode. The analytical optimum is 2 N, not the template's
1 N. The accepted design is independently re-evaluated, the chosen case is
settled at those same fitted parameters, and its actual displacement/velocity
state is transferred to the existing transient renderer. All case loads release
**before mechanics sample 37**, without a state reset. The sample horizon includes
the final short callback. The example is synthetic, not a calibrated instrument.

The optional playback flags require a WAV path, an exact case name, duration in
mechanics samples, a release sample inside that duration, and an explicit pressure
scale. Block size alone defaults to 512 and may be 1 through 65536. Duration is
bounded at 28800000 samples. All output is at the source model's mechanics rate:
there is no implicit resampling, gain, normalization, padding or filter tail.
The source pressure transfers and damping are unchanged. A source with all-zero
acoustic transfers produces silence, even when its mechanical fitting succeeds.

The existing fitting command, its constraints, stop reasons and nominal output
remain unchanged without playback flags. Combining playback with `--scenarios`
is currently refused rather than silently choosing a realization. Native API
callers may explicitly supply any admitted decision point and selected case.
No dependency or file schema is added.

## State and physical work

`EquilibriumDesign::playback_case` is implemented in
`design::forward::playback`. It reuses the design evaluator's parameter binding
and selected-case preload, including shared spring/contact fields and
case-specific actuator-force variables. It returns `DesignPlayback`, whose
`into_renderer` transfers the complete state to `ScheduledRenderer` and the
existing PCM stream APIs. Its `CaseForceEvent` supports separate releases and
later re-excitation of each load in the selected experiment; two loads at the
same attachment remain independent. The command offers the common release-all
case. It does not prescribe a hammer impulse or synthesize a force waveform.

One candidate preparation and one selected-case preload are charged. The shared
preparation still requires capacity for a full case family even though playback
solves only its selected case. The command reserves that extra preparation within
`--evaluations`, in addition to the existing final objective/adjoint audit; it
requires at least three evaluations. Counts in the output include that work.
Transient work is separately bounded by the explicit sample horizon, admitted
component/contact caps and existing per-step nonlinear budgets. Event compilation
uses at most 262144 modal controls and 16777216 projection terms, rather than
reinterpreting the static template's unused event budgets.

A budget-stopped fit can be played back, but remains labelled unconverged. The
`playback` object identifies the actual case, rate, duration, release, scale,
clipping and initial stored energy. A static objective or constraint certificate
does **not** establish transient audio accuracy, identify damping/radiation, or
prove that static response limits remain satisfied after releasing the loads.
The original mechanical and contact safety/admission limits remain active.

## Files and failures

The destination is exclusively created after input, fit, audit and playback
construction succeed. Existing files are never replaced, including the race
between the early existence check and file creation. Runtime physics or I/O
failure may leave an incomplete WAV; the command then emits an error, no success
JSON, and does not claim that the file is finalized. Inputs are never modified.
The output writer retains its bounded block buffers and original PCM scaling.

Focused native checks:

```bash
cargo test -p fs-couple --test equilibrium_playback
cargo test -p fs-ascent --features equilibrium-design \
  --bin equilibrium_fit --test equilibrium_playback_cli
```

Nine Rust tests were added across the native and command increments. They were
not executed in the authoring environment because Cargo/rustc are unavailable.
Independent matrix-exponential and load-projection checks are not a substitute
for running these tests or the repository's deterministic numerical code.
