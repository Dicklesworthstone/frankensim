# Render your own reduced structural model

```bash
cargo run -p fs-couple --bin music_render -- \
  modal crates/fs-couple/examples/modal-ports.performance \
  /tmp/modal-ports.wav --block 37
```

Use new output paths: the command refuses an existing WAV or provenance sidecar.
The example declares two independent modal voices, three modes, and three force
ports. The first voice starts under static preload, then releases and changes
its two actuators independently. The second starts with a supplied modal velocity
and receives a later physical force. These are authored numerical examples, not
measured instruments or sourced material specimens.

The command executes the supplied model, not a fixed fixture. It uses the existing
mass-normalized modal runtime, physical-force projection, sample-accurate scheduler,
and streaming PCM16 exporter. Block size only partitions output work. Mode data,
force timings, duration, pressure scaling, and numerical limits belong to the file;
`--seconds`, `--full-scale-pa`, and `--schedule` are refused in `modal` mode rather
than overriding or silently ignoring them. Existing `reed`/`string` commands are
unchanged. The command requires 48000 Hz; it does not retune or resample input.

## File layout

ASCII whitespace separates fields. Record order and field counts are exact;
comments, blank records, unknown fields, and trailing records are not accepted.
Indices start at zero. Voice, mode, and port order define their shared basis.

```text
frankensim-modal-performance-v1
sample_rate_hz RATE
samples EXACT_OUTPUT_SAMPLE_COUNT
full_scale_pa FULL_SCALE_PRESSURE
limits NYQUIST_FRACTION MAX_Q MAX_V MAX_ENERGY_J MAX_PRESSURE_PA
compile_limits MAX_EXPANDED_CONTROLS MAX_PROJECTION_TERMS
voices VOICE_COUNT
voice retain-state|static-preload MODE_COUNT PORT_COUNT
mode OMEGA_RAD_S DAMPING_RATIO H_RE H_IM Q_INITIAL V_INITIAL
... exactly MODE_COUNT mode records ...
port INITIAL_FORCE_N SHAPE_0 ... SHAPE_N_MINUS_1
... exactly PORT_COUNT port records, then the remaining voices ...
events FORCE_EVENT_COUNT
force SAMPLE VOICE_INDEX PORT_INDEX FORCE_N
... exactly FORCE_EVENT_COUNT force records ...
```

Every mode is mass-normalized: Q has units m sqrt(kg), V has units m sqrt(kg)/s,
and energy is in joules. OMEGA is angular frequency in rad/s, not Hz or a MIDI note.
H_RE/H_IM are the real/imaginary components of pressure per modal velocity,
in Pa s/(m sqrt(kg)), under the existing exp(-i omega t) convention. They must come
from your declared radiation model; the importer does not manufacture gains.

Each port column contains mode-shape values in 1/sqrt(kg) at its physical actuator.
The generalized force is the ordered sum of each column times its force in newtons.
The same column recovers port velocity from modal velocity. Signed forces and
shape weights are legal. Zero force releases only that port. Same-sample changes
are grouped before projection, with the last assignment to the same port winning.
Events apply before their named sample and must lie in `[0, samples)`.

`retain-state` preserves all supplied Q/V values. `static-preload` requires zero
Q/V placeholders and initializes the actual state from the equilibrium under the
initial port forces. This models a load settled before the simulation window;
releasing it excites the existing resonator. It is not a hammer/contact model.

Numerical limits apply independently to each voice, not to the summed observer
channel. The existing encoder reports clipping after physical mixing; it does not
normalize the waveform. The input's compilation limits bound projected controls
and multiply-add terms, not elapsed time. The loader additionally caps bytes at
4 MiB, voices at 64, total modes at 4096, total port weights and events at 65536,
expanded controls at 262144, projection terms at 16777216 and samples at 28800000.
These are admission ceilings, not performance guarantees.

## Replay and scope

An independent v1 performance can opt into computed broadband receivers by
appending this block after its last force event:

```text
broadband VOICE_COUNT
transfer PRODUCER_ID DELAY_S PRODUCER_REL_ERROR ORDER ITERATIONS DELAY_HALF MAX_REL_ERROR TRAINING_COUNT HELD_OUT_COUNT
frequency OMEGA H0_RE H0_IM ... H_LAST_RE H_LAST_IM
... exactly TRAINING_COUNT complete frequency rows ...
withheld OMEGA H0_RE H0_IM ... H_LAST_RE H_LAST_IM
... exactly HELD_OUT_COUNT complete withheld rows, then the remaining voices ...
```

Each voice needs its own block in original voice order. Samples are computed
pressure per generalized modal velocity in Pa s/(m sqrt(kg)), under
`exp(-i omega t)`, including `DELAY_S`. Both grids must be increasing, positive,
disjoint, below Nyquist, and withheld points must lie inside the training band.
Declare the actual producer identity (for example an fs-bem run hash) and its
discretization error estimate separately from `MAX_REL_ERROR`; neither is inferred
or absorbed into the fitting budget. Both grids cap at 512 points. The complete
4096-mode basis is retained in batches of at most 256 input fields.

The receiver fits a continuous proper rational operator, removes the declared
pure delay before fitting, and advances the augmented mode/filter system exactly
under the original held forces. It uses no velocity differences or Tustin warp.
The fractional delay uses an explicit centered Blackman/sinc kernel with the
declared half-length, admitting only causal taps within the physical delay.
An integer delay is exact. Acoustic filter/delay history starts at zero; this is
bandlimited interpolation of sampled pressure, not exact within-step retiming of
the continuous held-force images. The receiver does not load the mechanical mode.

Admission measures the maximum relative error on the withheld computed operator,
including the actual delay kernel. It is a sampled estimate, not a rigorous
continuum bound. The sidecar retains `EstimateOnly`, producer uncertainty, fit
controls, actual fractional delay, and the observed error. Constant samples with
zero delay retain bitwise narrow-band equilibrium-relative observation. Coupled
or nonlinear v2-v5 trajectories refuse this independent-mode observer.

Use `scripts/e2e/music_lane.sh --fixture modal-broadband --out-dir NEW_PATH
--binary MUSIC_RENDER --qualification-log NATIVE_BROADBAND_TEST_STDOUT` for the
analytic sphere replay/refusal lane and the native three-axis qualification
receipts. Native qualification includes the exact held-sinusoid endpoint oracle;
an unheld continuous sinusoid is a different forcing problem.

The sidecar records the exact input-byte hash, mode/voice/event counts, output
scale, statistics, and the finalized WAV hash. Moving identical input to another
path does not change those identities. Editing whitespace does change the input
hash: it identifies bytes, not a canonical semantic representation.

This path makes supplied linear, mass-normalized structural images executable.
It does not certify their origin, eigenbasis normalization, truncation error,
material validity, or observer transfer accuracy. It is not a general CAD assembly
loader, nonlinear-contact instrument solver, MIDI synthesizer, stereo room model,
or real-time device callback. No runtime dependency was added.

Focused native checks:

```bash
cargo test -p fs-couple --lib render::schedule::force::file
cargo test -p fs-couple --bin music_render
cargo test -p fs-couple --test music_render_model
```

These tests were added with the feature. They were not executed in the authoring
environment, which did not contain a Rust toolchain.
