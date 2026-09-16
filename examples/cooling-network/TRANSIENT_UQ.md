# Uncertainty in transient cooling peaks

`cooling-network-uq` can propagate uncertain inputs through complete transient
cooling runs, including the existing nonlinear material, contact, adaptive-step,
fixed-repeat and cycle-rate thermostat paths. It calls the normal cooling
producer for every sample; it does not replace the transient with a steady
solve or an independent lumped model.

```bash
frankensim --json cooling-network-uq \
  examples/cooling-network/nonlinear-contact-pulse.json \
  examples/cooling-network/uq-nonlinear-contact-pulse.json \
  --checkpoint pulse-first.uqcp --max-new-samples 3

frankensim --json cooling-network-uq \
  examples/cooling-network/nonlinear-contact-pulse.json \
  examples/cooling-network/uq-nonlinear-contact-pulse.json \
  --resume pulse-first.uqcp --checkpoint pulse-complete.uqcp
```

The first command deliberately exits `BUDGET` with three completed trajectories,
not three timesteps. The supplied eight-sample example is a small illustrative
study, not evidence of a precise tail probability or validated hardware.

## Select the observable explicitly

Add this field to the UQ document, not the base cooling document:

```json
"qoi": {"kind": "transient-sampled-peak"}
```

The spatial objective still comes from the base's `objective`: mean wall,
maximum wall, selected vertices, or the maximum discrete solid temperature.
One observation is its maximum over the initial state and every accepted
endpoint of ONE completed trajectory. A low final temperature cannot hide a
higher earlier pulse. For `transient.repeat.cycles`, selection uses
`repeated_cycles.sampled_peak_objective_k`, not the last cycle's local peak.
The all-cycle result must report the declared fixed cycle count as complete.

A transient base without the explicit selector refuses. The default selector,
`steady-objective`, retains the existing steady behavior. Nested steady or
transient design searches refuse, as does `repeat.until_periodic`: this UQ mode
requires the same prescribed time horizon in every sample. A fixed number of
cycles, including its fixed controller policy, is supported.

Result and sequential-decision metadata identify the spatial objective,
sampled-endpoint time scope, cycle count and `one-completed-trajectory`
observation unit. A sample-dependent adaptive mesh in TIME can change the
endpoint locations; the probability concerns the declared numerical algorithm,
not a continuous-time peak bound. No individual timestep is counted as an
independent Monte Carlo observation.

## Sample active transient quantities

In addition to the existing fluid, inlet, coefficient and material targets,
these target forms are available (indices are zero-based):

```json
{"kind":"initial-temperature"}
{"kind":"volumetric-heat-capacity"}
{"kind":"element-heat-capacity","element":0}
{"kind":"interval-power-scale","interval":0}
{"kind":"interval-component-power","interval":1,"component":"chip"}
{"kind":"interval-fan-speed-ratio","interval":1}
```

Initial temperature is in kelvin and selects a declared uniform initial state;
it cannot replace a nodal initial field. It is independent of an inlet-air
uncertainty target. Heat capacity is volumetric, in J/(m3 K): the uniform target
requires uniform capacity, while the element target requires an existing
per-tetrahedron array. Neither converts the other's input mode.

An interval power scale multiplies the base sources through the normal workload
producer. An interval component target changes one named absolute watt value in
`component_powers_w`; it leaves the base power map and its fixed footprints
unchanged. These two workload modes are mutually exclusive. Zero scale or watts
is a valid off interval; negative or non-finite draws refuse rather than being
clipped, skipped or redrawn.

Use `interval-fan-speed-ratio` for transient fan uncertainty. The old base
`fan-speed-ratio` target is refused because interval speeds override it. The
base `component-power` target remains useful where an interval uses
`power_scale`, but refuses when all intervals provide absolute component watts.
Interval targets are checked against real schedule indices and input modes
before reserving a checkpoint or launching a child.

One parameter vector is drawn ONCE for a complete trajectory and reused in all
its fixed repeated cycles. This models uncertain fixed schedule parameters,
not a stochastic process redrawn at every timestep or cycle. Durations,
footprints, assignments, time-integration controls and controller rules remain
unchanged. Existing independent and explicitly joint-Gaussian dependence
policies apply across these parameter coordinates.

## Confidence and recovery

The existing `--compliance-probability`, `--confidence-alpha`, and
`--min-decision-samples` options apply to the event
`sampled trajectory peak <= temperature_limit_k`. That limit is declared in the
UQ request; it does not silently inherit the base trajectory's reporting limit.
All three statistical settings must be explicit, and the minimum must fit the
UQ sample cap. The eight-sample example generally cannot resolve a demanding
tail-probability target; increase the declared sample budget before starting
such a study. See `COOLING_UQ.md` for confidence assumptions and exit statuses.

Checkpointing retains completed scalar trajectory observations, not partial
fields inside a trajectory. Parent-timeout interruption discards that unfinished
trajectory and retries the same seed/ordinal from its declared initial state on
resume. Completed trajectories are not rerun. A genuine failed endpoint, child
budget refusal, invalid material evaluation or missing peak is a terminal model
failure; it invalidates the current checkpoint output rather than contributing
a partial peak or filtering the observation.

Observable identity, the entire base model/schedule, UQ plan, executable and
optional confidence policy remain bound by the existing checkpoint machinery.
Use trusted checkpoint inputs and fresh output paths. The same deterministic
runtime profile is required for exact replay.

## Focused product regressions

```bash
cargo test -p fs-cli --bin frankensim uq_command::model
cargo test -p fs-cli --test cooling_transient_uq
```

The actual-command tests compare zero-uncertainty observations with the direct
nonlinear contact-pulse peak, distinguish peak from final temperature, retain
all-cycle peaks, exercise off workloads, compare complete result/checkpoint
bytes after chunking and time exhaustion, and test sequential decisions and
refusals. Their presence is not a claim that Rust execution has been verified
on this checkout.

All results concern declared discrete-model temperatures and input probability
laws. They do not certify inter-step maxima, time/space discretization error,
missing material/physical-model uncertainty, experimental validity or safety.
