# Rare-event probabilities through actual cooling solves

`cooling-network-uq BASE UQ --importance EPSILON,SHIFT...` uses the existing
`fs_uq::ImportanceExecution` to evaluate the strict event
`QoI > temperature_limit_k`. The limit must be present in the UQ document;
equality is compliant. Every observation calls the original cooling product,
with its original geometry, materials, target lowering and child watchdog.

```sh
cargo run -p fs-cli --bin frankensim -- --json cooling-network-uq \
  examples/cooling-network/fan-correlated-hotspot.json \
  examples/cooling-network/uq-importance-hotspot.json --importance 0.2,3
```

The example declares Gaussian inlet-temperature scatter with mean 300 K and
standard deviation 0.25 K. It spends 64 cooling solves. The illustrative shift
of three is in independent standard-normal coordinates, NOT three kelvin.
Choose the proposal, threshold, sample count and model before seeing these
production responses. This example does not establish that the selected limit
is rare or that the proposal is effective for this particular cooling device.

Provide exactly one finite shift per declared parameter in declaration order.
The defensive target mass EPSILON must be in (0,1] with finite reciprocal.
The proposal is `epsilon*phi(z) + (1-epsilon)*phi(z-shift)` in independent germs;
the existing admitted Gaussian factor maps these to physical parameters. Explicit
independence or a joint-Gaussian law is required, including admitted singular
laws. Uniform marginals, unknown dependence and invalid shifts refuse before
physical evaluation. Input-domain failures are not clipped, replaced or redrawn.
Zero shift or epsilon one retains the ordinary target Monte Carlo sampling path.

## What the result means

`frankensim.cooling-network-uq.importance.v1` reports the original parameter
law, sampling key, actual proposal, paid calls and the strict failure threshold.
`estimated_failure_probability` is the mean of the target/proposal-weighted
failure indicators. It is NOT the raw frequency of failures under the proposal,
not normalized by the observed total weight, and not clipped into [0,1].
`failure_standard_error` and `failure_effective_samples` are descriptive
sampling diagnostics, not confidence bounds or valid stopping rules.
`mean_likelihood_weight` diagnoses importance weighting; its population value
under the ideal sampling model is one. `proposal_failure_count` is explicitly
a diagnostic of the proposal, not a target-law probability.

No physical temperature mean, quantiles or distribution are inferred from the
biased proposal responses. No observed failures or a zero empirical standard
error does not establish zero risk. A single shift can miss disconnected modes
or increase variance. Finite-generator grids, numerical likelihood evaluation,
solver error and physical-model validity are not bounded by this report.

Steady objectives and explicit `qoi.kind=transient-sampled-peak` share the same
adapter. One transient observation is the initial/accepted-endpoint maximum of
one COMPLETE declared trajectory, including fixed repeated cycles. It is not a
single time step or a continuous-time peak certificate.

The existing 10,000-call ceiling applies. Without recovery options, one soft
wall allowance includes setup, all physical solves and report construction.
Bounded synchronous work is not preemptible; an expired allowance produces
BUDGET and no partial probability. Physical refusal remains terminal even when
a deadline also expires. There is no fallback to ordinary Monte Carlo.

QMC, PCE, Sobol/Shapley sensitivity, mean controls, sequential compliance and
candidate selection remain incompatible. These conflicts refuse before input
reads or output reservation. The separate `--importance-pilot` mode learns a
proposal but does not yet have a durable combined checkpoint; its recovery flags
still refuse. See `PILOT_IMPORTANCE_UQ.md` for that mode.

## Recover every completed physical solve

Fixed `--importance` proposals now support the existing recovery options:

```sh
frankensim --json cooling-network-uq \
  examples/cooling-network/fan-correlated-hotspot.json \
  examples/cooling-network/uq-importance-hotspot.json --importance 0.2,3 \
  --checkpoint first.imp --max-new-samples 3

frankensim --json cooling-network-uq \
  examples/cooling-network/fan-correlated-hotspot.json \
  examples/cooling-network/uq-importance-hotspot.json --importance 0.2,3 \
  --resume first.imp --checkpoint completed.imp
```

The first command saves three responses and their original scaled likelihood
weights. It exits BUDGET with `importance.progress.v1`, not a probability report.
The second evaluates only the remaining 61 inputs. Its final report and
checkpoint match an uninterrupted run on the same deterministic runtime profile.
Repeated chunking retains the original lifetime count. Explicit chunk limits,
including zero, require a new checkpoint destination.

A successful save is staged, synchronized and atomically published through the
existing file transport after each accepted call. Source checkpoints are never
overwritten: every output path must be new, and using the resume path as output
refuses. Filesystem failures stop the command. A real physical/model failure
invalidates this invocation's output with a terminal diagnosis, so a failed
sample cannot be skipped by resuming its earlier successful prefix.

Each invocation has a fresh evaluation-time `wall_seconds` allowance after
admission/restoration; it never refunds paid calls. An unfinished cooling child
is killed and reaped by the existing watchdog and retries the same ordinal on
resume. An unfinished transient restarts from its declared initial state;
completed trajectories are not repeated. When all original observations are
already retained, resume returns the complete result with no further physical
calls even under a tiny new evaluation-time allowance.

A shortened run publishes only counts, next ordinal, proposal, observable and
checkpoint location. Neither weighted failure probability nor empirical error
is reported from that prefix. Without a destination, an interrupted resume
returns a budget diagnostic; its source file remains unchanged. The complete
result keeps the existing `importance.v1` schema and statistics unchanged.

The `FSIMP001` codec binds every plan field, original budget, exact ordered
proposal shifts, defensive mass, threshold, base-request bytes, target lowering,
executable and sampler/likelihood version. A changed proposal, law, model, seed,
threshold or allocation refuses before output reservation or physical work.
Only the per-invocation wall allowance and chunk length may change. Each record
retains 16 bytes (QoI and scaled weight); file reads are bounded independently
of the smaller scalar checkpoint formats. No other sampler's codec changes.

Checksums detect corruption, not authenticity or physical provenance. Restore
only trusted files: weight range checks do not prove that the producer evaluated
the right model or likelihood. Restoration recomputes failure flags and reports,
not physics, and does not provide a model-validation or cross-runtime guarantee.

Focused native checks:

```sh
cargo test -p fs-cli --bin frankensim uq_command::execute::importance
cargo test -p fs-cli --test cooling_importance
cargo test -p fs-uq --test product_importance_checkpoint
```

These tests exercise actual cooling calls, independent reconstruction through
the library, target-MC identity, complete transient observations, admission,
deadlines, terminal physics failures, process restarts, unchanged final report
bytes, immutable proposal bindings and protected checkpoint inputs. They require
native execution before claiming a Rust test pass or a physical-accuracy result.
