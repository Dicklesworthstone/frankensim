# UQ over actual coupled Cooling solves

`cooling-network-uq` runs Monte Carlo by executing the existing steady
`cooling-network` product once for every sampled parameter vector. Fixed-count
empirical reporting is the default; optional sequential compliance stops when
a confidence sequence resolves a declared probability target. Neither path
substitutes `fs-uq` unit-test callbacks for the cooling model.

```bash
cargo run -p fs-cli --bin frankensim -- --json cooling-network-uq \
  examples/cooling-network/fan-correlated-hotspot.json \
  examples/cooling-network/uq-fan-hotspot.json
```

The first argument is an ordinary **steady** cooling-network request. It must
have `objective.gradient=false` and may not contain transient or design controls.
The second document uses schema `frankensim.cooling-network-uq.v1` and declares a
fixed lifetime sample budget, deterministic seed, per-invocation evaluation-time
budget, optional temperature limit, explicit dependence, and uncertain parameters.

Supported parameter targets are:

* `air-density`
* `air-specific-heat`
* `inlet-temperature` with its external-inlet index
* `fan-speed-ratio`
* `surface-htc` for a surface that declares scalar `htc_w_m2_k`
* `component-power` for a named `solid.component_power` row

A sampled scalar `h` cannot override a correlation-derived convection surface.
Changing component power also recomputes the declared component total before the
sample is admitted, so the existing power-map audit remains active. Fan speed,
inlet temperature and all other sampled values pass through the ordinary cooling
parser and producer gates on every sample.

Distributions are `gaussian` (`mean`, `std_dev`) or `uniform` (`lo`, `hi`). A
uniform support must remain inside the target's basic physical domain. Gaussian
sampling is not silently truncated: if a sampled value leaves the physical or
model domain, that model evaluation refuses the **entire** UQ execution rather
than dropping the sample and biasing the distribution.

Dependence is explicit:

```json
{"kind":"independent"}
{"kind":"unknown"}
{"kind":"joint-gaussian","matrix":[[1,0.4],[0.4,1]]}
```

`unknown` with multiple random dimensions refuses before a cooling sample is
run. `joint-gaussian` uses the existing `fs-uq` numerical PSD admission and
requires Gaussian marginals. Correlations alone are not interpreted as a copula.

The fixed-count result reports empirical mean, standard deviation, p05/p50/p95,
observed min/max, descriptive standard error, and optional empirical probability
`temperature <= temperature_limit_k`. All are labeled **Estimated empirical
Monte Carlo**. This default `result.v1` does not infer a confidence sequence,
optional-stopping guarantee, model-form uncertainty bound, mesh-convergence
certificate, experimental validation, or native `.fsim`/ledger package.

Each sample is handed to a child invocation of the same `frankensim
--json cooling-network /dev/stdin` binary. This deliberately spends process
startup time to ensure the parser and physics product boundary are identical to a
normal cooling run. The present handoff therefore targets the repository's Unix
platforms (Linux/macOS). The parent feeds stdin and drains bounded child output
concurrently, so a blocked request write cannot disable its deadline watchdog.
The active child is killed and reaped when the invocation's evaluation budget
expires. Input parsing, executable hashing, and final filesystem synchronization
are outside that evaluation-time allowance; this is not a hard real-time bound
on process startup or I/O.

## Checkpoint and resume expensive runs

Add `--checkpoint NEW-PATH` to retain each completed sample using the existing
`fs-uq` checkpoint format. The destination must not exist. Each update is written
to a sibling staging file, synced, and atomically renamed over the output that
this invocation reserved. Unix directory synchronization persists the rename.
A staging/write error refuses the command rather than reporting saved work.
An abrupt stop can leave a `.pending` file; the previously published checkpoint
remains the last accepted prefix. Resume to a fresh destination rather than
removing or reusing that staging file automatically.

For example, pause the eight-sample example after three actual cooling solves:

```bash
cargo run -p fs-cli --bin frankensim -- --json cooling-network-uq \
  examples/cooling-network/fan-correlated-hotspot.json \
  examples/cooling-network/uq-fan-hotspot.json \
  --checkpoint first.uqcp --max-new-samples 3
```

This deliberately returns the `BUDGET` exit status and a
`frankensim.cooling-network-uq.progress.v1` document: three samples retained,
next ordinal three, and termination `sample-chunk`. It does **not** publish a
partial distribution or a compliance decision. `--max-new-samples 0` is allowed
for admitting and checkpointing an empty prefix; a chunk limit always requires
`--checkpoint` so paid work is not discarded.

Continue the same plan, without rerunning the three accepted solves:

```bash
cargo run -p fs-cli --bin frankensim -- --json cooling-network-uq \
  examples/cooling-network/fan-correlated-hotspot.json \
  examples/cooling-network/uq-fan-hotspot.json \
  --resume first.uqcp --checkpoint completed.uqcp
```

A wall-time interruption also retains a checkpoint and returns progress with
termination `wall-time-budget`. If the child was interrupted before delivering
its QoI, that sample is retried with the **same seed and ordinal** on resume.
Completed samples are not rerun; interrupted samples are never skipped or
replaced by zeros. A real child refusal or non-finite result is different: it
permanently refuses the execution and invalidates this invocation's output with
a retained failure diagnosis that the checkpoint decoder rejects.

The original sample count, parameter ordering, distributions, dependence,
threshold and seed remain immutable across resume. `wall_seconds` is a new
per-invocation evaluation allowance and may be changed without changing the
sample sequence; `--max-new-samples` also does not reset the lifetime count.
A completed checkpoint is terminal and returns the completed report without
additional cooling evaluations.

Resume binds the exact base-request bytes (even whitespace), parameter-target
lowering, and executable content, as well as every library plan field. A changed
base model, seed, plan, executable, corrupt payload, or excessive observation
count refuses before any child evaluation or new output reservation. The resume
input is never overwritten; use a new checkpoint path for further progress.
Without `--checkpoint`, the legacy two-argument command still refuses budget
exhaustion without publishing partial statistics.

Checkpoints must come from a trusted source. BLAKE3 detects accidental corruption
and identity mismatch; it does not authenticate the producer or prove that stored
observations came from the model. Do not replace the installation during a run;
Linux pins hashing and child launch to `/proc/self/exe`, while other platforms
use the current executable pathname. The same deterministic software/hardware
profile is still required for bitwise numerical replay.

For deterministic replay, the Philox sample ordinal is keyed to the UQ seed and
fixed plan. Running the same base model and UQ document produces the same sample
vectors; deterministic cooling mode then provides the same empirical reduction
on the same admitted software/hardware profile. The integration regression
compares both final result bytes and final checkpoint bytes for uninterrupted
execution versus a 3+2+3 split through real cooling child solves.

## Stop when a probability target is resolved

Declare all three policy values explicitly, before inspecting the run:

```bash
cargo run -p fs-cli --bin frankensim -- --json cooling-network-uq \
  examples/cooling-network/fan-correlated-hotspot.json \
  examples/cooling-network/uq-fan-compliance.json \
  --compliance-probability 0.95 \
  --confidence-alpha 0.05 \
  --min-decision-samples 32 \
  --checkpoint compliance-first.uqcp
```

This example asks whether **the declared numerical model's probability of its
objective being at most 303 K is at least 0.95**. The sample budget is 2048;
that is an upper limit, not a promise to evaluate every sample. The inputs are
illustrative, not measured hardware or experimentally validated probability laws.

After each accepted cooling solve, the driver calls the existing
`UqExecution::assess_compliance` on every retained indicator `QoI <= limit`.
The producer uses `fs-eproc::GaussianMixtureCs` with sigma=1/2 and rho=1.
After the minimum sample count, a lower bound at or above the required
probability yields `meets-probability-target`; an upper bound strictly below it
yields `below-probability-target`. Otherwise the decision is `indeterminate`.
The empirical success fraction alone is never the stopping criterion.

The response uses `frankensim.cooling-network-uq.compliance.v1`, distinct from
the fixed-count distribution result. It retains the actual sample count, target,
alpha, minimum count, empirical success fraction, probability confidence
sequence, parameter declarations, and model-only authority. No observations means
`null` probability estimates and bounds, not zero uncertainty. It deliberately
does not present stopped temperature means or quantiles as confidence bounds.

Both resolved decisions have status `decision-reached`, termination
`probability-target`, and a successful execution exit code. **A successful exit
is not a cooling-design approval:** consumers must inspect `decision` and its
scope. Exhausted samples, wall time, or an invocation chunk leave status
`inconclusive` and exit `BUDGET`, even when every planned sample was evaluated.
A model failure remains a refusal with no confidence result from a filtered prefix.

The confidence-sequence mathematics is time-uniform under its fixed-conditional-
mean sampling assumptions, for a fixed model, law, threshold and predeclared
alpha. This licenses looking after each sample; it does not control searches
across seeds, designs, confidence levels or probability targets. Numerical bounds
are not outward-rounded certificates. Geometry, mesh, solver and model error
remain outside this sampling statement; **no physical safety certification** is
created. See Howard, Ramdas, McAuliffe and Sekhon, *Time-uniform, nonparametric,
nonasymptotic confidence sequences*, Annals of Statistics (2021),
DOI 10.1214/20-AOS1991, and the existing `fs-uq` compliance producer's contract.

The current general-purpose Gaussian mixture can be conservative for extreme
reliability targets. A tight target may remain inconclusive at the product's
10,000-sample cap even when every observed sample complies. Neither the sample
budget nor an all-success empirical fraction relaxes the probability criterion.

### Resume with the same stopping policy

All three policy values and the confidence-method version are bound into the
checkpoint identity, in addition to the original model and UQ plan. Repeat them
on resume; changing or omitting any one refuses before evaluating a sample or
reserving a new output. The lifetime sample budget remains immutable.

```bash
cargo run -p fs-cli --bin frankensim -- --json cooling-network-uq \
  examples/cooling-network/fan-correlated-hotspot.json \
  examples/cooling-network/uq-fan-compliance.json \
  --compliance-probability 0.95 \
  --confidence-alpha 0.05 \
  --min-decision-samples 32 \
  --resume compliance-first.uqcp --checkpoint compliance-next.uqcp
```

A checkpoint already at its first resolved decision returns that decision
without another cooling evaluation. An unresolved prefix continues at its next
ordinal. The sequential integration regression compares uninterrupted execution
with 5+4+remaining chunks, including exact stopped-result and checkpoint bytes,
and verifies terminal replay with no time allowance for another child solve.
