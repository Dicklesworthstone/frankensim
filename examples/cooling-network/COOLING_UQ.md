# Empirical UQ over actual coupled Cooling solves

`cooling-network-uq` runs fixed-count Monte Carlo by executing the existing
steady `cooling-network` product once for every sampled parameter vector. It does
not substitute `fs-uq` unit-test callbacks for the cooling model.

```bash
cargo run -p fs-cli --bin frankensim -- --json cooling-network-uq \
  examples/cooling-network/fan-correlated-hotspot.json \
  examples/cooling-network/uq-fan-hotspot.json
```

The first argument is an ordinary **steady** cooling-network request. It must
have `objective.gradient=false` and may not contain transient or design controls.
The second document uses schema `frankensim.cooling-network-uq.v1` and declares a
fixed lifetime sample count, deterministic seed, per-invocation evaluation-time
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

The result reports empirical mean, standard deviation, p05/p50/p95, observed
min/max, descriptive standard error, and optional empirical probability
`temperature <= temperature_limit_k`. All are labeled **Estimated empirical
Monte Carlo**. No confidence sequence, optional-stopping guarantee, model-form
uncertainty bound, mesh-convergence certificate, experimental validation, or
native `.fsim`/ledger package is inferred.

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
