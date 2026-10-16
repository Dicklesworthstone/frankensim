# Replicated quasi-Monte Carlo cooling uncertainty

Use the same physical request and uncertainty input as `cooling-network-uq`,
with an explicit number of independently keyed Owen-scrambled Sobol replicates:

```sh
cargo run -p fs-cli --bin frankensim -- --json cooling-network-uq \
  examples/cooling-network/fan-correlated-hotspot.json \
  examples/cooling-network/uq-fan-hotspot.json --qmc-replicates 2
```

This example's eight samples become two complete four-point nets. For a larger
study, set `samples` to 1024 and use `--qmc-replicates 8` (128 points per net).
The total sample count must divide exactly by the replicate count; points per
net must be a power of two, at least two. Replicates must be in 2..=256 and the
existing command's 10,000-total-sample ceiling still applies. At most ten
parameters are admitted; there is no automatic Monte Carlo tail or fallback.

Every point runs the actual cooling child command. Steady objectives and
`qoi.kind=transient-sampled-peak` both reuse the existing physical request
lowering and completed-trajectory observation path. One transient observation
is one complete trajectory, not an individual time step or an early prefix.
Material, contact, fan, load, and transient uncertainty targets retain their
existing admission and probability-model restrictions.

The `frankensim.cooling-network-uq.qmc.v1` response includes `mean_k`,
`replicate_means_k`, `between_replicate_standard_error_k`, and, when a temperature
limit is supplied, `estimated_probability_of_compliance` and
`between_replicate_probability_standard_error`. Errors are computed from
replicate means/proportions, never by treating dependent points inside a net
as iid draws. The report does not attach a confidence interval or approval.

## Retain and resume expensive cooling samples

Add `--checkpoint NEW-PATH` to save every completed model evaluation, including
points inside an unfinished net. Chunking does not change the fixed integration
rule: the lifetime sample count and number of replicates stay immutable.

```sh
cargo run -p fs-cli --bin frankensim -- --json cooling-network-uq \
  examples/cooling-network/fan-correlated-hotspot.json \
  examples/cooling-network/uq-fan-hotspot.json --qmc-replicates 2 \
  --checkpoint first.qmc --max-new-samples 3

cargo run -p fs-cli --bin frankensim -- --json cooling-network-uq \
  examples/cooling-network/fan-correlated-hotspot.json \
  examples/cooling-network/uq-fan-hotspot.json --qmc-replicates 2 \
  --resume first.qmc --checkpoint completed.qmc
```

The first invocation exits BUDGET and returns
`frankensim.cooling-network-uq.qmc.progress.v1`: three accepted samples, zero
completed replicates, and next replicate/point ordinals 0/3. The second evaluates
only the remaining five points. Its completed result and checkpoint match an
uninterrupted run exactly on the same deterministic runtime profile.
`--max-new-samples 0` retains an empty prefix and requires a checkpoint destination,
as does every explicit chunk limit.

A wall-time interruption also returns progress and BUDGET when a checkpoint
output is supplied. An interrupted cooling child is killed and reaped; its
unfinished observation retries the same Sobol point on resume. A transient
sample restarts its entire declared trajectory. `wall_seconds` is a fresh
per-invocation evaluation allowance and may change on resume. A completed
checkpoint returns its final report without any further model evaluations.

Progress reports contain work counts, layout, next ordinals and checkpoint path.
They contain no mean, replicate error or compliance estimate, even if some nets
are complete: a value-dependent timeout cannot turn a shortened layout into a
completed fixed-budget inference. Without `--checkpoint`, a timeout exits BUDGET
without a partial distribution. Sequential compliance and candidate-selection
flags remain unavailable for QMC.

Recovery uses the existing bounded, atomic checkpoint file transport. Destinations
must be new; input checkpoints are never overwritten. Each successful update is
written to a new sibling staging file, synchronized, renamed over this invocation's
reserved output, and directory-synchronized on Unix. A filesystem error refuses
the command. A physical solver failure also refuses the run and replaces this
invocation's checkpoint with a terminal failure diagnosis, so a rejected point
cannot be silently dropped and the remaining samples reported as valid.

QMC checkpoint framing is distinct from MC. It binds exact base-request bytes,
parameter-target lowering, executable content, every plan field, total sample
budget, replicate count, points per replicate, and versioned sampler semantics.
Changed layouts, models, seeds and corrupted files refuse before evaluating a
sample or reserving a new output. The library's `QmcExecution::checkpoint` and
`QmcExecution::restore` expose the same durable recovery for other consumers.
Checkpoints must come from a trusted source: checksums detect corruption and
identity mismatch but do not authenticate the producer or prove the observations.

Sampling uses the midpoints of the existing 32-bit Sobol grid and the existing
approximate inverse-normal transform. Between-replicate errors do not bound
finite-grid bias, transform error, solver error, physical-model discrepancy,
or unsampled continuous-time thermal peaks. Zero replicate variation does not
prove an exact answer, nor justify optional-stopping inference.
