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

A wall-time interruption exits BUDGET without a partial distribution. A physical
solver error refuses the run rather than dropping the sample. This command
currently rejects `--checkpoint`, `--resume`, `--max-new-samples`, sequential
compliance flags, and candidate-selection flags with QMC. Its library producer,
`fs_uq::QmcExecution`, supports chunked execution and in-memory clone/resume;
durable QMC checkpoint transport is not implemented. The default Monte Carlo
path and its existing recovery/confidence behavior are unchanged.

Sampling uses the midpoints of the existing 32-bit Sobol grid and the existing
approximate inverse-normal transform. Between-replicate errors do not bound
finite-grid bias, transform error, solver error, physical-model discrepancy,
or unsampled continuous-time thermal peaks. Zero replicate variation does not
prove an exact answer, nor justify optional-stopping inference.
