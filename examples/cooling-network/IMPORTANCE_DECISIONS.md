# Confidence-based decisions from weighted cooling failures

Fixed `--importance EPSILON,SHIFT...` runs can use the existing explicit policy:

```sh
frankensim --json cooling-network-uq BASE.json UQ.json \
  --importance 0.2,3 \
  --compliance-probability 0.95 --confidence-alpha 0.05 \
  --min-decision-samples 32 --checkpoint first.imp --max-new-samples 20

frankensim --json cooling-network-uq BASE.json UQ.json \
  --importance 0.2,3 \
  --compliance-probability 0.95 --confidence-alpha 0.05 \
  --min-decision-samples 32 --resume first.imp --checkpoint continued.imp
```

The UQ request must declare a temperature limit and enough total samples for the
minimum. These settings, the model, and the proposal must be fixed before seeing
production results. The example is a syntax example, not a claim that 32 calls
resolve 95% compliance. The caller's lifetime sample limit is a cap, not an
allocation silently enlarged to get a favorable decision.

The event being decided is `P(numerical QoI <= temperature_limit_k) >= P_required`.
The existing `ImportanceExecution::assess_failure_probability` computes a
Gaussian-mixture confidence sequence for the complementary strict-failure event.
For a defensive proposal q=epsilon*p+(1-epsilon)*p_shift, each confidence
observation is epsilon*(p/q)*1_failure in [0,1]. The original library rescales the
result by epsilon. A proposal identical to the target uses its known unit-weight
indicator sequence. No second statistical bound is implemented by the command.

After the minimum, the run stops if the failure upper bound is at most
`1-P_required`, or its lower bound is strictly greater. The empirical proposal
failure frequency and descriptive importance standard error are NOT decision
rules. A small epsilon can make the conservative range-based confidence interval
very wide; rare-event variance reduction does not imply equally tight confidence
bounds. Zero observed failures cannot by itself certify low risk.

The new `frankensim.cooling-network-uq.importance-compliance.v1` response reports
both failure and complementary compliance confidence intervals, actual confidence
sample count, policy, proposal and physical work. `meets-probability-target` and
`below-probability-target` are resolved calculations and exit SUCCESS. Exhausted
sample, wall or chunk allowances without separation return `indeterminate` and
exit BUDGET, even when every planned sample completed. Means at a data-dependent
stopping time are not claimed unbiased. No stopped temperature distribution or
physical-approval claim is emitted.

Every accepted cooling call updates the original state before assessment. The
existing checkpoint codec binds the model and proposal; an additional versioned
binding freezes the stopping policy. Changed target, alpha, minimum, or removal
of the policy refuses before output reservation. A checkpoint at the stopping
ordinal returns the identical resolved report without new physics, including
under a tiny new evaluation allowance. Checkpoints remain trusted-input artifacts,
not authentication or proof that a physical solver ran. Recovery paths must be
new and input checkpoints are never overwritten.

A physical failure remains terminal and invalidates the current checkpoint.
A confidence-calculation error emits a diagnostic, never a fallback approval.
The per-invocation deadline limits new physical evaluations. The bounded
confidence replay and report are synchronous; they are not hard real-time work.
Conditional sampling assumptions, finite-generator/rounding limitations and
numerical-model error boundaries remain those of the existing confidence owner;
intervals are not outward-rounded certificates. Separate studies, models,
thresholds or policies need their own selection/multiplicity control.

Without the three policy flags, the existing fixed-count importance behavior and
checkpoint binding are unchanged. Other sampling/mean-control modes and design
candidate selection remain incompatible. Pilot-learned decisions are not yet
admitted in this increment.

Focused tests (native execution required):

```sh
cargo test -p fs-cli --bin frankensim uq_command::execute::importance
cargo test -p fs-cli --test cooling_importance_decision
```
