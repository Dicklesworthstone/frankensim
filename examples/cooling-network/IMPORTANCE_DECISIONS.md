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
candidate selection remain incompatible.

## Learned proposals: decisions use independent production only

The same three policy flags also work with `--importance-pilot`. Every declared
pilot must finish and its last proposal must freeze before any response can
enter the confidence sequence. The minimum counts production observations, not
the total physical call count. Admission refuses when the production allowance
left after all pilots is smaller than the minimum; it neither skips pilots nor
enlarges the original budget.

```sh
frankensim --json cooling-network-uq BASE.json UQ.json \
  --importance-pilot 2,8,0.2,0.25,1,8 \
  --compliance-probability 0.95 --confidence-alpha 0.05 \
  --min-decision-samples 32 --checkpoint learning.imp --max-new-samples 10
```

This syntax needs at least 48 total calls in the UQ document: 16 for the pilots
and at least 32 for production. It does not promise resolution in that budget.
Repeat the learning and confidence settings with `--resume learning.imp` and a
fresh checkpoint destination to continue. A prefix still inside a pilot reports
zero confidence samples and null confidence intervals, even when its observed
pilot failures are numerous. Only final-production outcomes may stop the run.

The `importance-pilot-compliance.v1` report separates pilot and production work,
retains learned proposals and phase seeds, and uses the same failure-CS owner as
the fixed-proposal command. Conditional on the completed pilots, production has
one fixed proposal, independent sampling stream, model and event; the original
bounded-mean assumptions apply to that production stream. Pilots are never
pooled into the confidence sequence or its minimum. The declared alpha is not
spent on separate pilot decisions because none are made.

The existing complete-study checkpoint retains learned history and the current
phase. Its additional policy identity prevents changing alpha, target or minimum
after inspecting production results. Restarting across a fit/production boundary
retains paid calls; restarting a resolved decision makes no further physical
calls. The original fixed-count pilot command, without confidence flags, remains
unchanged. No new sampler, proposal update, confidence formula or checkpoint
codec is introduced by this integration.

Focused tests (native execution required):

```sh
cargo test -p fs-cli --bin frankensim uq_command::execute::importance
cargo test -p fs-cli --test cooling_importance_decision
cargo test -p fs-cli --test cooling_importance_pilot
```
