# Select cooling designs under declared uncertainty

`cooling-network-uq` can select among a finite, predeclared family of fan-speed
or workload multipliers using the existing real cooling solver, uncertainty
sampler and anytime probability assessment. Nominal temperatures and empirical
success fractions alone do not select a design.

```bash
frankensim --json cooling-network-uq \
  examples/cooling-network/fan-hotspot.json \
  examples/cooling-network/uq-candidate-hotspot.json \
  --power-candidates 0.5,0.75,1 \
  --compliance-probability 0.9 --confidence-alpha 0.05 \
  --min-decision-samples 64 --checkpoint workload-candidates.uqcp
```

For fan selection, replace the candidate flag with
`--fan-speed-candidates 0.75,1,1.25`. For a nonlinear contact-cooling trajectory:

```bash
frankensim --json cooling-network-uq \
  examples/cooling-network/nonlinear-contact-pulse.json \
  examples/cooling-network/uq-candidate-pulse.json \
  --power-candidates 0.5,0.75,1 \
  --compliance-probability 0.9 --confidence-alpha 0.05 \
  --min-decision-samples 64 --checkpoint pulse-candidates.uqcp
```

These input distributions and engineering data are illustrative declarations,
not measured data. The commands may return unresolved budgets rather than a
selection. Checkpoint destinations must be new.

## Candidate semantics

Choose exactly one candidate control and provide 2 through 64 finite, strictly
increasing multipliers. Fan values must be positive; workload values can include
zero. Fan preference is the lowest multiplier; workload preference is the
highest. This is a discrete family, not a monotonicity assumption or an
interpolation bracket. Every encountered candidate is evaluated independently.

Each parameter vector is drawn by the existing sampler first. The candidate
multiplier is then applied to that sampled request. A sampled fan speed or load
is never overwritten with an unrelated nominal value. Parameter dependence is
unchanged. All candidates reuse the declared seed and draw ordinal, allowing
common random numbers without assuming independence between designs.

For steady fan requests the factor scales `hydraulics.fan.speed_ratio`. For
transients it scales every interval's active fan speed. Steady workload
selection requires `solid.component_power`: every component's watts is scaled
and the total recomputed. A transient factor instead scales every interval's
`power_scale` or absolute `component_powers_w` map, exactly once. Base component
watts are not also multiplied. Footprints, durations, material assignments,
contact pairing, initial fields and integration rules remain unchanged.

Every candidate runs actual hydraulics, convection, solid/contact physics and
physical acceptance checks. A transient observation is the peak of a COMPLETE
trajectory, including all fixed repeated cycles, not its final temperature.
Existing forward adaptive integration and fixed-count controller models remain
whatever the base request declares; they are not differentiated. Nested design
searches, variable periodic horizons and per-sample adjoints are refused.
Nonphysical samples, bad fan domains and solver failures are errors, not hot
observations, clipped values, skipped samples or redraw opportunities.

## Confidence and selection

All three compliance flags are required. `--confidence-alpha` is the error
budget for the ENTIRE predeclared candidate family. Each existing confidence
sequence receives `next_down(alpha / candidate_count)`. The union bound then
controls the family's mathematical confidence events across all inspected
sample counts under the existing sampler/confidence-sequence assumptions.
No independence between candidates is needed, including when their draws are
shared. Changing the family or choosing new seeds after seeing results is a
new analysis, not covered by this family's allocation.

Candidates are visited in preference order. Sampling stops at a candidate's
first resolved probability decision, subject to the explicit minimum count.
A candidate still unresolved at its lifetime cap remains unresolved; the next
candidate can nevertheless be evaluated. The first qualifying candidate ends
the search, but does not erase uncertainty about better alternatives.

The `frankensim.cooling-network-uq.design.v1` result distinguishes:

- `selected`: the candidate's lower confidence bound meets the required
  probability AND every more-preferred candidate's upper bound is below it.
- `no-qualified-candidate`: every predeclared candidate is below the probability
  target. This says nothing about unlisted or intermediate designs.
- `inconclusive`: preferred alternatives remain unresolved, or the computation
  exhausted its resources. `selected_multiplier` stays null. A nonnull
  `qualified_multiplier` identifies a qualifying candidate, not a resolved
  optimum while better alternatives remain uncertain.

The first two statuses exit successfully as computations. The second is NOT
an approved design. Inconclusive results use budget exit class 6. Each row
retains its evaluated count, decision, empirical fraction and confidence bounds;
unevaluated rows are explicit. No nominal field is fabricated for a stochastic
selection. Apply the reported multiplier to the SAME load/speed semantics.

## Budgets and recovery

`samples` is the lifetime cap PER candidate. Candidate count times this cap must
not exceed 10,000. The one `wall_seconds` allowance and `--max-new-samples` count
cover the entire invocation, not a fresh budget at each candidate. The latter
requires `--checkpoint` and counts completed model evaluations across candidates.

The existing atomic checkpoint writer retains every candidate's existing
`fs-uq` envelope in one bounded file after each accepted observation. Resume
binds the exact base request, executable, UQ plan, candidate list, control,
threshold and confidence policy; individual envelopes additionally bind their
candidate index. Changed families, swapped entries, corruption and truncation
refuse. Checksums are not authenticity: resume only trusted files.

```bash
# Continue the same candidate list and compliance flags, preserving their values.
frankensim --json cooling-network-uq \
  examples/cooling-network/fan-hotspot.json \
  examples/cooling-network/uq-candidate-hotspot.json \
  --power-candidates 0.5,0.75,1 \
  --compliance-probability 0.9 --confidence-alpha 0.05 \
  --min-decision-samples 64 \
  --resume workload-candidates.uqcp --checkpoint workload-continued.uqcp
```

Resume uses a fresh evaluation-time allowance but never resets sample ordinals
or caps. Interrupted child evaluations retry the same candidate and ordinal;
a trajectory restarts from its declared initial state. A genuine model failure
invalidates this invocation's output with its diagnosis. Existing files,
including the resume source, are never replaced. A resolved checkpoint returns
the same result without another child solve. This does not add native `.fsim`
ledger/report integration or recovery inside an individual trajectory.

## Verification boundary

The implementation adds ten Rust regressions, including four actual-command
tests for fan/workload selection, complete repeated-trajectory peaks, exact
cross-candidate checkpoint/result replay, non-unit candidate scaling, unresolved
budgets, terminal errors and changed-family refusal:

```bash
cargo test -p fs-cli --bin frankensim uq_command
cargo test -p fs-cli --test cooling_uq_design
```

Those Rust tests were not executed in the authoring environment, which had no
Rust toolchain. Independent Python checks exercised exact-rational allocation,
selection tables, abstract prefix replay, high-precision confidence radii, and
finite-horizon Bernoulli boundary probabilities. They are mathematical checks,
not Rust or physical-solver execution.

Confidence concerns the declared numerical model only. The existing confidence
implementation is not outward rounded. There is no physical/model-discrepancy,
mesh-error, continuous-time peak, experimental-validation or global-design
certificate. A narrow sampling interval cannot discharge those obligations.
