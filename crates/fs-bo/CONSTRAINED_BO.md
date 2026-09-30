# Noisy design optimization with outcome constraints

`constrained::minimize_constrained` connects the existing Gaussian processes,
heteroscedastic kernel learner and CMA-ES engine to a multi-output design loop.
Each physical callback returns an objective plus every declared constraint,
each with its own observation-noise variance. All values remain in their
original units; signal and noise variances use squared output units.

```bash
cargo run -p fs-bo --example constrained_cooling
cargo test -p fs-bo --lib constrained
cargo test -p fs-bo --all-targets
```

The example uses explicitly synthetic algebraic fan-power, temperature and
sound responses, not an actual thermal/acoustic solver. Its noiseless feasible
speed interval is [0.5, 5/7], while the unconstrained power minimum is at zero.
On success it makes six initial calls and four two-point batches: 14 complete
multi-output evaluations. Objective and temperature kernels learn after
initialization and completed batches 2 and 4, each with at most 24 likelihood
probes per refit. The sound kernel stays fixed. This implies at most 144
likelihood probes and nine ordinary conditioning fits; acquisition posterior
factorizations are separate work. No measured Rust example result is claimed.

## Joint feasibility rather than an objective penalty

For each posterior draw, the acquisition uses the best feasible utility

```text
U(S) = max({0} union {reference - objective(x): x in S and all constraints pass})
acquisition = average[U(baseline union candidates) - U(baseline)]
```

The finite `reference` is a caller-declared zero-utility objective value, not an
estimated incumbent or a bound inferred from observations. It defines the case
where no baseline point is feasible and gives no reward to objectives above
it. Pick a reference that covers objective values worth exploring: an overly
optimistic reference can suppress acquisition. This reference-capped utility
is not identical to ordinary uncapped q-NEI on every posterior draw.

Each output has an independent GP. Baseline and candidate values for that
output are sampled together, retaining their spatial correlation. Every
constraint must pass on the SAME draw at a point; neither a posterior-mean
feasibility filter nor the best raw noisy observation supplies the incumbent.
An infeasible low-objective observation cannot suppress feasible improvement.
Equal coordinates, including signed zero, share a single latent draw per
output. Repeated noisy measurements remain distinct likelihood observations.
Duplicated zero-noise constraints refuse instead of being repaired by adding
undeclared noise.

The bank is point-major, output-minor within each sample row: objective first,
then constraints. Row-wise prefixes preserve common random numbers while a
batch grows. The existing `joint_normal_bank` supplies Sobol columns and an
explicit Philox Monte Carlo tail beyond the embedded Sobol dimension ceiling.
The existing joint posterior's jitter and degenerate-covariance fallback are
inherited, not strengthened into a covariance certificate.

This follows the general sampled feasible-utility approach to outcome
constraints; see BoTorch's outcome-constraint discussion and the Letham et al.
noisy-constraint treatment for context. The finite reference, inclusive upper
bounds, independent-output assumption and deterministic bank grammar above
are the actual contract of this implementation.

- https://botorch.org/docs/next/constraints
- https://arxiv.org/abs/1706.07094

## Physical units, learning and recommendations

`ConstrainedBoConfig.search` supplies the initial objective kernel, prior mean,
box, seed, batch and acquisition-search limits. Each `OutcomeConstraint`
supplies its own kernel, prior mean and physical upper bound. The driver
subtracts the appropriate fixed mean from each observed outcome and bound.
For a lower limit, negate that outcome and its limit together; its variance
is unchanged. Output cross-covariance and correlated measurement noise are
not modeled.

`objective_learning` and each constraint's `learning` optionally declare an
independent `NoisyLearningConfig`. Each scheduled refit warm-starts from its
own previous selected kernel and retains the original observation variances.
Between refits, every output still conditions on ALL new observations. The
winning learner model is reused without an extra unbudgeted fit. Records name
the output, data prefix, kernel, resolved seed, local likelihood gain and work.
Likelihoods are only compared within the same data prefix. Exhaustion does not
imply convergence or trigger additional physical callbacks.

After initialization and each batch, the driver returns the lowest objective
posterior mean among evaluated points whose modeled joint latent feasibility
probability reaches `recommendation_probability`. The probability is the
product of independent Gaussian marginal probabilities using the existing
approximate normal CDF. `None` explicitly means that no evaluated point meets
that model policy. It is not a proof of global infeasibility. Recommendations
can exceed the acquisition's objective reference; the two policies are distinct.

**This is not safe Bayesian optimization.** The probability policy filters
recommendations, not all proposed evaluations. Candidate evaluations may
violate constraints. Hard physical safety requires separate authoritative
admission. Neither these GP probabilities nor the Monte Carlo acquisition are
frequentist confidence bounds, physical validation or feasibility certificates.
The driver is synchronous and dense; it does not add sparse-history scaling,
within-factorization/callback cancellation, correlated outputs, replicate
value-of-information, statistical stopping or cross-ISA guarantees.

## Verification boundary

Eight acquisition and eight closed-loop Rust regressions cover finite-draw
feasible utilities, simultaneous constraints, absent feasible incumbents,
duplicate identity, covariance, centering, complete-data conditioning,
recommendations, independent refit cadence, learning/replay and refusals.
These Rust tests and the example remain unexecuted in the authoring environment:
DSR, RCH, Cargo and rustc are absent. Independent Python references passed 2,000
finite-draw utility identities, and a sampled feasible Gaussian utility agreed
with its analytical value within 1.2e-7. Those equation checks do not establish
that the Rust sources compile or pass their tests.
