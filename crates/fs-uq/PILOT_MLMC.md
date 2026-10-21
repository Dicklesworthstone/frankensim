# Pilot-allocated multilevel physical simulation

`PilotMlmcExecution` chooses a fixed production allocation from a complete,
separate pilot. It wraps `MlmcExecution`; there is no second fine/coarse sampler,
statistics implementation, or model-call state machine.

```rust
use fs_uq::{CorrelationModel, ParameterUncertainty, PilotMlmcConfig,
    PilotMlmcExecution, PropagationMethod, UqPlan, UqStatus};

let plan = UqPlan::new("temperature at one second", PropagationMethod::MultilevelMonteCarlo, 8192)
    .with_parameter(ParameterUncertainty::uniform("initial temperature", 0.0, 1.0, "K"))
    .with_correlation(CorrelationModel::Independent);
let config = PilotMlmcConfig {
    pilot_samples_per_level: 16,
    physical_model_costs: vec![1.0, 2.0, 4.0],
    target_standard_error: 0.025,
    minimum_production_samples_per_level: 4,
    correction_std_dev_floor: 1e-4,
};
let mut execution = PilotMlmcExecution::new(&plan, config)?;
let report = execution.advance(usize::MAX, || false, |level, parameters| {
    let steps = 1usize << level;
    let factor = (1.0_f64 / (1.0 + 1.0 / steps as f64)).powi(steps as i32);
    Ok::<f64, &'static str>(parameters[0] * factor)
});
assert_eq!(report.status, UqStatus::Complete);
assert!(report.model_calls_evaluated <= plan.budget_max_samples);
assert!(report.estimate.is_some());
# Ok::<(), Box<dyn std::error::Error>>(())
```

The example is backward-Euler cooling, not an accuracy qualification of a
particular cooling product. Its exact finest-grid expectation is
`0.5 * (1 / 1.25)^4`; the continuum expectation is a different quantity.

## Allocation and evidence

Costs price one physical solve at each resolution. The allocator prices level
zero at `cost[0]` and every higher correction at `cost[l] + cost[l-1]`. Its
continuous sample counts follow `sqrt(V_l/C_l) * sum(sqrt(V_j*C_j)) / target^2`,
rounded up and bounded below by the declared production minimum. The positive
standard-deviation floor is expressed in QoI units and affects allocation only.
A zero pilot variance is not treated as proof that the model is deterministic.

The original plan budget covers both the full pilot and all production solves;
a higher-level correction always costs two calls. An allocation that exceeds
that cap refuses before production rather than silently reducing its counts.
The modeled-cost allocation is not an exact integer optimum under a separate
physical-call constraint: refusal is not a proof that every alternative
allocation is infeasible. Common cost rescaling and QoI unit changes preserve
the intended allocation. Floating-point borderline counts can differ by one.

Pilot and production use the canonical Philox sampler with different keys:
`pilot_seed = plan.seed`; `production_seed = plan.seed XOR 0x4d4c4d4350524f44`.
Both are reported. Production ordinals start at zero under the production key.
These addresses are disjoint within the run; independence remains the existing
pseudorandom-stream model. Cross-run callers must avoid reusing either key.
The selected production plan and allocation are available through read-only
accessors, so a standalone `MlmcExecution` reproduces the production result.

Only production observations enter `report.estimate` and `report.compliance`.
A complete pilot report remains visible as diagnostics. Production allocation
is immutable, even when its realized standard error misses the target.
`observed_sampling_target_met` reports that descriptive comparison separately;
`Complete` means the selected design completed, not that physical accuracy was
proved. Compliance uses physical threshold indicators but has no separate
allocation target. No bias extrapolation, level selection, confidence bound,
optional-stopping guarantee, or physical-model validation is supplied.

## Continuation and failures

`advance` spends a caller allowance across both phases. Zero allowance is a
no-op. Fine/coarse pairs retain identical physical inputs, and a paid fine
result is never repeated because its coarse partner was interrupted. Cloning
preserves both execution owners, their observations, and the selected design.
`advance_interruptible` treats `Ok(None)` as unfinished work and retries the
same physical call; that solver's interrupted internal cost remains its own.

Cancellation brackets allocation and every physical call. Resuming after the
last pilot call does not re-evaluate the pilot. Errors and nonfinite physical
outputs permanently refuse the experiment; no draw is skipped or replaced.
Long callbacks must implement their own cancellation and preserve a fixed
physical meaning throughout the entire experiment. There is no durable
checkpoint codec in this producer.

`tests/product_mlmc_allocation.rs` specifies analytical cost/count checks,
unit/cost rescaling, floors, numerical extremes and budget refusals.
`tests/product_mlmc_pilot.rs` specifies a cooling hierarchy, independent
production reconstruction, paired-call continuation, phase-boundary
cancellation, terminal failures, singular Gaussian admission and a deliberately
missed observed target. These Rust tests require native execution; independent
Python numerical checks are not substitutes for that execution.
