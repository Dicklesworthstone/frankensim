# Worst-case propagation of bounded parameter uncertainty

`EpistemicExecution` implements the method-specific `EpistemicBounding` path.
The single-level Monte Carlo `UqPropagator` still refuses that method: point
samples cannot enclose a continuum of parameter values.

## What the evaluator promises

The callback takes a closed parameter box in the original plan order and
returns `[lower, upper]` enclosing the value of one fixed real-valued model at
**every point** of that box. Singleton boxes use the same contract: their
response may still have nonzero enclosure width. For example, a caller can
convert the inputs to `fs_ivl::Interval` and evaluate an outward interval model.
A point solver needs its own valid error bound and parameter-range enclosure
before it can serve as this evaluator. Merely checking corners is not enough:
`1 - x*x` on `[-1, 1]` reaches its maximum at the interior point zero.

The producer preserves inclusion if the evaluator fulfills that contract. It
cannot establish an arbitrary callback's soundness, validate the physical
model, or mint Verified evidence. In particular, the returned values are not a
probability distribution, confidence interval, or empirical sample range.
No floating-point point estimate is promoted to an exact physical response.

The input domain is explicitly the full Cartesian product of the supplied
finite epistemic intervals. Unknown stochastic dependence does not become
probabilistic independence. Extra joint feasibility constraints are not modeled;
a retained parameter witness need not satisfy an unrepresented constraint.
Probability laws and correlation matrices refuse rather than being converted
into interval domains. Names and units are retained; dimensional compatibility
remains the caller's obligation. The original seed is retained but unused.

## Bounds, witnesses, and continuation

Each accepted leaf encloses its entire subdomain. Midpoint singleton evaluations
supply actual parameter witnesses. The minimum is bracketed between the least
leaf lower bound and the best witness upper bound; the maximum is bracketed
between the best witness lower bound and the greatest leaf upper bound.
These bound the infimum and supremum without assuming continuity or attainment.
`Complete` requires both brackets to meet `absolute_tolerance`, in QoI units.
The actual range of the response can be much wider than that tolerance.

Refinement selects the cell that could most improve either extremum. Axes are
chosen by the fewest previous halvings, not by comparing numerical lengths with
incompatible physical units. Priority arithmetic does not enter any bound.
Children cover the parent's two closed halves, and their enclosures are
intersected with the parent's. The parent stays live until BOTH children and
their midpoint observations are accepted. Thus a partial split cannot remove
unvisited possible extrema. Reported outer ranges never widen under a sound
fixed evaluator. Disjoint nested enclosures refuse explicitly.

Cells whose entire output ranges lie between the best minimum witness upper
bound and maximum witness lower bound cannot improve either extremum. They are
retired before requesting more live-cell storage. An aggregate hull retains
their complete inclusion, so pruning frees parameter-box storage without
discarding domain coverage. `live_cells` counts retained boxes; `retired_cells`
counts these resolved interior domains; `covered_cells` includes both. Witness
improvements can only strengthen the retirement criterion. Cancellation during
pruning preserves the accumulated coverage and deterministic continuation.

The lifetime plan budget counts both box and singleton calls, capped at 65536.
Storage holds at most 4096 live leaves over at most 64 parameters, plus two
prospective children and two witnesses. Zero-width axes are admitted; adjacent
binary endpoints are not split into a fictitious interior. Unrepresentable
width arithmetic becomes an unresolved width, never a spurious completion.

`advance` consumes an additional call allowance. `advance_interruptible` also
accepts `Ok(None)` for an unfinished callback, which is retried on the identical
box. Paid child and point results survive chunking, cancellation and cloning.
Cancellation is polled between calls, between cell-selection visits, and before
completion. The solver inside a long callback owns its own cancellation.
`Cancelled` may expose the earlier valid full coverage, unlike restart-only
geometric APIs. Refused executions expose no numerical bounds or witnesses.

Evaluation exhaustion, cell exhaustion and parameter-resolution exhaustion are
separate terminal statuses, not successful tolerance closure. Zero allowance
and calls on terminal executions do nothing. A callback error or malformed
interval is terminal and never becomes a zero response or a skipped subdomain.
There is no durable checkpoint codec, portable certificate, derivative-based
contraction, or implicit numerical-model verification.

## Decision-driven upper-limit screening

With `plan.compliance_threshold` set, every report distinguishes:

- `AllCompliant`: the full-domain output upper bound is at or below the limit.
- `AllViolating`: the full-domain output lower bound is strictly above it.
- `WitnessedViolation`: an actual parameter witness has its entire response
  enclosure strictly above the limit, while a universal verdict is unavailable.
- `Undetermined`: none of those implications is established; a loose upper bound
  above the limit is not by itself a violation.

`advance_until_compliance` and its interruptible counterpart stop immediately
when one of the first three decisions is established. A broad but sufficient
root enclosure may decide compliance after one call, without an unnecessary
midpoint call or tolerance-level extremum search. Equality is compliant for the
declared `QoI <= threshold` convention. No frequency or probability is inferred.

The status is `DecisionReached`, distinct from `Complete`. Ordinary `advance`
can resume that exact paid state to refine both extrema; it does not restart
the root or repeat a completed half-split. A missing threshold returns an error
without altering the execution. Cancellation precedes decision publication;
the cancelled report may still carry the earlier valid bounds and their logical
decision. Resumption can publish that decision without repeating physical work.
Conversely, tolerance completion can leave a near-threshold decision unresolved.

## Example

For a linear thermal response `T = ambient + power * resistance`, positive
parameter intervals give a simple outward enclosure:

```rust
use fs_uq::{EpistemicConfig, EpistemicExecution, ParameterUncertainty,
    PropagationMethod, UqPlan};

let plan = UqPlan::new("junction temperature", PropagationMethod::EpistemicBounding, 16384)
    .with_parameter(ParameterUncertainty::interval("ambient", 20.0, 25.0, "K"))
    .with_parameter(ParameterUncertainty::interval("power", 2.0, 4.0, "W"))
    .with_parameter(ParameterUncertainty::interval("resistance", 1.0, 2.0, "K/W"));
let mut execution = EpistemicExecution::new(&plan, EpistemicConfig {
    absolute_tolerance: 0.05,
    maximum_cells: 4096,
})?;
let report = execution.advance(usize::MAX, || false, |b| {
    let product_lo = (b[1][0] * b[2][0]).next_down();
    let product_hi = (b[1][1] * b[2][1]).next_up();
    Ok::<_, &'static str>([
        (b[0][0] + product_lo).next_down(),
        (b[0][1] + product_hi).next_up(),
    ])
});
let lower_extremum = report.minimum;
let upper_extremum = report.maximum;
# Ok::<(), &'static str>(())
```

The analytical extrema of this declared model are 22 K and 33 K. This example
qualifies neither a particular thermal solver nor an actual cooling product.
Focused regressions are in `tests/product_epistemic.rs`; they require native
Rust execution. Independent numerical reference checks are not Rust test passes.
