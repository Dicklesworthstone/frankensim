# Constrain sampled thermal hot spots

`ScanPeakResponse` selects a strictly increasing, nonempty list of mesh vertex
IDs and an inclusive range of retained endpoints. It uses kelvin for temperature
or J/kg for specific enthalpy. Supply all vertex IDs explicitly to cover the
whole nodal field; empty never means "everywhere".

For selected values `y` and an explicit positive smoothing scale `tau`, the
response is `max(y) + tau * ln(sum(exp((y - max(y)) / tau)))`. In exact arithmetic
it is at least the sampled maximum and at most `tau * ln(N)` above it. This uses
a sum, NOT an average: a log-mean-exponential can underestimate the hottest
sample. Reports keep the actual sampled maximum and added smoothing margin
separate from physical solver tolerances. Ordinary binary64 arithmetic is not
an outward-rounded interval certificate.

## Physical inverse-design API

`fs_conduction::bc::history::march::enthalpy::scanning::sensitivity::objective::responses::peaks`
exports `evaluate_scan_objective_with_peaks`, `ScanPeakResponse`, and
`ScanPeakEvaluation`. Existing point responses and the fitting objective share
the same forward tape. Each peak contributes one reverse sweep through its
last endpoint, even when it contains thousands of selected nodal/time samples.
The reverse pass injects local softmax seeds at all selected times and retains
earlier-heating, initial-enthalpy, moving-source, nonlinear-conductivity and
radiation derivatives. On a latent plateau, a temperature peak can have zero
sensitivity while an enthalpy peak still changes.

`fs_ascent::scan_design::constraints` exports `ScanPeakBound` and
`ConstrainedScanStudy::new_with_peak_bounds`. Pass ordinary control bounds,
energy caps, optional point-state bounds, peak bounds, a combined
`ScanResponseBudget`, a total `max_peak_samples`, and the dense KKT dimension cap.
Peak bounds are upper limits only; lower-bounding a maximum cannot establish
that every sample exceeds a minimum. Point-state lower bounds remain available.

For a two-endpoint, eight-vertex enthalpy cap:

```rust
let peak = ScanPeakBound {
    response: ScanPeakResponse {
        first_endpoint: 0,
        last_endpoint: 1,
        vertices: (0..8).collect(),
        observable: ScanObservable::SpecificEnthalpy,
        smoothing_scale: 0.05,
    },
    maximum: 155.0,
    scale: 10.0,
};
```

This single SQP row needs two additional physical transpose calls, not sixteen
separate reverse sweeps. The response budget counts objective plus point/peak
reverse calls; retained gradient counts include original source controls and
initial-state gradients. `max_peak_samples` bounds the sum of selected samples
across all peak responses before the forward trajectory. It is not a byte or
peak-allocation certificate. The original tape storage budget still applies.

All optimizer-coordinate derivatives use the existing signed-centre/log-width/
log-power chain. Unselected off values remain zero. `replay_with_constraints`
returns the objective and every constraint from one trajectory, including
`peak_values`, `sampled_peak_values`, and `peak_smoothing_gaps`. Input-only
`constraints_at` refuses attached peak limits rather than omitting them.

## Executable synthetic example

```bash
cargo run --manifest-path examples/scan_power_fit/Cargo.toml -- --peak-limits
cargo run --manifest-path examples/scan_power_fit/Cargo.toml -- --peak-limits 155 0.2
cargo test --manifest-path examples/scan_power_fit/Cargo.toml -- peak_limits
cargo test -p fs-conduction --test scan_response_rows
```

An off-centre source heats a unit cube initially on a 310 K melting plateau.
Targets request 152 and 157 J/kg at one relatively cold vertex. One peak bound
covers all eight vertices at both endpoints, so the fit cannot meet that target
by ignoring a hotter node elsewhere. CSV exports the fitted held power schedule;
comments report both envelope and actual sampled maximum. An unsuccessful KKT
or feasibility result exits unsuccessfully. The example also exercises budget
pause/resume, source export and final physical replay.

The limit and smoothing bias have different roles. A smaller tau reduces the
bias but can increase local optimization curvature; the implementation does not
alter tau during optimization. Refining the sample set changes N and the
function being constrained. Select tau and the spatial/time samples explicitly.

These are fixed-grid local design limits, not unsampled temporal peaks, a bound
on the continuum solution, experimental calibration, or machine safety. SQP may
evaluate infeasible intermediate candidates. Material/phase corners and physical
solver failures retain their original refusals; the peak adapter never widens
those gates or modifies the forward temperature/enthalpy fields.
