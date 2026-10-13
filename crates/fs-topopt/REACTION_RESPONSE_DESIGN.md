# Prescribed-motion force and displacement fitting

`reaction_response_sdf3` fits a declared axial actuator-force target together
with an integral displacement target on an embedded-supported slab:

```sh
cargo run -p fs-topopt --features cutfem-marquee --release \
  --example reaction_response_sdf3 -- 40 250000 0.0012
```

Use the repository's remote execution lane where available. Arguments are
accepted updates, cumulative Krylov iterations and target reaction. The example
prints accepted history, measured reaction, displacement integral, densities,
actual optimizer stop, constraint violation and consumed work. The default
reaction is a design objective, not experimentally measured data. Iteration
limits do not imply convergence; nonlinear volume feasibility is reported
separately from augmented-Lagrangian acceptance.

## Numerical reaction

Both `CutElasticity3` and `AdaptiveElasticity3` expose `embedded_reaction`:

```
R_h = integral_Gamma h dot [sigma_s(u)n - gamma_s(u-g)]
    = -b_h(s)^T u + sum_c s_c integral_Gamma_c gamma_c h dot g
```

The sign is traction exerted **on the solid**. A constant unit-vector mode
selects a force component. `h=e_i cross (x-origin)` selects a moment component.
Zero the mode on unobserved support patches. Only the retained embedded support
is included, not any additional strongly clamped box nodes. The mode, prescribed
motion and support law must be pure reference-configuration data; numerical
surface quadrature must resolve their spatial variation.

The returned state gradient is `-b_h(s)`. Its direct cell-scale gradient holds
u fixed. For a solved response, use both that direct derivative and the
primal/adjoint terms from changes in stiffness **and** prescribed-motion lifting.
Raw stress sampling without the numerical penalty is not the same functional.
The implementation reuses the existing lifting/penalty, traces and adaptive
constraint transpose. No extra stiffness model or reaction solver is introduced.

## Response fitting API

Keep displacement observations in `ResponseCase3::targets`. Supply one slice of
`ReactionTarget3` per case to `evaluate_responses_with_reactions` or
`ProjectedResponseStudy3::new_with_reactions`. Reaction-only cases may have an
empty displacement slice. Target and scale units follow the selected force or
moment. Each trial shares one preparation across the family and one aggregate
adjoint per case, including both kinds of target. Reaction values are retained
separately in `ResponseEvaluation3::reaction_responses`.

The original displacement-only constructors and arithmetic remain available.
The same projected state retains bounds, multiplier/spectral history, work and
accepted physical fields across `run` calls. Probes restore incoming scales;
only accepted optimizer states install new scales. No dense KKT system is added.

## Reaction-aware two-grid refinement

Append `--estimate` to the example to assess the retained mixed objective under
the remaining geometry and linear-work budgets. It reports coarse and enriched
loss, reaction-offset transfer, the complete correction and original-grid marks.
The globally enriched probe is not installed or reoptimized.

Programmatic callers use `estimate_response_enrichment_with_reactions`, passing
the same reference experiment laws and parallel reaction target slices used by
the fit. Each grid reintegrates observations and reactions with its own retained
rules. Physical stiffness is parent-inherited for estimation, not produced by
refiltering raw densities. All coarse fields and the complete mixed objective
are checked before preparing either grid. A failure returns no partial estimate
and restores incoming source scales. The existing identity/Jacobi/two-level/
multilevel backends prepare once per grid and share their actions across cases.

The secant coefficient `w*(R_f+R_c-2*target)/(2*scale^2)` reconstructs the exact
quadratic loss change in exact arithmetic, including at an exact coarse fit.
Reaction is affine in displacement, so the separate constant-term change must
be retained in addition to derivative-weighted residuals. Cancellation between
these terms can be large; the reported identity defect is scaled by its actual
terms and does not certify relative accuracy of a near-zero fitted loss.

Marking uses real hierarchical residual contributions. Reaction-offset and
volume-measure transfers remain separate rather than becoming invented local
indicators. Empty marks do not certify accuracy. The original displacement-only
entry point still rejects reaction-bearing results instead of dropping targets.

These numerical two-grid differences are not continuum-force bounds, physical
actuator energy, follower loads, shape derivatives or unique identification.
The focused native reaction workflow includes the relevant Rust tests; check
its individual test results separately from dependency-bootstrap status.

## Selective adaptive refitting

Use `--refit` instead of `--estimate` for one complete fit/estimate/mark/refit
cycle with the same mixed force and displacement targets:

```sh
cargo run -p fs-topopt --features cutfem-marquee --release \
  --example reaction_response_sdf3 -- 40 250000 0.0012 --refit
```

The update cap applies separately to the coarse and fine optimizer phases. The
single cumulative Krylov allowance and geometry control are NOT reset. The
uniformly enriched error probe is discarded; only the two marked original
leaves (or fewer, plus octree balancing) are selectively rebuilt for fitting.
The marking fraction and whether its 0.5 target was met are reported even when
the two-mark cap permits only partial coverage. No marks cause an explicit
refit refusal, not an accuracy or convergence claim.

The new candidate inherits **raw** densities and then rebuilds its own filter
and projection. Actual fine-grid volume is restored with the optimizer's
specified density floor before its first physical baseline. A fresh optimizer
starts at iteration zero; coarse gradients, multipliers, spectral steps and
primal/adjoint fields are not transferred into its history. Reference force
and displacement observations are reintegrated, and reaction modes are evaluated
with current fine material on every trial. This is a new discrete design
problem, unlike the fixed-material error probe.

Fine-phase CSV blocks are prefixed `refit_`. They contain the new baseline and
accepted history, force/displacement measurements, inherited and feasible-start
raw densities, and accepted raw/projected densities. Diagnostics separate the
coarse-to-new-baseline loss jump from fine-phase optimization. A late fit error
is propagated **after** emitting the last aligned accepted fine prefix. A
successful optimizer return is still checked separately for volume feasibility;
iteration/evaluation limits are never labeled as convergence.

For programmatic use, `fit_reference_responses_with_reactions` accepts reference
cases, one reaction-target slice per case, raw densities, projected options,
a phase step allowance, and the shared solve control.
`refit_reference_responses_with_reactions` additionally takes an explicitly
built candidate `CutDensityStudy3`, source accepted evaluation, and transfer-term
budget. It borrows the source immutably and returns a separate candidate with
its initial/accepted fields, inherited/starting densities, history, work, and
`fit.outcome`. Initialization failure returns `Err`; interruption after a valid
baseline remains in `fit.outcome` with the accepted prefix. Inspect both outcome
and volume feasibility before replacing the source study.

Reaction-only cases are supported. Empty reaction rows use the same legacy
arithmetic. Missing rows, excessive total observations, invalid targets and
absent embedded supports are rejected before reference-law callbacks or solves.
The displacement-only refit entry refuses reaction-bearing evidence rather
than silently dropping force objectives. Only raw-design provenance is checked
for transfer; cached coarse objectives do not enter the new fit. Keeping SDF,
support/observation laws, target values and units consistent remains a caller
obligation, not something inferred from callback pointers.

Seven additional tests in `reaction_response_sdf3::adaptive` exercise real
marked-grid refitting, nodal/reference optimizer equivalence, reaction-only
volume restoration, input admission, cancellation and work accounting,
source/evaluation consistency, and displacement-only parity. Native execution
must be checked in the reaction workflow; source checks alone are not a Rust
build or a numerical-validation result.

## Error-driven mixed-response campaigns

`--adapt-rounds` runs the reusable library campaign rather than the one-shot
refit. The count is a maximum of **2..=4 total retained grids**, including the
initial grid. Tolerances are independent, finite and nonnegative:

```sh
cargo run -p fs-topopt --features cutfem-marquee --release \
  --example reaction_response_sdf3 -- 40 250000 0.0012 \
  --adapt-rounds 3 --response-tol 0.001 --objective-tol 0.001 --volume-tol 0.0001
```

The same 250000-iteration Krylov allowance covers every fit, restoration, probe
and rejected trial. Geometry uses one cumulative quadrature control, a 4096-leaf
background cap, and maximum level five to reserve room for the final probe.
Per-grid updates, two-mark cap and 0.5 marking target remain separate controls;
partial marking coverage is reported. Default tolerances are the values in the
command above. Do not combine this mode with `--estimate` or `--refit`.

A stop requires all of: absolute squared-loss change divided by the fixed
optimizer objective scale; **each** displacement and force/moment change divided
by its observation scale; and absolute numerical volume-fraction change. The
largest displacement and reaction discrepancies are reported separately.
Even zero-weight observations are checked. Equal squared losses on two grids
can hide opposite observation errors; a small total-loss change alone never
passes this stopping policy. This is a numerical two-grid comparison, **not a
continuum error bound**, KKT stationarity test or unique parameter-recovery test.

The campaign probes every active source leaf, then refits only selected source
leaves plus balance closure. The probe inherits physical stiffness; the refit
inherits raw controls and rebuilds the original filter/SIMP law. Protected solid
and void cells inherit labels through the checked parent map. The best feasible
accepted endpoint is retained even when the final augmented-Lagrangian iterate
has a worse objective or violates the nonlinear volume constraint. Optimizer
stops and KKT residuals still describe that final AL iterate; `selected_latest`
reports whether it is also the retained feasible design.

The library owner is `sdf3::response::adaptive::AdaptiveResponseStudy3`. Construct
it with an existing study and matching background, immutable reference cases,
optional parallel reaction slices, starting raw controls, explicit options and
a solve control. Call `run(additional_refinements, build, control)` with the same
control and a geometry factory returning an elasticity backend for the requested
`Probe` or `Refit` tree. Capture the same bounded quadrature control in that
factory. Geometry, supports, material and pure experiment laws remain fixed.

`run` caches complete probes across pauses and clears them at every promotion.
It checks the final retained grid even when the refinement allowance is reached;
a zero-refit call can therefore spend work on a probe. `TwoGridToleranceMet`,
`RefinementLimit`, `Paused` and `NoRefinementSignal` are distinct outcomes. A
failed tolerance test with no residual marks is not silently converted into an
accuracy pass: offset and measure transfers are not invented local indicators.

Before promotion, failures leave the source intact and retain any completed
candidate in `unpromoted()`. Cancellation after the atomic promotion leaves the
complete new grid, fields and history aligned. `accepted()` returns the selected
feasible endpoint; `fit()` retains its baseline, latest AL point, best feasible
point and actual optimizer outcome. `probe()` only describes the current source.
The example prints promoted-round summaries and retained fields before returning
late errors; candidate failures never overwrite the source output.

Focused native targets are `adaptive_response_campaign` (nine integration
regressions), the campaign's two policy unit tests, and the reaction example's
three argument tests. They cover manual/automatic equivalence, segmented replay,
cancelling loss changes, material masks, promotion cancellation, factory refusal,
input admission and legacy parity. They are authored but not locally executed:
DSR, RCH and Cargo are absent in the editing environment. Use the repository's
DSR lane for native verification; lexical and source checks are not compilation.
