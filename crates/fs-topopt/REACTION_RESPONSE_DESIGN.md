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
Automatic reaction-aware adaptive refitting is not yet connected.

These numerical two-grid differences are not continuum-force bounds, physical
actuator energy, follower loads, shape derivatives or unique identification.
The focused native reaction workflow includes the relevant Rust tests; check
its individual test results separately from dependency-bootstrap status.
