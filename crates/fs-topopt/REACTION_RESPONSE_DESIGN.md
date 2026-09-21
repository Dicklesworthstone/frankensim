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

Reaction-aware DWR and adaptive reaction refitting are not supplied here.
`estimate_response_enrichment` rejects nonempty reaction observations instead
of returning a displacement-only error estimate for a mixed objective. This
numerical support-traction response is not a continuum-force bound, physical
actuator energy, follower load, shape derivative or unique material identification.
Native Rust compilation/tests remain unverified in the development environment;
independent NumPy checks are not execution of the Rust physics or optimizer.
