# Combined-load requirements on two physical grids

The finite-case refinement check is not enough for a simultaneous-load set
`f = F q, ||q||_2 <= 1`: passing every basis column leaves cross compliances
unchecked. The worst direction can also change when the physical grid changes.

`CutDensityStudy3::estimate_load_space_cap_enrichment` accepts an adaptive
study, an actually enriched `AdaptiveElasticity3`, the original reference load
laws, raw controls, explicit nominal SIMP parameters, robustness, an authored
`LoadSpaceCap3`, existing `GoalRefinementOptions3`, and a shared `SolveControl`.
It freshly solves the coarse ERODED design and inherits its physical stiffness
on the fine operator. It does not refilter fine raw controls or reoptimize them.

The existing reference-load DWR owner still performs every fine solve with the
same selected preconditioner and residual/identity checks. The new call retains
those complete fine basis fields long enough to reconstruct every cross
compliance. Old estimator APIs do not retain this extra family or change their
solve arithmetic. There are no extra pairwise-polarization state solves.

Both grids independently compute their full normalized compliance Gram matrix,
checked spectrum, maximizing coefficients, force, displacement and physical
work. `stable_and_feasible(smooth_tolerance, change_tolerance)` requires both
physical caps and smooth constraints to pass. It also checks

```text
||G_fine - G_coarse||_F <= change_tolerance,
G = F^T K_eroded^-1 F / authored_compliance_cap.
```

In exact arithmetic the Frobenius norm bounds the absolute normalized
compliance change for every unit coefficient vector. Unlike comparing two
maximum values or their diagonals, it detects changing adverse directions and
off-diagonal coupling. This computed numerical quantity is not an interval
enclosure, a remaining-error estimate, or a continuum certificate. A failed
check returns a complete assessment whose admission predicate is false, not a
silently relaxed requirement.

## Solid/void material

`estimate_binary_load_space_cap_enrichment` uses the same two-grid assessment
on a declared `BinaryComplianceOptions3` endpoint. It filters and thresholds
the coarse source, freshly solves its binary-eroded basis, and inherits the
actual solid/void stiffness multipliers. It never refilters/rethresholds the
fine grid or reuses the gray maximizing direction. The report retains the
coarse three-field masks and both grids' worst combined fields. The explicit
positive void stiffness and loads/supports in voids remain part of the model.

Both methods refuse unchanged active spaces. Reference laws, SDF, domain,
material, homogeneous supports and integration policies must represent the
same physical problem; inter-grid transfer checks compatibility, not SDF
identity. All physical work shares the supplied control. Every returned
outcome restores the borrowed study's original parameters and stiffness;
errors publish no partially checked load set. Geometry construction retains
its own caller-managed budget.

The retained basis-wise DWR masses are diagnostic, not a localization of error
for the entire load set. These APIs do not certify dilated-volume error,
connectivity, stress/buckling safety, feature sizes, zero-void behavior or global
optimality. They are library checks; existing CLI modes do not automatically
perform them or claim a two-grid pass.

```sh
cargo test --release -p fs-topopt --features cutfem-marquee \
  --lib sdf3::robust::envelope::load_space::design::limits::refinement
```
