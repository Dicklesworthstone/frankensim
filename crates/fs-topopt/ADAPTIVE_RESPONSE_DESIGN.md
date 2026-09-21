# Goal-driven response refitting

The existing motion-response example can now execute the whole
**fit → estimate → mark → selectively refine → restore volume → refit** loop.
It uses the same reference forces, prescribed boundary motions, response
observations and target values on every grid.

```sh
rch exec -- cargo run -p fs-topopt --features cutfem-marquee --release \
  --example motion_response_sdf3 -- 200 250000 --projected --adapt-rounds 2
```

`--adapt-rounds` counts the initial grid and admits 2 through 4 grids. It requires
explicit `--projected`; ordinary SQP, `--projected` without adaptation, and the
old `--estimate` path retain their existing behavior. `--level 1..3` selects the
initial grid. An additional `--estimate` checks the final retained design too.
Geometry caps (4 levels, 4096 leaves and the quadrature allowances) still apply;
an unavailable enrichment refuses rather than silently increasing those limits.

The first positional budget is accepted steps **per grid**. Optimizer callback
limits also apply per new fit. The second budget is cumulative Krylov work for
synthetic target generation, every fit, all restorations and all enriched
primal/adjoint solves. Neither the solve controller nor the geometry controller
is reset between rounds. Setup work is separately accumulated as before. These
are work categories, not equivalent flops or wall-clock measurements.

Four optional target scalars follow the two budgets. Without those scalars, the
initial grid generates disclosed synthetic targets once. Those values are frozen
for subsequent grids, not regenerated into an easier new fitting problem.
Synthetic fitting is not experimental calibration or unique material recovery.

## What gets transferred

A global one-level refinement is an **estimation probe only**. The existing
response DWR estimator selects up to two cells from the original coarse tree;
only those marks and balancing constraints define the next retained grid.
Signed corrections, achieved marking fraction and a possibly unmet 50% target
are printed explicitly. Zero signal is not an accuracy certificate.

The next fit inherits **raw densities only**, checks its own projected material
volume using its new filter and cut measures, and restores feasibility toward
the optimizer's declared density floor. It then integrates fresh nodal forces
and observations and solves a new baseline. Old displacements, adjoints,
gradients, multipliers and spectral optimizer state never become fine-grid data.
The baseline jump is reported separately from optimization; no monotonic
improvement is asserted across different discretizations.

## Library API and failure behavior

`CutDensityStudy3::fit_reference_responses` binds `ReferenceResponseCase3` laws to
one grid and runs the existing `ProjectedResponseStudy3`. A failure during
initialization returns an error without installing material. Once initialization
succeeds, the returned fit always retains the initial and last accepted fields;
a later optimizer error is in its `outcome`, not hidden as success.

`source.refit_reference_responses(candidate, accepted, cases, options, steps,
max_transfer_terms, control)` consumes a separately constructed candidate
pipeline. It never mutates the source. The result contains inherited raw values,
restored starting values, the candidate study and its complete fit outcome.
The caller must inspect that outcome and volume feasibility before promoting it.
Same-grid continuation with retained optimizer state still uses the original
`ProjectedResponseStudy3`; a refit intentionally starts fresh on a different grid.

The example promotes only completed, volume-feasible candidates. A refused
transfer, unavailable probe, exhausted solve budget, failed candidate fit or
infeasible final candidate leaves the previous design in place, prints its
retained grid keys, densities and responses, and exits nonzero. Partial candidate
work remains charged. An iteration-limited but feasible fit can be promoted;
its actual numerical stopping reason is printed, not called converged.

## Verification scope

Focused targets are `adaptive_response_sdf3` and the tests of the
`motion_response_sdf3` example, alongside the existing projected-response and
response-goal tests. Native compilation and execution are pending in the authoring
environment. Numerical two-grid estimates and finite quadrature are not rigorous
continuum bounds, measured physical validation or a mesh-independent performance
claim. Prescribed motion is on the fixed reference surface, not a follower-load
or shape derivative.
