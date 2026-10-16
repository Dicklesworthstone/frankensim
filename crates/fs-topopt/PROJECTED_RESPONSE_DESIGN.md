# Response design beyond the dense SQP size limit

The existing `ResponseDesignStudy3` remains the small-dense SQP path.
`ProjectedResponseStudy3` uses the new `fs-ascent::projected_al` engine for the
same response objective, one nonlinear projected-volume constraint, and exact
raw-density bounds. It does not allocate a dense Hessian, KKT system, or two
Jacobian rows per density. Optimizer vectors use O(n) storage; accepted scalar
history is bounded by the evaluation allowance. Physical fields and prepared
linear solvers retain their independent geometry and work limits.

```sh
# Original SQP behavior and original level-1 geometry:
rch exec -- cargo run -p fs-topopt --features cutfem-marquee --release \
  --example motion_response_sdf3 -- 40 250000
# Explicit projected AL; level 3 has 384 active densities on this slab:
rch exec -- cargo run -p fs-topopt --features cutfem-marquee --release \
  --example motion_response_sdf3 -- 200 250000 --projected --level 3
# The existing motion-aware response DWR is still available on an accepted fit:
rch exec -- cargo run -p fs-topopt --features cutfem-marquee --release \
  --example motion_response_sdf3 -- 200 250000 --projected --level 2 --estimate
```

Without four positional response targets the example generates disclosed
synthetic forward targets, charged to the same linear-solve budget. They are
not experimental observations or a uniquely identified material distribution.
Use DSR for the repository quality gate; RCH commands above are focused probes.

## Numerical model and stopping

Projected AL uses the existing Powell-Hestenes-Rockafellar inequality model,
with spectral projected-gradient inner steps and monotone Armijo backtracking.
It updates the multiplier only after the current box-constrained subproblem
meets its inner stationarity threshold. Spectral information is reset after a
multiplier or penalty change. There is no logistic density transformation and
no differentiation through the optimizer or linear solver.

This is a specialized box-plus-one-inequality optimizer, not general sparse
SQP or an implementation of every ALSPG algorithm. The scalar first-order
engine trades dense factorizations for potentially many additional physical
function/gradient evaluations. Removing the storage restriction does not
establish wall-clock speed, mesh-independent convergence, or a global optimum.

Bounds hold at all trial points. Nonlinear volume feasibility may fail at
accepted augmented-Lagrangian steps, and original objective values may increase
while volume violation is reduced. Check the actual stop, constraint violation,
and all KKT components. Stationarity is the box-projected mapping
`||rho - P_box(rho - grad L)||_inf`, not the unconstrained gradient norm.
Penalty, multiplier, evaluation and iteration limits never imply convergence.
A zero-step run returns cached state without changing multipliers or solving.

`objective_scale` divides objective/gradient only in optimizer coordinates;
physical objectives and responses remain unchanged. The caller must select
response scales and tolerances appropriate to the intended numerical problem.

## Physics and continuation

All samples go through the original `evaluate_responses`: current-density
prescribed-motion lifting, full Nitsche/ghost bilinear derivatives, unchanged
SIMP/projection/filter pullback, and shared primal/adjoint preparation. Each
trial restores incoming physical scales. Only a fully accepted optimizer point
installs new scales and its corresponding fields. Cancelled and failed probes
retain spent work. Complete dual updates and the spectral step survive split
runs; an interrupted line search may repeat probes, without refunding work.

The same accepted `ResponseEvaluation3` feeds response DWR. No motion work is
reinterpreted as external traction, reaction force, or actuator work. Native
Rust compilation and tests remain pending in this delivery environment; the
independent NumPy/SciPy development checks do not execute repository kernels.
