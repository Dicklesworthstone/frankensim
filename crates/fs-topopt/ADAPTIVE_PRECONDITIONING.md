# Adaptive optimization and enrichment preconditioning

`CutDensityStudy3` uses one numerical preparation per evaluated density and
shares it across independent loads. Goal refinement separately prepares once
for its enriched physical design. Geometry stays fixed inside each optimization
call; numeric factors are never reused after a density change. Bare library
backends retain identity preconditioning.

```sh
# Existing defaults: two-level enrichment, Jacobi optimization, initial level 1.
cargo run --locked -p fs-topopt --features cutfem-marquee --release \
  --example adaptive_elastic_sdf3 -- 2 3 250000
# Recursive correction for BOTH optimization and enriched goal solves:
cargo run --locked -p fs-topopt --features cutfem-marquee --release \
  --example adaptive_elastic_sdf3 -- 2 3 250000 multilevel multilevel 2
# Start at level 3 to exercise an enriched first coarse space above 512 DOFs.
cargo run --locked -p fs-topopt --features cutfem-marquee --release \
  --example adaptive_elastic_sdf3 -- 2 2 250000 multilevel multilevel 3
```

Use the repository's DSR/RCH execution lane where available. Arguments are
rounds, updates per round, cumulative outer iteration limit, enrichment mode,
optimization mode, and initial uniform level. Enrichment modes are identity,
Jacobi, two-level and multilevel; optimization modes are Jacobi, two-level and
multilevel. All size, geometry and numeric-work refusals remain explicit.

## Recursive correction

Construct `AdaptiveSolveSpace3::multilevel(fine, &coarser, options, checkpoint)`
from the existing CutFEM preconditioning module. The borrowed coarse geometries
are ordered nearest-to-farthest; construction retains their sparse Q1 transfer,
not their stiffness matrices. Every adjacent pair must satisfy the existing
box, reference material, clamp, active-support and refinement admission.
The example retains uniform lower-resolution spaces across refinement rounds;
other callers may supply valid locally refined nested spaces.

The finest action is additive: `D_f^-1 + P_f V P_f^T`. Intermediate `V` levels
use forward Gauss-Seidel, recursive residual correction, and the backward
adjoint sweep. Only the final level is factored densely. The finest operator
remains matrix-free; this is a hybrid recursive preconditioner, not a claim
that every level is matrix-free or a full finest-level multiplicative V-cycle.

The physics owner contracts the actual constrained element blocks and
current-density ghost jumps into `P_f^T A_f P_f`. It does not probe every fine
coordinate, assemble the full fine matrix, re-integrate geometry or substitute
a rediscretized coarse stiffness. Deeper matrices are sparse Galerkin products.
Intermediate dimensions can exceed 512; only the bottom factor has that hard
ceiling (default 192). Entry, transfer, level and accumulated-product caps bound
setup structures and arithmetic categories, not peak RSS or wall-clock time.

`prepare_with_work` rebuilds the diagonal, sparse matrices and bottom factor at
every evaluated density, including rejected trials. Prepared values borrow the
fine operator and are fixed throughout each outer CG solve. Symmetric smoothing
and a direct bottom inverse avoid a residual-dependent nonlinear inner solve.
Rank/positivity failures refuse without pivot shifts or silent fallback.

## Goal refinement

Select `GoalPreconditioner3::Multilevel` and call
`estimate_compliance_enrichment_with_coarse_levels`. The accepted grid is always
the first correction space for the enriched problem; pass only additional
geometries BELOW it, nearest first. With no extra levels, the original enrichment
entry point also admits a single sparse bottom when it fits the cap. Other
policies reject extra ladders instead of ignoring them.

The complete accepted coarse load family is checked before numerical setup.
Physical stiffness inheritance, independent fine solves, actual-field residual
checks and the DWR decomposition remain shared by all solve policies. Refinement
marks apply to the original coarse tree, not the globally enriched probe.
Two-grid differences are not continuum bounds, and zero marking signal is not
an accuracy certificate.

## Work and limits

`SolveWork` separates `linear_iterations`,
`preconditioner_operator_applications` and `preconditioner_galerkin_products`.
The latter counts admitted sparse upper-triangle summands, including discarded
and cancelled construction. Do not add those three numbers as equivalent work
units. Two-level fine applications and recursive local products are different
setup mechanisms. Sparse setup has a whole-hierarchy per-preparation budget;
outer iterations retain their existing cumulative campaign budget. Callbacks
can enforce an additional campaign-level policy from the cumulative counters.

Original identity, exact Jacobi and dense two-level paths remain available.
Recursive cycles allocate temporary vectors and are not preemptible within a
single apply. No allocation-free, hard-real-time, mesh-independent convergence,
asymptotic scalability or continuum certificate is claimed. Native Rust tests
remain pending; independent numerical references are not Rust execution.
