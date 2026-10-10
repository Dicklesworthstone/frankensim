# Multilevel solvers in native 3-D design studies

Native `fsim-sdf3-study` files can select the existing sparse recursive Galerkin
solver for optimization and goal-refinement solves. This connects a reusable
geometric hierarchy to each newly evaluated physical density, including the
stress adjoints and accepted/incumbent endpoint checks during stress recovery.

```sh
cargo run --release -p fs-cli --features sdf3-study --bin frankensim -- \
  --json study examples/marquee/bracket-3d-multilevel.fsim multilevel-study.db
```

The complete example has two independent body loads, two continuation stages,
and an initial level-2 background. Its material and load numbers are a numerical
demonstration, not an engineering specification or a physical validation.

## Explicit declaration

Append the following optional section after the optimizer and, when present,
after `design-regions`. Fields have the shown order; duplicates and unknown
fields are refused rather than ignored.

```lisp
(linear-solver
  :type multilevel
  :coarsest-level 0
  :max-transfer-terms 500000
  :max-matrix-entries 500000
  :max-galerkin-products 100000000
  :max-diagonal-contributions 100000000)
```

`coarsest-level` must be below `initial-level`. Uniform correction geometries
are built once, nearest first, through the SAME native domain/material/support
builder. With initial level 2 and coarsest level 0 this gives level 1 then level
0 beneath the optimization grid. An enriched goal solve also prepends the
CURRENT accepted grid to that ladder. The final compact dense factor is bounded
at 512 coordinates; intermediate grids are sparse and do not need to fit that
bottom limit. The existing initial-level, leaf, quadrature and fine-DOF envelopes
are unchanged. This is not an unlimited-grid mode.

Geometry-only corrections determine interpolation, never a replacement physical
stiffness. Numeric preparations use the CURRENT fine-density bulk, Nitsche and
ghost terms. One preparation is shared across a compliance evaluation's
independent right-hand sides. No factors from an older density are reused.
The existing true-residual and numerical gradient gates remain in force.

## Budgets, interruption and recovery

Transfer, sparse-matrix, Galerkin-product and diagonal-contribution limits govern
each bounded hierarchy construction/preparation, not the sum of all work over
an optimization. The existing report's `preconditioner_galerkin_products` and
`preconditioner_operator_applications` record cumulative setup work separately
from `linear_iterations`; these cost units are not interchangeable. Galerkin
products include unsuccessful preparations. The study's original wall and
Krylov controls still cover all evaluations and rejected trials.

Correction geometry consumes the original cumulative quadrature allowance.
Memory admission adds an explicit hierarchy/transfer/geometry envelope to the
existing field estimate; it is not measured peak RSS. Large declared setup caps
can therefore require a larger `memory-bytes` declaration even for a small grid.
A requested hierarchy that cannot resolve its coarse geometry, support, transfer
or factorization refuses instead of silently switching algorithms. A missing
coarse-space support cannot be repaired by substituting a box clamp.

Compliance resume retains its verified-stage-replay semantics: it reconstructs
the same solver policy and charges repeated setup and geometry again. Stress
resume reconstructs the hierarchy and verifies retained endpoints without
replaying prior optimization steps. Both policies bind the original canonical
source and executable, and preserve physical solid/void regions. Interrupted
proposals do not replace accepted fields or geometry.

Omitting `linear-solver` preserves the old native arithmetic: Jacobi optimization
and two-level goal enrichment. No automatic default switch is made. The
multilevel option is an available numerical method, not a guarantee of fewer
iterations, lower wall time, mesh-independent convergence, or continuum safety
on every geometry/material distribution.

## Focused regression commands

```sh
cargo test --release -p fs-topopt --features cutfem-marquee --lib sdf3::adaptive_continuation::multilevel_tests
cargo test --release -p fs-cli --features sdf3-study --lib study::elasticity::sdf3::solver
cargo test --release -p fs-cli --features sdf3-study --test study_multilevel_cli
```

The tests exercise independent equilibria and gradients, preparation reuse,
actual adaptive updates, rejected work, cancellation, and both recovery paths.
Test source is not an execution receipt.
