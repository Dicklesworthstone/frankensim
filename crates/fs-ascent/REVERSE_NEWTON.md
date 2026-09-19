# Matrix-free reverse Newton studies

`ReverseNewtonStudy` connects the live algebraic IR to the existing
Newton--Krylov trust-region engine. Users provide neither an objective callback
nor a Hessian callback. The accepted primal tape is reused by every Steihaug
Hessian-vector product, including the final quadratic-model check. There are no
coordinate perturbations and no dense Hessian or BFGS approximation.

```rust,ignore
use fs_ascent::{ReverseNewtonStudy, StopRule};
use fs_opt::{ReverseProblem, reverse::ReverseLimits};
let oracle = ReverseProblem::new(&problem, ReverseLimits {
    max_nodes: 100_000,
    max_scalar_slots: 1_000_000,
})?;
let mut study = ReverseNewtonStudy::new(&oracle, &initial_point, Some(&cx))?;
let report = study.run(&StopRule::GradNorm(1e-8), 100, 2_000, Some(&cx))?;
let checkpoint = study.clone();
let accepted = study.snapshot();
```

The final numeric run argument is a **cumulative Hessian-product ceiling**,
separate from the problem's cumulative objective-evaluation limit. Every
`StopRule::Budget` leaf further limits objective attempts, even inside `All`.
Initialization counts as one objective/gradient attempt. Products count when
attempted, including failed or cancelled products; no new primal is evaluated
for them. An authorized final call may complete and publish its iteration.
An exhausted objective budget takes priority over simultaneous convergence.
`HessianBudget` means another product would be required but is unfunded; it is
not a numerical stall or a convergence claim.

Each complete existing trust iteration is staged before publication. A
non-finite trial objective or gradient contracts the radius and records the
original refusal in `last_rejection()`. A failed Hessian, cancellation, or a
mid-Krylov Hessian budget stop discards the staged iteration and leaves the
accepted point, tape, gradient, radius and history unchanged. Spent objective
and Hessian work is not refunded. Raising a product allowance permits retry;
unfinished Krylov work is restarted, so such a resume can repeat products.
Splits at complete outer-iteration boundaries retain trajectory and counts.
Snapshots expose actual attempted counts, not the temporary kernel counters.

The oracle is borrowed immutably and the accepted primal is shared using Arc.
A checkpoint cannot silently resume with different variables, weights or
constraints. Budgets belong to each checkpoint; clones do not install a shared
cross-branch resource pool. AD sweeps poll Cx internally, but Packing, state
cloning and existing vector kernels remain whole phases. There is no hard
wall-clock cancellation bound or process-restart serialization. Staging copies
solver state and history; this is not an allocation-free implementation.

## Derivative API and scope

`ReverseEvaluation::hessian_vector_product` applies a weighted-root Hessian.
`ProblemEvaluation::objective_hessian_vector_product` applies the signed,
weighted objective Hessian. `lagrangian_hessian_vector_product` adds constraint
Hessians in declaration order with fixed multipliers. These reuse stored
primals with a directional forward and differentiated reverse sweep. Working
storage is three scalar-slot arrays, a node mask and the output, not n squared.
The reverse compiler's scalar/node caps are not a total solver-memory budget.

Products use **ambient point coordinates**. They do not differentiate the
manifold pullback and are not intrinsic Riemannian Hessians. The Newton study
therefore admits only unconstrained Euclidean variables. Use the existing
constrained study for constraints; no declared constraint is dropped. Kinks,
external PDE/UQ nodes and unsupported tags retain existing refusals. Active
primitive derivatives must remain finite. These are chain-rule derivatives
of the mathematical operators, not derivatives of floating-point rounding,
interval certificates, global-optimality proofs or a PDE-adjoint implementation.
Zero gradient alone does not prove positive curvature or escape an exact saddle.
Existing legacy trust, first-order and constrained paths remain unchanged.

## Runnable nonconvex example and checks

```sh
cargo run -p fs-ascent --example reverse_newton
cargo test -p fs-opt --test hessian_products --test reverse_problem
cargo test -p fs-ascent --test reverse_newton --example reverse_newton
cargo test -p fs-ascent --lib trust::tests
```

The 64-coordinate example uses nonnegative quartic double-well potentials and
nearest-neighbor quadratic coupling. Starting at 0.25 exercises genuine negative
curvature, then approaches the known uniform positive minimum. It is explicit
dimensionless algebra, not an identified material or mesh-based physics model.
Native regression targets cover analytic products, shared vectors, constraints,
symmetry, directions, second-order failures, domain recovery, strict dual work
budgets, checkpoint splits and real-Cx cancellation. They were added but not
run in the authoring environment because Cargo/rustc are unavailable. Executed
Python mathematical translations are not evidence of a successful Rust build.
