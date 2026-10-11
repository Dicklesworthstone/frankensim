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

## Preconditioning and inexact Newton steps

`run_preconditioned` uses a caller-supplied **fixed linear SPD inverse action**
inside each Krylov solve. The map may change with the accepted point between
outer iterations. It does not replace the true reverse-mode Hessian, change the
Euclidean trust ball or redefine the gradient-based outer stopping test.
Coordinates follow variable-declaration order; an inverse callback receives the
accepted packed point and an objective-normalized residual. A single fixed
binary normalization per inner solve controls extreme inverse-action units.
A positive residual product is checked locally; this does not certify global
SPD, symmetry or linearity of an arbitrary callback.

```rust,ignore
use fs_ascent::trust::KrylovOptions;
// Use the declared, positive stiffness diagonal of your model.
let mut inverse = |_: &[f64], residual: &[f64]| {
    residual.iter().zip(&stiffness_diagonal)
        .map(|(r, diagonal)| r / diagonal).collect()
};
// Optional inexact-subproblem policy. Defaults remain 1e-8 and at most 2*n
// Krylov products. This setting is retained when cloning the study.
study.set_krylov_options(KrylovOptions::new(1e-3, 20)?);
let report = study.run_preconditioned(&mut inverse,
    &StopRule::GradNorm(1e-8), 100, 2_000, 1_000, Some(&cx))?;
println!("inverse attempts: {}", report.preconditioner_evals);
```

The `20` above caps **Krylov products per outer iteration**, not total Hessian
work: a truncated step still needs one more Hessian product for its final model
check. It is an approximate Newton direction, and must pass the usual true
objective/gradient evaluation and model-agreement acceptance. This local cap
never borrows work from the cumulative `2_000` Hessian ceiling. Similarly, a
looser relative residual tolerance changes inner accuracy, not outer convergence.
Local caps work for both preconditioned and ordinary `run` calls.

The cumulative `1_000` limit counts inverse attempts, including malformed results,
contained callback panics, cancelled calls and discarded staged iterations.
`PreconditionerBudget` is a distinct resumable stop. Other inverse refusals return
`ReverseNewtonError::Preconditioner`; accepted state and tape remain intact.
`study.preconditioner_evals()` exposes attempted inverse work even after errors.
No inverse is evaluated after its budget, after residual convergence, merely to
prepare a direction beyond the configured local cap, or when the next Hessian
product is already unaffordable. Inverse callbacks are polled for cancellation
before and after invocation, but cannot be interrupted internally.

A solver checkpoint does not capture external mutable callback state. Replaying
it requires reproducing the same inverse map and inner policy. Changing the
policy spends no callbacks, refunds nothing and does not revive a numerical
stall. These controls do not make the solver allocation-free or install a total
memory/time budget.

The complete synthetic stiff-diagonal example is runnable with:

```sh
cargo run -p fs-ascent --example preconditioned_newton
```

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
budgets, checkpoint splits and real-Cx cancellation. Preconditioning regressions
also exercise exact-budget recovery after partial searches and inverse faults.
The focused production-source probe tests the numerical engine in unoptimized
and optimized builds; it does not replace native Cargo/Cx/IR integration or the
workspace DSR gate.
