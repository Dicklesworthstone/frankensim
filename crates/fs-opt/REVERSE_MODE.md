# Reverse-mode objectives and constraints

`fs_opt::reverse::ReverseProgram` compiles chosen scalar roots from the live,
sealed Problem IR. The solver-facing `fs_opt::ReverseProblem` adapter compiles
all objectives and constraints together, then applies the declared objective
weights and senses. Repeated gradient and Jacobian-transpose products reuse one
valuation; there is no coordinate perturbation or dense Jacobian construction.

```rust,ignore
use fs_opt::ReverseProblem;
use fs_opt::reverse::ReverseLimits;

let limits = ReverseLimits {
    max_nodes: 1 << 18,
    max_scalar_slots: 1 << 22,
};
// problem is a sealed fs_opt::Problem. Bind every declared variable, in
// declaration order, using its manifold's point-storage coordinates.
let oracle = ReverseProblem::new(&problem, limits)?;
let point = oracle.evaluate(&bindings)?;
let objective = point.objective_value();
let gradient = point.objective_gradient()?; // one vector per variable
let residuals = point.constraint_values();
let jt_lambda = point.constraint_pullback(&multipliers)?;
let lagrangian_gradient = point.lagrangian_gradient(&multipliers)?;
let parameter_gradient = point.objective_parameter_gradient()?;
```

## Meaning of the results

The minimized scalar objective is the sum of `weight * f` for minimization
entries and `-weight * f` for maximization entries, accumulated in declaration
order. `objective_values()` retains the raw individual values. A non-finite
weighted sum is refused with the objective-list index, rather than being
reported as a finite objective or a successful optimization result.

`constraint_values()` retains every original residual in declaration order,
including negative values of satisfied inequalities. Inspect the corresponding
`oracle.problem().constraints()` entries for `EqZero` or `LeZero` kinds. Neither
feasibility clipping nor penalty construction occurs in this adapter. Each
constraint pullback requires exactly one finite multiplier per constraint.
Multiplier sign feasibility belongs to the constrained solver.

The objective, constraint and Lagrangian gradient methods return either ambient
point-coordinate vectors or authoritative manifold-parameter vectors. SO(3)
uses four quaternion point components and three body-frame parameters. Sphere
and Stiefel pullbacks delegate to the existing manifold implementation; this
adapter introduces no independent retraction or projection formula.

The valuation owns its captured primal point through the underlying reverse
tape. Later modifications of the caller's bindings do not change retained
values or gradients. The immutable Problem and compiled program remain borrowed,
so a tape cannot be rebound to a different graph.

## Direct root-level use

The existing lower-level API remains unchanged:

```rust,ignore
use fs_opt::reverse::{ReverseLimits, ReverseProgram};

let limits = ReverseLimits { max_nodes: 1 << 18, max_scalar_slots: 1 << 22 };
let program = ReverseProgram::new(&problem, &[objective, constraint], limits)?;
let tape = program.evaluate(&bindings, None)?; // explicit cancellation opt-out
let root_values = tape.values();
let gradient = tape.pullback(&[1.0, 0.0], None)?;
let parameters = tape.parameter_pullback(&[1.0, multiplier], None)?;
```

The limits above are example values, not defaults. `max_nodes` bounds graph,
root and variable counts separately. `max_scalar_slots` bounds the underlying
program's scalar storage, including every declared variable's point. Each
valuation stores its primals; each pullback additionally allocates adjoints and
returned per-variable gradients. The adapter also retains root weights and
uses a combined seed vector for constraint/Lagrangian products. Allocations
are fallible; these caps are not byte-accurate memory or time certificates.

## Solver example and checks

The example supplies objective and constraint derivatives to the existing SQP
engine for a quadratic objective with affine equality and inequality constraints:

```sh
cargo run -p fs-ascent --example reverse_ir_optimize
cargo test -p fs-opt --test reverse_problem
cargo test -p fs-ascent --example reverse_ir_optimize
```

It checks the analytic optimum `(1.2, 0.8)` and the returned KKT residual. Its
`expect` calls bridge a closed, globally smooth fixture to the current
infallible solver callbacks; they are not a general panic-free solver adapter
for domain-limited expressions. Its problem evaluation budget is explicitly
unlimited and its SQP iteration count is bounded.

## Boundaries

This is a valuation/derivative API, not a solve loop: callers own and enforce
objective-evaluation budgets. Existing Study finite-difference trajectories,
checkpoint semantics, and default solver routing are unchanged.

`ReverseProblem` refuses missing objectives and chance, bilevel or multi-fidelity
tags instead of silently ignoring their semantics. The underlying compiler
rejects reachable kinks and unimplemented physics/UQ execution. It does not
manufacture PDE adjoints from availability metadata. Non-finite primal or
adjoint arithmetic remains a typed refusal.

Every problem-level evaluation and derivative operation has a `_cancellable`
variant. Wrapper root/seed/scalarization loops poll in bounded chunks and pass
the context to the existing numerical sweeps. The underlying compiler is
bracketed by cancellation checks as a whole phase, not interrupted internally.
Binding validation and public manifold operations retain the same whole-phase
boundaries as the underlying API. Errors and cancellation return no partial
valuation or derivative and do not mutate the compiled problem or bindings.

Gradients are chain-rule derivatives of the mathematical operators, not of
floating-point rounding. No interval certificate, cross-ISA bit identity,
Hessian product, generic nonlinear convergence guarantee, or PDE adjoint is
claimed by this adapter.
