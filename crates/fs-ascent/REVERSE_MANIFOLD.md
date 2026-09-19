# Optimize mixed manifold variables directly from the IR

`ReverseManifoldStudy` solves smooth unconstrained problems whose variables may
mix `Rn`, `Sphere`, `So3` and `Stiefel`. It consumes the existing `ReverseProblem`
without user callbacks or coordinate finite differences. Every sample uses one
shared primal evaluation and one weighted-objective reverse sweep. Objective
senses and weights remain exactly those declared in the problem.

```rust,ignore
use fs_ascent::{ReverseManifoldStudy, StopRule};
use fs_opt::{ReverseProblem, reverse::ReverseLimits};
let oracle = ReverseProblem::new(&problem, ReverseLimits {
    max_nodes: 100_000, max_scalar_slots: 1_000_000,
})?;
let mut study = ReverseManifoldStudy::new(&oracle, &initial_point, 7, Some(&cx))?;
let report = study.run(&StopRule::GradNorm(1e-7), 100, Some(&cx))?;
let checkpoint = study.clone();
let point = study.point();
let parameter_gradient = study.gradient();
```

## Geometry and the actual solve

The ordered `ProductManifold` uses variable declaration indices as factor IDs.
Its new `parameter_gradient`, `validate_parameter_tangent`, `retract_curve` and
`transport_parameter` operations delegate to the existing factor authority.
Point storage and optimization parameters are never interchangeable: SO(3)
stores four quaternion components but contributes three body-frame parameters.
Only quaternion starts adopt their canonical antipodal representative; other
valid initial points keep their original bits. SO(3) objectives should respect
`q ~ -q`; this driver cannot prove that an arbitrary graph does so.

L-BFGS directions are projected where appropriate, then the existing fallible
strong-Wolfe search evaluates the real retraction curve and its derivative.
Old curvature pairs are transported to the accepted point and readmitted using
positive, scale-relative curvature. The update uses the same cautious stretch
scaling as the single-manifold Riemannian engine. A numerical transport refusal
restarts optional memory and is exposed through `memory_restarts` and
`last_rejection()`. It does not corrupt an otherwise valid Wolfe step. Memory
zero is retracted steepest descent. No single-factor geometry is reimplemented.

The product uses the existing factor parameter metrics. Reported norms are in
these concatenated coordinates, not physical units or an automatic physical
preconditioner. Nondimensionalize mixed physical variables deliberately.

## Budgets, failures and continuation

The sealed problem's evaluation limit and every `StopRule::Budget` leaf are
hard cumulative sample-attempt ceilings, including beneath `All`. Initialization
costs one attempt. Rejected geometry, numerical trials and failed in-flight
samples also count; each attempts at most one primal/reverse pair. Budgets are
checked before every Wolfe probe including zoom. Budget has stop-attribution
priority and a final funded trial may still be accepted. Zero-step continuation
only inspects cached state. Clone retains the point, gradient, curvature,
history and counts; the borrowed oracle prevents switching problem meanings.
Budgets are per checkpoint, not a shared pool across cloned branches.

Invalid initial values are errors. Trial domain/nonfinite errors instead cause
backtracking and retain their original diagnostic. Cancellation and other
errors do not publish an unfinished step or refund attempted work. Transport
is staged before accepted-state publication. A cancelled search restarts on
resume and may repeat probes; completed-step splits preserve trajectory and
accounting. Cx polls run inside reverse sweeps and bracket geometry/memory
operations; those operations, Packing and allocation remain whole phases.
There is no bounded wall-clock cancellation or allocation-free claim.

Additional EqZero/LeZero constraints are refused, not dropped: use the separate
Euclidean constrained driver for those. Unsupported tags, kinks and external
physics/UQ execution retain the reverse compiler's refusals. This is not a
Riemannian Newton method, a process-restart serialization format, a certificate
of global optimality, or proof that a stationary point is a minimum. Tape caps
are not a bound on total optimizer memory. Legacy solvers are unchanged.

```sh
cargo run -p fs-ascent --example reverse_manifold
cargo test -p fs-opt --test product_differential
cargo test -p fs-ascent --test reverse_manifold
```

The example couples a position to a direction and simultaneously optimizes a
rotation and two-column frame: 15 stored coordinates and 14 parameters. Its
nonnegative dimensionless objective has a known zero minimum. The native tests
cover this solve, all 13 completed-step splits, logarithmic-domain recovery,
budget/retry boundaries, real-Cx cancellation, quaternion antipodes, zero memory
and typed refusals. These Rust targets have not been executed in the authoring
environment (Cargo/rustc unavailable). Independent NumPy reference checks are
not evidence that the Rust targets build or pass.
