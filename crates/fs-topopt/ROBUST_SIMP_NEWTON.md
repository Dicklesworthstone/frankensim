# Robust eroded/dilated design with native Newton derivatives

The existing `RobustPipeline` can now be used by the shared live-IR constrained
Newton path. `SimpComplianceGoal::new_robust` provides the low-level native
linearization; `SimpNewtonModel::new_robust` exposes named compliance and volume
views for `PhysicsBinding`, `LogisticBox`, and `AugmentedNewtonStudy`.

```rust,ignore
let physics = SimpNewtonModel::new_robust(
    &robust_pipeline, &elasticity, &independent_loads,
    &physical_cell_volumes, solve_options, lifetime_budget,
)?;
let raw_compliance = physics.compliance("eroded-compliance")?;
let raw_volume = physics.volume("dilated-material")?;
let compliance = coordinates.bind(&raw_compliance)?;
let volume = coordinates.bind(&raw_volume)?;
```

## Which physical design is being optimized?

For the nominal threshold `eta` and offset `delta = eta_offset`, compliance uses
`eta + delta` (erosion), while the volume view uses `eta - delta` (dilation).
Thresholds must remain in `[0,1]`; the offset must be finite and nonnegative.
`projection_thresholds()` reports the resolved pair. Zero offset is ordinary
SIMP. The older `SimpNewtonModel::new` keeps its original equal-threshold meaning.

The **volume constraint is dilated physical material**. It is neither nominal
volume nor SIMP-penalized stiffness. Unlike the older robust OC driver's adaptive
nominal-volume target, this Newton path has no implicit target adaptation: the
caller writes the dilated allowance explicitly in the IR. The example reports
all three volumes so this distinction remains visible.

One forward Helmholtz solve produces the common filtered density. Independent
loads are solved on eroded stiffness, and the compliance pullback includes the
full second derivative of the eroded projection/SIMP map. Volume uses its own
dilated projection slope and curvature:

```text
H_compliance v = F^T [m_e' H_modulus(m_e' Fv) + m_e'' Fv g_modulus]
H_volume v     = F^T [physical_volume_weights * h_d'' * Fv]
```

The existing logistic coordinate adapter adds its own full coordinate Hessian;
no finite-difference production derivative, dense design Hessian, or alternative
Newton/linear solver is introduced.

## Work, refusal and continuation

The two scalar views share one exact-point cache and one work owner. Fresh
preparation reserves `number_of_loads + 3` solves, exactly as for joint ordinary
SIMP: one forward filter, the independent primal loads, and two transpose
filters. Reading the other scalar at the same point does not repeat physics.
A compliance product reserves `positive_weight_loads + 2` solves; a nonlinear
volume product reserves two filters. Reservations are conservative, while
started solves and actual Krylov iterations are reported separately.

A failed new-point preparation preserves the preceding cache and charged work.
Older accepted tapes reconstruct their own point after a successful trial evicts
it. Increasing `set_budget` permits recovery without refunding attempts. Use a
real `Cx` for cancellation; `None` explicitly disables cooperative cancellation.

### Automatic immutable-stage continuation

`fs_ascent::newton_continuation::NewtonContinuation` runs an explicit schedule
of `NewtonStage` values, each borrowing its own `ReverseProblem`. It advances
only when the current stage meets its declared numerical KKT tolerance. Only
accepted coordinates cross the boundary: primal tapes, gradients, multipliers,
merit history and trust radius are initialized against the new model. A failed
or cancelled initialization leaves the preceding accepted stage available.

Create ONE `SimpNewtonResources` and clone its handle into every
`SimpNewtonModel::new_robust_with_resources` (or `new_with_resources` for ordinary
SIMP). The models keep separate exact-point caches but share all lifetime solve
reservations and actual Krylov work. Creating a stage does not reset the physical
allowance. `model.work()`, `model.budget()` and `set_budget()` then refer to that
shared account; original constructors still create independent accounts.

The continuation driver separately enforces one cumulative objective/Hessian/
multiplier budget. Every attempted stage initialization is conservatively charged
as one objective attempt, even when admission or a physical callback fails.
Zero-iteration calls do no solves or handoffs. Clone retains current progress;
provider caches and work accounts stay outside, and are shared rather than reset.
Use `active()`, `work()` and `completed_stages()` after a refusal, raise the
relevant absolute ceiling, and run again. No stalled stage is silently skipped.

Stage variables must have identical names, dimensions and Euclidean layouts.
Reuse the SAME coordinate chart across stages: the driver cannot infer a hidden
change of coordinate meaning inside a callback. Objectives belong to different
models; their stage reports are not one monotone descent history. The last
stage's KKT convergence is not a global or closed-box boundary certificate.

## Run and verify

```bash
cargo run --release -p fs-topopt --features newton-design --example robust_simp
cargo run --release -p fs-topopt --features newton-design --example robust_simp -- --continuation
cargo test --release -p fs-ascent --test newton_continuation
cargo test --release -p fs-topopt --features newton-design \
  --test robust_simp_newton --test simp_newton --test simp_hessian
```

Both example modes have 48 independent design cells, unequal tetrahedral cell
volumes, two independent loads, and a dilated-volume allowance of 0.6. The
unchanged default mode uses one fixed model. The continuation mode schedules
`(penal, beta, offset)` as `(1,0,0)`, `(2,1,0.05)`, `(3,2,0.1)` with KKT tolerances
`1e-4`, `3e-5`, `1e-5`; its fixed objective normalization is 100 joules. It reports
the actual reached stage and stopping reason, rather than calling every run
complete. Each published stage is checked by reassembly and separate solves
through the older first-order robust path. These verification solves are not
charged to the optimization provider's work counter.

The six continuation regressions cover exact split runs, combined budgets,
failed/cancelled stage entry, constrained KKT solutions and schema admission.
Native compilation, these tests and the examples were **not executed in the
authoring runtime**: `dsr quality --tool frankensim` returned exit 127; the
alternate DSR path, RCH, Cargo and rustc were unavailable. Repository DSR policy
remains authoritative. Earlier dense-reference derivative checks do not verify
this newly added orchestration or substitute for native execution.

## Scope

This is fixed-mesh, fixed-load, homogeneous-support linear-elastic SIMP design.
Eroded stiffness must stay in the native Hessian's smooth projected interval
`(1e-12,1)`. A logistic open box does not waive this restriction or provide
closed-box boundary KKT. Numerical three-field order is not a certified minimum
feature size, binary topology, manufacturing yield, global optimum, or continuum
certificate. Stress/shape derivatives, global byte admission, and serialized
process restart remain outside this implementation.
