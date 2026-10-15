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
Projection/penalization continuation must rebuild immutable physics bindings;
there is no automatic continuation schedule or budget reset in this API.

## Run and verify

```bash
cargo run --release -p fs-topopt --features newton-design --example robust_simp
cargo test --release -p fs-topopt --features newton-design \
  --test robust_simp_newton --test simp_newton --test simp_hessian
```

The example has 48 independent design cells, unequal tetrahedral cell volumes,
two independent loads, and a dilated-volume allowance of 0.6. It reassembles the
candidate through the older first-order robust path, checks agreement with the
Newton report, and prints actual termination, chart-coordinate KKT residuals,
all three volumes, and provider work. Its separate verification solves are not
included in the optimization provider's work counter.

The new native integration tests cover both objective/constraint orientations
against an analytically known uniform-density optimum, fresh physical checks,
split replay, shared work and failed-cache budget recovery. Native compilation,
these tests and the example were **not executed in the authoring runtime**:
`dsr quality --tool frankensim` returned exit 127, and the alternate DSR path,
RCH, Cargo and rustc were unavailable. Repository DSR policy remains authoritative.
Independent dense reference equations passed 72 derivative cases; those are not
substitutes for compiling or executing the Rust implementation.

## Scope

This is fixed-mesh, fixed-load, homogeneous-support linear-elastic SIMP design.
Eroded stiffness must stay in the native Hessian's smooth projected interval
`(1e-12,1)`. A logistic open box does not waive this restriction or provide
closed-box boundary KKT. Numerical three-field order is not a certified minimum
feature size, binary topology, manufacturing yield, global optimum, or continuum
certificate. Stress/shape derivatives, global byte admission, and serialized
process restart remain outside this implementation.
