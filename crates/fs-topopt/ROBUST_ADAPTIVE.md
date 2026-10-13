# Robust adaptive design: assess, selectively rebuild, reoptimize

```sh
cargo run --release -p fs-topopt --features cutfem-marquee \
  --example elastic_sdf3 -- --robust-adaptive 200 250000
```

The example runs the robust `(penal,beta)=(1,1),(2,2),(3,8)` schedule, then
performs up to two selective refinement/refit cycles at the final material
model. Each model/grid gets up to 200 accepted primal steps. The entire call
shares 250,000 Krylov iterations; geometry builds share one separate cumulative
quadrature allowance, including failed builds. The fixed objective conditioning
scale is explicitly 1000 in this dimensionless example. All policy values and
spent work are printed. The original `--robust-continuation` mode remains a
continuation plus diagnostic probe, with no fine-grid reoptimization.

## Reusable APIs

`CutDensityStudy3::refit_robust_reference` takes an explicitly constructed finer
operator, source raw controls, the new material parameters, reference load
laws and `RobustRefitOptions3`. The checked inter-grid map must enrich the
active space without losing source support. Raw controls and protected physical
solid/void labels inherit through that same map. The source filter radius is
preserved. Small declared raw floors are respected rather than raised to 1e-3.

The candidate rebuilds its filter and cut-volume measures, reintegrates every
independent reference load, restores its **dilated** material constraint,
checks eroded-compliance and dilated-volume gradients and runs a fresh robust
optimizer. No nodal force vector, solved displacement, old gradient or optimizer
multiplier is interpolated. The source is immutable; the returned candidate is
not automatically installed. Its report can retain an interrupted prefix.

`robust::adaptive::controlled_adaptive_robust_design` closes the repeated loop.
The initial grid receives the supplied material schedule; later grids refit the
final admitted model. Each eroded-compliance probe uses globally enriched
geometry with inherited **physical stiffness**. Its Dörfler marks are applied
to the original source tree, not to the already enriched probe tree. The new
design problem then inherits **raw controls** and refilters them. These are
separate numerical experiments: their compliance values do not establish
cross-grid optimization descent.

The caller's builder must preserve SDF, domain, base material, supports and
integration policy, enforce a cumulative geometry budget and honor its
checkpoint. Tree level/leaf limits and sparse transfer bounds remain enforced.
The driver installs a completed feasible candidate, its tree and its report
only after a publication checkpoint. Failures before it preserve the previous
grid; interruption after it preserves the complete new grid. Reports separate
installed grids, diagnostic probes and uninstalled candidate results. Nominal
operator scales are working state; use the retained point's eroded fields.

## Stop meanings

A refinement-count limit is **not** an accuracy or optimizer-convergence pass.
An empty marking signal also does not establish accuracy. Optional
`eroded_relative_tolerance` compares each positive-weight load's independently
solved coarse/fine compliance; signed changes across different loads cannot
cancel into a pass. With that gate enabled, the current terminal grid is probed
even when the refinement budget is exhausted. Without it, reaching the requested
refinement count stops without an unnecessary final diagnostic solve.

A compliance-tolerance stop is only a two-grid numerical observation. It does
not certify continuum error, dilated-volume discretization error, geometric
offsets, minimum feature size, safety or a global optimum.

## Native validation commands

```sh
cargo test --release -p fs-topopt --features cutfem-marquee --lib sdf3::robust::refit::tests
cargo test --release -p fs-topopt --features cutfem-marquee --lib sdf3::robust::adaptive::tests
cargo check --release -p fs-topopt --features cutfem-marquee --example elastic_sdf3
```

These commands have not run successfully in the editing environment, which
lacks Cargo/rustc. Independent finite-element analogues validate the transfer,
filter/projection gradients and constrained-refit formulation, not execution
of these Rust implementations or production CutFEM tests.
