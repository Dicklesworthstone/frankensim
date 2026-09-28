# Mesh-checked native stress constraints

The single-load `projected-stress` study now accepts an optional `:mesh-check`.
The original sampled von Mises limit and its explicit absolute allowance apply
on the optimization grid AND every declared finer grid. A coarse-grid-safe
candidate cannot be accepted when a finer-grid sample exceeds that same limit.

```sh
cargo run -p fs-cli --bin frankensim -- --json study \
  examples/marquee/bracket-mesh-stress-2d.fsim /tmp/stress-mesh.db
```

This is a numerical software example with explicit synthetic material/load
scales and a 20 Pa sample limit, not a calibrated bracket or an engineering
allowable. Neither baseline admission, mesh agreement nor descent is promised.
Inspect the actual status and measurements; do not raise a limit to conceal a
refinement-sensitive stress concentration.

Place this block after `:stress-tolerance-pa` and before any `:design-regions`:

```lisp
    :mesh-check (refinement
      :extra-levels 1
      :absolute-compliance-tolerance-j 0.000001
      :relative-compliance-tolerance 0.05
      :area-tolerance-m2 0.005)
```

All four fields are required. They retain the existing volume-mode mesh-check
meaning: one or two extra levels through level seven, adjacent-grid compliance
agreement and comparable material area. The stress limit is NOT inferred from
those tolerances. It remains the existing `:sampled-stress-limit-pa` plus
`:stress-tolerance-pa`, unchanged at every level. No stress-convergence rate or
continuous maximum is inferred from passing the sampled limits.

The accepted baseline and each otherwise admissible candidate are evaluated by
the original CutFEM stress evaluator. One solve at each level provides actual
compliance, material area, sampled stress, maximum location and probe count.
The same bilinear geometry is prolonged; no fine-grid area projection,
redistancing or interpolation of displacements alters the comparison. A
candidate must improve compliance on every matching grid and meet every grid's
stress limit before the existing projected owner changes its accepted state.

An unresolved or finer-grid-overstressed baseline returns `mesh-unresolved`
with the budget exit code and a durable unchanged endpoint. A rejected
candidate contracts the existing search; exhausted attempts return
`no-feasible-descent`, never convergence. The requested work uses the existing
lifetime wall budget and CG/stress-cell cancellation callbacks. Cancellation
discards unfinished checks, not prior accepted geometry or complete checks.
Assembly, prolongation, quadrature and ledger calls remain indivisible.

The ordinary `study --resume`, `report` and `package` interfaces retain the
workflow. `constraints.mesh_resolution` contains the policy and last complete
check, including baseline rungs and the accepted candidate's rungs. Earlier
checks remain in predecessor receipts. JSON retains each maximum's location,
probe count, field snapshot, area and compliance; HTML shows a per-grid stress
table. Terminal mesh-checked receipts return unchanged on resume, without new
physics. Report readers recheck the original gates rather than accepting an
unsupported headline pass. The package remains structural Estimated evidence.

Omitting `:mesh-check` preserves existing canonical declarations and numerical
behavior. Combining it with independent-load or stress-restoration mode is
explicitly refused for now; no load family is silently replaced by one load.
Volume-mode mesh checks and explicit coarse-to-fine volume studies remain
separate and unchanged. This increment does not provide automatic adaptivity,
a stress adjoint, continuous stress certification or physical validation.

Focused checks use the existing test lanes:

```sh
cargo test --release -p fs-topols --lib projected_stress::resolution
cargo test --release -p fs-cli --lib study::elasticity
cargo test --release -p fs-cli --test study_checkpoint_cli
```

Regression source is not a passing-test claim. Consult the actual native run
before treating a particular build as verified.
