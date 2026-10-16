# Adaptive minimum-material design under physical response limits

```sh
cargo run -p fs-topopt --features cutfem-marquee --release \
  --example minimum_material_sdf3 -- 80 250000 0.0012 0.003 \
  --adapt-rounds 3 --response-tol 0.001 --volume-tol 0.0001
```

This runs the existing tension/compression actuator example with repeated
selective refinement. `3` is the maximum total number of retained grids,
including the initial grid; the command accepts 2..=4. The 80-update allowance
applies per grid. One 250000-iteration Krylov budget and one bounded quadrature
control cover all fits, diagnostic probes and rejected candidates. The example
uses at most two marks per round, plus balancing, and reserves a final probe
level under the existing 4096-leaf cap. Legacy fit, `--estimate` and `--refit`
modes are unchanged and cannot be combined with `--adapt-rounds`.

A final acceptance requires BOTH the smooth constraint and every original
upper/lower limit to pass on a newly solved enriched physical model. It also
requires every normalized individual response change and the volume-fraction
change to meet their declared tolerances. Fine feasibility uses the same
numerical tolerance as the minimum-material optimizer's feasible archive;
loosening response-change tolerances cannot waive a violated force bound.
The example returns nonzero on any other outer stop, after exporting the last
retained design and, when available, all final source/probe residuals.

The enriched diagnostic inherits physical stiffness. Selective optimization
instead transfers raw density controls, rebuilds the unchanged filter/SIMP law,
inherits protected solid/void labels and solves a fresh baseline. Force windows
are not monotone in density: no volume-based scaling is used to invent feasible
responses. Only a successfully fitted candidate with an actual accepted feasible
incumbent may replace the source. Candidate failure preserves the source and
retains any completed candidate separately. Cross-grid volume decrease is not
asserted; source, baseline and selected volumes are reported separately.

Programmatic entry: `limits::refinement::adaptive::AdaptiveMinimumMaterialStudy3`
under `fs_topopt::sdf3::response`. Construct it with a matching background,
reference limit cases, raw design, explicit options and the original solve
control. `run(additional_refinements, factory, control)` accepts an elasticity
factory for `Probe` or `Refit` trees. The factory must use unchanged geometry,
material/support laws and the same bounded geometry control. Completed probes
are cached across pauses and invalidated on promotion. Even `run(0)` assesses
an unprobed source. Retry after a failed candidate spends additional work;
this is not an on-disk or mid-solve optimizer checkpoint.

`TwoGridLimitsSatisfied`, `RefinementLimit`, `Paused`, `NoRefinementSignal`,
`NoFeasibleSource` and `NoFeasibleCandidate` have different meanings. The
least-volume feasible incumbent is separate from the latest AL iterate and
its KKT report. Offset/measure changes are not fabricated into local indicators;
empty marks with a failed check are not an accuracy pass. These are numerical
two-grid tests, not continuum bounds, global optimality or manufacturing proof.
