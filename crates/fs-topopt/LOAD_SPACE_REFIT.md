# Refit simultaneous-load requirements on finer design grids

```sh
cargo run --release -p fs-topopt --features cutfem-marquee \
  --example worst_load_sdf3 -- --load-space --refine-material 1500 200 250000 0.02
```

Arguments are an authored physical compliance cap, accepted optimizer updates
per grid (default 200), one cumulative Krylov-iteration allowance (250000), and
a positive full normalized-Gram change tolerance (0.02). The numbers shown are
illustrative requirements in dimensionless units, not calibrated ratings or a
promise that this fixture has a feasible result. The load set remains the
original Y and correlated Y/Z columns: `f = Fq`, `||q||_2 <= 1`.

The command restores and optimizes an initial adaptive grid, freshly rechecks
its complete load set, and assesses that physical material on an enriched grid.
If both caps and the full-Gram change check pass, it finishes without refitting.
Otherwise it performs ONE refined-grid minimum-material solve and assesses the
NEW design on another enriched grid. The cap, load laws, material parameters,
filter radius and spectral temperature remain fixed throughout. Failure of the
last check returns nonzero; completing an optimizer allowance is not a pass.

This example refines every active background cell. The basis-wise DWR marking
masses are not treated as error localization for the entire load ellipsoid.
The reusable refit API also accepts checked selectively refined operators.
At most four geometry builds occur: initial, initial probe, refit, final probe.
They share one quadrature allowance, including failed builds. All filtering,
restoration, gradient audits, optimization and rechecks share one SolveControl.

## One numerical path for initial and refined designs

`solve_reference_load_space_cap` consumes a working adaptive density study,
integrates every reference law and runs existing full-spectrum restoration,
gradient checking, constrained optimization and a complete-basis recheck.

`refit_load_space_cap` borrows the source immutably and accepts a genuinely
finer operator, raw controls, explicit nominal parameters, the SAME reference
laws and `LoadSpaceCap3`, `LoadSpaceCapRefitOptions3`, and shared solve control.
It transfers raw controls and protected solid/void labels through the checked
parent map, preserves filter radius, and calls that same initial-grid path.
Fine filtering and cut-volume measures are rebuilt. Small declared raw floors
are not silently raised. No force vector, displacement, old worst direction,
gradient or optimizer state is transferred.

This refit differs from the two-grid diagnostic: the latter inherits physical
stiffness; the new optimization inherits and refilters raw controls. Their
compliances do not demonstrate cross-grid optimization descent. The caller
must preserve SDF, domain, material, supports and reference-load identity, and
own the geometry budget. Transfer checks compatibility, not geometric identity.

`LoadSpaceCapRefit3` owns the uninstalled candidate and separately retains its
restoration seed, accepted history, least-volume feasible incumbent, optimizer
outcome and fresh full-basis recheck. Resource limits, initialization failures
and late cancellation cannot silently replace the source or discard an admitted
prefix. `completed_and_rechecked()` requires a completed bounded optimizer run
and a fresh physically feasible check. It does not mean stationarity or accuracy.
The option/error types reuse the existing finite-case refit policy units; the
constraint itself always uses the full simultaneous-load spectrum.

The command exports grid-labeled histories, all three density fields, worst
load coefficients and actual force/displacement fields before later work can
fail. Restoration seeds are labeled as seeds, not optimized results. Each probe
prints its own pair of spectra, physical compliances, worst coefficients and
full normalized-Gram change. No passing individual basis case can authorize a
combined-load pass, and the final probe assesses the refitted material rather
than recycling the preceding coarse design's diagnostic.

This is numerical, linear-elastic, fixed-SDF evidence. It does not certify
continuum error, dilated-volume discretization error, stress/buckling safety,
connectivity, minimum feature sizes, binary-manufacturing behavior or global
optimality. Existing minimum-material and binary command routes are unchanged.

## Native checks

```sh
cargo test --release -p fs-topopt --features cutfem-marquee \
  --lib sdf3::robust::envelope::load_space::design::limits::refit
cargo test --release -p fs-topopt --features cutfem-marquee \
  --example worst_load_sdf3 load_space::refinement::tests
```

These new tests have not been executed in the editing environment because DSR,
its documented fallback, RCH and Cargo are absent. Independent body-fitted FEM
reference checks are not execution of the Rust/CutFEM implementation.
