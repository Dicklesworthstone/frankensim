# Minimum material with a simultaneous-load compliance requirement

```sh
cargo run --release -p fs-topopt --features cutfem-marquee \
  --example worst_load_sdf3 -- --load-space --minimum-material 1500 200 250000
```

The first value is an authored physical compliance cap, in this example's
explicit dimensionless units, NOT a baseline-derived rating. The remaining
values bound accepted optimizer updates and cumulative Krylov iterations.
The two load columns are the same Y and correlated Y/Z laws as `--load-space`.
The set is `f = F q`, `||q||_2 <= 1`, not a finite list or independently bounded
coefficient box. Passing each column separately does not establish this cap.

## Reused numerical path

`load_space::design::limits::LoadSpaceCap3` makes the requirement explicit.
Its fixed normalization equals the authored physical cap. The existing eroded
basis solves, cross-compliance Gram matrix and admitted numerical eigensolver
compute the entire load-set maximum. The minimum-material formulation uses
DILATED volume as objective and the smooth spectral maximum minus one as
constraint. The original volume-constrained constructor retains its original
objective, constraint, recovery and numerical loop; both formulations use that
same projected augmented-Lagrangian engine.

`LoadSpaceLimitedDesign3::best_feasible` selects the least-volume accepted
point meeting BOTH the numerical spectral maximum and the physical worst-force
work cap. Neither optimizer tolerance nor rounding of a normalized ratio may
waive a physical excess. Smoothing is conservative in exact arithmetic; it is
not a floating-point enclosure. Final KKT diagnostics belong to the final AL
iterate, not necessarily the retained feasible incumbent.

`limits::start::controlled_load_space_feasible_start` searches the segment
from incoming raw controls toward raw one. Every probe re-solves the full
load set; no worst direction is frozen. Caps, basis, material parameters and
protected physical regions remain fixed. The default is at most 16 bisections
after the incoming/full endpoint solves, with a 0.001 smooth reserve. A rejected
full endpoint is not a global infeasibility proof. Interruption retains any
fully checked feasible seed without calling it an optimized result.

The command prints the selected source, history, three density fields,
maximizing coefficients, physical cap comparison, worst force/displacement
and spectral diagnostics. All restoration, gradient, optimization and final
recheck work shares one physical budget. A fresh full-load-set recheck does
not modify AL history or counts. Solved fields are emitted before a late error.
Absence of a feasible seed or incumbent is never reported as success.

This is fixed-grid linear elasticity with homogeneous supports and an authored
zero-centered ellipsoid (at most 16 basis columns). It is not a probability or
risk model, rigorous continuum bound, stress/buckling qualification, geometric
length-scale guarantee, global optimum, or binary-material certification.

```sh
cargo test --release -p fs-topopt --features cutfem-marquee \
  --lib sdf3::robust::envelope::load_space
cargo test --release -p fs-topopt --features cutfem-marquee \
  --example worst_load_sdf3 load_space::material::tests
```
