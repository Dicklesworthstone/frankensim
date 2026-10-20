# Binary endpoints under every eroded compliance cap

```sh
cargo run --release -p fs-topopt --features cutfem-marquee \
  --example worst_load_sdf3 -- --binary-material 1000 1500 200 250000
```

This mode runs the existing fixed-grid cap-constrained minimum-material study,
then checks five explicit thresholds `[0.2, 0.35, 0.5, 0.65, 0.8]`. The first two
arguments are authored compliance requirements, not baseline-inferred ratings.
The remaining arguments bound gray optimizer updates and cumulative Krylov
iterations. Endpoint filtering, preconditioner preparation and all fresh load
solves share that same physical allowance; it is not reset after optimization.

## What changes when gray material becomes solid/void

One raw design is filtered once and projected at the declared eroded, nominal
and dilated thresholds. Each resulting physical field is made binary at the
SAME endpoint threshold, with all exact ties included. Protected solid/void
regions remain fixed in every realization. The masks must be pointwise ordered;
an ordering defect is refused, not repaired by changing the material.

For each candidate, actual eroded solid cells use unit reference stiffness;
eroded voids retain an explicitly declared positive stiffness (1e-3 in this
example). EVERY independent load is freshly solved on that material with the
existing preconditioner and true-residual gate. Gray displacements, gray
feasibility, and gray SIMP gradients are not reused. No derivative is claimed
for the discontinuous threshold operation.

`evaluate_binary_compliance_caps` returns the three binary masks and numerical
cut volumes, actual eroded stiffness multipliers, per-load compliances, the
smooth aggregate and fresh displacement fields. `feasible` requires every
physical compliance to satisfy its own cap directly; numerical tolerance belongs
only to the smooth aggregate. Both returned errors and successful evaluations
restore the source's incoming nominal scales and leave its parameters unchanged.

## Bounded selection

`scan_binary_compliance_caps` admits the entire threshold schedule before work,
then chooses the feasible candidate with least DILATED binary material. It also
checks a separate explicit dilated-volume limit. Equal-volume ties retain the
first authored schedule entry. It does not assume feasibility is monotone,
interpolate between candidates, or split ties to meet an exact volume.

An interrupted scan preserves completed summaries and a completed feasible
incumbent, if one exists. The example exports that binary prefix before returning
nonzero. A complete scan with no feasible candidate returns `best = None` and
the example exits nonzero; previously printed gray fields are not a successful
binary fallback. The minimum is over the declared candidates only, not all
thresholds or material layouts.

## Scope

This is a fixed-grid ersatz MATERIAL endpoint, not a boundary mesh, CAD or printer
file. Loads, supports, cut geometry and quadrature remain present in void cells;
there is no hidden load deletion or zero-void claim. Feasibility can depend on
the positive void stiffness. The eroded family is solved; nominal/dilated masks
are not silently advertised as separately verified load solves. There is no
continuum guarantee, connectivity proof, manufactured feature-size guarantee,
stress/buckling qualification, or global optimality claim.

The existing weighted-sum, worst-load, minimum-material and refinement modes
retain their original routes and numerical algorithms.

## Native checks

```sh
cargo test --release -p fs-topopt --features cutfem-marquee \
  --lib sdf3::robust::envelope::limits::binary
cargo test --release -p fs-topopt --features cutfem-marquee \
  --example worst_load_sdf3 binary::tests
```
