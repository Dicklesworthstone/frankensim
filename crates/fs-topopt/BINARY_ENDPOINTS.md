# Binary material endpoints with fresh physical checks

```sh
cargo run -p fs-topopt --features cutfem-marquee --release \
  --example minimum_material_sdf3 -- 80 250000 0.0012 0.003 --binary
```

The existing minimum-material fit is followed by a bounded sweep of projected
physical-density thresholds: 0.2, 0.35, 0.5, 0.65, 0.8. Every endpoint is solved
against the same independent load, prescribed-motion and observation laws.
The sweep shares the fit's remaining Krylov/preparation work budget. No extra
geometry is built. Do not combine `--binary` with the adaptive/refit modes.

`evaluate_binary_response_limits` accepts one explicit threshold and positive
void-stiffness floor. Source controls pass through the original filter and
projection once. Thresholded material does **not** pass through them again:
that would reintroduce gray stiffness. A projected value equal to the threshold
is solid; ties are not split by cell index to fabricate an exact volume.
Protected solid and void physical regions retain their labels. The exact
stiffness multipliers are 1 for solid cells and the declared floor for void.

The current-density prescribed-motion lifting and numerical reaction are
recomputed. Displacements, force/moment observations and original signed
residuals are freshly evaluated through the same primal implementation as the
differentiable limit evaluator. No gray gradient, adjoint, displacement or
feasibility result is reused. Incoming material is restored after successful
and failed assessments. `BinaryLimitEvaluation3` is deliberately NOT a raw
optimizer design and carries no differentiable density gradient.

`scan_binary_response_limits` accepts an explicit threshold list and
`BinaryScanOptions3`: family/sharpness policy, void floor, candidate-count cap,
solid-volume cap and normalized feasibility tolerance. It filters the source
once, then evaluates candidates in caller order. Force windows need not be
monotone in threshold; no bisection or interpolation of feasibility is used.
Duplicate masks are still solved. Equal-volume ties keep the earliest complete
candidate. Only the current and incumbent field families are retained, along
with scalar candidate summaries. A complete sweep with no feasible candidate
is not proof that every possible binary design is infeasible.

The least-volume completed candidate must meet the smooth upper aggregate,
EVERY original signed response limit, and the actual numerical solid-volume
cap. A later stop retains that candidate and preceding scalar summaries in the
report, with the error in `outcome` and all consumed work in `work`. Invalid
input or failure of the initial filter returns an error before any candidate
report exists. The example prints retained endpoint fields before propagating
late errors and returns nonzero when no threshold succeeds. Tight force windows
may legitimately reject all tested thresholds even after a feasible gray fit.

The example retains the gray study's stiffness floor, uses a volume cap of one,
and permits no normalized response violation. Its CSV includes masks in active
cell order, exact stiffness scales, original values/bounds/residuals and fresh
master-node displacement fields. It preserves the gray source operator.

This is a fixed-grid **ersatz-material** check, not geometry reconstruction:
void cells retain nonzero stiffness, reference forces and support patches are
not deleted, and cut quadrature is unchanged. Lowering the floor can change
responses or defeat the numerical solve. Passing at one floor does not prove
zero-void behavior, connectivity, continuum safety, printable feature sizes or
manufacturability. A body-fitted or newly cut removed-cell model needs its own
geometry/load/support admission and independent final solves. None is claimed
by this endpoint API.
