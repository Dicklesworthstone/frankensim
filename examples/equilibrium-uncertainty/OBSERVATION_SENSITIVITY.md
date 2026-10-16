# Which parameters do the displacement experiments actually resolve?

`equilibrium_sensitivity` reads the existing mechanical model and design files,
evaluates all displacement observations at one explicitly named parameter point,
and returns their analytic Jacobians plus a local numerical information analysis.
No fitting, parameter perturbation, sampling, or alternate mechanical solver is
introduced. Each load case uses its original preload and tangent preparation;
each target uses the original response adjoint. Physical and derivative admission
remain mandatory, including the contact-switch exclusion margin.

## Two perfect fits, different information

From the repository root:

```sh
cargo run -p fs-cli --bin equilibrium_sensitivity -- \
  examples/equilibrium-uncertainty/sensitivity-parallel.model \
  examples/equilibrium-uncertainty/sensitivity-parallel.fit \
  --rank-relative-tolerance 0.00001 --max-observations 16 --max-adjoints 16 \
  --point-x left 0 --point-x right 0
```

Two parallel springs restrain one mass. The observations are `F/(k_left+k_right)`
at two loads. Both targets are matched exactly at the declared point, but the
weighted Jacobian has rank one: increasing one stiffness and decreasing the
other by the same amount leaves every observation unchanged. More repeats or
more load amplitudes on this same configuration cannot separate the two springs.
The weak parameter direction is proportional to `[1,-1]` in these equal scales.

Replace `sensitivity-parallel` with `sensitivity-independent` in both filenames
to observe two independently supported masses instead. The two stiffnesses then
have separate observations: the analytical Jacobian has full rank two and
condition number two at this point. The objective and its gradient are zero in
both examples. An objective gradient alone cannot distinguish these situations.
These are analytical fixture properties; native Rust execution remains unverified.

## Output and numerical meaning

The JSON includes exact source hashes, variable names and physical affine scales,
every case/target identity, physical predictions, raw `d_displacement_d_x`, and
weighted residual derivatives. The weighted residual is
`sqrt(weight)*(displacement-target_m)/scale_m`. Zero-weight targets keep their raw
physical derivatives but contribute zero weighted information. Columns always
follow design-file variable order; named command arguments may be reordered.
The result also records completed primal cases, attempted observation adjoints,
original equilibrium/adjoint residuals, and the count of response constraints
that this calibration query did **not** assess.

`information.numerical_rank` describes this weighted Jacobian only, at the supplied
relative singular-value threshold. The existing `fs-la` cyclic Jacobi solver acts
on the scaled Gram matrix; the domain adapter checks eigen residuals and
orthogonality. Singular-value ratios and matching parameter directions are ordered
from weakest to strongest. `jacobian_scale` is the largest absolute weighted
Jacobian entry removed before forming products, not an inferred noise standard
deviation. No numerical Gram inverse or parameter covariance is fabricated.

Gram construction squares conditioning. This bounded reference path admits only
1–32 variables and relative thresholds in `[1e-5,0.5)`. It does not resolve arbitrarily
small singular values. A mode too close to the threshold relative to the numerical
rounding/residual screen leaves rank `null`; rank deficiency or ambiguity leaves
condition number `null`. The numerical screen is not an outward-rounded bound.
Directions within a degenerate subspace do not have a unique physical basis.

Full local numerical rank does **not** prove global/structural identifiability,
calibrated material parameters, physical validity, active-constraint identifiability,
or a confidence/posterior interval. Fitting weights are not automatically noise
precisions. Coordinate scaling is explicit and changes the conditioning question.
No existing identifiability authority records or maturity claims are upgraded.

## Reusable library API

`EquilibriumDesign::evaluate_observations` accepts the existing `DesignControl`
and a separate `ObservationControl`. Row capacity is at most 1024; adjoint attempts
consume a cumulative allowance that can be extended but never reset by a retry.
The complete family is admitted before physics. Returned errors never expose a
partial Jacobian or refund attempted work. The immutable result provides
`gauss_newton_product(v, gate)` for `J^T J v` without a square matrix or new solves;
this omits residual-weighted second derivatives and is not the full Hessian away
from zero residual. `information(relative_threshold, gate)` uses no further physics.
The 32-column dense eigensolve polls before/after, not within its existing sweeps.

Focused native checks:

```sh
rch exec -- cargo test -p fs-couple --test equilibrium_observations
rch exec -- cargo test -p fs-couple --lib observations::information
rch exec -- cargo test -p fs-cli --bin equilibrium_sensitivity --test equilibrium_sensitivity
```

The implementation environment lacks `rch`, `cargo`, and `rustc`. These native
checks, compilation, rustfmt and Clippy have not passed here. Independent
physical-coordinate and NumPy matrix controls do not execute repository Rust.
