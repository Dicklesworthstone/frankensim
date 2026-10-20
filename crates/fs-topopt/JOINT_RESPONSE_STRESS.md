# Minimum material with response requirements and a stress cap

```sh
cargo run -p fs-topopt --features cutfem-marquee --release \
  --example response_stress_sdf3 -- 80 500000 0.0012 0.008 0.1
```

Use the repository's DSR/RCH lane where available. Arguments are accepted
updates, cumulative Krylov allowance, minimum/maximum actuator reaction, and
stress cap. The fixed implicit slab uses E=1, nu=0 and prescribed axial strains
+/-0.02 with independent opposite body loads. All values are declared in
reference/dimensionless units, not calibrated experimental or allowable data.

A pure prescribed-displacement stress cap can reward removing material. This
problem also requires a force window and displacement ceiling in BOTH
independent experiments. It does not substitute a response-fitting objective
for a physical inequality or discard an inactive original bound.

## One joint numerical constraint

`CutDensityStudy3::evaluate_response_stress_limits` is implemented in
`response::limits::strength`. Pass the original
`ResponseLimitCase3` family and `ResponseStressOptions3`. Every experiment must
contain at least one original response limit; its identical external load and
optional embedded motion enter the stress calculation. Stress always selects
`StressMeasure3::SampledPeakBound`: the unweighted p-norm includes both relaxed
and physical stress at every retained sample of every experiment. There are
no load weights that can suppress a stress-constrained case.

For original response residuals c_i, the existing upper log-sum-exp is R.
Let S=stress_aggregate/stress_limit-1. Applying the SAME log-sum-exp sharpness
to [R,S] gives the upper aggregate of all c_i and S in exact arithmetic.
Its derivative is a two-component softmax combination of the original exact
response gradient and the stress gradient divided by the cap. Neither the
normalizing subtraction used by an average nor objective priority weights are
introduced. Original measurements, signed bounds, stress aggregate and both
separately evaluated maxima remain visible and are all checked for feasibility.

The response tangent includes direct reaction-material terms. The stress
tangent includes direct relaxed/physical stress terms. Both include stiffness
AND prescribed-motion lifting derivatives and the original filter/projection/
SIMP pullback. The stress p-norm and outer log-sum-exp are conservative smooth
approximations of their finite sets; a narrow feasible window may be excluded
by smoothing. None of these bounds is a continuum stress guarantee.

The reference implementation intentionally calls both existing evaluators.
That means two material preparations and two primal families, plus their
component adjoints. This is not a fused solve or a performance improvement.
All setup, filter, primal and adjoint work uses ONE caller-owned control,
including completed response work before failed stress evaluation. Both phases
restore incoming scales on success and returned errors. No partial combined
result is published; callbacks must be pure reference laws and panics are not
caught. The retained component adjoints are not mislabeled as joint adjoints.

## Joint material search

`strength::design::ResponseStressDesignStudy3` delegates bound projection,
line search, multipliers, penalties and spectral steps to the existing
linear-storage projected augmented-Lagrangian engine. Its objective is actual
projected material volume. `run` segments preserve search state and budgets;
zero steps use cached derivatives. No dense constraint Jacobian is allocated.

The original stress-design ersatz-turnover admission is reused for the raw
floor and every evaluated projected free cell. Protected physical voids retain
the same explicit exception. Other fixed material labels and density mapping
are unchanged. This does not prove stiffness coercivity or physical connectivity.

Accepted AL iterates can violate constraints and increase material volume.
`best_feasible` archives the least-volume accepted point that passes the joint
aggregate, all original response limits, and sampled stress checks separately.
`into_best_feasible` explicitly installs that incumbent and refuses when none
exists. A late cancellation keeps the completed accepted transition and its
matching material/fields. Latest-iterate KKT diagnostics do not become a
stationarity claim for an earlier selected incumbent.

The example prints accepted history, every original bound/value/residual,
per-case stress maxima and external-force work, plus raw/projected densities.
It emits retained evidence before returning late errors. It returns an error
when no feasible incumbent was accepted; an iteration limit with a feasible
incumbent is reported as that limit, not as convergence.

This increment is fixed-grid. Response-only checkpoint, binary-endpoint and
DWR APIs do not represent the joint objective; joint persistence/refinement
requires explicit future integration. There is no follower load, shape
sensitivity, global optimality, actuator-energy or manufacturing certificate.

## Verification scope

Ten real-API Rust regressions cover component parity/work, complete mixed-load
finite differences, hidden-limit refusal, protected regions, second-phase
failure, admission, material descent, split-run replay, accepted-state
cancellation and impossible constraints. They were authored but not executed
in the editing environment: DSR, its documented absolute path, RCH and Cargo
all failed to start because the executables were absent.

An independent NumPy/SciPy 1-D Nitsche rod reference passed 36 coordinate checks
over nine cases (maximum relative discrepancy 9.17e-8). A separate SLSQP solve
reduced volume fraction 0.69125 to 0.40846 while satisfying all original signed
limits and the sampled-stress aggregate. That is not execution of Rust, 3-D
CutFEM or the repository optimizer, and provides no performance evidence.
