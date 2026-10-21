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

The joint path uses the same component kernels under ONE density transaction
and ONE material preparation. It computes the filter/projection once, solves
each experiment once, and shares those fresh primal fields and the volume
gradient with the stress kernel. Primal fields must pass the original explicit
equilibrium-residual check first. No field or factor is cached across density
trials, and there is no public API accepting arbitrary cached equilibria.

Response and stress retain their distinct adjoints and material gradients.
Prescribed-motion lifting is still reintegrated for stress RHS admission and
its derivative; it is not replaced by a frozen force. Both component field
families remain independently owned in the result, so this does not claim
zero-copy or reduced peak-memory operation. Neither component is skipped when
its outer softmax weight is small.

For N experiments the number of linear solves is 3*N+4, rather than 4*N+6:
N primal, 2*N component adjoint, one density filter, two constraint pullbacks
and one volume pullback. Two experiments therefore start 10 solves instead of
14. Only one preparation's Galerkin work is charged. These are algorithmic
counts, not native timings, total arithmetic or a measured wall-clock speedup.

All setup and solve work uses ONE caller-owned control, including completed
response work before failed stress evaluation. The joint owner restores
incoming scales on success and returned errors. No partial combined result is
published; callbacks must be pure reference laws and panics are not caught.
Component adjoints remain component adjoints, not relabeled joint adjoints.

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

This increment remains fixed-grid. The joint owner now exposes in-memory checkpoint/restore for the real projected-AL state. It preserves multiplier, penalty, spectral schedule, evaluation accounting, accepted history and the separate least-volume feasible incumbent, but never trusts saved displacement, reaction, stress, adjoint or gradient fields. Restoration re-solves the accepted endpoint and a distinct incumbent under the caller-bound original geometry, experiment family and policy before installing saved search state.\n\nThose endpoint evaluations are reserved against the original optimizer allowance before physics starts, and all solve/preparation work remains on the supplied cumulative SolveControl. Changed bounds, stress cap, motion/load laws or material mapping that change retained evidence refuse restoration. Callback identity remains a caller obligation; matching endpoints do not prove arbitrary laws globally equivalent. No byte/disk format or adaptive-grid migration is introduced.\n\nBinary-endpoint checks are separate physical post-processing and DWR APIs still do not represent the joint objective. Adaptive joint refinement remains future integration. There is no follower load, shape sensitivity, global optimality, actuator-energy or manufacturing certificate.

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

Three additional regressions check exact shared-versus-standalone fields and
gradients for mixed actuator/force experiments, a strict solve allowance that
cannot pay for duplicate equilibria, and once-per-density two-level/recursive
setup accounting. The existing parity and cancellation tests also exercise
the shared path. All existing optimizer consumers use this path without a new
mode or changed stopping policy. Native execution of these new tests remains
unverified: DSR (including its documented absolute path), RCH and Cargo are
unavailable in the editing environment. Source comparisons verify that the
component arithmetic was preserved; they are not Rust compilation evidence.


Three focused checkpoint regressions cover resumed trajectory equivalence with charged endpoint reconstruction, changed-requirement rollback, and malformed feasible-incumbent metadata. They are authored source evidence until the native Rust target runs.
