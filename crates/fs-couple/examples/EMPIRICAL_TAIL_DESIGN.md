# Physical design against empirical tail loss

```bash
cargo run -p fs-ascent --features equilibrium-design --bin equilibrium_fit -- \
  crates/fs-couple/examples/equilibrium-tail.model \
  crates/fs-couple/examples/equilibrium-tail.fit \
  --scenarios crates/fs-couple/examples/equilibrium-tail.scenarios \
  --tolerance 1e-9
```

A scenario file can now select tail risk by placing one optional record immediately
below `frankensim-equilibrium-scenarios-v1`, before `scenarios COUNT`:

```text
risk empirical-cvar 0.4
```

Without this record, the existing worst-case objective and output are unchanged.
Alpha must be finite and strictly between zero and one. This selection explicitly
assigns equal empirical mass `1/COUNT` to every realization. It does NOT infer
probabilities from tolerance endpoints or labels. Repeated offsets with distinct
names intentionally retain their mass; deduplicating them would change the problem.
No nonuniform weights or out-of-sample distributional guarantee are implied.

## What is optimized

For losses `J_i(x)`, the command minimizes

```text
eta + sum(z_i) / (COUNT * (1 - alpha))
```

subject to `J_i(x) - eta - z_i <= 0`, `-z_i <= 0`, every original physical
constraint in every realization, and the existing common nominal parameter bounds.
This is the finite Rockafellar-Uryasev CVaR formulation. The existing SQP engine
owns all steps, line searches and KKT checks; the existing equilibrium and adjoint
implementations supply the physical rows. There is no differentiated sort, smooth
maximum, second optimizer, resampling or discarded failed scenario.

The threshold and excesses are optimizer coordinates, not physical parameters.
The expanded dense admission is `3*n + 1 + COUNT*(3 + c)`, where `n` is the number
of physical decisions and `c` is the number of declared physical requirements.
The original KKT dimension cap is not enlarged. Physical evaluation and case-solve
budgets are unchanged, including the reserved final complete-family re-solve.

The empirical tail can have fractional mass. With four scenarios and alpha=0.4,
its mass is 2.4 realizations, not a rounded choice of two or three. Alpha=0.9 with
only four equally weighted scenarios reduces to the maximum loss; higher alpha
cannot manufacture evidence about rarer events than the supplied set represents.

## Editable analytical example

The model is one mass-normalized elastic coordinate with natural frequency
2 rad/s and an actuator/observation shape of 1. Its static displacement is `u=F/4`
metres. The target is 1 metre, and the dimensionless loss is `0.5*(u-1)^2`.
Three realizations use the nominal load and one uses a +4 N offset. These are
synthetic, declared design cases, not measurements or an inferred manufacturing law.

An independent analytical reference gives:

| Selection | Nominal load | Selected risk | Maximum scenario loss |
|---|---:|---:|---:|
| Default worst case | 2 N | 0.125 | 0.125 |
| Empirical CVaR, alpha=0.4 | 7/3 N | 35/288 | 49/288 |

The tail objective improves by accepting a higher single worst loss. That is the
explicit requested tradeoff, not a guarantee that CVaR is always preferable.
Force, penetration, travel and support-reaction requirements do not participate
in that tradeoff: they remain required in EVERY scenario, even one with zero
loss-tail multiplier. The contact regressions exercise this distinction.

## Reading a result

CVaR runs use `frankensim-equilibrium-cvar-fit-v1` output. They report alpha,
equal scenario mass, fractional tail mass, threshold, optimized excess slacks,
all tail and original physical multipliers, every physical prediction, actual
worst loss, stop attribution and all four original KKT residuals.

`objective` is the optimizer's threshold-plus-slack value. `cvar_upper_bound`
recomputes `eta + sum(max(J_i-eta,0))/(COUNT*(1-alpha))` from the final re-solved
losses rather than trusting the slacks. This is a numerical upper bound on the
finite-set CVaR at the returned design, NOT an exact minimized order-statistic
risk statistic. A stopped/stalled threshold need not minimize it; a positive
`tail_violation` can make the optimizer's objective overoptimistic. The reported
stop/convergence flag is retained rather than relabelled. The canonical exact
order-statistic reporting service remains in `fs-robust`, unchanged.

The native API is `ScenarioProblem::new(...).with_cvar(alpha)`, consumed by the
same `ScenarioEquilibriumStudy`. `cvar_alpha()` identifies the fixed profile and
`cvar_upper_bound(evaluation, threshold)` provides the explicit recheck functional.
Cancellation and failed physics retain the same accepted numerical checkpoint and
spent-work accounting, including all excess variables.

## Checks and scope

```bash
cargo test -p fs-ascent --features equilibrium-design --lib equilibrium::scenarios
cargo test -p fs-ascent --features equilibrium-design --bin equilibrium_fit --test equilibrium_scenarios
```

Six core regressions and five command/parser regressions were added. They include
fractional-tail analytic optima, full Jacobians, contact requirements outside the
tail, dimension/budget admission, cancellation/split continuation, strict input,
relocation and a budget stop with a nonoptimal threshold. Native execution was
not possible in the authoring environment (`cargo`/`rustc` unavailable). Separate
NumPy/SciPy calculations validate the analytical examples and formulation, not the
Rust implementation or deterministic math library. No global optimum, continuous
uncertainty coverage, probability-of-failure, physical calibration or real-time
claim is made. Fixed-basis static contact and existing activity-margin restrictions
still apply. No dependencies, lockfile edges or physical input schemas changed.

Reference: R. T. Rockafellar and S. Uryasev, Optimization of Conditional Value-at-Risk,
Journal of Risk 2(3), 21-41 (2000), DOI 10.21314/JOR.2000.038.
