# Worst-load robust CutFEM design

```sh
cargo run --release -p fs-topopt --features cutfem-marquee \
  --example worst_load_sdf3 -- 200 250000
```

This fixed-grid example minimizes a smooth maximum of independent **eroded**
load compliances under the **dilated** material-volume cap. It uses the real
CutFEM filter/projection/SIMP chain and existing projected augmented-Lagrangian
optimizer. It does not replace the weighted-sum robust or adaptive commands.

## What is optimized

Choose a positive reference compliance `r_i` for every declared load. These
values have compliance units and are performance scales, not probabilities.
With `s_i = C_i / r_i`, temperature `t > 0`, and `m = max(s_i)`, the objective is

```text
J = m + t * ln(sum_i exp((s_i - m) / t))
dJ/dC_i = softmax_i(s/t) / r_i
```

In exact arithmetic, `m <= J <= m + t*ln(number_of_loads)`. The reported bias
bound is an algebraic smoothing statement, not a floating-point enclosure or
continuum certificate. Smaller temperatures approximate the maximum more
closely but sharpen the optimization problem. All references and temperature
remain fixed throughout a session. Duplicating a load changes the smooth
objective even when it leaves the exact maximum unchanged.

`LoadEnvelope3` admits the references and a caller-specified case-count limit.
`evaluate_robust_envelope` solves every load independently, computes its current
normalized compliance, and combines the exact self-adjoint contractions before
one filter-transpose pullback. Weights are recomputed at every design; they are
not frozen. Opposite forces cannot cancel as they would in a summed load.
The physical result retains each eroded displacement and all three density
fields. Its embedded weighted-sum compliance/gradient is a **uniform-mean
diagnostic**, not the objective and gradient in the envelope result.

## Constrained optimization and accepted results

`envelope::optimize::EnvelopeDesignStudy3` checks the baseline envelope and
volume gradients before starting the existing projected AL. Its single
inequality is `dilated_volume / volume_limit - 1 <= 0`. A reference is NOT an
additional hard per-load constraint: a returned design may exceed it. Inspect
the reported normalized maximum and every individual ratio.

The best feasible incumbent has the smallest **smooth envelope** among accepted
material-feasible designs. The actual maximum is reported separately; its
monotonic decrease is not guaranteed by a decrease in the smoothed objective.
Terminal KKT diagnostics belong to the final AL iterate, which may differ from
the retained incumbent. Refinement, material continuation and a changing
reference/temperature policy require a new appropriately defined problem; the
existing weighted-sum adaptive driver is not silently relabeled as minimax.

`run` can be segmented without resetting AL state, evaluation counts or the
caller's cumulative physical budget. Gradient probes, rejected trials and
interrupted evaluations all remain charged. Late cancellation after an
accepted step preserves its paired fields. This is in-memory continuation,
not a durable checkpoint or cross-process restart. The borrowed study's
incoming nominal stiffness remains unchanged on returned outcomes.

The example explicitly normalizes each case by its own initial eroded
compliance. Those initial solves are charged, and the references are printed
and fixed once. It exports accepted history, actual per-load ratios, the three
density fields and separate eroded displacement fields. It prints its retained
prefix before propagating a late stop. An iteration cap is not convergence.

## Scope and validation

The model requires independent density-invariant dead loads and homogeneous
supports on fixed geometry. It covers only the declared finite load family;
it does not implement CVaR, chance constraints, distributional uncertainty,
continuum guarantees, manufactured length scales or global optimality.

```sh
cargo test --release -p fs-topopt --features cutfem-marquee --lib sdf3::robust::envelope
cargo check --release -p fs-topopt --features cutfem-marquee --example worst_load_sdf3
```

Native execution was unavailable in the editing environment (Cargo/rustc
absent). Independent NumPy/SciPy gradient and constrained-design checks are
formulation evidence only, not execution of the Rust code or CutFEM tests.
