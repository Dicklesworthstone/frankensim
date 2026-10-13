# Analytic Newton for nonlinear percussion

`--analytic-newton` selects allocation-free analytic Newton in the existing
prepared Gonzalez time owner. It implies prepared nonlinear execution; it can
also accompany `--prepared-nonlinear`. Without either flag the reference path
is unchanged. `--prepared-nonlinear` alone retains finite-difference Newton.

```sh
# Actual curved shell, hysteretic stand felt and physical stick impact.
cargo run --release -p fs-couple --example percussion -- \
  splash 4096 --analytic-newton --strike-speed-m-s 3 > splash-analytic.csv

# Both stretching heads, reciprocal distributed air, two sticks and a muffler.
cargo run --release -p fs-couple --example percussion -- \
  drum-stretch 4096 --analytic-newton --cavity-modes \
  --strike-speed-m-s 4 --strike-position-m 0.06 0.01 \
  --second-stick-position-m -0.05 0.02 --second-stick-speed-m-s 2.5 \
  --muffler batter 0.08 0.01 0.4 > drum-analytic.csv

# The unchanged BEM microphone path consumes the same physical trajectories.
cargo run --release -p fs-couple --example percussion -- \
  splash-mic 4800 20 --analytic-newton > splash-analytic.wav
```

This is an execution choice, not a linearized instrument. fs-plate differentiates
its actual curved-shell metric strain and statically relaxed membrane energy,
including mixed-mode and geometric stiffness. fs-dcontact contributes each
active unilateral contact's tangent. fs-material supplies the conditioned felt
branch tangent; series Kelvin coordinates remain in the joint solve. Every
cavity and neck spring uses the same original overlap rows. Spatial viscous
mufflers still enter the resistance operator rather than the storage Hessian.
No geometry, material history, mode, contact or pressure transfer is duplicated.

fs-phs differentiates the **full Gonzalez discrete gradient**, including its
energy correction and roundoff guard. Using only half the midpoint Hessian
would differentiate a different equation. Each Jacobian column uses an
analytic Hessian-vector action rather than two complete finite-difference
residual probes. Momentum-only directions skip geometric potential work;
no numerical threshold drops a small physical coupling.

All original timestep, Newton, energy, slope, felt-densification and lifetime
budgets remain. Failed trials publish neither state nor history. Library users
can call `ImpactSystem::prepare_analytic()` or toggle
`PreparedImpactSystem::set_analytic_newton(bool)` without resetting a ringdown.
The `drum-modal`/snare image has its own joint contact solver and is not silently
converted. Existing vented-exterior-audio restrictions are unchanged.

The full Jacobian and dense fallback scratch remain. Radiation-loaded analytic
execution additionally condenses independent acoustic pairs as described below.
Neither change is a **measured speedup, real-time qualification or bandwidth/
convergence certificate**. Hard impacts and branch transitions can still exhaust
Newton; the original physical refusal is retained. These execution choices do
not change geometry, materials or the selected acoustic model and its limits.

Focused checks:

```sh
cargo test --release -p fs-phs --lib prepared::analytic
cargo test --release -p fs-couple --lib render::plate::impact::tangent
cargo test --release -p fs-couple --example percussion analytic -- --test-threads=1
```

## Radiation memory factorization

Loaded analytic percussion now eliminates independent acoustic pairs while
retaining the complete Gonzalez correction in a scalar border. All physical
states and histories remain. See `CONDENSED_RADIATION.md` for exact algebra,
dense fallback, comparison controls and the unmeasured performance boundary.

## Reuse geometric linearizations across Newton columns

Prepared analytic execution now assembles each shell or stretching-film Hessian
once at a trial displacement and reuses its columns. The original fs-plate
strain factors supply both material and geometric stiffness; negative curvature
and every mixed-mode coefficient remain. A condensed n-mode film takes
O(facets*n^2) work to assemble its full Hessian, rather than n separate
O(facets*n^2) directional evaluations. Ordinary shell assembly remains
O(facets*n^2), but evaluates strain/stress only once per facet across columns.
Neither operation count is a whole-step benchmark.

Reuse requires the **identical displacement bits of that body** and the same
immutable mechanical owner. No approximate cache key, cross-model reuse, frozen
previous-tick stiffness, modal deletion or state linearization is used. The full
nonlinear residual, Gonzalez correction, force derivatives and acceptance gates
are unchanged. A rejected trial can populate numerical scratch, never accepted
motion or material state. A changed displacement refreshes the affected body;
changing only momentum does not invalidate a potential Hessian.

Only the immutable shell/film potential is cached. Contact branches, felt
conditioning, Kelvin strain, material-relaxation memory, gas-film state and
radiation terms remain in the original live composed action. Thus unchanged
head displacement does not incorrectly freeze a changed felt unloading tangent.
Strings and other body types retain their original directional actions.

Preparation reserves sum(n_body^2) matrix entries plus displacement and
three-component strain-gradient scratch for the shell/film bodies. This is at
most the existing 256-mechanical-mode ceiling squared; it does not duplicate
an acoustic-memory Jacobian or a shell mesh. These routines allocate no memory
while stepping. The reference finite-difference path does not use the matrices.

`set_batched_body_tangents(false)` selects the previous directional evaluation
for numerical comparisons; changing it does not reset a strike, history or
clock and does not allocate. `body_tangent_work_counts()` reports body matrix
assemblies and matrix actions in the last solver attempt, not timing or accepted
physical work. For internal substeps it reports the last internal attempt.
Matrix assembly changes floating-point accumulation, so trajectories are not
promised bit-identical to the old directional path. Same-path retries and
repreparation retain deterministic operation order.

Focused regressions compare actual reduced shell/film actions, changed body
points, live felt history, Hertz-driven stretching-head trajectories, complete
energy balance and refused-step/reprepare replay. Native execution of these new
regressions remains required; independent mathematical checks are not Rust tests
or a real-time performance qualification.
