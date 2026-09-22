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
would differentiate a different equation. Each Jacobian column now uses an
analytic Hessian-vector action rather than two complete finite-difference
residual probes. Momentum-only directions skip geometric potential work;
no numerical threshold drops a small physical coupling.

All original timestep, Newton, energy, slope, felt-densification and lifetime
budgets remain. Failed trials publish neither state nor history. Library users
can call `ImpactSystem::prepare_analytic()` or toggle
`PreparedImpactSystem::set_analytic_newton(bool)` without resetting a ringdown.
The `drum-modal`/snare image has its own joint contact solver and is not silently
converted. Existing vented-exterior-audio restrictions are unchanged.

The Jacobian and LU are still dense. This removes a substantial class of
repeated nonlinear evaluations but is **not a measured speedup, real-time
qualification or bandwidth/convergence certificate**. Hard impacts and branch
transitions can still exhaust Newton; the original physical refusal is retained.
It does not supply missing measured geometry, damping, hand contact, full-band
cymbal modes or two-way radiation loading. These limits need separate work.

Focused checks:

```sh
cargo test --release -p fs-phs --lib prepared::analytic
cargo test --release -p fs-couple --lib render::plate::impact::tangent
cargo test --release -p fs-couple --example percussion analytic -- --test-threads=1
```
