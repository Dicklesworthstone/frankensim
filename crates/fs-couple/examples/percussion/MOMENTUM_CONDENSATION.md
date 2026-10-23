# Exact momentum elimination for nonlinear percussion

`--analytic-newton` now prepares a smaller dense solve even when the model
includes spatial cross-damping, nonlinear contact loss, and reciprocal
radiation memory. It retains every mechanical, material and acoustic state.
No new geometry, material preset, contact law or time integrator is introduced.
The ordinary finite-difference path is unchanged.

```sh
cargo run --release -p fs-couple --example percussion -- --drum-spec crates/fs-couple/examples/percussion/estimated_drum.fsd \
  snare 4096 --head-stretching --cavity-modes --analytic-newton \
  --impact-substeps 8 511

# Same loaded instrument, NOT a reduced pole bank or an output-only filter.
(set -C; cargo run --release -p fs-couple --example percussion -- --drum-spec crates/fs-couple/examples/percussion/estimated_drum.fsd \
  drum-stretch-mic 4800 20 --cavity-modes --radiation-feedback \
  --analytic-newton --impact-substeps 8 511 > loaded-drum.wav)
```

These are implemented invocation paths, not completed render or real-time
benchmark results. The examples remain estimated reference instruments, with
their existing bandwidth, mesh, slope, contact, energy and radiation-fit limits.

## Solve a smaller system without removing physics

With `m` canonical modes, `h` material memories and `p` radiation oscillators,
the full state has `2*m+h+2*p` scalar coordinates. The joint plan factors a dense
border of size `m+h+1`, plus `m` scalar pivots and `p` two-by-two pair factors.
The old loaded path retained both displacement and momentum in its `2*m+h+1`
border. An unloaded cross-damped or nonlinear-loss model previously required
full dense LU. Neither restriction is needed for the new kinematic pivots.

For example, 32 canonical modes, six memories and 64 acoustic oscillators keep
all 198 states, but use a 39-coordinate border instead of the pair-only border
of size 71. This is an algebraic dimension, NOT a measured speedup. Full
Jacobian assembly, Hessian actions, response construction, residual checks and
nonlinear iteration still cost work. Full dense fallback scratch remains
allocated; total memory use is not claimed smaller.

The Hamiltonian is `H(q,p,z)=p^T p/2+V(q,z)` in the existing mass-normalized
coordinates. The analytic Gonzalez Jacobian is `A=B+u*v^T`. A pivot uses the
**displacement equation row** and its **momentum unknown column**, rather than
assuming the momentum resistance matrix is diagonal. The momentum equations
remain in the dense border, with every spatial damper and contact-loss tangent.
The complete Gonzalez correction is retained using `y=v^T*delta_x`.

Radiation equations depend on these momentum unknowns. Their two-coordinate
blocks therefore follow the momentum pivots and include the lower coupling
terms in forward substitution. The combined eliminated block is checked as
block lower triangular, NOT declared block diagonal. Both directions of the
mechanical/acoustic reaction remain in the original Newton equation. There is
no delayed acoustic force, partitioned time step or omitted wire coupling.

Every removed increment is reconstructed before the full original-equation
check and physical acceptance. Hertz/Hunt-Crossley contact, nonlinear
shell/head stiffness, snare modes, cavity overlaps, Kelvin creep, conditioned
felt, hereditary material energy, pressure source rows and force ports are
unchanged. State addresses do not move; acoustic memories may occur before or
after hereditary material states and at odd scalar addresses.

## Selection, fallback and controls

Multi-mode loaded systems automatically use the joint kinematic/acoustic plan.
One-mode loaded systems retain the existing pair-only plan. Unloaded models
with at least two modes use kinematic pivots for nonlinear loss or cross-damping;
the simpler diagonal-resistance models retain their established scalar-leaf
traversal. This avoids needlessly changing their arithmetic.

Every analytic update verifies the ACTUAL base Jacobian structure. A nonzero
upper inter-block coefficient, even subnormal, is not dropped. A singular pivot,
nonfinite elimination, or failed original-equation backward error invokes full
LU of the SAME Jacobian. Bounded refinement uses the original residual before
fallback. No force, state, tolerance or time step is changed to admit a solve.

```rust,ignore
let mut prepared = system.prepare_analytic()?;
// Unloaded model:
prepared.set_momentum_condensation(false)?; // full dense comparison
prepared.set_momentum_condensation(true)?;  // restore appropriate momentum plan
// Radiation-loaded model (use these INSTEAD of the unloaded controls):
prepared.set_radiation_condensation(false)?; // full dense comparison
prepared.set_radiation_condensation(true)?;  // restore joint plan, or one-mode pairs
```

The controls allocate only during preparation/reconfiguration, preserve accepted
motion/history/time, and leave the previous plan intact on rejection. The
unloaded setter still refuses radiation-loaded models. Turning analytic Newton
off makes execution ignore either plan. Different traversal can change floating
point rounding and iteration counts; bit identity with dense Newton is not
promised. `condensed_newton_dimension()` reports the planned border including
the energy-correction scalar. `newton_linear_solve_counts()` reports condensed
and full-dense solves in the last solver call, not a whole substepped sample.

## Focused checks

Shared `fs-phs` regressions cover separate row/column maps, lower couplings,
singular uncorrected matrices, exact structural refusal, cancellation and
forced cross-damped quartic motion. `fs-couple` tests compare the same complete
loaded state against dense execution in both material-attachment orders, retain
Hertz rebound/felt memory, and exercise simultaneous spatial damping/contact
loss. Example regressions retain curved shells, stretching heads, active snare
contacts, cavity modes and the complete twenty-strand bank.

```sh
cargo test --release -p fs-phs --lib prepared::condensed
cargo test --release -p fs-couple --lib render::plate::impact
cargo test --release -p fs-couple --example percussion condensed_impact_tests
```

Native tests were authored/updated but could not execute in the implementation
environment: `cargo` and `rustc` were absent. Independent Python/SciPy checks
ran 108 block matrices and two 400-step nonlinear loaded/unloaded trajectories;
those validate equations, not the Rust implementation, full-shell audio,
allocation instrumentation or a real-time deadline.
