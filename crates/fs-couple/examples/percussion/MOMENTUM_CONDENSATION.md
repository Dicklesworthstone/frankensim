# Exact momentum elimination for nonlinear percussion

The existing `--analytic-newton` execution now selects momentum condensation
for eligible models without radiation memory. No new instrument, geometry,
contact law, material preset or time integrator is introduced. The ordinary
finite-difference path is unchanged.

```sh
cargo run --release -p fs-couple --example percussion -- \
  splash 4096 --analytic-newton > splash.csv
cargo run --release -p fs-couple --example percussion -- \
  drum 4096 --analytic-newton > drum.csv
```

These remain physical reference experiments, not calibrated drum sounds or
real-time-qualified rendering. The existing geometry, frequency windows,
force/slope/energy limits, acoustic model and timestep still apply.

## What is eliminated, and what is not

For separable kinetic energy, the mechanical Hamiltonian has the form
`H(q,p,z) = p^T p/2 + V(q,z)`, with material memory `z`. The analytic Gonzalez
Newton matrix is a base matrix plus the full rank-one energy correction.
Its independent momentum coordinates can be solved as scalar leaves. Every
nonlinear displacement and material-memory coordinate remains in the dense
border, together with one auxiliary coordinate for the Gonzalez correction.
All momenta are recovered before checking the original full equation.

With `m` canonical modes and `h` extra memory states, this changes the dense
factorization dimension from `2*m+h` to `m+h+1`, plus `m` scalar leaf solves.
For example, 32 modes and six memory states retain all 70 physical states but
use a 39-coordinate dense border. This is an algebraic dimension, not a measured
speedup: full Jacobian construction, Schur updates, residual checks, geometric
Hessians and nonlinear iterations still cost work. Full dense fallback scratch
also remains allocated.

Hertz penetration, shell/membrane geometric stiffness, cavity overlaps, Kelvin
creep, conditioned felt and the accepted-history/energy gate are unchanged.
No small coupling, physical mode, material history or radiated-pressure input
is removed to obtain the smaller factorization.

## Selection and fallback

Preparation selects the momentum plan only with at least two modes, no
radiation memory, no nonlinear dissipative port, and an exactly diagonal
momentum block of `J-R`. Thus a spatial muffler is not silently replaced with
independent modal damping. This is only a cold eligibility screen. The actual
analytic Jacobian is checked on **every** Newton update; any nonzero coupling
between selected leaves, including a subnormal coefficient, requires full LU.
Singular leaf factors or failed original-equation backward error also use full
LU. Bounded iterative refinement and transactional publication are retained.

Radiation-loaded systems keep their existing independent-acoustic-pair plan.
Mechanical momenta exchange power with these pairs, so combining both sets as
independent leaves would be incorrect. This change does not do so, and does not
remove radiation loading to make a model eligible.

Library users can select dense comparison execution without resetting a strike:

```rust,ignore
let mut prepared = system.prepare_analytic()?;
prepared.set_momentum_condensation(false)?; // unloaded model: original full LU
prepared.set_momentum_condensation(true)?;  // exact, checked momentum plan
```

This setter is cold: it allocates, preserves accepted state/history/time, and
refuses radiation-loaded models without changing their existing plan. Explicit
enabling on a cross-damped or nonlinear-loss model may attempt condensation,
but cannot bypass the exact structure check. Use the existing
`set_radiation_condensation` for loaded comparisons; it replaces the numerical
plan rather than silently combining selections. Switching analytic Newton off
ignores either plan. Changing numerical traversal can change rounding and
iteration counts; bit identity with full dense execution is not promised.

`condensed_newton_dimension()` reports the selected border (including the energy
correction), and `newton_linear_solve_counts()` distinguishes condensed solves
from full solves/fallbacks. Neither diagnostic measures elapsed time.

## Validation boundary

The shared owner is `fs-phs::StepWorkspace::set_condensed_leaves(pairs, scalars)`.
Its existing pair-only API and pair regressions remain. New regressions cover
mixed/permuted/scalar-only layouts, energy-correction borders, changed coupling,
subnormal refusal, cancellation, retry, and nonlinear work/loss. Percussion
regressions use actual reduced curved-shell and stretching-head owners, Hertz
contact, conditioned felt and Kelvin memory, reciprocal volume work, spatial
loss, and unchanged radiation selection.

```sh
cargo test --release -p fs-phs --lib prepared::condensed
cargo test --release -p fs-couple --lib render::plate::impact::prepared::momentum
```

Native compilation, all regression results and deadline benchmarks must be
checked in an equipped environment. Independent equation checks are not Rust
execution, measured instrument validation, allocation instrumentation or a
real-time performance certificate. Refer to [ANALYTIC.md](ANALYTIC.md) for the
unchanged tangent, nonlinear solve and acceptance contracts.
