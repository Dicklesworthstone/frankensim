# Keep radiation physics, eliminate the large acoustic factorization

`--radiation-feedback --analytic-newton` now prepares exact acoustic-pair
condensation automatically for the nonlinear percussion image. The existing
splash, drum, snare and paired-cymbal commands retain their admitted geometry,
force programs, contact/material laws, acoustic pole models and clocks. No
extra command-line flag or alternate physical model is required.

```sh
# Same physically loaded command; a smaller Newton factorization, not fewer poles.
(set -C; cargo run --release -p fs-couple --example percussion -- \
  drum-stretch-mic 4800 20 --cavity-modes \
  --radiation-feedback --analytic-newton --impact-substeps 8 511 \
  --microphone-right -0.08,0.05,0.35 > loaded-drum.wav)
```

This is a supported invocation, not a claim that this complete render has
finished, is calibrated, or can run in real time. Radiation-fit, power,
spatial-resolution, state-count and mechanical-validity refusals still apply.

## What gets smaller

For `b` original mechanical/material scalar states and `p` acoustic oscillators,
the previous analytic path factored a dense `(b+2p)` square matrix at every
Newton update. The new path factors a `(b+1)` square border and solves `p`
independent 2-by-2 systems through the existing `fs-la::LuWorkspace` owner.
It then recovers ALL `b+2p` state increments. The acoustic coordinates still
store energy, dissipate it, exert reciprocal reaction, and survive between
samples. They have not become a delayed external force or an output filter.

The extra border unknown is essential. Gonzalez's discrete-gradient energy
correction couples even otherwise independent acoustic pairs. Writing the
complete Jacobian as `A = B + u v^T`, the implementation introduces
`y = v^T delta_x`, solves the bordered system, and eliminates only independent
pairs in `B`. The full derivative of nonlinear contact loss, including its
state and effort dependence, is retained. This is algebraic elimination of
the existing Newton equation, not a new time integrator or midpoint substitute.

The factorization's dominant work changes from cubic in `b+2p` to cubic in
`b+1` plus pair/border work proportional to `p*(b+1)^2`. Full Jacobian assembly,
Hessian actions, residuals, energy checks and the original dense fallback
scratch remain. Memory use is not claimed smaller. The 132-state core test
asserts that all its Newton updates use a 5-by-5 border, while evolving and
checking the original 132-state equation against the full dense solve.

## Exact structure and safe fallback

Every update checks the pair-block structure with EXACT zero comparisons.
An arbitrarily small nonzero coupling between two eliminated blocks is not
removed. Unexpected structure, a singular leaf/border solve, nonfinite
elimination, or unresolved full-equation backward error uses the existing full
dense LU on the SAME complete analytic matrix. Before that fallback, the
solver can correct roundoff with at most two solves of the same condensed
system against the original full-equation residual, as described below. The fallback does not remove
radiation, change contact strength, relax an energy tolerance or skip time.

Using a scalar border also avoids requiring the entire uncorrected `B` to be
invertible. The regression suite includes a singular `B` with an invertible
corrected Jacobian. Both numerical paths still use the original nonlinear
iteration budget, residual gate and physical work/energy admission.

The elimination plan uses actual stored scalar addresses. Radiation can be
attached before OR after hereditary material memory; it need not be the final
state block or start at an even scalar index. The original mechanical force
ports, cavity addresses, felt histories and receiver source rows do not move.

## Library controls and work diagnostics

`ImpactSystem::prepare()` prepares the plan when radiation memory exists;
`prepare_analytic()` activates analytic Newton. Turning analytic Newton on or
off later retains accepted motion, material history and time. Finite-difference
calls ignore the plan and preserve their original arithmetic.

`PreparedImpactSystem::set_radiation_condensation(false)` selects the original
full dense analytic factorization for comparison; `true` rebuilds the plan.
Changing the plan allocates cold scratch and must be done outside a hard-real-
time callback. It does not reset physical state. With no load, or a completely
zero load that the existing owner omits, no condensation plan is installed.

`condensed_newton_dimension()` reports the planned border size when analytic
condensation is enabled. `newton_linear_solve_counts()` returns the condensed
and full dense solves from the most recent solver call. On a substepped image,
that is the last internal attempt, NOT the sum across the entire output tick.
The low-level equivalents are `StepWorkspace::set_condensed_pairs`,
`condensed_dimension`, and `linear_solve_counts`.

Preparation allocates; elimination adds no step-time allocation. All state and
history still pass through the existing accepted-step gate. Cancellation or a
failed complete output tick restores acoustic history together with mechanics,
felt, and hereditary memory, and consumes no staged player force.

## Verification boundary

Six core regressions cover full-matrix equivalence with permuted pairs, the
Gonzalez/nonlinear-loss derivative, a 132-state evolution, exact retry,
subnormal-coupling fallback and unchanged dense/finite-difference execution.
Four impact regressions include real unilateral contact, felt/Kelvin and
Maxwell history, both attachment orders, substep rollback and zero-load parity.
Two executable regressions use actual curved shells with two sticks and stand
felts, and stretching heads with enclosed air and material memory. Those two
use a declared passive test load to isolate the numerical comparison; existing
radiation-feedback tests separately exercise BEM fitting and pressure output.

```sh
cargo test --release -p fs-phs --lib condensed
cargo test --release -p fs-couple --lib render::plate::impact::prepared::condensed_tests
cargo test --release -p fs-couple --example percussion condensed -- --test-threads=1
```

Authored tests and matrix arithmetic alone are not a native passing result.
No wall-clock speedup, audio fidelity, real-time deadline or reduced-memory
claim is made without measurement. Stationary-reference radiation and all
previously declared cymbal/drum modeling limitations remain unchanged.


## Sparse directional assembly and numerically balanced elimination

Loaded analytic execution now uses both row and column traversal of the
CURRENT `J-R` operator. Rows still evaluate the original residual. Columns
apply the operator to an analytic Hessian direction, skipping only exactly
zero direction components. Every nonzero coefficient, including subnormal
couplings, is retained; an updated operator rebuilds both traversals before
use. Per-row sums still receive terms in increasing source-index order.
The unselected dense analytic and finite-difference arithmetic is unchanged.
Finite-difference calls do not build the column traversal.

This matters because the uncorrected Hessian of each independent acoustic
state has only one nonzero component. The old traversal multiplied every
operator edge by a mostly zero vector for every acoustic column. The new
traversal visits only edges incident on nonzero components. Nonlinear storage,
contact-loss tangents and the COMPLETE rank-one Gonzalez correction remain.
No declarations about a permanently diagonal Hessian are trusted: the actual
computed direction supplies its support on every call.

The 132-state regression has 326 operator edges. It compares 43,032 scalar
operator products for the original dense analytic assembly with 652 for the
new assembly, including the rank-one direction. This 66-fold reduction is
ONLY in that counted operation, not a measured overall speedup. Hessian-vector
callbacks, residual evaluation, full Jacobian storage, state/history checks,
Schur products and factorization have not disappeared.

The scalar border is also balanced using reciprocal powers of two on its
rank-one factors. Original factors remain untouched for the full-matrix
check and dense fallback. Scaling is discarded if any coefficient would not
round-trip exactly, so it cannot erase underflowed entries. A zero outer
product uses a dummy zero auxiliary equation instead of an arbitrarily large
irrelevant one. This changes no physical coordinate or tolerance.

Subtractive cancellation during leaf recovery can fail the original backward
error gate even when the complete matrix is well-conditioned. At most two
residual corrections reuse the SAME eliminated operator. Every candidate is
checked against the original full floating-point matrix before publication.
The error allowance remains `128 * epsilon * (n+1)` times the absolute equation
scale. Cancellation or an unsuccessful correction retains the original dense
fallback and never publishes a partial correction. No extra Newton iteration,
mechanical step, force evaluation, or history update is counted as a correction.

`PreparedImpactSystem::newton_flow_product_count()` (low-level:
`StepWorkspace::jacobian_flow_product_count()`) reports scalar operator products
in analytic Jacobian assembly during the most recent solver call. It excludes
Hessian, residual, factorization and finite-difference work; counts include
rejected Newton attempts and saturate at `usize::MAX`. On a substepped image
this is still the last internal attempt, not an output-tick sum.

The new scratch is allocated during preparation: an optional column index
buffer reserved for a same-size dense operator, plus linear-size balance,
residual, correction and flow vectors. It adds memory; it does not claim a
reduced-memory implementation. Refresh and execution stay within that capacity.

Five additional regressions cover extreme auxiliary scaling, exact preservation
of subnormal coefficients, correction/cancellation at a tiny leaf pivot, bitwise
column/row traversal parity after topology changes, and complete Jacobian
parity with the 132-state operation count. Independent arithmetic checks are
not Rust test execution or a real-time deadline measurement.
