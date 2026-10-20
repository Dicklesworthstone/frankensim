# Constrained heating schedules and scan geometry

`fs_ascent::scan_design::constraints::ConstrainedScanStudy` binds a
`ScanPowerStudy` to the existing fallible, resumable small-dense SQP engine.
The physical objective, full-history adjoint, source normalization, phase charts,
material/radiation feedback and physical acceptance criteria are unchanged.

```rust
use fs_ascent::scan_design::constraints::{
    ConstrainedScanStudy, ScanControlBound, ScanEnergyCap,
};

let mut fit = ConstrainedScanStudy::new(
    cx, study,
    vec![ScanControlBound {
        coordinate: 0, lower: Some(0.01), upper: Some(3.5), scale: 1.0,
    }],
    vec![ScanEnergyCap {
        start_time_s: 0.0, end_time_s: 2.0, maximum_j: 5.0, scale_j: 1.0,
    }],
    32,
)?;
let report = fit.advance_optimizer(cx, 1e-8, 120)?;
let constraints = fit.constraints_at(cx, fit.point())?;
let candidate_sources = fit.sources(cx)?;
```

Here `study` is the caller's previously constructed power-only or joint study;
`cx` is its real execution context. Construction spends one physical evaluation.
The final argument caps decision plus constraint dimension before any dense
constraint Jacobian or SQP matrix is allocated. It does not bound the underlying
thermal mesh or its separate retained-adjoint storage.

## Physical limits and tolerances

Control indices use the original point layout: selected powers first, followed
by selected geometry controls. Bounds use watts for powers and metres for center
values or half widths, not logarithms or dimensionless search coordinates.
Each declared bound may have a lower limit, an upper limit, or both. Two equal
limits are refused; leave a fixed variable unselected. Only selected authored
values are bounded: this does not constrain a different, unselected value or
impose a continuous trajectory envelope on a polynomial path.

Each inequality is scaled in its original physical units. A residual allowance
`tau` permits `tau * scale` watts/metres of violation for a bound, or
`tau * scale_j` joules for an energy cap. Report rows are lower then upper for
each bound in declaration order, followed by energy caps in declaration order.
`constraints_at` retains the actual bound values, total absorbed energy, signed
scaled residuals and row-major Jacobian, without another physical solve.

Joule caps include **every attached moving source**, including fixed sources
and unselected values. The input integral uses the original constant, ramp, or
held/linear-table interpolation, includes all knots and clamped tails, and
supports multiple independent time windows. The source API
`integrate_absorbed_power(cx, start, end, max_values)` also exposes the exact
piecewise integral and its authored-power weights in seconds directly.
Chebyshev power is refused by this integration path; no positivity or quadrature
policy for that representation is silently introduced.

These are domain-absorbed input joules. They are not electrical consumption,
incident optical energy, background volumetric heating or net storage after
boundary losses. A finite-step thermal solve may apply different temporal
quadrature to continuously varying inputs. The schedule integral is deliberately
not presented as that solver's delivered-energy audit.

## Optimization and continuation

SQP may visit constraint-infeasible but physically evaluable trial designs. It
does not clip sources and does not change the objective to hide constraint
violations. Native source-support, material, chart-smoothness, solve and energy
errors propagate. The limits are final-design constraints, not an actuator
controller or a restriction on every intermediate simulation. The caller's run
factory must still keep all non-source physical data unchanged across probes.

The wrapper owns its immutable limits and accepted SQP state. It cannot resume
that state against a replacement constraint set. Physical evaluation allowances
include failed simulations and persist through budget extensions. Cancellation
or callback failure preserves accepted optimizer fields; spent work is retained.
`extend_evaluation_budget` increases the cumulative allowance. `sources` exports
the current candidate without spending an evaluation; `replay` spends one to
reevaluate its original physical objective. Neither operation establishes
convergence or feasibility on its own.

Check `report.solution.converged`, the four native KKT residuals, and the physical
constraint report. Iteration/evaluation exhaustion or a stalled search is not
an infeasibility certificate. A converged constrained best fit can have a
nonzero target error; it is a local first-order result, not a global optimum.

## Executable synthetic cases

```bash
cargo run --manifest-path examples/scan_power_fit/Cargo.toml -- --limits
cargo run --manifest-path examples/scan_power_fit/Cargo.toml -- --limits 3 20
cargo test --manifest-path examples/scan_power_fit/Cargo.toml limits::tests
```

The unit-mass fixture starts at 150 J/kg on the 310 K latent plateau, with
enthalpy targets of 152 and 157 J/kg after two one-second heating intervals.
Without limits the analytical fit is 2 W followed by 5 W. With the default
3.5 W per-value ceiling and 5 J total ceiling, the constrained analytical fit
is 2 W followed by 3 W: the final target cannot be reached within that budget.
With a 3 W ceiling and loose 20 J budget the constrained fit is 3 W followed
by 3 W, not the 2/3 W obtained by clipping the unconstrained answer afterward.
The executable prints predictions, original objective, KKT residuals, actual
joules, violation and the exported held schedule; nonconvergence or an excessive
constraint residual produces an unsuccessful exit.

These are same-model numerical fixtures, not experimental calibration. Native
tests are authored, but were not executed in the delivery environment because
the Rust toolchain and remote build runners were absent.
