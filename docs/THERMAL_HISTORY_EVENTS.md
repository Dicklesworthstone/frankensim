# Event-aligned thermal and phase histories

`advance_to` retains its original right-continuous endpoint sampling. Use
`advance_event_step` to avoid skipping a short heater pulse or radiative
reservoir change. Each call advances **one** interval, bounded by the earliest
source, boundary, Robin-reference or reservoir knot and an explicit maximum
timestep. It may stop before the requested target:

```rust
while run.time_s() < target_s {
    let step = run.advance_event_step(cx, target_s, maximum_step_s, config)?;
    // step.start_time_s and step.end_time_s bound the accepted interval.
    // step.report is the existing solver's complete physical response.
}
```

This method exists on `ThermalHistoryRun`, `EnthalpyHistoryRun`,
`RadiativeThermalHistoryRun` and `RadiativeEnthalpyHistoryRun`. For either
radiative wrapper, the loop condition is `run.conduction().time_s() < target_s`.
Temperature runs return `HistoryEventStep<ThermalHistoryStep>`; phase and
radiative runs retain their original full solver report types. The constructors,
material assignments, phase charts, capacities, interfaces and attached
radiation declarations are unchanged.

## Which side of a jump?

For a held table, a knot's new value acts **after** that knot. The event step
therefore uses the exact *left limit* at its end. A heater that is on from
0.125 to 0.25 seconds contributes its full 0.125-second pulse, rather than the
off value at 0.25 erasing the pulse. No epsilon is subtracted from physical
time, and even adjacent representable knots remain distinct. Continuous ramps,
linearly interpolated tables and Chebyshev histories keep ordinary endpoint
sampling. Their supplied knots and clamp endpoints also constrain stepping.

At a prescribed-temperature jump the reported field is the pre-jump state.
The new prescribed value enters the following finite step, with the existing
boundary-node capacity reaction. This does not model an instantaneous jump
impulse. Enthalpy histories still refuse temperature Dirichlet conditions;
a plateau temperature alone does not select a latent enthalpy.

## Work, cancellation and limits

Each successful interval consumes one existing lifetime step-budget slot.
Cancellation, exhausted budgets and rejected physical solves leave the previous
accepted state unchanged. Retain the returned reports, supply a fresh `Cx`
after cancellation, and explicitly extend the lifetime budget when necessary.
There is no hidden unbounded loop or unbounded event-list allocation. The maximum
step must be positive, finite and large enough to advance represented time.

Held source work is integrated over the correct intervals. This is **not** an
exact solution of temperature-dependent convection, conduction or radiation.
Those still use the existing backward-Euler endpoint equations and physical
energy gates. Extra event stops change the integration grid, so they can change
the trajectory relative to an unsplit run. No temporal error certificate,
substep maximum-temperature bound, melt-flow model or moving geometry is added.

Regressions cover short source pulses, independently assembled pre-jump
boundaries, nonlinear radiation, analytic latent plateau energy, distinct knots,
step caps, and cancellation/budget/future-reservoir refusals. The accompanying
independent equation checks are not a substitute for executing native tests.
