# Finite-time cooling in the native project workflow

`cooling-transient.fsim` runs the seven-stage import, solve and retained-report
workflow with heat storage. It uses the existing `plate.stl` tetrahedron and
`aa6061.fsmcdpk` conductivity fixture. The solid starts at 293.15 K, receives
5 W for two seconds and exchanges heat through a fixed 10 W/(m² K) Robin
boundary at 293.15 K. The declared volumetric heat capacity is a **synthetic
100 J/(m³ K)** chosen to make this small example's time response visible; it
is not an AA6061 property claim. Each capacity retains its own source text.

From the repository root, with a built `frankensim` on `PATH`:

```bash
frankensim --json validate data/reference-project/cooling-transient.fsim
frankensim --json import data/reference-project/cooling-transient.fsim \
  data/reference-project/plate.stl transient.db --unit m --max-hole-edges 0
frankensim --json solve data/reference-project/cooling-transient.fsim transient.db \
  --materials data/reference-project/aa6061.fsmcdpk
```

Export the retained report or package using the run ID returned by `solve`:

```bash
frankensim --json report <run-id> transient.db
frankensim --json package <run-id> transient.db
```

The conduction stage uses backward Euler on four 0.5 s steps, then starts
again from the declared initial field on eight 0.25 s steps. Its twelve-step
cap includes both trajectories. The fine trajectory supplies the published
field and temperature maximum **at 2 s**. Every accepted step retains its
time, final-region maximum, stored energy change, net heat input, recomputed
linear residual and storage-minus-input energy residual. The energy gate is
1 µJ per step. The common endpoint energy report shows `storage_w` separately
from numerical `closure_w`.

The report's temporal estimate is `1.25 * abs(fine_final_k - coarse_final_k)`
for the actual selected final-region maximum, with assumed order one. This
is Estimated temporal error. It does not establish spatial convergence,
observed order, a maximum between time nodes, or a complete uncertainty budget.
If both grids agree exactly in floating point, the estimate stays unavailable
instead of claiming zero error. The temperature-limit margin is a final-time
nominal comparison; missing engineering error terms keep the decision
indeterminate.

For another native project, add `:transient (transient ...)` inside its
`conduction` declaration, providing an initial absolute temperature, horizon,
maximum coarse step, combined step cap, energy tolerance and one sourced
volumetric capacity per conduction region. All quantities use coherent SI.
For static power, the horizon and maximum step imply
`N = ceil(horizon / max-step)`; admission requires
`3*N <= max-steps <= 10000`. With scheduled power, the same limit applies to
the sum of the step counts between successive switches. Capacities and every
other input are part of canonical project and run identity. Existing v10
projects migrate with storage absent, so they retain their steady behavior;
v11 transient projects migrate with their existing static power unchanged.

This product slice supports fixed capacity, constant or temperature-dependent
conductivity from the bound material card, Dirichlet/Neumann/Robin boundaries
and matching finite contact. A temperature-dependent curve is evaluated at
each trial endpoint, including its actual conductivity derivative in the
Newton/FGMRES tangent; the solver never substitutes conductivity at the old
temperature. Each step starts from unchanged physical history until its
nonlinear residual and energy balance both pass. Initial
temperature may differ from prescribed boundary temperature; the first step
includes the discrete boundary storage reaction. An insulated heated body is
admitted because positive heat capacity anchors its finite-time equation.

For temperature-dependent conductivity, `solver.tolerance-rel` controls the
relative nonlinear residual against that endpoint's initial residual. The
absolute residual floor is 1% of the declared energy tolerance divided by the
square root of the vertex count. The existing numerical policy limits each
endpoint to 32 Newton updates and 24 backtracks per update; the same linear
iteration budget covers all of that endpoint's Newton corrections. These
controls, actual updates, backtracks and joule residuals are retained in the
`transient.nonlinear` object and each time-step row. An exhausted numerical
budget or material-validity boundary refuses the unfinished stage. Constant
conductivity retains the linear path and records `nonlinear: null`.

Radiation, natural or coupled airflow convection, spatial ladder/adaptive
studies and steady adjoint requests are not admitted
together with this declaration. Heat capacity remains temperature independent;
latent heat and phase changes require a different storage law. Stage
cancellation retains the preceding
completed pipeline stages; an unfinished conduction stage restarts from its
declared initial state.

## Pulse and duty-cycle heating

Schema v12 admits optional `:power-schedules` inside the transient declaration.
For example, this schedule deposits 5 J through a 20 W pulse, then switches
the solid off for the rest of the two-second window:

```lisp
:power-schedules
(power-schedules
  (schedule :region "solid" :source "declared 20 W startup pulse"
    :steps (steps
      (step :until 0.25s :watts 20.0kg·m^2·s^-3)
      (step :until 2.0s :watts 0.0kg·m^2·s^-3))))
```

Each `watts` value is the absolute delivered power from the preceding switch
(or time zero) until `until`. It replaces all static delivered power for that
volume region; the static row's duty factor is not applied a second time.
Unscheduled volume regions and surface heat inputs retain their static values.
Every schedule needs a source, strictly increasing positive times and an exact
final time equal to the transient horizon. Negative power and surface schedules
are refused.

Both time grids land on every switch, including short pulses between their
old uniform sampling points. For the example above with `max-step = 0.5 s`,
the two segments need `ceil(0.25/0.5) + ceil(1.75/0.5) = 5` coarse steps and
10 fine steps, so declare at least `max-steps = 15`. Actual step lengths and
source watts are retained. For schedules the uniform `coarse_step_s` and
`fine_step_s` fields are null; the `*_max_step_s` fields and each row's `dt_s`
describe the grid. The `workload` object includes the full sourced schedule
and integrated scheduled input in joules.

`cooling-pulsed.fsim` is the complete runnable project for this pulse. Use it
in the same validate/import/solve commands above with a fresh `pulsed.db`.
The reported temperature limit still applies at the final time, which can be
lower than the temperature reached during the pulse.
