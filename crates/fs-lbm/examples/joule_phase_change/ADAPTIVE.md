# Material-sensitive adaptive heating and melting

Use the existing material files and optional electrical-interface file with:

```sh
cargo run -p fs-lbm --features implicit-enthalpy-coupling \
  --example joule_phase_change -- \
  crates/fs-lbm/examples/joule_phase_change/transport_phase_dependent.txt \
  10 60 0.5 crates/fs-lbm/examples/joule_phase_change/contact_low.txt \
  --adaptive 1e-10 0.03 1e-7
```

The final values explicitly specify the numerical Jacobian stencil in **cell
joules**, the full run's **per-cell temporal-indicator budget in joules**, and
the **minimum retry span in seconds**. The earlier `0.5` is the maximum
comparison span, not a promised accepted timestep. Change only `contact_low.txt`
to `contact_high.txt` to change the interface material data. The data format and
physical setup are the same as in explicit and fixed implicit modes. These
fixtures are synthetic, not measurements of named metals or alloys.

The 60 observation intervals divide the per-cell budget equally. Within each
interval, the library compares one full backward-Euler step with two half steps
and measures `mass * abs(h_two_halves - h_full)`. It apportions that interval's
allowance by comparison duration and accumulates the accepted indicators. This
includes latent storage even when both trajectories have the same temperature.
The accepted state is the conservative two-half-step solution; it is **not** a
Richardson-extrapolated state. Indicator sums are estimates, not certified global
error bounds. Refine timestep, stencil and spatial resolution independently.

Rejected comparisons and failed nonlinear trials spend work but contribute no
physical heat. Only accepted fine steps count toward stored heat and source work.
Current-solve work is bounded over the entire observation interval, including
coarse solves, fine solves, numerical derivative probes and discarded attempts.
The original limit is 20,001 current solves per observation interval; the
adaptive mode also limits it to 20,000 implicit attempts. Minimum span and either
work limit are refusals, not permission to relax physical or accuracy checks.

Only nonlinear nonconvergence and excessive temporal indicators trigger retries.
Unsupported material data, numerical stencils outside the chart, cancellation
and physical energy failures are returned unchanged or with explicit numerical
context. There is no material extrapolation or fallback to the explicit solver.
The nonlinear Jacobian/error adapter is shared with the finite-circuit coupling
path through `fs-solver::finite_difference`; neither consumer owns another
Jacobian implementation.

The power pulse ends exactly at 30 seconds in this invocation. Every observation
interval lies on one side of shutoff. Stored enthalpy and phase continue evolving
after power is removed. Library callers must likewise split at drive changes and
scheduled events. Target-time events remain pending until processed explicitly.

CSV retains all original columns and appends the implicit diagnostics plus
`max_interval_cell_temporal_indicator_J` and `interval_rejections`. In adaptive
mode, `newton_iterations` accumulates **accepted fine-trajectory** iterations;
`current_solves` includes **all** trial work. `max_interval_cell_residual_J` is
the largest physical cell defect on accepted halves in the latest interval,
not the defect of just the last half. No incomplete interval emits an output row,
but earlier accepted stdout rows may exist if a later interval refuses.

Without a mode flag, the original explicit path is used. `--implicit STENCIL_J`
retains fixed-step closed feedback. Both implicit modes require the feature;
a build without it refuses rather than silently selecting another algorithm.
This remains stationary Ohmic transport and fixed-reference-mass thermal phase
evolution, not moving melt, induction or a mechanical-history model.
