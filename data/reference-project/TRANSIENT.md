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
The horizon and maximum step imply `N = ceil(horizon / max-step)`; admission
requires `3*N <= max-steps <= 10000`. Capacities and every other input are part
of canonical project and run identity. Existing v10 projects migrate with
storage absent, so they retain their steady behavior.

This product slice supports fixed scalar conductivity, fixed capacity,
Dirichlet/Neumann/Robin boundaries and matching finite contact. Initial
temperature may differ from prescribed boundary temperature; the first step
includes the discrete boundary storage reaction. An insulated heated body is
admitted because positive heat capacity anchors its finite-time equation.
Time-dependent workloads, radiation, natural or coupled airflow convection,
spatial ladder/adaptive studies and steady adjoint requests are not admitted
together with this declaration. Stage cancellation retains the preceding
completed pipeline stages; an unfinished conduction stage restarts from its
declared initial state.
