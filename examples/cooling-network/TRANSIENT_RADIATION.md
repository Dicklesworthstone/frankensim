# Radiative cooling during workload pulses

```bash
frankensim --json cooling-network \
  examples/cooling-network/radiative-contact-pulse.json

frankensim --json cooling-network-uq \
  examples/cooling-network/radiative-contact-pulse.json \
  examples/cooling-network/uq-radiative-pulse.json \
  --checkpoint radiative-pulse.uqcp
```

The checkpoint destination must be new. These are illustrative inputs, not
measured hardware data. The example retains the existing 20 W localized P1
source, two k(T) materials, heterogeneous constant heat capacities, matching
contact, and fan-driven bypass mixing. A 30-second pulse is followed by
120 seconds of cooling, using 2-second backward-Euler steps.

This extends the initial steady-only implementation described in
`RADIATIVE_COOLING.md`: the same root `radiation` declaration can now accompany
`transient`. The physical patch model is unchanged. Each named cooling surface
also radiates to a fixed, large isothermal black reservoir with view factor one.
Constant emissivity and surroundings temperature remain fixed throughout the
trajectory. The patch mean drives its secant Robin coefficient and fourth-power
heat law; this is not pointwise integration of T(x)^4 or enclosure radiosity.

## New-temperature radiation, not a frozen old-temperature correction

Every endpoint solves the existing discrete storage equation together with
nonlinear material conduction, convection and radiation. At fixed air references,
the radiation loop recomputes its coefficient from the new surface means and
calls the existing linear or nonlinear backward-Euler solver. Both the raw
radiation-temperature residual and independently recomputed nonlinear/applied
heat discrepancy must pass. The temperature and watt tolerances are explicit.

Every radiation and air iteration uses exactly the SAME previous accepted
physical temperature field. Iterating a boundary law never advances time or
adds another storage increment. Material k(T) still requires `transient.nonlinear`;
its Newton policy, capacity matrix, contact operator and energy checks are not
replaced. Constant-material requests keep their linear inner solver.

## Separate mechanisms and acceptance

The air callback receives only convective heat. The combined FEM Robin report
contains convection plus applied radiation, and its independently accumulated
region sums must agree with the whole-boundary total. Accepted endpoints also
recompute the nonlinear radiative heat from their actual surface means.

The accepted timestep obeys the original absolute energy gate on

```
delta_stored_j - dt * (source_w - air_heat_gain_w - radiative_heat_w)
```

Radiation is positive outward. Hot surroundings can produce negative radiation
loss and increase stored energy even with zero workload; neither sign is clipped.
Only accepted endpoints contribute to integrated heat. Adaptive full-step trials,
rejected half-step pairs and intermediate coupling/radiation fields do not enter
physical history. Their completed FEM evaluations do count as computational work.

The original wall deadline covers every nested solve. `radiation.max_iterations`
is the maximum number of inner FEM evaluations per air-reference evaluation.
Exhaustion refuses without publishing a partial trajectory. This is not a new
per-iteration time allowance or an enlarged linear/coupling budget.

## Output and repeated schedules

Every noninitial `transient.history` row adds `radiative_heat_w` when radiation
is enabled. `transient.radiative_energy_loss_j` is the backward-Euler sum over
accepted endpoints of that cycle. `air_energy_gain_j` remains air heat alone.
`energy_residual_j` includes both mechanisms. `total_solid_solves` and
`forward_solid_solves` count actual nested FEM evaluations, not just air sweeps.

The root `radiation` object describes the FINAL accepted endpoint: patch means,
emissivities, surroundings, secant coefficients, applied/recomputed powers and
maximum mismatch. Its powers are not integrated window energies. For repeated
runs the final `transient` object still describes the last cycle in local time;
`repeated_cycles.radiative_energy_loss_j` sums ALL completed cycles, and each
cycle summary retains its own radiative energy. The global energy gate includes
this sum and storage still telescopes from the original initial field.

Ordinary fixed and adaptive timesteps, fixed repeated cycles, existing forward
periodic/controller execution, and derivative-free transient sizing all use the
same endpoint producer. Radiative adjoints remain unavailable and explicitly
refuse, including gradient-guided sizing. Mesh studies and steady nested design
searches with radiation retain their prior refusals. No initial-state radiation
energy is invented before the first timestep.

## Uncertainty

Select `qoi.kind=transient-sampled-peak` in the UQ request. Existing radiation
emissivity and surroundings-temperature targets now vary complete trajectories,
including a declared fixed number of repeated cycles. One draw is held for the
whole trajectory, not redrawn at every timestep. Existing chunk/checkpoint and
family-confidence mechanisms remain unchanged. Interrupted trajectories provide
no observation; resume retries the same ordinal from its original initial state.
Invalid physical samples remain terminal failures, never filtered observations.

Peaks include the initial state and accepted endpoints only. This does not bound
an inter-step maximum, discretization error, radiation-model discrepancy or
physical reliability. The small sample example is an empirical propagation,
not evidence that sixteen trajectories resolve a stringent probability target.

## Focused checks

```bash
cargo test -p fs-cli --test cooling_transient_radiation
cargo test -p fs-cli --test cooling_radiation
cargo test -p fs-cli --test cooling_uq_radiation
```

Nine new actual-command Rust tests cover independently manufactured hot/cold
reservoir endpoints, nonlinear contact pulses, signed/window heat accounting,
repeated versus explicitly unrolled schedules, accepted adaptive half-steps and
rejected trials, default-path and patch-order replay, missing-policy/derivative/
budget refusals, and exact transient-UQ checkpoint replay. They were not executed
in the authoring environment, which had no Rust toolchain.

Independent NumPy P1 calculations compared twelve constant/nonlinear/contact
endpoints against a separate nested-iteration mathematical mirror. The maximum
temperature difference was 3.29e-11 K. Deliberately freezing radiation at the old
temperature changed an endpoint by up to 0.001506 K. Two analytically manufactured
uniform endpoints agreed within 1.60e-14 K; their adaptive controls required four
rejected trials each. These are reference calculations, not Rust execution.

For the full declared pulse, the independent reference peak is 306.135975 K at
30 seconds, versus 306.165103 K without radiation. Over 150 seconds, 600 J input
splits into approximately 365.260990 J stored, 35.939220 J air gain and
198.799790 J radiative loss. Two repeated cycles peak at 307.020212 K at
180 seconds and retain 402.705557 J of signed radiative loss over 300 seconds.
None of these values is an experimental validation or continuum error bound.
