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
conductivity from the bound material card, Dirichlet/Neumann/Robin boundaries,
matching finite contact, declared ambient radiation, natural convection and
quasi-steady airflow from the native fan network. Natural convection and
airflow are separate endpoint modes. A temperature-dependent curve is evaluated at
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
solid response to 32 Newton updates and 24 backtracks per update; the same
linear iteration budget covers all Newton corrections within that response.
Without convection or radiation coupling there is one response per endpoint.
Each airflow or natural-convection evaluation can repeat the solid response
through radiation's separate outer-iteration cap. The receipt discloses the
complete product of these
bounded allowances for one physical endpoint. These controls, actual updates,
backtracks and joule residuals are retained in the `transient.nonlinear`
object and each time-step row. An exhausted numerical budget or
material-validity boundary refuses the unfinished stage. Constant conductivity
uses the linear solid response and records `nonlinear: null`.

Spatial ladder/adaptive studies, steady adjoint requests and simultaneous
natural and airflow-convection boundaries are not admitted together with this
declaration. Airflow is quasi-steady; its own heat storage and time-varying fan
dynamics are not modeled. Heat capacity remains temperature independent;
latent heat and phase changes use the explicit total-enthalpy storage below. Stage
cancellation retains the preceding
completed pipeline stages; an unfinished conduction stage restarts from its
declared initial state.

## Total enthalpy and latent heat

`cooling-enthalpy.fsim` exercises the spatial enthalpy solver on the same
tetrahedron, with an insulated boundary and 5 W of internal heating for two
seconds. It supplies a **synthetic** 1 kg/m³ reference density and a declared
enthalpy chart with a 293.35 K latent plateau from 20 to 80 J/kg. These are
educational model inputs, not properties of the AA6061 conductivity fixture.
The source text records that distinction in the project itself.

```bash
frankensim --json validate data/reference-project/cooling-enthalpy.fsim
frankensim --json import data/reference-project/cooling-enthalpy.fsim \
  data/reference-project/plate.stl enthalpy.db --unit m --max-hole-edges 0
frankensim --json solve data/reference-project/cooling-enthalpy.fsim enthalpy.db \
  --materials data/reference-project/aa6061.fsmcdpk
```

The uniform initial enthalpy is zero at 293.15 K. The tetrahedron's volume is
1/6 m³, so the declared density gives a reference mass of 1/6 kg. The
independent energy balance predicts a 10 J increase after two seconds,
specific enthalpy 60 J/kg, temperature 293.35 K and liquid mass fraction
2/3. The spatial solver computes that uniform limit and retains the nodal
enthalpy and phase fields plus accepted-step energy checks. Temperature
agreement between the two grids on a plateau is not proof of phase or
enthalpy convergence.

Schema v13 adds `:enthalpy` inside `transient`, alongside an explicitly empty
`:capacities (capacities)`. Every region needs a source, attributed card
identity, reference density, initial specific enthalpy and an ordered chart.
Specific enthalpy includes latent heat once; a separate capacity cannot be
active with this storage model. Temperature alone cannot select the initial
phase fraction on a plateau, which is why initial enthalpy is mandatory.
The chart supports all-solid, all-liquid or full solid-to-liquid data without
extrapolation. Newton iterations and line-search halvings are explicit bounded
controls. Existing v12 projects migrate with their capacities unchanged.

This native mode supports fixed reference geometry, temperature-dependent
conductivity, regional power schedules, matching finite-resistance contact,
prescribed Neumann/Robin boundaries and sourced ambient radiation. Prescribed
temperatures, natural/airflow convection and adjoint requests currently refuse
for enthalpy storage. The declared curves and card identities create no
material-query receipts, expansion, phase kinetics or experimental validation.

### Latent storage with ambient radiation

`cooling-enthalpy-radiation.fsim` combines the same explicit synthetic phase
chart with convection, the sourced 0.85 gray-emissivity fixture, and a 20 W
startup pulse lasting 0.25 s. The pulse delivers 5 J before the power switches
off. Both grids land on that switch and follow the remaining cooling window;
the final stored energy accounts for the heat lost to the exterior.

```bash
frankensim --json validate data/reference-project/cooling-enthalpy-radiation.fsim
frankensim --json import data/reference-project/cooling-enthalpy-radiation.fsim \
  data/reference-project/plate.stl enthalpy-radiation.db --unit m --max-hole-edges 0
frankensim --json solve data/reference-project/cooling-enthalpy-radiation.fsim \
  enthalpy-radiation.db --materials data/reference-project/aa6061.fsmcdpk \
  --materials data/reference-project/gray-surface.fsmcdpk
```

Every radiation trial solves from the same accepted old enthalpy. A step must
pass both the complete nonlinear residual and the physical storage balance
using actual radiation, not just a small change in a trial coefficient. The
enthalpy trajectory's common `net_input_w`, `energy_residual_j` and cumulative
energy fields use the physical balance, including signed nonlinear radiation.
Each step's `radiation` object retains the applied and nonlinear heat, original
convection, physical residual and work. The normal final radiation report keeps
the exact emissivity card evidence and convection/radiation decomposition.

`transient.radiation` counts all outer trials, Newton updates and Krylov
iterations against their complete multiplied budgets. Detailed Newton/linear
histories describe the final inner solve of each accepted step, which is the
history exposed by the shared owner. The reservoir model uses patch-mean
temperature; it does not infer enclosure exchange or pointwise T(x)^4
integration. The independent regular-tetrahedron regressions cover both
radiative heat loss and heating from a hotter reservoir, plus latent-plateau
energy and phase accounting.

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

## Ambient radiation during a pulse

`cooling-radiative-pulse.fsim` adds the existing sourced 0.85 gray-emissivity
surface card to the complete pulse example. The surface exchanges radiation
with a 293.15 K reservoir alongside its prescribed convection. Supply both
immutable material packs:

```bash
frankensim --json validate data/reference-project/cooling-radiative-pulse.fsim
frankensim --json import data/reference-project/cooling-radiative-pulse.fsim \
  data/reference-project/plate.stl radiative-pulse.db --unit m --max-hole-edges 0
frankensim --json solve data/reference-project/cooling-radiative-pulse.fsim \
  radiative-pulse.db --materials data/reference-project/aa6061.fsmcdpk \
  --materials data/reference-project/gray-surface.fsmcdpk
```

The declaration uses the existing `:radiation (radiation ...)` block inside
`conduction`; no new emissivity or radiation schema is needed. The shared
gray-patch law uses the actual area-mean endpoint temperature and a fixed
reservoir. The same accepted old temperature is held fixed through every
radiative trial and every inner conductivity correction. Only a complete
endpoint can advance physical time.

A small radiation temperature change alone cannot accept a step. The patch
heat mismatch, actual implicit residual and physical energy balance must all
pass their declared gates. The inner solid residual tolerance is tightened
to reserve room for the radiative coupling error. Actual surface temperatures
must remain within the emissivity card's validity domain, including when
the reservoir supplies heat to the solid.

The `radiation` object in each time-step row retains applied radiation watts,
nonlinear radiation watts, convection watts, physical joule residual and its
threshold, physical energy residual, and the prescribed-temperature reaction
with boundary storage. Existing step `net_input_w` and `energy_residual_j`
retain the frozen-secant operator's balance; the explicitly named physical
fields check the full nonlinear endpoint. Their distinction prevents a
frozen radiative coefficient from being reported as an exact nonlinear law.

`transient.radiation` reports coarse, fine and total radiative trial counts,
the maximum physical energy residual and the complete per-endpoint Krylov
allowance. A solid response shares one linear budget across all of its Newton
corrections; at most `radiation.max-iterations` responses occur per outer
airflow or natural-convection evaluation (or per endpoint without either).
The original complete `conduction.radiation` receipt describes the final fine
endpoint and retains the immutable surface-card evidence. Numerical
cancellation remains cooperative; wall-time checks occur between numerical
operations and do not promise an intra-kernel deadline.

This example and its capacity are synthetic. The result is an Estimated
final-time temperature with a coarse/fine temporal comparison, without a
continuous-time peak, enclosure-radiation, phase-change or transient-adjoint
claim.

## Fan-cooled startup pulses

`examples/heatsink-fan/heatsink-fan-pulse.fsim` adds a native transient
declaration to the existing finned-heatsink and fan-network project. The
20 W pulse lasts 0.25 s and is followed by 1.75 s with regional power off.
Both grids land on the switch: five coarse steps and ten fine steps share
the declared fifteen-step cap. The project declares a **synthetic
100000 J/(m³ K)** solid capacity, with an explicit source saying that this is
not an AA6061 property claim. The existing fan curve is illustrative as well.

From the repository root:

```bash
frankensim --json validate examples/heatsink-fan/heatsink-fan-pulse.fsim
frankensim --json import examples/heatsink-fan/heatsink-fan-pulse.fsim \
  examples/heatsink-fan/heatsink.stl fan-pulse.db --unit m --max-hole-edges 0
frankensim --json solve examples/heatsink-fan/heatsink-fan-pulse.fsim \
  fan-pulse.db --materials data/reference-project/aa6061.fsmcdpk
```

The existing `airflow-convection` laws supply each branch, ordered segment,
inlet temperature and convection correlation. The flow-network stage supplies
the operating point; retained exterior areas and the correlation card supply
the conductance. Air warms along each declared path at every thermal endpoint.
Independent branches retain their own inlet and mass flow while sharing one
solid field. The operating point stays fixed through the trajectory: only the
solid stores heat, and neither air storage nor fan acceleration is inferred.

Every air evaluation and inner radiation/Newton solve sees the same accepted
old solid temperature. The shared IQN-ILS solver has at most 100 airflow
evaluations per endpoint. It uses the existing branch-local heat gates and
tightens its raw reference-temperature criterion to at most
`0.01 * energy_tolerance_j / (dt_s * sum(hA))`. A final independent joule gate
uses actual heat carried away by every air branch, off-path convection,
nonlinear ambient radiation, prescribed boundary reaction and solid storage.
Radiation exchanges heat with its own reservoir; it is excluded from the
air's heat gain.

Each step's `conjugate` object retains the standard exchange receipt,
`applied_references_k`, all inner work counts and the full coupled energy
residual. The published solid is exactly the last counted response at those
applied references. The air marcher also reports its updated references;
their difference is bounded by the retained `reference_tolerance_k`.
No additional solve is performed after convergence.

`transient.conjugate` aggregates both grids' airflow evaluations and solid
responses, discloses the full endpoint Krylov allowance, and records the fine
trajectory's integrated coupled heat balance. `conduction.conjugate` describes
the final fine endpoint. Existing applied-operator energy and radiation
physical-energy fields keep their meanings; `coupled_energy_residual_j`
separately checks actual air transport. Final-time temperature limits and
the coarse/fine temporal estimate retain the same scope as the other native
transient examples.

## Passive cooling after a startup pulse

`cooling-passive-pulse.fsim` replaces the radiative pulse example's fixed
convection coefficient with the existing Churchill–Chu vertical-plate law.
It declares a 0.1 m characteristic length and 293.15 K ambient, retaining
the gray-emissivity surface and its separate 293.15 K radiation reservoir.
The geometry, heat capacity and vertical-plate idealization are synthetic;
this example is not experimental validation of passive cooling.

```bash
frankensim --json validate data/reference-project/cooling-passive-pulse.fsim
frankensim --json import data/reference-project/cooling-passive-pulse.fsim \
  data/reference-project/plate.stl passive-pulse.db --unit m --max-hole-edges 0
frankensim --json solve data/reference-project/cooling-passive-pulse.fsim \
  passive-pulse.db --materials data/reference-project/aa6061.fsmcdpk \
  --materials data/reference-project/gray-surface.fsmcdpk
```

Use the existing boundary declaration:

```lisp
(natural-convection :target "solid"
  :characteristic-length 0.1m
  :ambient-temperature 293.15K
  :correlation "convection.churchill-chu-vertical-plate")
```

At every backward-Euler endpoint, the area-mean wall temperature determines
the actual Rayleigh number and convection coefficient. Heated and cooled
walls retain their signed heat exchange. Natural convection composes with
scheduled power, prescribed-temperature storage reactions, matching contact,
endpoint `k(T)` and ambient radiation. All coefficient, radiation and Newton
trials hold the same accepted old field fixed.

The existing initial coefficient uses a wall guess 10 K above ambient.
That is only a starting coefficient, so a powered solid may start exactly at
ambient. The solved candidate must satisfy the actual law. Exact equilibrium
and actual Rayleigh numbers outside the card's `[0.1, 1e12]` domain retain
their named refusals; no near-equilibrium extension or clipped coefficient is
introduced. An unpowered equilibrium or a cooling tail that leaves the
card's domain therefore refuses the unfinished conduction stage.

A small coefficient update cannot accept an endpoint by itself. Publication
requires all three gates: relative coefficient change at most `1e-10`, a
full physical free-residual norm below its fixed target, and the original
declared joule energy tolerance. The residual target is fixed from the first
fully evaluated candidate, using the requested relative tolerance plus a
one-percent energy allowance divided by the square root of the free-degree
count. Subsequent inner solves tighten their allowances to reserve room
for the outer law discrepancy. Total energy balance alone cannot conceal
opposing local heat-equation errors.

The natural loop permits at most 80 coefficient evaluations per endpoint.
The admitted work product includes every inner radiation response and every
Newton/Krylov iteration; direct physical endpoint audits have a separate
160-evaluation cap. A successful endpoint uses `2 * iterations - 1` direct
audits. An exhausted limit or invalid candidate publishes no partial field.

Each step's `natural` object retains applied and physical coefficients,
the physical free-residual norm and target, physical convection/radiation,
the prescribed-temperature reaction, the physical energy residual, the
existing natural-convection card receipt and all nested work counts.
The published field and common energy report belong to the last counted
solid response at the retained applied coefficients. The physical-law
audit is reported separately; convergence does not trigger an extra solve.

`transient.natural` aggregates both grids' coefficient iterations, solid
responses and direct endpoint audits, together with the maximum accepted
physical energy residual and the fine trajectory's integrated heat balance.
`conduction.natural` retains the final fine endpoint's physical card receipt.
The fifteen-step pulse grid and final-time temperature/QoI scope are the
same as in the other pulse examples.
