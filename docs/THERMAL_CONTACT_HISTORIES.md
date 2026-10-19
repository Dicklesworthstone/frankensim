# Material contact feedback in scheduled phase-change runs

`bc::history::march::enthalpy::contact::ContactEnthalpyHistoryRun` connects an
unadvanced `EnthalpyHistoryRun` to the complete `TemperatureContacts` set. Both
must use the identical `ConductionMesh`. Constant interfaces belong in that
contact set alongside the temperature-dependent ones; the base history's
`ThermalHistoryProblem::interfaces` must be `None`. A second set is refused,
never silently ignored or double-counted.

Build the phase-history run with its original material assignments, specific
enthalpies, source and boundary schedules. Construct `TemperatureContact` laws
from ordered interface cards, explicit query coordinates and face pairs, then
bind the complete set with `TemperatureContacts::new`. Attach it with
`ContactEnthalpyHistoryRun::new(cx, run, contacts, coupling)`. Optional sourced
radiation is attached once with `with_radiation(cx, patches, config)` before
physical stepping. The existing `scheduled_melting` and
`temperature_contact_phase` examples show the underlying phase/material inputs.

Three stepping methods retain the same contact laws and optional radiation:

- `advance_to` uses the explicit right-continuous endpoint.
- `advance_event_step` stops at the next load, boundary or reservoir knot and
  preserves held pulses with exact interval-end left limits.
- `advance_adaptive_step` compares a full interval with two half intervals.
  Temperature, specific enthalpy and liquid fraction must all pass their local
  discrepancy gates. It accepts the actual half-step trajectory without
  extrapolating fields or phase.

Every original contact iteration queries resistance at the candidate interface
mean temperature. Each trial starts from the same input enthalpy; source queries
and radiation iterations do not advance physical time. The accepted endpoint
is re-evaluated with its own contact resistance and radiation before its local
residual and total energy are admitted. Internal contact exchanges are not
external heat sources. `ContactEnthalpyStep` preserves applied and physical
contact states, selected-source receipts, cumulative work and the complete
physical energy audit. Use its `physical` report for coupled heat accounting.

The wrapper exposes only immutable state/contact access. Errors and cancellation
publish no partial state; the adaptive pair is atomic and consumes two physical
step slots. `extend_step_budget` preserves materials, contacts and radiation.
Endpoint work is bounded by the original contact/radiation/Newton policies;
adaptive `endpoint_solves` counts top-level coupled endpoint calls, each subject
to those nested bounds. It is not a count of every internal Krylov iteration.

Geometry, reference mass, contact pressure, source selection and nonthermal
query coordinates remain fixed. This does not add opening/separation mechanics,
contact-area evolution, melt flow, enclosure exchange or a global time-error
certificate. The new driver reuses the original physical solvers and their
refusals; it does not relax their gates when a temporal trial fails.
