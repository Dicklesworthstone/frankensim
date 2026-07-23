# CONTRACT: fs-airflow

> Status: ACTIVE for bead `frankensim-extreal-program-f85xj.5.5`.

## Purpose and layer

`fs-airflow` is the L3 low-cost enclosure-airflow rung and the dependency-safe
consumer seam for the thermal vertical's decision-facing QoIs. It retains typed fan
pressure/volume-flow data, validates monotone interpolation, applies bounded
fan-law speed scaling, composes quadratic loss elements in series and parallel,
requires an explicit leakage branch, solves the fan/system operating point, and
combines that record with a steady `fs-conduction` solution without making either
lower-layer solver depend on the other.

Runtime dependencies are `fs-blake3`, `fs-conduction`, `fs-convection`,
`fs-evidence`, `fs-ivl`, `fs-math`, `fs-qty`, and `fs-regime`. Results flow
outward as evidence-bearing typed quantities and as a Re/Pr handoff to the
existing convection rung.

## Public types and semantics

- `FanCurve` owns typed `(VolumetricFlowRate, Pressure)` points, source identity,
  tolerance authority, a stall boundary, and the validity range for RPM affinity
  scaling.
- `FanBank` represents identical fans in series or parallel. At speed ratio
  `s`, flow scales as `s` and pressure as `s^2`; series pressure adds and
  parallel flow adds.
- `LossElement` uses a typed `LossResistance` in Pa/(m^3/s)^2. `LossNetwork`
  composes quadratic elements recursively in series and parallel.
- `EnclosureNetwork` cannot be constructed without a distinct
  `LeakageElement`; leakage is not an implicit constant.
- `solve_operating_point` returns an `OperatingPoint` with an interval-Newton
  unique-root bracket for the nominal model, weaker physical uncertainty
  estimates, per-terminal flows, and nominal leakage fraction.
- `OperatingPoint::correlation_handoff` converts one branch flow to typed mean
  velocity, computes Reynolds through role-tagged `fs-regime` dimensions, and
  produces `fs-convection::CorrelationInputs` without discarding evidence.
- `qoi::extract_thermal_qois` consumes a `ConductionSolution`, its exact
  `ConductionMesh`, an `OperatingPoint`, declared junction/surface regions, a
  cited fan-efficiency interval, and a cited maximum-temperature requirement.
  It emits the five E05.10 families: deterministic maximum junction
  temperature, pressure drop, fan input power, surface mean/uniformity, and
  thermal margin. An absent requirement refuses; there is no default limit.
- Every emitted `ThermalQoi<T>` carries both the existing `Evidence<T>` view and
  an `EngineeringUncertaintyBudget` with exactly one term for roundoff,
  solver/algebraic, discretization, geometry, parameters, boundary conditions,
  model form, and measurement. A missing propagation theorem/receipt is a
  named `Unknown`, never `Negligible` or zero.
- `ThermalQoiSet::audit_operating_envelope` is the mandatory final E05.10
  product gate. It requires exactly one consumed-card declaration for each of
  the seven emitted records, derives each incoming color from the actual
  `Evidence` receipt, runs one all-card/all-point `fs-regime` audit, applies
  each exact receipt to the matching eight-term budget, and returns the audit
  beside the updated values. Missing, duplicate, or foreign QoI declarations
  refuse; callers cannot supply the pre-audit color. Fully in-domain admission
  leaves the QoI set byte-for-byte unchanged, while any partial/out-of-domain
  envelope makes the affected model-form term explicitly Unknown under the
  exact receipt identity. Overrides remain acknowledgements only.
- `FanCurve::model_card` exposes the exact fan-law card already consumed by
  airflow evidence so product admission can audit the retained base-curve flow
  and speed-ratio axes without synthesizing a parallel card definition.
- The E05.10 battery also consumes the public Dittus-Boelter correlation card
  alongside that fan card. Its adversarial variant changes only Reynolds
  number, proving the correlation receipt is the sole demotion cause while the
  actual fan card remains in-domain.
- `ThermalQoiSet::audit_operating_envelope_with_cards` accepts owner-neutral
  `RegimeAuditCard` identity/validity projections for matdb and other native
  card schemas. The existing evidence-`ModelCard` method is an exact wrapper;
  parity is executable, so widening the card source cannot change receipts or
  budget demotions.
- The operating-point pressure/flow envelope populates the conditional
  boundary-condition term for pressure and power. The cited total-efficiency
  interval populates the fan-power parameter term. Both remain accompanied by
  an `Unknown` model-form term, so finite conditional bands cannot upgrade the
  synthetic/quadratic airflow model to validated product authority.
- Surface mean uses exact P1 triangle integration (`area * vertex mean`).
  Spread is the selected surface-vertex range. The reported standard deviation
  is the area-weighted dispersion of face means; it is explicitly not an exact
  integral of the pointwise squared P1 field.

## Invariants

1. Fan flow points increase strictly and pressure never increases.
2. Fan-law scaling refuses outside its caller-declared speed-ratio validity
   range.
3. The interval below the explicit stall boundary is non-admissible; the solver
   refuses an intersection there.
4. Every loss resistance is finite, positive, typed, and carries source and
   uncertainty authority.
5. A complete operating-point result has exactly one `Certified` interval root
   and no `Possible` root boxes for the nominal declared model.
6. Manufacturer/loss/leakage uncertainty is attached as model-form `Estimate`
   evidence. It is never relabeled as a rigorous physical enclosure.
7. Terminal branch order follows deterministic network traversal, and all
   provenance hashes are order-stable and bind the complete fan curve, source,
   tolerance, fan-bank configuration, recursive network topology, loss data,
   and explicit leakage identity.
8. Junction and surface declarations canonicalize index order and reject
   duplicates. Equal junction maxima choose the lowest canonical vertex index.
9. Every thermal QoI budget has exactly eight terms. Widening a valid upstream
   pressure/flow or efficiency interval cannot shrink the corresponding
   conditional term; changing a requirement or efficiency authority rebinds
   the affected QoI identity even when the nominal scalar is unchanged. The
   temperature-QoI identity binds canonical tetrahedral connectivity and every
   physical vertex coordinate, so geometry-only changes also rebind it.
10. Raw temperature extrema, surface summaries, and margin remain
    `NumericalKind::NoClaim` until an admitted DWR/refinement-to-QoI map exists.
    A conduction residual measured in watts is not converted into kelvin by
    dimensional wishful thinking.

## Error model

`AirflowError` refuses malformed curves, invalid tolerance or speed domains,
empty network groups, zero/invalid resistances, stall operation, absent curve
intersections, incomplete/ambiguous root searches, unknown branches, and bad
convection-handoff inputs. Caller input does not panic.

## Determinism class

Interpolation and network traversal have fixed order. Square roots use
`fs-math::det`; the numerical root uses deterministic `fs-ivl` subdivision.
Results are intended to be bit-stable on the same ISA. Cross-ISA G5 evidence is
not yet retained.

## Cancellation behavior

Curve evaluation and each finite network reduction are bounded scalar work.
The interval search has an explicit 65,536-box ceiling and returns a structured
refusal rather than running without bound. No asupersync cancellation poll is
required for this rung.

## Unsafe boundary

None. Workspace unsafe-code denial applies.

## Feature flags

None.

## Conformance tests

- G0 fan interpolation, monotonicity refusal, and series/parallel resistance
  composition;
- G0 identical-fan series/parallel affinity identities;
- explicit stall refusal and a sign-changing, unique interval root bracket;
- three declared fan speeds and leakage-resistance sensitivity;
- typed branch velocity/Reynolds handoff into `fs-convection`;
- semantic-identity separation when uncertainty authority changes without
  changing the nominal operating point.
- G0 E05.10 fixture emitting all five QoI families and seven records (the
  uniformity family has mean, spread, and face-mean standard deviation), each
  with a complete eight-term budget and term-by-term provenance rendering;
- deterministic region-order/tie-break equivalence, missing-requirement and
  malformed-region refusals, G3 upstream-envelope widening monotonicity, and
  source-only identity rebinding for fan power and margin.

## No-claim boundaries

- Fan points and tolerances are caller-supplied. Synthetic fixtures prove API
  and algebra behavior, not a manufacturer product's performance.
- AMCA/ASHRAE fan laws justify similarity scaling, not arbitrary operation:
  the curve's declared speed-ratio domain remains binding.
- Piecewise-linear interpolation and quadratic loss coefficients do not model
  swirl, recirculation, acoustic interaction, thermal buoyancy, compressibility,
  fouling, or installation system effects.
- The nominal root bracket certifies only the declared mathematical model.
  Tolerance-propagated flow, pressure, branch splits, and Reynolds values remain
  `Estimate`, not validated hardware envelopes.
- Parallel fans are identical and equally loaded. Unequal curves, unstable
  parallel operation, active control, fan-fan interference, and transient
  startup remain outside this slice.
- No retained manufacturer table, wind-tunnel corpus, CFD comparison, or
  experimental enclosure validation exists; there is no L4 or product claim.
- The thermal QoI consumer does not close E05.10's external validation, DWR,
  mesh-refinement, sensor, or naked-scalar lint obligations. It emits the rich
  budget and now enforces final validity intersection when the orchestrator
  supplies the complete consumed-card registry, card-use map, and operating
  envelope. Completeness/authenticity of those supplied authorities remains an
  orchestration and package/checker responsibility; the broader E05.10 bead
  remains open.
- The conditional pressure/flow and efficiency intervals are only as sound as
  their caller-declared source envelopes and the stated quadratic model. The
  always-explicit unknown model-form term prevents their interpretation as a
  whole-product uncertainty bound.
- Fan power means input power under the declared total efficiency,
  `Delta p * Q / eta_total`. Motor/controller transients, reactive power,
  acoustic power, and installation effects are outside this slice.
