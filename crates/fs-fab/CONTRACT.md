# CONTRACT: fs-fab

Manufacturing, fabrication & code compliance as its own layer: optimization
without fabrication semantics produces fantasy artifacts.

## Purpose and layer

Layer L4 (constraints/optimization). Depends only on `fs-evidence` (the `Color`
for Evidence-typed cost/carbon). Pure, deterministic.

## Public types and semantics

- `FabConstraint { name, kind, differentiability, cert, sense, limit, units }` —
  a one-sided scalar constraint. `margin(value)` (`>= 0` ⇒ satisfied),
  `satisfied`, `margin_gradient` (`±1` differentiable / `None` discrete),
  `repair(value)` (a `Repair` to the limit when violated).
- `ConstraintKind` — `Fabrication(Process)` or `Code(standard)`;
  `Differentiability` (Differentiable / Subdifferentiable / Discrete);
  `CertAvailability`; `Sense` (AtLeast / AtMost).
- Catalog constructors: `overhang_angle`, `min_feature_size`,
  `cnc_tool_radius`, `draft_angle`, `bolt_spacing_aisc` (AISC-360 §J3.3 = 2⅔d),
  `member_length_transport`, `rebar_spacing_aci` (ACI-318 §25.2.1 = max(25, d)).
- `evaluate(&FabConstraint, value) -> ConstraintResult`; `check_all(&[(FabConstraint,
  f64)]) -> FabReport { results, feasible, violations }` (detection +
  localization + repair).
- `Estimate { mean, rel_std }` with `std` + `color` (estimated); `process_cost`,
  `embodied_carbon` — Evidence-typed modeled quantities.

## Invariants

- A design with ANY violated constraint is `feasible = false` — no fantasy
  artifact passes.
- Repair suggestions move the feature exactly to the limit; satisfied
  constraints yield no repair.
- Differentiable constraints expose a `±1` margin gradient matching finite
  differences; discrete ones expose `None`.
- The AISC/ACI catalog limits encode the published rules.
- Cost/carbon estimates are `Color::Estimated` (modeled quantities).

## Error model

Total functions; no panics.

## Determinism class

Fully deterministic: margins, reports, and estimates are pure functions of the
inputs.

## Cancellation behavior

None (synchronous pure functions).

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/fab.rs` (7 cases): constraint families detect violations + margins;
AISC/ACI rules encode published limits; violations carry repair suggestions;
differentiable constraints pass a gradient gate (discrete has none); the report
localizes violations + prevents fantasy artifacts; cost/carbon are
Evidence-typed estimates; determinism.

## No-claim boundaries

- v1 models each constraint as a ONE-SIDED SCALAR feature check with metadata
  (differentiability / cert / kind) + a linear repair; the full geometric
  evaluators (overhang via surface-normal fields, tool reachability via
  clearance/visibility, drainability/no-enclosed-voids via topo-persistence)
  consume geom-queries and are staged.
- Sheet-metal, composite-layup, casting parting-surface, weld access,
  assembly-sequence precedence, and full ACI development-length / cover
  families are staged with their interfaces defined here.
- The catalog encodes representative AISC-360 / ACI-318 limits; a complete,
  edition-pinned rule set is a data-expansion.
- Cost/carbon unit rates are SUPPLIED parameters; sourcing them from a costing
  database is a downstream integration.
