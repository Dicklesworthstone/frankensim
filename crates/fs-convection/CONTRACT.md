# CONTRACT: fs-convection

> Status: ACTIVE for bead `frankensim-extreal-program-f85xj.5.2`.

## Purpose and layer

`fs-convection` is the L3 low-cost heat-transfer correlation rung. It maps
declared dimensionless regime points to Nusselt numbers and typed convective
coefficients without letting a formula escape its model card, source,
validity box, or model-form discrepancy.

It is intentionally separate from `fs-matdb`: the material database owns
immutable data and no executable closures. It is also separate from
`fs-conduction`: the conduction solver consumes a Robin coefficient but does
not choose or evaluate a convection model.

Runtime dependencies are `fs-conduction`, `fs-evidence`, `fs-math`, and
`fs-qty`. The dependency is one-way: the correlation rung lowers a typed,
evidence-paired boundary into conduction; conduction does not select a
correlation.

## Public types and semantics

- `CorrelationId` names eleven implemented relations: circular-duct 3.66 and
  4.36 limits; a thermally developing circular-duct relation; rectangular-
  duct constant-temperature and constant-flux limits; Dittus-Boelter;
  Gnielinski; laminar and turbulent average flat-plate relations;
  Churchill-Bernstein cylinder crossflow; and Churchill-Chu vertical-plate
  natural convection.
- `CorrelationCard` combines an `fs-evidence::ModelCard`, bibliographic
  provenance, and `DiscrepancyBasis`.
- `CorrelationInputs` makes Re, Pr, L/Dh, aspect ratio, Ra, and the
  heating/cooling convention explicit. Pe is deterministically derived from
  Re and Pr when both exist.
- `evaluate` returns `NusseltEvaluation` only after every shared
  `ValidityDomain` axis is present and inside its inclusive bounds.
- `NusseltEvaluation::heat_transfer_coefficient` returns
  `Evidence<HeatTransferCoefficient>` with `h = Nu k/L`. `k` is typed W/(m K),
  `L` is typed metres, and `h` is typed W/(m² K).
- `NusseltEvaluation::robin_boundary` returns a `CorrelationRobinBoundary`
  that owns both the evidence-bearing coefficient and the exact
  `fs-conduction::ThermalBc` row lowered from it. Private fields prevent the
  pair from drifting.

## Invariants

1. No unitless heat-transfer coefficient leaves the crate. The only public
   conversion returns `Evidence<HeatTransferCoefficient>`.
2. Missing, non-finite, non-positive, or out-of-domain dimensionless inputs
   refuse; there is no silent extrapolation.
3. Every successful value carries exactly the selected model card, its
   validity domain, assumptions, discrepancy, and deterministic evaluation
   provenance.
4. Card order and validity diagnostics are deterministic.
5. The 3.66/4.36 and Shah-London rectangular rows are analytic ideal limits.
   Their zero discrepancy denotes no empirical fit residual under the stated
   idealization; it is not a zero-error claim for hardware.
6. Other v1 discrepancy bands are conservative engineering allowances,
   explicitly labeled as such. They are not fabricated source-published
   confidence intervals and cannot earn a validation color.

## Error model

`CorrelationError` is total and teaching:

- `InvalidGroup` retains the axis and exact rejected bits;
- `OutOfDomain` lists every missing or violated axis with value and inclusive
  range;
- `InvalidDimensionalInput` rejects non-positive/non-finite `k` or `L`;
- `NonFiniteResult` refuses arithmetic overflow or a non-positive output.

No library path panics for caller input.

## Determinism class

Formula powers, logarithms, and square roots use `fs-math::det`; polynomial
evaluation has a fixed Horner tree. Evaluation provenance binds the stable
card id, direction convention, sorted group names, and exact float bits.

## Cancellation behavior

Each evaluation is bounded O(number of validity axes), currently at most four,
plus straight-line scalar arithmetic. No cancellation poll is required.

## Unsafe boundary

None. Workspace unsafe-code denial applies.

## Feature flags

None.

## Conformance tests

- catalog completeness, source presence, and shared-domain integrity;
- Level-A 3.66/4.36 limiting values and rectangular square/parallel-plate
  endpoints;
- per-formula frozen spot values;
- inclusive boundary acceptance plus missing/outside/non-finite refusals;
- G3 unit-rescaling invariance of the dimensioned `Nu k/L` conversion;
- evidence attachment and non-certification of empirical predictions;
- three-flow-rate heatsink slab integration through
  `fs-conduction::ThermalBc::robin`, including solve summaries and monotone
  heat removal.

## No-claim boundaries

- Source citations identify formula authority; the repository does not yet
  retain licensed copies of every source table or a cross-code validation
  dataset for these correlations.
- `EngineeringAllowance` is a declared design band, not a statistical
  confidence interval and not L4 experimental validation.
- The formula arithmetic is an `Estimate`, not an outward-rounded numerical
  enclosure.
- The library does not select among competing valid correlations, blend
  transition regimes, solve a boundary layer, compute pressure drop, or model
  fan operating points.
- Plate-fin behavior is represented only by smooth rectangular-channel
  limiting rows; interrupted-fin, louver, offset-strip, and full array effects
  remain outside v1.
- Cylinder crossflow is an isolated-pin baseline. Tube-bank and heatsink-array
  interference require separate cards and validation.
- The Robin integration test proves the coefficient-to-conduction seam on a
  small deterministic fixture. It is not conjugate CFD and does not promote
  the correlation prediction beyond its model evidence.
