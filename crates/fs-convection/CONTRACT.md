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

- `CorrelationId` names twelve implemented relations: circular-duct 3.66 and
  4.36 limits; a thermally developing circular-duct relation; rectangular-
  duct constant-temperature and constant-flux limits; one simultaneously
  developing constant-temperature rectangular-duct table slice at
  `alpha*=0.5`, `Pr=0.72`; Dittus-Boelter; Gnielinski; laminar and
  leading-edge-corrected mixed average flat-plate relations;
  Churchill-Bernstein cylinder crossflow; and Churchill-Chu vertical-plate
  natural convection.
- `CorrelationCard` combines an `fs-evidence::ModelCard`, bibliographic
  provenance, and `DiscrepancyBasis`.
- `CorrelationInputs` makes Re, Pr, L/Dh, aspect ratio, Ra, and the
  heating/cooling convention explicit. Pe is deterministically derived from
  Re and Pr when both exist; Gz is deterministically derived as
  `Re Pr/(L/Dh)` when all three inputs exist.
- `evaluate` returns `NusseltEvaluation` only after every shared
  `ValidityDomain` axis is present and inside its inclusive bounds.
- `NusseltEvaluation::heat_transfer_coefficient` returns
  `Evidence<HeatTransferCoefficient>` with `h = Nu k/L`. `k` is typed W/(m K),
  `L` is typed metres, and `h` is typed W/(m² K).
- `NusseltEvaluation::robin_boundary` returns a `CorrelationRobinBoundary`
  that owns both the evidence-bearing coefficient and the exact
  `fs-conduction::ThermalBc` row lowered from it. Private fields prevent the
  held pair from drifting, but the row accessor exposes the downstream raw
  coherent-SI representation for interoperability.

## Invariants

1. The standalone coefficient conversion returns
   `Evidence<HeatTransferCoefficient>`. Robin lowering is an explicit
   interoperability boundary: its accessor exposes the downstream
   `ThermalBc` raw coherent-SI row, which has no standalone `fs-convection`
   authority if cloned or detached from `CorrelationRobinBoundary`.
2. Missing, non-finite, non-positive, or out-of-domain dimensionless inputs
   refuse; there is no silent extrapolation.
3. Every successful value carries exactly the selected model card, its
   validity domain, assumptions, discrepancy, and deterministic evaluation
   provenance.
4. Card order and validity diagnostics are deterministic.
5. The 3.66/4.36 and fully developed Shah-London rectangular rows are
   analytic ideal limits. Their zero discrepancy denotes no empirical fit
   residual under the stated idealization; it is not a zero-error claim for
   hardware.
6. The developing rectangular card executes only the Chapter VII, Table 52
   `alpha*=0.5`, `Pr=0.72` row for `Gz` in `[10, 220]`. Its `Gz` interval
   `[0.01, 10]` is a declared linear bridge between the analytic CWT limit at
   zero and the first source row; it is not represented as a source-published
   curve.
7. Other v1 discrepancy bands, including that developing table card, are
   conservative engineering allowances, explicitly labeled as such. They are
   not fabricated source-published confidence intervals and cannot earn a
   validation color.

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
plus straight-line scalar arithmetic or one bounded 14-knot lookup. No
cancellation poll is required.

## Unsafe boundary

None. Workspace unsafe-code denial applies.

## Feature flags

None.

## Conformance tests

- catalog completeness, source presence, and shared-domain integrity;
- direct `fs-vvreg` execution bindings for both Level-A circular-duct limits,
  with a complete two-row family partition, JSON comparison verdicts, and
  frozen rectangular-duct square values;
- per-formula frozen spot values;
- independent test-owned transcriptions of four Shah-London Table 52 points,
  including one interpolation midpoint, plus exact source metadata checks;
- narrow developing-rectangular refusal checks for `Gz`, `Pr`, `Re`, and
  aspect ratio;
- limiting-behavior checks for the nonconstant cards, each asserting a
  relation derivable from the cards' own published coefficients rather than
  from this implementation's output:
  - the developing-flow Hausen card converges to the fully-developed
    constant-wall-temperature limit as the Graetz number vanishes, evaluated
    against `CircularDuctLaminarCwt` at points admitted by both cards
    (`L_over_Dh` in `[50, 1000]`), with a strictly positive gap, the analytic
    bracket `gap <= 0.0668 Gz`, and first-order halving of the residual;
  - the developing rectangular card approaches the matching aspect-0.5 CWT
    limit through the declared lower-`Gz` bridge, and a fixed-geometry
    Reynolds perturbation changes Nu and evaluation provenance;
  - the two flat-plate cards are continuous across `Re = 5e5`, which is
    simultaneously the inclusive upper bound of the laminar card and the
    inclusive lower bound of the mixed card; the relative disagreement is
    below 0.1% and is exactly Prandtl-independent, and the raw difference at
    `Pr = 1` pins the rounding residue of the published `871` constant;
  - the Dittus-Boelter direction exponents degenerate bitwise at `Pr = 1`
    while provenance still records the declared direction;
  - the Gnielinski correction denominator collapses exactly at `Pr = 1`,
    reducing the relation to its Prandtl-free form;
  - the Churchill-Chu `Ra -> 0` intercept, which is outside the admitted
    domain and unreachable through `evaluate`, is recovered by affine
    extrapolation in `Ra^(1/6)` from two in-domain points;
- inclusive boundary acceptance plus missing/outside/non-finite refusals;
- G3 unit-rescaling invariance of the dimensioned `Nu k/L` conversion;
- evidence attachment and non-certification of empirical predictions;
- three-flow-rate heatsink slab integration through
  `fs-conduction::ThermalBc::robin`, including solve summaries and monotone
  heat removal.

## No-claim boundaries

- Source citations identify formula authority; the repository does not retain
  a licensed copy of Shah-London or a cross-code validation dataset. The
  developing rectangular card manually encodes only its declared Table 52
  numeric slice.
- The two Level-A duct limits are resolved from `fs-vvreg` while the tests
  execute, but no comparison receipt or machine fingerprint is persisted into
  the corpus. This is a test-time binding, not a registry authority promotion.
- Those two rows compare a card that returns its literal unconditionally
  against an equal registry literal, so they exercise no arithmetic. They
  establish that the domain gate admits the point and that the two crates
  agree on the constant; they are not evidence of limiting behavior. The
  limiting-behavior authority for `3.66` is the Hausen convergence check,
  which is the only test that fails when the developing-flow relation drifts
  away from that constant.
- Except for the separately labeled Table 52 spot checks, the per-formula
  frozen spot values are outputs of this implementation, not source-published
  table entries. They detect drift; they do not validate any card against its
  source. No broader external or published-table oracle exists.
- The limiting-behavior checks verify internal consistency between formulas
  and their own published coefficients. Agreement between two cards at a
  shared domain edge constrains the constants that relate them; it is not
  experimental validation and does not upgrade either card's evidence colour.
- `EngineeringAllowance` is a declared design band, not a statistical
  confidence interval and not L4 experimental validation.
- The formula arithmetic is an `Estimate`, not an outward-rounded numerical
  enclosure.
- The library does not select among competing valid correlations, blend
  transition regimes, solve a boundary layer, compute pressure drop, or model
  fan operating points.
- Plate-fin behavior is represented only by smooth rectangular-channel
  limiting rows and one narrow simultaneously developing source-table slice.
  Interrupted-fin, louver, offset-strip, shroud, bypass, and full-array effects
  remain outside v1.
- The fully developed rectangular-duct cards admit aspect ratios from 0.001
  through 1.0; the developing table card admits exactly 0.5. The exact
  parallel-plate limit at aspect ratio zero is outside those domains and is
  not executed; current tests freeze the admitted square value at 1.0.
  Unlike the Churchill-Chu intercept, this endpoint is not recovered by
  extrapolation: the shape function is a quintic in the aspect ratio, so no
  two-point affine fit reaches it, and no such test is claimed here.
- The Churchill-Bernstein `Re -> 0` conduction limit of `0.3` is neither
  admitted nor approachable: the `Pe >= 0.2` floor forces `Pr >= 0.2` at
  `Re = 1`, where the card still returns roughly twice that value. No test
  constrains it.
- Cylinder crossflow is an isolated-pin baseline. Tube-bank and heatsink-array
  interference require separate cards and validation.
- The Robin integration test proves the coefficient-to-conduction seam on a
  small deterministic fixture. It is not conjugate CFD and does not promote
  the correlation prediction beyond its model evidence. A downstream
  `ThermalBc` row cloned through the public accessor is detachable raw
  coherent-SI data and does not independently carry that evidence.
