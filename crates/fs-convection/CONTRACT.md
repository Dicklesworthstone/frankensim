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

- `CorrelationId` names twenty-one implemented relations: circular-duct 3.66 and
  4.36 limits; a thermally developing circular-duct relation; rectangular-
  duct constant-temperature and constant-flux limits; one simultaneously
  developing constant-temperature rectangular-duct table slice at
  `alpha*=0.5`, `Pr=0.72`; Dittus-Boelter; Gnielinski; laminar and
  leading-edge-corrected turbulent average flat-plate relations;
  Churchill-Bernstein cylinder crossflow; and Churchill-Chu vertical-plate
  natural convection; Churchill-Chu horizontal-cylinder natural convection;
  and the stable, unstable-laminar and unstable-turbulent horizontal-plate
  rows of ASHRAE 2025 Fundamentals, Chapter 4, Table 9 (T9.7--T9.9); and
  Manning--Qureshi aiding laminar mixed convection over an isothermal vertical
  plate, including its pure forced and natural boundary-layer limits; and
  Churchill isolated-sphere natural convection, including its conductive limit;
  and Whitaker forced convection over an isolated sphere with explicit
  bulk/surface viscosity correction; and Hollands horizontal enclosed air-layer
  transfer, including gas conduction and stable/subcritical limits; and Hollands
  inclined air-gap transfer with the explicit source inclination correction.
- `CorrelationCard` combines an `fs-evidence::ModelCard`, bibliographic
  provenance, and `DiscrepancyBasis`.
- `CorrelationInputs` makes Re, Pr, L/Dh, aspect ratio, Ra, bulk/surface
  viscosity ratio, inclination in degrees, and the
  heating/cooling convention explicit. Pe is deterministically derived from
  Re and Pr when both exist; Gz is deterministically derived as
  `Re Pr/(L/Dh)` when all three inputs exist. `mixed(Re, Ra, Pr)` supplies both
  flow drivers; their magnitudes presume the caller has established aiding
  vertical alignment. The dimensional owner checks the actual signs.
  `horizontal_air_layer(signed_Ra, Pr, gap_over_min_span)` retains the enclosure's
  stability sign and its whole-plate aspect ratio instead of an external-flow
  length or a mesh-cell width. `inclined_air_layer(Ra, Pr, gap_over_min_span,
  inclination_degrees)` retains full-gravity Ra and the physical inclination
  separately; the source formula applies its own gravity projection.
- `evaluate` returns `NusseltEvaluation` only after every shared
  `ValidityDomain` axis is present and inside its inclusive bounds.
- `NusseltEvaluation::heat_transfer_coefficient` returns
  `Evidence<HeatTransferCoefficient>` with `h = Nu k/L`. `k` is typed W/(m K),
  `L` is typed metres, and `h` is typed W/(m² K).
- `EvaluationSourceRegion` makes the executed authority region machine
  readable. Closed-form rows report `CitedFormula`; the developing rectangular
  card reports `TableSliceInterpolation` at and above its first retained knot
  and `DeclaredEngineeringBridge` below it.
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
2. Missing, non-finite, incorrectly signed, or out-of-domain dimensionless
   inputs refuse; there is no silent extrapolation. The aiding mixed card
   admits zero Re or Ra (and corresponding derived Pe/Gz); at least one flow
   driver must be positive. The isolated sphere admits Ra=0 as the physical
   steady conduction limit Nu=2. The horizontal enclosed air layer retains
   finite signed Ra, with stable and zero-driving total conduction Nu=1.
   The inclined air layer admits nonnegative Ra and angles in [0, 60] degrees,
   retaining Nu=1 below the projected onset. Other supplied groups remain positive.
3. Every successful value carries exactly the selected model card, its
   validity domain, assumptions, discrepancy, and deterministic evaluation
   provenance.
4. Card order and validity diagnostics are deterministic.
5. The 3.66/4.36 and fully developed Shah-London rectangular rows are
   analytic ideal limits. Their zero discrepancy denotes no empirical fit
   residual under the stated idealization; it is not a zero-error claim for
   hardware.
6. For `Gz` in `[10, 220]`, the developing rectangular card linearly
   interpolates in `Gz` between fourteen retained Chapter VII, Table 52
   `alpha*=0.5`, `Pr=0.72` knots. It does not interpolate in `log(Gz)` or the
   source's `x*` coordinate. Its `Gz` interval `[0.01, 10)` is a declared
   linear-in-`Gz` engineering bridge between the analytic CWT limit at zero
   and the first retained table knot; it is not represented as a
   source-published curve. `EvaluationSourceRegion` distinguishes those paths,
   with `Gz=10` belonging to the table slice. At `Gz=15`, the linear-in-`Gz`
   chord is about 0.9% below a power-law interpolation through the two
   bracketing knots; that low bias is an interpolation-choice observation, not
   a rigorous error bound. The piecewise-linear interpolant is only C0 and its
   derivative jumps at knots.
7. Other v1 discrepancy bands, including that developing table card, are
   conservative engineering allowances, explicitly labeled as such. They are
   not fabricated source-published confidence intervals and cannot earn a
   validation color.

## Horizontal enclosed air layer

`HollandsHorizontalAirLayer` implements the horizontal reduction of
Hollands, Raithby and Unny, report COO-2597-2 (1976), §4 Eq. (4-1),
[doi:10.2172/7125623](https://doi.org/10.2172/7125623), reproducing the original
[1976 paper](https://doi.org/10.1115/1.3450517). For unstable layers it uses
`Nu = 1 + 1.44 max(0, 1-1708/Ra) + max(0, (Ra/5830)^(1/3)-1)`.
Stable horizontal stratification, zero driving and subcritical positive Ra
retain total `Nu=1`; gas conduction must not be added again. Radiation is separate.

The signed group uses `g beta (T_lower-T_upper) gap³/(nu alpha)`, with gas
properties at the mean of the two actual isothermal wall temperatures.
The unstable branch admits Ra through `1e5`; finite negative Ra describes the
analytic stable conduction limit. Applying the report's thin-layer assumption
to **both** spans requires `gap/minimum_span < 0.1`. The numerical Pr interval
`[0.70, 0.72]` is a conservative implementation policy around the source's air
value near `0.71`, **not a source-published validation interval**. Pr alone does
not establish fluid identity. Density/property variation, quasi-steady gas
storage neglect, and whole-plate geometry remain the dimensional owner's
explicit physical admissions. Inclination, ventilation and wall shear are absent.

Analytic signed-Ra partials are zero in the stable and subcritical branches,
including Ra=0. The Ra derivative is unavailable at the source's kinks `1708`
and `5830`, or at the unstable domain endpoint. Unrelated interior partials
remain available. G0/G3 tests in `tests/enclosed_air_layer.rs` compare independent
high-precision arithmetic, stability, source limits, continuity and typed scaling.

## Inclined enclosed air layer

`HollandsInclinedAirLayer` implements the inclined equation reproduced in the
[EnergyPlus 8.3 Engineering Reference, Heat Transfer Coefficients](https://bigladdersoftware.com/epx/docs/8-3/engineering-reference/solar-collectors.html#heat-transfer-coefficients),
citing [Hollands et al. (1976)](https://doi.org/10.1115/1.3450517). With inclination
`theta` measured from horizontal, `R = Ra cos(theta)`, and `x+ = max(0, x)`, it uses
`Nu = 1 + 1.44 (1 - 1708 sin(1.8 theta)^1.6/R) (1 - 1708/R)+
+ ((R/5830)^(1/3) - 1)+`. The first factor is not clipped; it is nonnegative
whenever its multiplied onset is active in this admitted domain. At zero or
subcritical projected Ra the implementation returns exactly one before dividing.
Replacing only Ra with its normal gravity component would omit a source term.

This card requires nonnegative full-gravity Ra through `1e5` and angle in
`[0, 60]` **degrees**. The 60-degree ceiling is an explicit conservative
implementation subset of the original measurements. Air identity, mean-wall
properties, full normal gap, `gap/minimum_span < 0.1`, Pr admission `[0.70, 0.72]`,
and quasi-steady storage assumptions follow the horizontal card's policy. The
coefficient includes total gas conduction plus convection; radiation remains
separate. Heating from above at nonzero inclination refuses: taking an absolute
Ra or swapping wall temperatures would falsely change the physical orientation.
Vertical cavities, non-air fluids, ventilation and wall shear remain unsupported.

Zero inclination shares the original horizontal arithmetic bit for bit; the
signed horizontal card retains its existing source identity and stable branch.
Analytic partials use the full-gravity Ra and angle **per degree**, including both
source inclination terms. Ra and angle partials are unavailable at projected
branch switches `R=1708` and `R=5830`; validity endpoints omit only the partial
that crosses them. Unrelated admitted partials stay available. The numerical
leaf and `tests/inclined_air_layer.rs` compare independent high-precision values
and derivatives, source/domain refusals, horizontal compatibility and typed
conductance scaling. These are implementation checks, not experimental validation.

## Error model

`CorrelationError` is total and teaching:

- `InvalidGroup` retains the axis and exact rejected bits;
- `OutOfDomain` lists every missing or violated axis with value and inclusive
  range;
- `InvalidDimensionalInput` rejects non-positive/non-finite `k` or `L`;
- `NonFiniteResult` refuses arithmetic overflow or a non-positive output.

No library path panics for caller input.

## Determinism class

Formula powers, logarithms, trigonometry and square roots use `fs-math::det`; polynomial
evaluation has a fixed Horner tree. Evaluation provenance binds the stable
card id, direction convention, sorted group names, and exact float bits. The
source-region marker is deterministically derived from the already-bound card
id and exact `Gz` bits, so adding it does not perturb unrelated cards'
provenance identities.

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
- G0/G3 aiding mixed-plate checks of independent 60-decimal-digit source-formula
  arithmetic, strict laminar bounds, continuous pure limits, power-sum
  homogeneity down to tiny drivers, and analytic derivatives including the
  Prandtl-dependent mixing exponent and a representable derivative whose
  intermediate power weight underflows;
- G0/G3 forced-sphere checks of independent 70-decimal-digit source-formula
  arithmetic, explicit viscosity correction, original data-domain refusals,
  analytic input partials and dimensional sensitivity propagation;
- direct `fs-vvreg` execution bindings for both Level-A circular-duct limits,
  with a complete two-row family partition, JSON comparison verdicts, and
  frozen rectangular-duct square values;
- per-formula frozen spot values;
- independent test-owned transcriptions of three of the fourteen retained
  Shah-London Table 52 knots (`Gz=10,20,100`), one explicitly
  implementation-owned interpolation midpoint, plus exact source metadata and
  bridge-disclosure checks;
- machine-readable source-region checks on both sides of `Gz=10`, including
  the exact boundary and the analytic companion card;
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

- `WhitakerSphere` implements the average isolated-isothermal-sphere relation
  in [Whitaker (1972), Eq. (33), Table 1](https://doi.org/10.1002/aic.690180219).
  It admits the source data ranges inclusively: `Re in [3.5,76000]`,
  `Pr in [0.71,380]`, and `mu_bulk_over_surface in [1,3.2]`. The caller must
  supply the viscosity ratio explicitly. All other transport properties,
  including conductivity in `h=Nu k/D`, use the bulk/free-stream temperature;
  the wall temperature supplies only the surface viscosity. Film-property
  substitution, silently assuming equal viscosities, and clamping the ratio
  change the source model and are not allowed. Heating a gas with a hotter
  sphere commonly produces a ratio below one and therefore refuses this card.
  The source's steady conductive term Nu=2 is retained, but zero flow remains
  outside the fitted data domain and is not extrapolated. No buoyancy blend,
  evaporation, rarefaction, confinement or particle-array effects are supplied.
  The 30% engineering allowance reflects the paper's maximum observed scatter;
  it is neither a statistical confidence interval nor independent validation.
  Analytic Re, Pr and viscosity-ratio partials are available only when the
  corresponding input does not move a validity endpoint.

- `ManningQureshiAidingVerticalPlate` implements the **average**, aiding-only
  isothermal vertical-plate correlation in Manning and Qureshi, Journal of Heat
  Transfer 116(3), 777--780 (1994), doi:10.1115/1.2910940, Eqs. (4), (8), (10),
  and (11). Its constant-property Boussinesq boundary-layer assumptions remain
  visible on the card. The source's strict laminar bounds `Re < 5e5` and
  `Ra < 1e9` are encoded as inclusive bounds at the immediately preceding
  floating-point numbers. `Pr in [0.003,100]` is a declared implementation
  policy spanning source comparison points, not a continuously validated
  experimental range. An adapter's film averaging is a declared frozen-property
  approximation; the source does not mandate that sampling temperature. No
  opposing-flow difference formula, transition blend, finite-width correction,
  stagnant conductive floor, or experimental confidence interval is supplied.
  The pure limits remain boundary-layer models and do not resolve
  diffusion-dominated or creeping-flow heat transfer. At Re=Ra=0 the card
  refuses; thermal owners can handle zero temperature difference as equilibrium.
  The formula cancels `Pr^(1/4) (Ra/Pr)^(1/4)` and normalizes component powers
  before summation, avoiding a division by Re or a singular mixed parameter
  `Gr/Re^2`. Input partials include the derivative of `n=3.5 Pr^0.075` and omit
  only directions that move a validity endpoint, including each zero driver.
  Formula tests establish arithmetic and limit behavior, not external validation.

- `NaturalConvectionGeometry` selects vertical plates, horizontal cylinders,
  isolated spheres, or the stable/unstable horizontal face from signed
  `g beta (Tw-Ta)`. Characteristic lengths are height for the vertical plate,
  diameter for the cylinder and sphere, and full planform area/perimeter for
  the horizontal plate. The horizontal plate uses film properties and retains the
  handbook's Rayleigh ranges. Its source does not specify a universal Prandtl
  interval; omission of that axis is not all-fluid experimental validation.
  The cylinder's Prandtl box is explicitly an implementation policy. Numeric
  range endpoints are admitted inclusively. The unstable plate selects the
  laminar card at Ra=8e6 and turbulent above, retaining the approximately 4.46%
  empirical jump. Each fixed card has local input partials; no derivative or
  continuity is claimed across selection. Signed-expansion orientation follows
  the Boussinesq force direction, not experimental validation of anomalous
  liquids. No low-Rayleigh extrapolation, enclosure, mixed-convection blend,
  or confidence interval is invented. The new tests check exact power-law
  values, orientation/sign metamorphisms, source domains and derivatives, and
  the cylinder against the handbook's approximate 733 W heater calculation
  within a declared 3% calculation tolerance; this is not measured validation.

- `ChurchillSphere` supplies the **average** coefficient for an isolated,
  isothermal sphere in an unbounded quiescent single-phase fluid, using its
  diameter in both Ra and `h=Nu k/D`. The formula and domain follow
  [COMSOL 6.4 Heat Transfer Theory, External Natural Convection, Sphere,
  Eq. (4-194)](https://doc.comsol.com/6.4/doc/com.comsol.help.heat/heat_ug_theory.07.100.html):
  `Nu=2+0.589 Ra^(1/4)/(1+(0.469/Pr)^(9/16))^(4/9)`, `Ra<=1e11`, `Pr>=0.7`.
  The same relation is attributed to Churchill (1990) by ASHRAE 2025
  Fundamentals, Chapter 4, Table 9 (T9.11). No upper Prandtl bound is published;
  the finite-number domain uses `f64::MAX`, without claiming validation in
  every fluid. Transport properties are sampled at the mean wall/ambient
  temperature; an adapter's beta-based Rayleigh number is its declared
  Boussinesq approximation, not COMSOL's finite-density-difference definition.
  Unlike the plate and cylinder cards, the sphere admits Ra=0: Nu=2 gives the
  exact steady spherical conduction limit `G=4 pi k r` to an unbounded fluid.
  This conductive component is already included and must not be added again
  as a separate parallel path. The model supplies no transient fluid diffusion,
  local angular coefficient, confinement, particle-array correction, or mixed
  forced/natural blend. The 20% discrepancy is an engineering allowance, not
  source-published statistical validation. `tests/natural_sphere.rs` checks
  independent high-precision arithmetic, diameter-based conductive scaling,
  signed-buoyancy symmetry, range admission and analytic input partials. At
  Ra=0 the singular Ra partial is absent; at positive tiny Ra its derivative
  is retained even when the buoyant increment rounds away in Nu.

- Source citations identify formula authority; the repository does not retain
  a licensed copy of Shah-London or a cross-code validation dataset. The
  developing rectangular card manually encodes only its declared Table 52
  numeric slice.
- Only 3 of the developing card's 14 retained numeric knots are independently
  repeated in the test suite. The other eleven values and the regular `Gz`
  knot grid are implementation data without an in-repository source copy.
  The grid may be a re-gridding of source `x*` values; no direct-transcription
  claim is made for its abscissae or the eleven unchecked ordinates.
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
- Apart from the explicit geometry/stability selector above, the library does
  not select among competing valid correlations, blend
  transition regimes, solve a boundary layer, compute pressure drop, or model
  fan operating points.
- Plate-fin behavior is represented only by smooth rectangular-channel
  limiting rows, one narrow simultaneously developing table slice, and its
  separately marked lower-`Gz` engineering bridge.
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


## Compartment-fire correlation rung and public comparison data

`fire` reuses the existing typed quantities, source provenance, and validity
owner for Heskestad flame height, virtual origin, plume entrainment and
centerline estimates; MQH hot-layer temperature; Thomas/Babrauskas flashover
heat release; and Kawagoe ventilation capacity. Heskestad buoyancy uses actual
supplied ambient density, heat capacity, temperature and gravity. The flaming
entrainment branch is explicit and joins the above-flame relation. Exactly
zero fire power is dormant; invalid empirical regimes are refused. Total and
convective heat release are distinct. MQH implements its published kW-based
coefficient with explicit SI conversions and enclosing area excluding openings.
Kawagoe requires explicit heat release per mass of air and efficiency and does
not replace a species/oxygen balance or a pressure-resolved vent calculation.

Every empirical result retains its source card, qualitative applicability,
and unquantified model discrepancy. No universal finite validity box,
zero uncertainty, calibration claim, or fire-safety certificate is invented.
Six focused tests cover hand-calculated plume/density behavior, branch joining,
dormant/refusal cases, flashover/ventilation units, metadata, and all eight
published NIST SP 1169 Table 2.1 MQH examples (0.006 K tolerance).

`fire::validation` computes allocation-free log-ratio sample statistics over
positive model/observation pairs. Unknown experimental standard uncertainty
leaves separated model scatter and adjusted FDS bias absent. With explicitly
supported experimental relative sigma it applies the pinned FDS guide's
small-uncertainty approximation and refuses separation when its scatter
constraint fails. Three tests independently check these formulas/refusals.
Temperature comparisons use positive rises over a stated ambient baseline.

`data/reference/compartment_fire` preserves six Steckler 1982 cases with
pinned upstream commits and per-file Git blob hashes, original profile and
reduced-layer data, original CFAST input decks, and their intact NIST notice.
Authored CFAST wall properties are distinguished from the experimental
assembly; reduced-layer measurement standard uncertainty remains unknown.
Dalmarnock is excluded. The runnable `compartment_fire_correlations` example
requires both the retained cases CSV and an explicit effective wall coefficient.
Its h=10 W/(m² K) illustration produces geometric upper-layer-rise ratio
1.314403236 and combined log scatter 0.086888962; this is an uncalibrated
comparison, not a validation pass. No parameter is fitted to the observations.

## Finite developing circular-duct sections

`developing::evaluate_circular_duct_section` integrates the existing Hausen
mean relation over explicit heated-entrance distances in diameters. It retains
the original admitted upstream/downstream `NusseltEvaluation`s and returns a
distinct `CircularDuctSection`; the ordinary correlation catalog is unchanged.
The original whole-length card and section evaluator share one formula owner.

With `s=x/D` and `P=Re Pr`, cumulative transfer is
`s Nu_mean(s) = 3.66 s + 0.0668 P s^(2/3)/(s^(2/3)+0.04 P^(2/3))`.
The section mean is its difference divided by section length. A rationalized
cube-root quotient avoids subtracting nearly equal cumulative values for short
sections. At the exact heated entrance, cumulative transfer has the analytic
zero limit; the first section reproduces the original downstream mean. Every
positive upstream point and downstream point must satisfy the unchanged card's
Re, Pr and L/D support. No small-distance floor or domain extrapolation exists.

The result's positive dimensional coefficient uses `h=Nu_section k/D` and
retains its original mean-card receipts. It is derived held-property transfer,
not a new fitted correlation or a validated local discrepancy bound. The mean
card's engineering allowance does not automatically bound a difference of
cumulative predictions. No new sensitivity or numerical enclosure is claimed.
Subdivision additivity assumes the same fluid properties, diameter, velocity
and wall state throughout the compared intervals. Local source/temperature
updates remain a separate held-property approximation.

Three regression sources cover cumulative additivity versus entrance resets,
the analytic differential limit for adjacent floating-point lengths, and
original domain/dimensional refusal. Source/API and scalar algebra review were
completed; Rust parsing, compilation and native numerical execution were
unavailable in the delivery environment. These tests are not reported executed.

## Finite external flat-plate strips

`developing::evaluate_flat_plate_section` evaluates the heat-transfer coefficient
over physical leading-edge positions `[a,b]`, given `Re_b`, `Pr` and `a/b`.
It retains the same original mean card at both endpoints and uses
`Nu_strip = h b/k = (Nu_mean(b)-Nu_mean(a))/(1-a/b)`. The characteristic length
for dimensional conversion remains **b**, not the strip width. Thus splitting
a physical plate preserves total heat without starting a new boundary layer
at every surface cell. The original mean and section paths share one arithmetic
owner; existing mean-card arithmetic is unchanged.

Laminar strips use `Nu_mean(b)/(1+sqrt(a/b))`. Turbulent strips cancel the common
871 entrance term analytically and evaluate the power difference with `expm1`
and the deterministic logarithm, preserving the local limit even for adjacent
representable positions. A strip beginning at the exact leading edge retains
the original whole-plate result, coefficient evidence and Robin row exactly.
Every positive upstream point must satisfy the SAME original card as the
downstream point. The turbulent whole-plate mean already includes its original
laminar entrance; a later strip must lie entirely in the supported turbulent
range. No cross-card transition blend, upstream clamp or zero-distance floor is
introduced, and the rounded cards' transition mismatch is not hidden.

Positive-strip coefficients retain both original source evaluations, the held
fluid and isothermal-wall assumptions, and a provenance binding to the physical
fraction, conductivity and downstream distance. Their model discrepancy is
explicitly unbounded: the mean-card engineering allowance is not a validated
local or section allowance. Numerical evidence remains an estimate, with no
new enclosure, derivative or physical validation claim. The source-bearing
coefficient lowers directly into the existing private `CorrelationRobinBoundary`
pair. Heat additivity assumes the same properties, speed and wall state across
compared intervals; varying local film or wall temperatures remains a separate
held-property approximation, not a resolved boundary-layer solution.

Four native scalar tests cover independent analytic means and local limits,
unequal physical subdivisions, streamwise decay and cancellation in thin strips.
`tests/flat_plate_sections.rs` covers original endpoint/domain refusal,
dimensional heat additivity, exact first-strip compatibility, source/evidence
retention and Robin lowering. Full crate integration execution still requires
the repository's constellation dependencies.
