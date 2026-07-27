# CONTRACT: fs-airflow

> Status: ACTIVE for beads `frankensim-extreal-program-f85xj.5.5` and
> `frankensim-extreal-program-f85xj.5.7` (conjugate coupling).

## Purpose and layer

`fs-airflow` is the L3 low-cost enclosure-airflow rung and the dependency-safe
consumer seam for the thermal vertical's decision-facing QoIs. It retains typed fan
pressure/volume-flow data, validates monotone interpolation, applies bounded
fan-law speed scaling, composes quadratic loss elements in series and parallel,
requires an explicit leakage branch, solves the fan/system operating point, and
combines that record with a steady `fs-conduction` solution without making either
lower-layer solver depend on the other.

The `conjugate` module closes that combination into a two-way exchange. It owns
the stream-wise air state the flow-network types do not carry, and the
partitioned fixed point between a solid conduction solve and that air path.
`fs-couple` has no coupling-loop driver — its only iteration entry points are
hardwired to a private added-mass fixture map, and its contract assigns the
driver to the consumer — so this crate is that consumer. Interface quasi-Newton
(IQN-ILS) does not exist in `fs-couple` and is not called.

Runtime dependencies are `fs-blake3`, `fs-conduction`, `fs-convection`,
`fs-couple`, `fs-evidence`, `fs-exec`, `fs-ivl`, `fs-ladder`, `fs-math`,
`fs-qty`, and `fs-regime`. Results flow outward as evidence-bearing typed
quantities, as a Re/Pr handoff to the existing convection rung, and as the
CHT ladder's correlation-rung transfer.

## Public types and semantics

- `FanCurve` owns typed `(VolumetricFlowRate, Pressure)` points, source identity,
  tolerance authority, a stall boundary, and the validity range for RPM affinity
  scaling.
- `FanBank` represents identical fans in series or parallel. At speed ratio
  `s`, flow scales as `s` and pressure as `s^2`; series pressure adds and
  parallel flow adds.
- `fan_law_source` exposes the scaling authority as structured
  `SourceProvenance`: AMCA Publication 201-02 (R2011), section 6.5.1, rather
  than the AMCA 210 / ASHRAE 51 laboratory test-method standard.
- `LossElement` uses a typed `LossResistance` in Pa/(m^3/s)^2. `LossNetwork`
  composes quadratic elements recursively in series and parallel. A loss
  element remains cardless unless its owner explicitly attaches a finite,
  nonempty `ValidityDomain`; `with_regime_validity` then binds the card name to
  `airflow.loss.<element>` and the card version to the retained source
  identifier. `LossNetwork::regime_audit_cards` and the enclosure wrapper
  expose those exact cards without silently deduplicating ambiguous names.
- `sharp_edged_orifice_loss` lowers a declared open area and typed air
  density to a `LossElement` through the sharp-edged thin-orifice model
  `R = ρ / (2 Cd² A²)`, with `Cd` the retained plateau constant
  `SHARP_EDGED_ORIFICE_CD` (0.61) cited through
  `ORIFICE_CD_SOURCE_CITATION`/`ORIFICE_CD_SOURCE_IDENTIFIER`. The declared
  `ORIFICE_RESISTANCE_UNCERTAINTY_REL` (0.15, engineering allowance) covers
  the sourced `Cd ∈ [0.60, 0.65]` large-opening spread, and every produced
  element carries a `loss_reynolds` validity card over
  `ORIFICE_PLATEAU_REYNOLDS` so the plateau assumption is audited after the
  solve. Non-finite or non-positive area or density refuses as
  `InvalidOrificeInput`.
- `EnclosureNetwork` cannot be constructed without a distinct
  `LeakageElement`; leakage is not an implicit constant.
- `solve_operating_point` returns an `OperatingPoint` with an interval-Newton
  unique-root bracket for the nominal model, weaker physical uncertainty
  estimates, per-terminal flows, and nominal leakage fraction. Flow bounds use
  the low-fan/high-resistance and high-fan/low-resistance corners. Pressure
  bounds separately solve the low-fan/low-resistance and
  high-fan/high-resistance corners, then evaluate each root with the same
  resistance that produced it.
- `OperatingPoint::correlation_handoff` converts one branch flow to typed mean
  velocity, computes Reynolds through role-tagged `fs-regime` dimensions, and
  produces `fs-convection::CorrelationInputs` without discarding evidence.
- The test-only plate-fin fixture composes that public seam through three
  Shah-London rectangular-channel cards and the real `fs-conduction` solve:
  one narrow simultaneously developing CWT table slice drives fin flanks,
  while the fully developed CWT/CHF limits provide companion checks and the
  CHF base-floor row.
  Its base plate and three fins are separate tetrahedral blocks joined by
  three named, card-backed bondlines; internal fin flanks and channel floors
  remain separately named Robin traces. Channel area, wetted perimeter,
  hydraulic diameter, aspect ratio, and length ratio all derive from the same
  mesh constants.
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
   uncertainty authority. Coefficient uncertainty alone never creates a
   validated-domain card. An attached loss domain has at least one named,
   finite axis and an exact retained source-version identity.
5. A complete operating-point result has exactly one `Certified` interval root
   and no `Possible` root boxes for the nominal declared model.
6. Manufacturer/loss/leakage uncertainty is attached as model-form `Estimate`
   evidence. It is never relabeled as a rigorous physical enclosure.
7. Terminal branch order follows deterministic network traversal, and all
   provenance hashes are order-stable and bind the complete fan curve, source,
   tolerance, fan-bank configuration, recursive network topology, loss data
   — INCLUDING each element's declared `regime_validity`, axis names and
   bounds together, with absence encoded as a state distinct from any declared
   domain — and explicit leakage identity. Equivalently: two `LossElement`
   values that compare unequal never share an operating identity, which
   `loss_elements_that_compare_unequal_never_share_an_identity` checks
   exhaustively over the constructible variants. The `regime_validity` clause
   is called out by name because omitting it silently falsified this invariant
   (bead `frankensim-yq435`): that field is the sole input to
   `regime_audit_cards`, so two networks differing only there reach a reviewer
   as an admitted result and a demoted one while hashing identically.
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
11. Pressure uncertainty endpoints come from the independently solved
    low-fan/low-resistance and high-fan/high-resistance corners. A flow root is
    never recombined with the opposite resistance bound to manufacture a wider
    but physically unreachable band.
12. Every successfully admitted `AirPath` and every published march or solve
    result keeps its represented conjugate reductions inside the finite `f64`
    domain. Path admission checks the derived positive capacity rate, segment
    conductance, NTU, effectiveness terms, and total wetted area. The march
    checks its temperature difference before products, every segment heat, and
    the running air-heat total. The coupled loop checks solid totals, the
    air/solid difference, each area-weighted residual contribution, the running
    weighted sum, and its final division. Final assembly independently
    rechecks every region difference and both heat-rate totals before it can
    publish a balance. Zero remains legal for signed heat rates, differences,
    and residuals; it refuses where the model mathematically requires a
    strictly positive denominator or exchange factor. The separately
    documented infallible decomposition attachment is a poisoned diagnostic,
    not a published finite-result guarantee.

## Error model

`AirflowError` refuses malformed curves, invalid tolerance or speed domains,
empty network groups, zero/invalid resistances, stall operation, absent curve
intersections, incomplete/ambiguous root searches, unknown branches, and bad
convection-handoff inputs. Loss-card attachment additionally refuses an
unconstrained/unusable validity box or malformed projected identity. Caller
input does not panic. Conjugate path aggregates that must be positive refuse as
`InvalidConjugateInput` with their exact field; signed runtime arithmetic
refuses as `NonFiniteCoupling` with the exact producing stage. A finite row
cannot disappear into a non-finite multi-row total.

## Determinism class

Interpolation and network traversal have fixed order. Square roots use
`fs-math::det`; the numerical root uses deterministic `fs-ivl` subdivision.
Repeated execution is checked against the complete public `OperatingPoint`
artifact, with explicit bit checks on the nominal values and certified root
endpoints. Results are intended to be bit-stable on the same ISA. The solver
has no internal parallel scheduler, so worker-count invariance is not a
separate execution mode. Cross-ISA G5 evidence is not yet retained.

## Cancellation behavior

Curve evaluation and each finite network reduction are bounded scalar work.
The interval search has an explicit 65,536-box ceiling and returns a structured
refusal rather than running without bound. No asupersync cancellation poll is
required for that rung.

`conjugate::solve_conjugate` DOES poll: every outer iteration opens with
`Cx::checkpoint`, so a cancelled exchange stops at an iteration boundary and
returns `AirflowError::Cancelled` carrying the index of the iteration that had
not yet run. Each outer iteration is also a resume boundary:
`ConjugateIteration::reference_temperatures_k` is the complete state that
iteration was solved against, `AirflowError::Cancelled` carries that vector so
a caller can resume from the refusal itself without having retained a history,
and `solve_conjugate_from` restarts from either. The poll happens at iteration
entry AND after the exchange that meets the criterion, so a cancellation
requested during the final solid solve refuses instead of publishing success.
A resumed run restarts iteration numbering, the energy audit, and the history
at zero: it is a fresh exchange seeded from a retained state, not a
continuation of the original record.
Under `Relaxation::Fixed` the resumed tail reproduces the uninterrupted run
bitwise; under `Relaxation::Aitken` the scalar relaxer's own omega history is
not carried across the resume, so the tail is a valid continuation and not a
bitwise replay.

## Unsafe boundary

None. Workspace unsafe-code denial applies.

## Feature flags

None.

## Conformance tests

- G0 fan interpolation, monotonicity refusal, and series/parallel resistance
  composition;
- G0 identical-fan series/parallel affinity identities;
- explicit stall refusal and a sign-changing, unique interval root bracket;
- three declared fan speeds, complete-artifact repeat replay, and
  leakage-resistance sensitivity;
- clause-addressed fan-law authority and a pressure-corner falsifier that
  rejects cross-pairing each solved flow with the opposite resistance;
- typed branch velocity/Reynolds handoff into `fs-convection`;
- G0/G3 plate-fin end-to-end wiring at three fan-solved speed ratios:
  card-backed base-to-fin bondlines, separately attributed fin-flank and
  base-floor Robin heat, an explicit channel-mean bulk-temperature closure,
  a wrong-hydraulic-diameter falsifier, and structured refusal of a
  fan-solved Reynolds point above the laminar-card ceiling. At fixed geometry
  and fluid properties, the developing source-table row makes fin-flank Nu
  and scalar `h` increase with each solved Reynolds point, while the companion
  fully developed limits correctly keep scalar `h` invariant. Evaluation
  provenance binds every changed point;
- semantic-identity separation when uncertainty authority changes without
  changing the nominal operating point.
- G0 sharp-edged-orifice vent lowering: closed-form resistance derivation
  with an independent operation order, allowance coverage of both sourced
  `Cd` endpoints, non-physical-input and empty-name refusals, plateau regime
  card identity/axis/bounds, and a fully vent-lowered enclosure (two parallel
  orifice vents plus an orifice leakage seam) solving to a certified
  operating point with the leakage branch reported.
- actual loss-card projection from an explicitly validated `LossElement`,
  including refusal of unconstrained pseudo-validity and an isolated
  out-of-domain loss Reynolds perturbation that demotes all seven E05.10
  records and rebinds their model-form terms to the exact receipts;
- G0 E05.10 fixture emitting all five QoI families and seven records (the
  uniformity family has mean, spread, and face-mean standard deviation), each
  with a complete eight-term budget and term-by-term provenance rendering;
- deterministic region-order/tie-break equivalence, missing-requirement and
  malformed-region refusals, G3 upstream-envelope widening monotonicity, and
  source-only identity rebinding for fan power and margin.
- G1 conjugate air path: the single-segment march against the closed-form
  heated channel `T_w - (T_w - T_in)e^(-NTU)` evaluated through an independent
  `f64::exp` path; the defining identity `h A (T_w - T_ref,eff) == Q` across
  four decades of `h` and three of area; and physical-interval containment up
  to `NTU ~ 25`.
- G3 segment-refinement invariance: splitting one uniform-wall channel into
  `N` equal segments leaves the outlet temperature and total heat unchanged for
  `N` up to 256, EXACTLY, because `e^(-NTU) = (e^(-NTU/N))^N`. The rejected
  arithmetic-mean model is built inside the test and shown to fail the same
  check with a defect that grows with `NTU` (sub-millikelvin at `NTU ~ 0.6`,
  above 10 K at `NTU ~ 20`), so the invariance is demonstrably sensitive rather
  than merely satisfied.
- G1 coupled fixed point against a manufactured solution: over a lumped solid
  the exchange has closed forms for both the outlet air temperature
  (`T_in + P/(m_dot c_p)`, a global energy statement) and the solid temperature
  (`T_in + P/(m_dot c_p (1 - e^(-sum NTU)))`); neither is evaluated by the
  driver.
- convergence battery: monotone residual contraction, Aitken reaching the same
  fixed point, and a typed non-convergence refusal that withholds the reached
  temperatures.
- G4 drills: cancellation at an outer-iteration boundary with the solid side
  proven not to re-run, and checkpoint resume reproducing the uninterrupted
  tail bitwise.
- G0 derived-arithmetic refusals: positive mass-flow and heat-capacity inputs
  whose product underflows, segment conductance/NTU range failures, overflowing
  total wetted area, per-segment and multi-segment air heat, multi-region solid
  heat, the air/solid difference, both area-weighted residual stages, and every
  independently recomputed converged-balance reduction retain exact diagnostic
  attribution.
- flux-balance sensitivity: the per-region audit closes to 2.4e-11 W on a 1.2 W
  fixture AND is shown to open under an injected dropped-face fault that the
  converged temperatures alone cannot reveal.
- G0/G3 conjugate end-to-end (`tests/conjugate_e2e.rs`): a solved fan operating
  point drives a validity-gated `CircularDuctLaminarCwt` card into a real
  `fs-conduction` FEM solve at `Re ~ 854`, `NTU ~ 2.35`, converging in 32 outer
  iterations; the air carries exactly the dissipated 1.2 W, the wall
  temperature rises monotonically downstream with a 14.9 K stream-wise tilt, a
  faster fan lowers peak, tilt, and air rise together, and the exchange replays
  bitwise. A companion test solves the SAME bar under the old one-way declared
  ambient and shows the coupled peak is 8.4 K hotter, so the coupling is not a
  refinement of the previous model but a materially different answer.
- CHT transfer: `restrict o prolongate = identity` on the coarse space, and the
  property `Refine1d` cannot state — prolongating a wall state onto the refined
  path preserves the air-side outlet and heat rate exactly at refinement
  factors 2, 5, and 16.

## No-claim boundaries

- Fan points and tolerances are caller-supplied. Synthetic fixtures prove API
  and algebra behavior, not a manufacturer product's performance.
- AMCA Publication 201-02 (R2011), section 6.5.1, justifies the retained
  speed-scaling relation, not arbitrary operation: the curve's declared
  speed-ratio domain remains binding.
- Repeat replay proves the current scalar implementation is stable on the
  exercised same-ISA path. It does not claim an internal worker-count mode that
  the crate does not implement, or cross-ISA equality.
- The plate-fin fixture is synthetic and one-way coupled. It does not solve
  channel momentum or fluid energy: its mean bulk-air reference temperature
  is the midpoint implied by declared base power, constant density/heat
  capacity, and the solved branch flow. It is neither conjugate CFD nor
  manufacturer or experimental validation. The `conjugate` module supersedes
  that closure for callers who want a two-way exchange; the plate-fin fixture
  is deliberately left as-is so the one-way and coupled answers stay
  separately inspectable.
- The convergence criterion is TWO gates, not one. A reference-temperature
  residual in kelvin cannot bound an interface heat rate in watts without a
  conductance, so the interface-imbalance gate is a separate declared
  obligation and a run that meets the temperature criterion while its worst
  per-region imbalance exceeds it refuses with
  `AirflowError::ConjugateBalanceUnclosed`. The dropped-face fixture proves the
  gap is real: it reaches its temperature fixed point and is caught only by the
  watt gate.
- The watt gate is HYBRID, not a constant. The admitted imbalance is
  `max(balance_tolerance_w, balance_relative_tolerance * scale)` with
  `scale = max_j max(|Q_solid,j|, |Q_air,j|)` from the response being judged,
  and the refusal reports that effective threshold. An absolute constant alone
  is scale-dependent in both directions — a microwatt-scale wiring fault hides
  under a watt-scale constant (measured: the dropped-face fault at 2e-6 W
  produces ~2.5e-7 W of signal), and kilowatt-scale floating-point residue can
  spuriously refuse — and a mW-vs-W unit mistake shrinks signal and gate
  together, which is exactly when the gate is needed most.
  `balance_relative_tolerance` admits `[0, 1)`; `0` restores a purely absolute
  gate. The gate still cannot certify anything smaller than the floor on an
  interface whose true scale is zero.
- Relaxation factors are ADMITTED, not merely finite. `Fixed` requires
  `0 < ω ≤ 2` and `Aitken` requires `0 < ω_init ≤ ω_max`; an inert `ω = 0`
  or beyond-over-relaxation factor refuses at admission as
  `InvalidConjugateInput` instead of stalling for `max_iterations` solid
  solves and blaming convergence. A reference vector the relaxation itself
  overflows refuses as `NonFiniteCoupling` at the "relaxed reference
  temperature" stage, attributed to the relaxation rather than to the next
  solid solve. `ω_max` is otherwise uncapped and is the caller's own risk
  budget.
- `SolidRegionState::mean_reference_temperature_k` is an OPTIONAL wiring
  claim. `Some` is held to the reference vector the driver sent (within
  1e-9 relative, kelvin-denominated) and refuses as
  `ReferenceTemperatureMismatch`; `None` makes no claim and skips the check —
  it is a no-claim boundary, not a wiring proof. `from_robin_flux` always
  claims, because `fs_conduction::RobinFlux` carries the value for free. The
  watt gate cannot see a reference-wiring error smaller than `hA` times its
  tolerance; this check is denominated directly in kelvin.
- `SegmentRefinementTransfer::{restrict, prolongate}` POISON a state whose
  arity disagrees with the bound path: the infallible `fs_ladder::Transfer`
  contract cannot refuse, so the mismatch returns all-NaN at the correct
  output arity rather than silently truncating the downstream-most entries or
  restricting a partial block with full-block weights. Downstream response
  gates refuse non-finite walls.
- "Exact" in this crate's conjugate documentation means exact in REAL
  arithmetic. Measured floating-point bands on the current fixtures: segment
  refinement holds the outlet to 1e-10 K and the heat rate to 1e-8 W up to
  N = 256; the transfer round trip lands within 1e-9 K rather than bit-exact,
  because the restriction evaluates a ratio of geometric sums; the coupled
  per-region balance closes to 2.4e-11 W on a 1.2 W fixture.
- `SegmentRefinementTransfer::restrict` is the DOWNSTREAM-WEIGHTED mean
  `W = (1-a) * sum_i a^(f-1-i) W_i / (1 - a^f)` with `a = e^(-NTU/f)`, not a
  plain mean. A plain mean is outlet-preserving only on the uniform blocks
  `prolongate` emits, so a round-trip test cannot distinguish them; on a
  nonuniform fine state the plain mean is wrong by 15.65 K at NTU 2 over four
  sub-segments. The nonuniform check in `tests/conjugate.rs` is the one that
  discriminates.
- The conjugate exchange is NOT CFD. Air is a stream-wise chain of well-mixed
  1-D segments: no lateral mixing, recirculation, buoyancy, heating-driven flow
  redistribution, or momentum coupling back to the operating point. The
  `RANS` and `LES` rungs remain declarations.
- `h` is FROZEN across the outer loop. The coupling variable is the reference
  temperature vector alone; air properties are not re-evaluated at the drifting
  film or bulk temperature, so a temperature-dependent `h` is outside the
  model. A caller who wants that must re-run the driver with a new `AirPath`.
- Relaxation is SCALAR. `fs_couple::AitkenRelaxation` is a scalar delta-squared
  relaxer; the driver projects the vector residual onto one area-weighted
  scalar and applies a single omega to the whole vector. This is not a vector
  interface accelerator. On this seam the composite gain is
  `(1 - eps/NTU)` times the solid's gain, both strictly below one, so plain
  staggering already contracts — the retained Aitken path is wiring, and no
  claim is made that relaxation rescues a divergent case, because the battery
  contains none.
- The per-region flux balance is a WIRING falsifier, not a conservation proof.
  Once `T_ref,eff` is defined so that `h A (T_w - T_ref) == Q`, the per-region
  balance is an algebraic identity at the fixed point for a uniform `h`; it
  catches a dropped face, a mis-bound region, or a wrong mass flow, and it
  catches nothing about physics. The independent cross-check is
  `decomposition_residual_w`, which compares the declared regions' sum against
  `fs-conduction`'s own whole-domain `robin_out_w` loop and so detects a Robin
  face owned by no region or counted twice. That is also wiring. The physical
  content lives in the closed-form and refinement-invariance checks above.
- `ConjugateSolution::with_decomposition_cross_check` is an INFALLIBLE
  compatibility surface. It cannot return a structured refusal without a
  public API change. A non-finite whole-domain Robin total or subtraction
  poisons `decomposition_residual_w` with `NaN`, which withholds a finite
  diagnostic but does not identify which operand was invalid. Callers needing
  exact refusal attribution must validate both operands and ensure their
  subtraction is representable before attaching the total; this method alone
  is not an admission gate.
- `AirSegment::ntu(capacity)` is also an INFALLIBLE compatibility helper. It
  performs unchecked raw arithmetic for an arbitrary caller-supplied capacity
  and can return zero or a non-finite value. The checked `AirPath` admission
  and march paths validate capacity, conductance, NTU, and effectiveness before
  using them; callers must not treat the standalone helper as an admission
  receipt.
- The seam is typed, not ledgered. `conjugate::SEAM_PORT_KIND` declares the
  `fs-couple` port kind and its effort dimension is checked against
  `Temperature::DIMS`, and `fs_couple::EnergyAudit` records every exchange's
  imbalance with NaN poisoning. The window-balance ledger path
  (`BoundaryTemperatureReference`, `WindowEvidenceRef`) is NOT wired, so no
  ledgered entropy accounting is claimed.
- `SegmentRefinementTransfer`'s fine side is a REFINED CORRELATION STATE, not a
  RANS field. It defines the correlation rung's own refinement semantics and
  the state shape a RANS rung would have to accept; it is not evidence that a
  RANS rung exists. The `RANS` and `LES` relative-cost hints in the CHT ladder
  remain unmeasured declarations, because those rungs do not exist to measure.
- The conjugate fixtures carry no maturity registration, retained corpus
  receipt, machine fingerprint, or experimental comparison. No capability
  level may be inferred from them.
- The fin-flank row is the narrow Shah-London Chapter VII, Table 52 slice for
  smooth simultaneously developing rectangular flow at `alpha*=0.5` and
  `Pr=0.72`. Its lower-`Gz` bridge is a declared engineering interpolation,
  not a published curve. It does not model a commercial plate-fin array,
  shroud/bypass flow, conjugate heat transfer, or interrupted-fin effects.
- The companion rectangular CWT/CHF rows are fully developed analytic limits.
  Reynolds is a validity axis but not a formula variable for those rows, so
  their `Nu` and scalar `h` stay invariant when only fan speed changes.
  Evaluation provenance still binds the changed Reynolds point, and a
  sufficiently high fan-solved point refuses rather than extrapolating.
- Summing named Robin heat rates back to the whole-domain Robin total is only
  a boundary-wiring falsifier; both reductions use the same face quadrature.
  The wrong-`D_h` case likewise proves the geometry handoff is load-bearing,
  not that any hydraulic diameter produces validated hardware physics.
- Piecewise-linear interpolation and quadratic loss coefficients do not model
  swirl, recirculation, acoustic interaction, thermal buoyancy, compressibility,
  fouling, or installation system effects.
- A `LossElement` with no attached validity remains an estimate with no regime
  card. An attached domain is only the source owner's stated operating box; it
  does not validate the quadratic loss law or turn its uncertainty into a
  rigorous enclosure.
- The sharp-edged-orifice lowering models a clean thin sharp-edged opening
  only. It is not a discharge model for louvered, screened, filtered, or
  ducted vents, claims no compressibility or installation effects, and its
  declared allowance is an engineering spread over sourced `Cd` values, not a
  published confidence interval. Density is a caller declaration; deriving it
  from an operating envelope is the lowering stage's obligation, and no
  retained orifice measurement corpus exists.
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
