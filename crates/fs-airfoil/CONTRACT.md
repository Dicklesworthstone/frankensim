# fs-airfoil — CONTRACT

Bead: frankensim-wf-root-guzez.5.1.1 (E4.0a, Wright Flyer program).
Spec: COMPREHENSIVE_PLAN_FOR_REAL_TIME_WRIGHT_FLYER_SIM_WITH_FRANKENSIM.md §5.2.1
(ROUND 6 steady state). Conventions: frame-conventions-v1 + geometry-conventions-v1
(frozen, E1.4/E1.9).

## Purpose and layer

L2. Generic airfoil section machinery consumed by wing, canard, rudder AND
propeller crates (sections live at L2 precisely so fs-wing/fs-airscrew need no
L3→L3 edge). E4.0a ships: section geometry with provenance binding, analytic
thin-airfoil and separated flat-plate baselines, exact wind↔body decomposition
+ moment-reference transfer, and the admission domain. E4.0b adds
provenance-bound coefficient tables + regime-partitioned constrained B-spline
residuals; E4.0c adds coherent-draw uncertainty, OOD refusals over fitted
domains, and indicial kernels.

## Public types and semantics

- `SectionGeometry { chord_m, camber_ratio, dossier_record, digitization_class }`
  — geometry is provenance-bound (dataset re-expression rule): an empty
  dossier record refuses at admission.
- `SectionCoefficients { cl, cd, cm_quarter }` — wind-axis; `cm` positive
  nose-up per moment-signs-v1.
- `NormalAxial { cn, ca, cm_ref, x_ref_over_c }` — the through-stall
  representation; explicit moment reference station.
- `thin_airfoil(α, f/c, log10 Re)` — classical exact parabolic-camber results:
  cl = 2π(α + 2f), cm_c/4 = −πf, cd = 0 (inviscid baseline).
- `flat_plate_separated(α, log10 Re)` — cn = CD90·sin α with CD90 = 1.98,
  ca = 0, cp walking c/4 → c/2; the low envelope of post-stall.
- `wind_to_body` / `body_to_wind` — exact rotation pair; `transfer_moment` —
  cm_B = cm_A + cn·(x_B − x_A).
- `fit::BsplineAxis` — clamped uniform cubic axis (n_coef = 1 ⇒ degenerate/
  constant axis); `fit::ResidualSurface` — tensor-product surface over
  (α, log Re, δ) with declared `DiffConstraint`s, penalized-LS `fit` (dense
  Cholesky, deterministic ridge) and FAIL-CLOSED constraint verification
  (`fit-constraint-violated`, never silent projection);
  `fit::verify_regime_continuity` — C⁰ face check between abutting patches.
- `table::CoefficientTable` — provenance-bound (dossier record + validated
  `ConventionBlock` against the frozen E1.4 ids), SurfaceKind-separated
  (Wing/Canard/Rudder/Prop), regime patches tiling α; `eval` refuses
  outside the FITTED partition (`alpha-outside-table`) — distinct from the
  global admission domain; `eval_strict` additionally refuses when the
  covering patch's log Re / δ box is exceeded
  (`query-outside-fitted-domain`, box stated) instead of the spline clamp.
- `uncertainty::UncertainSurface` — mean + low-rank coefficient modes
  (≤ 16); `realize(realization_id)` is the ONLY uncertain-query path
  (coherent draw; per-query independent intervals are impossible by
  construction); weights are a pure fs-blake3 function of the id under
  `org.frankensim.fs-airfoil.coef-realization.v1`.
- `indicial::IndicialKernel` — two-pole φ(s) = 1 − a₁e^(−b₁s) − a₂e^(−b₂s)
  with registered constants `WAGNER_JONES` (φ(0) = 0.5) and
  `KUSSNER_2POLE` (ψ(0) = 0); `IndicialState` advances by the EXACT
  diagonal exponential (sub-step composition is bit-tight);
  `reduced_time_increment` implements the CHORDWISE clock
  ds = 2·U_conv·dt/c — freezes at U_conv = 0, REFUSES reversed flow
  (`indicial-flow-reversed`), never |U|.

## Invariants

- Angles in radians everywhere (units rule, frame-conventions-v1).
- Baselines are pure functions: no state, no allocation on the query path
  (geometry admission allocates only in the refusal branch).
- Thin-airfoil odd symmetry at zero camber is bit-exact; flat-plate cn/cm are
  odd in α; decomposition round-trips to 1e-14.
- Analytic results match classical closed forms to 1e-15 at the tested points.

## Error model

Typed `Refusal { code, message, ranked_repairs }` — never panics on the query
path. Codes: `non-finite-input`, `alpha-outside-domain`,
`reynolds-outside-domain`, `chord-outside-domain`, `camber-outside-domain`,
`provenance-missing`; fit/table layer: `axis-domain-invalid`,
`axis-coef-count-invalid`, `constraint-axis-invalid`, `insufficient-samples`,
`fit-normal-equations-singular`, `fit-constraint-violated`,
`regime-boundary-mismatch`, `regime-boundary-discontinuity`,
`convention-block-mismatch`, `convention-block-missing`, `table-empty`,
`alpha-outside-table`, `query-outside-fitted-domain`,
`uncertainty-modes-invalid`, `realization-id-empty`,
`kernel-params-invalid`, `reduced-time-increment-invalid`,
`timestep-invalid`, `indicial-flow-reversed`. Applicability-domain refusals STATE the admitted domain.
Caps are tested at cap AND cap+1 (workspace law).

## Determinism class

Deterministic: trig routed through `fs_math::det` (libm doctrine); the
analytic-polar golden digest is pinned under
`org.frankensim.fs-airfoil.analytic-polar.v1` and the golden-bump protocol
applies.

## Cancellation behavior

All entry points are synchronous pure functions; nothing to cancel.

## Unsafe boundary

`unsafe_code = "forbid"` via workspace lints; no unsafe anywhere.

## Feature flags

None. (E4.0b/c additions stay feature-free; fidelity is data-driven, not
cfg-driven.)

## Conformance tests

`tests/fit_battery.rs`: basis partition-of-unity + Greville linear
reproduction; synthetic 3-D fit round-trip (off-grid < 1e-8); constraint
falsifier (non-monotone data under a monotone constraint refuses, and fits
without it); sample caps at n and n−1; regime-continuity twins; convention/
provenance falsifiers; pinned fit golden f25f5e76.

`tests/uncertainty_indicial_battery.rs`: coherent-draw law (same id →
bit-identical surface; two-point mode-structure prediction; anonymous-draw
+ mode-shape falsifiers; mode caps at 16 and 17); strict OOD twins (box
edge admitted, next float refused with the box stated); Wagner/Küssner
exact references (φ(0), ψ(0), monotonicity, closed-form tracking to 1e-13,
sub-step composition to 1e-14); chordwise-clock freeze/refusal battery;
pinned Wagner-trace golden 7897e5d7.

`tests/analytic_battery.rs`: thin-airfoil exact classical results (zero-lift
angle, 4πf at α = 0, 2π slope by central difference, α-independent cm_c/4,
bit-exact odd symmetry); flat-plate shape (CD90 anchor at ±90°, cp at
mid-chord, odd symmetry, separated-below-attached check); decomposition
round-trips + moment-transfer inversion + cm-about-cp vanishing; refusals at
cap and cap+1 with domain-stating messages; provenance-gate falsifier; pinned
polar golden. JSONL receipts per case.

## No-claim boundaries

- Analytic baselines are CLASSICAL MODELS, not Wright measurements; they make
  no Wright-specific claim. Wright section data arrives via E4.0b tables bound
  to dossier records.
- cd = 0 in the thin-airfoil baseline is inviscid by construction; nothing
  here estimates profile drag.
- The flat-plate post-stall shape is a declared low envelope; the
  a2-synthesized-stall record's prohibition (no Wright-specific deep-stall
  validation) binds every consumer.
- No aircraft-level claim of any kind: hinge moments, planform effects, and
  biplane interference live in fs-wing (plan §5.2 ownership).
