# CONTRACT: fs-thrust-e2e

CertQD-Thrust — a screened quality-diversity discovery campaign for
self-propelling point-vortex thrusters. An L6 (HELM) capstone that orchestrates
physics × illumination × fidelity-management × evidence colors × provenance.

Nothing in this crate is `Verified`. The drift comes from unchecked RK4 and the
only trust signal is an impulse-conservation residual on a different functional,
so every published drift is `Estimated` and the impulse test is called a SCREEN,
not a certificate (bead `frankensim-extreal-program-f85xj.2.30`).

## Purpose and layer

Layer L6 (HELM orchestration). Depends downward only: `fs-vpm` (L3 physics),
`fs-archive` (L4 QD), `fs-surrogate` (L4 certify-or-escalate), `fs-evidence`
(L1 colors), `fs-report` (L6 notebook). Ambition tag `[F]` (frontier synthesis;
the physics is a 2-D inviscid smoke tier).

## Public types and semantics

- `Design { gamma, d, l, ratio }` — a four-vortex quadrupole thruster; `vortices`
  (leading `±Γ` dipole + trailing `±Γ·ratio` dipole, total circulation 0),
  `budget` (`Σ|Γᵢ|`), `descriptor` (`[budget, length]`).
- `SimResult { drift, impulse_error }`; `simulate_thrust(design, steps, dt,
  core)` — RK4 vortex sim → net `x`-drift of the vortex mean + the
  linear-impulse-conservation error.
- `CampaignBudget` — full/short horizons, `dt`, core, bins, conformal `alpha`,
  `decision_tol`, `conserve_tol`, seed (the Five Explicits).
- `design_grid()` — the deterministic design sweep.
- `calibration_designs()` — the 8 designs whose short-vs-full residuals fit the
  conformal band; `calibration_support() -> CalibrationSupport` — their PER-AXIS
  support (`gamma`, `d`, `l`, `ratio`, `budget`), which is the surrogate's
  declared validity domain. `CalibrationSupport::contains(&Design)` is the exact
  predicate the campaign uses; it revalidates finiteness because `Design` is
  publicly constructible. A pinned axis has a degenerate support (`lo == hi`),
  and that is a claim, not an oversight: the residuals say nothing about how the
  short-vs-full gap moves along an axis the calibration never varied.
- `run_campaign(&CampaignBudget) -> CampaignReport` — the whole campaign;
  `CampaignReport` carries coverage, QD-score, best design + drift, screened vs
  unscreened elite tallies, full/short sim counts, steps spent vs all-full, the
  conformal band, the conservation-screened drift hull, the campaign color rank,
  and the content-addressed lab-notebook Markdown. `AtlasEntry` carries the
  elite's `conservation_screened` flag and its `rank` (always `Estimated`).

## Invariants

- PHYSICS: a four-vortex thruster with a leading dipole self-propels in `+x`
  (`drift > 0`); a converged inviscid sim conserves the exact linear impulse
  `I = (Σ Γᵢ yᵢ, −Σ Γᵢ xᵢ)` to a small error.
- CONSERVATION SCREEN (not a certificate): an escalated full sim PASSES the
  screen iff drift, impulse error, impulse scale, and `conserve_tol` are finite;
  the error, scale, and tolerance satisfy error ≥ 0, scale > 0, and tolerance
  ≥ 0; the relative error `‖ΔI‖ / Σ|Γᵢ|` is within the inclusive tolerance; and
  `drift ± impulse_error` is finite and ordered. Malformed state fails closed to
  an infinite-dispersion no-claim. A passing sim is colored
  `Estimated{estimator: "vpm-full-impulse-conserving", dispersion:
  impulse_error}`; a leaking one `Estimated{"vpm-full-nonconserving", …}`; a
  surrogate estimate `Estimated{"vpm-short-surrogate", band_half_width}`. NO
  path mints `Verified`. Conservation of `I` is not an error bound on the
  mean-`x` drift (for a zero-total-circulation quadrupole the drift is not a
  component of `I`), and `fs_vpm::simulate` is unchecked RK4 with no step-size
  control or outward rounding, so no executable enclosure exists to publish. The
  previous `Verified{drift ± max(impulse_error, 1e-9)}` band asserted a ±1e-9
  certificate whose width was chosen by the code, not derived.
- VALIDITY DOMAIN: the surrogate serves a design only if EVERY calibrated axis
  contains it — `gamma ∈ [1.0, 1.4]`, `d ∈ [0.7, 0.7]`, `l ∈ [0.6, 1.8]`,
  `ratio ∈ [0.4, 1.0]`, `budget ∈ [2.8, 5.6]` — not merely the 2-D `(budget,
  length)` descriptor hull. `d` is pinned by the calibration set, so no other
  transverse spacing is in-domain: a dipole self-advects at `~Γ/(2πd)`, making
  `d` a first-order driver of the drift the surrogate extrapolates, and it is
  invisible to `Design::descriptor`. Under the descriptor hull alone, 63 of 84
  served designs sat at spacings no calibration residual ever saw and each was
  handed the `d = 0.7` band as its uncertainty. Of the 192 swept designs, 18 are
  now in-domain (including the 8 calibration designs) and 174 escalate.
- CERTIFY-OR-ESCALATE COST: the short surrogate is used only for in-domain
  designs when the conformal band is within `decision_tol`; everything else
  escalates. Whether that is cheaper than the naive all-full sweep is
  arithmetic, not a promise: the 8 PAIRED calibration sims are charged
  unconditionally, so `steps_spent < steps_all_full` iff
  `short_sims·(full_steps − short_steps) > 8·(short_steps + full_steps)`. At the
  default budget that is `18·340 = 6120 > 8·460 = 3680`, a 3.2% saving; the
  campaign is CHEAPER only while `short_steps/full_steps < 10/26 ≈ 0.385`, and
  the reported `step_savings` metric is signed and goes negative above that
  ratio. No claim is made that surrogate-served answers match all-full answers:
  they do not, which is exactly why the campaign rank is `Estimated`. The
  default eight-residual calibration uses binary-exact `alpha = 0.125`, safely
  above the `1/(n+1)` minimum required for the retained eighth order statistic;
  it never relies on rank clamping.
- NO LAUNDERING: the campaign-level color rank is the weakest elite color
  (`min` over `ColorRank`), which is `Estimated`. The
  `conservation_screened_drift_hull` is a plain endpoint hull of the screened
  elites' `drift ± impulse_error` bands — it is NOT an `IntervalOp::Hull`
  `compose` of certified intervals, because there are no certified intervals
  here, and it is not an error bound on any drift.
- DETERMINISM: no RNG; a fixed design grid + fixed physics ⇒ the notebook content
  hash and all metrics are bit-stable across runs.

## Error model

Total functions; `run_campaign` never panics on the default/grid inputs
(`conformal_band` is fed a non-empty calibration residual set by construction).

## Determinism class

Fully deterministic (G5): the sweep, sims, archive, colors, and notebook are pure
functions of `CampaignBudget`.

## Cancellation behavior

None here (a synchronous batch campaign); the production path would poll `Cx` at
sim-tile boundaries.

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/thrust.rs` (6 cases): a four-vortex thruster self-propels and conserves
impulse; the campaign illuminates a screened diverse family (coverage/QD-score,
best drift > 0, both fidelities used, the calibration-repayment identity behind
the step saving, screened+unscreened tally, screened drift hull, no-laundering
rank, content-addressed notebook); no elite claims an interval certificate from
the impulse screen (regression, bead `.2.30`); every surrogate-served design is
inside the calibration support (regression, bead `.2.29`); the campaign is
deterministic (identical content hash + metrics across runs); and a public `+∞`
conservation tolerance cannot mint any screened elite. The in-source battery
additionally pins malformed public `SimResult` fields, invalid scales/tolerances,
band overflow refusal, valid inclusive boundary behavior, the exact bead-`.2.30`
witness (`drift 2.0`, `impulse_error 0.0` is `Estimated`, never
`Verified{2.0 ± 1e-9}`), the per-axis support values and their rejection of
off-calibration and non-finite genes, and admissibility of the default conformal
order statistic.

Representative run (default budget): 28 niches, coverage 0.438, QD-score 116.74,
best drift ≈ 9.005 (`gamma 1.8, d 0.4, l 0.6, ratio 1.0`), 28 screened + 0
unscreened elites, 174 full / 18 surrogate sims, 74,360 of 76,800 steps
(≈3.2% saving vs all-full), band half-width ≈ 0.7639, conservation-screened
drift hull ≈ [1.338, 9.005], campaign rank `Estimated`.

The numbers moved when the two defects above were fixed, and the movement is the
point. Before: 108 full / 84 surrogate sims and ≈32% savings — bought by serving
63 designs off-calibration in `d`; and 19 "Verified" + 9 Estimated elites with a
"certified drift envelope ≈ [1.3, 9.0]" — bought by reading an impulse residual
as a drift bound. No elite is surrogate-served at the default budget now: the
18 in-domain designs all sit at `d = 0.7` and never out-drift the tighter-spaced
ones, so the surrogate's honest domain contains no niche winner.

## No-claim boundaries

- The physics is `fs-vpm`'s 2-D INVISCID point-vortex core (no viscosity, no
  free surface, no body); "self-propulsion" is the vortex mean-position drift of
  a zero-total-circulation cluster, not a solid swimmer. A hybrid BEM+VPM body,
  viscous PSE, and 3-D filaments are the fuller physics, staged in `fs-vpm`.
- The surrogate is a linear time-extrapolation of a short-horizon sim; a POD/
  neural operator surrogate is `fs-surrogate`'s fuller deliverable.
- No drift here is bounded. Producing a real enclosure needs validated
  integration (interval/Taylor-model RK4, or step-doubling with outward
  rounding) in `fs-vpm`; until that exists the campaign can screen runs and
  report residuals, and that is all it claims.
- The calibration set pins `d` and samples two `gamma`, two `l`, and two `ratio`
  values, so the honest validity domain is a thin slab of the design space and
  the surrogate serves 18 of 192 designs. Widening the domain requires widening
  the CALIBRATION — new paired sims that actually vary `d` — not widening the
  test.
- The design sweep is a regular grid; a QD variation/emitter loop (CMA-ME) over
  the archive is the fuller illumination.
- The lab notebook is `fs-report`'s v0 (deterministic Markdown + reproducing IR);
  ledger persistence and semantic design diffs are downstream integrations.
