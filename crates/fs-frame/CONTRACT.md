# fs-frame CONTRACT

## Purpose and layer

Layer: L6 (HELM). Flagship 2 (plan §15.2, bead mye.3): the
SEISMIC-MINIMAL building frame, SMOKE TIER — minimum material with
certified fragility, run end-to-end through crates that each carry
their own battery: fs-truss (layout LP + sizing), fs-solid/fs-material
(fiber hysteresis), fs-scenario (Kanai–Tajimi ensembles), fs-eproc
(anytime-valid confidence sequences), fs-uq (MLMC levels), and fs-robust
(canonical empirical CVaR risk algebra).

## Public types and semantics

- `layout::layout_and_size(..., cx) -> Result<LayoutReport, LayoutError>`:
  admits physical parameters and a strictly increasing positive catalog, then
  constructs the immutable ground grid and support/load case through
  `fs-truss`'s bounded, cancellation-aware APIs. It assembles the LP only after
  exact shape, sparse-resource, numerical, and surviving-load checks. The
  PDHG layout LP is solved at σ_y = 1 (yield stress scales only the
  objective, never equilibrium — a 250 MPa σ_y measurably stalled the
  primal-dual scaling, objective separation stuck at 1.0), relative
  primal/dual objective separation + equilibrium residual are diagnostics, and
  physical returned-iterate volume is rescaled on report. Independently, the
  shared `fs-truss-e2e` promotion gate binds an exact-feasibility repair and
  outward-checked dual to the LP/settings/iterates; verified optimum endpoints
  are outward-divided by physical σ_y. It then calls
  fs-truss `size_and_snap` (Euler floors, catalog up-snap,
  mandatory post-prune equilibrium refit, member code rows).
- `history::StoryFrame`: single-story, two fiber-hinge columns —
  concentrated plasticity: drift x → hinge curvature x/(h·l_p), the
  TRUE Mander/Menegotto–Pinto fiber section returns the moment, story
  shear V = 2M/h. Newmark average acceleration with Newton on the
  fiber tangent, one commit per step. SI units (probed k₀ ≈ 4.4e7
  N/m, V_y ≈ 4.8e5 N; the 280 t default mass gives T ≈ 0.5 s).
- `history::GroundMotion` + `StoryFrame::run_checked` bind a borrowed
  acceleration record to explicit SI units, its sample interval, and sample/Newton
  limits. The checked path returns paired relative-displacement and actual
  path-dependent fiber restoring-shear histories plus their peak absolute
  values, maximum final equilibrium residual, and maximum Newton work. The
  restoring channel excludes the viscous contribution and is not a full
  support reaction. Samples have step-end semantics: a table with an explicit
  `t = 0` row must handle that row before constructing the advancing sequence.
  The path validates the whole record, reserves both output channels, and
  integrates a staged clone; a refusal publishes neither response nor committed
  hinge state. The original `run(&[f64], dt)` remains the compatibility surface
  for the synthetic smoke studies and preserves its fixed 30-correction behavior.
  `frame-008` additionally pins the PEER/OpenSees Gallery Peknold El Centro 1940
  NS record (1,559 samples, 0.02 s, values in g), converts by exactly
  9.80665 m/s² per g, and produces a bounded diagnostic response. Its source
  revision, upstream/repository hashes, terminal-LF normalization, license,
  external model, and comparison authority are retained beside the fixture.
- `fragility::e_stopped_fragility` → `FragilityReport`: exceedance
  P(peak drift ratio > limit) over an fs-scenario ensemble, estimated
  by an fs-eproc Gaussian-mixture confidence sequence (σ = ½ is the
  HARD sub-Gaussian bound for indicators, not a plug-in); the study
  stops itself when the radius is decision-grade — validity AT the
  stopping time is the CS's construction, not a correction. An fs-uq
  MLMC report over dt levels rides along as level-design evidence. The
  smoke-tier entry point accepts only a non-empty Kanai–Tajimi
  ground-acceleration ensemble; wind spectra and material-parameter bands are
  refused before their untyped realization vectors can reach the integrator.
- `cvar::empirical_cvar` and `cvar::cvar` are direct re-exports of the
  canonical `fs-robust` report and scalar surfaces, respectively. The report
  retains deterministic VaR/minimizer and fractional-boundary metadata; empty
  or non-finite losses and beta outside `(0,1)` are structured refusals.
- `cvar::cvar_mass_min` → `CvarDesign`: Rockafellar–Uryasev empirical
  CVaR from `fs-robust`, bisection on the section scale (monotone at smoke
  scale), catalog UP-snap with an independent CVaR re-check.
  `ensemble_cvar` exposes the monotonicity probe. These smoke-tier
  orchestrators generate their own finite, non-empty loss sets and treat a
  canonical risk-algebra refusal as an internal programmer-contract defect.

## Invariants

1. Layout LP diagnostics: objective separation 3.4e-7, equilibrium residual
   4.1e-7 on the smoke fixture; returned-iterate volume positive and
   physically rescaled (frame-001). Those diagnostics do not confer proof
   strength. A separate Neumann-repaired exactly feasible primal and outward-
   feasible scaled dual produce the finite physical `Color::Verified` optimum
   interval; numerical proof unavailability remains `Estimated`.
2. Sizing: post-prune equilibrium refit 1.8e-13, every member code
   row passes post-snap (frame-002).
3. Dynamics: elastic runs do not ratchet over 10× duration (Newmark
   average acceleration, fiber tangent Newton); yielding cycles
   dissipate positive hysteretic work through the fibers at 3.3%
   peak drift (frame-003). The checked surface reproduces the compatibility
   displacement path bit-for-bit on an admitted fixture, retains one finite
   restoring-shear value per displacement, reports peaks from those exact
   channels, and leaves the frame bit-replayable from its pre-call state after
   a nonlinear refusal.
4. Fragility: the CS at the DATA-DEPENDENT stop covers the fixed-N
   reference (p_ref 0.105 inside 0.098 ± 0.120 after 163/200
   members); the threshold discriminates (16/200 exceedances); the
   e-stop saving vs fixed-N is measured and LEDGERED (18% on the
   smoke fixture) (frame-004).
5. CVaR: monotone in the section scale (0.408 → 0.0147 across the
   catalog range); bisection + snap yields a feasible design under
   the independent re-check (frame-005).
6. Replay: bitwise-identical reruns; budget exhaustion reports
   honest indecision (no early stop claimed); infeasible CVaR limits
   fire the diagnostic instead of returning a design (frame-006).
7. Recorded-input integrity: the normalized Peknold fixture has 16,032 bytes,
   FNV-1a identity `9eda012ee82cb084`, 1,559 finite samples, 0.02 s spacing,
   and 0.31882 g peak magnitude. Its checked response publishes 1,558 step-end
   displacement/restoring-shear pairs after explicitly designating source row 0
   as initial forcing, or returns a structured refusal (frame-008).

## Error model

`layout_and_size` returns `LayoutError::Construction` for malformed geometry,
physical parameters, catalog, resource excess, allocation refusal, or observed
cancellation, and `LayoutError::Solver` for rejected PDHG/diagnostic state. It
returns `LayoutError::Certificate` for malformed proof state, allocation
failure, or cancellation during receipt construction/binding. A sound numerical
proof refusal is published only as Estimated, never as finite bounds. It does
not publish a partial layout. Direct empirical-CVaR calls return
`RobustError`; smoke-tier orchestration contracts still panic with teaching
messages when internally generated losses violate that contract, an ensemble
spec is malformed, or a CVaR study is infeasible (the drill gates the
diagnostic). Statistical outputs carry their own uncertainty: the CS radius and
stopping state ARE the answer's quality statement.

`StoryFrame::run_checked` returns `HistoryError` for an empty or non-finite
record, invalid SI time step or story parameter, sample-budget excess, invalid
Newton limits/tolerances, response-allocation failure, non-finite dynamic state,
or exhaustion of the per-sample correction limit. Convergence currently means
both the absolute Newton displacement correction and the absolute dynamic-
equilibrium residual are below their caller-supplied finite, positive
tolerances, including a final post-correction residual check. The checked path
is all-or-error with respect to both the returned histories and the frame's
committed state; the compatibility `run` retains the older unchecked fixed-cap
contract.

## Determinism class

Bit-deterministic per platform: Philox-streamed ensembles
(fs-scenario), fixed iteration orders, deterministic solvers.
frame-006 pins bitwise study replay; frame-003 additionally pins checked versus
compatibility Newmark displacement bits and post-refusal frame replay.

## Cancellation behavior

Ground construction and LP assembly poll the explicit `Cx` at deterministic
bounded strides and return structured cancellation before publication. The
outward certificate stage polls the same context through repair, verification,
identity binding, and publication. Later fixed solver/dynamics loops remain
synchronously bounded by iteration/member budgets. Checked history integration
also has explicit sample and per-sample Newton limits and stages all mutable
fiber state until success; it is synchronously bounded but does not yet accept
an explicit `Cx`. The e-stop is itself the anytime-cancellation story: stopping
at ANY member count leaves a valid interval.

## Unsafe boundary

`#![deny(unsafe_code)]` via workspace lints; no capsules.

## Feature flags

None (the smoke tier ships enabled; heavier tiers will gate).

## Conformance tests

`tests/battery.rs`: frame-001 LP diagnostics and physical outward optimum bounds;
frame-002 sizing code
rows plus pre-cancelled construction refusal; frame-003 elastic stability,
hysteretic dissipation, bounded/atomic displacement + restoring-shear response
history, and admission refusals;
frame-004 e-stopped fragility coverage + ledgered savings; frame-005
CVaR monotonicity + design; frame-006 replay, infeasibility, and structured
canonical-CVaR refusal drills; frame-007 wrong-physics and empty-ensemble
refusals at the realization boundary; frame-008 pinned El Centro record,
provenance/license lint, checked diagnostic response, and diagnostic-only
OpenSees displacement rounding envelope.

## No-claim boundaries

- SMOKE TIER geometry: one story, two identical fiber-hinge columns.
  Distributed-plasticity frames (fs-solid `ForceBasedElement`
  columns), multi-story assemblies, and joint modeling are recorded
  successors.
- One source- and license-pinned Peknold El Centro NS record now exercises the
  checked integrator, but it is a diagnostic fixture, not a recorded-motion
  suite or spectral-matching lane. The OpenSees Example 3 peak roof displacement
  is retained with its source-rounding interval, but its distributed-column
  model is materially non-equivalent to the shipped concentrated-hinge fixture,
  so the interval is not applied as a FrankenSim acceptance gate. Example 3
  publishes no total dynamic support reaction; restoring shear is not mislabeled
  as that missing oracle. Admitted model mapping, a support-reaction artifact,
  and fs-vvreg authority remain pending.
- Newmark average acceleration ships; the fs-time VARIATIONAL
  integrator swap (the plan's long-duration drift story) is a named
  successor — the 10×-duration stability gate stands in.
- The Michell-continuum catalogue comparison remains LEDGERED PENDING
  exactly as in the fs-truss contract; `:oracle (michell :tol 0.08)`
  lands with the fs-fab spec constants.
- Million-member ground structures, SOCP layout, trust-region
  Newton–Krylov multi-variable sizing, and arclength global-buckling
  sweeps: perf/full-tier scope.
- MLMC here is a LEVEL-DESIGN REPORT, not the estimator of record
  (the CS is); driving the fragility estimate itself through MLMC
  with e-stopping per level is full-tier scope.
- No `explain()` artifact chain yet — the diagnostic and evidence fields are
  the auditable record; the fs-ledger integration is staged with the
  study-program (Appendix C) runner.
