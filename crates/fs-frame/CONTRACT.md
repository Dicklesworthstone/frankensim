# fs-frame CONTRACT

## Purpose and layer

Layer: L6 (HELM). Flagship 2 (plan §15.2, bead mye.3): the
SEISMIC-MINIMAL building frame, SMOKE TIER — minimum material with
certified fragility, run end-to-end through crates that each carry
their own battery: fs-truss (layout LP + sizing), fs-solid/fs-material
(fiber hysteresis), fs-scenario (Kanai–Tajimi ensembles), fs-eproc
(anytime-valid confidence sequences), fs-uq (MLMC levels).

## Public types and semantics

- `layout::layout_and_size` → `LayoutReport`: ground-structure grid,
  PDHG layout LP solved at σ_y = 1 (yield stress scales only the
  objective, never equilibrium — a 250 MPa σ_y measurably stalled the
  primal-dual scaling, gap stuck at 1.0), duality gap + equilibrium
  residual as the CERTIFICATE, physical volume rescaled on report;
  then fs-truss `size_and_snap` (Euler floors, catalog up-snap,
  mandatory post-prune equilibrium refit, member code rows).
- `history::StoryFrame`: single-story, two fiber-hinge columns —
  concentrated plasticity: drift x → hinge curvature x/(h·l_p), the
  TRUE Mander/Menegotto–Pinto fiber section returns the moment, story
  shear V = 2M/h. Newmark average acceleration with Newton on the
  fiber tangent, one commit per step. SI units (probed k₀ ≈ 4.4e7
  N/m, V_y ≈ 4.8e5 N; the 280 t default mass gives T ≈ 0.5 s).
- `fragility::e_stopped_fragility` → `FragilityReport`: exceedance
  P(peak drift ratio > limit) over an fs-scenario ensemble, estimated
  by an fs-eproc Gaussian-mixture confidence sequence (σ = ½ is the
  HARD sub-Gaussian bound for indicators, not a plug-in); the study
  stops itself when the radius is decision-grade — validity AT the
  stopping time is the CS's construction, not a correction. An fs-uq
  MLMC report over dt levels rides along as level-design evidence.
- `cvar::cvar_mass_min` → `CvarDesign`: Rockafellar–Uryasev empirical
  CVaR (exact for empirical measures: t* = β-quantile), bisection on
  the section scale (monotone at smoke scale), catalog UP-snap with
  an independent CVaR re-check. `ensemble_cvar` exposes the
  monotonicity probe.

## Invariants

1. Layout LP certificate: duality gap 3.4e-7, equilibrium residual
   4.1e-7 on the smoke fixture; volume positive and physically
   rescaled (frame-001).
2. Sizing: post-prune equilibrium refit 1.8e-13, every member code
   row passes post-snap (frame-002).
3. Dynamics: elastic runs do not ratchet over 10× duration (Newmark
   average acceleration, fiber tangent Newton); yielding cycles
   dissipate positive hysteretic work through the fibers at 3.3%
   peak drift (frame-003).
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

## Error model

Programmer contracts panic with teaching messages (ensemble spec
defects, infeasible CVaR studies — the drill gates the diagnostic).
Statistical outputs carry their own uncertainty: the CS radius and
stopping state ARE the answer's quality statement.

## Determinism class

Bit-deterministic per platform: Philox-streamed ensembles
(fs-scenario), fixed iteration orders, deterministic solvers.
frame-006 pins bitwise replay.

## Cancellation behavior

Bounded synchronous loops (member budgets, iteration caps); the
e-stop is itself the anytime-cancellation story: stopping at ANY
member count leaves a valid interval.

## Unsafe boundary

`#![deny(unsafe_code)]` via workspace lints; no capsules.

## Feature flags

None (the smoke tier ships enabled; heavier tiers will gate).

## Conformance tests

`tests/battery.rs`: frame-001 LP certificate; frame-002 sizing code
rows; frame-003 elastic stability + hysteretic dissipation;
frame-004 e-stopped fragility coverage + ledgered savings; frame-005
CVaR monotonicity + design; frame-006 replay + drills.

## No-claim boundaries

- SMOKE TIER geometry: one story, two identical fiber-hinge columns.
  Distributed-plasticity frames (fs-solid `ForceBasedElement`
  columns), multi-story assemblies, and joint modeling are recorded
  successors.
- Motions are SYNTHETIC Kanai–Tajimi only; recorded-motion suites and
  spectral matching are staged with fs-scenario's data lanes.
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
- No `explain()` artifact chain yet — the certificate fields are the
  auditable record; the fs-ledger integration is staged with the
  study-program (Appendix C) runner.
