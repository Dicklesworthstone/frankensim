# fs-orbit — Periodic-orbit service (HB + shooting + continuation)

The first implementation of the fs-vmanifest I09 periodic-machine
orbit slot (music bead `frankensim-music-v8-root-3ez8g.11.1`). Music
is the first consumer, not the owner.

## Purpose and layer

Layer L2. General periodic-orbit machinery over the ISLAND-PLUS-
LINEAR-PORT problem shape: a nonlinear part evaluated pointwise in
time (reed law, device card, cubic spring) and a linear part
specified per harmonic in the frequency domain (`s I` for an ODE, a
TMM impedance for an acoustic load). Serves music lock questions,
flutter boundaries, and thermoacoustics through one problem-shaped
API — never instrument-shaped.

## Public types and semantics

- `OrbitProblem` — `dim`, `island(t, x, out)` (the AFT-sampled
  nonlinearity), `port(s)` (row-major `d x d` complex operator,
  default `s I`), `autonomous()`.
- `HbAnchor` — `Forced { omega }`, `Autonomous { omega_guess }`
  (phase anchored by `Im X_1[0] = 0`), `Backbone { amplitude,
  omega_guess }` (conservative backbone point closed by an APPENDED
  amplitude-norm equation `|X_1[0]| = a/2`; every balance row stays
  enforced and the phase stays free, so `omega` is
  phase-independent).
- `HbBudget` / `solve_hb` / `solve_hb_seeded` -> `HbOrbit` (omega,
  canonical coefficients, per-iteration residual trace, verified residual,
  final verification sample count, measured sampling defect).
- `ShootBudget` / `solve_shooting` -> `ShootOrbit` (period, orbit
  point, Floquet multipliers from the forward-differenced monodromy
  through the fs-la complex eigensolver).
- `ContinuableProblem` / `ContinuationBudget` / `continue_branch` —
  TRUE pseudo-arclength (bordered corrector on coefficients and the
  parameter jointly), so parameter folds are traversed, not stalled
  at.

## Invariants

- Harmonic truncation `N` is fixed, disclosed structure (X-Struct).
  AFT starts with `m = (4N+4).next_power_of_two()` samples and uses
  prepared det-routed trigonometric sums. Basis tables and state/RHS
  buffers are reused through Newton and finite-difference Jacobian
  evaluations; the DC port is evaluated once per residual.
- A converged candidate is re-evaluated on a finer `2m+1` grid.
  Both its residual and the norm of the change in retained balance
  coefficients must meet the same tolerance and coefficient scaling.
  Otherwise Newton resumes on the finer grid. There are at most four
  refinements; the fourth is a verification grid. Persistent disagreement
  returns `SamplingUnresolved`, and refinements also consume the existing
  Newton iteration budget.
- This sampling gate also applies to every continuation corrector.
  Accepted points retain their actual residual trace and final verified
  residual (including the arclength equation), never a fabricated zero.
  `aft_samples` records the final verification grid and
  `sampling_defect` records the measured retained-coefficient change.
- Canonical unknown ordering: DC block, then per-harmonic
  `Re`/`Im` blocks, state-major; the packing is the seed format.
- A converged residual below the relative tolerance is the ONLY
  accepted discrete-balance evidence. Autonomous/backbone HB solves
  refuse equilibrium by name (`TrivialCollapse`), checking the full
  non-DC spectrum rather than the fundamental alone. The existing
  absolute autonomous collapse threshold is 1e-9 in packed state units.
  Forced HB solves admit DC, zero, arbitrarily small, and higher-harmonic
  responses: their period is prescribed, not inferred from a cycle.
  Equilibrium deflation applies only to autonomous/backbone anchors.
  Shooting retains its separate equilibrium check.
- Conservative families carry a parity/phase nullspace: the Newton
  step regularizes a singular factorization with a FIXED Tikhonov
  jitter (deterministic); the converged residual gate is unchanged.

## Error model

Typed `OrbitError`: `BadParameter`, `NewtonStalled` (with the full
per-iteration residual trace), `SamplingUnresolved` (finest sample count,
residual and measured defect), `SingularJacobian`,
`ContinuationExhausted`, `TrivialCollapse`, `TorusSuspected` (a
non-trivial unit-circle Floquet pair with nonzero angle:
quasi-periodic drift is a NAMED no-claim in v1, never a wrong
answer), forwarded `Eigen`.

## Determinism class

Deterministic: bounded deterministic grid refinement, fixed iteration caps and step policies,
det-routed transcendentals, no RNG, no time dependence; bitwise
repeatability is a conformance case.

## Cancellation behavior

Pure synchronous computation; the heaviest unit is one dense LU per
Newton iteration. No `Cx` integration.

## Unsafe boundary

`#![deny(unsafe_code)]` via workspace lints; no capsules.

## Feature flags

None.

## Conformance tests

`tests/orbit_conformance.rs` (ob-001..ob-006) — oracles with no
music in them: the conservative Duffing backbone vs the fs-nlmodal
pinned first-order law with the band authored at the NAMED
second-order Lindstedt coefficient `(15/256) eps^2 a^4` (measured:
the HB deviation IS that term); van der Pol at mu = 1 against TWO
independent published sources (Amore arXiv:2111.12198
T = 6.663286859323130 — matched to nine digits — and
Gasull–Giacomini–Grau arXiv:1602.00113 T ~ 6.6632866) plus the
weak-mu perturbation laws; pseudo-arclength continuation TRAVERSING
both folds of the forced-Duffing S-curve with fold frequencies
within 0.2% of the independent first-harmonic scalar law; shooting
cross-validating the same vdP orbit to 1e-12 relative period with
Floquet classification (trivial multiplier 1 + transverse 8.6e-4);
refusals by name including `TorusSuspected` on a constructed
incommensurate two-center fixture; bitwise determinism.

The G1/G3 analytic regressions in `tests/harmonic_balance.rs` compare
higher-harmonic forced oscillator coefficients (including their phase and
velocity) with the closed-form transfer function, exercise small and DC
forced responses, and retain autonomous-equilibrium refusal. A seventh-order
nonlinear spring exposes the old eight-point alias (0.9968112456 amplitude
instead of the retained Galerkin answer 1.0); an independent shifted
65-point quadrature checks the corrected coefficients. Further regressions
verify truthful continuation residuals against the analytic scalar balance
and a bounded named refusal for persistently unresolved forcing.

## No-claim boundaries

- Quasi-periodic orbits (tori) are a NAMED refusal, not a claim.
- The continuation corrector is forced-anchor v1; autonomous-branch
  continuation (period as a continuation unknown) is a named
  successor.
- Collocation is described by the I09 slot but not implemented here;
  shooting is the second method.
- No a-priori truncation-error bound: `N` remains disclosed structure.
  Sampling refinement improves the retained Galerkin balance, not the
  number of harmonics or the continuum model. Agreement between finite
  grids is a measured convergence check, not a rigorous universal aliasing
  bound: sufficiently high frequencies can alias on both grids. Explicit
  input bandwidth must be admitted by domain adapters; shooting is an
  independent cross-check only for problems with an equivalent time-domain
  formulation.
- The monodromy is forward-differenced, not a consistent tangent;
  the I09 lattice's consistent-tangent element is a successor.
- Grazing/hybrid (non-smooth) islands are out of scope in v1.
