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
  default `s I`), `autonomous()`, and optional `time_domain_rhs()`.
  The latter supplies a `TimeDomainRhs<Self>` callback for the explicit
  derivative `x' = f(t, x)` in the same state and time coordinates.
  It defaults to unavailable; a plain ODE explicitly supplies
  `Some(Self::island)`, while a mass-weighted problem supplies `M^-1 g`.
- `HbAnchor` — `Forced { omega }`, `Autonomous { omega_guess }`
  (phase anchored by `Im X_1[0] = 0`), `Backbone { amplitude,
  omega_guess }` (conservative backbone point closed by an APPENDED
  amplitude-norm equation `|X_1[0]| = a/2`; every balance row stays
  enforced and the phase stays free, so `omega` is
  phase-independent). The authored `a` is a positive finite physical peak
  amplitude, with a representable half-amplitude Fourier coefficient; its
  relative tolerance must be finite and strictly between zero and one.
- `HbBudget` / `solve_hb` / `solve_hb_seeded` -> `HbOrbit` (omega,
  canonical coefficients, per-iteration residual trace, verified residual,
  final verification sample count, measured sampling defect).
- `HbScales` / `solve_hb_scaled_seeded` -> `ScaledHbOrbit` conditions this
  same Newton owner with separate physical state, equation and frequency
  magnitudes. The seed and returned orbit remain in original units. Acceptance
  bounds the absolute normalized Parseval RMS defect and independent sampling
  change; large static preloads never relax it. Both scaled diagnostics are
  retained alongside the original physical residual trace. This is the entry
  for mixed-unit or very small physical responses; the original unscaled
  entry points retain their existing coefficient-relative convergence rule.
- `HbRefinementBudget` / `solve_hb_adaptive` /
  `solve_hb_adaptive_seeded` -> `AdaptiveHbOrbit`. These additive entry
  points retain the fixed-order Newton budget, start at its harmonic order,
  and warm-start successively larger orders up to an explicit maximum.
  Per-equation residual tolerances and per-state waveform-change tolerances
  have the corresponding physical units; frequency-change tolerance is relative.
  `HbRefinementCheck` records each solved order, assessed band and sample
  count, total/retained/omitted/sampling ratios, and actual inter-order changes.
- `ShootBudget` / `solve_shooting` -> `ShootOrbit` (period, orbit
  point, Floquet multipliers from the forward-differenced monodromy
  through the fs-la complex eigensolver). Admission requires the explicit
  time-domain callback; a frequency-only problem returns
  `TimeDomainUnavailable` before evaluating its island.
- `ContinuableProblem` / `ContinuationBudget` / `continue_branch` —
  TRUE pseudo-arclength (bordered corrector on coefficients and the
  parameter jointly), so parameter folds are traversed, not stalled
  at. This original entry is forced-anchor only.
- `autonomous::continue_autonomous_branch` adds complete autonomous
  coefficient/frequency/parameter continuation with a full-waveform integral
  phase condition, independent physical scales and explicit work/payload caps.
  It returns `AutonomousBranch`, with exactly two corrected seeds plus the
  declared number of accepted additional points, not a parameter endpoint.

## Autonomous cycle continuation

`AutonomousScales` separately conditions original states, balance rows,
frequency and parameter. State coefficients use Parseval weights, so exact-zero
harmonic padding does not change the arclength metric. Large DC offsets do not
relax the absolute scaled residual tolerance or non-DC oscillation floor.

The corrector retains every original HB row and adds an integral phase row
from all non-DC state harmonics, plus a fixed-parameter seed row or arclength
plane. It works when the first state has no oscillatory component. Frequency
and physical parameter are independent unknowns. A phase-projected secant
predicts later points; central differences, native LU, residual-decreasing
backtracking and bounded step halving/growth correct them. No damping,
equilibrium deflation or diagonal jitter is introduced on this path.

Original AFT buffers and trigonometric tables are reused within each corrector.
The existing independent odd-grid sampling test applies to every accepted
point. Its `HbOrbit` retains physical coefficients/frequency and unscaled HB
residual/trace/defect; separate fields retain the combined dimensionless
residual, phase/arclength rows and scaled sampling discrepancy.

Dimensions, full retained output/AFT/dense payload and total residual calls
are bounded before the corresponding work. Cancellation is polled between
residuals, Jacobian columns and steps; one native residual or LU remains atomic.
On refusal the mutable family parameter returns to the last accepted value,
or the first seed parameter before acceptance. No partial branch is published.
Callbacks must be deterministic and may not conceal evolving state.

The authored `tests/autonomous_continuation.rs` checks an analytical fold of
cycles with a static first coordinate, amplitude-dependent frequency, common
phase shifts, independent state/equation unit changes, zero harmonic padding,
pre-source resource refusals, equilibrium refusal and parameter restoration.
Its analytical family is lambda=(r^2-1)^2, omega=1+0.2*r^2. These are finite
numerical regression claims, not existence, primitive-period, branch-identity,
stability, fold-location or infinite-tail certificates. Extra-neutral
conservative families and hybrid cycles remain outside this entry's claim.

## Invariants

- The original `solve_hb` / `solve_hb_seeded` and continuation APIs use
  fixed, disclosed harmonic truncation `N` (X-Struct).
  AFT starts with `m = (4N+4).next_power_of_two()` samples and uses
  prepared det-routed trigonometric sums. Basis tables and state/RHS
  buffers are reused through Newton and finite-difference Jacobian
  evaluations; the DC port is evaluated once per residual.
- A converged candidate is re-evaluated on a finer `2m+1` grid.
  Both its residual and the norm of the change in retained balance
  coefficients must meet the same tolerance and coefficient scaling
  (explicit absolute physical normalization for the autonomous branch above).
  Otherwise Newton resumes on the finer grid. There are at most four
  refinements; the fourth is a verification grid. Persistent disagreement
  returns `SamplingUnresolved`, and refinements also consume the existing
  Newton iteration budget.
- This sampling gate also applies to every continuation corrector.
  Accepted points retain their actual residual trace and final verified
  residual, never a fabricated zero. The original forced corrector includes
  arclength in that norm; the autonomous result separates physical HB norms
  from its normalized phase/arclength diagnostics.
  `aft_samples` records the final verification grid and
  `sampling_defect` records the measured retained-coefficient change.
- The adaptive entry points grow solved orders as `min(2N+1, max_harmonics)`.
  Each candidate is independently assessed through harmonic `2N+1`, including
  harmonics above the maximum solved order when the cap is reached. Omitted
  state coefficients are exactly zero, so their equation defect is `-G_n`;
  no undefined or unnecessary higher-frequency port is evaluated.
  Assessment quadrature uses the same bounded finer-grid policy, starting no
  coarser than the candidate's verified sampling.
- Adaptive acceptance requires each full measured per-equation RMS defect,
  combining retained and omitted bands, to fit its single physical tolerance.
  Separate band diagnostics do not each receive the whole error budget.
  Quadrature disagreements must also fit the physical tolerances,
  plus waveform and frequency agreement between two actually solved orders.
  Parseval weights are one for DC and two for positive-frequency coefficients.
  Scaling is per row/state, so a large DC bias or another coordinate's units
  cannot dilute a small but relevant harmonic. The maximum must leave room
  for at least one larger solve. A retained-row failure refuses immediately
  rather than increasing order to compensate for a loose Newton tolerance.
- Warm starts retain every solved coefficient and move the solved frequency
  to its new canonical slot. Forced and autonomous phase gauges stay unchanged;
  comparison of conservative backbone shapes aligns only their free phase
  using the explicitly pinned first-state fundamental. Frequencies and
  unknown-frequency seeds must remain positive and finite.
- Canonical unknown ordering: DC block, then per-harmonic
  `Re`/`Im` blocks, state-major; the packing is the seed format.
- The explicitly scaled entry uses one row/column unit for DC and that unit
  divided by sqrt(2) for each positive-frequency real/imaginary coefficient.
  Thus Euclidean norms in the conditioned chart are real-waveform RMS norms.
  Equation units are independent of state units, as required by general
  mechanical/acoustic/electrical frequency ports. Newton differences use each
  coordinate's declared physical magnitude; the dense solve is conditioned on
  both sides, and physical steps are restored before original-model evaluation.
  The backbone's extra row is its relative fundamental-amplitude error. No
  amplitude or frequency is silently normalized in the returned waveform.
  Invalid scales, physical seeds and budgets refuse before island/port calls.
- Backbone acceptance independently bounds `|2|X_1[0]|/a - 1|` by the
  requested tolerance. Other state coordinates and DC offsets cannot relax
  this relative amplitude gate. The appended Newton row is the unsquared
  norm difference `|X_1[0]|-a/2`, with its analytic norm gradient; this avoids
  squared-unit scaling and differencing a norm with a perturbation larger
  than the physical oscillation. The dynamical residual's coefficient scale
  uses only the non-DC coefficients, with the prescribed amplitude as its
  floor instead of one SI unit, so static preloads cannot hide an incorrect
  frequency or unresolved oscillatory balance.
  Backbone state-column differences use each coordinate's retained
  oscillation scale (the target supplies a zero-coordinate fallback);
  the frequency-column and non-backbone differencing rules are unchanged.
  Deflation scales only the physical equations, consistently in the Jacobian,
  right-hand side and line-search merit; it does not scale the amplitude row.
- A converged residual below the declared tolerance is the ONLY
  accepted discrete-balance evidence. Unanchored autonomous HB solves
  refuse equilibrium by name (`TrivialCollapse`), checking the full
  non-DC spectrum rather than the fundamental alone. A backbone's independently
  resolved prescribed positive amplitude excludes equilibrium without an
  arbitrary absolute displacement floor. The original
  absolute autonomous HB collapse threshold is 1e-9 in packed state units;
  the new branch uses its explicit scaled RMS `minimum_oscillation`.
  Forced HB solves admit DC, zero, arbitrarily small, and higher-harmonic
  responses: their period is prescribed, not inferred from a cycle.
  Equilibrium deflation applies only to the original autonomous/backbone anchors.
  Shooting applies its derivative-based equilibrium check only to
  autonomous problems: a forced response can have zero velocity at an
  ordinary turning point.
- Every RK4 stage, shooting Jacobian, monodromy evaluation, and autonomous
  equilibrium check uses the admitted time-domain callback. Shooting never
  substitutes `island` for a general linear port's realized dynamics.
- In the original HB solver, conservative families carry a parity/phase
  nullspace: its Newton step regularizes a singular factorization with FIXED
  Tikhonov jitter; the converged residual gate is unchanged. Autonomous
  parameter continuation does not use that regularization.

## Error model

Typed `OrbitError`: `BadParameter`, `NewtonStalled` (with the full
per-iteration residual trace), `SamplingUnresolved` (finest sample count,
residual and measured defect), `HarmonicsUnresolved` (last finite-band check,
including whether retained equations or omitted/order-change evidence failed),
`SingularJacobian`,
`ContinuationExhausted`, `TrivialCollapse`, `TimeDomainUnavailable`,
`TorusSuspected` (a
non-trivial unit-circle Floquet pair with nonzero angle:
quasi-periodic drift is a NAMED no-claim in v1, never a wrong
answer), forwarded `Eigen`.

`autonomous::AutonomousContinuationError` adds typed argument, cancellation
and residual-evaluation-limit refusals while retaining native numerical errors.

## Determinism class

Deterministic: bounded deterministic grid refinement, fixed iteration caps and step policies,
det-routed transcendentals, no RNG, no time dependence; bitwise
repeatability is a conformance case.

## Cancellation behavior

Synchronous computation; no direct `Cx` dependency. Original HB/shooting and
forced continuation do not poll cancellation. The new autonomous entry accepts
an explicit cancellation callback; an adapter may supply its original `Cx`
checkpoint. Native AFT residuals and dense LU factorizations remain atomic.

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
incommensurate two-center fixture; bitwise determinism. Additional G1/G3
regressions refuse a frequency-only port before island evaluation, compare
HB and shooting for `2 x' = -x + cos(t)` against the exact phase-zero value
0.2 and Floquet multiplier `exp(-pi)`, and retain a forced cosine response
whose initial derivative vanishes.

The G1/G3 analytic regressions in `tests/harmonic_balance.rs` compare
higher-harmonic forced oscillator coefficients (including their phase and
velocity) with the closed-form transfer function, exercise small and DC
forced responses, and retain autonomous-equilibrium refusal. A seventh-order
nonlinear spring exposes the old eight-point alias (0.9968112456 amplitude
instead of the retained Galerkin answer 1.0); an independent shifted
65-point quadrature checks the corrected coefficients. Further regressions
verify truthful continuation residuals against the analytic scalar balance
and a bounded named refusal for persistently unresolved forcing.

The physically scaled path additionally recovers the exact fundamental and
generated third harmonic of a nonlinear cubic cascade under independent state,
equation and clock-unit changes, including picometre-scale coordinates. A
1e12 unrelated static coordinate cannot hide its missing third harmonic; an
insufficient Newton budget refuses. Other G1/G3 cases retain prescribed
backbone amplitude/frequency under free phase, refine the seventh-order alias
in small physical units, refuse invalid scales before source calls, and replay
the accepted coefficients and scaled residual bit-for-bit on the same toolchain.

Backbone regressions additionally correct a seed with twice the prescribed
micrometre or picometre displacement, under phase shifts and optional
deflation, without rejecting the resolved small motion as equilibrium. A
large unrelated static coordinate cannot hide an unresolved amplitude or an
incorrect frequency; adequate iterations recover both.
A nonlinear Duffing coordinate-rescaling test compares the complete retained
waveform and frequency after changing displacement scale by 1e-6; invalid
zero, negative, nonfinite and unrepresentable-half amplitudes, and relative
tolerances outside `(0,1)`, refuse before model callbacks. These verify the
discrete amplitude constraint and scale
handling, not continuous-orbit existence or a global nonlinear convergence
guarantee.

Adaptive-order regressions use the exactly soluble nonlinear cascade
`x' = -x + cos(t), y' = -y + x^3`: the third harmonic is absent from an
N=1 solution even when its retained balance is converged. The adaptive solve
checks orders 1, 3, and 7 against exact complex coefficients, and requires a
second actual solve after the missing harmonic enters. A capped N=2 solve
detects the omitted third harmonic without querying ports outside its domain;
rescaling an equation and adding DC bias do not mask that defect. Further
cases verify backbone-frequency warm starts, retained-residual refusal under
a deliberately loose Newton tolerance, admission before model callbacks, and
refusal when separately sub-budget retained/omitted bands exceed the single
combined physical RMS budget.

## No-claim boundaries

- Quasi-periodic orbits (tori) are a NAMED refusal, not a claim.
- The original `continue_branch` remains forced-anchor. The additive autonomous
  entry solves frequency independently but is fixed-order and does not certify
  folds, stability, primitive periods, global branches or extra-neutral families.
- Collocation is described by the I09 slot but not implemented here;
  shooting is the second method.
- No a-priori truncation-error bound: the returned `N` remains disclosed
  structure. Fixed-order sampling refinement improves quadrature; the
  optional adaptive-order entry points also compare finite omitted bands
  and successive solved waveforms. Neither supplies a rigorous infinite-tail
  bound, an existence proof, or a conditioning-dependent forward-error bound.
  Agreement between finite grids is a measured convergence check, not a
  universal aliasing bound: sufficiently high frequencies can alias on both
  grids or lie outside the checked spectral band. Explicit input bandwidth
  must be admitted by domain adapters. Continuation remains fixed-order;
  adaptive order does not add relative-periodic, hybrid, matrix-free, or
  resume capabilities. Shooting is an independent cross-check
  only for problems with an equivalent time-domain
  formulation. The explicit callback is a model-supplied realization, not
  a proof of that equivalence. Frequency responses are never automatically
  inverted or assigned a state-space realization by the shooting solver.
- The monodromy is forward-differenced, not a consistent tangent;
  the I09 lattice's consistent-tangent element is a successor.
- Grazing/hybrid (non-smooth) islands are out of scope in v1.
