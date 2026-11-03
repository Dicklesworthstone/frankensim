# fs-wing — CONTRACT

Bead: frankensim-wf-root-guzez.5.3.1 (E4.2-i, Wright Flyer program).
Spec: COMPREHENSIVE_PLAN §5.2 (ROUND 6 steady state). Conventions:
frame-conventions-v1 (frd body axes).

## Purpose and layer

L3. Coupled multisurface lifting-surface aerodynamics. E4.2-i ships the
WeissingerLLinear machinery: multisurface horseshoe-panel layout (both
wings, both canard planes, verticals; ≥2 chordwise rows where hinge
moments matter), influence assembly at the three-quarter-chord control
points, one deterministic dense LU with a condition estimate per solve,
Kutta–Joukowsky strip forces per surface. E4.2-ii adds the nonlinear
section closure (fs-airfoil datasets vs induced alpha, warm-started
safeguarded iteration, branch identity).

## Public types and semantics

- `Panel { surface, a, b, ctrl, normal, width }` — bound vortex a→b at the
  quarter chord, control point at three-quarter chord, unit normal.
- `flat_surface(...)` — rectangular n_span×n_chord layout builder.
- `solve_weissinger_linear(panels, freestream, rho) -> SolveReport
  { gamma, condition_est, surface_lift_n, total_lift_n, total_force_n }`.

- `images::CertifiedGround` + `solve_weissinger_ground` — the exact
  vortex-image system across the certified flat plane (mirrored corners,
  a/b swapped, mirrored-stream trailing legs); admissible ONLY inside the
  FlatnessCertificate bands (cross-linked with fs-flyer::prelaunch;
  `ground-uncertified` otherwise); `aircraft-below-ground` refusal.
- `nonlinear::solve_nonlinear` — the production force path: strip-level
  section closure vs TRAILING-SHEET induced alpha, warm-started
  safeguarded relaxation, branch identity, enforced operator-reuse rule.

### Complete circulatory force vectors

`SolveReport::total_force_n` and `NonlinearReport::total_force_n` expose all
three body-frame components [N] of the existing Kutta–Joukowsky reduction,
`sum(rho * gamma_j * (velocity_j cross (b_j - a_j)))`. Linear free-air,
certified-ground and nonlinear solves accumulate the streamwise and lateral
components alongside the existing lift reduction. `total_lift_n` and
`surface_lift_n` retain their negative body-z meaning; they are not wind-axis
lift magnitudes. The new vector's z component is assigned directly from
`-total_lift_n`, preserving that exact projection and the existing scalar
multiplication and accumulation order.

This makes vertical-fin sideforce and canted-wing loads available to linear
solve consumers without reconstructing circulation forces. It adds no
allocation or extra influence solve. The vector includes only the represented
circulatory force; profile drag, thrust and weight belong to their existing
owners, and this report adds no separate induced-drag model. Each solve retains
its existing choice of incident velocity, wash and circulation closure.

The nonlinear section equation still measures incidence and trailing-sheet
induction in the global body xz plane. Reporting its full force vector does
not make that section closure objective for arbitrarily oriented or vertical
lifting surfaces. The linear rotation fixtures do not assert such a nonlinear
model extension. Ground-image rotation invariance is asserted only for
rotations preserving the certified horizontal plane.

`tests/force_vector_battery.rs` supplies physical regressions for signed
vertical/canted linear loading, rigid-rotation covariance of circulation and
force, force perpendicularity to uniform incident flow, ground-plane-preserving
rotation, and the analytic `q*S*cl` force projection under explicit local flow
and wash for a constant section coefficient of either sign.

### Circulatory spatial wrench and local-flow hinge loads

`wrench::CirculatoryWrenchAccumulator` streams accepted panel circulation and
the corresponding local incident velocity into an `fs_ga::Wrench`. Density
and moment origin are fixed at construction; all geometry, velocities and the
origin share one Cartesian frame. Each panel contributes
`F = rho * gamma * (V cross (B - A))` at the bound midpoint and
`M = (midpoint - origin) cross F`. The full bound segment retains swept,
canted and vertical loads and all three moments. Uniform or per-strip flow
uses the same API, without allocation on success, ownership reconstruction or
an additional aerodynamic solve. Insertion order and native Vec3 arithmetic
are retained.

Malformed geometry, nonfinite circulation/flow and unrepresentable panel or
accumulated wrenches refuse with `circulatory-wrench-invalid`; a failed add
leaves the accumulator unchanged. These are circulatory loads only: profile
and induced-drag models, section couples and apparent mass retain their own
owners. This force assembly does not extend the nonlinear section-plane model.

`hinge::hinge_load_with_incident` accepts one local velocity per panel in the
full geometry, retaining the coupled solution's wash, gust and rigid-motion
effects on selected hinge surfaces. It projects each panel's checked wrench
onto the hinge axis, then adds the existing section couples as free vectors.
The uniform entry point shares this owner while retaining its original scalar
projection/sum order and frozen golden.

Three wrench tests establish per-panel zero incident-flow work, six-component
bit parity with the previous native propulsion loop, shifted-origin moments,
rigid-frame/coadjoint transport and atomic refusal. The seven-case hinge
battery includes nonlinear local-wash/gust loading, an independent panel
quadrature, a live uniform-flow counterexample, origin shifts and uniform
bit equality. These and the unchanged nonlinear/continuation/oscillation
batteries pass in a 29-test native numerical-source run; this is not a full
workspace or frontend qualification.

### Spatially varying strip inflow

`nonlinear::solve_nonlinear_with_incident` accepts one body-frame incident
relative-air velocity per strip. Axial wash adds to the local x component before
incidence and chord-plane speed are computed; the same resulting local velocity
enters circulation and Kutta–Joukowsky lift. Spanwise flow does not increase the
2-D section speed. Nonfinite or nonpositive chordwise flow and incomplete strip
ownership refuse. The reference freestream still defines common straight
trailing legs and linear chordwise row ratios: this is quasi-steady strip gust
loading, without chord-resolved unsteady gust memory.

Both nonlinear entry points cache immutable strip/panel Biot–Savart coefficients
once per solve and reuse iteration scratch, retaining the original product/sum
order. `nonlinear_battery` adds uniform exact-equivalence, mirrored nonuniform
loading/local lift, and wash-vector/refusal tests. Existing goldens are unchanged.

`NonlinearReport::alpha_eff_rad` retains each strip's effective incidence at the
accepted circulation. It reuses the prepared trailing-sheet coefficients once
after convergence, with the same local flow and geometric twist as the lift
closure. Section drag/moment consumers can therefore use the accepted incidence
without a second solve, repeated Biot–Savart construction, or ambiguous inversion
of a stalled lift coefficient. G0 checks invert a linear section law from the
returned circulation under nonuniform inflow and twist and verify that the
reported incidence includes a nonzero induction correction.

### Noncontracting alternating strip modes

The nonlinear owner also detects two successive residual-vector reversals
(angle greater than 120 degrees) without contraction in the Euclidean norm of
the unrelaxed circulation mismatch. Shared scaling makes this comparison finite;
an allowance of 64 machine epsilons recognizes a roundoff-level plateau. This
allowance is exclusively a safeguard diagnostic and never changes the acceptance
tolerance. The scalar relative residual alone cannot identify this defect: a
two-cycle may alternate between large and small denominators, while a gradually
growing negative mode may stay below the existing 25% growth threshold.

On detection, the solver halves its current relaxation and reevaluates the same
state. It uses the existing `OMEGA_MIN=0.02`, `MAX_ITERATIONS=600` (including
repeated evaluations), `GAMMA_TOL=1e-10`, section closure and starting state.
There is no Newton/global branch search, extrapolation, cold restart or linear
fallback. The additional scratch is one residual vector; reductions are fixed
order. Ordinary successful updates retain their prior arithmetic and the
original nonlinear golden remains unchanged. Nonfinite circulation mismatch
and an exhausted relaxation floor return `nonlinear-did-not-converge`.

`tests/nonlinear_oscillation.rs` freezes two actual low-speed Huffman launch
wing evaluations exposed by the corrected actuator-disk wash. The unchanged
owner exhausted 600 evaluations at residuals 3.220064062700418 and
6.113375193275809. The safeguard resolves the same attached solutions in 422
and 28 evaluations. Independent reconstruction from the returned circulation,
the complete public horseshoe field and a separate finite-bound-segment
subtraction measures final relative mismatches `8.095e-13` and `4.729e-11`,
respectively, both below the unchanged tolerance. Replay and warm continuation
retain that branch. These are discrete nonlinear-solver checks, not evidence of
global uniqueness, aerodynamic validity outside the section model, or full
flight calibration.

### Circulation sign changes and divergence detection

The two-consecutive-growth guard compares the Euclidean norm of the dimensional
circulation mismatch, using the same overflow-safe scaled squared norms as the
oscillation guard. Its 25% growth threshold, relaxation bounds, iteration budget
and final componentwise relative tolerance are unchanged. The convergence
denominator `|Gamma|+1e-9` is state dependent: a strip approaching zero circulation
can increase the relative residual even while its actual equation mismatch
contracts. Such a crossing must not trigger progressively smaller relaxation.

The third `nonlinear_oscillation` fixture freezes the real distributed wing and
canard flow from assisted Historical(3) flight at tick 974. The old guard made
five unnecessary cuts to `omega=0.021875` and exhausted 600 evaluations at
`9.460544473899003e-7`. The corrected guard solves the same attached branch in
162 evaluations; independent public-field reconstruction measures final relative
mismatch `7.657e-11`, below the original `1e-10` gate. Signed opposing wing loads,
bitwise replay and warm continuation are checked. The existing nonlinear golden
and both low-speed alternating-mode fixtures remain unchanged. In the actual
native SimLoop, the assisted flight now reaches physical ground contact at tick
1102 (8.575 seconds airborne, three undulations); the fixed-control fixture's
outer prop/airframe coupling refusal at tick 553 remains explicit.

## Invariants

- No scalar biplane factor anywhere: interference EMERGES from the
  influence matrix (battery: Munk-class gap trend).
- Deterministic assembly + LU (fixed pivot tie rule): mirror-symmetric
  inputs give mirror-equal circulations to 1e-12 relative.
- Every solve reports a condition ESTIMATE (honestly labeled estimate).

## Error model

Typed `Refusal`: `panel-count-invalid` (cap AND cap+1), `panel-invalid`,
`freestream-invalid`, `influence-singular`, `influence-ill-conditioned`,
`layout-invalid`.

## Determinism class

Deterministic: det:: sqrt in norms, fixed elimination order, fixed pivot
ties; golden pinned under `org.frankensim.fs-wing.e42i-golden.v1`.

## Cancellation behavior

Synchronous pure functions; nothing to cancel.

## Unsafe boundary

Workspace `deny(unsafe_code)`; none.

## Feature flags

None.

## Conformance tests

`tests/weissinger_battery.rs`: monoplane CL_alpha within 12% of the
lifting-line trend at AR 6; the EMERGENT biplane gap effect (1903 gap
factor < 0.97, deficit shrinking with gap, approaching 1 from below);
bitwise mirror symmetry; canard-wing coupling live (sign + measurable
wing-loading change); refusal caps; pinned golden over the 1903
five-surface layout at the Dec-17 condition.

## No-claim boundaries

- WeissingerLLinear is an EXACT FIXTURE / admission-selected debug mode —
  never the production force path (E4.2-ii's nonlinear closure) and never
  an automatic fallback after nonlinear failure (plan law).
- Flat-plate panels here: camber/thickness enter through the E4.2-ii
  section closure, not panel geometry.
- No unsteady effects, no ground images (E4.4a), no wake models (E4.7):
  trailing legs are straight semi-infinite lines.
- The condition number is an ESTIMATE (single-probe Hager class).

## Explicit circulation continuation under changing local flow

`nonlinear::continuation::solve_with_incident` adds a caller-owned continuation
lane for the existing strip equations. `WingContinuation` stores accepted panel
circulation, effective strip incidence, layout identity and measured solver/domain
diagnostics. It is numerical history, not air inertia, unsteady circulation physics
or a delayed load. The existing cold/Picard APIs and their arithmetic remain in
place; malformed raw warm vectors now refuse before indexing. The existing
nonlinear golden remains `2797bcda52ee7d52ff0cb8283c640fe463d1f5996bd4c3e4a5d02162823e1fe0`.

For strip total `G_s`, the exact residual is
`F_s = (G_s + 0.5*V_s*c_s*cl(alpha_s))/(abs(G_s)+1e-9)`, where `alpha_s` uses the
original trailing-leg induction, local incident velocity, vector wash and twist.
The same original linear row ratios scatter strip totals, and the original
Kutta--Joukowsky owner evaluates forces with the actual local section velocities.
The final accepted circulation is reevaluated and every residual component must
be finite and strictly below the original `1e-10` threshold. Solver convergence
alone is insufficient. The new lane reports regimes and incidence at that same
accepted circulation.

Newton--Krylov, FGMRES and Armijo backtracking are the existing `fs-solver` owners.
The adapter supplies a matrix-free Jacobian action: analytic row-scatter and
downwash derivatives, a declared centered numerical section slope, and the
derivative of `abs(G_s)+1e-9` itself. A finite-difference check differentiates the
complete normalized residual with signed circulation and signed row weights,
away from nonsmooth points. This numerical section tangent is not an exact
adjoint or a differentiability claim at a constitutive kink or zero circulation.
The source method is C. T. Kelley, *Iterative Methods for Linear and Nonlinear
Equations*, SIAM (1995), chapters
[6](https://doi.org/10.1137/1.9781611970944.ch6) and
[8](https://doi.org/10.1137/1.9781611970944.ch8). No unconditional convergence
guarantee follows from selecting this method.

`ContinuationControls` declares finite positive incidence and incidence-change
bounds, a signed row-sum condition cap, minimum incidence-probe distance from
straight trailing filaments in chord units, section-difference stencil and work
budgets. Absolute incidence must remain below pi/2. A section-slope stencil must
be representable at the current incidence and remain in that absolute domain.
Regime labels are recorded, not treated as mathematical branch boundaries: a
smooth attached-to-blended constitutive transition is explicitly tested and
accepted. Incidence continuity is a local admission rule, not proof of uniqueness
or physical validity of a post-stall root. The multiple-root and spanwise
oscillation risks are documented, for example, by Wu, Gopalarathnam and Kim,
[AIAA 2005-6169](https://ntrs.nasa.gov/api/citations/20050207439/downloads/20050207439.pdf?attachment=true);
this lane does not implement that paper's decambering model.

Signed chord-row ratios remain legal. A multirow strip refuses if
`sum(abs(g_linear))/abs(sum(g_linear))` exceeds its declared cap, is nonfinite,
or its total is at the original `1e-14 m^2/s` equal-share fallback floor. No ratio
clipping or equal-share fallback is admitted by this lane. Near-filament cases
refuse rather than altering the vortex kernel. Geometry/strip ownership changes
invalidate the anchor; current freestream, twist, incident flow and wash can
change, with the existing influence-operator rebuild rule still enforced.

Work is finite: at most 64 Newton attempts, 600 residual evaluations and 8192
Jacobian products may be declared, with a restart bounded by 32 strip unknowns
and four inner cycles per outer attempt. Section queries are bounded by
`n_strips*(residual_evaluations + 3*jacobian_products)`. Geometry coefficients and
scratch buffers are prepared once per solve. Nonfinite section sources and budget
exhaustion refuse. An out-of-domain line-search trial is a rejected trial that
native Armijo may backtrack from; it does not latch a terminal source failure.
An actual native-Armijo test exercises that retry and independently checks the
resulting scalar equation. No approximate answer or alternate cold root is
returned on failure, and a borrowed accepted anchor is unchanged.

Validation includes two private Jacobian/Armijo tests and four external
continuation tests: an independent dense section-equation and force oracle with
nonuniform flow and wash; a smooth actual stall-blend transition; stale-layout,
incidence, cancellation and filament admission; and bounded source/solver
failures, caps, malformed warm vectors and density admission. The existing
10-case nonlinear, two-case oscillation, six-case Weissinger and six-case hinge
batteries also pass. These are focused native runs of the actual wing source and
native solver numerics, not a claim that the complete workspace build passed.
