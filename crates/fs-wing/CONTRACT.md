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
  { gamma, condition_est, surface_lift_n, total_lift_n }`.

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
