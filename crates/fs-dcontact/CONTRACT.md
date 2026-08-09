# fs-dcontact — Contract

## Purpose and layer

Layer L3. Distributed unilateral contact (bead
frankensim-fsim-distributed-contact-q6nmy): one-sided power-law
penalty potentials over collocation lines/profiles — string-fretboard
rattle, sitar/tanpura jawari bridges, snare wires, reed lays — as
fs-phs storage, so the Gonzalez discrete-gradient stepper makes
collisions energy-exact by construction (the Bilbao-Chatziioannou
energy-consistent collision doctrine through the pHS machinery; no
LCP, no bespoke scheme).

## Public types and semantics

- `Obstacle` — collocation matrix `Phi[i][k]`, per-point gaps and
  quadrature weights, shared `(K, alpha)` power law with a PROVENANCE
  string (logged, never invented; matdb lookup deferred until packs
  carrying contact-law parameters exist — no fake wiring). Admission
  refuses NaN/non-finite entries, negative weights/stiffness, and
  `alpha < 1` (the force law loses C^1 there) by name.
- `ContactStorage` — wraps any fs-phs `Storage` with
  `Phi_c = sum_i w_i K/(alpha+1) [p_i]_+^(alpha+1)`,
  `p_i = (Phi q)_i - c_i`; exact analytic gradient; `probe()` reports
  active points, max penetration (the authored-ceiling FLAG for
  stiffness inadequacy — deliberately not a refusal), and contact
  energy.
- `polyline_heights` — linear-interpolated obstacle profiles with
  typed refusals (unsorted, short, out-of-span).
- `string_collocation` — mass-normalized sine-mode collocation for
  fixed-fixed strings.
- Stepping, the iteration budget, and the disclosed solver residual
  are fs-phs's (`step`, `StepRecord::{newton_iters,
  solver_residual}`, `NewtonStalled`): ONE frozen solver law.

## Invariants

1. Polyline gap geometry matches analytic interpolation to 1e-14.
2. Bouncing analytic fixture: pure-potential restitution is 1 within
   1e-3, and max penetration matches the closed form
   `E = w K/(alpha+1) p^(alpha+1)` within 2% plus the bounded gravity
   correction.
3. String-fret rattle with >= 3 contact events conserves H to 1e-8
   relative over 20k steps; the explicit-integrator mutation
   (symplectic Euler on the same field) visibly grows energy by
   >= 1e3x the discrete-gradient drift or diverges.
4. Dropped one-sidedness is detected as ATTRACTION: a separated state
   feels zero force from the true storage and a large spurious force
   from the two-sided mutant.
5. Iteration budget holds across a 3-level velocity sweep (max
   newton_iters <= 50, histogram logged); the stall path is the TYPED
   `NewtonStalled` (executed at absurd stiffness x step size).
6. Collocation refinement converges (8 -> 32 -> 128 point contact
   energy differences contract).
7. Jawari casebook: a graded bridge profile grazes the swinging
   string and pumps >= 5% (measured 75%) of the energy into the high
   band, against the CLEAN-TERMINATION control (the plain string,
   which preserves the pluck's band split exactly). Band energies
   logged as JSON. Executed design lesson: a tiny-gap rattling point
   is itself a collision exciter, not a "hard bridge" control.
8. Determinism: repeated runs are bitwise identical.

## Error model

Typed `DContactError` (shape, parameter); solver refusals pass
through as fs-phs `PhsError`. No silent degradation.

## Determinism class

Deterministic: fixed collocation, fs-phs stepping, no RNG or time.

## Cancellation behavior

Synchronous; cost is the implicit step (fs-phs Newton with contact
scan per gradient evaluation). No `Cx` integration
(workspace `frankensim-ccmn`).

## Unsafe boundary

`#![deny(unsafe_code)]` via workspace lints; no capsules.

## Feature flags

None.

## Conformance tests

`tests/dcontact.rs` (9): polyline geometry + refusals; bouncing
restitution/penetration closed form; fret-rattle conservation +
explicit mutation; one-sidedness mutation; iteration budget sweep +
typed stall; collocation refinement; jawari vs clean termination
(band-energy JSON); bitwise determinism; typed parameter refusals.

## No-claim boundaries

- NORMAL contact only: tangential friction (bowing) is its own future
  bead (stated in the bead's polish round).
- No contact-internal viscous loss (Hunt-Crossley damping inside the
  collision): the bare potential's restitution is exactly 1 and
  losses enter through the modal damping `R`; the lossy contact law
  is a recorded follow-up (it needs state-dependent `R`, outside the
  constant-R pHS form).
- `alpha < 2` laws carry an unbounded contact Hessian at the boundary
  (`d2f ~ p^(alpha-2)`); the tight-graze regime stalled the
  FD-Jacobian Newton at `alpha = 1.5` (typed refusal, executed) —
  stiff grazing fixtures should use `alpha >= 2`, recorded as
  guidance, and a globalized Newton is a potential fs-phs follow-up.
- Contact-law `(K, alpha)` materials data: no packs exist yet; the
  provenance field is the wiring point when they do.
- Obstacles are static (no moving frets/bridges) in v1.
