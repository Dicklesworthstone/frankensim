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
  physical quadrature measures, shared `(K, alpha)` power law with a PROVENANCE
  string (logged, never invented; receipted construction via
  `Obstacle::from_receipt` from fs-matdb contact packs). Fields are
  PRIVATE so admission is mandatory (accessors + a documented
  `from_raw_parts` trust escape for mutation batteries). Admission
  refuses NaN/non-finite entries, negative weights/stiffness, and
  `alpha < 1` by name. Smoothness, stated precisely: the POTENTIAL is
  C^1 for all `alpha >= 1` (what Gonzalez needs); the force is C^1
  only for `alpha > 1`, and the contact Hessian is unbounded for
  `1 < alpha < 2` and discontinuous at `alpha` in {1, 2}.
  A point has measure 1; line and area sites use m and m2. `K` is per unit
  measure, so a line law has units N/m^(alpha+1). Measures are never normalized
  to a probability distribution. In particular, splitting a uniform physical
  line into 8/16/32 cells preserves its total force within 1e-12 relative.
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
   >= 1e3x the discrete-gradient drift or diverges. HONEST SCOPE
   (review-corrected): under Gonzalez stepping the correction term
   forces dg.dx = dH for ANY gradient, so conservation certifies the
   INTEGRATOR, not the contact gradient — the gradient's own oracle
   is invariant 10.
4. Dropped one-sidedness is detected as ATTRACTION: a separated state
   feels zero force from the true storage and a large spurious force
   from the two-sided mutant.
5. Iteration budget: an INTERIOR bound (max newton_iters <= 24,
   measured <= 8) holds across a 3-level velocity sweep with the
   histogram logged — asserting the fs-phs cap itself would be
   vacuous, since a successful step cannot exceed it
   (review-corrected); the stall path is the TYPED `NewtonStalled`
   (executed at absurd stiffness x step size).
6. Collocation refinement converges (8 -> 32 -> 128 point contact
   energy differences contract).
7. Jawari casebook: a graded bridge profile grazes the swinging
   string and pumps >= 5% (measured 75%) of the energy into the high
   band, against the CLEAN-TERMINATION control (the plain string,
   which preserves the pluck's band split exactly). Band energies
   logged as JSON. Executed design lesson: a tiny-gap rattling point
   is itself a collision exciter, not a "hard bridge" control.
8. Determinism: repeated runs are bitwise identical.
9. `polyline_heights` and `string_collocation` refuse non-finite and
   non-physical inputs by name (NaN slipped ordering comparisons —
   executed).
10. THE contact-gradient oracle: central finite differences of the
    coded `H` match `gradient()` to 1e-5 relative on a multi-point,
    non-uniform-weight obstacle with a MIXED active/inactive contact
    set, and `probe().contact_energy` equals the Hamiltonian split
    exactly.

## Error model

Typed `DContactError` (shape, parameter); solver refusals pass
through as fs-phs `PhsError`. No silent degradation.
Both obstacle and storage `dissipative_modal_forces` now return `Result`.
Velocity/collocation/state mismatches refuse even when chi is zero, and
nonfinite state or derived forces refuse before a force vector is returned.
`Constitutive { owner, field }` preserves a higher-level constitutive refusal
without making this crate depend on its owner. The loss formula and local
non-attractive unloading saturation are unchanged.

## Contact parameter derivatives

`Obstacle::resistance_parameter_vjp_into` differentiates the original total
positive resisting contact force with respect to `[K, chi, gap_0, ...]`, with
exponent, collocation and physical measures fixed. It contracts a generalized
force seed without allocating or dividing by coefficients, so `K=0` and `chi=0`
retain their admitted-side sensitivities. Open contact and strictly capped
unloading select zero derivatives; cap equality selects the existing
unsaturated one-sided derivative. Shape, malformed-owner and nonfinite-input
refusals leave output unchanged; arithmetic refusal invalidates scratch.
`tests/contact_parameters.rs` checks these local derivatives against the
original storage/loss forces, including nonuniform multipoint measures and the
inactive/capped branches. This is a local constitutive parameter derivative,
not a source-card derivative or an identification/uniqueness claim.

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

M4.1-W1 adds G0 physical-measure refinement and velocity-shape refusal;
G1 Hertz contact duration against Johnson section 11.4 on three timestep rungs
(finest relative error <=0.5%, observed order >=1.8); and quantitative
Hunt--Crossley restitution against its exact logarithmic relation for
chi*v in {0.01,0.1,0.5} (absolute error <=1e-3). The latter uses the original
implicit dissipative port, never a lagged explicit force. These are named
test obligations pending central native execution, not physical validation.

The HC chi*v=0.5 fine rung exposed a numerical two-cycle in the original
Gonzalez numerator `H(b)-H(a)-gradH(mid).(b-a)`: cancellation of O(H) endpoint
energies contaminated its much smaller cubic remainder, then the nonlinear
loss's effort derivative amplified that noise in the residual. ContactStorage
now supplies `Storage::midpoint_energy_remainder` when its inner storage borrows
its current exact quadratic matrix or supplies its own complete stable remainder.
No quadratic is inferred from sampled gradients, a cached classification or an
inactive contact. The TOTAL unilateral storage remains nonlinear and does not
acquire the exact-quadratic Newton-floor allowance.

Active power-law remainders use a centered binomial series only where every
successive term contracts by at most 1/4; through degree 33 its remaining tail is
below 2^-63 of the leading cubic term. This evaluates the same power-law energy
remainder to roundoff, including affine-collocation rounding corrections. Entry,
separation and large increments retain the original unilateral potential.
Unknown inner storage retains the original endpoint-subtraction path. Both the
canonical residual and analytic Newton use the same numerator; solver requests,
floors, update limits, conservative/loss laws and physical gates are unchanged.
This may change contact trajectories' final rounding bits, not their physical
parameters; existing music replay/baseline checks remain required centrally.

Focused G0 tests compare Hertz remainders with independent double-double
sqrt/multiply at three increment scales, check unilateral entry/separation, and
replay the central failed HC pre-step state under its original 1e-13 relative,
zero absolute, 50-update request. An opaque-inner mutation retains the old
refusal; the stable source must converge below that SAME request without a
raised floor and retain energy closure. Native verification remains pending.

`tests/dcontact.rs` (10): polyline geometry + refusals; bouncing
restitution/penetration closed form; fret-rattle conservation +
explicit mutation; one-sidedness mutation (force swept over separated
depths); iteration budget sweep + typed stall; collocation refinement
(0.1 gate, measured 0.0099); FD gradient oracle + probe consistency;
jawari vs clean termination (band-energy JSON); bitwise determinism;
typed parameter refusals.

## No-claim boundaries

- NORMAL contact only: tangential friction (bowing) is its own future
  bead (stated in the bead's polish round).
- Hunt–Crossley internal loss is a **dissipative port force**
  `χ K [p]_+^α ṗ` ([`ContactStorage::dissipative_modal_forces`]), not
  a term in `H` and not a state-dependent `R`. `χ = 0` (the default)
  keeps elastic restitution 1. Tangential friction is still
  `fs-tribo`, not this crate.
- `1 < alpha < 2` laws carry an unbounded contact Hessian at the
  boundary (at `alpha = 1` the Hessian coefficient vanishes and the
  FORCE has a C^0 kink instead — review-corrected); the tight-graze
  regime stalled the FD-Jacobian Newton at `alpha = 1.5` (typed
  refusal, executed) — stiff grazing fixtures should use
  `alpha >= 2`, and a globalized Newton is a potential fs-phs
  follow-up.
- Contact-law `(K, alpha)` materials data: WIRED (music bead
  3ez8g.13.1). fs-matdb's `ContactLawCard`/`ContactReceipt` carry
  the pair + geometry context, identification, force/velocity
  validity, and the graze advisory; `Obstacle::from_receipt` builds
  from the receipt verbatim and formats the provenance from the pack
  identity. Absent packs refuse upstream BY NAME (no default K
  anywhere on the path); out-of-validity lookups refuse in fs-matdb
  rather than extrapolate. The jawari fixture's migration to the
  receipted card is proven PROVENANCE-ONLY (bitwise-identical
  2000-step trajectory, `contact_pack_migration_is_provenance_only`).
- Obstacles are static (no moving frets/bridges) in v1.
