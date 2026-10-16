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
  string (logged, never invented; receipted construction via
  `Obstacle::from_receipt` from fs-matdb contact packs). Fields are
  PRIVATE so admission is mandatory (accessors + a documented
  `from_raw_parts` trust escape for mutation batteries). Admission
  refuses NaN/non-finite entries, negative weights/stiffness, and
  `alpha < 1` by name. Smoothness, stated precisely: the POTENTIAL is
  C^1 for all `alpha >= 1` (what Gonzalez needs); the force is C^1
  only for `alpha > 1`, and the contact Hessian is unbounded for
  `1 < alpha < 2` and discontinuous at `alpha` in {1, 2}.
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
