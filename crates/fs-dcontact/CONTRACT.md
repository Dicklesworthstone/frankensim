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
  only for `alpha > 1`. The potential Hessian (contact stiffness) jumps at
  `alpha = 1` and is continuous at entry for `alpha > 1`; its derivative is
  unbounded for `1 < alpha < 2` and discontinuous at `alpha = 2`.
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
- `ContactStorage::continuous_tangent` retains a declared inner continuous
  tangent through boxing and composes the exact elastic contact Hessian action
  without additional allocation. Momentum and internal-memory rows retain their original
  inner law. At alpha=1 and exactly zero penetration, a nonzero normal direction
  refuses because the two-sided continuous derivative does not exist; directions
  tangent to the surface remain valid. The existing `hessian_vector_with`
  Newton helper retains its inactive-side convention. Hunt--Crossley loss is
  still separately owned and is not included in this elastic tangent.
- `polyline_heights` — linear-interpolated obstacle profiles with
  typed refusals (unsorted, short, out-of-span).
- `string_collocation` — mass-normalized sine-mode collocation for
  fixed-fixed strings.
- Stepping, the iteration budget, and the disclosed solver residual
  are fs-phs's (`step`, `StepRecord::{newton_iters,
  solver_residual}`, `NewtonStalled`): ONE frozen solver law.

## Momentum-consistent distributed contact work

For an exactly symmetric declared quadratic inner storage, or an inner storage
with its own `ConstitutiveDiscreteGradient`, `ContactStorage` now supplies a
componentwise discrete gradient. Each original scalar contact potential uses
the existing stable `OpeningContactStep` secant and is pulled back through the
same collocation row as its analytic force. This keeps contact effort out of
momentum and unrelated internal coordinates. A row that annihilates a rigid
translation exchanges equal and opposite discrete forces, so an internal pair
preserves total momentum under the canonical pHS update. Its contact impulse
is independent of a common Galilean boost, subject to represented geometry.

The previous full-state Gonzalez correction could conserve total energy while
adding a net internal force proportional to center-of-mass displacement. The
new composition preserves the original Hamiltonian, physical quadrature,
unilateral entry/separation, symmetry in endpoints and second-order consistency.
Hunt--Crossley loss remains the separately owned dissipative port. An arbitrary
nonlinear inner storage without a declared rule, and a merely tolerated
asymmetric quadratic matrix, retain the original Gonzalez path and acquire no
new momentum claim.

An allocation-free endpoint Jacobian differentiates this actual composed rule.
It delegates an owned inner action or uses the exact quadratic half-Hessian,
then adds each contact's scalar secant derivative through its collocation row.
A centered series prevents cancellation at nearly coincident active endpoints;
the coincident zero-penetration alpha=1 kink uses the inactive Newton side.
An inner rule without an analytic action retains the original probe fallback.
Inner coordinate scales are retained for those probes.
`tests/contact_momentum.rs` checks equal/opposite force, unchanged kinetic and
tail efforts, physical work, Galilean/origin invariance, independent cubic
force/tangent identities, contact crossings, and inner-rule composition.

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

`opening::AffineContactStep` also retains each point's Hunt--Crossley loss
increment before adding it to the conservative reaction. Its heat is the
original capped loss times local opening velocity, never the subtraction of
two nearly equal total forces. Thus weak representable dissipation survives
even when adding its force to a much larger elastic reaction rounds away the
increment. Pointwise unloading remains non-attractive before mixed-sign
collocation projection. G1 tests check both velocity signs at weak loss and
simultaneous capped unloading/loading with cancelling generalized elastic forces.

## Restitution-based loss calibration

`hunt_crossley_loss_from_restitution(e, approach_speed_m_s)` and
`Obstacle::with_restitution_loss` convert a measured normal restitution and
positive reference speed into the existing fixed Hunt–Crossley `chi` [s/m].
They invert the exact isolated-impact relation
`x - ln(1+x) = -e*x - ln(1-e*x)`, `x = chi * approach_speed`, rather than a
weak-loss approximation. Mass, positive stiffness and the power-law exponent
cancel from this relation. The setup-time solve is allocation-free and uses a
factored series near the elastic limit to retain weak loss. No calibration
work is added to the force or timestep path.

`e=1` gives zero loss; `e=0` requires infinite loss and refuses, as do invalid
speeds and unrepresentable coefficients. The supplied restitution applies to
an initially uncompressed, isolated normal impact at the reference speed;
external loading, other speeds, coupled modes and simultaneous contacts need
not reproduce it. This does not impose a restitution impulse or change the
obstacle's caller-owned provenance. Tests compare the inverse with the exact
forward relation and exercise the original implicit contact port at several
restitutions, masses, stiffnesses and speeds.

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
(workspace `frankensim-ccmn`). The one exception is `ContactInterface`
(feature `component`): its reaction search polls the kernel's `Cx` once per
trial reaction, and a cancelled search refuses as
`DContactError::Reaction { what: "cancelled", .. }`.

## Unsafe boundary

`#![deny(unsafe_code)]` via workspace lints; no capsules.

## Feature flags

`component` (off by default): `component::ContactInterface`, a normal
contact between two bodies as a Component of the composition kernel. It adds
the optional dependencies `fs-ports` and `fs-exec`; nothing else in the crate
changes with it.

## Contact interface on the kernel's load line (feature `component`)

- Snapshots. `ContactInterface` has the state codec of fs-ports'
  `Component`, so a system with one can be written to a durable snapshot
  and resumed. The state is two words, the closure and the reaction. It is
  written, and read back, only as a state the interface holds: both values
  finite and the potential at the closure finite. The identity is the
  obstacle's own, then every budget of the solve, then which ends read an
  end input; the ports are in a system's own description.
  `Obstacle::identity_words` (in `lib.rs`, with or without the feature) is
  that obstacle identity: the point count; the count and the bits of the
  collocation rows, the gaps and the quadrature weights; the stiffness, the
  exponent and the internal loss; and the provenance as its length in bytes
  and its bytes, eight to a word, so that it cannot run into what follows.
  The provenance is a label and no part of the force law, and it is bound
  because a changed source is a changed declaration. The obstacle, the
  interface and the solve are each written out by a pattern with no rest,
  so a new field does not compile until it is bound. Three tests in
  `tests/component.rs` (claims written before the first run): the state is
  read back apart and in contact, and no words, one word, three words, a
  closure that is not finite or whose potential is not, and a reaction that
  is not finite are refused and not written; an obstacle's words are equal
  for a rebuild and different for another collocation row, a second point,
  a second mode, another gap, weight, stiffness or exponent, an internal
  loss, another provenance and a provenance that continues the first; the
  interface's words are equal for a rebuild and different for another
  obstacle, each of the five budgets and each of three choices of ends.
  Not tested here: a snapshot through a kernel (fs-cli's scene tests go
  through fs-system).

- `ContactInterface::new(obstacle, ContactSolve)` admits one one-point
  `Obstacle` on the unit opening map `[-1]` (the map fs-couple's modal impact
  owner admits) with a non-empty provenance and nonnegative Hunt–Crossley
  loss. Ports `a` and `b` are mechanical force–velocity power ports
  (`FlowIn`, width 1). Positive loss adds `heat` after them (`EffortIn`,
  width 1): positive temperature in [K], negative generated heat rate out
  [W]. Zero loss retains the original two-port interface and conservative
  arithmetic. Its state `ContactState` is the closure `x` [m] (end `a`
  towards end `b`) and the last reaction.
- A step moves the closure by `dt` times the ends' mean closure rate. The
  reaction is the `OpeningContactStep` secant of the potential between the
  start and end closures, so a lossless contact stores what the ends give
  it: `stored_change` is the potential's change and `numerical_loss` is the
  port work minus it, with `dissipated` zero.
- With loss, the existing `AffineContactStep` evaluates the original
  one-point obstacle at openings `-x0`, `-x1` and opening velocity
  `-(x1-x0)/dt`. The elastic secant is `E`, the damping coefficient is
  `D = chi E`, and the total non-attractive reaction is
  `E + max(D (x1-x0)/dt, -E)`. The native response retains the signed loss
  increment before adding it to `E`, so representable weak heat survives
  rounding in the total force. Its dissipated power leaves through `heat`;
  the receipt records `dt` times that power and entropy generation `Q/T`.
  Stored energy remains the original conservative potential. Numerical
  loss is settled mechanical work minus its change minus exported heat,
  including the residual and arithmetic rounding; physical loss is not
  silently counted as numerical loss. All derived heat, entropy and energy
  values are checked before writing a candidate or outputs.
- On the load line the end closure is `x0 + dt (released_a - released_b) -
  S R` with `S = dt (c_aa - c_ab - c_ba + c_bb)`. The reaction is found by
  the bracket bisection of fs-couple's `bracket_reaction`, with the same
  unfused closure arithmetic and the same acceptance tolerance
  `force_absolute_tolerance_n + force_relative_tolerance max(|R|, |law|)`.
  Negative relative compliance refuses. For nonnegative compliance, the
  elastic secant and nonnegative loss factor are monotone in final closure,
  so the reaction at zero trial reaction bounds the same scalar root.
  With zero compliance the step is `propose`'s, to the bit.
- `ContactSolve` carries the fields and admission of fs-couple's
  `ModalContactConfig`: `max_iterations` in `1..=128`, positive finite force
  and penetration ceilings and tolerances, relative tolerance below one. The
  search refuses by name when the bracket passes `maximum_force_n`, when the
  iterations run out and when the bisection runs out of representable
  midpoints. A step that ends above `maximum_penetration_m` refuses.
- `reading_ends([a, b])` adds signal inputs `a.end` and `b.end` for
  neighbours that tell the end of their flow; the law uses only the mean
  closure rate, so the end flows are checked finite and not otherwise used.

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

`tests/component.rs` (`required-features = ["component"]`): a mass on a
Hertz wall in an fs-system graph, against the closed forms
`p_max = ((alpha+1) m v^2 / (2K))^(1/(alpha+1))` and
`t_c = 2 I p_max / v` (`I` computed by quadrature in the test), at 400 and
800 steps per contact time. The exit speed is within the energy the
search's force tolerance can leave over the penetration; the peak lies
between the closed form and one sampling interval's depth below it; the
contact time is within `4/n^2` relative (a cap about four times the
measured coefficient, not a derived bound), and the observed order is in
`[1.8, 2.2]`. A linear law settles where its closed-form reaction on a
load line puts it; the kernel's conformance kit passes apart, entering and
inside the contact; end inputs join told neighbours and leave the step
unchanged; refusals by name (shape, budgets, both ceilings, iteration
exhaustion, lengths, NaN, zero `dt`).

Its `component/dissipation.rs` module covers the closed-form linear
Hunt–Crossley load-line root, weak heat that rounds away in total force,
capped non-attractive unloading, heat and end-input ordering, bitwise
zero-loss compatibility, and invalid temperatures or entropy overflow
refusing before publication. A two-body impact checks momentum,
restitution and conversion of lost mechanical energy into a connected
thermal owner's energy; checkpoint replay and cancellation/retry check
that generated heat is neither lost nor duplicated.

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
- `ContactInterface` (feature `component`) remains one-point and does not
  shorten a step to the instant of touch or reach distributed obstacles.
  Hunt–Crossley loss is fixed by its obstacle; the heat-port temperature
  determines entropy generation, without temperature-dependent contact
  softening, thermal expansion or inferred material validity. Its heat
  outlet requires a connected thermal owner or an explicit reservoir.
  Its tests do not compare it with fs-couple's modal impact owner; that
  comparison belongs to the drum's kernel graph.

## Shared unilateral power-potential remainder

`power_potential_midpoint_remainder` exposes the same cancellation-stable
power-law numerator already used by ContactStorage. Callers supply admitted
finite arguments, exponent at least one and nonnegative coefficient, and
retain ownership of their physical domain and affine work map. It evaluates
the original positive-part potential; it neither supplies a force law nor
advances history. The existing contact implementation and arithmetic are
unchanged. The original felt impact owner reuses this primitive for its
conditioned and virgin power-law segments, with separate compensated affine
normalization and original material history.
