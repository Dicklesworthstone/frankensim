# CONTRACT: fs-tribo

> Status: `frankensim-ext-tribo-dry-baseline-tgbj` baseline. This is not an
> admitted material-data result or an Euler-disc validation result.

## Purpose and layer

`fs-tribo` is an L3 dry-contact constitutive leaf. It provides typed dry
friction, elastic Hertz reference forms, scalar rolling/contour loss laws,
explicit dissipated-work accounting, and caller-owned Archard wear state. It
does not depend on contact-solver internals: the consuming solver owns contact
complementarity, tangent reactions, and residual embedding.

## Public types and semantics

- `InterfaceSystemRef` names ordered surfaces, history, medium, and a nonblank
  caller `source_id`. `InputAuthority` is a caller-declared ceiling
  (`CallerDeclared`, `SyntheticFixture`, or `Estimated`), never an admission or
  verification receipt. Every constitutive response returns this provenance.
- `ContactFrame` normalizes a finite contact normal. `TangentialSlip` accepts
  only a finite velocity with no material normal component; it refuses rather
  than silently projects a closing/separating velocity into friction.
- `FrictionLaw` provides Coulomb, velocity-dependent, and Stribeck rungs. At
  zero slip it reports static capacity only; it never invents a stick reaction.
- `HertzSpherePlane` and `HertzCylinderPlane` are G1 elastic closed forms,
  expressed in SI. They use caller-supplied radius/modulus and retain the
  caller ceiling; a mathematical closed form does not validate its inputs.
- `ConstantRollingMoment` uses `M = -sign(omega) N a` (N m) and
  `P = N a |omega|` (W). `ConstantContourForce` uses
  `F = -sign(v) F_c` (N) and `P = F_c |v|` (W). They are separate mechanisms;
  callers must not combine them without independent source support.
- `HeatPartition`, `DissipationStep`, `WorkLedger`, and `WearState` have
  private invariant fields and checked getters. `ArchardLaw` computes
  `dV = k F ds/H` in m3 and commits only a validated candidate state.

## Invariants

- All public construction and evaluation paths refuse missing identity,
  non-dry medium, negative/non-finite input, and every non-finite derived
  candidate from finite extreme inputs.
- Tangential traction and scalar resistance oppose their declared relative
  rates; output dissipated power/work is finite and non-negative.
- Heat shares are finite, non-negative, and sum to one. A dissipation step's
  channels are finite, non-negative, and close to its total before a ledger can
  mutate. Ledger and wear candidate totals are checked before assignment.
- No mutable material table, hidden state, or unordered material-pair lookup
  exists. Ordered interface and history identities remain explicit.

## Error model

`TriboError` is a total refusal surface for missing identity, non-dry media,
invalid/overflowing physical input, non-finite vectors, normal slip, invalid
partitions, and forged/invalid dissipation states. A refusal makes no partial
state change to `WorkLedger` or `WearState`.

## Determinism class

Pure scalar operations are deterministic for the same input and ISA. Norms use
the platform `hypot` sequence to avoid avoidable intermediate overflow. No
cross-ISA bit-stability claim is made.

## Cancellation behavior

This bounded scalar leaf has no asynchronous work and no cancellation scope.
Callers scheduling many contact points own cancellation between calls.

## Unsafe boundary

None. The crate forbids unsafe code.

## Feature flags

None.

## Conformance tests

Inline G0/G1/G3 tests cover identity/dry/tangent refusals; analytic rigid-block
equilibrium; independent numerical Hertz values and force-pressure-radius
cross-relations; reversal/scaling/provenance; resistance sign; partition
closure; forged negative/non-finite work rejection; and rollback on work/wear
overflow. Test coefficients are explicitly synthetic fixtures.

## No-claim boundaries

This crate does not mint material admission, calibration, InterfaceSystemCard
query receipts, flash-temperature or thermal-port solutions, roughness,
adhesion, plasticity, finite-patch partial slip, lubrication/EHL, contact
geometry evolution, wear geometry updates, stop time, Euler-disc ranking,
one-millimetre optimum, or experimental/video correspondence. A caller may
record an estimated or synthetic result, but must not promote it based on this
crate alone.
