# CONTRACT: fs-tribo

> Status: baseline implementation for `frankensim-ext-tribo-dry-baseline-tgbj`; it is not an
> admitted material-data or Euler-disc result.

## Purpose and layer

`fs-tribo` is an L3 dry-contact constitutive leaf. It supplies typed friction, Hertz reference
responses, optional rolling/contour loss laws, heat/work accounting, and Archard wear evolution.
It does not depend on contact-solver internals. A solver owns contact complementarity and embeds
the returned traction, state, tangent, and accounting records into its residual.

## Public semantics

- `InterfaceSystemRef` names an **ordered** surface/system/history identity. The caller may bind it
  to an admitted material card or explicitly declare a synthetic fixture. The latter remains
  `Synthetic` and cannot become material authority.
- `FrictionLaw` provides Coulomb, affine velocity-dependent, and Stribeck-plus-viscous rungs.
  At zero slip the law reports a stick capacity rather than fabricating a tangential reaction;
  contact owns that solve.
- `HertzSpherePlane` and `HertzCylinderPlane` are G1 analytic reference forms. They require
  `AnalyticReference` authority plus caller-supplied finite positive reduced modulus and radius.
  They are not plastic, rough-surface, adhesive, layered, or finite-patch claims.
- `ResistanceLaw` is the generic rolling/contour loss interface. The supplied constant-moment and
  constant-contour-force laws dissipate non-negative power; their coefficients are caller-declared
  and never calibrated here.
- `DissipationStep` conserves declared frictional work into surface-A heat, surface-B heat, and an
  explicit other channel. `WorkLedger` accumulates only non-negative declared work.
- `ArchardLaw` evolves a caller-owned `WearState` by `dV = k F ds / H`. Hardness appears only in
  this wear law; no ductility, plastic-contact, damage, fracture, or life claim is inferred.

## Invariants

- Every public constructor and evaluation refuses non-finite, negative, missing, or incompatible
  inputs. A dry law refuses a wet or undeclared-medium context.
- The interface order and history are non-empty, identity-bearing strings. Coefficients require an
  explicit authority source; no material-pair lookup or default coefficient exists.
- Traction and resistance oppose their corresponding velocity. Reported dissipated power and work
  are non-negative. Heat partitions close to one within floating-point roundoff.
- Evaluation is pure apart from an explicit caller-owned `WearState` or `WorkLedger`; neither
  material data nor hidden global state is mutated.

## Error model and no-claim boundary

`TriboError` is a total refusal surface. `UnsupportedAuthority`, `SyntheticInput`, and
`NotDryInterface` are ordinary outcomes, not fallback modes. The crate carries no material-card
adapter until the seed interface data is admitted. It has no calibration, temperature evolution,
roughness, mixed/EHL lubrication, flash-temperature solution, wear-volume geometry update,
contact-patch solution, plasticity, or Euler-disc prediction claim.

## Determinism, cancellation, unsafe

The functions are deterministic pure scalar arithmetic on one ISA; no cross-ISA bit claim is made.
There is no asynchronous work or cancellation scope in this bounded leaf. There are no unsafe
blocks.

## Evidence

Inline tests cover G0 admission/refusal and accounting laws, G1 Hertz closed forms and rigid-block
stick/slip threshold, plus G3 scaling, reversal, partition, and replay metamorphics. All numerical
coefficients used by tests are declared `SyntheticFixture` inputs only.
