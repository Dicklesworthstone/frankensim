# Optimize physical responses, not only parameter bounds

```bash
cargo run -p fs-ascent --features equilibrium-design --bin equilibrium_fit -- \
  crates/fs-couple/examples/equilibrium-design.model \
  crates/fs-couple/examples/equilibrium-response-limits.fit
```

This example varies a physical actuator load to approach a receiver displacement
while limiting contact force, indentation, support reaction and receiver travel.
It reuses the existing supported 40-gram mass and quadratic contact model.
The initial 1.5 N load violates the response limits but is mechanically valid:
it is a legitimate SQP starting point, not a failed simulation.

An independent physical-coordinate reference predicts a final load of about
0.90777709 N and a normal contact force of 0.8 N. The actuator's 0.1..4 N design
bounds are inactive. The force constraint, rather than a parameter bound, sets
the optimum. These numbers are synthetic reduced-model predictions, not native
Rust execution, measured material properties or a global-optimality claim.

## Additive input section

The existing design-v1 file may end after its variable bindings, exactly as
before, or append one complete section:

```text
constraint_limits 4
constraints 4
constraint normal-cap 0 contact-force 0 at-most 0.8 1
constraint indentation-cap 0 contact-penetration 0 at-most 0.00012 0.0001
constraint support-floor 0 spring-force 0 at-least -0.2 1
constraint travel-cap 0 displacement 1 0 at-most 0.00009 0.0001
```

A row is `constraint NAME CASE QUANTITY LOCATION SENSE BOUND SCALE`.
`displacement` uses `COMPONENT PORT` for its location, referring to the original
model's complete mass-normalized port map. The other quantities take one index
in the original spring or normal-contact order. Case indices are zero-based.
Names must be unique. Counts are explicit and hard-capped at 64. Every field is
required, and unknown/trailing rows refuse rather than disappearing.

| Quantity | Physical value | Bound and scale units |
|---|---|---|
| `displacement` | Signed attachment displacement | m |
| `spring-force` | Signed reaction on the left: `-k*(left-right-rest)` | N |
| `contact-force` | Nonnegative compressive normal reaction | N |
| `contact-penetration` | Positive closure beyond the gap, zero if separated | m |

The senses are `at-most`, `at-least`, and `equal`. Upper and lower constraints
produce `(value-bound)/scale <= 0` and `(bound-value)/scale <= 0`, respectively.
Equalities produce `(value-bound)/scale == 0`. Scale must be finite and strictly
positive; it normalizes the optimizer residual, not the physical constitutive
law or an inferred measurement uncertainty. No residual is squared or clipped
before SQP sees it. To constrain absolute signed force, declare both faces.

## Derivatives and numerical admission

`EquilibriumDesign::with_constraints` attaches immutable requirements to the
same native problem used by `EquilibriumStudy`. Each row reuses its case's
already-solved equilibrium and prepared tangent; at most one extra adjoint is
needed per row. There is no extra primal solve per parameter or constraint.
Reaction derivatives retain BOTH the implicit displacement effect and direct
stiffness/rest/gap/weight dependence. Shared physical fields sum all applicable
terms, and actuator derivatives retain physical-newton scaling.

The existing physics ceilings remain hard refusals. Design bounds should sit
inside those ceilings: SQP may evaluate a design-constraint violation, but it
cannot ask the simulator to bypass force, energy or penetration admission.
All normal contacts must still be separated from their activity switches by the
original derivative margin. There is no smooth-gradient claim across a switch.

Each adjoint and pullback retains the existing sensitivity query budget; the
constraint count bounds their repetition. The dense SQP cap now includes
`3*variables + physical_constraints`, with both parameter-box faces included.
Repeated evaluations retain the existing cumulative evaluation/case ceilings.
A failed or cancelled evaluation returns no partial row family.

The command's `constraints` JSON array reports each physical value, unit, bound,
scale, signed residual, violation, adjoint residual and normalized multiplier.
Equalities use the equality multiplier family; physical inequalities follow
parameter-box multipliers. KKT residuals include all rows. The reserved final
physical re-solve must reproduce the accepted constraints and Jacobians too.
A budget stop can retain positive violations and is not relabelled convergence.
With no constraint section, input hash domains and command output remain unchanged.

## Focused native checks

```bash
cargo test -p fs-couple --test equilibrium_constraints
cargo test -p fs-couple --lib render::schedule::force::file::design
cargo test -p fs-ascent --features equilibrium-design --lib equilibrium::constraints_tests
cargo test -p fs-ascent --features equilibrium-design --bin equilibrium_fit --test equilibrium_fit
```

These commands were attempted in the authoring environment, which has no Cargo
or rustc. Native compilation, tests, rustfmt and Clippy remain unverified.
Independent numerical checks exercise the derivation, not Rust or its
deterministic math implementation. No dependencies or lockfile edges changed.
