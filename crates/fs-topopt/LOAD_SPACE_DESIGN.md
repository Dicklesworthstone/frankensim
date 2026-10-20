# Worst simultaneous-load design

```sh
cargo run --release -p fs-topopt --features cutfem-marquee \
  --example worst_load_sdf3 -- --load-space 200 250000
```

The example minimizes the worst response over a **continuous set of simultaneous
loads**, not just the listed basis columns taken separately. Its two columns are
an authored Y body load and an authored correlated Y/Z body load. The command
allows 200 accepted optimizer updates and one cumulative 250,000 Krylov-iteration
budget. A charged initial solve fixes a single compliance normalization scale
for the entire set. Existing finite-case, cap-constrained and binary modes do
not change.

## Load set and objective

`envelope::load_space::LoadSpace3` defines a fixed scale and spectral temperature
for the supplied density-independent basis F:

```text
f(q) = F q,       ||q||_2 <= 1
G = F^T K_eroded^-1 F / reference
maximum = largest_eigenvalue(G)
J = temperature * log(sum_k exp(eigenvalue_k(G) / temperature))
```

The columns already include physical amplitudes and correlations. They are not
probabilities or mutually exclusive scenarios. This is a zero-centered ellipsoid
(possibly singular), not a box of independently bounded coefficients, and not an
ellipsoid around a nonzero nominal load. A change of orthonormal basis preserves
the set; adding duplicate columns generally does not. The implementation admits
at most 16 columns to bound the existing dense eigensolver's indivisible work.

Every basis load is solved on the same eroded material, including the existing
support, ghost and hanging-node treatment. ALL cross compliances are retained.
The small symmetric eigensolver is reused, with checks on reciprocity before
averaging numerical cross terms, eigenpair residual, orthogonality and positive
semidefiniteness to the declared numerical gate. Eigenvalues are not clipped.
The worst force and its displacement are linear combinations of the independently
solved basis fields. They are exported with the maximizing coefficients.

The smooth objective differentiates the entire spectral function rather than
choosing one eigenvector's gradient. This handles repeated eigenvalues: the
worst direction may be nonunique while the smooth objective remains differentiable.
In exact arithmetic `maximum <= J <= maximum + temperature*ln(dimension)`.
That algebraic statement is not a floating-point enclosure or continuum bound.
The exact discrete SIMP/projection/filter derivative is computed from weighted
modal stiffness contractions, with no finite-difference optimization.

## Constrained optimization

`load_space::design::LoadSpaceDesignStudy3` reuses `EnvelopeDesignOptions3` for
the raw box, dilated-volume cap, gradient-check policy and existing projected
augmented-Lagrangian solver. The baseline spectral objective and dilated volume
must pass numerical gradient checks. A `run` segment retains multipliers,
spectral step, accepted fields and cumulative budgets. The best feasible
incumbent minimizes the smooth objective; the actual worst compliance is also
reported and its decrease is not implied by every smooth-objective decrease.

`recover_feasible` uses the existing strict dilated-volume restoration and a
fresh load-space evaluation. It returns a separate candidate, without changing
AL state, history or optimizer counters. The example selects it only when it
beats the feasible accepted incumbent, labels its source and emits solved
fields before propagating a late failure. It never promotes an iteration limit
to an optimizer-convergence claim.

Only fixed-grid linear elasticity, homogeneous supports and fixed dead-load
columns are covered. There is no probability/risk model, nonzero load offset,
load-dependent stiffness, stress/buckling safety, geometric length-scale claim,
continuum certificate or global-optimality guarantee.

```sh
cargo test --release -p fs-topopt --features cutfem-marquee \
  --lib sdf3::robust::envelope::load_space
cargo test --release -p fs-topopt --features cutfem-marquee \
  --example worst_load_sdf3 load_space::tests
```
