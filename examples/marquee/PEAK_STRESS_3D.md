# Sampled-peak stress constraints in native 3-D design studies

`sampled-peak-limited-simp` minimizes projected material volume under a bound
on the **retained numerical stress samples**. It uses the existing CutFEM
operator, exact discrete adjoint, physical density map and projected
augmented-Lagrangian optimizer. The existing `stress-limited-simp` type keeps
its original normalized volume/load-weighted average; it is not silently changed.

```sh
cargo run --release -p fs-cli --features sdf3-study --bin frankensim -- \
  --json study examples/marquee/bracket-3d-peak-stress.fsim peak-study.db
```

The example declares numerical demonstration material and loading values, not
an engineering material specification or a physically validated design. It uses
two independent equilibria on the same fixed domain and an explicit p = 16.
A successful process invocation or satisfied discrete constraint does not establish
mesh convergence, manufacturability, global optimality or continuum safety.

## What is constrained

For each retained bulk quadrature point in each independent case, define

```text
r        = projected physical density in the point's cell
k(r)     = e_min + (1 - e_min) * r^penal
s_ref    = von_mises(C_reference : strain(u_case))
s_relax  = r^q * s_ref
s_phys   = k(r) * s_ref
A        = (sum_over_all_cases_and_points(s_relax^p + s_phys^p))^(1/p)
constraint = A / stress_limit_pa - 1 <= numerical_tolerance
```

There is **no division by volume, quadrature weight, number of samples, or
load weight** in A. A small-volume hotspot and a case with a small or zero
weight cannot disappear into an average. Both stress terms are intentional:
qp relaxation retains the declared low-density strength penalty, while the
physical term prevents a prescribed zero-density region from erasing its
modeled ersatz stress. The old density-floor/turnover admission still applies
to optimizable material; the new functional does not loosen that protection.

If M is the largest relaxed or physical stress sample and N is the number of
retained bulk points summed over cases, then in exact arithmetic

```text
M <= A <= (2*N)^(1/p) * M.
```

The evaluator also checks A against its independently computed numerical
relaxed and physical maxima before returning a result. It refuses a failed
bound instead of clipping the value with an inconsistent derivative. This is
not outward-rounded certification: the stresses, equilibria and aggregate are
numerical. In particular, no maximum between samples or in the continuum is
bounded. The accepted feasibility tolerance permits A up to the declared limit
multiplied by `(1 + tolerance)`.

Conservatism depends on N and p. More samples or duplicated cases can increase
A even if the physical maximum is unchanged. Increasing p tightens the bound
but concentrates sensitivities. Keep the same retained sampling within each
optimization; independently assess discretization sensitivity before drawing
engineering conclusions. The existing admitted range is 2 <= p <= 64.

## Input and load semantics

Select the new type in the existing stress optimizer block; all other required
fields, order, SI units and work controls remain the same:

```lisp
(optimizer
  :type sampled-peak-limited-simp
  :initial-density 0.75
  :filter-radius-m 0.15
  :max-updates 80
  :density-floor 0.05
  :stress-limit-pa 16.0
  :relaxation-power 1.0
  :aggregation-power 16.0
  :max-stress-points 50000
  :max-evaluations 2000
  :max-backtracks 40
  :tolerance 0.000002
  :penal 3.0
  :beta 0.0)
```

Every declared body/pressure/traction case is solved independently. Weights
remain explicit family metadata and are reported in normalized form, but do
not affect this constraint. Peak mode admits weights in [0,1], requiring at
least one positive weight for the family declaration. Zero-weight cases still
consume their original point and solve budgets and contribute to the adjoint.
The original normalized-average and compliance modes retain positive-only
native weight admission. Nonzero imposed displacement is still outside this
fixed-reference-load stress interface.

Physical solid/void regions, constructive domains, embedded supports, and the
explicit multilevel solver use the existing shared builders. A physical
zero-density region retains `e_min` stiffness, loads and quadrature. Use
constructive subtraction for an actual empty bore or cavity; a density label
is not geometric removal.

## Optimization, output and recovery

The full derivative includes the direct relaxed **and** physical density terms,
the equilibrium adjoint, and the same filter/projection pullback. The native
initial gradient gate checks the selected functional before optimization. One
current-density preparation is shared across all primal and adjoint cases;
there is no extra equilibrium per physical/relaxed stress component.

Only accepted states replace the current design. Accepted optimizer steps can
be infeasible; the least-volume accepted feasible incumbent is retained and
exported separately. Cancellation, point exhaustion and rejected work preserve
that existing ownership. Point limits count each retained point once, not each
of its two stress terms; the implementation combines the pair without doubling
the point-family allocation.

The report identifies `unweighted-sampled-peak-bound`, marks the finite sampled
maximum constraint explicitly, and keeps continuum/unqualified maximum-stress
authority false. It includes per-case maxima, point count and feasibility
policy. Reports and evidence packages use retained fields without a new solve.

The canonical optimizer type changes the study identity. Stress resume restores
the accepted optimizer state and re-solves its accepted and distinct incumbent
endpoints using the same original functional, loads, geometry, regions and
solver under the original work allowances. It does not replay earlier optimizer
steps or authorize switching from a peak problem to an averaged problem.

Focused tests:

```sh
cargo test --release -p fs-topopt --features cutfem-marquee --lib sdf3::stress
cargo test --release -p fs-cli --features sdf3-study --lib study::elasticity::sdf3::stress::peak_tests
```

These exercise real independent cut equilibria, full-chain gradients, hidden
hot cases, physical ersatz stresses, accepted material updates, recovery and
interruption. Test source is not an execution receipt.
