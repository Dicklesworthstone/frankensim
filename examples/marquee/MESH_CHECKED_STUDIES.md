# Mesh-checked native elasticity optimization

The existing `projected-volume` mode can require measured finer-grid agreement
before accepting a design update. This addresses coarse-grid-only improvements
without changing the optimizer, geometry representation, loads or material model.

```sh
cargo run -p fs-cli --bin frankensim -- --json study \
  examples/marquee/bracket-mesh-checked-2d.fsim /tmp/mesh-checked.db
```

The supplied numerical example may stop with `mesh-unresolved`; it is not a
fixture whose tolerances are chosen to guarantee success. Inspect the actual
retained grid differences. A finer optimization lattice is a separately funded
new study, not a silently substituted geometry or a reset resume budget.

Add this optional block after `:cg-poll-iters`, before any `:design-regions`:

```lisp
    :mesh-check (refinement
      :extra-levels 1
      :absolute-compliance-tolerance-j 0.000001
      :relative-compliance-tolerance 0.05
      :area-tolerance-m2 0.005)
```

All four fields are required. One or two extra dyadic levels are admitted,
through level seven. Each adjacent compliance difference must be no larger
than the absolute allowance plus the relative allowance times finer-grid
compliance. Adjacent area changes and finer area versus the original target
must meet the declared area allowance. These are separate from the stricter
area-projection gate on the optimization grid.

The actual accepted baseline is solved on every declared grid. Each otherwise
admissible coarse candidate is independently solved on the same grids and
must reproduce the requested strict compliance decrease on all of them. Its
area must also remain comparable with the baseline on each grid. Failure
contracts the existing candidate search; an unresolved baseline stops before
any proposal. Prolongation retains the bilinear field up to interpolation
rounding. There is no fine-grid area projection or redistancing to manufacture
a favorable comparison, and no interpolated displacement is called a solve.

`constraints.mesh_resolution` in the normal ledger receipt and report records
the policy and last complete search. It includes actual baseline rungs and,
when accepted, candidate rungs, with level, compliance, area and field snapshot.
Earlier checks remain in the predecessor receipt chain. Interrupted work does
not replace the last complete check or accepted geometry. The ordinary
`study --resume`, `report` and `package` commands remain the interfaces.

`mesh-unresolved` returns the budget exit code and a durable unchanged endpoint,
not `completed`. A failed bounded candidate search remains `no-feasible-descent`.
A terminal receipt returns unchanged on resume. Extra numerical work consumes
the same lifetime wall budget and shares CG/phase cancellation. Assembly,
prolongation and individual quadrature operations remain indivisible.

This is Estimated, observed grid sensitivity. Passing may miss a common bias:
it is not a guaranteed continuum error bound, observed-order theorem, continuous
stress certificate, physical validation or optimum. The optimization lattice
is not automatically adapted. Omitting the block preserves the existing
canonical source and unguarded projected-volume numerical path.

Library entry points are `fs_topols::resolution::assess_compliance_resolution_controlled`
and `ProjectedOptimizer::advance_one_resolution_controlled`. Focused tests run
through the existing topology and native-study lanes; execution results must
be checked separately from the presence of authored regression tests.
