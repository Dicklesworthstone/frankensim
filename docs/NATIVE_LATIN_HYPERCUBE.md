# Native Latin-hypercube cooling studies

`frankensim study` accepts `:method latin-hypercube` in a version-1 native
uncertainty study. It calls `fs_uq::LhsExecution`; it does not wrap Monte Carlo
samples in an LHS label or introduce another physical solver.

```bash
frankensim --json study data/reference-project/cooling-uncertainty-lhs.fsim study.db
```

The tracked example declares three independent hypercubes with ten strata per
coordinate, totaling thirty native import/solve/QoI pipelines. Counts need not
be powers of two. The product must exactly equal `:samples`, whose native limit
remains 256. Marginals are the native study's explicit uniform laws; dependence
must be `independent`. The existing parameter bindings, material domains,
geometry quarantine, physical solve limits and advisory decision gate apply.

## Interrupted runs

```bash
frankensim --json study data/reference-project/cooling-uncertainty-lhs.fsim study.db --budget 7
frankensim --json study --resume "$RUN" study.db
frankensim --json report "$COMPLETED_RUN" study.db
frankensim --json package "$COMPLETED_RUN" study.db
```

`RUN` and `COMPLETED_RUN` are the corresponding `run` strings returned by the
commands. The first command intentionally returns the budget exit code with a
retained seven-sample prefix. Resume uses the sealed base project, source
geometry and card packs in the ledger; it does not need the original files.
The exact canonical study, seed, replicate layout, solver inputs and accepted
sample ordinals remain bound. An unfinished model evaluation retries the same
point; a failed evaluation terminates rather than deleting an inconvenient
observation. Completed studies remain terminal.

## Reports and scope

The result's `method` is `latin-hypercube`. Its `lhs` object reports complete
replicates, partial-replicate work, replicate means, the mean temperature and
empirical numerical compliance fraction. Only complete hypercubes enter those
estimates. A partial hypercube remains paid, retained work. Before one complete
hypercube the estimate is null; before two complete hypercubes the descriptive
standard error is null. The ordinary iid Monte Carlo `statistics` field stays
null. HTML and evidence packages identify LHS explicitly.

The standard error measures between-replicate variation, not independent-point
variation inside one hypercube. It is not a confidence interval, an
optional-stopping guarantee, or a bound on numerical, finite-grid or physical
model error. Child engineering uncertainty budgets and verdicts are unchanged.
Gaussian copulas, version-2 sequential Bernoulli policies, and version-3 mean
controls are refused for this method, not silently reinterpreted.

Focused regressions: `cargo test -p fs-project uncertainty::lhs` and
`cargo test -p fs-cli --test native_lhs`. These cover non-power-of-two layouts,
exclusive admission, exact sampler coordinates, retained physical QoIs,
partial-hypercube statistics and source-independent checkpoint resumption.
