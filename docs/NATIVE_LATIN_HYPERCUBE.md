# Native Latin-hypercube cooling studies

`frankensim study` accepts `:method latin-hypercube` in a version-1 native
uncertainty study or a version-3 study with a frozen mean control.
It calls `fs_uq::LhsExecution`; it does not wrap Monte Carlo
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

## Calibrated mean estimates

To estimate the same mean with a fixed linear control, change the study to
version 3 and declare one of these policies alongside its existing LHS layout:

```lisp
:version 3
:mean-control (coordinate-secant :max-solves 4)
```

This example allows two support-endpoint solves for each of two varying inputs.
Alternatively, `:mean-control (nominal-adjoint :max-solves 1)` computes the
supported native derivatives together at the declared physical marginal means.
These calibration solves occur before the first probability observation and
consume the same wall allowance and per-invocation `--budget` as ordinary solves.
The explicit `:max-solves` is their additional lifetime cap; `:samples` still
counts probability observations only. An unsupported derivative or failed probe
refuses calibration before sampling.

The separately reported `mean_control.estimate` averages
`temperature - gradient · (parameters - marginal_means)`. Its descriptive error
uses complete independent hypercubes, named by
`statistical_unit: "complete-independent-latin-hypercubes"`. The original
observations, `lhs` mean, and numerical compliance estimate remain unchanged.
Calibration is not selected or refitted using random observations; it can reduce
or increase variance. The report retains the measured variance ratio, including
values above one. Resume retains completed probes and frozen coefficients, and
continues partial hypercubes from their original ordinals without the source files.

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
Gaussian copulas and version-2 sequential Bernoulli policies are refused for this
method, not silently reinterpreted. Controlled means have the same limits: they
are separate estimators, not adjusted physical temperatures or compliance events.

Focused regressions: `cargo test -p fs-uq product_lhs`,
`cargo test -p fs-project uncertainty::lhs` and
`cargo test -p fs-cli --test native_lhs`. These cover non-power-of-two layouts,
exclusive admission, exact sampler coordinates, retained physical QoIs,
both calibration methods, partial-hypercube statistics and source-independent
checkpoint resumption, including interruption during calibration.
