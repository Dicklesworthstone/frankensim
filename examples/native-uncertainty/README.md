# Probability studies of a native cooling project

This example varies the heat input and ambient temperature of the existing
native cooling reference project. Each sample imports the real tetrahedral STL,
resolves the AA6061 material card, executes the staged cooling solve, and reads
the retained maximum-temperature QoI. The study declares its probability laws
and their independence explicitly; engineering ranges and material tolerances
do not become probability distributions automatically.

Run from the repository root:

```bash
cargo run -p fs-cli --bin frankensim -- --json study examples/native-uncertainty/study.fsim native-uq.db
```

The completed result contains a `study-…` run identifier and retained report,
HTML, observations, checkpoint, and evidence-package hashes. Export the report
and package using that exact identifier:

```bash
cargo run -p fs-cli --bin frankensim -- --json report study-RECEIPT_HASH native-uq.db
cargo run -p fs-cli --bin frankensim -- --json package study-RECEIPT_HASH native-uq.db
```

The report lists the sampled parameters, child run and project identities,
retained QoI receipt, and maximum temperature for every accepted sample. After
all eight samples complete, it also reports the mean temperature, sample
standard deviation, and empirical fraction at or below the native requirement's
limit after its explicit margin is applied.

## Pause and resume under a sample budget

`--budget` limits the number of additional native samples in this invocation.
For example, start with two:

```bash
cargo run -p fs-cli --bin frankensim -- --json study examples/native-uncertainty/study.fsim partial-uq.db --budget 2
```

This returns `budget-exhausted` and exit code 6, with a retained partial run.
Use its returned identifier to finish:

```bash
cargo run -p fs-cli --bin frankensim -- --json study --resume study-RECEIPT_HASH partial-uq.db
```

`--budget 0` retains the admitted study and assets without executing a sample.
Resume uses the retained assets, so the original project, mesh, and card files
need not remain at their initial paths. Partial runs retain their observations
but publish no completed fixed-count statistics. A refused child solve stops
the study and cannot silently disappear from its sample plan.

## Stop when a probability target is resolved

The version-2 `compliance.fsim` example declares a probability decision before
sampling: a required pass probability of 0.5, significance level 0.05, at least
eight observations, and a lifetime cap of 32 native solves. Run it with:

```bash
cargo run -p fs-cli --bin frankensim -- --json study examples/native-uncertainty/compliance.fsim compliance-uq.db
```

Every completed native QoI contributes exactly the binary event
`temperature-max <= limit - margin`. The fixed Beta(1/2,1/2) Bernoulli likelihood
mixture yields a probability confidence sequence. Once the minimum count is
met, the study stops if its lower bound reaches the required probability or its
upper bound falls below it. The result is `decision-reached`, exit 0, with a
`meets-probability-target` or `below-probability-target` decision. Reaching the
sample cap with an interval that still straddles the target returns
`budget-exhausted`, exit 6, and `indeterminate`.

The report's `compliance` object contains the empirical frequency, both
confidence endpoints, policy, sample count, and decision. Even an all-pass
sample retains a nonzero interval width. Policy runs keep `statistics: null`:
a sample stopped for this decision does not become a fixed-count temperature
distribution. The same sample-budget, resume, report and package commands
apply. Resume keeps the original retained probability policy; editing the
source file does not change the question for an existing run.

## Interpretation

These are advisory estimates under the declared independent uniform input
model. Version 1's empirical pass fraction is a fixed-count descriptive
statistic. Version 2 supplies a time-uniform probability confidence sequence
under the stated Bernoulli sampling assumptions. Its deterministic floating-point
calculation is not an outward-rounded certificate. Neither route supplies
compliance signoff, certifies continuum error, validates the material model,
or replaces the native QoI's engineering uncertainty budget. The unchanged
base project's solver, material domains, requirements, and per-solve budgets
still apply to every sample.

## Fan speed through the native coupled cooling model

`fan-speed.fsim` varies the existing finned heatsink's fan speed, power and
inlet temperature. Its three uniform laws and independence are illustrative
study declarations, not measured manufacturing distributions. Run it with:

```bash
cargo run -p fs-cli --bin frankensim -- --json study examples/native-uncertainty/fan-speed.fsim fan-speed-uq.db
```

The `fan-speed-ratio` target names a `fan-system` bank by its exact `:id`:

```lisp
(uniform :name "fan-speed" :target fan-speed-ratio :entity "heatsink-bank"
  :low 0.6 :high 0.8)
```

Both bounds are positive dimensionless numbers inside that bank's existing
`:speed-domain-lo`/`:speed-domain-hi`. Each sample sets the **absolute speed
ratio relative to the bank's source curve**. For example, a sampled 0.8 sets
`:speed-ratio 0.8` even when the base project has `:speed-ratio 0.7`. The curve,
its source and tolerance, rated point, fan count, arrangement, topology and
other banks retain their original declarations. The native flow solver
applies the existing fan affinity law and resolves the operating point; the
heatsink's `airflow-convection` boundary then recomputes heat transfer and the
coupled solid/air temperature. A boundary with a fixed declared convection
coefficient keeps that coefficient and need not show a thermal response to
fan speed.

Missing banks, duplicate identities, invalid fan systems, nonpositive speeds,
and support extending beyond the declared speed domain refuse before any
sample runs. Fan speed does not expand a convection card's operating regime:
if a sampled physical solve leaves its admitted regime, the study retains the
refusal. The same sample-budget and source-independent resume commands apply.
Contact resistance remains selected from immutable interface cards and is not
a mutable uncertainty target.

## Replicated randomized Sobol quadrature

The version-1 `qmc.fsim` example uses the same native fan-speed, heat-input and
air-inlet parameters with four independently scrambled Sobol nets of four
points each. Every point still runs the ordinary native cooling pipeline.

```bash
cargo run -p fs-cli --bin frankensim -- --json study examples/native-uncertainty/qmc.fsim qmc-uq.db
```

Select the method and its entire layout before sampling:

```lisp
:version 1
:samples 16
:method quasi-monte-carlo
:qmc (owen-scrambled-sobol :replicates 4 :samples-per-replicate 4)
```

The sample count must equal the replicate count times the points per
replicate. There must be at least two replicates, with at least two points per
replicate; only the point count must be a power of two. The native limit
remains 256 full solves, and Sobol quadrature supports at most ten declared
parameters. Zero-width uniform laws remain valid. The parser refuses missing
or inconsistent layouts, excess dimensions, a QMC layout on Monte Carlo, and
QMC combined with the version-2 Bernoulli stopping policy.

The report keeps `statistics: null` and adds a `qmc` object with the exact
layout, completed replicate count, unfinished point count, retained replicate
means, temperature mean and probability of the numerical pass event.
Standard errors come from variation between complete independently scrambled
nets. Dependent points inside a net are never treated as independent
Bernoulli trials or used to calculate the Monte Carlo standard error.

The existing `--budget` and resume commands work at any point, including
inside a net. An unfinished net is retained paid work but is excluded from
estimates. One complete net supplies a mean with an unavailable standard
error; two or more supply the descriptive between-replicate error. Resume
uses the same retained assets, layout, scramble keys and next point, and
produces the same completed observations, report, package and checkpoint
bytes as uninterrupted execution. Refused physical samples terminate the
study and suppress every quadrature estimate.

Randomized quadrature may reduce integration error for a given solve budget,
but its reported standard errors are descriptive, not confidence intervals
or permission to stop when a result looks favorable. They do not bound
finite-grid bias, native numerical errors or physical-model errors; zero
variation between replicates does not prove an exact result.

## Dependent inputs with an explicit joint law

`dependent-inputs.fsim` runs the native heatsink with dependent fan speed,
power and inlet temperature, reusing the same uniform physical supports and
replicated Sobol layout. Its dependence is an illustrative declared model,
not an estimate fitted to manufacturing data:

```bash
cargo run -p fs-cli --bin frankensim -- --json study examples/native-uncertainty/dependent-inputs.fsim dependent-uq.db
```

The existing `:correlation` field accepts either `independent` or a Gaussian
copula with an explicit latent-normal correlation matrix:

```lisp
:correlation (gaussian-copula
  :latent-correlation ((1 0.4 0.2) (0.4 1 0.5) (0.2 0.5 1)))
```

Rows and columns follow the exact parameter declaration order: fan speed,
power, inlet temperature in this example. The sampler draws correlated
standard normals, maps them through the normal CDF, then maps each quantile
into that parameter's declared uniform support. The supplied correlations
describe the **latent normals**, not Pearson correlations of the physical
uniform inputs. A bare correlation matrix does not specify the joint law and
is refused. The report retains the model name, matrix, coordinate meaning and
parameter order beside the sampled physical inputs and child-run identities.

The numerical owner admits finite, symmetric, positive-semidefinite matrices
with a unit diagonal and entries in [-1, 1] before creating the study ledger.
Consistent singular matrices are allowed, including perfect positive or
negative dependence. Constant marginals remain exact constants and still
occupy their declared matrix coordinate. The factorization and normal CDF
are floating-point calculations, not distributional or PSD certificates.
No draw is clipped, rejected, fitted or replaced to force agreement.

The same joint law works with fixed-count Monte Carlo and with version-2
Monte Carlo probability decisions. Dependence is within each input vector;
Monte Carlo vectors remain independently drawn, so their raw native pass
indicators can feed the existing Bernoulli confidence sequence. QMC still
uses complete independent replicates for descriptive errors and refuses
the Bernoulli stopping policy. JSON, HTML and evidence-package statements
identify the declared law in all three modes.

Sample budgets and source-independent resume retain the original matrix,
physical supports, sampler and next sample. Editing the source's dependence
declaration does not change an existing run. Each physical sample must still
satisfy the unchanged native fan, material and convection domains; a refused
physical solve terminates propagation rather than being discarded.
