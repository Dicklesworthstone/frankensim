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
