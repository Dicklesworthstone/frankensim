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

## Interpretation

These are advisory estimates under the declared independent uniform input
model. The empirical pass fraction is a fixed-count descriptive statistic,
not a confidence sequence or compliance signoff. It does not certify continuum
error, validate the material model, or replace the native QoI's engineering
uncertainty budget. The unchanged base project's solver, material domains,
requirements, and per-solve budgets still apply to every sample.
