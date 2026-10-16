# Uncertainty through an authored mechanical equilibrium

`equilibrium_uq` reads the same model and design files as `equilibrium_fit`.
It samples declared uncertainty and executes the existing `fs-couple`
forward equilibrium owner for every independent load case at every draw. It does
not optimize, integrate artificial warm-up dynamics, or introduce another solver.
The output is one selected physical displacement in metres, with a one-sided
compliance event `displacement_m <= limit_m` (not absolute displacement).

Run the supplied one-mode example:

```sh
cargo run -p fs-cli --bin equilibrium_uq -- \
  examples/equilibrium-uncertainty/linear.model \
  examples/equilibrium-uncertainty/linear.fit \
  --method mc --samples 1024 --seed 73 \
  --case applied-force --target 0 --limit-m 0.0078125 --independent \
  --uniform-x force-N -1 1
```

Every variable requires exactly one `--uniform-x NAME LOWER UPPER` or
`--fixed-x NAME VALUE`. These are **dimensionless design coordinates**, not newtons
or metres: physical values remain `reference + scale*x`, as declared by each
binding in the design file. Here `F=0.5+0.5*x` N, so the requested uncertainty is
uniform force in `[0,1]` N. Sharing across physical fields is preserved. Argument
order does not change the sampler's variable order. Input files are never edited.

The example's modal stiffness and port map give `u=F/64` m. Its population mean
is `1/128` m, standard deviation `1/(64*sqrt(12))` m, and the declared compliance
probability is `1/2`. These are synthetic analytical controls, not experimental
validation of a physical specimen.

All probability-law declarations, their endpoint decoding, case/target selection,
and resource bounds are admitted before sampling. Independence is never inferred
from a box. This command currently admits only independent uniform or fixed
variables, not arbitrary copulas, Gaussian tails, or epistemic intervals. A box
is a uniform probability law only because the caller explicitly declares it so.

The original sample cap and `samples * load_cases` case-solve cap remain hard
ceilings. Each draw uses `EquilibriumDesign::evaluate_forward`: original force,
penetration, full-network residual and stored-energy checks apply, but no tangent
or adjoint is built. Contact activity margins restrict derivatives, not these
forward observations; the inverse-design derivative API retains its exclusions.
Any failed model evaluation terminates the run with no partial JSON statistics;
no draw is replaced, clipped to a domain, or discarded. No new dynamics, transient
response, contact-switch derivative, memory-allocation or real-time guarantee is
claimed. Without a confidence policy the CLI executes the complete finite plan;
with the optional MC policy below it may stop on a completed prefix. It has no
disk resume interface.

Output includes both source hashes, seed, laws, physical coordinate scaling,
actual work counts, mean, sample standard deviation, descriptive standard error,
and empirical compliance. Schema v2 labels the physics path `primal-only` and
limits compliance to the selected displacement. Any declared response constraints
not assessed by this query are counted explicitly; they are not deemed feasible. Evidence stays **Estimated**. Neither standard error
nor a finite observed range establishes a confidence bound or validates physical,
material, geometric or numerical model error. Optional-stopping validity applies
only to the MC compliance bounds explicitly requested below, not these descriptive
statistics or the QMC standard errors.

Native checks to execute in the project's Rust/constellation environment:

```sh
rch exec -- cargo test -p fs-cli --bin equilibrium_uq --test equilibrium_uq
```

These Rust checks have not been run in the authoring environment; it lacks `rch`,
`cargo`, and `rustc`.

## Replicated randomized QMC

To use the existing Owen-scrambled Sobol owner instead of Monte Carlo, replace
`--method mc` above with `--method rqmc --replicates 4`. The total `--samples`
ceiling must equal `replicates * points_per_net`, with 2..=256 replicates and a
power-of-two net size of at least two. QMC admits at most ten declared variables;
there is no silent Monte Carlo tail or dropped remainder. `mc` rejects a
`--replicates` option rather than ignoring it.

This path calls the upstream `QmcExecution` directly, not a new sampler. Each
complete net has its own scramble; observations are never passed to the MC
confidence-sequence code. The reported displacement-mean and compliance standard
errors are computed across complete-net means and proportions. They do not imply
confidence intervals, exactness at zero observed variation, or bounds on finite
digital-grid bias. `displacement_std_dev_m` is null for QMC: a between-net standard
error is not the distribution's physical standard deviation. No stopping rule
based on these descriptive standard errors is implemented.

The coupled-contact example from `crates/fs-couple/examples/equilibrium-design.*`
also works with this command. Declare `--uniform-x support-N-per-m 0 0.5`,
`--uniform-x contact-N-per-m2 0 0.8`, and `--uniform-x gap-m 0 0.5`, select
`--case load-1.6N --target 1 --limit-m 0.00014`, and keep `--independent` explicit.
That varies physical spring stiffness in [400,600] N/m, quadratic contact
coefficient in [1e8,1.8e8] N/m², and gap in [0.0001,0.0002] m. It consumes the
computed receiver displacement, not the optimization objective or a synthetic
sound amplitude. All three load cases are solved independently for each draw.

## Forward contact transitions

`contact-onset.model` and `contact-onset.fit` describe a physical spring plus a
unilateral linear contact. The design deliberately gives zero tangent/query
budget and an exclusion margin of 0.02 m: it cannot support the adjoint path,
but the original primal force and energy checks can still admit observations.

```sh
cargo run -p fs-cli --bin equilibrium_uq -- \
  examples/equilibrium-uncertainty/contact-onset.model \
  examples/equilibrium-uncertainty/contact-onset.fit \
  --method rqmc --replicates 4 --samples 1024 --seed 73 \
  --case applied-force --target 0 --limit-m 0.015625 --independent \
  --uniform-x force-N -0.5 0.5
```

Here `F=1+x` N is uniform in [0.5,1.5]. The contact begins at F=1 N and
q=1/64 m. Below onset q=F/64; above it q=(F+2)/192. The exact population mean
is 11/768 m and the declared compliance probability is 1/2. Both pieces describe
the same continuous displacement, but their one-sided derivatives differ.
`--fixed-x force-N 0` evaluates onset itself without fabricating a derivative.
This example remains a synthetic reduced-model control, not specimen validation.

## Decide a fixed compliance requirement with Monte Carlo

Declare both `--require-probability P --confidence-alpha A` to ask whether
`P(displacement_m <= limit_m) >= P`. Both values must lie strictly inside (0,1);
alpha must have a finite reciprocal. The required probability, event, model,
law and alpha are fixed before samples are observed.

```sh
cargo run -p fs-cli --bin equilibrium_uq -- \
  examples/equilibrium-uncertainty/contact-onset.model \
  examples/equilibrium-uncertainty/contact-onset.fit \
  --method mc --samples 4096 --seed 73 \
  --case applied-force --target 0 --limit-m 0.02 --independent \
  --uniform-x force-N -0.5 0.5 \
  --require-probability 0.9 --confidence-alpha 0.05
```

This uses **the existing `UqExecution::assess_compliance` owner**, not a new
confidence-sequence formula. It assesses the retained indicators at sample
counts 2,4,8,... and the final sample cap. A lower bound at or above the required
probability yields `compliance-satisfied`; an upper bound strictly below it yields
`compliance-violated`. If neither separation occurs, the budget stops with an
`inconclusive` decision, even when every observed draw passed. Neither outcome
uses the sample standard error to declare certainty.

The JSON separates `stop_reason` from the sampler's `status`. An early decision
has `status: "budget-truncated"` because the original planned cap was not filled;
a full-cap report has `status: "complete"`, which alone does not mean compliance.
`planned_samples` and actual `samples`/`case_solves` remain distinct. The
`compliance_assessment` object reports decision, requirement, alpha, interval,
method and sample count; it is null when the policy was not requested.
A computed violation is a successful command result, not a process error:
consumers must inspect the decision. Physical failures or interruption produce
an error and no partial JSON or recycled confidence result.

Time-uniform coverage is inherited only under the original fixed-conditional-mean,
bounded-indicator sampling assumptions for the declared numerical model. The
floating-point interval is not an outward-rounded certificate and does not bound
physical-model, material, discretization or finite-digital-grid error. Estimates
of mean and dispersion at the data-dependent stop remain descriptive. Searching
across models, seeds, thresholds or policies needs separate multiplicity control.
Evidence remains **Estimated**, never physical validation.

`rqmc` refuses these confidence flags before reading inputs or solving physics:
within-net points must not be treated as independent MC observations. Its existing
complete-net standard errors are unchanged and remain descriptive.
