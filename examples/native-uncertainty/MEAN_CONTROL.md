# Native cooling mean controls

Version 3 of `fsim-uncertainty-study` connects the existing mean-control
estimators to actual native import/solve/QoI runs. It does not require a user
to supply a gradient, and it never substitutes a linear model for random
physical evaluations.

```lisp
:version 3
:mean-control (coordinate-secant :max-solves 4)
```

The method performs two predeclared calibration solves per nonconstant input:
its low endpoint, then its high endpoint, with other inputs held at their
analytic marginal means. The resulting full-support secant coefficients are
frozen before the first random observation. Constant inputs need no probe.
`max-solves` must cover the complete plan and cannot exceed 64; insufficient
allowances refuse before the ledger opens. Every probe is separately validated
as an ordinary native project. These are **secants, not adjoints**, a numerical
local-derivative certificate, or a promise of variance reduction.

The probe model includes everything the ordinary child solve includes: its
geometry, actual conduction field, material cards, fan operating point,
conjugate exchange and declared radiation. A physically inadmissible probe
refuses the requested study rather than clipping the point, skipping the input,
or silently switching to an independent or uncontrolled model. For a singular
copula, coordinate probes may lie off the joint probability support. They are
deterministic calibration queries within the declared marginal bounds, not
probability observations. Their failure is still terminal.

## Execute and resume

```bash
cargo run --release -p fs-cli --bin frankensim -- --json study \
  examples/native-uncertainty/mean-control.fsim mean-control.db --budget 2

# Use the exact study-... run identifier returned above.
cargo run --release -p fs-cli --bin frankensim -- --json study \
  --resume study-RECEIPT_HASH mean-control.db --budget 16
```

For **version 3 only**, `--budget N` allows N further completed native
calibration/sample evaluations. Zero performs no physics. The first command
stops after two of the example's four calibration probes and returns the usual
budget exit status with a resumable receipt. `:samples` remains the lifetime
random-sample cap, separate from `:mean-control :max-solves`. All work shares
`:wall-time`; elapsed work remains charged across resume. Interrupted native
calls do not consume a completed ordinal; their partial native stages remain
resumable. Versions 1 and 2 retain sample-only invocation-budget semantics.

Completed calibration child runs, exact input vectors, QoIs and coefficients
are retained separately from the raw sample observations. Resume reconstructs
inputs from the retained model, verifies each probe against the sealed native
child, and checks the coefficient bits before continuing. Original source
paths are not needed. A completed native child is reused, not re-solved, and
a failed calibration is not resumable as a successful sample prefix. The
ledger's hashes and seals bind retained contents; they do not authenticate an
external producer or validate the physical model.

## Read the result

`mean_control.calibration.probes` are calibration runs, **not samples**.
`observations`, `statistics`, raw QMC results, empirical compliance and all child
engineering verdicts keep their original meaning. The new
`mean_control.estimate` reports the separately adjusted mean
`average(Y - g dot (X - E[X]))`, raw and controlled descriptive standard errors,
and the observed variance ratio. Ratios above one are retained. A missing or
failed estimate is never represented as zero uncertainty.

The coefficients are in physical parameter units, including Gaussian-copula
inputs. Fixed-count Monte Carlo and replicated randomized QMC are supported.
QMC uses complete equal-sized scrambles as the statistical units, not dependent
individual net points. Incomplete nets remain retained work and are excluded
from both mean estimates. Version 3 cannot combine these fixed-count mean
estimates with the version-2 optional-stopping compliance policy.

The usual `report` and `package` commands export the retained JSON, HTML and
Estimated mean claim. The extra claim is only issued for a completed study with
a representable standard error. Calibration does not bound solver, quadrature,
finite-grid, normal-transform or physical-model error. Automatic nominal-adjoint
extraction remains a distinct, cheaper future coefficient producer; this lane
makes the native mean-control workflow usable without mislabeling its method.
