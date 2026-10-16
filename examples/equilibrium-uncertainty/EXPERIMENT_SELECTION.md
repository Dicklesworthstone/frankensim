# Choose complementary mechanical experiments

`equilibrium_sensitivity` can now select a limited number of complete load
cases from an authored candidate family. It uses the actual local physical
observation adjoints, rather than objective gradients (which can vanish at a
perfect fit). Every target in a selected case is retained as one indivisible
experiment. Neither model nor design file is rewritten.

```sh
cargo run -p fs-cli --bin equilibrium_sensitivity -- \
  examples/equilibrium-uncertainty/sensitivity-independent.model \
  examples/equilibrium-uncertainty/experiment-candidates.fit \
  --rank-relative-tolerance 0.00001 --max-observations 16 --max-adjoints 16 \
  --point-x left 0 --point-x right 0 \
  --select-cases 2 --design-ridge 1 --max-design-factorizations 16
```

The example uses two independently supported masses and three possible load
cases. At the declared point the weighted displacement-Jacobian rows are
`[-10,0]`, `[-9,0]`, and `[0,-4]`. All displacement targets fit exactly. Ranking
only by row norm would choose the first two loads and leave the second support
unobserved. The regularized information criterion instead chooses `strong-left`
and `independent-right`. Their score is `ln(101*17) = 7.448333860897476`, compared
with `ln(182) = 5.204006687076795` for the two redundant left-side loads. These
are independent analytical reference values, not a native-test pass claim.

## The decision being made

For the selected set S, the score is

```text
log det(I + sum_{g in S} J_g^T J_g / ridge_precision)
```

`J_g` contains all of case g's weighted residual derivatives in the declared
scaled parameter coordinates. The residual is
`sqrt(weight)*(displacement-target_m)/scale_m`, exactly as in the original
sensitivity command. Zero weights contribute no information but do not hide
raw physical observation derivatives. The strictly positive ridge makes the
score defined before enough experiments have been selected to span every
parameter. It is an explicit numerical design preference, NOT an inferred
noise variance or a certified Bayesian prior. Changing coordinate scaling,
weights or ridge changes the question being optimized.

The selector greedily tests every unused complete case against the current
selection, then adds the best resolved gain. Numerically unresolved ties retain
declaration order. It is not exhaustive subset optimization, a global-optimum
certificate, or a guarantee of identifiable parameters. Correlated measurement
noise across cases is not modeled by this additive score. No physical experiment
is run or authorized, and no new noise, safety or instrument assumptions are made.

The bounded selector lives in `fs-uq::experimental_design`; `fs-la` remains the
Cholesky owner and `fs-couple` remains the mechanical/adjoint owner. Selection
reuses the complete candidate Jacobian: all candidate cases must first be solved,
even those not chosen. The selected-case budget is NOT the computational budget
for candidate characterization. No additional physics calls are used for scoring.

## Output and budgets

All three selection controls are required together. `--select-cases` is a maximum
number, not a promise to fill the set with zero-information cases. The existing
command output is unchanged when the controls are absent. With them, a separate
`experiment_selection` object records case names/indices, exact observation row
indices, marginal and cumulative scores, ridge, requested/actual cardinality,
actual factorizations and termination. Top-level `information` still describes
ALL candidates, not the selected set; it does not include ridge information.

`selection-limit` means the requested cardinality was reached;
`factorization-budget` means another whole candidate-comparison round could not
fit; `no-resolved-gain` means no gain exceeded the numerical screen. A partial
round never installs a candidate merely because earlier entries were visited.
The score screen is `128*EPSILON*(1+abs(old)+abs(new)+parameters)`, not an
outward-rounded error bound or a physical tolerance. Extremely ill-conditioned
or unrepresentable regularized matrices refuse through the numerical owner.

There are at most 32 parameters, 64 groups, 1024 total rows and 4096 score
factorizations. Admission occurs before physical work when input dimensions and
controls suffice. Original force, energy, contact-margin, adjoint and complete
family checks remain intact; a failed candidate is not silently discarded to
produce a seemingly better experiment set. Constraints remain explicitly
unassessed by this command. Cancellation returns no partial selection. The small
dense factorization polls before and after, not inside its existing kernel.

Focused native tests:

```sh
cargo test -p fs-uq --test experimental_design
cargo test -p fs-cli --bin equilibrium_sensitivity \
  --test equilibrium_sensitivity --test equilibrium_experiment_selection
```

The local authoring environment lacks Cargo/Rust; these commands were not run
successfully there. Independent NumPy/SVD comparisons validate reference
mathematics, not execution of the Rust implementation. The existing uncertainty
workflow includes these focused tests without changing bootstrap or compiler
admission.
