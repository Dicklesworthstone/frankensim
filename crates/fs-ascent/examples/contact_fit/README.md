# Fit shared physical parameters against independent load experiments

```bash
cargo run -p fs-ascent --example contact_fit
cargo run -p fs-ascent --example contact_fit -- \
  --data observations.txt --scale-m 0.0002 --iterations 128 --evaluations 256
```

This is a working-path reference consumer of the **existing** fallible L-BFGS or SQP,
nonlinear preload and implicit adjoint implementations. It does not introduce
another optimizer. The native executable and tests have not been run in the
authoring environment; the numerical references described below are independent.

The example fixes a simple physical rig: a supported 0.04 kg translating mass
presses against one 100 rad/s mass-normalized elastic coordinate through a
quadratic normal contact. The unknowns are support stiffness in N/m, contact
coefficient in N/m² and contact gap in metres. Mass, basis, observation geometry,
exponent and damping remain fixed. The mass's actuator/observation shape is
`1/sqrt(0.04) = 5`; no unit shape is substituted for it.

Each optional input row is one **independent** experiment:

```text
0.8 0.0003190067261188236 0.00006085959643287059
1.6 0.0004212400757625386 0.00013472559545424769
2.5 0.0005284197010700035 0.0002182948179357998
```

The columns are `load_N mass_displacement_m receiver_displacement_m`. There is
no header, ignored column or inferred unit conversion. Input is bounded at
8192 bytes and 32 complete finite rows. These particular rows are disclosed
**synthetic** targets from an independent quadratic-equilibrium formula, with
planted parameters 600 N/m, 1.8e8 N/m² and 0.2 mm. They are not specimen data.
Without `--data`, the example uses those targets; supplied data are labelled
`provided-displacements`, never automatically promoted to validated measurements.

## What the loop actually does

The reusable domain entry point is
`fs_couple::render::schedule::force::coupled::equilibrium::sensitivity::objective::design::EquilibriumDesign`.
It accepts other caller-authored modal/free-mass networks, springs, normal
contacts, physical load cases, displacement targets and named decision fields.
It has no dependency on the optimizer layer.

For every trial, the evaluator decodes all physical parameters and rebuilds
fresh zero-state candidate networks. Each case is independently settled under
its own loads. The existing force-balance/activity admission is rechecked, then
one existing adjoint per case produces all selected parameter derivatives.
Objective and derivatives are summed across cases in fixed order. Shared
fields receive the same parameter and contribute additively to its gradient.
A case-specific actuator variable contributes only to that experiment, retaining
its physical per-newton derivative before decision scaling.

Supported fields are spring stiffness/rest extension, contact stiffness/gap/
quadrature weight, and physical actuator force. A decision can share several
fields only when they have the same parameter kind; shared contact stiffnesses
must also have the same exponent. Duplicate assignments and incompatible units
refuse. Modified contact-law values are labelled inverse-design candidates, not
the original calibrated coefficients. Bases and shape maps are not re-extracted.

The default `--solver lbfgs` sends this callback to `LbfgsState::try_new` and
`try_run`. The only recoverable trial barrier is a decoded parameter outside its
declared domain, using the engine's existing positive-infinity rejection rule.
Physical solve failures, activity-margin refusals, cancellation and exhausted
work limits propagate as errors; they are not disguised as unfavorable designs.
The zero placeholder gradient for an infinite barrier is never accepted.

## Bounds, accounting and continuation

Decisions use `physical = reference + scale * decision`, with initial decisions
zero. Defaults are:

| Parameter | Reference | Scale | Admitted physical interval |
|---|---:|---:|---:|
| Support stiffness, N/m | 400 | 400 | 200 to 1000 |
| Quadratic contact coefficient, N/m² | 1e8 | 1e8 | 2e7 to 4e8 |
| Gap, m | 0.0001 | 0.0002 | 0.00001 to 0.0005 |

For L-BFGS these bounds define an admissible domain, **not** a box-constrained
optimizer. A boundary solution can stall rather than satisfying a constrained
KKT test. The SQP route below instead treats every bound as an inequality and
reports its multiplier. The default synthetic optimum lies strictly inside the
domain. Neither route claims global optimality or uniqueness/identifiability.

`--scale-m` is an explicit objective normalization, not an estimated noise level.
L-BFGS uses a dimensionless decision-gradient infinity norm of 1e-8. SQP uses
all four KKT residuals at 1e-8, not the unconstrained gradient norm.
`--evaluations` includes initialization, rejected domain trials, failed attempts
that start work, and a reserved final physical re-solve. Its admitted range is
2 through 4096. `--iterations` permits 0 through 512 additional iterations.
Per-evaluation case visits have a separate cumulative ceiling; the original
primal/adjoint setup, iteration, force, energy and activity limits remain intact.

The in-process `FitSession` example retains the original L-BFGS checkpoint and
physical work counters. Cloning at an accepted-iteration boundary preserves
continuation. Cancellation before a new run leaves it unchanged. An interrupted
line search retains spent work and resumes by restarting that search; no
bit-identical in-flight replay or durable serialized checkpoint is promised.

The final JSON names the actual stop reason, decoded physical values, fresh
objective/predictions, maximum observation discrepancy and adjoint residual,
alongside iterations, total evaluations and attempted cases. A budget stop is
not reported as convergence. Final reporting does not trust only a cached
optimizer objective; it re-solves every load case at the accepted parameters.

## Box-constrained physical fitting with SQP

```bash
cargo run -p fs-ascent --features equilibrium-design --example contact_fit -- \
  --solver sqp --iterations 256 --evaluations 1024
cargo run -p fs-ascent --features equilibrium-design --example contact_fit -- \
  --solver sqp --iterations 256 --evaluations 1024 \
  --data crates/fs-ascent/examples/contact_fit/bounded-observations.txt
```

`--solver sqp` uses the reusable `fs_ascent::EquilibriumStudy`, not another
example-local optimizer. The optional `equilibrium-design` feature enables an
L4-to-L3 dependency on the existing physics owner. Without the feature the
option refuses before reading observations or performing physical work; it
never silently falls back to a different engine. L-BFGS remains the default.

`EquilibriumStudy` supplies the exact summed physical adjoint gradient and
signed box residuals to the existing small-dense `SqpState`. Each physical
residual is divided by its declared decision scale, so the bound Jacobian
contains exact +/-1 rows. The dense admission cap includes decisions plus both
faces of each box: nine for this three-variable rig. Out-of-domain trials are
explicitly unavailable, not clamped or assigned an invented objective/gradient.
Original primal, adjoint and activity-margin refusals still stop the solve.

The reusable study borrows an immutable problem and an exclusive caller-owned
`DesignControl`. Read-only `optimizer()` and `accepted()` refer to the same
accepted point and complete case family, even after a failed search. Accepted
step boundaries can be resumed without reevaluation; spent work is not reset.
A physical budget can only be extended, not refunded. This borrowed study is
not a cloneable independent work allowance or a serialized checkpoint.

The SQP JSON additionally reports stationarity, primal feasibility, dual
feasibility and complementarity, the raw gradient norm, bound multipliers and
all case predictions in metres. Multipliers are ordered lower/upper support,
lower/upper contact coefficient, then lower/upper gap, in **dimensionless
decision coordinates**. A nonzero raw gradient at an active bound is compatible
with constrained stationarity. `converged` follows the actual stop reason;
`kkt_within_tolerance` separately exposes the numerical residual test. The
reserved final physical re-solve must reproduce the retained complete case
evidence; a mismatch refuses output rather than publishing a stale design.

`bounded-observations.txt` contains disclosed **synthetic**, not measured, data
from support=600 N/m, contact coefficient=8e8 N/m² and gap=0.2 mm. Its contact
coefficient exceeds this rig's allowed 4e8 N/m². An independent SciPy SLSQP
reference gives approximately 616.223 N/m, 4e8 N/m² and 0.183839 mm, with a
nonzero objective of 0.000552698 and an active upper-contact multiplier. These
numbers do not execute the repository's SQP or certify experimental validity.

## Verification and remaining scope

Six `fs-couple` tests cover analytic multi-case gradients, shared-field sums,
physical load units, late case failure, work/cancellation replay and invalid
bindings. Four tests in this example cover parameter recovery with the actual
L-BFGS interface, split-run equivalence, cancellation, budgeted final audit, and
strict external data/options. An additional selection test preserves the
original default and rejects invalid/repeated solver flags. Four opt-in example
tests cover nonlinear SQP recovery, a genuinely active physical bound, budgeted
final audits, and typed physical refusal. Six reusable-study tests cover both
active box faces, independent cases, replay, cancellation and cumulative work.

```bash
cargo test -p fs-couple --test equilibrium_design
cargo test -p fs-ascent --example contact_fit
cargo test -p fs-ascent --features equilibrium-design --test equilibrium_study
cargo test -p fs-ascent --features equilibrium-design --example contact_fit
```

`fs-couple` is an optional normal dependency for the reusable study and remains
an existing dev-dependency for the example; `fs-dcontact` remains dev-only. The
default optimizer runtime dependency graph is unchanged. Cargo/rustc are unavailable in the authoring environment, so
native compilation, Rust tests, formatting, Clippy and lockfile regeneration
remain unverified. Run ordinary Cargo before a `--locked` check to refresh the
manifest's dev-dependency edges; no external numerical runtime was added.

The independent NumPy/SciPy check differentiates the physical-coordinate
quadratic equations and re-solves perturbed cases; it does not execute Rust or
`fs-math`. Its synthetic recovery is not experimental material identification.
This static objective cannot identify damping, transient impact response,
frequency-dependent radiation or uncertainty. Near contact switches, the
existing exclusion-distance rule still refuses a smooth derivative.
