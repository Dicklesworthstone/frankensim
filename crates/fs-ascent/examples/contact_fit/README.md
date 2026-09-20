# Fit shared physical parameters against independent load experiments

```bash
cargo run -p fs-ascent --example contact_fit
cargo run -p fs-ascent --example contact_fit -- \
  --data observations.txt --scale-m 0.0002 --iterations 128 --evaluations 256
```

This is a working-path reference consumer of the **existing** fallible L-BFGS,
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

The example sends this value/gradient callback to `LbfgsState::try_new` and
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

These bounds define an admissible domain, **not** a projected/box-constrained
optimizer. A boundary solution can stall rather than satisfying a constrained
KKT test. The synthetic optimum lies strictly inside the domain. No global
optimality or uniqueness/identifiability claim is made for arbitrary data.

`--scale-m` is an explicit objective normalization, not an estimated noise level.
The stop rule is a dimensionless decision-gradient infinity norm of 1e-8.
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

## Verification and remaining scope

Six `fs-couple` tests cover analytic multi-case gradients, shared-field sums,
physical load units, late case failure, work/cancellation replay and invalid
bindings. Four tests in this example cover parameter recovery with the actual
L-BFGS interface, split-run equivalence, cancellation, budgeted final audit, and
strict external data/options.

```bash
cargo test -p fs-couple --test equilibrium_design
cargo test -p fs-ascent --example contact_fit
```

The two additional `fs-ascent` dependencies are **dev-only existing workspace
crates**, needed by this executable/test consumer; normal optimizer dependencies
are unchanged. Cargo/rustc are unavailable in the authoring environment, so
native compilation, Rust tests, formatting, Clippy and lockfile regeneration
remain unverified. Run ordinary Cargo before a `--locked` check to refresh the
manifest's dev-dependency edges; no external numerical runtime was added.

The independent NumPy/SciPy check differentiates the physical-coordinate
quadratic equations and re-solves perturbed cases; it does not execute Rust or
`fs-math`. Its synthetic recovery is not experimental material identification.
This static objective cannot identify damping, transient impact response,
frequency-dependent radiation or uncertainty. Near contact switches, the
existing exclusion-distance rule still refuses a smooth derivative.
