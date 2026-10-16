# Design against declared physical tolerance scenarios

```bash
cargo run -p fs-ascent --features equilibrium-design --bin equilibrium_fit -- \
  crates/fs-couple/examples/equilibrium-design.model \
  crates/fs-couple/examples/equilibrium-response-limits.fit \
  --scenarios crates/fs-couple/examples/equilibrium-response-tolerances.scenarios \
  --iterations 128 --evaluations 256 --tolerance 1e-8
```

The original model/design files still describe the physical assembly, load cases,
objectives, variable units and response constraints. The additional input declares
**fixed additive offsets** to their design parameters. Nothing samples a distribution
or assigns probabilities, and the nominal command without `--scenarios` is unchanged.

The example varies the actual applied load by -0.1 N, 0 N and +0.1 N. Every realization
must satisfy all four original response limits, including the 0.8 N normal-force cap.
The worst displacement error occurs at low load, but the high-load realization sets
the force limit. A nominal-only design does not account for that distinction.

An independent physical-coordinate reference predicts nominal load 0.80777708764 N,
versus 0.90777708764 N without tolerances. The high-load contact force is 0.8 N; the
worst normalized displacement objective is about 0.9632655963. These are reference
calculations, **not** results from executing the unverified Rust binary.

## Scenario file

```text
frankensim-equilibrium-scenarios-v1
scenarios 3
scenario low-load
offset load-N -0.1
scenario nominal
offset load-N 0
scenario high-load
offset load-N 0.1
```

Every scenario must name **every design variable exactly once**, including explicit
zero offsets. Offset rows may be reordered because they bind by variable name, not
by column position. Names refer to `variable` declarations in the design file, not
to individual fields. Shared variables retain their original shared-field mapping.
Offsets use each variable's physical units: newtons for load, metres for gap, and
the declared stiffness units for the corresponding coefficient. They are not
percentages. The implementation converts offset/scale to a fixed decision shift.
Normal floating-point evaluation applies; this is not interval uncertainty arithmetic.

Input is limited to 64 KiB, 1 through 32 scenarios, and the original maximum of
128 design variables. Missing/extra rows, duplicate names or offsets, unknown
variables, nonfinite offsets, and unrepresentable shifts refuse. The common nominal
domain is the intersection of the original parameter bounds and all shifted bounds;
empty or zero-width intersections refuse. Nothing clips a realized parameter.

## Worst-case objective and constraints

The native `fs_ascent::equilibrium::scenarios` module supplies `ScenarioProblem`
and `ScenarioEquilibriumStudy`. They use the existing physical oracle, analytic
adjoints, response-row adapter and SQP engine. The formulation is:

```text
minimize t
subject to J_s(x + delta_s) - t <= 0         for every scenario s
           c_s(x + delta_s) <= 0 or = 0     for every original response constraint
           nominal x lies in the common parameter domain
```

`t` is an auxiliary objective bound, not a physical parameter. Keeping one smooth
objective inequality per scenario handles ties without manufacturing a derivative
of `max(J_s)`. No scenario can disappear because its objective is not currently
worst: its physical constraints are still enforced. Equalities are required in
every realization and can make the finite problem infeasible. A stalled local
search is not an infeasibility certificate.

A zero-offset scenario is the original physical evaluation. Nominal physical
constraints are enforced only when a zero-offset realization is explicitly supplied;
there is no unrequested extra solve. The original mechanical force/energy limits,
contact-switch derivative exclusions, and adjoint residual checks remain hard
refusals, distinct from optimizable response violations.

## Work, continuation and output

`--evaluations` remains a **total physical-evaluation** limit. One full scenario
family consumes up to S physical evaluations, each visiting every original load
case. At least 2*S evaluations must be available for initialization and the reserved
final re-solve. SQP gets at most floor(evaluations/S)-1 callback attempts; this is
conservative when an out-of-domain trial performs no physical solves. Failed work
is retained and no partial family is reported. `--max-kkt-dimension` includes all
3*n + 1 + S*(1+C) decision/constraint dimensions. Dense SQP is a small-problem path.

In-process studies retain accepted SQP and physical states together. Cancellation
before a segment preserves them; an interrupted search restarts without refunding
spent work. Final `recheck` evaluates every realization again and requires exact
reproduction of the accepted physical evidence. No disk checkpoint is added.

Scenario-mode JSON has its own explicit scope and includes the actual worst
objective separately from t and its positive violation. It retains common-bound,
epigraph and per-scenario response multipliers, realized physical parameters,
source offsets, case predictions, response violations, complete KKT residuals and
both work counts. Model/design hashes keep their existing meanings; scenario
semantics are retained as names and physical offset values, not claimed to share
those hashes. The original nominal JSON is unchanged.

This is local optimization over a **finite, authored** set. No probability of
failure, continuous tolerance-box guarantee, unseen-scenario protection, confidence
interval, global optimality, calibration or dynamic robustness is established.
Contact-active branches are re-admitted at each realization; derivatives across
switches and frictional/static design remain outside the source oracle's scope.

```bash
cargo test -p fs-ascent --features equilibrium-design --lib equilibrium::scenarios
cargo test -p fs-ascent --features equilibrium-design --bin equilibrium_fit --test equilibrium_scenarios
cargo test -p fs-ascent --features equilibrium-design --test equilibrium_fit
```

Seven library, two parsing and four actual-command tests were added. They have not
been executed in the authoring environment, which lacks Cargo/rustc. Independent
NumPy/SciPy references and source checks do not establish native build success.
No dependency, feature or lockfile edges changed.
