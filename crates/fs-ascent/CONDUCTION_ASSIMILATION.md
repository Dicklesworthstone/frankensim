# Spatial temperature reconstruction

`conduction-assimilation` connects the existing tetrahedral P1 conduction
producer to `transient::variational::{WeakConstraintWindow, joint::JointWindow}`
and their existing L-BFGS studies. Enable it explicitly; it also enables
`transient-design`. The default optimizer dependency graph is unchanged.

## Executable CSV path

```sh
# A labeled synthetic example using the production finite-element solver.
cargo run --release -p fs-ascent --features conduction-assimilation \
  --example conduction_assimilate

# Inspect the exact vertex indices and coordinates of the example's mesh.
cargo run --release -p fs-ascent --features conduction-assimilation \
  --example conduction_assimilate -- --nodes

# Reconstruct with four physical steps per observation interval.
cargo run --release -p fs-ascent --features conduction-assimilation \
  --example conduction_assimilate -- readings.csv --model-sigma 0.03 --substeps 4
```

The header is `time_s,node,temperature_k,sigma_k`. Node IDs are zero-based mesh
vertices, not free-DOF indices. Rows must be ordered by acquisition time; a
repeated node/time pair is rejected. Independent measurement scales must be
positive and finite. Input is capped at 64 KiB, 256 readings, eight positive
observation times, and ten seconds. Missing readings are omitted, not replaced
by zero. There is no interpolation or nearest-time substitution.

The example uses a **declared**, not experimentally validated, 0.12 x 0.04 x
0.04 m slab: 45 vertices, 96 tetrahedra, conductivity 2 W/(m K), volumetric heat
capacity 1000 J/(m3 K), 300 K prescribed at x=0, and insulated other faces.
Its source is `amplitude * x/0.12` W/m3. It jointly fits the source amplitude,
one sensor bias, and 36 free temperatures at every time point. Source amplitude
is an unbounded signed density coordinate; negative values represent a sink.
It is not a total-watt parameter or identification of an arbitrary source field.

The initial background is 300 K with 0.1 K independent scales. The amplitude
prior is 1000 +/- 2000 W/m3; the bias prior is 0 +/- 0.2 K. These are fixed
Gaussian standard deviations. `--model-sigma` is an independent **endpoint
increment** scale in kelvin, not white-noise intensity or a solver tolerance.
The example prints loss components, fitted parameters, and all reconstructed
nodal fields. It exits unsuccessfully when the gradient tolerance is not met.
No posterior standard deviations or physical-validation claims are emitted.

## Solver resolution is separate from inference resolution

`--substeps N` accepts integers from 1 through 32 (default 1). It divides each
observation interval into N production backward-Euler steps, without adding
observations, reconstruction controls, or model-error priors. The synthetic
case still has four knots and 146 controls at every refinement. The endpoint
model-error scale is not divided by N, multiplied by dt, or applied at each
fine step. Refining the numerical forecast can change the fitted parameters;
it does not change the declared statistical tradeoff.

The default one-step mode retains its original numerical path without extra
replay. Refined mode retains small endpoint witnesses and recomputes the fine
trajectory through binary checkpoints, keeping one inner PDE linearization
live instead of all N sparse Jacobians. Counts in the output distinguish
native forward steps, reverse replay steps, and optimization controls. This
trades additional solves for bounded retained trajectory memory.

Without a CSV, synthetic readings are generated using the selected fine clock.
That is a manufactured self-consistency demonstration, not independent physical
validation. Comparing resolutions for supplied readings is useful, but the
option itself provides neither an adaptive error estimate nor a temporal-error
certificate. Unrepresentable intermediate times are rejected rather than
silently skipped.

## Library composition

Implement `ConductionWindowModel` with a borrowed `BackwardEuler` capacity
object and an endpoint `ConductionProblem` for every interval. Construct
`ConductionWindowPolicy` with the same mesh, fixed prescribed boundary values,
explicit time grid, `Cx`, and numerical/mesh caps. Build window priors and
controls in `policy.free_vertices()` order; use `gather_field` and `expand_field`
for checked conversion to and from physical full fields. The observation loss
must apply the same map. A reading at a prescribed vertex can still inform a
sensor-bias parameter, but cannot change that vertex's prescribed temperature.

Use `WeakConstraintStudy::new` with the policy for a fixed physical model, or
`JointWindowStudy::new` with a `ParameterFamily` for shared parameters. Override
`parameter_pullback` to include parameter effects through physical inputs.
`StepLinearization` already provides source-density and capacity contractions;
the default callback explicitly declares parameter-independent physical inputs.
Direct sensor derivatives belong to `WindowObjective::parameter_partials`.

A base `ConductionWindowPolicy` interval is one backward-Euler step. For several
steps per reconstruction interval, create a `SubstepGrid` from the coarse
measurement times and construct the base conduction policy on `grid.fine_times()`.
Wrap that policy with `CheckpointedIntervals::new(base, grid, budget)` from
`transient::variational::intervals::substeps`. Pass the wrapper to the same
`WeakConstraintStudy` or `JointWindowStudy`; build the window and its priors on
`grid.knot_times()`, not on every solver substep. `SubstepGrid::new` also admits
an explicit nonuniform fine grid with declared knot indices.

With the wrapper, `problem(k)` and `parameter_pullback(k, ...)` receive the
GLOBAL FINE-STEP index. Source/boundary schedules must supply those same indices
on replay. Do not reinterpret k as the observation interval or evaluate a
fresh random load during a repeated step. `SubstepBudget` bounds state/parameter
dimensions, parked outer checkpoints and inner record calls per reverse sweep;
each conduction solve retains its own work limits. The example admits at most
256 fine steps, six parked outer states, and 512 reverse fine-record calls per
coarse interval. Its optimizer evaluation count is not a count of PDE solves.

The adapter retains the actual transient Jacobian, supports the producer's
optional nonlinear k(T) and matching contacts, and never differentiates the
Krylov/Newton stopping decisions. In-solve cancellation uses the supplied `Cx`
gate; the extra window callback is checked between physical operations.
Changed replay endpoints or failed solves/derivatives publish no result.
Endpoint witnesses are diagnostics, not certificates of supplied derivatives.

Every model forecast passes the existing energy/residual gates. A reconstructed
history includes separately penalized model-error increments and need not itself
conserve energy. Geometry, prescribed-temperature values and covariance scales
remain fixed. Air-network feedback, radiation sensitivities, covariance learning,
certified derivatives and parameter identifiability are outside this adapter.

```sh
cargo test --release -p fs-ascent --features conduction-assimilation \
  --lib 'conduction_assimilation::'
cargo test --release -p fs-ascent --features conduction-assimilation \
  --example conduction_assimilate
```
