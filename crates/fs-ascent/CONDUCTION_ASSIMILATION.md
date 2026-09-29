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

# Reconstruct from observations of that declared slab.
cargo run --release -p fs-ascent --features conduction-assimilation \
  --example conduction_assimilate -- readings.csv --model-sigma 0.03
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

Each interval is **one** production backward-Euler step. Refine the supplied
time grid for temporal accuracy. The adapter retains the actual transient
Jacobian, supports the producer's optional nonlinear k(T) and matching contacts,
and never differentiates the Krylov/Newton stopping decisions. In-solve
cancellation uses the supplied `Cx` gate; the extra window callback is checked
between physical operations. A failed solve or derivative publishes no result.

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
