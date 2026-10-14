# One material state drives gas-spring motion and acoustics

This executable connects `ClosedGasParcel::mechanical_load` to the existing
`fs-mbd::RigidBodyIntegrator`. Its Newton solve iterates the midpoint gas force
and rigid-body motion together. The gas force tangent uses the same `rho*c^2`
as sound propagation. Final motion changes volume; the common gas state then
re-evaluates pressure, temperature, sound speed, impedance and transport.
No second gas equation, species registry or rigid-body integrator is introduced.

From this directory:

```sh
cargo run --release -- mixture-a.tsv 500 0.0002 > mixture-a.csv
cargo run --release -- mixture-b.tsv 500 0.0002 > mixture-b.csv
cargo test
```

Or pass this directory's `Cargo.toml` with `--manifest-path` and the mixture
file's path from the working directory. This small standalone example workspace
keeps the cross-owner dependency on `fs-mbd` out of the constitutive kernel and
does not change the main workspace lockfile or production dependency graph.

The two input files differ ONLY in mole fractions. No rebuild or physics-code
change is needed to switch them. Seven whitespace-separated fields per species:
identity, mole fraction, molar mass [kg/mol], constant `cp/R`, linear `cp/R`
slope [1/K], viscosity at 200 K [Pa s], conductivity at 200 K [W/(m K)].
The loader lowers these **synthetic** parameters into the existing NASA-9 and
material-card models. Transport curves double from 200 K to 2000 K. This is an
example input format, not a replacement production material schema. A/B are
not measured real gases. No empirical accuracy or named-gas fidelity is claimed.

The piston is 0.5 kg, its area 0.01 m², initial chamber volume 0.001 m³,
initial gas 400 K / 100000 Pa, exterior 100000 Pa and initial +x velocity
0.2 m/s. The center-applied force preserves this initialized one-axis motion;
there is no collision or joint-constraint solver here. Geometry, loading and
initial conditions are explicit example inputs, not reconstructed historical data.

CSV records motion, all gas observations, restoring stiffness, force residual
and total-energy drift. Energy includes gas internal energy, rigid-body kinetic
energy and the constant-pressure exterior's `p_ext V` contribution. The implicit
midpoint discretization is not exactly energy conserving for nonlinear gas;
its drift and motion must converge with timestep refinement. Four Rust tests
cover refinement, composition swaps, failed-trial atomicity and malformed input.
Native execution was unavailable in the authoring environment; independent
Python equations showed refinement ratios 4.09 and 3.66 for the two fixtures,
with final energy defects below 5e-7 J over 0.1 s at dt=0.001 s. Those results
do not constitute a Rust test pass.

Frozen composition, reversible adiabatic uniform gas only. No combustion,
relaxation, mass flow, condensation, spatial acoustic wave or audio render is
claimed. Source-domain failures and the bounded 32-iteration coupling budget
return errors rather than clamping the gas or publishing a partial body update.
