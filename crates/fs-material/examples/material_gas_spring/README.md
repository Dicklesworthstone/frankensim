# One material state drives heated gas-spring motion and acoustics

This executable connects the existing source-backed `ClosedGasParcel` and
`fs-mbd::RigidBodyIntegrator`. Its Newton solve iterates gas force and body motion
together. `midpoint_heat_volume` now resolves **simultaneous heat and volume work**
in conservative internal-energy/volume coordinates. The exact same midpoint
pressure supplies the piston force and gas work. The material owner supplies the
discrete pressure derivative; the example does not implement another gas law,
caloric inversion or rigid-body integrator.

From this directory:

```sh
cargo run --release -- mixture-a.tsv 500 0.0002 > mixture-a.csv
cargo run --release -- mixture-b.tsv 500 0.0002 > mixture-b.csv
# Same model and geometry, plus an explicit 100 W source during motion:
cargo run --release -- mixture-a.tsv 500 0.0002 100 > heated-a.csv
cargo run --release -- mixture-b.tsv 500 0.0002 100 > heated-b.csv
cargo test
```

The optional fourth argument is signed constant **watts**, integrated over each
actual step before calling the material API in joules. Omitted means zero heat.
A failed thermal or mechanical trial, invalid source domain, exhausted bounded
solve or late cancellation publishes neither accepted owner. The existing
no-source entry point and first twelve CSV columns remain; explicit heat input
adds accumulated supplied heat, step gas work, step energy residual and step gas
entropy change. `total_energy_drift_j` subtracts accumulated supplied heat.

Or pass this directory's `Cargo.toml` with `--manifest-path` and the mixture
file's path from the working directory. This standalone example workspace
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
0.2 m/s. The centered load preserves this initialized one-axis motion;
there is no collision or joint-constraint solver here. Geometry, loading and
initial conditions are explicit example inputs, not reconstructed historical data.

## Energy and material coupling

The same immutable material model supplies temperature-dependent heat capacity,
pressure, sound speed, impedance and transport at the midpoint AND endpoint.
With nonlinear heat capacity, the midpoint temperature is NOT the mean of the
endpoint temperatures. The Newton derivative is the discrete material pressure
response in final volume, not the static isentropic stiffness. The CSV stiffness
observation remains the endpoint's static adiabatic `rho*c^2*area^2/volume`.

Each accepted step checks, in absolute joules:

`delta(gas internal energy) + delta(body kinetic energy) + p_ext*delta(volume) = supplied heat`.

The original rigid-body step supplies actual kinetic-energy change; the material
step supplies actual stored energy and pressure work. No post-step correction
rescales momentum or gas temperature to force agreement. `advance_heated` accepts
whole-step heat, so tests turn the source off without resetting gas or motion.

Six driver tests plus the existing input test cover motion refinement, data-only
swaps, work balance, source-off history and atomic failure/cancellation. Six more
library tests cover the constitutive update and its analytic discrete derivatives.
Native tests were **not run** in the authoring environment: `cargo` and `rustc`
are unavailable. Independent Python equations ran 13 trajectories (2,896 steps),
including a 6 J pulse and source shutoff. Maximum total-energy discrepancy was
below 5e-13 J; heated displacement refinement ratios were 4.008 and 4.058 for the
two mixtures. These are independent equation checks, not Rust test results.

## Limits

Frozen composition, spatially uniform gas and explicitly prescribed heat only.
Midpoint is not an exact finite isentrope: reported gas entropy change is NOT
entropy production or a certificate of nonnegative entropy production. Motion
and thermal trajectories need timestep refinement. No combustion, mass flow,
condensation, molecular relaxation, spatial acoustic wave or audio render is
claimed. The material temperature brackets, original source support and
32-iteration mechanical coupling budget remain hard refusals, not clamps.
