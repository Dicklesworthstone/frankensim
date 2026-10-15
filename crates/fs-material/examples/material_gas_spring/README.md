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
# Heat transfer from a 500 K reservoir through a material gas layer:
# transfer area 0.02 m^2 and conduction distance 0.001 m
cargo run --release -- mixture-a.tsv 200 0.0005 wall 500 0.02 0.001 > reservoir-a.csv
cargo test
```

The optional fourth argument is signed constant **watts**, integrated over each
actual step before calling the material API in joules. Omitted means zero heat.
Alternatively, `wall temperature_K area_m2 distance_m` selects material-driven
exchange with a constant-temperature reservoir. The source uses the actual
midpoint material, not the initial or previous-step conductivity.
A failed thermal or mechanical trial, invalid source domain, exhausted bounded
solve or late cancellation publishes neither accepted owner. The existing
no-source entry point and first twelve CSV columns remain; explicit heat input
adds accumulated supplied heat, step gas work, step energy residual and step gas
entropy change. `total_energy_drift_j` subtracts accumulated supplied heat.
Reservoir mode additionally appends `step_source_residual_j`; the original
12-column no-source and 16-column prescribed-power CSV layouts are unchanged.

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

## Material-driven conductive reservoir

The explicit lumped gas-layer source is

`Q = dt * k(midpoint material) * transfer_area / distance * (reservoir_T - midpoint_T)`.

The transfer area is independent of piston area. The existing species material
cards and mixture owner supply `k`; no air preset, new property table or
conductivity interpolation is added. The source, gas heat/work and rigid-body
force are iterated together from the same accepted state. Force updates include
the material owner's discrete `d(midpoint_pressure)/dQ` and
`d(midpoint_pressure)/d(final_volume)` responses. Both the force residual and
whole-step heat-source residual must converge before either owner is published.
The source law is checked separately from first-law energy accounting; a good
energy balance alone cannot admit the wrong heat flux.

Only changing the input conductivity columns therefore changes absorbed heat,
temperature, pressure, piston motion and the shared sound-speed observation.
A colder reservoir extracts heat; zero transfer area explicitly isolates it.
Source shutdown does not reset mechanical or thermodynamic history. The finite
32-iteration solve can refuse stiff thermal exchange; no convergence for all
conductances or timesteps is claimed, and no fallback freezes the source law.

Four additional Rust tests cover the heat-source parser and local property
response, conductivity-only changes across the full motion trajectory, cooling,
isolation parity and cancellation/refusal without partial publication. They are
added but **not executed** here because the Rust toolchain is unavailable.
Independent standard-library Python equations exercised six reservoir
trajectories (1,200 steps). Maximum cumulative first-law discrepancy was below
`8e-10 J`, and the largest accepted heat-source mismatch was below `3e-11 J`.
These are synthetic numerical checks, not native test passes or measured gas
validation.

## Limits

Frozen composition and spatially uniform gas only. The reservoir option is a
lumped conductive gas layer connected to an infinite heat-capacity reservoir;
its distance and area are caller assumptions. It is not a resolved wall PDE,
convection correlation, finite wall thermal mass or wall-material model. Thermal
conductivity is sampled at the gas midpoint, not integrated across a spatial
layer with variable properties. The reservoir's absolute temperature may differ
from the source-supported gas temperatures; every actual gas material evaluation
still must be inside its source domain.
Midpoint is not an exact finite isentrope: reported gas entropy change is NOT
entropy production or a certificate of nonnegative entropy production. Motion
and thermal trajectories need timestep refinement. No combustion, mass flow,
condensation, molecular relaxation, spatial acoustic wave or audio render is
claimed. The material temperature brackets, original source support and
32-iteration mechanical coupling budget remain hard refusals, not clamps.
