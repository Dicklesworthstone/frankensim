# Material convection through finite electrical, thermal and gas storage

This executable uses `convection::finite_gas::network::step_candidate`. Two
finite source-backed gas chambers exchange heat through three convective
patches on two material cells. The cells also conduct heat to each other and
receive Joule heat from the original finite capacitor/resistor network. Gas
heat and pressure work drive the existing independent inertial walls.

The convection adapter does not implement another circuit, thermal, gas or
momentum solver. It supplies material-resolved conductance to the existing
multi-contact residual. Every nonlinear trial shares one material-network
sample, and every gas contact sees that chamber's actual gas-energy midpoint.
The returned candidate contains ALL finite cells, charge/flux and chamber
endpoints. Internal heat and Joule heat are already included; do not add them
again when consuming the result.

From this directory:

```sh
cargo run --release -- 0.75 1 24 > mixture-a.csv
cargo run --release -- 0.25 1 24 > mixture-b.csv
cargo test
```

The arguments are first chamber's A mole fraction, prescribed relative speed
[m/s], and steps. The second chamber retains 50/50 A/B composition. Geometry,
material laws and the 0.00005 s clock are otherwise identical. The original
source-pack and circuit/thermal manufactured fixtures are reused explicitly;
A/B and the cells are NOT measured real materials. The CSV includes both cell
temperatures, charge, gas temperature/pressure, wall velocity, sound speed and
an independent total-energy discrepancy. No audio signal is synthesized here.

The standalone example workspace has only local path dependencies. It keeps
cross-owner fixture construction out of the production airflow dependency
graph and does not edit the root lockfile. It is not an automatically discovered
root-workspace example; use its manifest path from other directories:

```sh
cargo test --manifest-path crates/fs-airflow/examples/material_convective_network/Cargo.toml
```

## Public adapter

`ConvectiveContact` binds an original `Contact` (ID, chamber, cell) to a
`ForcedGasSurface`. `step_candidate` takes the existing `MaterialNetwork`,
thermal history, finite circuit, explicit per-branch conductivity/resistivity
choices, original chambers, patch list and original work/accuracy controls.

All patches across all chambers must fit each cell's exposed area after
existing contact/reservoir coverage. Duplicate IDs, unknown owners, invalid
geometry, double-counted authored ambient convection and unsupported flow
correlations refuse. Thermal conductivity is required only where the original
material network needs it; a lumped convective patch does not invent a solid
conductivity requirement. Sources are pinned by the original material-network
preparation, while their values remain live at each proposed temperature.

Opposite signed contact transfers remain separate even when their sum is
zero. A shared cell receives both debits, and each chamber is advanced once
per trial with its total heat. Boundary work, circuit supply, cell enthalpy,
gas internal energy and wall kinetic energy remain under the original joint
energy check. No state is accepted after a partial solve or late cancellation.

## Verification and limits

Five authored tests cover explicit original-solver equivalence, exact source
points, patch-order independence, repeated composition feedback into charge
and motion, opposite transfers on one cell, aggregate coverage, invalid
bindings, failure/cancellation and closed-face reduction.

Native Rust compilation and tests remain UNVERIFIED in the authoring
environment: cargo/rustc were unavailable. Independent Python equations ran
128 manufactured updates and their reordered counterparts plus 48 history
steps. Maximum step energy discrepancy was 2.68e-13 J; reordered results were
identical. Changing first-chamber A fraction 0.75 -> 0.25 changed the reference
24-step cell temperature 300.29717 -> 300.32200 K and wall velocity
-0.00073863 -> -0.00112049 m/s, together with charge and sound speed. These are
independent equation checks, not execution of this Rust binary or experiments.

This retains the original first-order thermal-endpoint/gas-midpoint scheme,
uniform cells/gases, fixed cell geometry and independent rigid walls. Relative
speed is prescribed, not a solved fan or boundary layer. Correlations use bulk
gas properties within their admitted domains. No gas mass transfer, chemistry,
radiation, jet thrust, viscous heating or pump-power accounting is added by this
adapter. Those processes are not implied by a successful energy check. Source
events require explicit interval boundaries, and temporal/spatial validity
requires separate assessment.
