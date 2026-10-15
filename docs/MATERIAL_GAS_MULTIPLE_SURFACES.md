# Several material walls, one gas state

`GasPiston::material_network_surfaces_candidate` extends the existing finite
circuit / thermal network / gas coupling to multiple distinct wall surfaces.
Every surface retains its temperature, material source and signed heat. All
surfaces share ONE gas energy midpoint and ONE pressure-work/momentum solve.
The existing single-surface APIs remain available. Both scalar and multi-contact
boundaries use the concurrently landed shared `boundary::heat` implementation;
this network integration does not replace it or add another gas heat solver.

This matters even at exactly zero net gas heat. A hot wall can lose heat while
a cold wall gains the same heat through a gas whose energy remains unchanged.
Replacing that situation with an insulated gas, or returning only net heat,
loses the physical evolution of both walls.

## Existing equations and ownership

At each proposed thermal-network endpoint:

- The original material network resolves every resistor and thermal link.
- Each wall supplies its own endpoint temperature and material conductivity.
- `thermal_ports_candidate` solves the sum of surface heats together with the
  existing gas/piston equations. Each surface uses
  `Q_i = dt * G_i(gas_mid, wall_i) * (T_wall_i - T_gas_mid)`.
- `network_surfaces_candidate` debits each evaluated `Q_i` from its original
  cell, then the existing `circuit_enthalpy` residual closes the network.

The actual per-surface heats are retained even when their sum is zero. The
scalar root's admitted net residual is NOT artificially distributed between
surfaces. Consequently the final combined balance also checks the difference
between net applied gas heat and summed evaluated/debited surface heat.

No circuit charge, flux, thermal enthalpy, gas state, volume or piston momentum
is changed by a trial. A fresh final network evaluation retains all matching
material and boundary records. The candidate returns all existing owners
for publication together; returned Joule and interface heat are already applied.

## Data-only material binding

```rust,ignore
use fs_couple::cavity_phs::piston::thermal::boundary::network::{
    NetworkGasSolve,
    material::GasConductionGeometry,
    surfaces::{NetworkThermalPort, material::{MaterialGasSurface, MaterialSurfacesDrive}},
};

let surfaces = [
    MaterialGasSurface {
        port: NetworkThermalPort { id: 100, cell: hot_wall_cell },
        geometry: GasConductionGeometry {
            area_m2: hot_area, gas_distance_m: gas_length_a, cell_distance_m: wall_length_a,
        },
    },
    MaterialGasSurface {
        port: NetworkThermalPort { id: 200, cell: cold_wall_cell },
        geometry: GasConductionGeometry {
            area_m2: cold_area, gas_distance_m: gas_length_b, cell_distance_m: wall_length_b,
        },
    },
];
let proposal = gas.material_network_surfaces_candidate(
    &network, &thermal, target_s, &circuit, &electrical_ports,
    MaterialSurfacesDrive { surfaces: &surfaces, input: piston_input,
        piston: piston_controls, exchange: net_heat_controls },
    network_gas_controls, &cell_energy_guess_j, || cancelled(),
)?;
let (next_thermal, next_circuit, next_gas, report) = proposal.into_parts();
// report.ports and report.destinations have the same increasing surface-ID order.
// Publish all three histories together. Do NOT add the reported heat again.
```

The adapter consumes the existing `MaterialNetwork`; no second material
registry or handwritten property callback is needed. Every surface uses
`MaterialConductionLayers`, preserving the original material card, phase chart,
property key and selection. Its body conductivity result must equal that of
the same cell in the current network evaluation. Gas conductivity comes from
the common gas midpoint, not a nominal atmosphere.

Distinct surfaces may share a cell when the uniform-cell approximation is
justified. Their heat debits are added to that cell exactly once. Surface IDs
must be unique and are sorted for order-independent evaluation. At most 256
surfaces are admitted; every surface query and all gas/thermal/circuit work
retain explicit bounded iteration and cancellation checks.

The SUM of gas patch areas per cell must fit the exposed area left after its
existing thermal links and reservoir contact. Ambient convection on contacted
cells must be zero so the gas patch is not counted again. The scalar area
check is not a geometric proof of non-overlap: actual non-overlapping patch
placement remains a caller obligation. Missing conductivity refuses even for
zero area, balanced heat, or an otherwise insulated surface.

## Scope and verification

The thermal network uses backward Euler with endpoint material properties;
the gas and moving wall retain their midpoint scheme. The combined scheme is
generally first order. Conservation and residual admission do not certify time
accuracy, uniqueness, passivity or experimental fidelity. Split at physical
source/topology events and refine timesteps independently.

This is one spatially uniform, closed, frozen-composition gas with one
translating inertial wall, fixed circuit topology/C/L and fixed thermal cells.
It does not model mass exchange, combustion, spatial gas temperature gradients,
fluid boundary layers, radiation, mechanical deformation of the thermal walls,
or multiple mechanically linked pistons. Current acoustic properties remain
available on the common gas state; an existing renderer is not automatically
recompiled or restarted.

Four additional Rust integration tests cover balanced finite hot/cold walls,
source records, input-order invariance, conductivity-only circuit/thermal/gas
trajectories, per-cell area sums, surface splitting, missing sources and late
failure/cancellation rollback. Existing upstream gas-port tests are unchanged.
Native execution was unavailable during authoring. Independent Python reference
equations check the same physical balances but are NOT a Rust test pass.

```sh
cargo test -p fs-couple --lib cavity_phs::piston::thermal::boundary:: -- --test-threads=1
```

The underlying series-layer and parallel-path heat balance is described in
MIT Unified Engineering, section 16.4, *Thermal Resistance Circuits*:
https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node118.html
