# Material circuit heat drives two gases and one shared mechanism

The shared-piston path now connects to the existing finite circuit/enthalpy
network. This closes the circuit-network gap left by `MATERIAL_SHARED_PISTON.md`:
the two gases do not require separate independently moving walls, and network
cells are not replaced by copied infinite-temperature reservoirs.

## One physical state for every owner

`MaterialDoubleActingGasPiston::material_network_candidate` combines the original
material-bound mechanism with `circuit_enthalpy::material::MaterialNetwork`.
The body's `ResolvedBody` supplies canonical inertia and retains regional density
and geometry sources. `MaterialCircuit` owns charges, fluxes, branch currents and
Joule heat. `DepositionSimulation` owns finite cell enthalpy, phase laws, links,
authored reservoirs and pending events. Each `ClosedGasParcel` owns its original
species, caloric law, EOS, transport and acoustic state.

For every proposed thermal endpoint the original network Newton solve resolves
resistors and thermal conductivity, then the paired gas-boundary solve resolves
both signed heat transfers and the single piston motion. Every proposal starts
from the same accepted histories. Accepted output contains all owners together;
no returned heat should be applied again.

```rust
// After constructing the existing source-backed owners and explicit controls:
let candidate = mechanism.material_network_candidate(
    &materials, &thermal, target_s, &circuit, &electrical_ports,
    drive, controls, contact_geometry, &initial_cell_energy_increments_j,
    || false,
)?;
let (thermal, circuit, mechanism, report) = candidate.into_parts();
```

`drive` is `boundary::double_acting::network::DoubleNetworkDrive`. Its two cell
IDs follow [primary gas, opposing gas] order. It includes the original paired
piston inputs/controls and both heat-exchange bounds. `controls` uses the existing
`boundary::network::NetworkGasSolve`. The two contact geometries use the original
`GasConductionGeometry`; they are thermal patches, not electrical current paths
or alternative piston face areas. Initial guesses are joule increments in sorted
thermal CellId order, not face or resistor order.

Both ports may name one thermal cell. Their signed heats are grouped before one
cell-source debit; their patch areas must fit its combined exposed area. Each
contact reuses the original same-card conductivity resolver. Missing data is an
error even for zero area. Ambient convection on occupied gas-contact cells is
refused rather than counted a second time. Other cells retain native boundaries.

## Independent energy balance

For a cell c, the original residual is

```text
Delta H_c = Joule_c + dt*(external power + net conductive power)_c
            - sum(Q_i for gas ports attached to c)
Q_i = dt G_i(T_cell_endpoint - T_gas_i_midpoint)
```

Both gases share the original reciprocal pressure-work/momentum update. The final
combined check uses one mechanical kinetic-energy change:

```text
Delta E_circuit + sum(Delta H_cells) + Delta U_a + Delta U_b + Delta K
    = electrical port work + external thermal heat + prescribed gas heats
      + additional applied force work
```

The thermal owner evaluates external power from its original rate law and the
pre-transfer source map through `external_step_power_w`. It is not reconstructed
by relabeling missing gas/Joule transfers as external supply. The adapter also
checks actual post-insertion source debits, per-cell residuals, separate gas heat
residuals, all nested budgets and final cancellation.

## Data changes and limits

Changing supported conductivity/resistivity/caloric cards changes current, cell
heat storage, gas pressure, common motion and updated gas acoustic properties
without new solver code. Changing the initial mechanical density card changes
inertia through the original body compiler. These are different initialized
material scenarios, not unaccounted live changes to mass or chemical energy.

Thermal rates use endpoint temperatures with circuit/gas midpoint mechanics, so
the coupled time discretization is generally first order. Tight nonlinear
residuals do not establish time accuracy, uniqueness or unconditional stability.
The integration retains fixed circuit C/L/topology, fixed thermal geometry and
caloric laws, uniform closed frozen-composition gases and one rigid translation.
Split intervals at source/events as required by the existing network. Mechanical
body material state remains frozen; network heating does not silently deform or
melt that body. No combustion, valves, species/mass transfer, conductor-motion
work, moving-grid conduction or emitted audio is implemented here.

The existing Rust piston target contains material-card feedback, same-cell
routing/coverage, numeric/material mechanism parity and atomic refusal tests.
The shared thermal target contains direct native rate-law and cancellation
regressions. Independent Python equation checks also compare the nested model
with a five-variable solve and test timestep refinement; they are synthetic
references, not execution of Rust or empirical material validation. This advances
MR07/MR20/MR43 without claiming those larger beads complete.
