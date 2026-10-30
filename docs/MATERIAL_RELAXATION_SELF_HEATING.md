# Material relaxation heats its original finite layers

`air_path::cylinder::frame::modal::relaxation::temperature::heat::FrameRelaxationThermalSystem`
closes the rate-only feedback loop on the existing complete mechanical owner:

**original layer internal strain -> integrated material heat -> finite layer
enthalpy and temperature -> original source relaxation time -> continuing motion.**

No second mechanical, thermal or material-history integrator is introduced.
The source retains its immutable card, selected claims, reference geometry,
layer inventory and every internal strain. Fluid momentum and foreign bodies'
material histories remain part of the complete reacting assembly.

## Original-source heat, not a sharing rule

`coefficients::heat::MaterialHeatIntegral` integrates each selected group of
original Maxwell source squares `(sqrt(k) b.q - z)^2 / tau` through
`fs_phs::linear_exact::PreparedSourceLoss`. It retains all signed coupling,
spring preload and the actual damping columns. The complete original modal,
connection and material coefficients determine the trajectory. This is not
endpoint power times duration, a residual-energy correction or a fraction of
total damping. Constrained material histories can dissipate even without motion.

The frame adapter forms disjoint groups for each original layer and independent
prefix/interior/suffix groups for foreign material histories. The sum is checked against
the native material-loss account. Only the bound layer's own heat enters its
cell. Foreign material loss, exterior-fluid/connection loss and authored modal
loss remain explicit environment sinks. Numerical dissipation is never heat.

The generic observer is also usable without a thermal model:

```rust,ignore
use fs_couple::render::schedule::force::coupled::coefficients::heat::MaterialHeatBudget;
let heat = mechanical.prepare_material_heat(
    &[core_memory_range, coating_memory_range], observation_budget, &gate,
)?;
let joules = heat.evaluate(&mechanical, &held_native_forces, &gate)?;
// This reads the accepted START state and does not advance or deposit anything.
```

A coefficient signature rejects stale rates, stiffness, damping, preload or
layout. Warm evaluations reuse the cached integral; they use bounded allocated
scratch and O(groups * state_dimension^2) contractions. They do not run an
exponential. The observer's doubled full-state dimension must fit 256, with
explicit aggregate group/byte/work budgets. This does not raise the native
solver's independent limits. External acoustic radiation, DG and nonlinear
contact trajectories are not accepted as if held endpoint forces described
their full paths.

## Bind the actual caloric inventory

```rust,ignore
use fs_couple::air_path::cylinder::frame::modal::relaxation::temperature::heat::
    FrameRelaxationThermalSystem;
use fs_couple::frame_vibration::thermal::FrameLayerCell;

// `material` is the original FrameRelaxationBinding; `mechanical` is the
// complete CoupledModalSystem returned by its material attachment.
let mut coupled = FrameRelaxationThermalSystem::new(
    material, mechanical, &thermal, &original_material_network,
    &[
        FrameLayerCell { element: 0, layer: 0, cell: core_cell },
        FrameLayerCell { element: 0, layer: 1, cell: coating_cell },
    ],
    observation_budget, &gate,
)?;
let report = coupled.step(
    &mut thermal, &held_native_forces, thermal_controls,
    &cell_energy_guesses_j, &gate,
)?;
// Existing FrameFluidPlacement maps these SAME components to physical motion.
placement.reconstruct_into(
    coupled.mechanical().components(), &mut displacement, &mut velocity, &gate,
)?;
```

Each selected layer needs exactly one distinct whole solid caloric cell. The
original card/context, density, physical volume, mass and initial temperature
must match. No artificial heat capacity, heat fraction or missing temperature
is supplied. Unselected source-network cells may conduct heat through the
existing network. The caloric law must remain fully solid with the original
reference density. Whole-layer thermal lumping is an explicit application
assumption, not a spatial-gradient or Biot-number certificate.

## One acceptance boundary and a stated thermal split

An interval advances exact frozen-start-rate mechanics. Its independent layer
heat drives the existing finite-enthalpy/Newton residual, including original
conduction and reservoirs. Every trial uses the original accepted mechanical
history; it does not repeatedly step that history. At the accepted caloric
endpoint, the same pinned claims supply tau(T), and native rate preparation
completes BEFORE joint acceptance. Mechanical states, finite cells, material
snapshots, coefficients and sample clock publish together. Source, budget,
late thermal/native preparation or cancellation refusal publishes none of them.

This is the original first-order caloric/coefficient split. Halve the shared
clock to assess its error independently of structural/fluid spectral refinement.
Endpoint rates act on the NEXT interval. Rate changes preserve stored energy
because stiffness, mass and internal-strain coordinates are unchanged.
Temperature-dependent storage, geometry, expansion, interlayer slip, plasticity
and phase change require other constitutive owners; they are not inferred here.

Changing tau rebuilds the original native propagator. The heat cache rebuilds
before the next attempted step; unchanged numerical rates reuse it. Structural
modes, fluid quadrature and physical histories are never regenerated. No native
throughput or real-time performance improvement is claimed.

## Focused qualification

```sh
rch exec -- cargo test -p fs-couple --lib render::schedule::force::coupled::coefficients::heat
rch exec -- cargo test -p fs-couple --lib air_path::cylinder::frame::modal::relaxation::temperature::heat
```

Three generic heat and three finite-layer regression groups are authored. They
cover ordinary/prepared native material-loss parity, analytic constrained loss,
stale-cache refusal, per-layer self-heating and source-rate updates, continuing
foreign history, missing/duplicate receivers and late caloric rollback/retry.
Native Rust compilation, test execution and rustfmt remain UNEXECUTED in the
editing executor: Rust/Cargo/rustfmt/RCH/DSR are absent and toolchain hosts cannot
be resolved. Source review and patch checks are not native qualification.

Independent NumPy/SciPy reference equations ran six source-group cases over 24
intervals and six layered dry/wet thermal cases at four time resolutions. The
largest heat integral discrepancy versus DOP853 was 2.06e-16 J. The maximum
combined-energy residual in the self-heating references was 1.37e-17 J; the
foreign constrained-history discrepancy was 8.68e-18. For the synthetic wet
beam, the two layers reached approximately 301.761 K and 300.0435 K after twelve
samples from 300 K, using their distinct integrated heat and original masses.
Refinement reduced error against the continuous coupled ODE, approximately
halving the temperature error per halved step in the rate-feedback cases.
These are manufactured equations, NOT execution of the Rust source, measured
material validation or a general physical/temporal error certificate.

## Several heated bodies on the same trajectory

`FrameRelaxationThermalSystem::new_shared` accepts original bindings and their
whole-layer cell destinations in declaration order. The existing `new` delegates
to this path with one binding; there is not a separate single-body thermal solver.
Element and layer numbers remain local to each frame, while receiving `CellId`s
must be distinct across the complete set of bound bodies.

```rust,ignore
let mut coupled = FrameRelaxationThermalSystem::new_shared(
    vec![(receiver_material, receiver_layer_cells),
         (driver_material, driver_layer_cells)],
    mechanical, &thermal, &original_shared_material_network,
    observation_budget, &gate,
)?;
let report = coupled.step(
    &mut thermal, &held_complete_native_forces,
    thermal_controls, &all_cell_energy_guesses_j, &gate,
)?;
let receiver_heat = &report.layers[report.material_layer_ranges[0].clone()];
let driver_heat = &report.layers[report.material_layer_ranges[1].clone()];
```

`materials()` retains all original bindings in that same declaration order;
`material()` remains the first binding for compatibility. `material_layer_ranges()`
and the report's matching ranges disambiguate equal local element/layer IDs.
The common generic fluid source type can be `&dyn FluidTransport` when the
original models borrow different gas/liquid types. No averaged medium replaces
those retained original sources.

Memory blocks may be separated by other bodies' histories and supplied in the
opposite order from native construction. Preparation admits disjoint ranges,
sorts only heat addresses, and independently integrates all unbound gaps.
It never sends a foreign constrained history's loss to a neighboring layer.
Overlapping histories, repeated whole-cell inventory, missing receivers and
insufficient aggregate group budgets refuse before the owner is returned.

Mechanical connections and cross-body thermal conduction retain their original
owners. The source network declares each conductive link's original material
conductivities, endpoints, area and half-cell distances. Changing only the
receiving body's conductivity can change its temperature and relaxation, then
its reaction on the driving body. No heat transfer is inferred from a mechanical
spring, and no mechanical joint is inferred from a thermal link.

Each interval uses one complete mechanical step and one shared caloric residual.
All bodies' endpoint source queries finish before ONE ordered native relaxation-
time replacement. Every mechanical history, source snapshot and thermal cell
publishes together. A late refusal cannot leave the driver at new conditions
and the receiver at old ones. The source-only API
`temperature::shared::set_layer_temperatures` uses the same batch proposal for
prescribed multi-body conditions; the original single-source API delegates to it.

The material heat cache remains aggregate and reuses the original full generator.
One source query per binding per thermal trial is still required; changed times
still rebuild the native operator and invalidate its heat cache. This is not a
measured throughput improvement, an increase in solver capacity, or a stronger
order of thermal integration. Fixed storage/geometry, fully solid caloric laws,
first-order splitting and whole-layer thermal assumptions are unchanged.

Three shared-condition and five shared-heat Rust regression groups are added.
They cover reordered source blocks, native-condition parity, all-body rollback,
disconnected controls, separate layer caloric balances, interior foreign memory,
data-only receiver conductivity and aggregate admission. Existing single-body
regressions exercise the same implementation. Native compilation, execution and
rustfmt remain UNEXECUTED because this editing executor has no Rust/RCH/DSR;
its toolchain download host could not be resolved. Exact source/patch checks do
not replace native qualification.

Independent dimensional equations ran four piecewise-condition cases and six
shared-heating cases at four timestep resolutions. Maximum state difference
against DOP853 in the prescribed cases was 6.81e-16 (normalized), with no boundary
storage jump. The shared-heating maximum combined-energy discrepancy was
2.31e-17 J; foreign internal histories agreed with constrained analytic decay
within 1.83e-17. Changing only the receiver conductivity from zero to 1e-4 W/(m K)
changed the driver's displacement by approximately -1.0664e-8 m in the synthetic
fluid-loaded case. The receiving core warmed from approximately 300.000204 K
to 300.343735 K. Refinement reduced errors against continuous coupled equations.
These deliberately manufactured material/geometry/caloric values are NOT real-
material predictions, native Rust execution or experimental validation.
