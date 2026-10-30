# Prepared source-viscosity feedback for material frames

`MaterialFluidFrame::prepare` retains the original full Hamiltonian, loaded
modal inertia, source sections, and every fluid-memory momentum. It uses the
existing native `StepWorkspace` with analytic Hessian/resistance actions and
exact independent-memory condensation. No modes or fluid coordinates are
removed. The factored resistance callback evaluates each station velocity
once, then projects its accumulated reaction once. Native algebra retains
its cancellation, residual, iteration and safe full-LU fallback gates.

Preparation allocates. Native stepping and its physical callbacks allocate
nothing after preparation; returned reports and material evaluations can
allocate. No elapsed-time, real-time, Rust benchmark or full-workspace
qualification is claimed. Complete source-frame state/dense limits still
apply, and the native solver additionally owns its bounded workspace.

`LiquidViscosityFrame::new` binds the exact original `LiquidTemperatureLaw`
for every loaded member. `step_at_temperatures` evaluates those pinned claims
at prescribed temperatures and changes the same physical resistances.
A validity-wide scalar constant-density source is mandatory: matching two
density samples is not enough. At fixed density/radius and dimensionless
spectrum, c_j/lambda_j is viscosity-independent. The original kinetic storage
and all momenta therefore remain unchanged when viscosity varies.
Dynamic and kinematic source-viscosity conventions use the existing material
owner; unsupported source domains and unowned rate/frequency laws refuse.

```rust,ignore
let mut prepared = model.prepare(&gate)?;
let next = prepared.step_candidate(&state, &nodal_forces, dt, energy_budget_j, &gate)?;
state = next.state; // only after the enclosing transaction admits the proposal

let mut viscous = LiquidViscosityFrame::new(&model, original_liquid_laws, &gate)?;
let next = viscous.step_at_temperatures(
    &state, &nodal_forces, &prescribed_temperatures_k, dt, energy_budget_j, &gate,
)?;
```

Prescribed temperature is NOT a solved finite thermal history. In particular,
assigning every unconfined Stokes loss to an arbitrary finite liquid volume
would not resolve spatial fluid heating. This path makes no such claim.
Variable density, thermal expansion, phase change, finite boundaries and
spatially varying fluid temperature require their corresponding physical owners.

Five native Rust regressions cover dense/condensed parity, source-only
viscosity changes, cold-hot operator comparison, density-source refusal and
retry after source/cancellation failures. They are authored but UNEXECUTED:
Rust, Cargo, RCH and DSR are unavailable in this editing environment.
Independent NumPy/SciPy matrix equations passed ten one-/four-mode cases over
viscosity ratios 0.01 through 100. Maximum normalized resistance-action
mismatch was 6.66e-17, dense/Schur state discrepancy 2.78e-17, and work residual
6.35e-17 J. These reference calculations do not execute or compile Rust.

## Joint finite-solid heating with the complete fluid trajectory

`frame_vibration::thermal::fluid::FrameFluidThermalSystem` joins that prepared
mechanical owner to the ORIGINAL finite-cell enthalpy/Newton solver. It reuses
the existing `FrameLayerCell` and whole-layer source binding, including original
card and caloric-law identity, mass, volume, phase and physical query context.
All solid density and elastic storage coefficients must be source-declared
constants; no thermal expansion or variable-stiffness work is silently omitted.

Every thermal trial starts from the same accepted mechanical state, ALL fluid
momenta and accepted cell enthalpy. The original section temperature laws
resolve fresh axial/bending viscosities at each proposed layer temperature.
Those viscosities determine the original physical port resistances in the SAME
mechanical interval. Actual discrete-effort port losses are then partitioned
by the fresh per-layer eta*A and eta*I integrals and deposited once into their
actual caloric cells. Zero reference viscosity can become positive; the code
uses actual source resistance, never a division by its old value.

The native thermal network supplies conduction, boundary sources and cooling.
It retains all existing event-order, source, nonlinear-budget and entropy
rules. The joint energy gate includes mechanical/fluids' stored energy, finite
solid enthalpy, actual nodal work, external thermal supply and the explicit
exported fluid loss. The returned entropy account covers the finite thermal
network only, not a global second-law claim for the external fluid.

The unconfined exterior remains at its original, externally held source state.
Its viscous loss is exported to that environment, NOT deposited into the solids
or an invented finite fluid bath. This owner and the prescribed liquid-temperature
path above are distinct: it does not solve the exterior's spatial temperature.

```rust,ignore
use fs_couple::frame_vibration::thermal::fluid::FrameFluidThermalSystem;

let mut coupled = FrameFluidThermalSystem::new(
    &model, initial_complete_state, &thermal, &material_network,
    &original_layer_cells, &gate,
)?;
let report = coupled.step(
    &mut thermal, target_s, &incremental_nodal_forces_and_moments,
    thermal_controls, &initial_cell_energy_increments_j, &gate,
)?;
// Both histories advanced together; on Err, both retain their accepted state.
let motion = model.modal_states(coupled.state())?;
```

The cell model assumes uniform temperature in each physical layer; its thermal
lumping must be justified independently for the application. Endpoint viscosity
with backward-Euler thermal evolution is first order in general. Refine time,
structural modes and fluid spectrum separately, and cut intervals at source or
load discontinuities. Prepared trials reuse the full storage, structural modes,
loaded inverse mass, fluid quadrature and native scratch; they do not rebuild
these on every Newton evaluation. Returned reports and source queries allocate.

Six additional native regressions cover complete heat/work closure, data-only
thermal-viscosity feedback in actual motion, independently heated concentric
layers, late source/budget refusal and valid retry, whole-layer inventory refusal,
and zero-to-positive source viscosity. Run through the project build owner:

```sh
rch exec -- cargo test -p fs-couple --lib air_path::cylinder::frame::viscosity
rch exec -- cargo test -p fs-couple --lib frame_vibration::material::thermal::fluid
```

These Rust tests and formatting remain UNEXECUTED in this editing environment.
Independent physical-coordinate NumPy/SciPy calculations passed six thermal
cases, a two-layer fresh-heat-partition case and three refinements against a
continuous ODE reference. Maximum one-step combined energy discrepancy was
4.77e-13 J. The two-layer core heat share changed from 0.25 to about 0.37445;
using the cold share would route heat incorrectly. Halving the thermal step
reduced the observed temperature error by about two, consistent with the stated
first-order thermal composition. These are equation checks, not execution of
Rust, empirical accuracy estimates, performance measurements or workspace
qualification. Variable elastic energy/density, phase transitions, changing
geometry and spatial exterior heating remain outside this adapter.
