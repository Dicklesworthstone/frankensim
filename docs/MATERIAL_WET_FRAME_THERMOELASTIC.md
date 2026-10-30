# Thermoelastic frames with the original fluid memory

`fs_couple::frame_vibration::thermoelastic::fluid::FrameFluidThermoelasticSystem`
connects the existing `MaterialFluidFrame` to the original whole-layer
thermoelastic, finite-enthalpy and optional electrical transaction. Temperature
can now change the solid's Young's modulus, integrated axial free expansion and
viscosity while the surrounding fluid's full added inertia and viscous-memory
momenta remain in the same mechanical trajectory.

The constant loaded inertia is factored once: `M = L L^T`, `q = L^-T r`.
Physical force and observation columns use the same `L^-1` map. This retains
cross-mode terms; it is not a frequency fit or diagonal added-mass approximation.
The original source factory first checks the complete structural basis and
reconstructs stiffness from its layer laws. Those laws and all material ports
then undergo the same congruence. Solid elasticity is included exactly once.

```rust,ignore
use fs_couple::frame_vibration::thermoelastic::fluid::{
    FrameFluidThermoelasticConfig, FrameFluidThermoelasticSystem,
};

let mut owner = FrameFluidThermoelasticSystem::new(
    &fluid_frame, &initial_complete_state, &thermal, &material_network,
    &layer_cells, preparation, &cancel_gate,
)?;
let report = owner.step(
    &mut thermal, &held_nodal_forces_and_moments,
    thermal_controls, &cell_energy_guesses_j, &cancel_gate,
)?;
owner.reconstruct_into(&mut displacement, &mut velocity, &cancel_gate)?;
```

`preparation` is an explicit `FrameFluidThermoelasticConfig`: sample rate,
modal/coupling/prepared budgets, projection work cap, physical strain envelope
and independent initial-energy tolerance. The source frame must retain ALL
structural modes, every member's material viscosity (explicit zero is valid),
zero frozen preload, and its whole-layer caloric cells. Solid reference density
is fixed; E(T) is scalar or globally affine, and thermal expansion uses the
existing scalar source law. Missing data is not treated as zero.

Use `modal_states()` for original dry-basis coordinates, or `reconstruct_into`
for physical nodal observations. `mechanical()` exposes the actual normalized
native state, which MUST NOT be passed directly to the dry basis. The original
`MaterialFluidFrame` remains immutable reference data, not the live stiffness.
`physical_state()` preserves its canonical motion/fluid-momentum layout for
inspection, but resuming it with that old frozen stiffness would change physics.

Each spectral fluid momentum gains a native free position gauge with no spring,
observation or external force. These gauges are retained and counted: there are
`m + k` native modal coordinates for `m` structural modes and `k` fluid memories,
not just `m` modes. All source ports and native setup terms are budgeted before
allocation. Physical force projection uses reusable scratch. Structural spectra,
loaded-inertia factors and cylinder quadrature are not rebuilt per thermal trial.
Native propagator/coefficient updates still have their own costs; no runtime
speedup or real-time qualification is claimed.

The ORIGINAL thermoelastic owner advances exact frozen-start-temperature
mechanics, then its finite-cell residual accounts for source `Psi + T*S`,
reversible heat, material dissipation and endpoint constitutive replacement.
This remains a first-order thermal split. Stiffness/expansion changes are not
free energy injections, and conservative spring work is not frictional heat.
Exterior-fluid loss remains an explicit environment sink, never duplicate solid
heating or an invented finite fluid bath. Failure/cancellation publishes neither
thermal state, mechanical state, source coefficients nor sample clock.

`with_material_circuit` and `step_electrical` reuse the original same-cell R(T),
finite electrical storage and Joule-heating owner. Its source network must borrow
the original cell/link slices. Electrical heat can drive expansion/softening in
this same fluid-loaded history. Fixed electrical geometry/C/L remain explicit;
no motional EMF, magnetic force or evolving wire dimensions are inferred here.

Reference geometry and exterior transport are fixed. This is small-strain
linear-reference thermoelasticity with local circular-member Stokes loading,
not geometric buckling, moving-boundary CFD, variable-density inertia, radial
fluid-domain expansion, melting, interacting wakes or acoustic radiation.
Thermal lumping and modal/spectral/time refinement remain application obligations.

## Focused checks

```sh
rch exec -- cargo test -p fs-couple --lib air_path::cylinder::frame::modal
rch exec -- cargo test -p fs-couple --lib frame_vibration::material::thermal::thermoelastic::fluid
```

Two coordinate-map and five coupled regressions are authored: source E/alpha
feedback, independent internal-energy and nodal-work checks, retained fluid
history, dry-owner parity, late failure/retry and original electrical heating.
Native Rust compilation, test execution and formatting remain UNEXECUTED in the
editing environment (Rust/Cargo/RCH/DSR unavailable).

Independent NumPy/SciPy equations passed eight coordinate cases and eight
wet/dry thermoelastic cases. Maximum coordinate state discrepancy was 2.78e-17;
maximum cumulative combined-energy discrepancy was 7.65e-14 J. Three time
refinements against continuous ODEs gave approximately 2x error reduction per
halving, consistent with the first-order split. The electrical reference uses
equivalent prescribed Joule power, not execution of the Kirchhoff solver.
These synthetic references do not compile/run Rust or validate measured alloys.
