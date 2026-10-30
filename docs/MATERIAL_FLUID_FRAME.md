# Material-defined frames in a gas or liquid

`fs_couple::air_path::cylinder::frame::MaterialFluidFrame` connects the existing
source-defined `fs_solid::frame2d::Frame2d` to the existing causal cylinder-flow
operator. A bound member's original material section supplies its outer radius;
its density, elastic response and optional solid viscosities continue to drive
the original dry frame. The exterior source supplies the same checked
`FluidTransport` view used by other fluid consumers. There are no material-name
branches, independent mass/stiffness overrides, synthetic gas records or fitted
per-mode fluid damping coefficients.

## Use the existing source objects

Prepare the frame with `add_material_element` or
`add_material_viscous_element`, then provide `FrameFluidBinding { element, fluid }`
for each wholly immersed circular member. The fluid can be an existing resolved
liquid, resolved gas/mixture, moist-air owner, or an explicitly authored transport
tuple. Its source object remains inspectable through `members()[i].flow().fluid()`.
Different fluid source types can be used in different preparations; each binding
vector has one Rust source type. The original frame is immutably borrowed.

```rust,ignore
use fs_couple::air_path::cylinder::frame::{
    FrameFluidBinding, FrameFluidBudget, MaterialFluidFrame,
};

let coupled = MaterialFluidFrame::new(
    &frame,
    retained_modes,
    vec![FrameFluidBinding { element: 0, fluid: resolved_exterior }],
    cylinder_resolution,
    dry_frame_budget,
    FrameFluidBudget {
        max_states: 512,
        max_dense_entries: 512 * 512,
        max_projection_terms: 100_000,
        mass_inverse_tolerance: 1e-10,
    },
    &cancel_gate,
)?;
let mut state = coupled.state_from_motion(&initial_modal_motion, &initial_fluid_velocities)?;
let proposal = coupled.step_candidate(
    &state, &incremental_nodal_forces_and_moments, duration_s,
    energy_tolerance_j, &cancel_gate,
)?;
// Publish only after any enclosing transaction has accepted this proposal.
state = proposal.native.x;
let modal_motion = coupled.modal_states(&state)?;
coupled.basis().reconstruct_into(
    &modal_motion, &mut nodal_displacement, &mut nodal_velocity, || false,
)?;
```

The complete state count is `2*retained_modes + 4*loaded_members*intervals`.
Every fluid velocity must be supplied explicitly at initialization, in member,
Gauss-point, then spectral-branch order. Retain the complete returned state on
subsequent steps; reinitializing fluid velocities would erase physical history.
The original basis converts loaded momentum back to velocity with the full
added-mass inverse, rather than treating momentum as the old dry modal velocity.

## Physical account and boundaries

Four positive Gauss points integrate the original cubic-Hermite transverse
traces. The full projected added inertia retains cross-mode terms. Each viscous
fluid-memory branch has its own momentum and equal/opposite resistance forces.
Original axial/bending material viscosity retains its full projected matrix.
The native `fs-phs` solver advances the combined system; this adds no integrator.

`FrameFluidStep` separates original material-port losses from each immersed
member's fluid losses and checks them against stored-energy change and actual
nodal work. None of those losses has yet been deposited into a thermal model.
Budgets include all fluid memory before allocation. Duplicate loading, absent
source sections, invalid source clocks, failed work gates and cancellation
refuse without mutating the source frame or supplied state. Native dense solves
are synchronous; polling around them does not promise bounded inner-solve latency.

This is frozen, small-amplitude planar vibration with stationary supports and
stationary far-field Newtonian fluid. Members are prismatic circular sections;
only exterior transverse local Stokes loading is included. End effects, axial
and rotary drag, interacting wakes, confinement, cavitation, thermal evolution,
internal fluid in hollow members, moving supports and acoustic radiation remain
outside this model. Modal and fluid-spectrum truncation must be refined; a
small energy residual does not establish trajectory or experimental accuracy.
Empty fluid bindings retain the original dry Galerkin model and material loss.

## Focused qualification

```sh
rch exec -- cargo test -p fs-solid --lib transverse_trace_tests
rch exec -- cargo test -p fs-couple --lib air_path::cylinder::frame
```

Five trace regressions and six coupled-frame regressions are authored. Rust,
Cargo, DSR and RCH were unavailable in the editing environment, so these native
tests, compilation and formatting are UNEXECUTED. Independent NumPy/SciPy
reference equations passed nine rotated/scaled trace cases, partial-immersion
cross-mass checks, complex impedance elimination, a 100-step work balance and
three time-refinement levels against a continuous matrix exponential. Synthetic
data-only solid-density/modulus and fluid-density/viscosity changes altered the
reference displacement. Those equation checks are not execution of the Rust
implementation, cross-platform qualification or measured-material validation.
