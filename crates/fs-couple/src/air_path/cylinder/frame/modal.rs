# Complete material frames in shared mechanical assemblies

`MaterialFluidFrame::coupled_component` exports the original frozen frame,
regional added inertia, every fluid-memory momentum, dry elastic storage and
source material viscosity to the existing `CoupledModalSystem`. This is the
complete mechanical component, unlike `modal_parts`, which intentionally leaves
solid laws to a separate constitutive owner. No caloric cells are needed for
frozen-material assembly, and no new time integrator is added.

```rust,ignore
let component = model.coupled_component(
    &initial_complete_state, sample_rate_hz, modal_budget,
    coupling_config, &gate,
)?;
let (placement, mut system) = component.into_system(prepared_budget, &gate)?;
let structural_force = placement.project_nodal_force(&nodal_force, &gate)?;
let mut force = vec![0.0; system.mode_count()];
force[..structural_force.len()].copy_from_slice(&structural_force);
let report = system.step(&force, &gate)?;
placement.reconstruct_into(
    system.components(), &mut nodal_displacement, &mut nodal_velocity, &gate,
)?;
```

For a larger native assembly, `into_parts(component_start, connection_start,
&gate)` returns a placement, ALL native components and ALL connections. Append
them at those exact addresses before adding external connectors. The original
native constructor admits the complete assembly, including prepared propagation,
contact or radiation attachments under their existing contracts. Stepping bodies
separately and exchanging endpoint forces is NOT equivalent to this simultaneous
composition. Material/source objects remain immutably borrowed.

The existing constant full-inertia map `q=L^-T r`, `M=L L^T`, applies to every
spring, material resistance, fluid traction and physical force column. Dry
stiffness is represented by its original positive modal spring factors; native
structural coordinates have zero intrinsic stiffness, so there is no duplicate
elastic energy. Cross-mode added mass and damping remain intact. A nonzero
initial fluid velocity remains present; native fluid position gauges are kept
and counted, not repeatedly reset.

`placement.connections()` identifies disjoint elastic, solid-damping and
per-fluid-region ranges. Elastic work is not heat. Exterior dissipation is not
deposited into finite solid cells. External assembly connectors have their own
separate accounts. These indices can address the native port ledger without
merging regions that happen to share an original structural element.

Use `placement.physical_state()` to recover original dry modal displacement,
loaded momentum and every fluid momentum, or `reconstruct_into()` for nodal
observations. Raw native coordinates are NOT dry-basis velocities. Placement
checks retain complete dimensions, clock and free-coordinate semantics; they do
not authorize changing source connections or swapping a different same-sized
body into the exported addresses. Physical output buffers are scratch on error.
The native owner retains transactional stepping and cancellation.

This remains frozen, small-amplitude linear reference geometry with explicit
retained modes and local transverse circular-cylinder exterior flow. It does
not introduce moving immersion, interacting wakes, fluid heat evolution, phase
changes, measured acoustic pressure or live E(T). The existing thermoelastic
owner remains the appropriate path for its admitted temperature-dependent
source laws. No native throughput improvement or real-time claim is made.

## Focused qualification

```sh
rch exec -- cargo test -p fs-couple --lib air_path::cylinder::frame::modal::component
```

Five new native test groups cover the original complete pHS trajectory and
independent work/loss (prepared and ordinary native paths), exact source/layout
preservation at relocated addresses, dry-frame parity, data-only substitutions,
and failed-step/cancellation/retry preservation. Native Rust compilation,
execution and rustfmt are UNVERIFIED in this executor: Rust/Cargo/DSR/RCH are
absent and the toolchain download host cannot be resolved.

Independent NumPy/SciPy Euler-Bernoulli/Stokes equations passed eight one-/four-
mode cases with 24 continuing held-force intervals each. Maximum normalized
complete-state discrepancy was 8.89e-17; the maximum scaled generator and energy
map discrepancies were 2.13e-16 and 1.20e-16. The cumulative independent
work/heat residual was at most 5.38e-19 J. These are synthetic equation checks,
not native Rust execution, measured-material validation or performance evidence.
