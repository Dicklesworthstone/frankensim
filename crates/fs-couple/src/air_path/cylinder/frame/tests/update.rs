//! Continuing physical memory, independent kinetic work and atomic source changes.
use super::*;

fn awake(sim: &MaterialFluidFrame<'_, FluidTransportState>) -> Vec<f64> {
    let motion: Vec<_> = (0..sim.basis().modes().len()).map(|i| ModalAcousticState {
        displacement_m_sqrt_kg: (i+1) as f64*1e-5,
        velocity_m_sqrt_kg_per_s: (i+1) as f64*-2e-3,
    }).collect();
    let fluid: Vec<_> = (0..sim.memory_masses_kg().len())
        .map(|i| (i+1) as f64*3e-5).collect();
    sim.state_from_motion(&motion, &fluid).unwrap()
}

fn exact(a: &[f64], b: &[f64]) {
    assert_eq!(a.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
        b.iter().map(|x| x.to_bits()).collect::<Vec<_>>());
}

#[test]
fn same_source_intervention_preserves_history_and_native_continuation_exactly() {
    let frame = frame(1000.0, 2e6, 500.0, true);
    let mut sim = model(&frame, 4, 1, 400.0, 0.13);
    let state = awake(&sim);
    let before = state.clone();
    let gate = CancelGate::new_clock_free();
    let force = [0.0; 9];
    let expected = sim.step_candidate(&state, &force, 0.001, 1e-9, &gate).unwrap();
    let basis_address = sim.basis() as *const FrameVibrationBasis;
    let update = sim.rebind_fluids(&state, vec![fluid(400.0, 0.13)], 0.0, &gate).unwrap();
    assert!(!update.operator_changed);
    assert_eq!(update.parameter_work_j, 0.0);
    assert_eq!(update.energy_residual_j, 0.0);
    assert_eq!(update.stored_after_j, update.stored_before_j);
    assert_eq!(basis_address, sim.basis() as *const FrameVibrationBasis);
    exact(&state, &before);
    let actual = sim.step_candidate(&state, &force, 0.001, 1e-9, &gate).unwrap();
    exact(&actual.native.x, &expected.native.x);
    exact(&actual.material_port_dissipation_j, &expected.material_port_dissipation_j);
    exact(&actual.fluid_member_dissipation_j, &expected.fluid_member_dissipation_j);
}

#[test]
fn viscosity_changes_continuing_diffusion_without_changing_any_inertia_or_storage() {
    let frame = frame(1000.0, 2e6, 500.0, true);
    let mut sim = model(&frame, 4, 1, 400.0, 0.13);
    let state = awake(&sim);
    let mass = sim.modal_mass_matrix().to_vec();
    let memory = sim.memory_masses_kg().to_vec();
    let q = sim.system().storage().exact_quadratic_matrix().unwrap().to_vec();
    let gate = CancelGate::new_clock_free();
    let previous = sim.step_candidate(&state, &[0.0; 9], 0.001, 1e-9, &gate).unwrap();
    // Non-binary viscosity ratio exposes spurious inertia from rounded c/lambda.
    let update = sim.rebind_fluids(&state, vec![fluid(400.0, 0.217)], 0.0, &gate).unwrap();
    assert!(update.operator_changed);
    assert_eq!(update.parameter_work_j, 0.0);
    assert_eq!(update.energy_residual_j, 0.0);
    assert_eq!(update.stored_after_j, update.stored_before_j);
    exact(&mass, sim.modal_mass_matrix());
    exact(&memory, sim.memory_masses_kg());
    exact(&q, sim.system().storage().exact_quadratic_matrix().unwrap());
    assert_eq!(sim.members()[0].flow().fluid(), &fluid(400.0, 0.217));
    let next = sim.step_candidate(&state, &[0.0; 9], 0.001, 1e-9, &gate).unwrap();
    assert!(next.native.x.iter().zip(&previous.native.x).any(|(a,b)| (a-b).abs() > 1e-12));
    assert!(next.fluid_member_dissipation_j[0] > 0.0);
    assert_ne!(next.fluid_member_dissipation_j, previous.fluid_member_dissipation_j);
    near(next.stored_after_j-next.stored_before_j
        +next.fluid_member_dissipation_j.iter().sum::<f64>()
        +next.material_port_dissipation_j.iter().sum::<f64>(), 0.0, 1e-10);
}

#[test]
fn density_intervention_has_independent_scalar_kinetic_work_and_is_reversible() {
    let frame = frame(1000.0, 2e6, 500.0, false);
    let mut sim = model(&frame, 1, 0, 400.0, 0.13);
    let state = awake(&sim);
    let memory = sim.memory_masses_kg().to_vec();
    // The identical Hermite solid/fluid mass shape gives M=1+rho_fluid/rho_solid.
    let expected = 0.5*state[1]*state[1]*(1.0/1.8-1.0/1.4)
        +state[2..].iter().zip(&memory)
            .map(|(p,m)| 0.5*p*p*(1.0/(2.0*m)-1.0/m)).sum::<f64>();
    let gate = CancelGate::new_clock_free();
    let update = sim.rebind_fluids(&state, vec![fluid(800.0, 0.13)], 1.0, &gate).unwrap();
    near(sim.modal_mass_matrix()[0], 1.8, 2e-10);
    assert!(update.parameter_work_j < 0.0);
    assert!((update.parameter_work_j-expected).abs() < 1e-9*expected.abs());
    for (old, new) in memory.iter().zip(sim.memory_masses_kg()) {
        assert_eq!(*new, 2.0*old);
    }
    near(update.stored_after_j-update.stored_before_j, expected, 1e-12);
    let reverse = sim.rebind_fluids(&state, vec![fluid(400.0, 0.13)], 1.0, &gate).unwrap();
    assert_eq!(reverse.parameter_work_j, -update.parameter_work_j);
    assert_eq!(reverse.stored_after_j, update.stored_before_j);
    exact(&memory, sim.memory_masses_kg());
}

fn two_mode_kinetic(mass: &[f64], p: &[f64]) -> f64 {
    // Independent closed-form inverse, not the assembly's Cholesky factor.
    let determinant = mass[0]*mass[3]-mass[1]*mass[2];
    0.5*(mass[3]*p[0]*p[0]-(mass[1]+mass[2])*p[0]*p[1]+mass[0]*p[1]*p[1])/determinant
}

#[test]
fn density_refresh_keeps_full_cross_mode_mass_and_every_regional_memory() {
    let frame = frame(1000.0, 2e6, 500.0, true);
    let gate = CancelGate::new_clock_free();
    let mut sim = MaterialFluidFrame::new_with_regions(&frame, 2, vec![
        FrameFluidRegion { element: 0, span: [0.1, 0.8], fluid: fluid(400.0, 0.13) },
        FrameFluidRegion { element: 1, span: [0.2, 0.9], fluid: fluid(800.0, 0.21) },
    ], resolution(), frame_budget(), budget(), &gate).unwrap();
    let state = awake(&sim);
    let original_state = state.clone();
    let old_mass = sim.modal_mass_matrix().to_vec();
    let old_memory = sim.memory_masses_kg().to_vec();
    assert!(old_mass[1].abs() > 1e-3);
    let update = sim.rebind_fluids(&state, vec![fluid(1200.0, 0.19), fluid(200.0, 0.33)],
        1.0, &gate).unwrap();
    let expected = two_mode_kinetic(sim.modal_mass_matrix(), &state[2..4])
        -two_mode_kinetic(&old_mass, &state[2..4])
        +state[4..].iter().zip(&old_memory).enumerate().map(|(i,(p,m))| {
            let ratio = if i < 12 { 3.0 } else { 0.25 };
            0.5*p*p*(1.0/(ratio*m)-1.0/m)
        }).sum::<f64>();
    assert!((update.parameter_work_j-expected).abs() < 1e-10*expected.abs());
    exact(&state, &original_state);
    assert_eq!(sim.members()[0].span(), [0.1, 0.8]);
    assert_eq!(sim.members()[1].span(), [0.2, 0.9]);
    assert_eq!(sim.members()[0].flow().fluid(), &fluid(1200.0, 0.19));
    assert_eq!(sim.members()[1].flow().fluid(), &fluid(200.0, 0.33));
    // Independently rebuild the same resolved physical model and compare the
    // actual coupled operator, not independently retuned modal frequencies.
    let rebuilt = MaterialFluidFrame::new_with_regions(&frame, 2, vec![
        FrameFluidRegion { element: 0, span: [0.1, 0.8], fluid: fluid(1200.0, 0.19) },
        FrameFluidRegion { element: 1, span: [0.2, 0.9], fluid: fluid(200.0, 0.33) },
    ], resolution(), frame_budget(), budget(), &gate).unwrap();
    exact(sim.modal_mass_matrix(), rebuilt.modal_mass_matrix());
    exact(sim.memory_masses_kg(), rebuilt.memory_masses_kg());
    exact(sim.system().structure().1, rebuilt.system().structure().1);
    let step = sim.step_candidate(&state, &[0.0; 9], 0.001, 1e-9, &gate).unwrap();
    near(step.stored_after_j-update.stored_before_j
        +step.material_port_dissipation_j.iter().sum::<f64>()
        +step.fluid_member_dissipation_j.iter().sum::<f64>(), update.parameter_work_j, 1e-10);
}

#[test]
fn work_refusal_cancellation_and_invalid_input_preserve_all_owners() {
    let frame = frame(1000.0, 2e6, 500.0, true);
    let mut sim = model(&frame, 4, 1, 400.0, 0.13);
    let state = awake(&sim);
    let q = sim.system().storage().exact_quadratic_matrix().unwrap().to_vec();
    let r = sim.system().structure().1.to_vec();
    let mass = sim.modal_mass_matrix().to_vec();
    let memory = sim.memory_masses_kg().to_vec();
    let gate = CancelGate::new_clock_free();
    assert!(matches!(sim.rebind_fluids(&state, vec![fluid(800.0, 0.217)], 0.0, &gate),
        Err(FrameFluidError::ParameterWork { .. })));
    assert!(sim.rebind_fluids(&state, vec![], 1.0, &gate).is_err());
    assert!(sim.rebind_fluids(&state, vec![fluid(800.0, 0.217)], f64::NAN, &gate).is_err());
    assert!(sim.rebind_fluids(&state, vec![fluid(800.0, 0.217)], -1.0, &gate).is_err());
    assert!(sim.rebind_fluids(&state[..3], vec![fluid(800.0, 0.217)], 1.0, &gate).is_err());
    gate.request();
    assert!(matches!(sim.rebind_fluids(&state, vec![fluid(800.0, 0.217)], 1.0, &gate),
        Err(FrameFluidError::Cancelled)));
    exact(&q, sim.system().storage().exact_quadratic_matrix().unwrap());
    exact(&r, sim.system().structure().1);
    exact(&mass, sim.modal_mass_matrix());
    exact(&memory, sim.memory_masses_kg());
    assert_eq!(sim.members()[0].flow().fluid(), &fluid(400.0, 0.13));
}

#[test]
fn complete_source_snapshot_survives_even_when_mechanical_coefficients_agree() {
    let frame = frame(1000.0, 2e6, 500.0, false);
    let mut sim = model(&frame, 1, 0, 400.0, 0.13);
    let state = awake(&sim);
    let replacement = FluidTransportState::try_new(400.0, 5100.0, 0.67, 0.13).unwrap();
    let update = sim.rebind_fluids(&state, vec![replacement], 0.0,
        &CancelGate::new_clock_free()).unwrap();
    assert!(!update.operator_changed);
    assert_eq!(sim.members()[0].flow().fluid(), &replacement);
    assert_eq!(sim.members()[0].flow().transport(), &replacement);
}

#[test]
fn atmospheric_pressure_reaches_live_frame_inertia_and_rejects_incoherent_gas() {
    use fs_material::gas::{GasSpec, GasState};
    let frame = frame(1000.0, 2e6, 500.0, false);
    let gate = CancelGate::new_clock_free();
    let air = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 100_000.0).unwrap();
    let thin_air = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 70_000.0).unwrap();
    let mut sim = MaterialFluidFrame::new(&frame, 1,
        vec![FrameFluidBinding { element: 0, fluid: air }],
        resolution(), frame_budget(), budget(), &gate).unwrap();
    let state = sim.state_from_motion(&[ModalAcousticState {
        displacement_m_sqrt_kg: 1e-5, velocity_m_sqrt_kg_per_s: 2e-3,
    }], &[3e-4; 12]).unwrap();
    let original_mass = sim.modal_mass_matrix()[0];
    let update = sim.rebind_fluids(&state, vec![thin_air], 1.0, &gate).unwrap();
    assert!(update.operator_changed);
    assert!(update.parameter_work_j > 0.0);
    assert!(sim.modal_mass_matrix()[0] < original_mass);
    assert_eq!(sim.members()[0].flow().fluid(), &thin_air);
    let q = sim.system().storage().exact_quadratic_matrix().unwrap().to_vec();
    let mut incoherent = thin_air;
    incoherent.density *= 2.0;
    assert!(sim.rebind_fluids(&state, vec![incoherent], 1.0, &gate).is_err());
    assert_eq!(sim.members()[0].flow().fluid(), &thin_air);
    exact(&q, sim.system().storage().exact_quadratic_matrix().unwrap());
}

#[test]
fn source_query_cancellation_cannot_publish_a_partial_material_intervention() {
    use std::{cell::Cell, rc::Rc};
    struct ObservedFluid {
        transport: FluidTransportState,
        queries: Rc<Cell<usize>>,
        cancel_on_query: bool,
        gate: Rc<CancelGate>,
    }
    impl FluidTransport for ObservedFluid {
        fn transport_state(&self) -> core::result::Result<FluidTransportState, fs_material::MaterialError> {
            self.queries.set(self.queries.get()+1);
            if self.cancel_on_query { self.gate.request(); }
            Ok(self.transport)
        }
    }
    let frame = frame(1000.0, 2e6, 500.0, false);
    let gate = Rc::new(CancelGate::new_clock_free());
    let original_queries = Rc::new(Cell::new(0));
    let mut sim = MaterialFluidFrame::new(&frame, 1,
        vec![FrameFluidBinding { element: 0, fluid: ObservedFluid {
            transport: fluid(400.0, 0.13), queries: original_queries.clone(),
            cancel_on_query: false, gate: gate.clone(),
        } }], resolution(), frame_budget(), budget(), &gate).unwrap();
    let state = sim.state_from_motion(&[ModalAcousticState::default()], &[3e-4; 12]).unwrap();
    let q = sim.system().storage().exact_quadratic_matrix().unwrap().to_vec();
    let replacement_queries = Rc::new(Cell::new(0));
    assert!(matches!(sim.rebind_fluids(&state, vec![ObservedFluid {
        transport: fluid(800.0, 0.217), queries: replacement_queries.clone(),
        cancel_on_query: true, gate: gate.clone(),
    }], 1.0, &gate), Err(FrameFluidError::Cancelled)));
    assert_eq!(original_queries.get(), 1);
    assert_eq!(replacement_queries.get(), 1);
    assert_eq!(sim.members()[0].flow().transport(), &fluid(400.0, 0.13));
    exact(&q, sim.system().storage().exact_quadratic_matrix().unwrap());
}
