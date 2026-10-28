use fs_couple::frame_vibration::{FrameVibrationBasis, FrameVibrationBudget, FrameVibrationError};
use fs_couple::modal_acoustic_time::{ModalAcousticState, ModalAcousticTimeBudget};
use fs_solid::frame2d::{Frame2d, FrameMassFormulation, Geometry, SectionModel};

fn budget() -> FrameVibrationBudget {
    FrameVibrationBudget { max_frame_dofs: 12, max_dense_entries: 144 }
}
fn cantilever(axial_fixed: bool) -> Frame2d {
    let mut frame = Frame2d::new(Geometry::Linear);
    frame.add_node(0.0, 0.0);
    frame.add_node(1.0, 0.0);
    frame.fix(0, [true; 3]);
    frame.fix(1, [axial_fixed, false, false]);
    frame.add_mass(0, 3.0);
    frame.add_mass(1, 2.0);
    frame.add_element(0, 1, 3, 0.0, &|| SectionModel::Elastic { ea: 200.0, ei: 8.0 }).unwrap();
    frame
}
fn basis(frame: &Frame2d, modes: usize) -> FrameVibrationBasis {
    FrameVibrationBasis::from_frame(frame, modes, budget(), || false).unwrap()
}
fn close(a: f64, b: f64) {
    assert!((a-b).abs() < 2e-10*(1.0+a.abs()+b.abs()), "{a:e} != {b:e}");
}

#[test]
fn condensation_retains_rotations_and_nodal_moment_work() {
    let frame = cantilever(true);
    let basis = basis(&frame, 1);
    let mode = &basis.modes()[0];
    close(mode.eigenvalue_per_s2, 12.0);
    close(mode.displacement[5]/mode.displacement[4], 1.5);
    let load = [0.0, 0.0, 0.0, 0.0, 4.0, 2.0];
    let mut force = [0.0];
    basis.project_force_into(&load, [0.0, 0.0], &mut force, || false).unwrap();
    close(force[0], 7.0*mode.displacement[4]);
    let state = [ModalAcousticState { displacement_m_sqrt_kg: 0.2, velocity_m_sqrt_kg_per_s: -0.3 }];
    let (mut u, mut v) = ([0.0; 6], [0.0; 6]);
    basis.reconstruct_into(&state, &mut u, &mut v, || false).unwrap();
    close(load.iter().zip(v).map(|(f,v)| f*v).sum(), force[0]*state[0].velocity_m_sqrt_kg_per_s);
    assert_eq!(&u[..3], &[0.0; 3]);
}

#[test]
fn supported_mass_is_not_fabricated_as_modal_participation() {
    let basis = basis(&cantilever(true), 1);
    close(basis.base_influence_mass_kg()[1], 5.0);
    close(basis.retained_effective_mass_kg()[1], 2.0);
    close(basis.retained_effective_mass_kg()[0], 0.0);
}

#[test]
fn exact_shared_time_owner_reproduces_physical_step_response() {
    let basis = basis(&cantilever(true), 1);
    let mut model = basis.time_model(1000, &[0.0], ModalAcousticTimeBudget::audible_reference()).unwrap();
    let mut force = [0.0];
    basis.project_force_into(&[0.0, 0.0, 0.0, 0.0, 4.0, 0.0], [0.0; 2], &mut force, || false).unwrap();
    let dt = 0.125;
    let frame = model.step_duration(&force, dt).unwrap();
    let (mut u, mut v) = ([0.0; 6], [0.0; 6]);
    basis.reconstruct_into(model.states(), &mut u, &mut v, || false).unwrap();
    let omega = 12.0_f64.sqrt();
    close(u[4], 4.0/24.0*(1.0-(omega*dt).cos()));
    close(v[4], 4.0/24.0*omega*(omega*dt).sin());
    close(frame.input_work_j, 4.0*u[4]);
    close(frame.total_modal_energy_j, frame.input_work_j);
    assert_eq!(frame.observer_pressure_pa, 0.0);
}

#[test]
fn base_acceleration_is_relative_inertial_force_not_a_nodal_force() {
    let basis = basis(&cantilever(true), 1);
    let mut force = [0.0];
    basis.project_force_into(&[0.0; 6], [7.0, 0.4], &mut force, || false).unwrap();
    close(force[0], -0.4*basis.modes()[0].base_participation[1]);
    let mut physical = [0.0];
    basis.project_force_into(&[0.0, 0.0, 0.0, 0.0, -0.8, 0.0], [0.0; 2], &mut physical, || false).unwrap();
    close(force[0], physical[0]);
}

#[test]
fn original_consistent_mass_drives_both_base_projection_paths() {
    let mut frame = Frame2d::with_mass_formulation(Geometry::Linear, FrameMassFormulation::ConsistentReference);
    frame.add_node(0.0, 0.0);
    frame.add_node(0.6, 0.8);
    frame.fix(0, [true; 3]);
    frame.add_element(0, 1, 3, 3.0, &|| SectionModel::Elastic { ea: 200.0, ei: 8.0 }).unwrap();
    let basis = basis(&frame, 3);
    let acceleration = [0.2, -0.4];
    let uniform = [0.2, -0.4, 0.0, 0.2, -0.4, 0.0];
    let load: Vec<f64> = frame.inertial_forces(&uniform).unwrap().iter().map(|f| -f).collect();
    let (mut base, mut nodal) = ([0.0; 3], [0.0; 3]);
    basis.project_force_into(&[0.0; 6], acceleration, &mut base, || false).unwrap();
    basis.project_force_into(&load, [0.0; 2], &mut nodal, || false).unwrap();
    for i in 0..3 { close(base[i], nodal[i]); }
}

#[test]
fn explicit_truncation_keeps_inventory_and_source_mode_count() {
    let frame = cantilever(false);
    let basis = basis(&frame, 1);
    assert_eq!(basis.available_mode_count(), 2);
    assert_eq!(basis.modes().len(), 1);
    close(basis.retained_effective_mass_kg()[0], 0.0);
    close(basis.retained_effective_mass_kg()[1], 2.0);
    close(basis.base_influence_mass_kg()[0], 5.0);
    assert!(FrameVibrationBasis::from_frame(&frame, 0, budget(), || false).is_err());
    assert!(FrameVibrationBasis::from_frame(&frame, 3, budget(), || false).is_err());
}

#[test]
fn attachment_keeps_reference_preload_and_the_same_reciprocal_column() {
    let mut frame = cantilever(true);
    frame.static_load(&[0.0, 0.0, 0.0, 0.0, 4.0, 0.0], 1).unwrap();
    let basis = basis(&frame, 1);
    let attachment = basis.translation_attachment(7, 1, [0.0, 1.0]).unwrap();
    assert_eq!(attachment.port.component, 7);
    close(attachment.reference_displacement_m, 4.0/24.0);
    close(attachment.port.shapes[0], basis.modes()[0].displacement[4]);
    let (mut u, mut v) = ([0.0; 6], [0.0; 6]);
    basis.reconstruct_into(&[ModalAcousticState::default()], &mut u, &mut v, || false).unwrap();
    assert_eq!(u.as_slice(), frame.displacements());
    assert!(basis.translation_attachment(0, 2, [0.0, 1.0]).is_err());
    assert!(basis.translation_attachment(0, 1, [0.0, 2.0]).is_err());
}

#[test]
fn shared_nyquist_and_damping_admission_are_not_bypassed() {
    let basis = basis(&cantilever(true), 1);
    let budget = ModalAcousticTimeBudget::audible_reference();
    assert!(basis.time_model(1, &[0.0], budget).is_err());
    assert!(basis.time_model(1000, &[-0.01], budget).is_err());
    assert!(basis.time_model(1000, &[], budget).is_err());
}

#[test]
fn bad_inputs_and_cancellation_do_not_mutate_source_or_runtime() {
    let frame = cantilever(true);
    let before = frame.displacements().to_vec();
    assert!(matches!(FrameVibrationBasis::from_frame(&frame, 1, budget(), || true), Err(FrameVibrationError::Cancelled)));
    assert!(FrameVibrationBasis::from_frame(&frame, 1, FrameVibrationBudget { max_frame_dofs: 6, max_dense_entries: 35 }, || false).is_err());
    let basis = basis(&frame, 1);
    let mut force = [17.0];
    assert!(basis.project_force_into(&[f64::NAN; 6], [0.0; 2], &mut force, || false).is_err());
    assert_eq!(force, [17.0]);
    assert!(matches!(basis.project_force_into(&[0.0; 6], [0.0; 2], &mut force, || true), Err(FrameVibrationError::Cancelled)));
    assert_eq!(frame.displacements(), before.as_slice());
}
