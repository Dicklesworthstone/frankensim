//! Independent constitutive participants share one physical state and budget.
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_couple::render::schedule::force::coupled::{
    CoupledModalSystem, ModalColumnConnection, ModalCouplingConfig,
};
use fs_exec::CancelGate;
use fs_math::c64::C64;
use fs_phs::RelaxationBranch;

fn config() -> ModalCouplingConfig {
    ModalCouplingConfig {
        max_modes: 8, max_connections: 8, max_setup_terms: 100_000,
        nyquist_guard_fraction: 0.9, maximum_total_energy_j: 1.0,
        maximum_abs_pressure_pa: 1e6, maximum_abs_connection_force_n: 1e6,
        solve_relative_tolerance: 1e-11, energy_absolute_tolerance_j: 1e-12,
        energy_relative_tolerance: 1e-9,
    }
}
fn model(omega: f64, q: f64) -> ModalAcousticTimeModel {
    let mut m = ModalAcousticTimeModel::try_new(48_000, vec![ModalAcousticMode {
        angular_frequency_rad_s: omega, damping_ratio: 0.01,
        pressure_per_modal_velocity: C64::ZERO,
    }], ModalAcousticTimeBudget::audible_reference()).unwrap();
    m.restore_states(&[ModalAcousticState {
        displacement_m_sqrt_kg: q, velocity_m_sqrt_kg_per_s: 0.0,
    }]).unwrap();
    m
}
fn arm(projection: Vec<f64>, stiffness: f64, tau: f64) -> RelaxationBranch {
    RelaxationBranch { projection, stiffness, relaxation_time_s: tau }
}
fn pair(config: ModalCouplingConfig) -> CoupledModalSystem {
    CoupledModalSystem::new_with_columns(vec![model(400.0, 1e-4), model(700.0, -2e-4)],
        vec![ModalColumnConnection { column: vec![1.0, -0.3], stiffness_n_m: 2e4,
            damping_n_s_m: 0.2, rest_extension_m: 0.0 }], config,
        &CancelGate::new_clock_free()).unwrap()
}
fn one(config: ModalCouplingConfig, q: f64) -> CoupledModalSystem {
    CoupledModalSystem::new(vec![model(100.0, q)], vec![], config,
        &CancelGate::new_clock_free()).unwrap()
}
fn same_state(a: &CoupledModalSystem, b: &CoupledModalSystem) {
    assert_eq!(a.material_memory(), b.material_memory());
    for (x,y) in a.components().iter().zip(b.components()) {
        assert_eq!(x.states(),y.states());
    }
}

#[test]
fn separately_bound_materials_equal_one_complete_coupled_operator() {
    let gate = CancelGate::new_clock_free();
    let a = arm(vec![1.0,0.0,0.0,0.0],9e4,0.001);
    let b = arm(vec![0.2,0.0,-1.0,0.0],5e4,0.004);
    let mut batch = pair(config()).with_material_relaxation(
        vec![a.clone(),b.clone()],Some(vec![0.007,-0.004]),&gate).unwrap();
    let first = pair(config()).with_material_relaxation(
        vec![a],Some(vec![0.007]),&gate).unwrap();
    let mut composed = first.with_additional_material_relaxation(
        vec![b],Some(vec![-0.004]),&gate).unwrap();
    assert_eq!(composed.material_memory(),&[0.007,-0.004]);
    let difference = (composed.linear_frequency_upper_bound_rad_s()
        -batch.linear_frequency_upper_bound_rad_s()).abs();
    assert!(difference < 1e-11);
    same_state(&batch,&composed);
    for tick in 0..32 {
        let force = if tick<16 { [0.02,-0.01] } else { [0.0,0.0] };
        let x = batch.step(&force).unwrap().clone();
        let y = composed.step(&force).unwrap().clone();
        same_state(&batch,&composed);
        assert_eq!(x.material_energy_j,y.material_energy_j);
        assert_eq!(x.material_dissipation_j,y.material_dissipation_j);
        assert_eq!(x.external_work_j,y.external_work_j);
        assert!(y.energy_residual_j.abs()<=y.energy_tolerance_j);
    }
}

#[test]
fn relaxed_addition_keeps_old_unrelaxed_history_and_zero_addition_is_identity() {
    let gate = CancelGate::new_clock_free();
    let first = pair(config()).with_material_relaxation(
        vec![arm(vec![1.0,0.0,0.0,0.0],9e4,0.001)],Some(vec![0.0]),&gate).unwrap();
    let before = first.total_energy_j().unwrap();
    let mut composed = first.with_additional_material_relaxation(
        vec![arm(vec![0.0,0.0,1.0,0.0],4e4,0.003)],None,&gate).unwrap();
    assert_eq!(composed.material_memory(),&[0.0,-0.04]);
    assert_eq!(before,composed.total_energy_j().unwrap());
    let history = composed.material_memory().to_vec();
    let bound = composed.linear_frequency_upper_bound_rad_s();
    composed = composed.with_additional_material_relaxation(vec![],None,&gate).unwrap();
    assert_eq!(history,composed.material_memory());
    assert_eq!(bound,composed.linear_frequency_upper_bound_rad_s());
    composed.step(&[0.0,0.0]).unwrap();
    assert!(composed.with_additional_material_relaxation(vec![],None,&gate).is_err());
}

#[test]
fn composition_admits_the_total_not_each_material_in_isolation() {
    let gate = CancelGate::new_clock_free();
    let mut limits = config(); limits.max_setup_terms = 50;
    let first = pair(limits).with_material_relaxation(
        vec![arm(vec![1.0,0.0,0.0,0.0],1e4,0.01)],None,&gate).unwrap();
    assert!(first.with_additional_material_relaxation(
        vec![arm(vec![0.0,0.0,1.0,0.0],1e4,0.01)],None,&gate).is_err());
    let first = one(config(),0.0).with_material_relaxation(
        vec![arm(vec![1.0,0.0],1e10,0.01)],None,&gate).unwrap();
    assert!(first.with_additional_material_relaxation(
        vec![arm(vec![1.0,0.0],1e10,0.01)],None,&gate).is_err());
    let mut limits = config(); limits.maximum_total_energy_j = 0.008;
    let first = one(limits,1e-4).with_material_relaxation(
        vec![arm(vec![1.0,0.0],1e6,0.01)],Some(vec![0.0]),&gate).unwrap();
    assert!(first.with_additional_material_relaxation(
        vec![arm(vec![1.0,0.0],1e6,0.01)],Some(vec![0.0]),&gate).is_err());
    let first = one(config(),0.0).with_material_relaxation(
        vec![arm(vec![1.0,0.0],1e4,0.01)],None,&gate).unwrap();
    assert!(first.with_additional_material_relaxation(
        vec![arm(vec![1.0,0.0],1e4,0.01)],Some(vec![]),&gate).is_err());
}
