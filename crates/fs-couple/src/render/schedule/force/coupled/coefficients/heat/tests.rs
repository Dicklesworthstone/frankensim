use super::*;
use crate::modal_acoustic_time::{ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel};
use crate::render::schedule::force::coupled::coefficients::rates::MaterialRelaxationTime;
use fs_math::c64::C64;

fn config() -> ModalCouplingConfig {
    ModalCouplingConfig { max_modes: 8, max_connections: 8, max_setup_terms: 1_000_000,
        nyquist_guard_fraction: 0.9, maximum_total_energy_j: 100.0,
        maximum_abs_pressure_pa: 1e6, maximum_abs_connection_force_n: 1e6,
        solve_relative_tolerance: 1e-11, energy_absolute_tolerance_j: 1e-11,
        energy_relative_tolerance: 1e-9 }
}
fn budget() -> MaterialHeatBudget {
    MaterialHeatBudget { max_groups: 8, max_dimension: 256,
        max_dense_bytes: 16*1024*1024, max_multiply_terms: 1_000_000_000 }
}
fn model(q: f64,v: f64,omega: f64) -> ModalAcousticTimeModel {
    let mut m = ModalAcousticTimeModel::try_new_with_free_coordinates(1000,
        vec![ModalAcousticMode { angular_frequency_rad_s: omega, damping_ratio: 0.02,
            pressure_per_modal_velocity: C64::ZERO }],ModalAcousticTimeBudget::audible_reference()).unwrap();
    m.restore_states(&[ModalAcousticState { displacement_m_sqrt_kg:q,velocity_m_sqrt_kg_per_s:v }]).unwrap();
    m
}
fn network(prepared: bool) -> CoupledModalSystem {
    let gate=CancelGate::new_clock_free();
    let preparation=if prepared { ModalPreparation::Prepared(PreparedModalBudget {
        max_modes:8,max_dense_bytes:16*1024*1024,max_multiply_terms:1_000_000_000 }) }
        else { ModalPreparation::Ordinary };
    CoupledModalSystem::new_with_columns_preparation(vec![model(0.01,0.03,20.0),model(-0.02,-0.04,30.0)],
        vec![ModalColumnConnection { column:vec![1.0,-0.7], stiffness_n_m:70.0,
            damping_n_s_m:0.5,rest_extension_m:0.003 }], &[],config(),preparation,&gate).unwrap()
        .with_material_relaxation(vec![
            RelaxationBranch { projection:vec![1.0,0.0,-0.4,0.0],stiffness:100.0,relaxation_time_s:0.02 },
            RelaxationBranch { projection:vec![0.5,0.0,1.0,0.0],stiffness:40.0,relaxation_time_s:0.06 },
            RelaxationBranch { projection:vec![0.0;4],stiffness:30.0,relaxation_time_s:0.04 },
        ],Some(vec![0.007,-0.004,0.02]),&gate).unwrap()
}
fn near(a:f64,b:f64) { assert!((a-b).abs()<=1e-11+1e-8*a.abs().max(b.abs()),"{a:e} != {b:e}"); }

#[test]
fn original_group_heat_matches_complete_native_loss_and_analytic_constrained_history() {
    let gate=CancelGate::new_clock_free();
    for prepared in [false,true] {
        let mut system=network(prepared);
        let groups=system.prepare_material_heat(&[0..1,1..2,2..3,3..3],budget(),&gate).unwrap();
        let together=system.prepare_material_heat(&[0..3],budget(),&gate).unwrap();
        let mut cumulative=0.0;
        for _ in 0..24 {
            let memory=system.material_memory().to_vec();
            let sample=system.samples_rendered();
            let input=[0.02,-0.01];
            let heat=groups.evaluate(&system,&input,&gate).unwrap();
            assert_eq!(system.material_memory(),memory);
            assert_eq!(system.samples_rendered(),sample);
            assert!(heat.iter().all(|h| *h>=0.0));assert_eq!(heat[3],0.0);
            let expected=0.5*memory[2]*memory[2]*(-(-2.0*0.001_f64/0.04).exp_m1());
            near(heat[2],expected);
            near(heat.iter().sum(),together.evaluate(&system,&input,&gate).unwrap()[0]);
            let frame=system.step_under_gate(&input,&gate).unwrap().clone();
            near(heat.iter().sum(),frame.material_dissipation_j);
            assert!(frame.connection_dissipation_j>0.0 && frame.component_dissipation_j>0.0);
            cumulative+=heat[2];
        }
        near(cumulative,0.5*0.02_f64.powi(2)*(-(-2.0*0.024_f64/0.04).exp_m1()));
    }
}

#[test]
fn stale_rate_and_spring_caches_refuse_then_reprepare_without_erasing_history() {
    let gate=CancelGate::new_clock_free();let mut system=network(false);
    let old=system.prepare_material_heat(&[0..3],budget(),&gate).unwrap();
    system.step_under_gate(&[0.02,-0.01],&gate).unwrap();
    assert!(old.is_current(&system,&gate).unwrap());
    let memory=system.material_memory().to_vec();
    system.set_material_relaxation_times(&[MaterialRelaxationTime {branch:1,relaxation_time_s:0.03}],&gate).unwrap();
    assert_eq!(system.material_memory(),memory);
    assert!(!old.is_current(&system,&gate).unwrap());
    assert!(old.evaluate(&system,&[0.02,-0.01],&gate).is_err());
    let fresh=system.prepare_material_heat(&[0..3],budget(),&gate).unwrap();
    let heat=fresh.evaluate(&system,&[0.02,-0.01],&gate).unwrap()[0];
    near(heat,system.step_under_gate(&[0.02,-0.01],&gate).unwrap().material_dissipation_j);
    let mut law=system.connection_coefficients(0).unwrap();law.rest_extension_m=0.004;
    system.set_connection_coefficients(&[law],&gate).unwrap();
    assert!(!fresh.is_current(&system,&gate).unwrap());
}

#[test]
fn overlap_budget_force_and_cancellation_refuse_without_changing_accepted_state() {
    let gate=CancelGate::new_clock_free();let system=network(false);
    let memory=system.material_memory().to_vec();let energy=system.total_energy_j().unwrap();
    for ranges in [vec![0..2,1..3],vec![2..1],vec![0..4],vec![]] {
        assert!(system.prepare_material_heat(&ranges,budget(),&gate).is_err());
    }
    for bad in [MaterialHeatBudget {max_groups:0,..budget()},MaterialHeatBudget {max_dimension:1,..budget()},
        MaterialHeatBudget {max_dense_bytes:1,..budget()},MaterialHeatBudget {max_multiply_terms:1,..budget()}] {
        assert!(system.prepare_material_heat(&[0..3],bad,&gate).is_err());
    }
    let cache=system.prepare_material_heat(&[0..3],budget(),&gate).unwrap();
    assert!(cache.evaluate(&system,&[0.0],&gate).is_err());
    assert!(cache.evaluate(&system,&[f64::NAN,0.0],&gate).is_err());
    let cancelled=CancelGate::new_clock_free();cancelled.request();
    assert!(cache.evaluate(&system,&[0.0;2],&cancelled).is_err());
    assert!(system.prepare_material_heat(&[0..3],budget(),&cancelled).is_err());
    assert_eq!(system.material_memory(),memory);assert_eq!(system.total_energy_j().unwrap(),energy);
    assert_eq!(system.samples_rendered(),0);
    assert!(cache.evaluate(&system,&[0.0;2],&gate).unwrap()[0]>0.0);
}
