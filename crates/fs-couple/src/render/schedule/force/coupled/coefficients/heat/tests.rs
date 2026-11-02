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
fn staged_contact_input_uses_original_group_integrals_without_replaying_or_double_preload() {
    let gate = CancelGate::new_clock_free();
    for prepared in [false, true] {
        let mut system = network(prepared);
        let groups = system.prepare_material_heat(&[0..1, 1..3], budget(), &gate).unwrap();
        let memory = system.material_memory().to_vec();
        let states = system.components().iter().flat_map(|m| m.states().iter().copied()).collect::<Vec<_>>();
        assert!(groups.evaluate_staged(&system, &gate).is_err());
        // The converged held contact/support effort differs from the requested
        // free-body force. Nonzero spring preload must enter exactly once.
        let held = [-0.3, 0.6];
        let expected = groups.evaluate(&system, &held, &gate).unwrap();
        let free = groups.evaluate(&system, &[0.0; 2], &gate).unwrap();
        system.stage_inner(&held, Some(&gate)).unwrap();
        let staged = groups.evaluate_staged(&system, &gate).unwrap();
        for (&actual, &expected) in staged.iter().zip(&expected) { near(actual, expected); }
        near(staged.iter().sum(), system.staged_frame.material_dissipation_j);
        assert!(staged.iter().zip(&free).any(|(a,b)| (a-b).abs() > 1e-12));
        assert_eq!(system.material_memory(), memory);
        assert_eq!(system.components().iter().flat_map(|m| m.states().iter().copied()).collect::<Vec<_>>(), states);
        assert_eq!(system.samples_rendered(), 0);
        system.stage_inner(&held, Some(&gate)).unwrap();
        assert_eq!(groups.evaluate_staged(&system, &gate).unwrap(), staged);
        // A refused trial cannot leave a previous candidate's heat observable.
        assert!(system.stage_inner(&[f64::NAN, 0.0], Some(&gate)).is_err());
        assert!(groups.evaluate_staged(&system, &gate).is_err());
        system.stage_inner(&held, Some(&gate)).unwrap();
        system.publish_staged();
        assert_eq!(system.samples_rendered(), 1);
        assert!(groups.evaluate_staged(&system, &gate).is_err());
    }
}

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
fn material_viscosity_and_memory_match_separate_native_loss_accounts() {
    let gate = CancelGate::new_clock_free();
    // The fixture's independent native component resistances are 2*.02*20=.8
    // and 2*.02*30=1.2. Two physical regions contribute .3 and .5 to the first;
    // a third owns the second. The joint dashpot remains outside material heat.
    let sources = [
        MaterialViscousHeatSource::from_material_region(0, 0, vec![(0, 0.3)]),
        MaterialViscousHeatSource::from_material_region(0, 1, vec![(0, 0.5)]),
        MaterialViscousHeatSource::from_material_region(1, 0, vec![(1, 1.2)]),
    ];
    for prepared in [false, true] {
        let mut system = network(prepared);
        let mut channel = system.prepare_material_heat_with_viscosity(
            &[0..1, 1..3], &sources, budget(), &gate).unwrap();
        for sample in 0..3 {
            if sample == 1 {
                let memory = system.material_memory().to_vec();
                let rates = [MaterialRelaxationTime { branch: 1, relaxation_time_s: 0.03 }];
                system.set_material_relaxation_times(&rates, &gate).unwrap();
                channel = channel.prepare_after_material_rates(&system, &rates, &gate).unwrap();
                assert_eq!(system.material_memory(), memory);
            }
            let input = [-0.3, 0.6];
            let memory = system.material_memory().to_vec();
            let energy = system.total_energy_j().unwrap();
            let heat = channel.evaluate(&system, &input, &gate).unwrap();
            assert_eq!(heat.len(), 5);
            assert_eq!(system.material_memory(), memory);
            assert_eq!(system.total_energy_j().unwrap(), energy);
            system.stage_inner(&input, Some(&gate)).unwrap();
            assert_eq!(channel.evaluate_staged(&system, &gate).unwrap(), heat);
            near(heat[..2].iter().sum(), system.staged_frame.material_dissipation_j);
            near(heat[2..].iter().sum(), system.staged_frame.component_dissipation_j);
            assert!(heat.iter().all(|h| *h > 0.0));
            assert!(system.staged_frame.connection_dissipation_j > 0.0);
            system.publish_staged();
        }
        assert!(system.admit_material_viscosity([&sources[0], &sources[0]], &gate).is_err());
        let excess = MaterialViscousHeatSource::from_material_region(0, 2, vec![(0, 0.1)]);
        assert!(system.admit_material_viscosity(sources.iter().chain([&excess]), &gate).is_err());
        let other_component = MaterialViscousHeatSource::from_material_region(0, 3, vec![(1, 0.1)]);
        assert!(system.admit_material_viscosity([&other_component], &gate).is_err());
        // Pure Kelvin bodies require no fictitious hereditary branch or state.
        let only_viscosity = system.prepare_material_heat_with_viscosity(
            &[], &sources, budget(), &gate).unwrap();
        let heat = only_viscosity.evaluate(&system, &[0.0; 2], &gate).unwrap();
        near(heat.iter().sum(), system.step_under_gate(&[0.0; 2], &gate).unwrap().component_dissipation_j);
    }
    let mut kelvin = CoupledModalSystem::new(vec![model(0.01, 0.03, 20.0)],
        vec![], config(), &gate).unwrap();
    let source = kelvin.prepare_material_heat_with_viscosity(&[], &sources[..2], budget(), &gate).unwrap();
    assert!(kelvin.material_memory().is_empty());
    kelvin.stage_inner(&[0.0], Some(&gate)).unwrap();
    let heat: f64 = source.evaluate_staged(&kelvin, &gate).unwrap().iter().sum();
    // Independent closed-form damped oscillator energy, with no material
    // history, force input, connection or alternative integrator in the oracle.
    let (omega, decay, time) = (20.0_f64, 0.4_f64, 0.001_f64);
    let wd = (omega*omega - decay*decay).sqrt();
    let (sin, cos) = (wd*time).sin_cos();
    let q = (-decay*time).exp() * (0.01*cos + (0.03+decay*0.01)/wd*sin);
    let v = (-decay*time).exp() * (0.03*cos - (decay*0.03+omega*omega*0.01)/wd*sin);
    near(heat, 0.5*(0.03_f64.powi(2)+(omega*0.01).powi(2)-v*v-(omega*q).powi(2)));
    near(heat, kelvin.staged_frame.component_dissipation_j);
    assert_eq!(kelvin.staged_frame.material_dissipation_j, 0.0);
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
