//! Moving velocity/force maps retain the original state and joint work solve.
use fs_couple::modal_acoustic_time::{ModalAcousticMode, ModalAcousticState,
    ModalAcousticTimeBudget, ModalAcousticTimeModel};
use fs_couple::render::schedule::force::coupled::{CoupledModalSystem, ModalAttachment,
    ModalColumnConnection, ModalCouplingConfig};
use fs_couple::render::schedule::force::coupled::contact::{ModalContact, ModalContactConfig};
use fs_couple::render::schedule::force::coupled::contact::friction::{FrictionDrive,
    FrictionModalSystem, ModalFrictionConfig, ModalFrictionPort};
use fs_couple::render::schedule::force::coupled::contact::friction::multiple::{
    MultiFrictionConfig, MultiFrictionModalSystem};
use fs_couple::render::schedule::force::coupled::contact::friction::multiple::surfaces::geometry::{
    FrictionGeometry, FrictionGeometryConfig};
use fs_couple::render::schedule::force::coupled::contact::friction::multiple::obstacles::{
    FrictionObstacleConfig, geometry::ObstacleGeometryConfig};
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;
use fs_math::c64::C64;
use fs_phs::RelaxationBranch;
use fs_tribo::FrictionLaw;

fn budget() -> FrictionGeometryConfig { FrictionGeometryConfig { max_update_terms: 100_000 } }
fn map(component: usize, shape: f64) -> ModalAttachment { ModalAttachment { component, shapes: vec![shape] } }
fn model(omega: f64, q: f64) -> ModalAcousticTimeModel {
    let mut model = ModalAcousticTimeModel::try_new(48000, vec![ModalAcousticMode {
        angular_frequency_rad_s: omega, damping_ratio: 0.01,
        pressure_per_modal_velocity: C64::new(1.0, 0.0),
    }], ModalAcousticTimeBudget::audible_reference()).unwrap();
    model.restore_states(&[ModalAcousticState { displacement_m_sqrt_kg: q,
        velocity_m_sqrt_kg_per_s: 0.01 }]).unwrap(); model
}
fn build(scales: &[f64], memory: bool, stop: u8) -> MultiFrictionModalSystem {
    let gate = CancelGate::new_clock_free();
    let mut network = CoupledModalSystem::new_with_columns(
        vec![model(500.0,1e-4), model(900.0,-2e-4)],
        vec![ModalColumnConnection { column: vec![1.0,-0.3], stiffness_n_m: 4e4,
            damping_n_s_m: 2.0, rest_extension_m: 1e-5 }],
        ModalCouplingConfig { max_modes: 8, max_connections: 8, max_setup_terms: 100_000,
            nyquist_guard_fraction: 0.9, maximum_total_energy_j: 100.0,
            maximum_abs_pressure_pa: 1e6, maximum_abs_connection_force_n: 1e6,
            solve_relative_tolerance: 1e-11, energy_absolute_tolerance_j: 1e-12,
            energy_relative_tolerance: 1e-8 }, &gate).unwrap();
    if memory {
        network = network.with_material_relaxation(vec![RelaxationBranch {
            projection: vec![1.0,0.0,-0.4,0.0], stiffness: 3e4, relaxation_time_s: 0.002,
        }], Some(vec![0.007]), &gate).unwrap();
    }
    let config = ModalFrictionConfig { max_iterations: 96, maximum_force_n: 100.0,
        maximum_slip_speed_m_s: 10.0, force_absolute_tolerance_n: 1e-12,
        force_relative_tolerance: 1e-10, sticking_tolerance_m: 1e-12 };
    let port = |i: usize| ModalFrictionPort { left: map(i,scales[i]), right: None,
        law: FrictionLaw::Hyperbolic { static_mu: 0.8, kinetic_mu: 0.4, characteristic_speed: 0.04 } };
    let core = FrictionModalSystem::new(network,port(0),config,&gate).unwrap();
    let mut joint = core.with_additional_ports((1..scales.len()).map(|i|(port(i),config)).collect(),
        MultiFrictionConfig { max_ports: 4, max_sweeps: 128, max_setup_terms: 100_000 }, &gate).unwrap();
    if stop != 0 {
        let contact = ModalContact { left: map(0,1.0), right: map(0,0.0),
            law: Obstacle::new(vec![-1.0],1,1,vec![0.0],vec![1.0],5000.0,1.0,
                "synthetic relocation test stop".into()).unwrap() };
        joint = joint.with_obstacles(vec![(contact, ModalContactConfig {
            max_iterations: 96, maximum_force_n: 10.0, maximum_penetration_m: 0.01,
            force_absolute_tolerance_n: 1e-12, force_relative_tolerance: 1e-10,
        })], FrictionObstacleConfig { max_contacts: 1, max_setup_terms: 100_000 }, &gate).unwrap();
        if stop == 2 {
            joint = joint.with_obstacle_geometry_motion(ObstacleGeometryConfig { max_update_terms: 100_000 }, &gate).unwrap();
        }
    }
    joint
}
fn geometry(scales: &[f64]) -> Vec<FrictionGeometry> {
    scales.iter().enumerate().map(|(i,&s)|FrictionGeometry { left:map(i,s),right:None }).collect()
}
fn drives(count: usize) -> Vec<FrictionDrive> {
    (0..count).map(|i|FrictionDrive { speed_m_s: if i==0 {0.3} else {-0.2}, normal_force_n:0.01 }).collect()
}
fn states(s: &MultiFrictionModalSystem) -> Vec<ModalAcousticState> {
    s.components().iter().flat_map(|m|m.states().iter().copied()).collect()
}
fn near(a:f64,b:f64) { assert!((a-b).abs()<=2e-11*a.abs().max(b.abs())+2e-14,"{a:e} != {b:e}"); }

#[test]
fn fixed_maps_use_the_original_single_and_joint_paths_bitwise() {
    let gate=CancelGate::new_clock_free();
    for scales in [vec![1.0],vec![1.0,0.8]] {
        let mut ordinary=build(&scales,true,0); let mut held=build(&scales,true,0);
        for _ in 0..32 {
            let a=ordinary.step(&[0.02,-0.01],&drives(scales.len())).unwrap().clone();
            let b=held.step_with_friction_geometry(&[0.02,-0.01],&drives(scales.len()),
                &geometry(&scales),None,None,budget(),&gate).unwrap().clone();
            assert_eq!(a,b); assert_eq!(states(&ordinary),states(&held));
            assert_eq!(ordinary.material_memory(),held.material_memory());
        }
    }
}

#[test]
fn relocated_maps_match_original_construction_with_material_memory_and_normal_reactions() {
    let gate=CancelGate::new_clock_free();
    for stop in [0,1,2] {
        let mut actual=build(&[1.0,0.8],true,stop);
        let mut expected=build(&[0.65,1.1],true,stop);
        assert_eq!(actual.material_memory(),&[0.007]);
        let branches=actual.material_relaxation_branches().to_vec();
        for tick in 0..32 {
            let a=actual.step_with_friction_geometry(&[0.02,-0.01],&drives(2),
                &geometry(&[0.65,1.1]),None,None,budget(),&gate).unwrap().clone();
            let b=expected.step(&[0.02,-0.01],&drives(2)).unwrap().clone();
            near(a.observer_pressure_pa,b.observer_pressure_pa);
            for (x,y) in states(&actual).iter().zip(states(&expected)) {
                near(x.displacement_m_sqrt_kg,y.displacement_m_sqrt_kg);
                near(x.velocity_m_sqrt_kg_per_s,y.velocity_m_sqrt_kg_per_s);
            }
            near(actual.material_memory()[0],expected.material_memory()[0]);
            assert_eq!(actual.material_relaxation_branches(),&branches);
            assert!(a.energy_residual_j.abs()<=a.energy_tolerance_j);
            if stop!=0 && tick==0 { assert!(actual.last_obstacle_frames().unwrap()[0].normal_force_n>0.0); }
        }
    }
}

#[test]
fn failed_trial_restores_port_maps_memory_clock_and_mixed_contact_responses() {
    let gate=CancelGate::new_clock_free();
    for stop in [0,1,2] {
        let mut actual=build(&[1.0,0.8],true,stop); let mut expected=build(&[1.0,0.8],true,stop);
        for _ in 0..8 { actual.step(&[0.0;2],&drives(2)).unwrap(); expected.step(&[0.0;2],&drives(2)).unwrap(); }
        let before=states(&actual); let memory=actual.material_memory().to_vec();
        let frame=actual.last_frame().unwrap().clone();
        // Wrong external dimension fails AFTER the candidate maps/cache swap.
        assert!(actual.step_with_friction_geometry(&[0.0],&drives(2),&geometry(&[0.7,1.2]),
            None,None,budget(),&gate).is_err());
        assert_eq!(states(&actual),before); assert_eq!(actual.material_memory(),memory);
        assert_eq!(actual.last_frame(),Some(&frame)); assert_eq!(actual.samples_rendered(),8);
        assert_eq!(actual.port(0).unwrap().left,map(0,1.0));
        assert_eq!(actual.port(1).unwrap().left,map(1,0.8));
        let a=actual.step(&[0.0;2],&drives(2)).unwrap().clone();
        let b=expected.step(&[0.0;2],&drives(2)).unwrap().clone();
        assert_eq!(a,b); assert_eq!(actual.material_memory(),expected.material_memory());
    }
}

#[test]
fn invalid_geometry_and_work_limits_refuse_before_publication() {
    let gate=CancelGate::new_clock_free(); let mut system=build(&[1.0],true,0);
    let initial=states(&system);
    let mut wrong=geometry(&[0.8]); wrong[0].left.component=1;
    assert!(system.step_with_friction_geometry(&[0.0;2],&drives(1),&wrong,None,None,budget(),&gate).is_err());
    for shapes in [vec![],vec![f64::NAN],vec![0.0]] {
        let g=vec![FrictionGeometry {left:ModalAttachment{component:0,shapes},right:None}];
        assert!(system.step_with_friction_geometry(&[0.0;2],&drives(1),&g,None,None,budget(),&gate).is_err());
    }
    assert!(system.step_with_friction_geometry(&[0.0;2],&drives(1),&geometry(&[0.8]),None,None,
        FrictionGeometryConfig{max_update_terms:1},&gate).is_err());
    assert_eq!(states(&system),initial); assert_eq!(system.samples_rendered(),0);
    assert_eq!(system.material_memory(),&[0.007]);
}
