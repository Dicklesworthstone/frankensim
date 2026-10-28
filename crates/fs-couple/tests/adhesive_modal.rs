//! Frozen adhesive patches exchange physical work with the native finite bodies.
use fs_couple::modal_acoustic_time::{ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel};
use fs_couple::render::schedule::force::coupled::{CoupledModalSystem, ModalAttachment, ModalCouplingConfig};
use fs_couple::render::schedule::force::coupled::contact::friction::{FrictionDrive,FrictionModalSystem,ModalFrictionConfig,ModalFrictionRegime};
use fs_couple::render::schedule::force::coupled::contact::friction::adhesive::AdhesiveModalBinding;
use fs_exec::CancelGate;
use fs_math::c64::C64;
use fs_tribo::{FrictionLaw};
use fs_tribo::adhesive::{AdhesiveContactState,AdhesiveFrictionLaw,AdhesiveTangentialContact};
const MASS:[f64;2]=[0.002,0.01];
const DT:f64=1.0/48000.0;
fn mechanics(initial_speed:f64, area:f64)->FrictionModalSystem {
    mechanics_with_pressure(initial_speed, area, 0.0,
        FrictionLaw::Coulomb {static_mu:0.0,kinetic_mu:0.0})
}
fn mechanics_with_pressure(initial_speed:f64, area:f64, normal:f64, pressure:FrictionLaw)->FrictionModalSystem {
    let models=MASS.iter().enumerate().map(|(i,m)| {
        let mut model=ModalAcousticTimeModel::try_new_with_free_coordinates(48000,
            vec![ModalAcousticMode {angular_frequency_rad_s:0.0,damping_ratio:0.0,pressure_per_modal_velocity:C64::ZERO}],
            ModalAcousticTimeBudget::audible_reference()).unwrap();
        model.restore_states(&[ModalAcousticState { displacement_m_sqrt_kg:0.0,
            velocity_m_sqrt_kg_per_s:if i==0 {initial_speed*m.sqrt()} else {0.0}}]).unwrap(); model
    }).collect();
    let network=CoupledModalSystem::new(models,vec![],ModalCouplingConfig {
        max_modes:2,max_connections:1,max_setup_terms:4096,nyquist_guard_fraction:0.9,
        maximum_total_energy_j:1.0,maximum_abs_pressure_pa:1e6,maximum_abs_connection_force_n:100.0,
        solve_relative_tolerance:1e-11,energy_absolute_tolerance_j:1e-12,energy_relative_tolerance:1e-9},&CancelGate::new()).unwrap();
    let attachment=|i:usize|ModalAttachment { component:i,shapes:vec![1.0/MASS[i].sqrt()]};
    let law=AdhesiveFrictionLaw::new(pressure,1e6,0.5e6).unwrap();
    let binding=AdhesiveModalBinding::new(AdhesiveTangentialContact::new(law,AdhesiveContactState::new(area,normal).unwrap()).unwrap()).unwrap();
    FrictionModalSystem::new_with_adhesion(network,attachment(0),Some(attachment(1)),binding,
        ModalFrictionConfig {max_iterations:80,maximum_force_n:100.0,maximum_slip_speed_m_s:10.0,
            force_absolute_tolerance_n:1e-12,force_relative_tolerance:1e-11,sticking_tolerance_m:1e-12},&CancelGate::new()).unwrap()
}
fn drive()->FrictionDrive {FrictionDrive {speed_m_s:0.0,normal_force_n:0.0}}
fn states(s:&FrictionModalSystem)->Vec<ModalAcousticState> {s.components().iter().flat_map(|m|m.states().iter().copied()).collect()}
fn momentum(s:&FrictionModalSystem)->f64 {s.components().iter().enumerate().map(|(i,m)|m.states()[0].velocity_m_sqrt_kg_per_s*MASS[i].sqrt()).sum()}
fn close(a:f64,b:f64,t:f64) {assert!((a-b).abs()<=t,"{a:e} != {b:e}");}
#[test]
fn adhesive_static_reaction_accelerates_both_real_bodies_without_heat() {
    let mut s=mechanics(0.0,1e-6);let f=0.1;
    let frame=*s.step(&[f/MASS[0].sqrt(),0.0],drive()).unwrap();
    assert_eq!(frame.regime,ModalFrictionRegime::Sticking);
    close(frame.traction_n,-f*MASS[1]/(MASS[0]+MASS[1]),1e-12);
    assert_eq!(frame.friction_dissipation_j,0.0);close(frame.slip_distance_m,0.0,1e-12);
    let state=states(&s);close(state[0].displacement_m_sqrt_kg/MASS[0].sqrt(),state[1].displacement_m_sqrt_kg/MASS[1].sqrt(),1e-14);
    close(momentum(&s),f*DT,1e-16);assert!(frame.energy_residual_j.abs()<=frame.energy_tolerance_j);
}
#[test]
fn zero_repulsion_sliding_preserves_momentum_and_closes_native_energy() {
    let mut s=mechanics(1.0,1e-6);let initial=s.total_energy_j().unwrap();let p=momentum(&s);let mut heat=0.0;
    for _ in 0..100 {
        let frame=*s.step(&[0.0,0.0],drive()).unwrap();assert_eq!(frame.regime,ModalFrictionRegime::Sliding);
        close(frame.traction_n,-0.5,1e-12);close(frame.friction_dissipation_j,frame.traction_n*frame.slip_distance_m,1e-18);
        assert!(frame.friction_dissipation_j>0.0);assert!(frame.energy_residual_j.abs()<=frame.energy_tolerance_j);
        heat+=frame.friction_dissipation_j;close(momentum(&s),p,1e-14);
    }
    assert_eq!(s.samples_rendered(),100);close(initial-s.total_energy_j().unwrap(),heat,1e-13);
}
#[test]
fn detached_patch_is_released_and_failed_trials_leave_native_state_unchanged() {
    let mut detached=mechanics(1.0,0.0);let energy=detached.total_energy_j().unwrap();
    let frame=*detached.step(&[0.0,0.0],drive()).unwrap();assert_eq!(frame.regime,ModalFrictionRegime::Released);
    assert_eq!(frame.traction_n,0.0);assert_eq!(frame.friction_dissipation_j,0.0);close(detached.total_energy_j().unwrap(),energy,1e-14);
    let mut s=mechanics(1.0,1e-6);s.step(&[0.0,0.0],drive()).unwrap();
    let before=states(&s);let last=s.last_frame().copied();let clock=s.samples_rendered();
    assert!(s.step(&[0.0,0.0],FrictionDrive {normal_force_n:1.0,..drive()}).is_err());
    assert_eq!(states(&s),before);assert_eq!(s.samples_rendered(),clock);assert_eq!(s.last_frame().copied(),last);
    let gate=CancelGate::new();gate.request();assert!(s.step_under_gate(&[0.0,0.0],drive(),&gate).is_err());
    assert_eq!(states(&s),before);assert_eq!(s.last_frame().copied(),last);
}

#[test]
fn adhesive_weakening_laws_continue_sliding_inside_static_cone_at_zero_repulsion() {
    for pressure in [
        FrictionLaw::Stribeck {static_mu:1.2,kinetic_mu:0.35,characteristic_speed:0.05,viscous_per_speed:0.0},
        FrictionLaw::Hyperbolic {static_mu:1.2,kinetic_mu:0.35,characteristic_speed:0.05},
    ] {
        let mut system=mechanics_with_pressure(0.018,1e-6,0.0,pressure);
        let initial=system.total_energy_j().unwrap();let initial_momentum=momentum(&system);let mut loss=0.0;
        for step in 0..3 {
            let frame=*system.step(&[0.0,0.0],drive()).unwrap();
            assert_eq!(frame.regime,ModalFrictionRegime::Sliding);
            close(frame.traction_n,-0.5,1e-12);
            assert!(frame.friction_dissipation_j>0.0);
            assert!(frame.energy_residual_j.abs()<=frame.energy_tolerance_j);
            close(momentum(&system),initial_momentum,1e-16);
            loss+=frame.friction_dissipation_j;
            if step==2 {
                let roots=frame.roots.unwrap();
                assert_eq!((roots.found,roots.selected),(2,1),"a moving junction must not be captured by its static cone");
            }
        }
        close(initial-system.total_energy_j().unwrap(),loss,1e-17);
        let before=states(&system);let last=system.last_frame().copied();
        assert!(system.step(&[0.0,0.0],FrictionDrive {normal_force_n:1.0,..drive()}).is_err());
        assert_eq!(states(&system),before);assert_eq!(system.last_frame().copied(),last);
    }
}

#[test]
fn nonzero_repulsive_load_uses_actual_adhesive_law_in_native_work_balance() {
    for pressure in [
        FrictionLaw::Stribeck {static_mu:1.2,kinetic_mu:0.35,characteristic_speed:0.05,viscous_per_speed:0.1},
        FrictionLaw::Hyperbolic {static_mu:1.2,kinetic_mu:0.35,characteristic_speed:0.05},
    ] {
        let normal=0.3;
        let mut system=mechanics_with_pressure(1.0,1e-6,normal,pressure);
        let initial=system.total_energy_j().unwrap();let initial_momentum=momentum(&system);let mut loss=0.0;
        for _ in 0..40 {
            let frame=*system.step(&[0.0,0.0],FrictionDrive {normal_force_n:normal,..drive()}).unwrap();
            assert_eq!(frame.regime,ModalFrictionRegime::Sliding);
            let binding=system.adhesive_binding().unwrap().tangential();
            let expected=binding.law().sliding_capacity_n(binding.contact(),frame.slip_distance_m.abs()/DT).unwrap();
            close(frame.traction_n,-expected,1e-11);
            close(frame.friction_dissipation_j,frame.traction_n*frame.slip_distance_m,1e-18);
            assert!(frame.energy_residual_j.abs()<=frame.energy_tolerance_j);
            close(momentum(&system),initial_momentum,1e-14);
            loss+=frame.friction_dissipation_j;
        }
        close(initial-system.total_energy_j().unwrap(),loss,1e-13);
    }
}
