//! Synthetic interface domains test physical wiring, not measured rosin data.
use fs_couple::modal_acoustic_time::{ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel};
use fs_couple::render::schedule::force::coupled::{CoupledModalSystem, ModalAttachment, ModalCouplingConfig};
use fs_couple::render::schedule::force::coupled::contact::friction::{FrictionDrive, FrictionModalSystem, ModalFrictionConfig, ModalFrictionError, ModalFrictionPort, ModalFrictionRegime};
use fs_couple::render::schedule::force::coupled::contact::friction::interface::DryFrictionBinding;
use fs_exec::CancelGate;
use fs_math::c64::C64;
use fs_tribo::{ApplicabilityRange as Range, DryFrictionApplicability, DryInterfaceSystemCard, FrictionLaw, InputAuthority, InterfaceMedium, InterfaceSystemRef, TriboError};

fn law() -> FrictionLaw { FrictionLaw::Stribeck { static_mu: 0.8, kinetic_mu: 0.4,
    characteristic_speed: 0.04, viscous_per_speed: 0.03 } }
fn binding(law: FrictionLaw, pressure: (f64,f64), speed: (f64,f64)) -> DryFrictionBinding {
    DryFrictionBinding::new(DryInterfaceSystemCard::new(
        InterfaceSystemRef::new("synthetic/string->bow", "synthetic/dry-history", "synthetic/test", InputAuthority::SyntheticFixture, InterfaceMedium::Dry).unwrap(),
        law, DryFrictionApplicability::new(Range::new(280.0,320.0).unwrap(),
            Range::new(pressure.0,pressure.1).unwrap(),Range::new(speed.0,speed.1).unwrap()).unwrap()).unwrap(),
        293.15,0.001).unwrap()
}
fn limits() -> ModalFrictionConfig { ModalFrictionConfig { max_iterations:96,
    maximum_force_n:100.0,maximum_slip_speed_m_s:10.0,force_absolute_tolerance_n:1e-11,
    force_relative_tolerance:1e-10,sticking_tolerance_m:1e-12 } }
fn network(omega: f64) -> CoupledModalSystem {
    let model = ModalAcousticTimeModel::try_new(48_000,vec![ModalAcousticMode {
        angular_frequency_rad_s:omega,damping_ratio:0.0,
        pressure_per_modal_velocity:C64::new(1.0,0.0), // authored synthetic readout
    }],ModalAcousticTimeBudget::audible_reference()).unwrap();
    CoupledModalSystem::new(vec![model],vec![],ModalCouplingConfig {
        max_modes:8,max_connections:4,max_setup_terms:4096,nyquist_guard_fraction:0.9,
        maximum_total_energy_j:100.0,maximum_abs_pressure_pa:1e6,maximum_abs_connection_force_n:1e6,
        solve_relative_tolerance:1e-11,energy_absolute_tolerance_j:1e-12,energy_relative_tolerance:1e-8,
    },&CancelGate::new_clock_free()).unwrap()
}
fn port() -> ModalAttachment { ModalAttachment {component:0,shapes:vec![1.0]} }
fn raw(omega: f64, law: FrictionLaw) -> FrictionModalSystem {
    FrictionModalSystem::new(network(omega),ModalFrictionPort {left:port(),right:None,law},limits(),&CancelGate::new_clock_free()).unwrap()
}
fn bound(omega: f64, binding: DryFrictionBinding) -> FrictionModalSystem {
    FrictionModalSystem::new_with_interface(network(omega),port(),None,binding,limits(),&CancelGate::new_clock_free()).unwrap()
}
fn drive(speed_m_s:f64,normal_force_n:f64) -> FrictionDrive {FrictionDrive {speed_m_s,normal_force_n}}
fn state(s:&FrictionModalSystem)->Vec<ModalAcousticState> {s.components()[0].states().to_vec()}

#[test]
fn complete_interface_law_preserves_original_motion_work_and_velocity_strengthening() {
    for law in [law(),FrictionLaw::Coulomb {static_mu:0.8,kinetic_mu:0.4}] {
        let source=binding(law.clone(),(1.0,1e6),(0.0,10.0));
        let mut actual=bound(800.0,source.clone()); let mut expected=raw(800.0,law);
        assert_eq!(actual.interface(),Some(&source));
        for i in 0..256 {
            let d=drive(if i<128 {0.3} else {-0.3},if i<192 {1.0} else {0.0});
            let a=*actual.step(&[0.02],d).unwrap(); let b=*expected.step(&[0.02],d).unwrap();
            assert_eq!(a,b); assert_eq!(state(&actual),state(&expected));
            assert!(a.energy_residual_j.abs()<=a.energy_tolerance_j);
        }
        assert!(actual.total_energy_j().unwrap()>0.0);
    }
}

#[test]
fn actual_average_slip_cannot_bypass_a_positive_source_speed_floor() {
    let narrow=binding(law(),(1.0,1e6),(0.001,1.0));
    let mut actual=bound(0.0,narrow);
    let mut expected=raw(0.0,law());
    // An external load supplies the large inertial force. Friction can stick:
    // v0=0, prescribed u=.01, v1=.02, hence both endpoint slips are .01,
    // while the work-conjugate average slip is ZERO and outside this card.
    let d=drive(0.01,1.0);
    let f=*expected.step(&[960.0],d).unwrap();
    assert_eq!(f.regime,ModalFrictionRegime::Sticking);
    assert!(matches!(actual.step(&[960.0],d),Err(ModalFrictionError::Law(TriboError::OutsideApplicability {field:"slip_speed_mps",..}))));
    assert_eq!(actual.samples_rendered(),0);
    assert_eq!(actual.components()[0].states()[0],ModalAcousticState::default());
}

#[test]
fn endpoint_domain_refusal_is_atomic_even_when_initial_and_average_rates_pass() {
    let mut actual=bound(0.0,binding(law(),(0.0,1e6),(0.0,0.02)));
    let mut expected=raw(0.0,law());
    let d=drive(0.0,0.001);
    let frame=*expected.step(&[1440.0],d).unwrap();
    assert!((frame.slip_distance_m/expected.sample_period_s()).abs()<0.02);
    assert!(expected.components()[0].states()[0].velocity_m_sqrt_kg_per_s>0.02);
    let before=state(&actual);
    assert!(actual.step(&[1440.0],d).is_err());
    assert_eq!(actual.samples_rendered(),0); assert_eq!(state(&actual),before);
    assert!(actual.last_frame().is_none());
    assert_eq!(actual.step(&[0.0],drive(0.0,0.0)).unwrap().network_energy_j,0.0);
}

#[test]
fn nominal_pressure_and_release_use_actual_load_without_silencing_vibration() {
    let source=binding(law(),(900.0,1100.0),(0.0,1.0));
    let mut actual=bound(800.0,source.clone());
    actual.step(&[0.0],drive(0.3,1.0)).unwrap();
    let before=state(&actual); let frame=*actual.last_frame().unwrap();
    assert!(matches!(actual.step(&[0.0],drive(0.3,2.0)),Err(ModalFrictionError::Law(TriboError::OutsideApplicability {field:"contact_pressure_pa",..}))));
    assert_eq!(state(&actual),before); assert_eq!(actual.last_frame(),Some(&frame));
    // Zero load bypasses pressure and speed admission: a lifted bow is not contact.
    let release=*actual.step(&[0.0],drive(100.0,0.0)).unwrap();
    assert_eq!(release.regime,ModalFrictionRegime::Released);
    assert_eq!(release.traction_n,0.0); assert!(release.network_energy_j>0.0);
    assert_eq!(actual.interface(),Some(&source));
}

#[test]
fn source_temperature_geometry_mismatched_law_and_rebinding_refuse() {
    let source=binding(law(),(1.0,1e6),(0.0,10.0));
    for area in [0.0,-1.0,f64::NAN,f64::INFINITY] {
        assert!(DryFrictionBinding::new(source.card().clone(),293.15,area).is_err());
    }
    assert!(DryFrictionBinding::new(source.card().clone(),321.0,0.001).is_err());
    let zero_pressure = binding(law(),(0.0,1e6),(0.0,10.0));
    let huge = DryFrictionBinding::new(zero_pressure.card().clone(),293.15,f64::MAX).unwrap();
    let mut unrepresentable = bound(800.0,huge);
    assert!(unrepresentable.step(&[0.0],drive(0.0,f64::MIN_POSITIVE)).is_err());
    assert_eq!(unrepresentable.samples_rendered(),0);
    assert!(raw(800.0,FrictionLaw::Coulomb {static_mu:0.8,kinetic_mu:0.4}).with_interface(source.clone()).is_err());
    assert!(bound(800.0,source.clone()).with_interface(source.clone()).is_err());
    let mut original=raw(800.0,law()); original.step(&[0.0],drive(0.3,1.0)).unwrap();
    let previous=state(&original); let frame=*original.last_frame().unwrap();
    let mut attached=original.with_interface(source).unwrap();
    assert_eq!(state(&attached),previous); assert_eq!(attached.last_frame(),Some(&frame));
    let stop=CancelGate::new_clock_free();stop.request();
    assert!(attached.step_under_gate(&[0.0],drive(0.3,1.0),&stop).is_err());
    assert_eq!(state(&attached),previous); assert_eq!(attached.last_frame(),Some(&frame));
}

use fs_couple::render::schedule::force::coupled::contact::friction::multiple::{MultiFrictionConfig,MultiFrictionModalSystem};
fn multi(source: Option<DryFrictionBinding>) -> MultiFrictionModalSystem {
    let original=raw(800.0,law());
    let original=if let Some(source)=source {original.with_interface(source).unwrap()} else {original};
    original.with_additional_ports(vec![(ModalFrictionPort {
        left:port(),right:None,law:law(),
    },limits())],MultiFrictionConfig {max_ports:4,max_sweeps:128,max_setup_terms:4096},
        &CancelGate::new_clock_free()).unwrap()
}

#[test]
fn every_joint_contact_retains_its_own_domain_and_a_late_refusal_publishes_nothing() {
    let wide=binding(law(),(1.0,1e6),(0.0,10.0));
    let tight=binding(law(),(900.0,1100.0),(0.0,0.02));
    let mut actual=multi(Some(wide.clone())).with_interface(1,tight.clone()).unwrap();
    assert_eq!(actual.interface(0),Some(&wide));assert_eq!(actual.interface(1),Some(&tight));
    let release=[drive(0.0,0.0);2];
    actual.step(&[0.0],&release).unwrap();
    let before=actual.components()[0].states().to_vec();let frame=actual.last_frame().unwrap().clone();
    // Port 1 initially passes at zero slip. External force drives its endpoint
    // outside the domain after the JOINT solve; port 0 must not publish either.
    assert!(actual.step(&[1440.0],&[drive(0.3,1.0),drive(0.0,1.0)]).is_err());
    assert_eq!(actual.components()[0].states(),before);assert_eq!(actual.last_frame(),Some(&frame));
    assert_eq!(actual.samples_rendered(),1);
    actual.step(&[0.0],&release).unwrap();assert_eq!(actual.samples_rendered(),2);
    // Independent N/area check for the added port, not only the original core.
    assert!(actual.step(&[0.0],&[drive(0.0,0.0),drive(0.0,2.0)]).is_err());
    assert_eq!(actual.samples_rendered(),2);
}

#[test]
fn interface_checks_preserve_joint_and_singleton_arithmetic_exactly() {
    let source=binding(law(),(1.0,1e6),(0.0,10.0));
    let mut a=multi(Some(source.clone())).with_interface(1,source.clone()).unwrap();
    let mut b=multi(None);
    for i in 0..128 {
        let ds=[drive(0.3,1.0),drive(-0.2,if i<64 {0.5} else {0.0})];
        assert_eq!(a.step(&[0.0],&ds).unwrap(),b.step(&[0.0],&ds).unwrap());
        assert_eq!(a.components()[0].states(),b.components()[0].states());
    }
    let joint=MultiFrictionConfig {max_ports:1,max_sweeps:1,max_setup_terms:4096};
    let mut scalar=bound(800.0,source.clone());
    let mut one=raw(800.0,law()).with_additional_ports(vec![],joint,&CancelGate::new_clock_free()).unwrap()
        .with_interface(0,source).unwrap();
    for _ in 0..32 {
        let s=*scalar.step(&[0.0],drive(0.3,1.0)).unwrap();
        let f=one.step(&[0.0],&[drive(0.3,1.0)]).unwrap();
        assert_eq!(f.observer_pressure_pa,s.observer_pressure_pa);
        assert_eq!(f.ports[0].traction_n,s.traction_n);
        assert_eq!(one.components()[0].states(),scalar.components()[0].states());
    }
}
