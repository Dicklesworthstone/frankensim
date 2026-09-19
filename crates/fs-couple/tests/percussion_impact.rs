use fs_couple::render::plate::impact::{ImpactBody,BodyPotential,ImpactSystem,ImpactConfig,VolumeSpring};
use fs_couple::render::plate::impact::felt::{FeltPad,KelvinBranch};
use fs_couple::render::plate::impact::striker::{RadiusStation,StrikerProperties};
use fs_couple::modal_acoustic_time::ModalAcousticState;
use fs_material::fiber::WoolFelt;
use fs_exec::CancelGate;
use fs_dcontact::Obstacle;
fn config()->ImpactConfig {ImpactConfig{dt_s:2e-6,max_steps:2000,maximum_energy_j:1.0,
    energy_absolute_tolerance_j:1e-10,energy_relative_tolerance:1e-7,maximum_generalized_force:1e4}}
fn fixture(felt:bool,contact:bool)->ImpactSystem {
    let (stick,weight)=ImpactBody::free_mass(0.02,-0.0002,0.8).unwrap();
    let body=ImpactBody{potential:BodyPotential::Linear(vec![2.0*std::f64::consts::PI*500.0]),
        initial:vec![ModalAcousticState::default()],damping_per_s:vec![1.0]};
    let law=Obstacle::new(vec![weight,-2.0],1,2,vec![0.0],vec![1.0],if contact {2e7}else{0.0},1.5,
        "synthetic impact regression".into()).unwrap();
    let pads=if felt {vec![FeltPad{area_m2:0.001,thickness_m:0.006,precompression_m:0.0006,
        weights:vec![0.0,2.0],law:WoolFelt::new(30000.0,0.2,2.2,3.0,0.15,0.7).unwrap(),
        prior_maximum_strain:0.1,creep:vec![KelvinBranch{stiffness_n_m:1500.0,viscosity_n_s_m:8.0}]}]}else{vec![]};
    ImpactSystem::new(vec![stick,body],vec![law],pads,vec![],config()).unwrap()
}
#[test]
fn launched_mass_strikes_excites_and_rebounds_without_an_authored_force() {
    let mut a=fixture(true,true);let mut uncoupled=fixture(false,false);let gate=CancelGate::new_clock_free();
    let initial=a.stored_energy_j();let mut loss=0.0;let mut receiver=0.0_f64;let mut returned=false;
    for _ in 0..1200 {
        let frame=a.step(&[0.0;2],&gate).unwrap();loss+=frame.dissipated_energy_j;
        assert!(frame.balance_residual_j.abs()<1e-8);
        receiver=receiver.max(a.state()[2].abs());returned|=a.state()[1]<0.0;
        uncoupled.step(&[0.0;2],&gate).unwrap();assert_eq!(uncoupled.state()[2],0.0);
    }
    assert!(receiver>1e-6 && returned);
    assert!(loss>1e-6);assert!((a.stored_energy_j()-initial+loss).abs()<1e-7);
    assert!(a.state()[4].abs()>0.0,"felt recovery must retain a physical state");
}
#[test]
fn cancellation_bad_inputs_and_shorter_partitions_preserve_full_history() {
    let mut a=fixture(true,true);let mut b=fixture(true,true);let gate=CancelGate::new_clock_free();
    for _ in 0..140 {a.step(&[0.0;2],&gate).unwrap();b.step(&[0.0;2],&gate).unwrap();}
    let state=a.state().to_vec();let h=a.felt_history(0);let cancel=CancelGate::new_clock_free();cancel.request();
    assert!(a.step(&[0.0;2],&cancel).is_err());
    assert!(a.step(&[f64::NAN,0.0],&gate).is_err());assert!(a.step(&[0.0],&gate).is_err());
    assert_eq!(state,a.state());assert_eq!(h,a.felt_history(0));assert_eq!(a.samples(),140);
    for count in [3,9,1,27] {for _ in 0..count {a.step(&[0.0;2],&gate).unwrap();b.step(&[0.0;2],&gate).unwrap();}}
    assert_eq!(a.state().iter().map(|v|v.to_bits()).collect::<Vec<_>>(),b.state().iter().map(|v|v.to_bits()).collect::<Vec<_>>());
    assert_eq!(a.felt_history(0),b.felt_history(0));
}
#[test]
fn sealed_volume_couples_two_heads_without_direct_force_on_the_receiver() {
    let body=|q|ImpactBody{potential:BodyPotential::Linear(vec![800.0]),initial:vec![ModalAcousticState{
        displacement_m_sqrt_kg:q,velocity_m_sqrt_kg_per_s:0.0}],damping_per_s:vec![0.0]};
    let mut a=ImpactSystem::new(vec![body(1e-5),body(0.0)],vec![],vec![],
        vec![VolumeSpring{bulk_modulus_pa:1.4e5,volume_m3:0.01,areas:vec![0.1,-0.1]}],config()).unwrap();
    let initial=a.stored_energy_j();let gate=CancelGate::new_clock_free();
    for _ in 0..100 {a.step(&[0.0;2],&gate).unwrap();}
    assert!(a.state()[2].abs()>1e-9);
    assert!((a.stored_energy_j()-initial).abs()<1e-10);
}
#[test]
fn profile_inertia_matches_a_uniform_cylinder_and_changes_with_grip() {
    let (length,radius,rho)=(0.4,0.007,800.0);
    let profile=[RadiusStation{position_m:0.0,radius_m:radius},RadiusStation{position_m:length,radius_m:radius}];
    let s=StrikerProperties::from_profile(&profile,rho,0.0,length).unwrap();
    let m=rho*std::f64::consts::PI*radius*radius*length;
    let inertia=m*(length*length/3.0+radius*radius/4.0);
    assert!((s.mass_kg-m).abs()<1e-15);assert!((s.pivot_inertia_kg_m2-inertia).abs()<1e-15);
    assert!((s.center_of_mass_m-length/2.0).abs()<1e-15);
    let grip=StrikerProperties::from_profile(&profile,rho,0.12,length).unwrap();
    assert_ne!(s.contact_effective_mass_kg.to_bits(),grip.contact_effective_mass_kg.to_bits());
    assert!(StrikerProperties::from_profile(&profile,rho,length,length).is_err());
}
