//! Static contact balance and continuation through the authentic dynamic owner.
use fs_couple::modal_acoustic_time::{ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel};
use fs_couple::render::schedule::force::coupled::{CoupledModalSystem, ModalAttachment, ModalConnection, ModalCouplingConfig, ModalCouplingError};
use fs_couple::render::schedule::force::coupled::contact::{ContactModalSystem, ModalContact, ModalContactConfig};
use fs_couple::render::schedule::force::coupled::contact::multiple::{MultiContactConfig, MultiContactModalSystem};
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;
use fs_math::c64::C64;

fn config() -> ModalCouplingConfig {
    ModalCouplingConfig { max_modes:16,max_connections:4,max_setup_terms:4096,nyquist_guard_fraction:0.9,
        maximum_total_energy_j:1000.0,maximum_abs_pressure_pa:1e6,maximum_abs_connection_force_n:1e6,
        solve_relative_tolerance:1e-11,energy_absolute_tolerance_j:1e-10,energy_relative_tolerance:1e-8 }
}
fn set() -> MultiContactConfig { MultiContactConfig { max_contacts:8,max_sweeps:128,max_setup_terms:4096 } }
fn limits() -> ModalContactConfig { ModalContactConfig { max_iterations:96,maximum_force_n:1e6,
    maximum_penetration_m:0.1,force_absolute_tolerance_n:1e-10,force_relative_tolerance:1e-11 } }
fn image(rate:u32, omega:f64, budget:ModalAcousticTimeBudget) -> ModalAcousticTimeModel {
    ModalAcousticTimeModel::try_new(rate,vec![ModalAcousticMode {angular_frequency_rad_s:omega,
        damping_ratio:0.01,pressure_per_modal_velocity:C64::new(1.0,0.0)}],budget).unwrap()
}
fn network(rate:u32, omegas:&[f64]) -> CoupledModalSystem {
    CoupledModalSystem::new(omegas.iter().map(|w|image(rate,*w,ModalAcousticTimeBudget::audible_reference())).collect(),
        vec![],config(),&CancelGate::new()).unwrap()
}
fn contact(left:usize,right:usize,k:f64,alpha:f64,gap:f64,chi:f64) -> (ModalContact,ModalContactConfig) {
    (ModalContact {left:ModalAttachment {component:left,shapes:vec![1.0]},right:ModalAttachment {component:right,shapes:vec![1.0]},
        law:Obstacle::new(vec![-1.0],1,1,vec![gap],vec![1.0],k,alpha,"authored preload test".into()).unwrap()
            .with_internal_loss(chi).unwrap()},limits())
}
fn states(s:&CoupledModalSystem)->Vec<ModalAcousticState> {s.components().iter().flat_map(|m|m.states()).copied().collect()}

#[test]
fn quadratic_contact_preload_matches_closed_form_not_independent_modal_equilibria() {
    let mut s=network(48000,&[800.0,800.0]);let contacts=vec![contact(0,1,1e8,2.0,0.0001,0.3)];
    let energy=s.initialize_contact_equilibrium(&[400.0,0.0],&contacts,set(),&CancelGate::new()).unwrap();
    let compliance=2.0/800.0_f64.powi(2);let free=400.0/800.0_f64.powi(2)-0.0001;
    let penetration=2.0*free/(1.0+(1.0+4.0*compliance*1e8*free).sqrt());
    let reaction=1e8*penetration*penetration;
    let q=states(&s);assert!((q[0].displacement_m_sqrt_kg-(400.0-reaction)/640000.0).abs()<1e-13);
    assert!((q[1].displacement_m_sqrt_kg-reaction/640000.0).abs()<1e-13);
    assert!(q.iter().all(|x|x.velocity_m_sqrt_kg_per_s==0.0));assert_eq!(s.samples_rendered(),0);
    assert!(s.last_frame().is_none());
    let expected=0.5*640000.0*(q[0].displacement_m_sqrt_kg.powi(2)+q[1].displacement_m_sqrt_kg.powi(2))
        +1e8*penetration.powi(3)/3.0;
    assert!((energy-expected).abs()<1e-10);
    let (point,c)=contacts[0].clone();let mut dynamic=ContactModalSystem::new(s,point,c,&CancelGate::new()).unwrap();
    for _ in 0..32 { assert!(dynamic.step(&[400.0,0.0]).unwrap().observer_pressure_pa.abs()<1e-7); }
    let mut peak=0.0_f64;
    for _ in 0..128 {peak=peak.max(dynamic.step(&[0.0,0.0]).unwrap().observer_pressure_pa.abs());}
    assert!(peak>0.01,"load release must retain and excite stored displacement");
}

#[test]
fn shared_contact_equilibrium_balances_all_bodies_and_is_sample_rate_and_loss_independent() {
    let loads=[400.0,0.0,-300.0];let omega=[800.0,1000.0,1200.0];let mut reference=None;
    for rate in [24000,48000,96000] {
        for chi in [0.0,20.0] {
            let mut s=network(rate,&omega);let contacts=vec![contact(0,1,1e8,2.0,0.0001,chi),contact(1,2,8e7,2.0,0.00015,chi)];
            s.initialize_contact_equilibrium(&loads,&contacts,set(),&CancelGate::new()).unwrap();let q=states(&s);
            let r0=1e8*(q[0].displacement_m_sqrt_kg-q[1].displacement_m_sqrt_kg-0.0001).max(0.0).powi(2);
            let r1=8e7*(q[1].displacement_m_sqrt_kg-q[2].displacement_m_sqrt_kg-0.00015).max(0.0).powi(2);
            assert!(r0>0.0 && r1>0.0);
            for (i,force) in [loads[0]-r0,r0-r1,loads[2]+r1].iter().enumerate() {
                assert!((omega[i]*omega[i]*q[i].displacement_m_sqrt_kg-force).abs()<1e-7);
            }
            if let Some(old)=&reference {assert_eq!(&q,old);} else {reference=Some(q);}
        }
    }
}

#[test]
fn bilateral_rest_offsets_and_inactive_contacts_do_not_change_the_static_owner() {
    let make=||CoupledModalSystem::new(vec![image(48000,800.0,ModalAcousticTimeBudget::audible_reference()),
        image(48000,1100.0,ModalAcousticTimeBudget::audible_reference())],vec![ModalConnection {
            left:ModalAttachment {component:0,shapes:vec![1.0]},right:ModalAttachment {component:1,shapes:vec![1.0]},
            stiffness_n_m:3e5,damping_n_s_m:200.0,rest_extension_m:0.0003}],config(),&CancelGate::new()).unwrap();
    let mut a=make();let mut b=make();let gate=CancelGate::new();
    let energy=a.initialize_static_equilibrium(&[100.0,-50.0],&gate).unwrap();
    let actual=b.initialize_contact_equilibrium(&[100.0,-50.0],&[contact(0,1,1e8,2.0,1.0,0.1)],set(),&gate).unwrap();
    assert_eq!(states(&a),states(&b));assert_eq!(energy.to_bits(),actual.to_bits());
}

#[test]
fn final_contact_restraint_can_admit_a_load_whose_free_prediction_exceeds_state_limits() {
    let budget=ModalAcousticTimeBudget {maximum_abs_displacement_m_sqrt_kg:0.0001,..ModalAcousticTimeBudget::audible_reference()};
    let mut s=CoupledModalSystem::new(vec![image(48000,100.0,budget)],vec![],config(),&CancelGate::new()).unwrap();
    let (mut point,c)=contact(0,0,1e8,1.0,0.0,0.0);point.right.shapes[0]=0.0;
    assert!(s.initialize_static_equilibrium(&[1000.0],&CancelGate::new()).is_err());
    let energy=s.initialize_contact_equilibrium(&[1000.0],&[(point,c)],set(),&CancelGate::new()).unwrap();
    assert!(energy>0.0);assert!((states(&s)[0].displacement_m_sqrt_kg-1000.0/(1e4+1e8)).abs()<1e-13);
}

#[test]
fn initially_overlapping_geometry_is_admitted_against_the_solved_not_placeholder_penetration() {
    let mut s=network(48000,&[100.0]);let (mut point,mut c)=contact(0,0,1e8,1.0,-0.01,0.0);
    point.right.shapes[0]=0.0;c.maximum_penetration_m=1e-4;
    s.initialize_contact_equilibrium(&[0.0],&[(point.clone(),c)],set(),&CancelGate::new()).unwrap();
    let q=states(&s)[0].displacement_m_sqrt_kg;assert!((q+1e6/(1e4+1e8)).abs()<1e-12);
    assert!(q+0.01<1e-4);
    assert!(ContactModalSystem::new(s,point,c,&CancelGate::new()).is_ok());
}

#[test]
fn force_penetration_energy_and_work_refusals_publish_nothing_and_allow_retry() {
    let contacts=vec![contact(0,1,1e8,2.0,0.0001,0.0)];let gate=CancelGate::new();
    let mut s=network(48000,&[800.0,800.0]);let before=states(&s);
    for which in 0..4 {
        let mut points=contacts.clone();let mut c=set();
        match which {0=>points[0].1.maximum_force_n=1.0,1=>points[0].1.maximum_penetration_m=1e-8,
            2=>points[0].1.max_iterations=1,_=>c.max_setup_terms=1}
        assert!(s.initialize_contact_equilibrium(&[400.0,0.0],&points,c,&gate).is_err());
        assert_eq!(states(&s),before);assert_eq!(s.samples_rendered(),0);assert!(s.last_frame().is_none());
    }
    let mut c=config();c.maximum_total_energy_j=0.114;
    // Network-only storage is below 0.114 J; contact storage must still count.
    let mut bounded=CoupledModalSystem::new(vec![image(48000,800.0,ModalAcousticTimeBudget::audible_reference()),
        image(48000,800.0,ModalAcousticTimeBudget::audible_reference())],vec![],c,&gate).unwrap();
    assert!(bounded.initialize_contact_equilibrium(&[400.0,0.0],&contacts,set(),&gate).is_err());
    assert!(states(&bounded).iter().all(|x|x.displacement_m_sqrt_kg==0.0));
    let cancelled=CancelGate::new();cancelled.request();
    assert!(matches!(s.initialize_contact_equilibrium(&[400.0,0.0],&contacts,set(),&cancelled),Err(ModalCouplingError::Cancelled)));
    assert_eq!(states(&s),before);
    s.initialize_contact_equilibrium(&[400.0,0.0],&contacts,set(),&gate).unwrap();
    let accepted=states(&s);assert!(s.initialize_contact_equilibrium(&[400.0,0.0],&contacts,set(),&gate).is_err());
    assert_eq!(states(&s),accepted);
}

#[test]
fn multiple_preloaded_contacts_continue_without_a_startup_impact() {
    let mut s=network(48000,&[800.0,1000.0,1200.0]);let contacts=vec![contact(0,1,1e8,2.0,0.0001,0.1),contact(1,2,8e7,2.0,0.00015,0.2)];
    let loads=[400.0,0.0,-300.0];let gate=CancelGate::new();
    let energy=s.initialize_contact_equilibrium(&loads,&contacts,set(),&gate).unwrap();
    let mut dynamic=MultiContactModalSystem::new(s,contacts,set(),&gate).unwrap();
    assert!((dynamic.total_energy_j().unwrap()-energy).abs()<1e-12);
    for _ in 0..37 {assert!(dynamic.step(&loads).unwrap().observer_pressure_pa.abs()<1e-7);}
    let mut peak=0.0_f64;for _ in 0..128 {peak=peak.max(dynamic.step(&[0.0;3]).unwrap().observer_pressure_pa.abs());}
    assert!(peak>0.01);
}
