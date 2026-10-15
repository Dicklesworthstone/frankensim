//! Source-only substitutions in the real joint contact/momentum/energy solve.
#[path="../../fs-contact/tests/source_point/common.rs"]
mod common;
use common::*;
use fs_contact::normal_patch::source::{PointContactGeometry, SourcePointContact};
use fs_couple::render::plate::impact::{ImpactBody, ImpactConfig, ImpactSystem};
use fs_couple::render::plate::impact::material::MaterialContact;
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;
use fs_ga::Vec3;
use fs_matdb::MaterialCard;
use fs_material::state_point::MaterialPropertySelection;
use fs_mbd_material::{RegionBinding, VolumeMoments, resolve_body};

fn config() -> ImpactConfig { ImpactConfig { dt_s:2e-6,max_steps:2500,maximum_energy_j:1.0,
    energy_absolute_tolerance_j:1e-12,energy_relative_tolerance:1e-7,maximum_generalized_force:1.0 } }
fn mass(card: &MaterialCard) -> f64 {
    resolve_body(&[RegionBinding { key:"specimen",material:card,point:&query(),
        geometry:VolumeMoments::cuboid(Vec3::new(0.01,0.01,0.01)).unwrap(),
        selection:MaterialPropertySelection::SingleClaimOnly }],1,||false).unwrap().inertia().mass()
}
fn pair(rho: f64, young: f64, chi: Option<f64>) -> (Vec<ImpactBody>,MaterialContact,[f64;2]) {
    let a=solid("A",rho,young);let b=solid("B",2000.0,2e7);
    let ma=mass(&a);let mb=mass(&b);
    let (a_body,wa)=ImpactBody::free_mass(ma,-2e-5,0.04).unwrap();
    let (b_body,wb)=ImpactBody::free_mass(mb,0.0,0.0).unwrap();
    let source=SourcePointContact::new(bind(&a,&b,chi,0.0),
        PointContactGeometry::SpherePlane {effective_radius_m:0.1},envelope()).unwrap();
    let contact=MaterialContact::new("interface".into(),source,vec![wa,-wb],0.0).unwrap();
    (vec![a_body,b_body],contact,[ma,mb])
}
fn system(rho:f64,young:f64,chi:Option<f64>) -> (ImpactSystem,[f64;2]) {
    let (bodies,contact,masses)=pair(rho,young,chi);
    (ImpactSystem::new_with_material_contacts(bodies,vec![contact],vec![],vec![],config()).unwrap(),masses)
}
fn states_close(a:&[f64],b:&[f64]) {
    assert_eq!(a.len(),b.len());
    for (x,y) in a.iter().zip(b) {assert!((x-y).abs()<1e-8,"{x} != {y}");}
}

#[test]
fn sourced_contact_uses_the_existing_potential_and_dissipative_runtime() {
    let (bodies,contact,_)=pair(1000.0,1e7,Some(0.4));
    let raw=Obstacle::new(contact.weights().to_vec(),1,2,vec![contact.gap_m()],vec![1.0],
        contact.source().stiffness_n_per_m_3_2(),1.5,"explicit comparison only".into()).unwrap()
        .with_internal_loss(contact.source().dissipation_s_per_m()).unwrap();
    let mut original=ImpactSystem::new(bodies.clone(),vec![raw],vec![],vec![],config()).unwrap();
    let retained=contact.source().model().identity();
    let mut sourced=ImpactSystem::new_with_material_contacts(bodies,vec![contact],vec![],vec![],config()).unwrap();
    let gate=CancelGate::new_clock_free();let mut peak=0.0_f64;
    for _ in 0..1000 {
        let a=original.step(&[0.0;2],&gate).unwrap();let b=sourced.step(&[0.0;2],&gate).unwrap();
        assert_eq!(original.state(),sourced.state());
        assert_eq!(a.stored_energy_j,b.stored_energy_j);
        assert_eq!(a.dissipated_energy_j,b.dissipated_energy_j);
        let r=sourced.material_contact_observation(0).unwrap();peak=peak.max(r.normal_force_n);
        assert!(r.normal_force_n>=0.0 && r.dissipated_power_w>=0.0);
    }
    assert!(peak>0.05);
    assert_eq!(sourced.material_contacts()[0].source().model().identity(),retained);
}

#[test]
fn reference_prepared_and_analytic_keep_one_source_and_one_energy_balance() {
    let (mut reference,m)=system(1000.0,1e7,Some(0.4));
    let mut fd=system(1000.0,1e7,Some(0.4)).0.prepare().unwrap();
    let mut analytic=system(1000.0,1e7,Some(0.4)).0.prepare_analytic().unwrap();
    let gate=CancelGate::new_clock_free();let initial=reference.stored_energy_j();let mut loss=0.0;
    for _ in 0..2000 {
        let f=reference.step(&[0.0;2],&gate).unwrap();
        fd.step(&[0.0;2],&gate).unwrap();analytic.step(&[0.0;2],&gate).unwrap();
        assert_eq!(reference.state(),fd.state());states_close(reference.state(),analytic.state());
        loss+=f.dissipated_energy_j;assert!(f.balance_residual_j.abs()<1e-10);
    }
    assert!(loss>1e-9);
    assert!((reference.stored_energy_j()+loss-initial).abs()<1e-9);
    assert!(reference.state()[1]<0.0 && reference.state()[3]>0.0);
    let p=m[0].sqrt()*reference.state()[1]+m[1].sqrt()*reference.state()[3];
    assert!((p-m[0]*0.04).abs()<1e-8);
    let before=analytic.state().to_vec();let sources=analytic.material_contacts().to_vec();
    let resumed=analytic.into_reference();assert_eq!(resumed.state(),before);
    assert_eq!(resumed.material_contacts(),sources);
}

#[test]
fn stiffness_density_and_loss_are_independent_source_only_substitutions() {
    let (base,m)=system(1000.0,1e7,None);let initial=base.state().to_vec();
    let mut base=base.prepare_analytic().unwrap();
    let mut stiff=system(1000.0,3e7,None).0.prepare_analytic().unwrap();
    assert_eq!(initial,stiff.state());
    let (heavy,hm)=system(2000.0,1e7,None);close(hm[0]/m[0],2.0);
    close(heavy.material_contacts()[0].source().stiffness_n_per_m_3_2(),
        base.material_contacts()[0].source().stiffness_n_per_m_3_2());
    let mut heavy=heavy.prepare_analytic().unwrap();
    let mut lossy=system(1000.0,1e7,Some(0.4)).0.prepare_analytic().unwrap();
    let gate=CancelGate::new_clock_free();let mut peaks=[0.0_f64;2];
    for _ in 0..2000 {
        base.step(&[0.0;2],&gate).unwrap();stiff.step(&[0.0;2],&gate).unwrap();
        heavy.step(&[0.0;2],&gate).unwrap();lossy.step(&[0.0;2],&gate).unwrap();
        peaks[0]=peaks[0].max(base.material_contact_observation(0).unwrap().normal_force_n);
        peaks[1]=peaks[1].max(stiff.material_contact_observation(0).unwrap().normal_force_n);
    }
    assert!(peaks[1]>1.2*peaks[0]);
    assert!(heavy.state()[3]/hm[1].sqrt()>base.state()[3]/m[1].sqrt()+0.005);
    assert!(lossy.stored_energy_j()<base.stored_energy_j()-1e-9);
}

#[test]
fn failed_solver_trial_and_cancellation_do_not_advance_source_state_or_clock() {
    let mut retry=system(1000.0,1e7,Some(0.4)).0.prepare_analytic().unwrap();
    let mut clean=system(1000.0,1e7,Some(0.4)).0.prepare_analytic().unwrap();
    let gate=CancelGate::new_clock_free();
    for _ in 0..400 {retry.step(&[0.0;2],&gate).unwrap();clean.step(&[0.0;2],&gate).unwrap();}
    let before=retry.state().to_vec();let sources=retry.material_contacts().to_vec();let sample=retry.samples();
    retry.set_iteration_limit(0).unwrap();assert!(retry.step(&[0.0;2],&gate).is_err());
    assert_eq!(retry.state(),before);assert_eq!(retry.samples(),sample);
    let cancelled=CancelGate::new_clock_free();cancelled.request();
    assert!(retry.step(&[0.0;2],&cancelled).is_err());assert!(retry.step(&[f64::NAN,0.0],&gate).is_err());
    assert_eq!(retry.state(),before);assert_eq!(retry.material_contacts(),sources);
    retry.set_iteration_limit(50).unwrap();
    retry.step(&[0.0;2],&gate).unwrap();clean.step(&[0.0;2],&gate).unwrap();
    assert_eq!(retry.state(),clean.state());
}

#[test]
fn incomplete_maps_duplicate_keys_and_initial_source_excursions_refuse() {
    let (bodies,c,_)=pair(1000.0,1e7,None);
    let build=|contacts|ImpactSystem::new_with_material_contacts(bodies.clone(),contacts,vec![],vec![],config());
    assert!(build(vec![c.clone(),c.clone()]).is_err());
    let short=MaterialContact::new("short".into(),c.source().clone(),vec![1.0],0.0).unwrap();
    assert!(build(vec![short]).is_err());
    let crushed=MaterialContact::new("crushed".into(),c.source().clone(),c.weights().to_vec(),-0.001).unwrap();
    assert!(build(vec![crushed]).is_err());
    assert!(MaterialContact::new("zero".into(),c.source().clone(),vec![0.0;2],0.0).is_err());
    let mut fast=bodies;fast[0].initial[0].velocity_m_sqrt_kg_per_s=1.0;
    assert!(ImpactSystem::new_with_material_contacts(fast,vec![c],vec![],vec![],config()).is_err());
}

#[test]
fn material_rate_gate_rejects_candidate_and_preserves_the_accepted_transaction() {
    for analytic in [false,true] {
        let (bodies,c,_)=pair(1000.0,1e7,None);
        let mut envelope=c.source().envelope();envelope.maximum_absolute_rate_m_per_s=0.041;
        let source=SourcePointContact::new(c.source().model().clone(),c.source().geometry(),envelope).unwrap();
        let contact=MaterialContact::new("bounded".into(),source,c.weights().to_vec(),0.0).unwrap();
        let s=ImpactSystem::new_with_material_contacts(bodies,vec![contact],vec![],vec![],config()).unwrap();
        let mut p=if analytic{s.prepare_analytic().unwrap()}else{s.prepare().unwrap()};
        let initial=p.state().to_vec();let gate=CancelGate::new_clock_free();
        // One N/sqrt(kg) is within the force budget. After a few dry steps it
        // exceeds the explicit relative-speed box; no failed endpoint publishes.
        let mut refused=false;
        for _ in 0..100 {
            let before=p.state().to_vec();let sample=p.samples();
            if p.step(&[1.0,0.0],&gate).is_err() {
                assert_eq!(p.state(),before);assert_eq!(p.samples(),sample);refused=true;break;
            }
        }
        assert!(refused);assert_ne!(p.state(),initial);
    }
}

#[path="material_impact/pair.rs"]
mod pair_tests;
