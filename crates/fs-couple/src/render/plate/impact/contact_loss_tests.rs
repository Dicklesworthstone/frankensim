//! Joint-contact regressions; authored mechanics are not a calibrated specimen.
use super::*;
use fs_phs::QuadraticStorage;

#[test]
fn contact_loss_port_matches_existing_law_is_passive_and_has_exact_tangent() {
    let ob=Obstacle::new(vec![2.0,-3.0],1,2,vec![0.01],vec![0.5],40.0,1.5,
        "authored two-body contact".into()).unwrap().with_internal_loss(0.4).unwrap();
    let mut q=vec![0.0;25];for i in 0..5 {q[5*i+i]=1.0;}
    let contact=ContactStorage::new(Box::new(QuadraticStorage::new(q,5).unwrap()),2,vec![ob.clone()]).unwrap();
    let x=[0.02,0.0,-0.01,0.0,0.7];let dx=[0.1,0.2,0.03,0.4,0.8];let de=[0.2,-0.1,0.3,0.2,0.9];
    for e in [[1.0,0.2,2.0,-0.1,3.0],[1.0,-4.0,2.0,1.0,3.0]] {
        let (mut out,mut tangent)=([0.0;5],[0.0;5]);
        assert!(contact.dissipative_flow_into(&x,&e,&mut out));
        let expected=ob.dissipative_modal_forces(2,&x,&[e[1],e[3]]);
        assert_eq!(out[1],-expected[0]);assert_eq!(out[3],-expected[1]);
        assert_eq!([out[0],out[2],out[4]],[0.0;3]);
        assert!(out.iter().zip(e).map(|(f,e)|f*e).sum::<f64>()>0.0);
        assert!(contact.dissipative_flow_tangent_into(&x,&e,&dx,&de,&mut tangent));
        let h=1e-6;let mut xp=x;let mut xm=x;let mut ep=e;let mut em=e;
        for i in 0..5 {xp[i]+=h*dx[i];xm[i]-=h*dx[i];ep[i]+=h*de[i];em[i]-=h*de[i];}
        let (mut plus,mut minus)=([0.0;5],[0.0;5]);
        assert!(contact.dissipative_flow_into(&xp,&ep,&mut plus));
        assert!(contact.dissipative_flow_into(&xm,&em,&mut minus));
        for i in 0..5 {assert!((tangent[i]-(plus[i]-minus[i])/(2.0*h)).abs()<1e-8);}
    }
    let mut out=[123.0;5];assert!(!contact.dissipative_flow_into(&x[..3],&x,&mut out));assert_eq!(out,[123.0;5]);
    let mut bad=x;bad[4]=f64::NAN;assert!(!contact.dissipative_flow_into(&bad,&x,&mut out));
    let separated=[-1.0,0.0,0.0,0.0,0.0];
    assert!(contact.dissipative_flow_into(&separated,&x,&mut out));assert_eq!(out,[0.0;5]);
}

fn config()->ImpactConfig {ImpactConfig {dt_s:2e-6,max_steps:2000,maximum_energy_j:1.0,
    energy_absolute_tolerance_j:1e-10,energy_relative_tolerance:1e-7,maximum_generalized_force:1e4}}
fn impact(chi:f64)->ImpactSystem {
    let (stick,w)=ImpactBody::free_mass(0.02,-0.0002,0.8).unwrap();
    let target=ImpactBody {potential:BodyPotential::Linear(vec![2.0*core::f64::consts::PI*500.0]),
        initial:vec![ModalAcousticState::default()],damping_per_s:vec![0.0]};
    let ob=Obstacle::new(vec![w,-2.0],1,2,vec![0.0],vec![1.0],2e7,1.5,
        "authored dissipative Hertz impact".into()).unwrap().with_internal_loss(chi).unwrap();
    ImpactSystem::new(vec![stick,target],vec![ob],vec![],vec![],config()).unwrap()
}

#[test]
fn lossy_contact_rebounds_in_one_joint_solve_with_reference_and_analytic_parity() {
    let gate=CancelGate::new_clock_free();let mut reference=impact(0.8);
    let mut fd=impact(0.8).prepare().unwrap();let mut analytic=impact(0.8).prepare_analytic().unwrap();
    let mut elastic=impact(0.0).prepare_analytic().unwrap();let initial=fd.stored_energy_j();
    let mut loss=0.0;let mut moved=false;let mut rebound=false;
    for _ in 0..1200 {
        let a=reference.step(&[0.0;2],&gate).unwrap();let b=fd.step(&[0.0;2],&gate).unwrap();
        let c=analytic.step(&[0.0;2],&gate).unwrap();elastic.step(&[0.0;2],&gate).unwrap();
        assert_eq!(reference.state(),fd.state());assert_eq!(a.sample,b.sample);assert_eq!(b.sample,c.sample);
        for (x,y) in fd.state().iter().zip(analytic.state()) {assert!((x-y).abs()<1e-8);}
        assert!(b.dissipated_energy_j>=0.0);assert!(b.balance_residual_j.abs()<1e-8);
        loss+=b.dissipated_energy_j;moved|=fd.state()[2].abs()>1e-6;rebound|=fd.state()[1]<0.0;
    }
    assert!(moved && rebound && loss>1e-5);
    assert!((fd.stored_energy_j()-initial+loss).abs()<1e-7);
    assert!(fd.stored_energy_j()<elastic.stored_energy_j()-1e-5);
}

#[test]
fn lossy_contact_refusals_and_preparation_preserve_exact_state_and_clock() {
    let gate=CancelGate::new_clock_free();let mut retried=impact(0.8).prepare_analytic().unwrap();
    let mut clean=impact(0.8).prepare_analytic().unwrap();
    for _ in 0..160 {retried.step(&[0.0;2],&gate).unwrap();clean.step(&[0.0;2],&gate).unwrap();}
    let before=retried.state().to_vec();let sample=retried.samples();
    retried.set_iteration_limit(0).unwrap();assert!(retried.step(&[0.0;2],&gate).is_err());
    assert_eq!(retried.state(),before);assert_eq!(retried.samples(),sample);
    retried.set_iteration_limit(50).unwrap();
    let cancel=CancelGate::new_clock_free();cancel.request();assert!(retried.step(&[0.0;2],&cancel).is_err());
    assert!(retried.step(&[f64::NAN,0.0],&gate).is_err());assert_eq!(retried.state(),before);
    retried.step(&[0.0;2],&gate).unwrap();clean.step(&[0.0;2],&gate).unwrap();assert_eq!(retried.state(),clean.state());
    let before=retried.state().to_vec();let sample=retried.samples();
    let mut resumed=retried.into_reference().prepare_analytic().unwrap();assert_eq!(resumed.state(),before);
    assert_eq!(resumed.samples(),sample);resumed.step(&[0.0;2],&gate).unwrap();clean.step(&[0.0;2],&gate).unwrap();
    assert_eq!(resumed.state(),clean.state());
}

#[test]
fn distributed_wire_sized_nonlinear_basis_is_retained_and_bounded() {
    let body=|n|ImpactBody {potential:BodyPotential::Linear(vec![0.0;n]),
        initial:vec![ModalAcousticState::default();n],damping_per_s:vec![0.0;n]};
    let system=ImpactSystem::new(vec![body(200)],vec![],vec![],vec![],config()).unwrap();
    assert_eq!(system.state().len(),400);
    assert!(ImpactSystem::new(vec![body(MAX_IMPACT_MODES+1)],vec![],vec![],vec![],config()).is_err());
}
