//! Reuse the existing actual mesh/eigensolve fixtures and physical test seams.
use super::*;
use super::super::BodyTangents;

fn compare(s:&ImpactSystem,cache:&mut BodyTangents,x:&[f64],d:&[f64]) {
    let mut direct=vec![0.0;x.len()];let mut batched=vec![0.0;x.len()];
    assert!(s.hessian_vector(x,d,&mut direct));
    assert!(s.hessian_vector_cached(x,d,&mut batched,cache));
    let scale=direct.iter().chain(&batched).fold(1.0_f64,|a,b|a.max(b.abs()));
    for (a,b) in direct.iter().zip(&batched) {assert!((a-b).abs()<3e-12*scale,"{a:e} != {b:e}");}
}

fn sweep(s:&ImpactSystem,cache:&mut BodyTangents,x:&[f64]) {
    for col in 0..x.len() {let mut d=vec![0.0;x.len()];d[col]=1.0;compare(s,cache,x,&d);}
}

#[test]
fn g0_g3_real_shell_and_film_assemble_once_per_body_point_not_once_per_column() {
    let a=body(shell());let split=a.potential.count();let b=body(membrane());let total=split+b.potential.count();
    let mut cfg=config();cfg.dt_s=1e-9; // Admits every supplied shell mode, no band deletion.
    let s=ImpactSystem::new(vec![a,b],vec![],vec![],vec![],cfg).unwrap();
    let mut cache=BodyTangents::new(Rc::clone(&s.mechanical));
    let mut x:Vec<_>=(0..2*total).map(|i|if i%2==0 {1e-5*((i+1) as f64).sin()}else{0.03}).collect();
    sweep(&s,&mut cache,&x);assert_eq!(cache.counts(),(2,total));
    // Momentum is deliberately absent from this potential's cache key.
    for i in 0..total {x[2*i+1]*=-2.0;}
    sweep(&s,&mut cache,&x);assert_eq!(cache.counts(),(2,2*total));
    // Change only the second body: the first body's exact point is reusable.
    x[2*split]+=1e-6;
    sweep(&s,&mut cache,&x);assert_eq!(cache.counts(),(3,3*total));
    let d:Vec<_>=(0..x.len()).map(|i|(i as f64+0.7).cos()).collect();compare(&s,&mut cache,&x,&d);
    assert_eq!(cache.counts().0,3);
    let accepted=s.state().to_vec();let mut malformed=x.clone();malformed[0]=f64::NAN;
    assert!(!s.hessian_vector_cached(&malformed,&d,&mut vec![0.0;x.len()],&mut cache));
    let mut overflow=x.clone();overflow[2*split]=1e200;
    assert!(!s.hessian_vector_cached(&overflow,&d,&mut vec![0.0;x.len()],&mut cache));
    compare(&s,&mut cache,&x,&d);assert_eq!(s.state(),accepted);assert_eq!(s.samples(),0);
    // A different owner with an identical shape cannot reuse another model's matrix.
    let other=ImpactSystem::new(vec![body(BodyPotential::Linear(vec![800.0;total]))],
        vec![],vec![],vec![],cfg).unwrap();
    let mut sentinel=vec![123.0;x.len()];
    assert!(!other.hessian_vector_cached(&x,&d,&mut sentinel,&mut cache));
    assert!(sentinel.iter().all(|v|*v==123.0));
}

#[test]
fn g0_cached_body_curvature_does_not_freeze_contact_felt_or_kelvin_history() {
    let membrane=body(membrane());let n=membrane.potential.count();
    let mut weights=vec![0.0;n];weights[0]=1.0;
    let contact=Obstacle::new(weights.clone(),1,n,vec![0.0],vec![1.0],2e7,1.5,
        "batched nonlinear material regression".into()).unwrap();
    let s=ImpactSystem::new(vec![membrane],vec![contact],vec![pad(weights,0.0012)],vec![],config()).unwrap();
    let mut cache=BodyTangents::new(Rc::clone(&s.mechanical));
    let mut x=s.state().to_vec();x[0]=1e-4;
    let mut d=vec![0.0;x.len()];d[0]=1.0;d[2*n]=0.3;
    compare(&s,&mut cache,&x,&d);let builds=cache.counts().0;
    let mut before=vec![0.0;x.len()];assert!(s.hessian_vector_cached(&x,&d,&mut before,&mut cache));
    // An explicit test-only history change at IDENTICAL body displacement.
    // The live material owner must still determine the unloading tangent.
    let previous=s.histories.borrow()[0].clone();
    s.histories.borrow_mut()[0]=s.pads[0].spec.law.update_state(0.4,&previous);
    compare(&s,&mut cache,&x,&d);
    let mut after=before.clone();assert!(s.hessian_vector_cached(&x,&d,&mut after,&mut cache));
    assert_eq!(cache.counts().0,builds);
    assert!(before.iter().zip(&after).any(|(a,b)|(a-b).abs()>1e-4));
    // Kelvin memory also changes pad strain without changing the body point.
    x[2*n]=0.02;compare(&s,&mut cache,&x,&d);assert_eq!(cache.counts().0,builds);
    // Different displacement detaches contact and forces a fresh body matrix.
    x[0]=-1e-4;compare(&s,&mut cache,&x,&d);assert_eq!(cache.counts().0,builds+1);
}

fn membrane_strike()->ImpactSystem {
    let film=body(membrane());let n=1+film.potential.count();
    let (stick,w)=ImpactBody::free_mass(0.02,-0.0002,0.8).unwrap();
    let mut row=vec![0.0;n];row[0]=w;row[1]=-2.0;
    let contact=Obstacle::new(row,1,n,vec![0.0],vec![1.0],2e7,1.5,
        "struck geometric film; synthetic contact and felt declarations".into()).unwrap();
    let mut weights=vec![0.0;n];weights[1]=2.0;
    let mut areas=vec![0.0;n];areas[1]=0.02;
    ImpactSystem::new(vec![stick,film],vec![contact],vec![pad(weights,0.0006)],
        vec![VolumeSpring {bulk_modulus_pa:1.4e5,volume_m3:0.01,areas}],config()).unwrap()
}

#[test]
fn g3_actual_struck_film_keeps_motion_energy_and_felt_under_batched_analytic_newton() {
    let mut direct=membrane_strike().prepare_analytic().unwrap();direct.set_batched_body_tangents(false);
    let mut batched=membrane_strike().prepare_analytic().unwrap();let gate=CancelGate::new_clock_free();
    let initial=batched.stored_energy_j();let mut net=0.0;let mut assemblies=0;let mut actions=0;
    let mut stretch=0.0_f64;let mut motion=0.0_f64;
    for tick in 0..600 {
        let mut force=[0.0;4];if tick<40 {force[0]=0.02;}
        let old=direct.step(&force,&gate).unwrap();let new=batched.step(&force,&gate).unwrap();
        assert_eq!(old.sample,new.sample);assert_eq!(old.time_s.to_bits(),new.time_s.to_bits());
        assert!(new.balance_residual_j.abs()<1e-8);
        net+=new.supplied_work_j-new.dissipated_energy_j;
        for (a,b) in direct.state().iter().zip(batched.state()) {
            assert!((a-b).abs()<3e-6*a.abs().max(b.abs()).max(1e-4),"tick {tick}: {a:e} != {b:e}");
        }
        let counts=batched.body_tangent_work_counts();assemblies+=counts.0;actions+=counts.1;
        assert_eq!(direct.body_tangent_work_counts(),(0,0));
        stretch=stretch.max(batched.membrane_observation(1).unwrap().stretching_energy_j);
        motion=motion.max(batched.state()[2].abs());
    }
    assert!(assemblies>0 && actions>=3*assemblies,"one film matrix serves its three Newton displacement columns");
    assert!(stretch>0.0 && motion>1e-7);
    assert!(batched.state().last().unwrap().abs()>0.0,"Kelvin history remains dynamical");
    assert!((batched.stored_energy_j()-initial-net).abs()<1e-7);
    assert!((batched.felt_observation(0).unwrap().0-direct.felt_observation(0).unwrap().0).abs()<1e-7);
}

#[test]
fn g4_batched_refusal_switching_and_repreparation_preserve_exact_retry() {
    let mut trial=membrane_strike().prepare_analytic().unwrap();let mut clean=membrane_strike().prepare_analytic().unwrap();
    let gate=CancelGate::new_clock_free();
    for _ in 0..140 {trial.step(&[0.0;4],&gate).unwrap();clean.step(&[0.0;4],&gate).unwrap();}
    let before=trial.state().to_vec();let history=trial.felt_history(0);
    let cancel=CancelGate::new_clock_free();cancel.request();
    assert!(matches!(trial.step(&[0.0;4],&cancel),Err(ImpactError::Cancelled)));
    trial.set_iteration_limit(0).unwrap();assert!(trial.step(&[0.0;4],&gate).is_err());
    trial.set_batched_body_tangents(false);trial.set_batched_body_tangents(true);
    assert_eq!(trial.state(),before);assert_eq!(trial.felt_history(0),history);assert_eq!(trial.samples(),140);
    trial.set_iteration_limit(50).unwrap();
    // Reprepare drops every numerical cache but preserves the authoritative state.
    let mut trial=trial.into_reference().prepare_analytic().unwrap();
    assert_eq!(trial.state(),before);assert_eq!(trial.felt_history(0),history);
    for _ in 0..60 {
        trial.step(&[0.0;4],&gate).unwrap();clean.step(&[0.0;4],&gate).unwrap();
        assert!(trial.state().iter().zip(clean.state()).all(|(a,b)|a.to_bits()==b.to_bits()));
        assert_eq!(trial.felt_history(0),clean.felt_history(0));
    }
}
