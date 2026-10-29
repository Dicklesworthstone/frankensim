use super::*;
use crate::visco::GeneralizedMaxwell;

fn diagonal(value:f64)->Tangent6 {
    core::array::from_fn(|i|core::array::from_fn(|j|if i==j {value}else{0.0}))
}
fn close(a:f64,b:f64,tolerance:f64) {
    assert!((a-b).abs()<=tolerance*(1.0+a.abs().max(b.abs())),"{a} != {b}");
}
fn model()->TensorProny {
    let mut a=[[0.0;6];6];a[0]=[2.0,1.0,-1.0,0.5,0.0,0.0];
    let mut b=[[0.0;6];6];b[0]=[0.0,3.0,1.0,-0.5,0.0,1.0];b[1]=[1.0,0.0,0.0,1.0,1.0,0.0];
    TensorProny::new_mandel(diagonal(10.0),vec![
        TensorPronyBranch::from_gram_rows(a,0.01).unwrap(),
        TensorPronyBranch::from_gram_rows(b,0.03).unwrap()],0.1).unwrap()
}

#[test]
fn diagonal_tensor_matches_the_existing_scalar_prony_history() {
    let scalar=GeneralizedMaxwell::new(2e9,vec![(3e9,0.02),(1e9,0.2)]).unwrap();
    let branches=scalar.terms.iter().map(|&(e,tau)| {
        let mut b=[[0.0;6];6];b[0][0]=det::sqrt(e);
        TensorPronyBranch::from_gram_rows(b,tau).unwrap()
    }).collect();
    let law=TensorProny::new_mandel(diagonal(scalar.e_inf),branches,0.1).unwrap();
    let mut tensor=law.history();let mut state=scalar.state();
    for (eps,dt) in [(0.002,0.0),(-0.001,1e-12),(0.003,0.03),(-0.004,0.2),(0.0,2.0),(0.0,10.0)] {
        let expected=scalar.step(&mut state,eps,dt);
        let actual=tensor.advance([eps,0.0,0.0,0.0,0.0,0.0],dt,||false).unwrap();
        close(actual.stress[0],expected,2e-12);
        assert!(actual.stress[1..].iter().all(|s|*s==0.0));
        close(tensor.work(),state.work,2e-12);close(tensor.dissipated(),state.dissipated,2e-12);
        let stored=0.5*scalar.e_inf*eps*eps+state.q.iter().zip(&scalar.terms)
            .map(|(q,(e,_))|0.5*e*q*q).sum::<f64>();
        close(tensor.stored_energy(),stored,2e-12);
        close(tensor.work(),tensor.stored_energy()+tensor.dissipated(),2e-12);
    }
    for omega in [0.0,1.0,50.0,1e4] {
        let (a,b)=law.modulus_mandel(omega).unwrap();let (c,d)=scalar.modulus(omega);
        close(a[0][0],c,2e-14);close(b[0][0],d,2e-14);
    }
}

#[test]
fn noncommuting_rank_deficient_branches_retain_coupled_passive_memory() {
    let law=model();let a=law.branches()[0].stiffness_mandel();let b=law.branches()[1].stiffness_mandel();
    let ab=(0..6).map(|k|a[0][k]*b[k][1]).sum::<f64>();
    let ba=(0..6).map(|k|b[0][k]*a[k][1]).sum::<f64>();
    assert!((ab-ba).abs()>1.0,"fixture must not share one eigentensor basis");
    let mut history=law.history();
    for step in 0..100 {
        let eps=core::array::from_fn(|i|0.003*det::sin((step*(i+1)) as f64*0.17));
        let report=history.advance(eps,0.003,||false).unwrap();
        assert!(report.dissipated>=0.0);assert!(report.stored_energy>=0.0);
        close(report.balance_residual,0.0,2e-13);
        close(history.work(),history.stored_energy()+history.dissipated(),2e-12);
    }
    let (_,loss)=law.modulus_mandel(50.0).unwrap();
    assert!(loss[0][2]<0.0,"negative cross terms must not be clipped");
    for i in 0..6 {for j in 0..6 {assert_eq!(loss[i][j],loss[j][i]);}}
    for seed in 0..32 {
        let v=core::array::from_fn(|i|det::sin((seed*6+i) as f64));
        assert!(dot(&v,&apply(&loss,&v))>=-1e-13);
    }
}

#[test]
fn shear_convention_and_algorithmic_tangent_match_physical_stress() {
    let elastic=crate::IsotropicElastic::new(12e9,0.25,0.1).unwrap();
    let law=TensorProny::from_tensor_shear(elastic.stiffness(),vec![],0.1).unwrap();
    let eps=[0.0,0.0,0.0,0.001,0.0,0.0];
    let report=law.history().advance(eps,0.0,||false).unwrap();
    let mu=elastic.lame().1;
    close(report.stress[3],2.0*mu*eps[3],2e-14);
    close(report.stored_energy,2.0*mu*eps[3]*eps[3],2e-14);
    assert_eq!(report.dissipated,0.0);
    let law=model();let mut initial=law.history();
    initial.advance([0.002,-0.003,0.001,0.004,0.0,0.001],0.0,||false).unwrap();
    for dt in [0.0,1e-10,0.02,1.0] {
        let tangent=law.tangent(dt).unwrap();let eps=[0.001;6];let h=1e-7;
        for col in 0..6 {
            let mut plus=eps;plus[col]+=h;let mut minus=eps;minus[col]-=h;
            let a=initial.clone().advance(plus,dt,||false).unwrap();
            let b=initial.clone().advance(minus,dt,||false).unwrap();
            for row in 0..6 {close((a.stress[row]-b.stress[row])/(2.0*h),tangent[row][col],2e-9);}
        }
    }
}

#[test]
fn harmonic_cycle_work_converges_to_the_full_loss_tensor() {
    let law=model();let omega=50.0;let amplitude=[0.002,-0.001,0.0,0.001,0.0,0.0005];
    let (_,loss)=law.modulus_mandel(omega).unwrap();let v=mandel(&amplitude);
    let expected=core::f64::consts::PI*dot(&v,&apply(&loss,&v));
    assert!(expected>0.0);
    let period=core::f64::consts::TAU/omega;let mut previous=f64::INFINITY;
    for samples in [64,128,256] {
        let mut history=law.history();let mut cycle_work=0.0;
        for k in 1..=12*samples {
            let phase=core::f64::consts::TAU*(k%samples) as f64/samples as f64;
            let strain=amplitude.map(|a|a*det::sin(phase));
            let report=history.advance(strain,period/samples as f64,||false).unwrap();
            if k>11*samples {cycle_work+=report.work;}
        }
        let error=(cycle_work/expected-1.0).abs();
        assert!(error<0.3*previous,"strain-ramp refinement: {error} vs {previous}");
        previous=error;
    }
    assert!(previous<2e-4);
}

#[test]
fn cancellation_and_refused_advances_are_atomic_and_retryable() {
    let law=model();let mut initial=law.history();
    initial.advance([0.001;6],0.02,||false).unwrap();
    let expected=initial.clone().advance([0.002;6],0.03,||false).unwrap();
    for stop in 1..=law.branches().len()+2 {
        let mut actual=initial.clone();let mut polls=0;
        assert_eq!(actual.advance([0.002;6],0.03,||{polls+=1;polls==stop}),Err(TensorPronyError::Cancelled));
        assert_eq!(actual,initial);
        assert_eq!(actual.advance([0.002;6],0.03,||false).unwrap(),expected);
    }
    for (eps,dt) in [([f64::NAN;6],0.1),([0.0;6],-0.1),([1.0;6],0.1),([0.0;6],f64::INFINITY)] {
        let mut actual=initial.clone();assert!(actual.advance(eps,dt,||false).is_err());assert_eq!(actual,initial);
    }
}

#[test]
fn material_admission_and_extreme_relaxation_scales_fail_closed() {
    for c in [diagonal(0.0),diagonal(-1.0),diagonal(f64::NAN)] {
        assert!(TensorProny::new_mandel(c,vec![],0.1).is_err());
    }
    let mut c=diagonal(1.0);c[0][1]=1.0;
    assert!(TensorProny::new_mandel(c,vec![],0.1).is_err());
    assert!(TensorPronyBranch::from_gram_rows(diagonal(1e308),1.0).is_err());
    assert!(TensorPronyBranch::from_gram_rows(diagonal(1.0),0.0).is_err());
    let branch=TensorPronyBranch::from_gram_rows(diagonal(1.0),1e308).unwrap();
    assert!(TensorProny::new_mandel(diagonal(1.0),vec![branch.clone();MAX_BRANCHES+1],0.1).is_err());
    let law=TensorProny::new_mandel(diagonal(1.0),vec![branch],0.1).unwrap();
    let (storage,loss)=law.modulus_mandel(1e308).unwrap();assert_eq!(storage,diagonal(2.0));assert_eq!(loss,diagonal(0.0));
    for omega in [-1.0,f64::NAN,f64::INFINITY] {assert!(law.modulus_mandel(omega).is_err());}
    let branch=TensorPronyBranch::from_gram_rows(diagonal(1.0),f64::MIN_POSITIVE).unwrap();
    let law=TensorProny::new_mandel(diagonal(1.0),vec![branch],0.1).unwrap();
    let mut history=law.history();history.advance([0.001;6],0.0,||false).unwrap();
    let report=history.advance([0.001;6],f64::MAX,||false).unwrap();
    close(report.dissipated,report.stored_energy,2e-14);
    assert_eq!(report.work,0.0);close(report.balance_residual,0.0,2e-14);
}

#[test]
fn original_mandel_and_tensor_shear_histories_share_source_strain_admission() {
    let law = model();
    let strain = [0.0, 0.0, 0.0, 0.07, 0.0, 0.0];
    let mandel = [0.0, 0.0, 0.0, 0.07 * core::f64::consts::SQRT_2, 0.0, 0.0];
    let mut original = law.state();
    let mut factored = law.history();
    let stress = law.step(&mut original, mandel, 0.01).unwrap();
    let other = factored.advance(strain, 0.01, || false).unwrap();
    close(stress[3], other.stress[3] * core::f64::consts::SQRT_2, 2e-12);
    let original_before = original.clone();
    let factored_before = factored.clone();
    for value in [0.071, f64::NAN, f64::INFINITY] {
        let tensor = [0.0, 0.0, 0.0, value, 0.0, 0.0];
        let mandel = [0.0, 0.0, 0.0, value * core::f64::consts::SQRT_2, 0.0, 0.0];
        assert!(law.step(&mut original, mandel, 0.01).is_err());
        assert!(factored.advance(tensor, 0.01, || false).is_err());
        assert_eq!(original, original_before);
        assert_eq!(factored, factored_before);
    }
}
