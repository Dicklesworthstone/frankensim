use super::*;
use crate::{PortHamiltonian, QuadraticStorage, step};

fn budget() -> PortExchangeBudget {
    PortExchangeBudget {max_left:32,max_right:1024,max_setup_terms:70_000_000,maximum_dt_coupling:1.0}
}
fn near(a: &[f64], b: &[f64]) {
    assert_eq!(a.len(),b.len());
    for (&x,&y) in a.iter().zip(b) {assert!((x-y).abs()<2e-12,"{x:e} != {y:e}");}
}
fn mul(a: &[f64], rows: usize, cols: usize, x: &[f64]) -> Vec<f64> {
    (0..rows).map(|i|(0..cols).map(|j|a[i*cols+j]*x[j]).sum()).collect()
}

#[test]
fn exact_rank_one_rotation_retains_the_unconnected_complement_for_both_shapes() {
    let dt=0.03;let angle: f64=5.0*dt;let (s,c)=(angle.sin(),angle.cos());
    // A wide L uses the transposed SVD; the physical signs must not transpose.
    let mut wide=PreparedPortExchange::new(2,1,&[3.0,-4.0],dt,budget()).unwrap();
    let (mut x,mut y)=([0.1,0.2],[0.3]);let projected=0.6*x[0]-0.8*x[1];
    let new=c*projected-s*y[0];let expected=[x[0]+0.6*(new-projected),x[1]-0.8*(new-projected)];
    let expected_y=[s*projected+c*y[0]];
    assert!(wide.apply(&mut x,&mut y).unwrap().energy_defect.abs()<1e-15);
    near(&x,&expected);near(&y,&expected_y);
    let mut tall=PreparedPortExchange::new(1,2,&[3.0,-4.0],dt,budget()).unwrap();
    let (mut x,mut y)=([0.3],[0.1,0.2]);let projected=0.6*y[0]-0.8*y[1];
    let new=s*x[0]+c*projected;
    let expected=[y[0]+0.6*(new-projected),y[1]-0.8*(new-projected)];
    let expected_x=[c*x[0]-s*projected];
    tall.apply(&mut x,&mut y).unwrap();near(&x,&expected_x);near(&y,&expected);
}

#[test]
fn full_quadratic_phs_refines_to_the_collective_continuous_flow() {
    let l=[0.7,-0.4,0.3,0.9,-0.2,0.5];let n=5;
    let mut j=vec![0.0;n*n];let mut q=vec![0.0;n*n];
    for i in 0..n {q[i*n+i]=1.0;}
    for row in 0..3 {for col in 0..2 {j[col*n+2+row]=-l[row*2+col];j[(2+row)*n+col]=l[row*2+col];}}
    let sys=PortHamiltonian::new(n,0,j,vec![0.0;n*n],vec![],Box::new(QuadraticStorage::new(q,n).unwrap())).unwrap();
    let x0=[0.2,-0.1,0.3,-0.4,0.1];let horizon=0.2;
    let mut exact=PreparedPortExchange::new(2,3,&l,horizon,budget()).unwrap();
    let (mut a,mut b)=([x0[0],x0[1]],[x0[2],x0[3],x0[4]]);exact.apply(&mut a,&mut b).unwrap();
    let expected:Vec<_>=a.into_iter().chain(b).collect();let mut errors=Vec::new();
    for count in [16,32,64] {
        let mut x=x0.to_vec();
        for _ in 0..count {x=step(&sys,&x,&[],horizon/count as f64).unwrap().x;}
        errors.push(x.iter().zip(&expected).map(|(x,y)|(x-y).abs()).fold(0.0_f64,f64::max));
    }
    assert!(errors[0]>1e-8 && errors[1]<0.27*errors[0] && errors[2]<0.27*errors[1],"{errors:?}");
}

#[test]
fn arbitrary_orthogonal_coordinate_changes_preserve_the_full_power_exchange() {
    let l=[1.0,-0.4,0.2,0.7,-0.3,0.6];let q=[0.6,-0.8,0.8,0.6];
    let direction=[1.0/3.0,2.0/3.0,2.0/3.0];
    let r:Vec<_>=(0..9).map(|i|(if i/3==i%3 {1.0}else{0.0})-2.0*direction[i/3]*direction[i%3]).collect();
    let mut transformed=vec![0.0;6];
    for i in 0..3 {for j in 0..2 {for a in 0..3 {for b in 0..2 {
        transformed[i*2+j]+=r[i*3+a]*l[a*2+b]*q[j*2+b];
    }}}}
    let (mut a,mut b)=([0.2,-0.3],[0.4,0.1,-0.2]);
    let mut ar=mul(&q,2,2,&a);let mut br=mul(&r,3,3,&b);
    let mut original=PreparedPortExchange::new(2,3,&l,0.05,budget()).unwrap();
    let mut rotated=PreparedPortExchange::new(2,3,&transformed,0.05,budget()).unwrap();
    for _ in 0..20 {original.apply(&mut a,&mut b).unwrap();rotated.apply(&mut ar,&mut br).unwrap();}
    near(&ar,&mul(&q,2,2,&a));near(&br,&mul(&r,3,3,&b));
}

#[test]
fn rank_deficient_and_zero_loads_keep_all_history_and_do_not_create_dissipation() {
    let mut rank_one=PreparedPortExchange::new(1,3,&[1.0,2.0,2.0],0.1,budget()).unwrap();
    let (mut a,mut b)=([0.0],[2.0,-1.0,0.0]);
    rank_one.apply(&mut a,&mut b).unwrap();assert_eq!(a,[0.0]);assert_eq!(b,[2.0,-1.0,0.0]);
    let mut deficient=PreparedPortExchange::new(2,3,&[1.0,2.0,2.0,4.0,0.0,0.0],0.01,budget()).unwrap();
    let (mut a,mut b)=([0.2,-0.1],[0.0,0.0,0.7]);
    deficient.apply(&mut a,&mut b).unwrap();near(&a,&[0.2,-0.1]);near(&b,&[0.0,0.0,0.7]);
    let mut zero=PreparedPortExchange::new(2,3,&[0.0;6],0.1,budget()).unwrap();
    let (mut a,mut b)=([-0.0,0.3],[0.2,-0.0,0.7]);let before=(a.map(f64::to_bits),b.map(f64::to_bits));
    assert_eq!(zero.apply(&mut a,&mut b).unwrap().energy_defect,0.0);
    assert_eq!((a.map(f64::to_bits),b.map(f64::to_bits)),before);
}

#[test]
fn every_cancellation_boundary_is_transactional_and_retries_exactly() {
    let l=[1.0,-0.4,0.2,0.7,-0.3,0.6];
    let mut work=PreparedPortExchange::new(2,3,&l,0.05,budget()).unwrap();
    let (a,b)=([0.2,-0.3],[0.4,0.1,-0.2]);let (mut expected_a,mut expected_b)=(a,b);let mut calls=0;
    let record=work.apply_controlled(&mut expected_a,&mut expected_b,||{calls+=1;false}).unwrap();
    for stop in 1..=calls {
        let (mut x,mut y)=(a,b);let mut count=0;
        assert_eq!(work.apply_controlled(&mut x,&mut y,||{count+=1;count==stop}),Err(PreparedStepError::Cancelled));
        assert_eq!(x,a);assert_eq!(y,b);
        assert_eq!(work.apply(&mut x,&mut y).unwrap(),record);assert_eq!(x,expected_a);assert_eq!(y,expected_b);
    }
    let (mut x,mut y)=(a,b);x[0]=f64::INFINITY;
    assert!(work.apply(&mut x,&mut y).is_err());assert_eq!(y,b);assert_eq!(x[0],f64::INFINITY);
}

#[test]
fn cold_work_rate_and_finite_limits_refuse_without_silent_rank_or_port_changes() {
    for dt in [0.0,-1.0,f64::NAN,10.0] {assert!(PreparedPortExchange::new(1,1,&[1.0],dt,budget()).is_err());}
    for bad in [f64::NAN,f64::INFINITY] {assert!(PreparedPortExchange::new(1,1,&[bad],0.01,budget()).is_err());}
    let mut b=budget();b.max_setup_terms=59;assert!(PreparedPortExchange::new(1,1,&[1.0],0.01,b).is_err());
    assert!(PreparedPortExchange::new(33,1,&[1.0;33],0.01,budget()).is_err());
    assert!(PreparedPortExchange::new(1,usize::MAX,&[],0.01,budget()).is_err());
    assert!(PreparedPortExchange::new(2,1,&[1.0],0.01,budget()).is_err());
    assert!(PreparedPortExchange::new(1,1,&[f64::from_bits(1)],0.1,budget()).is_err());
}
