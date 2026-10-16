use std::ops::ControlFlow;
use fs_solver::{CgState, CsrOp, LinearOp};
use fs_solver::op::two_level::{AdditiveTwoLevel, TwoLevelBudget, TwoLevelError};
use fs_sparse::{Coo, Csr, precond::{IdentityPrecond, Precond}};

fn fixture(nc: usize) -> (CsrOp, Csr) {
    let n = 2 * nc + 1;
    let mut a = Coo::new(n, n);
    for i in 0..n {
        a.push(i, i, 2.0);
        if i + 1 < n { a.push(i, i + 1, -1.0); a.push(i + 1, i, -1.0); }
    }
    let mut p = Coo::new(n, nc);
    for c in 0..nc {
        p.push(2*c, c, 0.5); p.push(2*c+1, c, 1.0); p.push(2*c+2, c, 0.5);
    }
    (CsrOp::symmetric(a.assemble()), p.assemble())
}
fn dot(a: &[f64], b: &[f64]) -> f64 { a.iter().zip(b).map(|(a,b)|a*b).sum() }
#[test]
fn additive_action_matches_an_independent_small_coarse_inverse() {
    let (op,p) = fixture(2);
    let pc = AdditiveTwoLevel::new(&op, &[0.5;5], p, TwoLevelBudget::default(),
        |_|ControlFlow::Continue(())).unwrap();
    // The coarse matrix is [[1,-1/2],[-1/2,1]], inverse (4/3)*[[1,1/2],[1/2,1]].
    let r = [1.0,-2.0,3.0,2.0,5.0];
    let rc = [0.5*r[0]+r[1]+0.5*r[2], 0.5*r[2]+r[3]+0.5*r[4]];
    let c = [(4.0/3.0)*(rc[0]+0.5*rc[1]), (4.0/3.0)*(0.5*rc[0]+rc[1])];
    let coarse = [0.5*c[0], c[0], 0.5*(c[0]+c[1]), c[1], 0.5*c[1]];
    let mut actual = [0.0;5]; pc.apply(&r,&mut actual);
    for i in 0..5 { assert!((actual[i]-(0.5*r[i]+coarse[i])).abs()<1e-12); }
    assert_eq!(pc.work().operator_applications,2);
}
#[test]
fn preconditioner_is_fixed_linear_symmetric_positive_and_replays() {
    let (op,p)=fixture(3);
    let pc=AdditiveTwoLevel::new(&op,&[0.5;7],p,TwoLevelBudget::default(),
        |_|ControlFlow::Continue(())).unwrap();
    let x=[1.0,-2.0,3.0,0.0,1.0,-0.5,2.0]; let y=[-1.0,1.0,0.0,2.0,-3.0,0.5,4.0];
    let mut bx=[0.0;7];let mut by=bx;let mut combo=bx;let mut replay=bx;
    pc.apply(&x,&mut bx);pc.apply(&y,&mut by);pc.apply(&x,&mut replay);
    let r=std::array::from_fn::<_,7,_>(|i|2.0*x[i]-3.0*y[i]);pc.apply(&r,&mut combo);
    assert_eq!(bx,replay);assert!(dot(&x,&bx)>0.0);
    assert!((dot(&x,&by)-dot(&y,&bx)).abs()<1e-11);
    for i in 0..7 {assert!((combo[i]-(2.0*bx[i]-3.0*by[i])).abs()<1e-11);}
}
#[test]
fn two_levels_reduce_outer_krylov_work_without_changing_equilibrium() {
    let (op,p)=fixture(31);let n=op.n();
    let rhs:Vec<f64>=(0..n).map(|i|((17*i)%31) as f64/31.0-0.5).collect();
    let pc=AdditiveTwoLevel::new(&op,&vec![0.5;n],p,TwoLevelBudget::default(),
        |_|ControlFlow::Continue(())).unwrap();
    let mut plain=CgState::new(&op,&IdentityPrecond,&rhs);
    assert!(plain.run(&op,&IdentityPrecond,1e-11,300).converged);
    let mut accelerated=CgState::new(&op,&pc,&rhs);
    assert!(accelerated.run(&op,&pc,1e-11,300).converged);
    assert!(accelerated.iters<plain.iters,"{} vs {}",accelerated.iters,plain.iters);
    let mut applied=vec![0.0;n];op.apply(&accelerated.x,&mut applied);
    let residual:Vec<f64>=applied.iter().zip(&rhs).map(|(a,b)|a-b).collect();
    assert!(dot(&residual,&residual)<1e-18*dot(&rhs,&rhs));
    for (x,y) in plain.x.iter().zip(&accelerated.x){assert!((x-y).abs()<1e-8);}
}
#[test]
fn setup_caps_cancellation_bad_diagonal_and_rank_failure_refuse() {
    let (op,p)=fixture(2);let mut calls=0;
    let result=AdditiveTwoLevel::new(&op,&[0.5;5],p,TwoLevelBudget{max_operator_applications:1,..Default::default()},
        |_|{calls+=1;ControlFlow::Continue(())});
    assert!(matches!(result,Err(TwoLevelError::Budget(_))));assert_eq!(calls,1);
    let (op,p)=fixture(2);let mut completed=0;
    let result=AdditiveTwoLevel::new(&op,&[0.5;5],p,TwoLevelBudget::default(),|w|{
        completed=w.operator_applications;
        if completed>0{ControlFlow::Break(())}else{ControlFlow::Continue(())}
    });
    assert!(matches!(result,Err(TwoLevelError::Cancelled)));assert_eq!(completed,1);
    let (op,p)=fixture(2);
    assert!(matches!(AdditiveTwoLevel::new(&op,&[0.0;5],p,TwoLevelBudget::default(),
        |_|ControlFlow::Continue(())),Err(TwoLevelError::Invalid(_))));
    let mut p=Coo::new(5,2);p.push(0,0,1.0);p.push(0,1,1.0);
    assert!(matches!(AdditiveTwoLevel::new(&op,&[0.5;5],p.assemble(),TwoLevelBudget::default(),
        |_|ControlFlow::Continue(())),Err(TwoLevelError::NotPositiveDefinite)));
}
