use super::*;
use crate::{PortHamiltonian, StepWorkspace, Storage};

fn matrix(base:&[f64],u:&[f64],v:&[f64])->Vec<f64> {
    let n=u.len();(0..n*n).map(|a|base[a]+u[a/n]*v[a%n]).collect()
}
fn fixture(n:usize,pairs:&[[usize;2]],singles:&[usize])->(Condensation,Vec<f64>,Vec<f64>) {
    let mut plan=Condensation::with_singletons(n,pairs,singles).unwrap();
    let mut base=vec![0.0;n*n];
    for row in 0..n {for col in 0..n {
        if plan.owner[row]==usize::MAX || plan.owner[col]==usize::MAX || plan.owner[row]==plan.owner[col] {
            base[row*n+col]=if row==col {3.0+row as f64*0.1}
                else {0.017*((row*7+col*11)%17) as f64-0.12};
        }
    }}
    for i in 0..n {plan.left[i]=0.023*(i as f64-2.0);plan.right[i]=0.031*((i*3)%11) as f64-0.1;}
    let rhs=(0..n).map(|i|0.2*(i as f64+1.0)).collect();
    (plan,base,rhs)
}

#[test]
fn mixed_noncontiguous_leaves_match_full_lu() {
    for pairs in 0..=5 {for singles in 1..=7 {
        let n=3+2*pairs+singles;
        let pair_indices:Vec<_>=(0..pairs).map(|i|[n-1-i,1+i]).collect();
        let singleton_indices:Vec<_>=(0..singles).map(|i|1+pairs+i).collect();
        let (mut plan,base,rhs)=fixture(n,&pair_indices,&singleton_indices);
        assert_eq!(plan.dimension(),4);
        assert_eq!(plan.responses.len(),(n-3)*5);
        let mut expected=vec![0.0;n];let mut actual=expected.clone();
        LuWorkspace::new(n).unwrap().solve_into(&matrix(&base,&plan.left,&plan.right),&rhs,&mut expected).unwrap();
        assert!(plan.solve(&base,&rhs,&mut actual,&mut ||Ok(())).unwrap());
        for (a,b) in actual.iter().zip(expected) {assert!((a-b).abs()<2e-12*b.abs().max(1.0));}
    }}
}

#[test]
fn scalar_only_and_all_eliminated_layouts_keep_the_rank_one_correction() {
    for n in [1,2,9,33] {
        let singles:Vec<_>=(0..n).rev().collect();
        let (mut plan,base,rhs)=fixture(n,&[],&singles);
        assert_eq!(plan.dimension(),1);
        let mut expected=vec![0.0;n];let mut actual=expected.clone();
        LuWorkspace::new(n).unwrap().solve_into(&matrix(&base,&plan.left,&plan.right),&rhs,&mut expected).unwrap();
        assert!(plan.solve(&base,&rhs,&mut actual,&mut ||Ok(())).unwrap());
        for (a,b) in actual.iter().zip(expected) {assert!((a-b).abs()<2e-12*b.abs().max(1.0));}
    }
    // The retained base is singular; only the full Gonzalez border makes it invertible.
    let mut plan=Condensation::with_singletons(3,&[],&[2,1]).unwrap();
    plan.left[0]=2.0;plan.right[0]=0.5;
    let mut answer=[17.0;3];
    assert!(plan.solve(&[0.,0.,0.,0.,2.,0.,0.,0.,4.],&[3.,4.,8.],&mut answer,&mut ||Ok(())).unwrap());
    assert_eq!(answer,[3.,2.,2.]);
}

#[test]
fn any_cross_leaf_coefficient_or_singular_scalar_refuses_without_publication() {
    let (mut plan,mut base,rhs)=fixture(6,&[[5,1]],&[4,2]);
    let mut output=[17.0;6];
    for (row,col) in [(1,4),(4,1),(4,2),(2,4)] {
        base[row*6+col]=f64::from_bits(1);
        assert!(!plan.solve(&base,&rhs,&mut output,&mut ||Ok(())).unwrap());
        assert_eq!(output,[17.0;6]);base[row*6+col]=0.0;
    }
    let saved=base[4*6+4];base[4*6+4]=0.0;
    assert!(!plan.solve(&base,&rhs,&mut output,&mut ||Ok(())).unwrap());assert_eq!(output,[17.0;6]);
    base[4*6+4]=saved;
    assert!(plan.solve(&base,&rhs,&mut output,&mut ||Ok(())).unwrap());
}

#[test]
fn every_mixed_solve_cancellation_boundary_preserves_output_and_retry() {
    let (mut first,base,rhs)=fixture(6,&[[5,1]],&[4,2]);
    let mut expected=[0.0;6];let mut polls=0;
    assert!(first.solve(&base,&rhs,&mut expected,&mut ||{polls+=1;Ok(())}).unwrap());
    for stop in 1..=polls {
        let (mut plan,_,_)=fixture(6,&[[5,1]],&[4,2]);let mut calls=0;let mut out=[17.0;6];
        assert_eq!(plan.solve(&base,&rhs,&mut out,&mut ||{
            calls+=1;if calls==stop {Err(PreparedStepError::Cancelled)} else {Ok(())}
        }).unwrap_err(),PreparedStepError::Cancelled);
        assert_eq!(out,[17.0;6]);
        assert!(plan.solve(&base,&rhs,&mut out,&mut ||Ok(())).unwrap());assert_eq!(out,expected);
    }
}

struct Quartic;
impl Storage for Quartic {
    fn hamiltonian(&self,x:&[f64])->f64 {
        let s=x[0]-x[2];
        (0..2).map(|i|0.5*((4.0+i as f64)*x[2*i].powi(2)+x[2*i+1].powi(2))).sum::<f64>()+2.5*s.powi(4)
    }
    fn gradient(&self,x:&[f64],out:&mut[f64]) {
        for i in 0..2 {out[2*i]=(4.0+i as f64)*x[2*i];out[2*i+1]=x[2*i+1];}
        let f=10.0*(x[0]-x[2]).powi(3);out[0]+=f;out[2]-=f;
    }
}
fn hessian(x:&[f64],d:&[f64],out:&mut[f64])->bool {
    for i in 0..2 {out[2*i]=(4.0+i as f64)*d[2*i];out[2*i+1]=d[2*i+1];}
    let f=30.0*(x[0]-x[2]).powi(2)*(d[0]-d[2]);out[0]+=f;out[2]-=f;true
}
fn system(drag:f64)->PortHamiltonian {
    let mut j=vec![0.0;16];let mut r=j.clone();let mut g=vec![0.0;4];g[1]=1.0;
    for i in 0..2 {let q=2*i;j[q*4+q+1]=1.0;j[(q+1)*4+q]=-1.0;r[(q+1)*4+q+1]=drag;}
    PortHamiltonian::new(4,1,j,r,g,Box::new(Quartic)).unwrap()
}

#[test]
fn coupled_nonlinear_motion_work_and_loss_match_dense_newton() {
    for drag in [0.0,0.3] {
        let sys=system(drag);let mut dense=StepWorkspace::new(&sys).unwrap();
        let mut reduced=StepWorkspace::new(&sys).unwrap();reduced.set_condensed_leaves(&[],&[1,3]).unwrap();
        assert_eq!(reduced.condensed_dimension(),Some(3));
        let mut x=[0.3,-0.6,-0.1,0.2];let mut reference=x;let initial=sys.hamiltonian(&x);
        let mut work=0.0;let mut loss=0.0;let mut solves=0;
        for step in 0..200 {
            let force=[if step<80 {0.2} else {0.0}];
            let (mut a,mut b)=([0.0;4],[0.0;4]);let (mut ya,mut yb)=([0.0;1],[0.0;1]);
            let record=reduced.step_into_analytic(&sys,&x,&force,0.003,&mut a,&mut ya,&hessian).unwrap();
            dense.step_into_analytic(&sys,&reference,&force,0.003,&mut b,&mut yb,&hessian).unwrap();
            for (a,b) in a.iter().zip(b) {assert!((a-b).abs()<2e-9);}
            assert!((ya[0]-yb[0]).abs()<2e-9);assert!(record.balance_residual().abs()<1e-10);
            assert_eq!(reduced.linear_solve_counts(),(record.newton_iters,0));
            assert!(record.dissipated>=0.0);solves+=record.newton_iters;
            work+=record.supplied;loss+=record.dissipated;x=a;reference=b;
        }
        assert!(solves>0);assert!((sys.hamiltonian(&x)-initial+loss-work).abs()<1e-9);
    }
}

#[test]
fn changed_momentum_coupling_uses_full_equation_not_independent_scalars() {
    let mut sys=system(0.3);let mut reduced=StepWorkspace::new(&sys).unwrap();
    reduced.set_condensed_leaves(&[],&[1,3]).unwrap();
    // A passive cross-velocity damper added after preparing the scalar layout.
    sys.r[1*4+3]=-0.1;sys.r[3*4+1]=-0.1;
    let x=[0.3,-0.6,-0.1,0.2];let (mut a,mut b)=([0.0;4],[0.0;4]);let (mut ya,mut yb)=([0.0;1],[0.0;1]);
    let mut dense=StepWorkspace::new(&sys).unwrap();
    let record=reduced.step_into_analytic(&sys,&x,&[0.2],0.003,&mut a,&mut ya,&hessian).unwrap();
    dense.step_into_analytic(&sys,&x,&[0.2],0.003,&mut b,&mut yb,&hessian).unwrap();
    assert_eq!(reduced.linear_solve_counts(),(0,record.newton_iters));
    for (a,b) in a.iter().zip(b) {assert!((a-b).abs()<1e-11);}
    // Finite differences remain the original path even with a scalar plan selected.
    reduced.step_into(&sys,&x,&[0.2],0.003,&mut a,&mut ya).unwrap();
    dense.step_into(&sys,&x,&[0.2],0.003,&mut b,&mut yb).unwrap();assert_eq!(a,b);assert_eq!(ya,yb);
}

#[test]
fn scalar_configuration_refusals_and_switches_leave_the_old_plan_usable() {
    for (pairs,singles) in [(vec![],vec![1,1]),(vec![[1,2]],vec![1]),(vec![],vec![4]),(vec![[0,4]],vec![1])] {
        assert!(Condensation::with_singletons(4,&pairs,&singles).is_err());
    }
    assert!(Condensation::with_singletons(4,&[],&[]).is_err());
    let sys=system(0.0);let mut work=StepWorkspace::new(&sys).unwrap();
    work.set_condensed_leaves(&[],&[1,3]).unwrap();
    assert!(work.set_condensed_leaves(&[[0,1]],&[1]).is_err());assert_eq!(work.condensed_dimension(),Some(3));
    let x=[0.3,-0.6,-0.1,0.2];let mut a=[17.0;4];let mut y=[19.0;1];
    assert_eq!(work.step_into_analytic_controlled(&sys,&x,&[0.2],0.003,&mut a,&mut y,&hessian,||true).unwrap_err(),PreparedStepError::Cancelled);
    assert_eq!(a,[17.0;4]);assert_eq!(y,[19.0;1]);
    work.step_into_analytic(&sys,&x,&[0.2],0.003,&mut a,&mut y,&hessian).unwrap();
    work.set_condensed_pairs(&[[0,1]]).unwrap();assert_eq!(work.condensed_dimension(),Some(3));
    work.set_condensed_leaves(&[],&[]).unwrap();assert_eq!(work.condensed_dimension(),None);
    let mut dense=StepWorkspace::new(&sys).unwrap();let mut expected=[0.0;4];let mut ey=[0.0;1];
    dense.step_into_analytic(&sys,&x,&[0.2],0.003,&mut expected,&mut ey,&hessian).unwrap();
    work.step_into_analytic(&sys,&x,&[0.2],0.003,&mut a,&mut y,&hessian).unwrap();assert_eq!(a,expected);assert_eq!(y,ey);
}

#[test]
fn scalar_refinement_reuses_factors_and_checks_original_coordinates() {
    let (mut plan,base,rhs)=fixture(6,&[[5,1]],&[4,2]);let mut out=[0.0;6];
    assert!(plan.solve(&base,&rhs,&mut out,&mut ||Ok(())).unwrap());let factors=plan.matrix.clone();
    plan.candidate[4]+=1e-5;
    assert_eq!(plan.check(&base,&rhs,&mut ||Ok(())).unwrap(),Some(false));
    assert!(plan.solve_rhs(&base,&mut ||Ok(())).unwrap());
    for i in 0..6 {plan.candidate[i]+=plan.correction[i];}
    assert_eq!(plan.check(&base,&rhs,&mut ||Ok(())).unwrap(),Some(true));assert_eq!(plan.matrix,factors);
    for (a,b) in plan.candidate.iter().zip(out) {assert!((a-b).abs()<1e-13);}
}
