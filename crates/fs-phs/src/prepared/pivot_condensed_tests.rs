//! Non-principal pivots: equation and unknown addresses are NOT interchangeable.
use super::*;
use crate::{PortHamiltonian, StepWorkspace, Storage};

fn full(plan:&Condensation,base:&[f64])->Vec<f64> {
    let n=plan.n;(0..n*n).map(|i|base[i]+plan.left[i/n]*plan.right[i%n]).collect()
}
fn fixture(modes:usize,air:usize,history:usize)->(Condensation,Vec<f64>,Vec<f64>) {
    let n=2*modes+2*air+history;
    let permute=|i:usize|(i+n/3)%n;
    let pivots:Vec<_>=(0..modes).map(|i|[permute(2*i),permute(2*i+1)]).collect();
    // History can put the first acoustic state at an odd absolute address.
    let pairs:Vec<_>=(0..air).map(|i|[permute(2*modes+history+2*i),permute(2*modes+history+2*i+1)]).collect();
    let mut plan=Condensation::with_pivots(n,&pivots,&pairs).unwrap();
    let mut base=vec![0.0;n*n];
    for (i,&row) in plan.retained_rows.iter().enumerate() {
        for (j,&col) in plan.retained.iter().enumerate() {
            base[row*n+col]=if i==j {3.0} else {0.004*((i*7+j*3)%11) as f64};
        }
    }
    for (p,leaf) in plan.leaves.iter().enumerate() {
        for (i,&row) in leaf.equations().iter().enumerate() {
            for (j,&col) in leaf.coordinates().iter().enumerate() {base[row*n+col]=if i==j {2.0}else{0.1};}
            for prior in &plan.leaves[..p] {for &col in prior.coordinates() {
                base[row*n+col]=0.003*((row+col*7)%13) as f64;
            }}
            for &col in &plan.retained {base[row*n+col]=0.006*((row*3+col)%7) as f64;}
        }
        for &row in &plan.retained_rows {for &col in leaf.coordinates() {
            base[row*n+col]=-0.008*((row+col*5)%17) as f64;
        }}
    }
    for i in 0..n {plan.left[i]=0.01*(i%5) as f64;plan.right[i]=-0.02*(i%7) as f64;}
    let rhs=(0..n).map(|i|0.1+i as f64*0.03).collect();
    (plan,base,rhs)
}

#[test]
fn ordered_nonprincipal_pivots_and_reciprocal_pairs_match_the_full_equation() {
    for modes in [1,2,4,8,16] {for air in [0,2,5] {for history in [0,3] {
        let (mut plan,base,rhs)=fixture(modes,air,history);let n=plan.n;
        assert_eq!(plan.dimension(),modes+history+1);
        assert_ne!(plan.retained,plan.retained_rows);
        let mut expected=vec![0.0;n];let mut actual=vec![0.0;n];
        LuWorkspace::new(n).unwrap().solve_into(&full(&plan,&base),&rhs,&mut expected).unwrap();
        assert!(plan.solve(&base,&rhs,&mut actual,&mut ||Ok(())).unwrap());
        for (a,b) in actual.iter().zip(expected) {assert!((a-b).abs()<2e-11);}
    }}}
}

#[test]
fn even_subnormal_upper_edges_and_singular_pivots_require_full_lu() {
    let (mut plan,mut base,rhs)=fixture(3,2,1);let n=plan.n;let mut out=vec![17.0;n];
    let row=plan.leaves[0].rows[0];let col=plan.leaves[1].indices[0];
    base[row*n+col]=f64::from_bits(1);
    assert!(!plan.solve(&base,&rhs,&mut out,&mut ||Ok(())).unwrap());assert_eq!(out,vec![17.0;n]);
    base[row*n+col]=0.0;
    assert!(plan.solve(&base,&rhs,&mut out,&mut ||Ok(())).unwrap());
    // A zero off-diagonal kinematic pivot is not evidence that A is singular.
    let mut plan=Condensation::with_pivots(2,&[[0,1]],&[]).unwrap();
    let mut out=[19.0;2];let identity=[1.0,0.0,0.0,1.0];
    assert!(!plan.solve(&identity,&[2.0,3.0],&mut out,&mut ||Ok(())).unwrap());assert_eq!(out,[19.0;2]);
    LuWorkspace::new(2).unwrap().solve_into(&identity,&[2.0,3.0],&mut out).unwrap();assert_eq!(out,[2.0,3.0]);
}

#[test]
fn singular_base_and_complete_row_column_permutation_keep_the_gonzalez_border() {
    let mut plan=Condensation::with_pivots(2,&[[0,1]],&[]).unwrap();
    plan.left[1]=1.0;plan.right[0]=1.0;
    let mut out=[0.0;2];
    assert!(plan.solve(&[0.0,1.0,0.0,0.0],&[3.0,2.0],&mut out,&mut ||Ok(())).unwrap());
    assert_eq!(out,[2.0,3.0]);
    let mut plan=Condensation::with_pivots(3,&[[0,2],[1,0],[2,1]],&[]).unwrap();
    assert_eq!(plan.dimension(),1);
    let base=[0.0,0.0,2.0,3.0,0.0,0.4,0.5,4.0,0.6];
    let mut out=[0.0;3];
    assert!(plan.solve(&base,&[6.0,4.2,10.3],&mut out,&mut ||Ok(())).unwrap());
    for (a,b) in out.iter().zip([1.0,2.0,3.0]) {assert!((a-b).abs()<1e-13);}
}

#[test]
fn cancellation_never_publishes_and_retry_reuses_all_triangular_scratch() {
    let (mut plan,base,rhs)=fixture(7,5,3);let n=plan.n;
    let mut expected=vec![0.0;n];let mut polls=0;
    assert!(plan.solve(&base,&rhs,&mut expected,&mut ||{polls+=1;Ok(())}).unwrap());
    for stop in [1,5,polls/2,polls-1,polls] {
        let mut count=0;let mut out=vec![23.0;n];
        let error=plan.solve(&base,&rhs,&mut out,&mut ||{
            count+=1;if count==stop {Err(PreparedStepError::Cancelled)}else{Ok(())}
        }).unwrap_err();
        assert_eq!(error,PreparedStepError::Cancelled);assert_eq!(out,vec![23.0;n]);
        assert!(plan.solve(&base,&rhs,&mut out,&mut ||Ok(())).unwrap());assert_eq!(out,expected);
    }
}

struct Quartic;
impl Storage for Quartic {
    fn hamiltonian(&self,x:&[f64])->f64 {
        0.5*(x[1]*x[1]+x[3]*x[3])+2.0*x[0]*x[0]+3.0*x[2]*x[2]+(x[0]-x[2]).powi(4)
    }
    fn gradient(&self,x:&[f64],out:&mut[f64]) {
        let t=4.0*(x[0]-x[2]).powi(3);
        out.copy_from_slice(&[4.0*x[0]+t,x[1],6.0*x[2]-t,x[3]]);
    }
}
fn hessian(x:&[f64],d:&[f64],out:&mut[f64])->bool {
    let t=12.0*(x[0]-x[2]).powi(2)*(d[0]-d[2]);
    out.copy_from_slice(&[4.0*d[0]+t,d[1],6.0*d[2]-t,d[3]]);true
}
fn coupled()->PortHamiltonian {
    let mut j=vec![0.0;16];let mut r=vec![0.0;16];let mut g=vec![0.0;4];
    for i in [0,2] {j[4*i+i+1]=1.0;j[4*(i+1)+i]=-1.0;}
    // One physical relative-velocity damper, not diagonal modal attenuation.
    r[5]=0.8;r[7]=-0.8;r[13]=-0.8;r[15]=0.8;g[1]=1.0;
    PortHamiltonian::new(4,1,j,r,g,Box::new(Quartic)).unwrap()
}

#[test]
fn cross_damped_nonlinear_steps_preserve_forced_energy_balance_and_dense_trajectory() {
    let sys=coupled();let mut fast=StepWorkspace::new(&sys).unwrap();let mut dense=StepWorkspace::new(&sys).unwrap();
    fast.set_condensed_pivots(&[[0,1],[2,3]],&[]).unwrap();assert_eq!(fast.condensed_dimension(),Some(3));
    let mut x=[0.3,-0.6,-0.1,0.2];let mut reference=x;
    let initial=sys.hamiltonian(&x);let mut net=0.0;let mut calls=0;
    for tick in 0..400 {
        let force=[if tick<100 {0.2}else{0.0}];let (mut a,mut b)=([0.0;4],[0.0;4]);let (mut ya,mut yb)=([0.0;1],[0.0;1]);
        let report=fast.step_into_analytic(&sys,&x,&force,0.003,&mut a,&mut ya,&hessian).unwrap();
        dense.step_into_analytic(&sys,&reference,&force,0.003,&mut b,&mut yb,&hessian).unwrap();
        for (a,b) in a.iter().zip(b) {assert!((a-b).abs()<1e-9);}
        assert!(report.dissipated>=0.0);net+=report.supplied-report.dissipated;
        assert_eq!(fast.linear_solve_counts().1,0);calls+=fast.linear_solve_counts().0;x=a;reference=b;
    }
    assert!(calls>0);assert!((sys.hamiltonian(&x)-initial-net).abs()<1e-9);
    // The default finite-difference image still ignores the numerical plan.
    let (mut a,mut b)=([0.0;4],[0.0;4]);let (mut ya,mut yb)=([0.0;1],[0.0;1]);
    fast.step_into(&sys,&x,&[0.1],0.003,&mut a,&mut ya).unwrap();
    dense.step_into(&sys,&x,&[0.1],0.003,&mut b,&mut yb).unwrap();assert_eq!(a,b);assert_eq!(ya,yb);
}

#[test]
fn invalid_equation_or_unknown_maps_preserve_the_previous_prepared_plan() {
    let sys=coupled();let mut workspace=StepWorkspace::new(&sys).unwrap();
    workspace.set_condensed_pivots(&[[0,1],[2,3]],&[]).unwrap();
    for (pivots,pairs) in [
        (vec![[0,1],[0,3]],vec![]),(vec![[0,1],[2,1]],vec![]),
        (vec![[4,1]],vec![]),(vec![[0,4]],vec![]),
        (vec![[0,1]],vec![[0,2]]),(vec![[0,1]],vec![[1,2]]),
    ] {
        assert!(workspace.set_condensed_pivots(&pivots,&pairs).is_err());
        assert_eq!(workspace.condensed_dimension(),Some(3));
    }
    let x=[0.2,0.1,0.0,-0.1];let mut out=[13.0;4];let mut y=[17.0;1];
    assert_eq!(workspace.step_into_analytic_controlled(&sys,&x,&[0.0],0.001,&mut out,&mut y,&hessian,||true)
        .unwrap_err(),PreparedStepError::Cancelled);
    assert_eq!(out,[13.0;4]);assert_eq!(y,[17.0;1]);
    workspace.step_into_analytic(&sys,&x,&[0.0],0.001,&mut out,&mut y,&hessian).unwrap();
    workspace.set_condensed_pivots(&[],&[]).unwrap();assert_eq!(workspace.condensed_dimension(),None);
}
