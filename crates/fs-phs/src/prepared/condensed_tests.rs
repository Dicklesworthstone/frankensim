use super::*;
use crate::{PortHamiltonian, StepWorkspace, Storage};
use super::super::dissipation::Dissipation;

fn full(base:&[f64],u:&[f64],v:&[f64])->Vec<f64> {
    let n=u.len();(0..n*n).map(|a|base[a]+u[a/n]*v[a%n]).collect()
}
#[test]
fn bordered_rank_one_solve_matches_full_lu_with_permuted_noncontiguous_pairs() {
    for count in 1..=12 {
        let n=3+2*count; let pairs:Vec<_>=(0..count).map(|p|[n-1-p,1+p]).collect();
        let mut plan=Condensation::new(n,&pairs).unwrap();
        assert_eq!(plan.dimension(),4);
        for seed in 0..8 {
            let mut base=vec![0.0;n*n];
            for row in 0..n { for col in 0..n {
                if plan.owner[row]==usize::MAX || plan.owner[col]==usize::MAX
                    || plan.owner[row]==plan.owner[col] {
                    base[row*n+col]=if row==col {3.0} else {0.04*((row*7+col*11+seed)%13) as f64-0.24};
                }
            }}
            for i in 0..n { plan.left[i]=0.1*(i as f64-2.0);plan.right[i]=0.03*((i+seed)%7) as f64-0.07; }
            let matrix=full(&base,&plan.left,&plan.right);
            let rhs:Vec<_>=(0..n).map(|i|0.2*(i as f64+1.0)).collect();
            let (mut expected,mut actual)=(vec![0.0;n],vec![0.0;n]);
            LuWorkspace::new(n).unwrap().solve_into(&matrix,&rhs,&mut expected).unwrap();
            assert!(plan.solve(&base,&rhs,&mut actual,&mut ||Ok(())).unwrap());
            for (a,b) in actual.iter().zip(expected) {assert!((a-b).abs()<1e-11);}
        }
    }
    // Full uncorrected base is singular, but the corrected matrix is not.
    // Eliminating a scalar border does not need Sherman--Morrison's base inverse.
    let mut plan=Condensation::new(3,&[[1,2]]).unwrap();
    plan.left[0]=1.0;plan.right[0]=2.0;
    let mut answer=[0.0;3];
    assert!(plan.solve(&[0.,0.,0.,0.,1.,0.,0.,0.,1.],&[4.,3.,5.],&mut answer,&mut ||Ok(())).unwrap());
    assert_eq!(answer,[2.,3.,5.]);
}

#[test]
fn unsuited_structure_and_singular_leaves_never_publish_an_approximation() {
    let mut plan=Condensation::new(5,&[[1,2],[3,4]]).unwrap();
    let mut base=vec![0.0;25];for i in 0..5 {base[i*5+i]=1.0;}
    let mut out=[17.;5];
    // Even the smallest subnormal coefficient is a real coupling.
    base[1*5+3]=f64::from_bits(1);
    assert!(!plan.solve(&base,&[1.;5],&mut out,&mut ||Ok(())).unwrap());assert_eq!(out,[17.;5]);
    base[1*5+3]=0.;base[1*5+1]=0.;
    assert!(!plan.solve(&base,&[1.;5],&mut out,&mut ||Ok(())).unwrap());assert_eq!(out,[17.;5]);
    base[1*5+1]=1.;let mut polls=0;
    assert_eq!(plan.solve(&base,&[1.;5],&mut out,&mut ||{polls+=1;if polls==4 {Err(PreparedStepError::Cancelled)}else{Ok(())}}).unwrap_err(),PreparedStepError::Cancelled);
    assert_eq!(out,[17.;5]);
    assert!(plan.solve(&base,&[1.;5],&mut out,&mut ||Ok(())).unwrap());assert_eq!(out,[1.;5]);
    for pairs in [vec![[1,1]],vec![[1,5]],vec![[1,2],[2,3]]] {assert!(Condensation::new(5,&pairs).is_err());}
}

struct Potential(usize);
impl Storage for Potential {
    fn hamiltonian(&self,x:&[f64])->f64 {
        10.0*x[0].powi(4)+(0..self.0/2).map(|i|0.5*((4.0+i as f64)*x[2*i].powi(2)+x[2*i+1].powi(2))).sum::<f64>()
    }
    fn gradient(&self,x:&[f64],out:&mut[f64]) {
        for i in 0..self.0/2 {out[2*i]=(4.0+i as f64)*x[2*i];out[2*i+1]=x[2*i+1];}
        out[0]+=40.0*x[0].powi(3);
    }
}
impl Potential {
    fn hessian(&self,x:&[f64],d:&[f64],out:&mut[f64])->bool {
        for i in 0..self.0/2 {out[2*i]=(4.0+i as f64)*d[2*i];out[2*i+1]=d[2*i+1];}
        out[0]+=120.0*x[0]*x[0]*d[0];true
    }
}
fn system(pairs:usize)->PortHamiltonian {
    let n=4+2*pairs;let mut j=vec![0.;n*n];let mut r=vec![0.;n*n];let mut g=vec![0.;n*2];
    for i in 0..n/2 {let a=2*i;j[a*n+a+1]=1.;j[(a+1)*n+a]=-1.;r[(a+1)*n+a+1]=0.1;}
    for p in 0..pairs {let a=5+2*p;let l=if p%2==0 {0.7}else{-0.4};j[n+a]=-l;j[a*n+1]=l;}
    g[2]=1.;g[3*2+1]=1.;
    PortHamiltonian::new(n,2,j,r,g,Box::new(Potential(n))).unwrap()
}
fn resist(x:&[f64],e:&[f64],out:&mut[f64])->bool {out.fill(0.);out[1]=(0.2+x[0]*x[0])*e[1];true}
fn tangent(x:&[f64],e:&[f64],dx:&[f64],de:&[f64],out:&mut[f64])->bool {
    out.fill(0.);out[1]=2.*x[0]*dx[0]*e[1]+(0.2+x[0]*x[0])*de[1];true
}
fn plan(work:&mut StepWorkspace,n:usize) {work.set_condensed_pairs(&(4..n).step_by(2).map(|i|[i,i+1]).collect::<Vec<_>>()).unwrap();}

#[test]
fn split_jacobian_keeps_every_energy_and_nonlinear_dissipation_derivative() {
    let sys=system(4);let n=sys.n;let action=|x:&[f64],d:&[f64],o:&mut[f64]|Potential(n).hessian(x,d,o);
    for loss in [false,true] { for zero_increment in [false,true] {
        let mut dense=StepWorkspace::new(&sys).unwrap();let mut sparse=StepWorkspace::new(&sys).unwrap();plan(&mut sparse,n);
        dense.flow.refresh(&sys).unwrap();sparse.flow.refresh(&sys).unwrap();
        let mut x0=vec![0.0;n];x0[0]=0.3;x0[1]=-0.6;x0[2]=-0.1;
        let mut x1=x0.clone();if !zero_increment {for i in 0..n {x1[i]+=0.005*(i as f64-2.0);}}
        dense.x.copy_from_slice(&x1);sparse.x.copy_from_slice(&x1);
        let diss=if loss {Some(Dissipation{force:&resist,tangent:Some(&tangent)})}else{None};
        dense.analytic_jacobian_into(&sys,&x0,0.01,&action,diss,&mut ||Ok(())).unwrap();
        sparse.analytic_jacobian_into(&sys,&x0,0.01,&action,diss,&mut ||Ok(())).unwrap();
        let c=sparse.condensation.as_ref().unwrap();let restored=full(&sparse.jacobian,&c.left,&c.right);
        for (a,b) in restored.iter().zip(&dense.jacobian) {assert!((a-b).abs()<1e-13);}
        if !zero_increment {assert!(c.right.iter().any(|x|x.abs()>1e-7));}
    }}
}

#[test]
fn full_nonlinear_steps_eliminate_128_linear_states_and_keep_work_and_loss() {
    let sys=system(64);let n=sys.n;let mut dense=StepWorkspace::new(&sys).unwrap();
    let mut sparse=StepWorkspace::new(&sys).unwrap();plan(&mut sparse,n);
    assert_eq!(n,132);assert_eq!(sparse.condensed_dimension(),Some(5));
    let mut x=vec![0.0;n];x[0]=0.3;x[1]=-0.6;x[2]=0.1;let mut z=x.clone();
    let (mut a,mut b)=(vec![0.0;n],vec![0.0;n]);let (mut ya,mut yb)=([0.0;2],[0.0;2]);
    let initial=sys.hamiltonian(&x);let mut work=0.;let mut loss=0.;let mut count=0;
    let action=|x:&[f64],d:&[f64],o:&mut[f64]|Potential(n).hessian(x,d,o);
    for step in 0..60 {
        let force=[if step<20 {0.7}else{0.},-0.03];let dt=0.003;
        let fast=sparse.step_into_dissipative_controlled(&sys,&x,&force,dt,&mut a,&mut ya,
            &resist,Some(&action),Some(&tangent),||false).unwrap();
        dense.step_into_dissipative_controlled(&sys,&z,&force,dt,&mut b,&mut yb,
            &resist,Some(&action),Some(&tangent),||false).unwrap();
        for (v,w) in a.iter().zip(&b) {assert!((v-w).abs()<2e-9);}
        assert!(fast.balance_residual().abs()<1e-10);assert!(fast.dissipated>=0.);
        assert!((fast.supplied-dt*(force[0]*ya[0]+force[1]*ya[1])).abs()<1e-15);
        assert_eq!(sparse.linear_solve_counts(),(fast.newton_iters,0));count+=fast.newton_iters;
        work+=fast.supplied;loss+=fast.dissipated;x.copy_from_slice(&a);z.copy_from_slice(&b);
    }
    assert!(count>0);assert!((sys.hamiltonian(&x)-initial+loss-work).abs()<1e-9);
}

#[test]
fn cancelled_condensed_step_retries_exactly_and_disabling_preserves_dense_path() {
    let sys=system(3);let n=sys.n;let mut work=StepWorkspace::new(&sys).unwrap();plan(&mut work,n);
    let mut x=vec![0.;n];x[0]=0.3;x[1]=-0.6;
    let mut out=vec![123.;n];let mut y=[456.;2];let mut polls=0;
    let action=|x:&[f64],d:&[f64],o:&mut[f64]|Potential(n).hessian(x,d,o);
    let error=work.step_into_analytic_controlled(&sys,&x,&[0.2,0.],0.003,&mut out,&mut y,&action,
        ||{polls+=1;polls==2*n+5}).unwrap_err();
    assert_eq!(error,PreparedStepError::Cancelled);assert_eq!(out,vec![123.;n]);assert_eq!(y,[456.;2]);
    let mut fresh=StepWorkspace::new(&sys).unwrap();plan(&mut fresh,n);
    let mut expected=vec![0.;n];let mut ye=[0.;2];
    work.step_into_analytic(&sys,&x,&[0.2,0.],0.003,&mut out,&mut y,&action).unwrap();
    fresh.step_into_analytic(&sys,&x,&[0.2,0.],0.003,&mut expected,&mut ye,&action).unwrap();assert_eq!(out,expected);assert_eq!(y,ye);
    assert!(work.set_condensed_pairs(&[[0,n]]).is_err());assert_eq!(work.condensed_dimension(),Some(5));
    work.set_condensed_pairs(&[]).unwrap();let mut dense=StepWorkspace::new(&sys).unwrap();
    work.step_into_analytic(&sys,&x,&[0.2,0.],0.003,&mut out,&mut y,&action).unwrap();
    dense.step_into_analytic(&sys,&x,&[0.2,0.],0.003,&mut expected,&mut ye,&action).unwrap();assert_eq!(out,expected);assert_eq!(y,ye);
}

#[test]
fn dense_fallback_retains_new_interpair_coupling_and_finite_difference_calls() {
    let mut sys=system(3);let n=sys.n;let mut work=StepWorkspace::new(&sys).unwrap();plan(&mut work,n);
    // Change operator topology after workspace construction. This must NOT be
    // discarded based on a previously valid pair plan.
    sys.j[5*n+7]=0.3;sys.j[7*n+5]=-0.3;
    let mut x=vec![0.;n];x[0]=0.3;x[1]=-0.6;
    let (mut a,mut b)=(vec![0.;n],vec![0.;n]);let (mut ya,mut yb)=([0.;2],[0.;2]);
    let mut dense=StepWorkspace::new(&sys).unwrap();
    let action=|x:&[f64],d:&[f64],o:&mut[f64]|Potential(n).hessian(x,d,o);
    let f=work.step_into_analytic(&sys,&x,&[0.2,0.],0.003,&mut a,&mut ya,&action).unwrap();
    dense.step_into_analytic(&sys,&x,&[0.2,0.],0.003,&mut b,&mut yb,&action).unwrap();
    assert_eq!(work.linear_solve_counts(),(0,f.newton_iters));for (a,b) in a.iter().zip(&b){assert!((a-b).abs()<1e-11);}
    work.step_into(&sys,&x,&[0.2,0.],0.003,&mut a,&mut ya).unwrap();
    dense.step_into(&sys,&x,&[0.2,0.],0.003,&mut b,&mut yb).unwrap();assert_eq!(a,b);assert_eq!(ya,yb);
}

#[test]
fn auxiliary_unit_scaling_does_not_change_the_original_equation_or_its_solution() {
    let n=9;let pairs=[[8,1],[7,2],[6,3]];
    let mut plan=Condensation::new(n,&pairs).unwrap();
    let mut base=vec![0.0;n*n];
    for i in 0..n {for j in 0..n {
        if plan.owner[i]==usize::MAX || plan.owner[j]==usize::MAX || plan.owner[i]==plan.owner[j] {
            base[i*n+j]=if i==j {2.0+i as f64} else {0.01*((i*5+j*3)%7) as f64-0.03};
        }
    }}
    let u:Vec<_>=(0..n).map(|i|0.1*(i as f64-3.0)).collect();
    let v:Vec<_>=(0..n).map(|i|0.07*(i as f64-4.0)).collect();
    let matrix=full(&base,&u,&v);
    let want:Vec<_>=(0..n).map(|i|0.2*(i as f64-2.0)).collect();
    let rhs:Vec<_>=(0..n).map(|i|(0..n).map(|j|matrix[i*n+j]*want[j]).sum()).collect();
    for scale in [1e-140,1e-70,1.0,1e70,1e140] {
        for i in 0..n {plan.left[i]=u[i]*scale;plan.right[i]=v[i]/scale;}
        let mut out=vec![91.0;n];
        assert!(plan.solve(&base,&rhs,&mut out,&mut ||Ok(())).unwrap(),"auxiliary scale {scale}");
        for (x,y) in out.iter().zip(&want) {assert!((x-y).abs()<2e-12);}
        assert!(plan.refinements<=MAX_REFINEMENTS);
    }
}

#[test]
fn refinement_solves_the_original_residual_with_the_same_factored_border() {
    let mut plan=Condensation::new(5,&[[1,2],[3,4]]).unwrap();
    let base=[3.0,0.2,-0.1,0.3,0.4, 0.1,2.0,0.2,0.0,0.0,
        -0.2,-0.1,4.0,0.0,0.0, 0.4,0.0,0.0,2.0,0.3, -0.1,0.0,0.0,-0.2,3.0];
    plan.left.copy_from_slice(&[1e-100,-2e-100,3e-100,0.0,4e-100]);
    plan.right.copy_from_slice(&[1e99,2e99,-1e99,3e99,-2e99]);
    let rhs=[0.2,0.1,-0.3,0.4,0.8];let mut answer=[0.0;5];
    assert!(plan.solve(&base,&rhs,&mut answer,&mut ||Ok(())).unwrap());
    let factored=plan.matrix.clone();let responses=plan.responses.clone();
    // A known inaccurate candidate must fail the unchanged componentwise gate.
    // Then solve its ORIGINAL-space residual, not the scaled Schur residual.
    plan.candidate[2]+=1e-6;
    assert_eq!(plan.check(&base,&rhs,&mut ||Ok(())).unwrap(),Some(false));
    assert!(plan.solve_rhs(&base,&mut ||Ok(())).unwrap());
    for i in 0..5 {plan.candidate[i]+=plan.correction[i];}
    assert_eq!(plan.check(&base,&rhs,&mut ||Ok(())).unwrap(),Some(true));
    assert_eq!(plan.matrix,factored);
    for row in 0..4 {for j in 0..plan.dimension() {
        assert_eq!(plan.responses[row*(plan.dimension()+1)+j],responses[row*(plan.dimension()+1)+j]);
    }}
    for (x,y) in plan.candidate.iter().zip(answer) {assert!((x-y).abs()<1e-14);}
}


#[test]
fn equilibration_that_would_underflow_a_nonzero_border_entry_uses_dense_fallback() {
    let mut plan=Condensation::new(4,&[[2,3]]).unwrap();
    let mut base=[0.0;16];for i in 0..4 {base[4*i+i]=1.0;}
    base[0]=1e300;base[1]=1e-100;
    let mut out=[17.0;4];
    assert!(!plan.solve(&base,&[1.0;4],&mut out,&mut ||Ok(())).unwrap());
    assert_eq!(out,[17.0;4]);
    // Restore resolvable units: a refused preparation did not poison reuse.
    base[0]=1.0;
    assert!(plan.solve(&base,&[1.0;4],&mut out,&mut ||Ok(())).unwrap());
    assert_eq!(out,[1.0;4]);
}

#[test]
fn scalar_border_scaling_preserves_extreme_rank_one_factors_and_zero_updates() {
    let n=5; let pairs=[[1,2],[3,4]];
    let mut base=vec![0.0;n*n]; for i in 0..n { base[i*n+i]=2.0+i as f64; }
    base[1]=-0.25; base[n]=0.5; base[3]=0.125; base[3*n]=-0.5;
    let u=[0.5,-0.25,0.125,0.25,-0.5]; let v=[0.25,0.5,-0.25,0.125,0.25];
    let rhs=[1.0,-2.0,3.0,-4.0,5.0]; let mut expected=[0.0;5];
    LuWorkspace::new(n).unwrap().solve_into(&full(&base,&u,&v),&rhs,&mut expected).unwrap();
    for shift in [-1000,-500,0,500,1000] {
        let mut plan=Condensation::new(n,&pairs).unwrap();
        for i in 0..n { plan.left[i]=scale_binary(u[i],shift); plan.right[i]=scale_binary(v[i],-shift); }
        let left=plan.left.clone(); let right=plan.right.clone(); let mut out=[99.0;5];
        assert!(plan.solve(&base,&rhs,&mut out,&mut ||Ok(())).unwrap());
        for (a,b) in out.iter().zip(expected) { assert!((a-b).abs()<2e-13); }
        assert_eq!(plan.left,left); assert_eq!(plan.right,right);
        for i in 0..n { for j in 0..n {
            assert_eq!(plan.scaled_left[i]*plan.scaled_right[j],left[i]*right[j]);
        }}
    }
    // An identically zero update must not require v^T*x to be representable.
    let mut plan=Condensation::new(n,&pairs).unwrap(); plan.right.fill(f64::MAX);
    let mut out=[0.0;5]; let huge_rhs=[1e100;5];
    assert!(plan.solve(&base,&huge_rhs,&mut out,&mut ||Ok(())).unwrap());
    LuWorkspace::new(n).unwrap().solve_into(&base,&huge_rhs,&mut expected).unwrap();
    for (a,b) in out.iter().zip(expected) { assert!((a-b).abs()<1e-14*b.abs()); }
}

#[test]
fn binary_balance_never_erases_nonzero_subnormal_couplings() {
    let mut plan=Condensation::new(3,&[[1,2]]).unwrap();
    plan.left=vec![1.0,f64::from_bits(1),-0.0]; plan.right=vec![1e-200,0.0,-1e-200];
    plan.balance();
    // Normalizing the large entries would lose left[1]; keep the original
    // factors rather than silently sparsifying the auxiliary equation.
    assert_eq!(plan.scaled_left,plan.left); assert_eq!(plan.scaled_right,plan.right);
    assert_eq!(plan.scaled_left[1].to_bits(),1);
    assert_eq!(plan.scaled_left[2].to_bits(),(-0.0_f64).to_bits());
    assert_eq!(exponent(f64::from_bits(1)),-1074);
    assert_eq!(exponent(f64::MIN_POSITIVE),-1022);
    assert_eq!(exponent(f64::MAX),1023);
    for x in [f64::from_bits(1),0.125,1.0,f64::MAX] {
        assert_eq!(scale_binary(x,0).to_bits(),x.to_bits());
    }
    // Exercise binary scaling beyond the exponent of a single finite factor.
    let mut extremes=Condensation::new(3,&[[1,2]]).unwrap();
    extremes.left[0]=f64::from_bits(1);extremes.right[0]=f64::MAX;extremes.balance();
    assert!(extremes.scaled_left[0].is_normal() && extremes.scaled_right[0].is_normal());
    assert_eq!(extremes.left[0]*extremes.right[0],extremes.scaled_left[0]*extremes.scaled_right[0]);
    let mut out=[17.0;3]; let mut calls=0;
    let base=[1.,0.,0.,0.,1.,0.,0.,0.,1.];
    assert_eq!(plan.solve(&base,&[0.25;3],&mut out,&mut || {
        calls+=1; if calls==3 {Err(PreparedStepError::Cancelled)} else {Ok(())}
    }).unwrap_err(),PreparedStepError::Cancelled);
    assert_eq!(out,[17.0;3]);
    assert!(plan.solve(&base,&[0.25;3],&mut out,&mut ||Ok(())).unwrap());
}


#[test]
fn residual_correction_recovers_a_well_conditioned_system_with_a_tiny_leaf_pivot() {
    let epsilon=1e-20_f64;
    let base=[1.,1.,0.,0.,0., 1.,epsilon,0.,0.,0., 0.,0.,1.,0.,0.,
        0.,0.,0.,2.,0., 0.,0.,0.,0.,3.];
    let rhs=[1.,2.,3.,4.,5.];
    // A is nonsingular and well-conditioned. Eliminating the tiny (not zero)
    // leaf pivot loses its recovered coordinate through subtractive cancellation.
    let naive_border=(1.-2./epsilon)/(1.-1./epsilon);
    let naive_leaf=2./epsilon-naive_border/epsilon;
    assert!((naive_border+naive_leaf-rhs[0]).abs()>0.5);
    let mut plan=Condensation::new(5,&[[1,2],[3,4]]).unwrap();
    // Locate the first residual-correction solve after the current factor
    // and equilibration traversal, rather than coupling to an old poll count.
    let mut before_correction=0;
    assert!(plan.factor(&base,&mut ||{before_correction+=1;Ok(())}).unwrap());
    plan.error.copy_from_slice(&rhs);
    assert!(plan.solve_rhs(&base,&mut ||{before_correction+=1;Ok(())}).unwrap());
    plan.candidate.copy_from_slice(&plan.correction);
    assert_eq!(plan.check(&base,&rhs,&mut ||{before_correction+=1;Ok(())}).unwrap(),Some(false));
    let mut answer=[17.;5];let mut polls=0;
    // Cancellation during correction must not publish the inaccurate candidate.
    assert_eq!(plan.solve(&base,&rhs,&mut answer,&mut || {
        polls+=1;if polls==before_correction+1 {Err(PreparedStepError::Cancelled)} else {Ok(())}
    }).unwrap_err(),PreparedStepError::Cancelled);
    assert_eq!(answer,[17.;5]);
    assert!(plan.solve(&base,&rhs,&mut answer,&mut ||Ok(())).unwrap());
    let mut expected=[0.;5];LuWorkspace::new(5).unwrap().solve_into(&base,&rhs,&mut expected).unwrap();
    for (a,b) in answer.iter().zip(expected) {assert!((a-b).abs()<1e-13);}
    let mut fresh=Condensation::new(5,&[[1,2],[3,4]]).unwrap();let mut retry=[0.;5];
    assert!(fresh.solve(&base,&rhs,&mut retry,&mut ||Ok(())).unwrap());assert_eq!(answer,retry);
}
