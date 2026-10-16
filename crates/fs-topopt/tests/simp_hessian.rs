//! Full raw-density Hessians, not a first-order filter around a stiffness HVP.
use std::ops::ControlFlow;
use fs_topopt::{DensityElasticity, DensityFilter, DesignPipeline, EvaluationStop, SimpParams,
    SolveBudget, SolveControl, SolveProgress};
use fs_topopt::elasticity::second_order::{ComplianceError, ComplianceSolveOptions};
use fs_topopt::pipeline::LoadCase;
use fs_topopt::simp_second_order::SimpComplianceGoal;

fn fixture(params: SimpParams) -> (DensityElasticity, DesignPipeline, Vec<f64>, Vec<f64>, Vec<f64>) {
    let (mesh, mut positions) = fs_feec::kuhn_cube(2);
    let tip: Vec<usize> = positions.iter().enumerate().filter(|(_,p)| p[0]==1.0).map(|(i,_)|i).collect();
    // Unequal cell volumes make substituting F for F^T detectable.
    for p in &mut positions {let [x,y,z]=*p;*p=[x*(1.0+0.25*z),y*(1.0+0.15*x),z*(1.0+0.1*x*y)];}
    let volumes = fs_feec::element_geometry(&mesh,&positions).vol_signed.iter().map(|v|v.abs()).collect();
    let el=DensityElasticity::new(&mesh,&positions,1.0,0.27,&|p|p[0]==0.0);
    let pipeline=DesignPipeline {filter:DensityFilter::new(&mesh,&positions,0.15),params};
    let mut a=vec![0.0;el.n()];let mut b=a.clone();
    for i in tip {a[3*i+2]=-0.2;b[3*i+1]=0.13*(1.0+positions[i][2]);}
    (el,pipeline,volumes,a,b)
}
fn options()->ComplianceSolveOptions {ComplianceSolveOptions {tolerance:1e-11,max_iterations:4000}}
fn close(a:&[f64],b:&[f64],tol:f64) {
    assert_eq!(a.len(),b.len());
    let err=a.iter().zip(b).map(|(a,b)|(a-b).abs()).fold(0.0,f64::max);
    let scale=b.iter().map(|b|b.abs()).fold(1.0,f64::max);
    assert!(err<tol*scale,"error {err:e}, scale {scale:e}, tolerance {tol:e}");
}
fn dot(a:&[f64],b:&[f64])->f64 {a.iter().zip(b).map(|(a,b)|a*b).sum()}

#[test]
fn full_simp_hessian_matches_legacy_gradient_differences_and_transpose_symmetry() {
    for (penal,beta) in [(1.0,0.0),(2.5,2.0),(3.0,6.0)] {
        let params=SimpParams {e_min:1e-3,penal,beta,eta:0.55};
        let (el,pipeline,_,a,b)=fixture(params);let n=el.cells();
        let (mut reference,_,_,_,_)=fixture(params);
        let rho:Vec<f64>=(0..n).map(|e|0.37+0.03*(e%7) as f64).collect();
        let v:Vec<f64>=(0..n).map(|e|0.02*((e%5) as f64-2.0)).collect();
        let w:Vec<f64>=(0..n).map(|e|0.03*((e%9) as f64-4.0)).collect();
        let loads=[LoadCase {force:&a,weight:0.4},LoadCase {force:&b,weight:1.2}];
        let mut callback=|_|ControlFlow::Continue(());let mut control=SolveControl::new(SolveBudget::default(),&mut callback);
        let goal=SimpComplianceGoal::new(&pipeline,&el,&rho,&loads,options(),&mut control).unwrap();
        assert_eq!(control.work().linear_solves,4,"one filter, two primals, one transpose");
        let old=pipeline.multi_load_compliance_and_gradient(&mut reference,&rho,&loads);
        close(&[goal.value()],&[old.compliance],1e-9);close(goal.gradient(),&old.gradient,2e-8);
        let before=control.work();let hv=goal.hessian_vector_product(&v,&mut control).unwrap();
        assert_eq!(control.work().linear_solves-before.linear_solves,4);
        let hw=goal.hessian_vector_product(&w,&mut control).unwrap();
        close(&[dot(&w,&hv)],&[dot(&v,&hw)],2e-8);
        let eps=5e-5;
        let plus:Vec<f64>=rho.iter().zip(&v).map(|(r,v)|r+eps*v).collect();
        let minus:Vec<f64>=rho.iter().zip(&v).map(|(r,v)|r-eps*v).collect();
        let gp=pipeline.multi_load_compliance_and_gradient(&mut reference,&plus,&loads);
        let gm=pipeline.multi_load_compliance_and_gradient(&mut reference,&minus,&loads);
        let fd:Vec<f64>=gp.gradient.iter().zip(&gm.gradient).map(|(a,b)|(a-b)/(2.0*eps)).collect();
        close(&hv,&fd,3e-6);
        close(&[(gp.compliance-gm.compliance)/(2.0*eps)],&[dot(goal.gradient(),&v)],3e-6);
        assert_eq!(el.moduli,vec![1.0;n],"the retained physical operator is immutable");
    }
}

#[test]
fn projected_volume_hessian_has_geometry_weights_and_no_simp_curvature() {
    let params=SimpParams {e_min:1e-3,penal:3.0,beta:4.0,eta:0.55};
    let (el,pipeline,volumes,a,_) = fixture(params);let n=el.cells();
    assert!(volumes.iter().any(|v|(v-volumes[0]).abs()>1e-5));
    let rho:Vec<f64>=(0..n).map(|e|0.37+0.03*(e%7) as f64).collect();
    let v:Vec<f64>=(0..n).map(|e|0.03+0.01*(e%5) as f64).collect();
    let mut callback=|_|ControlFlow::Continue(());let mut control=SolveControl::new(SolveBudget::default(),&mut callback);
    let goal=SimpComplianceGoal::new(&pipeline,&el,&rho,&[LoadCase {force:&a,weight:1.0}],options(),&mut control).unwrap();
    let (value,gradient)=goal.volume_value_gradient(&volumes,&mut control).unwrap();
    let (old_value,old_gradient)=pipeline.volume_and_gradient(&rho,&volumes);
    close(&[value],&[old_value],1e-12);close(&gradient,&old_gradient,1e-12);
    let before=control.work();let hv=goal.volume_hessian_vector_product(&volumes,&v,&mut control).unwrap();
    assert_eq!(control.work().linear_solves-before.linear_solves,2);
    assert!(hv.iter().any(|v|v.abs()>1e-7));
    let eps=1e-4;
    let plus:Vec<f64>=rho.iter().zip(&v).map(|(r,v)|r+eps*v).collect();
    let minus:Vec<f64>=rho.iter().zip(&v).map(|(r,v)|r-eps*v).collect();
    let (vp,gp)=pipeline.volume_and_gradient(&plus,&volumes);let (vm,gm)=pipeline.volume_and_gradient(&minus,&volumes);
    close(&hv,&gp.iter().zip(&gm).map(|(a,b)|(a-b)/(2.0*eps)).collect::<Vec<_>>(),1e-7);
    close(&[(vp-vm)/(2.0*eps)],&[dot(&gradient,&v)],1e-7);
}

#[test]
fn uniform_inverse_simp_formula_includes_material_curvature_and_linear_volume_is_zero() {
    let params=SimpParams {e_min:1e-3,penal:3.0,beta:0.0,eta:0.5};
    let (el,pipeline,volumes,a,_) = fixture(params);let n=el.cells();let rho=0.6_f64;
    let mut callback=|_|ControlFlow::Continue(());let mut control=SolveControl::new(SolveBudget::default(),&mut callback);
    let goal=SimpComplianceGoal::new(&pipeline,&el,&vec![rho;n],&[LoadCase {force:&a,weight:1.0}],options(),&mut control).unwrap();
    let modulus=params.e_min+(1.0-params.e_min)*rho*rho*rho;
    let first=3.0*(1.0-params.e_min)*rho*rho;let second=6.0*(1.0-params.e_min)*rho;
    close(&[goal.gradient().iter().sum()],&[-goal.value()*first/modulus],1e-9);
    let hv=goal.hessian_vector_product(&vec![1.0;n],&mut control).unwrap();
    let expected=goal.value()*(2.0*(first/modulus)*(first/modulus)-second/modulus);
    close(&[hv.iter().sum()],&[expected],1e-8);
    let before=control.work();
    assert_eq!(goal.volume_hessian_vector_product(&volumes,&vec![1.0;n],&mut control).unwrap(),vec![0.0;n]);
    assert_eq!(goal.hessian_vector_product(&vec![0.0;n],&mut control).unwrap(),vec![0.0;n]);
    assert_eq!(control.work(),before);
}

#[test]
fn refused_clamps_budget_and_cancellation_preserve_the_retained_smooth_goal() {
    let params=SimpParams {e_min:1e-3,penal:3.0,beta:2.0,eta:0.5};
    let (el,pipeline,_,a,_) = fixture(params);let n=el.cells();let loads=[LoadCase {force:&a,weight:1.0}];
    let mut callback=|_|ControlFlow::Continue(());let mut control=SolveControl::new(SolveBudget::default(),&mut callback);
    let goal=SimpComplianceGoal::new(&pipeline,&el,&vec![0.5;n],&loads,options(),&mut control).unwrap();
    let direction:Vec<f64>=(0..n).map(|e|0.1+0.02*(e%5) as f64).collect();
    let reference=goal.hessian_vector_product(&direction,&mut control).unwrap();
    assert!(matches!(SimpComplianceGoal::new(&pipeline,&el,&vec![0.0;n],&loads,options(),&mut control),
        Err(ComplianceError::InvalidInput(_))));
    let mut bounded_callback=|_|ControlFlow::Continue(());
    let mut bounded=SolveControl::new(SolveBudget {per_solve_iterations:4000,total_iterations:1},&mut bounded_callback);
    assert!(matches!(goal.hessian_vector_product(&direction,&mut bounded),
        Err(ComplianceError::Stopped(EvaluationStop::TotalBudget {..}))));
    assert_eq!(bounded.work().linear_iterations,1);
    let mut cancel=|p:SolveProgress| if p.stage=="compliance-tangent" {ControlFlow::Break(())} else {ControlFlow::Continue(())};
    let mut cancelled=SolveControl::new(SolveBudget::default(),&mut cancel);
    assert!(matches!(goal.hessian_vector_product(&direction,&mut cancelled),
        Err(ComplianceError::Stopped(EvaluationStop::Cancelled))));
    assert_eq!(cancelled.work().linear_solves,1);assert!(cancelled.work().linear_iterations>0);
    assert_eq!(goal.hessian_vector_product(&direction,&mut control).unwrap(),reference);
    let before=control.work();
    assert!(goal.hessian_vector_product(&[],&mut control).is_err());
    assert!(goal.hessian_vector_product(&vec![f64::NAN;n],&mut control).is_err());
    assert_eq!(control.work(),before);
}
