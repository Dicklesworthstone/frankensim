//! Native PDE second-order checks: no mock operator or differentiated solver.
mod support;
use fs_conduction::adjoint::ConductivityDesign;
use fs_conduction::{ConductionError, ConductionMesh, ConductionProblem, ConductivityModel,
    LinearConfig, ScalarField, ThermalBc, ThermalBoundary, ThermalBoundaryBuilder};
use fs_conduction::fixtures::{unit_cube,on_box_face};
use support::{with_cx,with_cancelled_cx};

struct Fixture { mesh: ConductionMesh,boundary: ThermalBoundary,material: ConductivityModel,source: ScalarField }
impl Fixture {
    fn new(robin: bool) -> Self {
        let (complex,positions)=unit_cube(2);
        let mesh=ConductionMesh::new(complex,positions).unwrap();
        let mut builder=ThermalBoundaryBuilder::new(&mesh).region("fixed",
            |f| on_box_face(f.centroid[2],0.0),ThermalBc::dirichlet(300.0).unwrap()).unwrap();
        if robin { builder=builder.region("robin",|f| on_box_face(f.centroid[0],1.0),
            ThermalBc::robin(8.0,290.0).unwrap()).unwrap(); }
        let boundary=builder.adiabatic_remainder().finish().unwrap();
        Self { mesh,boundary,material: ConductivityModel::constant_tensor(
            [[3.0,0.5,0.25],[0.5,2.0,0.3],[0.25,0.3,1.5]]).unwrap(),source: ScalarField::Uniform(4000.0) }
    }
    fn problem(&self) -> ConductionProblem<'_> {
        ConductionProblem { mesh: &self.mesh,boundary: &self.boundary,material: &self.material,
            source: &self.source,element_materials: None }
    }
    fn design(&self) -> ConductivityDesign<'_> {
        ConductivityDesign::new(self.problem(),LinearConfig { tolerance: 1e-12,max_iterations: 2000,restart: 30 }).unwrap()
    }
}
fn close(a: &[f64],b: &[f64],tol: f64) {
    assert_eq!(a.len(),b.len());
    let scale=b.iter().map(|x| x.abs()).fold(1.0_f64,f64::max);
    let error=a.iter().zip(b).map(|(a,b)| (a-b).abs()).fold(0.0_f64,f64::max);
    assert!(error<tol*scale,"max error {error}, reference scale {scale}");
}
fn dot(a: &[f64],b: &[f64])->f64 { a.iter().zip(b).map(|(a,b)| a*b).sum() }

#[test]
fn uniform_coefficient_scaling_has_exact_inverse_law_and_zero_boundary_tangent() {
    let f=Fixture::new(false);let design=f.design();let n=design.parameter_count();
    let rho:Vec<_>=(0..n).map(|i| 0.8+0.01*(i%17) as f64).collect();
    let weights=vec![1.0/design.dofs().n() as f64;design.dofs().n()];
    with_cx(|cx| {
        let goal=design.linearize_goal(cx,&rho,&weights).unwrap();
        let hv=goal.hessian_vector_product(cx,&rho).unwrap();
        close(&hv.values,&goal.gradient().iter().map(|g| -2.0*g).collect::<Vec<_>>(),2e-9);
        close(&hv.temperature_direction,&goal.solution().temperature.iter().map(|t| 300.0-t).collect::<Vec<_>>(),2e-9);
        assert!((dot(&rho,&hv.values)-2.0*(goal.value()-300.0)).abs()<1e-7);
        for &v in design.dofs().fixed() { assert_eq!(hv.temperature_direction[v],0.0); }
        assert!(hv.tangent.relative_residual<1e-12 && hv.incremental_adjoint.relative_residual<1e-12);
        assert!(goal.solution().primal_relative_residual<1e-12 && goal.dual_solve().relative_residual<1e-12);
    });
}

#[test]
fn anisotropic_robin_hessian_matches_independent_gradient_differences_and_symmetry() {
    let f=Fixture::new(true);let design=f.design();let n=design.parameter_count();
    let rho:Vec<_>=(0..n).map(|i| 0.7+0.07*(i%9) as f64).collect();
    let v:Vec<_>=(0..n).map(|i| ((i%7) as f64-3.0)/7.0).collect();
    let u:Vec<_>=(0..n).map(|i| ((i%11) as f64-5.0)/11.0).collect();
    let weights:Vec<_>=(0..design.dofs().n()).map(|i| if i%3==0 {-0.25} else {1.0}).collect();
    with_cx(|cx| {
        let goal=design.linearize_goal(cx,&rho,&weights).unwrap();
        let original=goal.gradient().to_vec();
        let hv=goal.hessian_vector_product(cx,&v).unwrap();let hu=goal.hessian_vector_product(cx,&u).unwrap();
        let h=2e-5;
        let plus:Vec<_>=rho.iter().zip(&v).map(|(r,v)| r+h*v).collect();
        let minus:Vec<_>=rho.iter().zip(&v).map(|(r,v)| r-h*v).collect();
        // Existing independent first-order IFT/GMRES route, not the new helper.
        let gp=design.gradient(cx,&plus,&weights).unwrap().0;
        let gm=design.gradient(cx,&minus,&weights).unwrap().0;
        let fd:Vec<_>=gp.iter().zip(&gm).map(|(p,m)| (p-m)/(2.0*h)).collect();
        close(&hv.values,&fd,2e-5);
        close(&[dot(&u,&hv.values)],&[dot(&v,&hu.values)],2e-8);
        let sum:Vec<_>=u.iter().zip(&v).map(|(a,b)| a+b).collect();
        let hs=goal.hessian_vector_product(cx,&sum).unwrap();
        close(&hs.values,&hu.values.iter().zip(&hv.values).map(|(a,b)| a+b).collect::<Vec<_>>(),2e-8);
        assert_eq!(goal.gradient(),original);
        assert_eq!(goal.hessian_vector_product(cx,&v).unwrap(),hv,"retained tape is repeatable");
    });
}

#[test]
fn zero_goals_directions_and_scaled_signed_weights_keep_physical_semantics() {
    let f=Fixture::new(true);let design=f.design();let rho=vec![1.0;design.parameter_count()];
    let weights=vec![1.0/design.dofs().n() as f64;design.dofs().n()];
    with_cx(|cx| {
        let goal=design.linearize_goal(cx,&rho,&weights).unwrap();
        let base=goal.hessian_vector_product(cx,&rho).unwrap();
        let zero=goal.hessian_vector_product(cx,&vec![0.0;rho.len()]).unwrap();
        assert!(zero.values.iter().chain(&zero.temperature_direction).all(|&x| x==0.0));
        assert_eq!(zero.tangent.iterations,0);assert_eq!(zero.incremental_adjoint.iterations,0);
        for scale in [-2.0,1e-120,1e120] {
            let weights:Vec<_>=weights.iter().map(|w| w*scale).collect();
            let other=design.linearize_goal(cx,&rho,&weights).unwrap();
            let hv=other.hessian_vector_product(cx,&rho).unwrap();
            close(&hv.values.iter().map(|v| v/scale).collect::<Vec<_>>(),&base.values,2e-8);
        }
        let empty=design.linearize_goal(cx,&rho,&vec![0.0;weights.len()]).unwrap();
        assert_eq!(empty.value(),0.0);assert_eq!(empty.dual_solve().iterations,0);
        assert!(empty.hessian_vector_product(cx,&rho).unwrap().values.iter().all(|&x| x==0.0));
    });
}

#[test]
fn invalid_inputs_and_exhausted_linear_solves_do_not_return_a_linearization() {
    let f=Fixture::new(true);let design=f.design();let rho=vec![1.0;design.parameter_count()];
    let weights=vec![1.0;design.dofs().n()];
    with_cx(|cx| {
        for bad in [vec![],vec![0.0;rho.len()],vec![f64::NAN;rho.len()]] {
            assert!(design.linearize_goal(cx,&bad,&weights).is_err());
        }
        for bad in [vec![],vec![f64::INFINITY;weights.len()]] {
            assert!(design.linearize_goal(cx,&rho,&bad).is_err());
        }
        let goal=design.linearize_goal(cx,&rho,&weights).unwrap();
        for bad in [vec![],vec![f64::NAN;rho.len()]] { assert!(goal.hessian_vector_product(cx,&bad).is_err()); }
        let limited=ConductivityDesign::new(f.problem(),LinearConfig { tolerance: 1e-12,max_iterations: 1,restart: 1 }).unwrap();
        assert!(matches!(limited.linearize_goal(cx,&rho,&weights),Err(ConductionError::LinearSolveFailed { .. })));
    });
}

#[test]
fn cancelled_preparation_and_action_preserve_the_immutable_linearization() {
    let f=Fixture::new(true);let design=f.design();let rho=vec![1.0;design.parameter_count()];
    let weights=vec![1.0;design.dofs().n()];
    with_cancelled_cx(|cx| assert!(matches!(design.linearize_goal(cx,&rho,&weights),Err(ConductionError::Cancelled { .. }))));
    with_cx(|cx| {
        let goal=design.linearize_goal(cx,&rho,&weights).unwrap();
        let before=goal.hessian_vector_product(cx,&rho).unwrap();
        with_cancelled_cx(|cx| assert!(matches!(goal.hessian_vector_product(cx,&rho),Err(ConductionError::Cancelled { .. }))));
        assert_eq!(goal.hessian_vector_product(cx,&rho).unwrap(),before);
    });
}
