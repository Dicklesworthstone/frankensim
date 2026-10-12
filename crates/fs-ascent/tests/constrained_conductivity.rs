#![cfg(feature = "conduction-design")]
//! Actual thermal constraints, independently known uniform-coefficient optimum.
use fs_ascent::conduction_design::ConductivityNewtonModel;
use fs_ascent::reverse_newton::{AugmentedNewtonStudy, AugmentedNewtonOptions, AugmentedNewtonBudget, AugmentedNewtonStop};
use fs_conduction::adjoint::ConductivityDesign;
use fs_conduction::{ConductionMesh, ConductionProblem, ConductivityModel, LinearConfig, ScalarField,
    ThermalBc, ThermalBoundaryBuilder, TEMPERATURE_DIMS};
use fs_conduction::fixtures::{unit_cube, on_box_face};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_opt::reverse::physics::{PhysicsBinding, PhysicsError, PhysicsHessian, PhysicsModel, PhysicsSample, PhysicsSignature};
use fs_opt::reverse::ReverseLimits;
use fs_opt::{ConstraintKind, Manifold, ProblemBuilder, ReverseProblem, Sense};
use fs_qty::Dims;

// A single uniform log multiplier expands onto ALL actual element parameters.
// This reduction is only a test fixture. Both first- and second-order products
// still execute the native P1 state/adjoint/sensitivity solves.
#[derive(Debug)]
struct Uniform<'a, 'd, 'm> { model: &'a ConductivityNewtonModel<'d, 'm>, n: usize }
impl PhysicsModel for Uniform<'_, '_, '_> {
    fn signature(&self) -> PhysicsSignature<'_> {
        PhysicsSignature { study: "uniform", parameter_count: 1,
            parameter_dims: Dims::NONE, residual_dims: TEMPERATURE_DIMS }
    }
    fn value_gradient(&self, x: &[f64], cx: Option<&Cx<'_>>) -> Result<PhysicsSample, PhysicsError> {
        if x.len()!=1 { return Err(PhysicsError::InvalidInput("uniform point length")); }
        let sample=self.model.value_gradient(&vec![x[0]; self.n], cx)?;
        Ok(PhysicsSample { value: sample.value, gradient: vec![sample.gradient.iter().sum()] })
    }
    fn second_order(&self) -> Option<&dyn PhysicsHessian> { Some(self) }
}
impl PhysicsHessian for Uniform<'_, '_, '_> {
    fn hessian_vector_product(&self, x: &[f64], v: &[f64], cx: Option<&Cx<'_>>) -> Result<Vec<f64>, PhysicsError> {
        if x.len()!=1 || v.len()!=1 { return Err(PhysicsError::InvalidInput("uniform direction length")); }
        Ok(vec![self.model.hessian_vector_product(&vec![x[0];self.n], &vec![v[0];self.n], cx)?.iter().sum()])
    }
}
fn with_design(test: impl FnOnce(&ConductivityDesign<'_>, &Cx<'_>)) {
    let (complex, positions)=unit_cube(2);let mesh=ConductionMesh::new(complex,positions).unwrap();
    let boundary=ThermalBoundaryBuilder::new(&mesh).region("base",|f|on_box_face(f.centroid[2],0.0),
        ThermalBc::dirichlet(300.0).unwrap()).unwrap().adiabatic_remainder().finish().unwrap();
    let material=ConductivityModel::constant_tensor([[3.0,0.5,0.25],[0.5,2.0,0.3],[0.25,0.3,1.5]]).unwrap();
    let source=ScalarField::Uniform(4000.0);
    let design=ConductivityDesign::new(ConductionProblem {mesh:&mesh,boundary:&boundary,material:&material,
        source:&source,element_materials:None},LinearConfig {tolerance:1e-12,max_iterations:2000,restart:30}).unwrap();
    let gate=CancelGate::new();let pool=fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx=Cx::new(&gate,arena,StreamKey {seed:1,kernel_id:81,tile:0,iteration:0},Budget::INFINITE,ExecMode::Deterministic);
        test(&design,&cx);
    });
}
fn work() -> AugmentedNewtonBudget { AugmentedNewtonBudget {objective_evals:300,hessian_products:1000,multiplier_updates:50} }

#[test]
fn conductivity_allowance_constrains_the_native_temperature_objective() {
    with_design(|design,cx| {
        let weights=vec![1.0/design.dofs().n() as f64;design.dofs().n()];
        let native=ConductivityNewtonModel::new(design,"native",weights.clone(),10000).unwrap();
        let model=Uniform {model:&native,n:design.parameter_count()};
        let initial=model.value_gradient(&[0.0],Some(cx)).unwrap().value;
        let mut b=ProblemBuilder::new();let v=b.var("uniform log rho",Manifold::Rn {dim:1},Dims::NONE).unwrap();
        let r=b.var_ref(v).unwrap();let x=b.component(r,0).unwrap();let rho=b.exp(x).unwrap();
        let node=b.pde_residual("uniform",v,true,TEMPERATURE_DIMS).unwrap();
        let base=b.konst(300.0,TEMPERATURE_DIMS).unwrap();let scale=b.konst(1000.0,TEMPERATURE_DIMS).unwrap();
        let rise=b.sub(node,base).unwrap();let f=b.div(rise,scale).unwrap();b.objective(f,Sense::Minimize,1.0).unwrap();
        let cap=b.konst(1.5,Dims::NONE).unwrap();let residual=b.sub(rho,cap).unwrap();
        b.constraint(residual,ConstraintKind::LeZero,"conductivity allowance").unwrap();let p=b.finish();
        let oracle=ReverseProblem::new_with_physics(&p,ReverseLimits {max_nodes:128,max_scalar_slots:4096},
            &[PhysicsBinding {node,model:&model}],None).unwrap();
        let mut study=AugmentedNewtonStudy::new(&oracle,&[0.0],AugmentedNewtonOptions::default(),Some(cx)).unwrap();
        let mut split=study.clone();let report=study.run(1e-5,150,work(),Some(cx)).unwrap();
        assert_eq!(report.stop,AugmentedNewtonStop::Converged,"{report:?}");assert!(report.kkt.within_tolerance(1e-5));
        let result_rho=fs_math::det::exp(report.x[0]);assert!((result_rho-1.5).abs()<2e-5);
        assert!(report.constraints[0]<1e-5 && report.multipliers[0]>0.0);
        let expected=(initial-300.0)/(1000.0*1.5);
        assert!((report.objective-expected).abs()<1e-5);
        let fresh=design.objective(cx,&vec![result_rho;design.parameter_count()],&weights).unwrap();
        assert!(((fresh-300.0)/1000.0-report.objective).abs()<1e-10);
        split.run(1e-5,4,work(),Some(cx)).unwrap();let replay=split.run(1e-5,146,work(),Some(cx)).unwrap();
        assert_eq!(replay,report,"shared physical caches must not change mathematical restart state");
    });
}

#[test]
fn a_native_temperature_constraint_drives_minimum_conductivity_with_budget_repair() {
    with_design(|design,cx| {
        let weights=vec![1.0/design.dofs().n() as f64;design.dofs().n()];
        let native=ConductivityNewtonModel::new(design,"native",weights.clone(),2).unwrap();
        let model=Uniform {model:&native,n:design.parameter_count()};
        let initial=model.value_gradient(&[0.0],Some(cx)).unwrap().value;
        let target=300.0+(initial-300.0)/1.5;
        let mut b=ProblemBuilder::new();let v=b.var("uniform log rho",Manifold::Rn {dim:1},Dims::NONE).unwrap();
        let r=b.var_ref(v).unwrap();let x=b.component(r,0).unwrap();let rho=b.exp(x).unwrap();
        b.objective(rho,Sense::Minimize,1.0).unwrap();
        let node=b.pde_residual("uniform",v,true,TEMPERATURE_DIMS).unwrap();
        let cap=b.konst(target,TEMPERATURE_DIMS).unwrap();let scale=b.konst(1000.0,TEMPERATURE_DIMS).unwrap();
        let difference=b.sub(node,cap).unwrap();let residual=b.div(difference,scale).unwrap();
        b.constraint(residual,ConstraintKind::LeZero,"temperature limit").unwrap();let p=b.finish();
        let oracle=ReverseProblem::new_with_physics(&p,ReverseLimits {max_nodes:128,max_scalar_slots:4096},
            &[PhysicsBinding {node,model:&model}],None).unwrap();
        let mut study=AugmentedNewtonStudy::new(&oracle,&[0.0],AugmentedNewtonOptions::default(),Some(cx)).unwrap();
        let before=study.inner().snapshot();let radius=study.inner().radius();
        assert!(study.run(1e-5,1,work(),Some(cx)).is_err());
        assert_eq!(study.inner().snapshot().x,before.x);assert_eq!(study.inner().radius(),radius);
        assert_eq!(native.reserved_solves(),2);assert_eq!(study.inner().snapshot().hv_evals,1);
        native.set_solve_limit(10000);
        let report=study.run(1e-5,200,work(),Some(cx)).unwrap();
        assert_eq!(report.stop,AugmentedNewtonStop::Converged,"{report:?}");
        assert!((report.objective-1.5).abs()<1e-4);assert!(report.kkt.within_tolerance(1e-5));
        let fresh=design.objective(cx,&vec![report.objective;design.parameter_count()],&weights).unwrap();
        assert!(fresh-target<0.01);assert!(((fresh-target)/1000.0-report.constraints[0]).abs()<1e-10);
        let tape=oracle.evaluate_cancellable(&[report.x.clone()],cx).unwrap();
        let g=tape.lagrangian_gradient_cancellable(&report.multipliers,cx).unwrap();
        assert!((g[0][0].abs()-report.kkt.stationarity).abs()<1e-10);
        assert!(native.reserved_solves()>2 && report.multiplier_updates>0);
    });
}
