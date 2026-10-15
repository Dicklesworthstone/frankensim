//! Real eroded-compliance/dilated-volume constraints through the shared Newton path.
use fs_ascent::reverse_newton::{AugmentedNewtonBudget, AugmentedNewtonOptions, AugmentedNewtonStop, AugmentedNewtonStudy};
use fs_exec::Cx;
use fs_opt::reverse::physics::{PhysicsBinding, PhysicsError, PhysicsHessian, PhysicsModel, PhysicsSample, PhysicsSignature};
use fs_opt::reverse::physics::box_map::LogisticBox;
use fs_opt::reverse::ReverseLimits;
use fs_opt::{ConstraintKind, Manifold, ProblemBuilder, ReverseProblem, Sense};
use fs_qty::Dims;
use fs_topopt::compliance_newton::{COMPLIANCE_DIMS, ComplianceNewtonBudget};
use fs_topopt::elasticity::second_order::ComplianceSolveOptions;
use fs_topopt::filter::heaviside;
use fs_topopt::pipeline::LoadCase;
use fs_topopt::simp_newton::SimpNewtonModel;
use fs_topopt::{DensityElasticity, DensityFilter, RobustPipeline, SimpParams};

fn fixture() -> (RobustPipeline, DensityElasticity, Vec<f64>, Vec<f64>) {
    let (mesh, mut positions) = fs_feec::kuhn_cube(2);
    for p in &mut positions { p[0] *= p[0]; }
    let volumes = fs_feec::element_geometry(&mesh, &positions).vol_signed.iter().map(|v| v.abs()).collect();
    let el = DensityElasticity::new(&mesh, &positions, 1.0, 0.29, &|p| p[0] == 0.0);
    let mut force = vec![0.0; el.n()];
    for (i,p) in positions.iter().enumerate() { if p[0] == 1.0 { force[3*i+2] = -0.2; } }
    let pipeline = RobustPipeline { filter: DensityFilter::new(&mesh, &positions, 0.15),
        params: SimpParams { e_min: 1e-3, penal: 3.0, beta: 2.0, eta: 0.5 }, eta_offset: 0.1 };
    (pipeline, el, volumes, force)
}
fn budget() -> ComplianceNewtonBudget {
    ComplianceNewtonBudget { solve_reservations: 100_000, linear_iterations: 5_000_000 }
}
fn optimizer_budget() -> AugmentedNewtonBudget {
    AugmentedNewtonBudget { objective_evals: 800, hessian_products: 2000, multiplier_updates: 100 }
}
fn options() -> ComplianceSolveOptions { ComplianceSolveOptions { tolerance: 1e-11, max_iterations: 4000 } }

// Test-only reduction to an analytically solvable uniform material design.
// Every value and derivative still uses the actual 48-cell native physics.
#[derive(Debug)]
struct Uniform<'a> { inner: &'a dyn PhysicsModel, cells: usize }
impl PhysicsModel for Uniform<'_> {
    fn signature(&self) -> PhysicsSignature<'_> {
        PhysicsSignature { parameter_count: 1, ..self.inner.signature() }
    }
    fn value_gradient(&self, x: &[f64], cx: Option<&Cx<'_>>) -> Result<PhysicsSample, PhysicsError> {
        if x.len() != 1 { return Err(PhysicsError::InvalidInput("uniform point")); }
        let sample = self.inner.value_gradient(&vec![x[0];self.cells], cx)?;
        Ok(PhysicsSample { value: sample.value, gradient: vec![sample.gradient.iter().sum()] })
    }
    fn second_order(&self) -> Option<&dyn PhysicsHessian> { Some(self) }
}
impl PhysicsHessian for Uniform<'_> {
    fn hessian_vector_product(&self, x: &[f64], v: &[f64], cx: Option<&Cx<'_>>) -> Result<Vec<f64>, PhysicsError> {
        if x.len() != 1 || v.len() != 1 { return Err(PhysicsError::InvalidInput("uniform direction")); }
        let hv = self.inner.second_order().unwrap().hessian_vector_product(&vec![x[0];self.cells], &vec![v[0];self.cells], cx)?;
        Ok(vec![hv.iter().sum()])
    }
}

#[test]
fn both_robust_constraint_orientations_reach_the_uniform_optimum_and_replay() {
    for minimum_material in [false, true] {
        let (pipeline, el, volumes, force) = fixture();
        let model = SimpNewtonModel::new_robust(&pipeline, &el,
            &[LoadCase { force: &force, weight: 1.0 }], &volumes, options(), budget()).unwrap();
        assert_eq!(model.projection_thresholds(), (0.6,0.4));
        let raw_c = model.compliance("eroded").unwrap(); let raw_v = model.volume("dilated").unwrap();
        let uc = Uniform { inner: &raw_c, cells: el.cells() }; let uv = Uniform { inner: &raw_v, cells: el.cells() };
        let chart = LogisticBox::new(&[0.05], &[0.95], Dims::NONE).unwrap();
        let c = chart.bind(&uc).unwrap(); let v = chart.bind(&uv).unwrap();
        let x0 = chart.encode(&[0.5], None).unwrap(); let baseline = c.value_gradient(&x0,None).unwrap().value;
        let target = 0.52; let p = pipeline.params;
        let modulus = |r| p.e_min+(1.0-p.e_min)*fs_math::det::pow(heaviside(r,p.beta,p.eta+pipeline.eta_offset),p.penal);
        let target_c = modulus(0.5)/modulus(target);
        let target_v = heaviside(target,p.beta,p.eta-pipeline.eta_offset);
        let mut b = ProblemBuilder::new(); let x = b.var("uniform bounded density",Manifold::Rn { dim:1 },Dims::NONE).unwrap();
        let cn = b.pde_residual("eroded",x,true,COMPLIANCE_DIMS).unwrap();
        let vn = b.pde_residual("dilated",x,true,Dims::NONE).unwrap();
        let scale = b.konst(baseline,COMPLIANCE_DIMS).unwrap(); let nc = b.div(cn,scale).unwrap();
        let (objective, constrained, cap) = if minimum_material { (vn,nc,target_c) } else { (nc,vn,target_v) };
        b.objective(objective,Sense::Minimize,1.0).unwrap(); let cap = b.konst(cap,Dims::NONE).unwrap();
        let residual = b.sub(constrained,cap).unwrap(); b.constraint(residual,ConstraintKind::LeZero,"robust allowance").unwrap();
        let problem = b.finish();
        let oracle = ReverseProblem::new_with_physics(&problem,ReverseLimits {max_nodes:128,max_scalar_slots:4096},
            &[PhysicsBinding {node:cn,model:&c},PhysicsBinding {node:vn,model:&v}],None).unwrap();
        let mut whole = AugmentedNewtonStudy::new(&oracle,&x0,AugmentedNewtonOptions::default(),None).unwrap();
        let mut split = whole.clone();
        let report = whole.run(1e-5,200,optimizer_budget(),None).unwrap();
        assert_eq!(report.stop,AugmentedNewtonStop::Converged,"{report:?}");
        assert!(report.kkt.within_tolerance(1e-5)); assert!(report.multipliers[0]>0.0);
        let density = chart.decode(&report.x,None).unwrap()[0]; assert!((density-target).abs()<1e-4);
        let (fresh_pipeline, mut fresh_el, fresh_volumes, fresh_force) = fixture();
        let (fresh_c, _) = fresh_pipeline.eroded_compliance_and_gradient(&mut fresh_el,&vec![density;el.cells()],&fresh_force);
        let fields = fresh_pipeline.three_fields(&vec![density;el.cells()]);
        let fresh_v: f64 = fields.dilated.iter().zip(&fresh_volumes).map(|(r,w)|r*w).sum::<f64>()/fresh_volumes.iter().sum::<f64>();
        let (objective, constraint) = if minimum_material { (fresh_v,fresh_c/baseline-target_c) }
            else { (fresh_c/baseline,fresh_v-target_v) };
        assert!((report.objective-objective).abs()<1e-7);
        assert!((report.constraints[0]-constraint).abs()<1e-7);
        split.run(1e-5,4,optimizer_budget(),None).unwrap();
        assert_eq!(split.run(1e-5,196,optimizer_budget(),None).unwrap(),report);
    }
}

#[test]
fn robust_views_share_work_and_failed_replacement_keeps_the_accepted_point() {
    let (pipeline, el, volumes, force) = fixture(); let n = el.cells();
    let model = SimpNewtonModel::new_robust(&pipeline, &el,
        &[LoadCase {force:&force,weight:1.0}],&volumes,options(),budget()).unwrap();
    let c = model.compliance("eroded").unwrap(); let v = model.volume("dilated").unwrap();
    let x = vec![0.5;n]; let c0 = c.value_gradient(&x,None).unwrap(); let v0 = v.value_gradient(&x,None).unwrap();
    assert_eq!(model.work().started_solves,4); assert_eq!(model.work().solve_reservations,4);
    let before = model.work();
    model.set_budget(ComplianceNewtonBudget {linear_iterations:before.linear_iterations,..budget()});
    assert!(matches!(v.value_gradient(&vec![0.52;n],None),Err(PhysicsError::Budget {..})));
    let failed = model.work(); assert_eq!(failed.solve_reservations,before.solve_reservations+4);
    assert_eq!(failed.linear_iterations,before.linear_iterations);
    assert_eq!(c.value_gradient(&x,None).unwrap(),c0); assert_eq!(v.value_gradient(&x,None).unwrap(),v0);
    assert_eq!(model.work(),failed);
    model.set_budget(budget()); let changed = v.value_gradient(&vec![0.52;n],None).unwrap();
    assert!(changed.value>v0.value); assert_eq!(model.work().solve_reservations,failed.solve_reservations+4);
    assert_eq!(c.value_gradient(&x,None).unwrap(),c0);
    let before = model.work();
    assert_eq!(v.hessian_vector_product(&x,&vec![0.0;n],None).unwrap(),vec![0.0;n]);
    assert_eq!(model.work(),before);
}
