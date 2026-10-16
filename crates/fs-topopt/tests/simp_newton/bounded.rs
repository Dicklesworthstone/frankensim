//! End-to-end constrained Newton using the real filtered tetrahedral model.
use super::*;
use fs_ascent::reverse_newton::{AugmentedNewtonBudget, AugmentedNewtonOptions, AugmentedNewtonStop, AugmentedNewtonStudy};
use fs_opt::reverse::physics::{PhysicsSample, PhysicsSignature};
use fs_opt::reverse::physics::box_map::LogisticBox;

// Only this test reduction ties cells together; every evaluation and derivative
// still runs the production 3-D filter, independent loads and structural solves.
#[derive(Debug)]
struct Uniform<'a> { native: &'a dyn PhysicsModel, n: usize }
impl PhysicsModel for Uniform<'_> {
    fn signature(&self) -> PhysicsSignature<'_> {
        PhysicsSignature { parameter_count: 1, ..self.native.signature() }
    }
    fn value_gradient(&self, x: &[f64], cx: Option<&Cx<'_>>) -> Result<PhysicsSample, PhysicsError> {
        if x.len() != 1 { return Err(PhysicsError::InvalidInput("uniform point shape")); }
        let sample = self.native.value_gradient(&vec![x[0];self.n], cx)?;
        Ok(PhysicsSample { value: sample.value, gradient: vec![sample.gradient.iter().sum()] })
    }
    fn second_order(&self) -> Option<&dyn PhysicsHessian> { Some(self) }
}
impl PhysicsHessian for Uniform<'_> {
    fn hessian_vector_product(&self, x: &[f64], v: &[f64], cx: Option<&Cx<'_>>) -> Result<Vec<f64>, PhysicsError> {
        if x.len() != 1 || v.len() != 1 { return Err(PhysicsError::InvalidInput("uniform direction shape")); }
        let hv = self.native.second_order().unwrap().hessian_vector_product(&vec![x[0];self.n], &vec![v[0];self.n], cx)?;
        Ok(vec![hv.iter().sum()])
    }
}
fn optimizer_budget() -> AugmentedNewtonBudget {
    AugmentedNewtonBudget { objective_evals: 600, hessian_products: 3000, multiplier_updates: 100 }
}

#[test]
fn bounded_heterogeneous_derivatives_keep_shared_physics_and_coordinate_curvature() {
    let (pipeline, el, a, b, volumes) = fixture(2.0); let n = el.cells();
    let native = SimpNewtonModel::new(&pipeline, &el,
        &[LoadCase { force: &a, weight: 0.7 }, LoadCase { force: &b, weight: 1.3 }],
        &volumes, options(), budget()).unwrap();
    let raw_c = native.compliance("c").unwrap(); let raw_v = native.volume("v").unwrap();
    let map = LogisticBox::new(&vec![0.05;n], &vec![0.95;n], Dims::NONE).unwrap();
    let c = map.bind(&raw_c).unwrap(); let v = map.bind(&raw_v).unwrap();
    let x: Vec<_> = (0..n).map(|i| 0.08*((i%9) as f64-4.0)).collect();
    let d: Vec<_> = (0..n).map(|i| 0.02*((i%5) as f64-2.0)).collect();
    let cs = c.value_gradient(&x, None).unwrap(); let work = native.work();
    let vs = v.value_gradient(&x, None).unwrap(); assert_eq!(native.work(), work);
    let ch = c.hessian_vector_product(&x, &d, None).unwrap();
    let vh = v.hessian_vector_product(&x, &d, None).unwrap();
    let eps = 1e-4;
    let plus: Vec<_> = x.iter().zip(&d).map(|(x,d)| x+eps*d).collect();
    let minus: Vec<_> = x.iter().zip(&d).map(|(x,d)| x-eps*d).collect();
    // Use the retained first-order pipeline, not the new Hessian implementation,
    // as the independent gradient oracle after the coordinate transformation.
    let (_, mut reference, _, _, _) = fixture(2.0);
    for (volume, expected) in [(false, &ch), (true, &vh)] {
        let mut transformed = Vec::new();
        for point in [&plus, &minus] {
            let rho = map.decode(point, None).unwrap();
            let gradient = if volume { pipeline.volume_and_gradient(&rho, &volumes).1 }
                else { pipeline.multi_load_compliance_and_gradient(&mut reference, &rho,
                    &[LoadCase { force: &a, weight: 0.7 }, LoadCase { force: &b, weight: 1.3 }]).gradient };
            transformed.push(gradient.iter().zip(&rho).map(|(g,r)| g*(r-0.05)*(0.95-r)/0.9).collect::<Vec<_>>());
        }
        let fd: Vec<_> = transformed[0].iter().zip(&transformed[1]).map(|(a,b)| (a-b)/(2.0*eps)).collect();
        close(expected, &fd, 5e-6);
    }
    // Evict through either view, then revisit the accepted point through both.
    v.value_gradient(&plus, None).unwrap();
    assert_eq!(c.value_gradient(&x, None).unwrap(), cs);
    assert_eq!(v.value_gradient(&x, None).unwrap(), vs);
    assert_eq!(c.hessian_vector_product(&x, &d, None).unwrap(), ch);
    let before = native.work();
    assert_eq!(v.hessian_vector_product(&x, &vec![0.0;n], None).unwrap(), vec![0.0;n]);
    assert_eq!(native.work(), before);
}

#[test]
fn native_simp_compliance_minimization_respects_physical_volume_and_split_resume() {
    with_cx(|_, cx| {
        let (pipeline, el, a, bload, volumes) = fixture(2.0); let n = el.cells();
        let loads = [LoadCase { force: &a, weight: 0.7 }, LoadCase { force: &bload, weight: 1.3 }];
        let native = SimpNewtonModel::new(&pipeline, &el, &loads, &volumes, options(), budget()).unwrap();
        let raw_c = native.compliance("c").unwrap(); let raw_v = native.volume("v").unwrap();
        let uc = Uniform { native: &raw_c, n }; let uv = Uniform { native: &raw_v, n };
        let map = LogisticBox::new(&[0.05], &[0.95], Dims::NONE).unwrap();
        let c = map.bind(&uc).unwrap(); let v = map.bind(&uv).unwrap();
        let start = map.encode(&[0.5], Some(cx)).unwrap();
        let base = c.value_gradient(&start, Some(cx)).unwrap().value;
        let mut builder = ProblemBuilder::new();
        let var = builder.var("bounded density", Manifold::Rn { dim: 1 }, Dims::NONE).unwrap();
        let cn = builder.pde_residual("c", var, true, COMPLIANCE_DIMS).unwrap();
        let vn = builder.pde_residual("v", var, true, Dims::NONE).unwrap();
        let scale = builder.konst(base, COMPLIANCE_DIMS).unwrap(); let objective = builder.div(cn, scale).unwrap();
        builder.objective(objective, Sense::Minimize, 1.0).unwrap();
        let cap = builder.konst(0.6, Dims::NONE).unwrap(); let residual = builder.sub(vn, cap).unwrap();
        builder.constraint(residual, ConstraintKind::LeZero, "physical material fraction").unwrap();
        let problem = builder.finish();
        let oracle = ReverseProblem::new_with_physics(&problem, ReverseLimits { max_nodes: 128, max_scalar_slots: 4096 },
            &[PhysicsBinding { node: cn, model: &c }, PhysicsBinding { node: vn, model: &v }], None).unwrap();
        let mut study = AugmentedNewtonStudy::new(&oracle, &start, AugmentedNewtonOptions::default(), Some(cx)).unwrap();
        let mut split = study.clone();
        let report = study.run(1e-5, 300, optimizer_budget(), Some(cx)).unwrap();
        assert_eq!(report.stop, AugmentedNewtonStop::Converged, "{report:?}");
        assert!(report.kkt.within_tolerance(1e-5));
        let rho = map.decode(&report.x, Some(cx)).unwrap();
        let (fresh_volume, _) = pipeline.volume_and_gradient(&vec![rho[0];n], &volumes);
        close(&[fresh_volume], &[0.6], 4e-5);
        let expected = (0.01+0.99*0.5_f64.powi(3))/(0.01+0.99*0.6_f64.powi(3));
        close(&[report.objective], &[expected], 1e-4);
        let (_, mut reference, _, _, _) = fixture(2.0);
        let fresh = pipeline.multi_load_compliance_and_gradient(&mut reference, &vec![rho[0];n], &loads);
        close(&[fresh.compliance/base], &[report.objective], 1e-8);
        close(&[fresh_volume-0.6], &report.constraints, 1e-11);
        assert!(report.multipliers[0] > 0.0);
        split.run(1e-5, 4, optimizer_budget(), Some(cx)).unwrap();
        assert_eq!(split.run(1e-5, 296, optimizer_budget(), Some(cx)).unwrap(), report);
    });
}

#[test]
fn minimum_material_uses_compliance_constraint_and_recovers_without_refunding_work() {
    with_cx(|_, cx| {
        let (pipeline, el, a, _, volumes) = fixture(2.0); let n = el.cells();
        let loads = [LoadCase { force: &a, weight: 1.0 }];
        let native = SimpNewtonModel::new(&pipeline, &el, &loads, &volumes, options(), budget()).unwrap();
        let raw_c = native.compliance("c").unwrap(); let raw_v = native.volume("v").unwrap();
        let uc = Uniform { native: &raw_c, n }; let uv = Uniform { native: &raw_v, n };
        let map = LogisticBox::new(&[0.05], &[0.95], Dims::NONE).unwrap();
        let c = map.bind(&uc).unwrap(); let v = map.bind(&uv).unwrap();
        let start = map.encode(&[0.5], Some(cx)).unwrap();
        let base = c.value_gradient(&start, Some(cx)).unwrap().value;
        let work = native.work();
        native.set_budget(ComplianceNewtonBudget { solve_reservations: work.solve_reservations, ..budget() });
        let limit = (0.01+0.99*0.5_f64.powi(3))/(0.01+0.99*0.6_f64.powi(3));
        let mut builder = ProblemBuilder::new();
        let var = builder.var("bounded density", Manifold::Rn { dim: 1 }, Dims::NONE).unwrap();
        let cn = builder.pde_residual("c", var, true, COMPLIANCE_DIMS).unwrap();
        let vn = builder.pde_residual("v", var, true, Dims::NONE).unwrap();
        builder.objective(vn, Sense::Minimize, 1.0).unwrap();
        let scale = builder.konst(base, COMPLIANCE_DIMS).unwrap(); let normalized = builder.div(cn, scale).unwrap();
        let cap = builder.konst(limit, Dims::NONE).unwrap(); let residual = builder.sub(normalized, cap).unwrap();
        builder.constraint(residual, ConstraintKind::LeZero, "compliance limit").unwrap();
        let problem = builder.finish();
        let oracle = ReverseProblem::new_with_physics(&problem, ReverseLimits { max_nodes: 128, max_scalar_slots: 4096 },
            &[PhysicsBinding { node: cn, model: &c }, PhysicsBinding { node: vn, model: &v }], None).unwrap();
        let mut study = AugmentedNewtonStudy::new(&oracle, &start, AugmentedNewtonOptions::default(), Some(cx)).unwrap();
        let accepted = study.inner().snapshot(); let radius = study.inner().radius();
        assert!(study.run(1e-5, 1, optimizer_budget(), Some(cx)).is_err());
        assert_eq!(study.inner().snapshot().x, accepted.x); assert_eq!(study.inner().radius(), radius);
        assert_eq!(study.inner().snapshot().hv_evals, accepted.hv_evals+1);
        assert_eq!(native.work(), work);
        native.set_budget(budget());
        let report = study.run(1e-5, 300, optimizer_budget(), Some(cx)).unwrap();
        assert_eq!(report.stop, AugmentedNewtonStop::Converged, "{report:?}");
        assert!(report.kkt.within_tolerance(1e-5));
        close(&[report.objective], &[0.6], 4e-5);
        let rho = map.decode(&report.x, Some(cx)).unwrap();
        let (_, mut reference, _, _, _) = fixture(2.0);
        let fresh = pipeline.multi_load_compliance_and_gradient(&mut reference, &vec![rho[0];n], &loads);
        close(&[fresh.compliance/base-limit], &report.constraints, 1e-8);
        let (volume, _) = pipeline.volume_and_gradient(&vec![rho[0];n], &volumes);
        close(&[volume], &[report.objective], 1e-11);
        assert!(native.work().solve_reservations > work.solve_reservations);
        let tape = oracle.evaluate_cancellable(&[report.x.clone()], cx).unwrap();
        let lagrangian = tape.lagrangian_gradient_cancellable(&report.multipliers, cx).unwrap();
        close(&[lagrangian[0][0].abs()], &[report.kkt.stationarity], 1e-10);
    });
}
