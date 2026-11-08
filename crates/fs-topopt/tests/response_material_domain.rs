//! Real response optimizers must backtrack physical strain-domain exits.
#![cfg(feature = "cutfem-marquee")]

use fs_ascent::projected_al::{ProjectedAlError, ProjectedAlOptions};
use fs_ascent::sqp::SqpError;
use fs_cutfem::elastic3::adaptive::AdaptiveElasticity3;
use fs_cutfem::elastic3::adaptive::enrichment::precondition::AdaptiveSolveSpace3;
use fs_cutfem::elastic3::{ElasticityFormulation3, ElasticityOptions3};
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::QuadratureControl3;
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_solver::LinearOp;
use fs_topopt::sdf3::response::limits::design::{MinimumMaterialOptions3, MinimumMaterialStudy3};
use fs_topopt::sdf3::response::limits::strength::design::{
    ResponseStressDesignOptions3, ResponseStressDesignStudy3,
};
use fs_topopt::sdf3::response::limits::strength::ResponseStressError3;
use fs_topopt::sdf3::response::limits::{
    BoundedResponse3, ResponseLimit3, ResponseLimitCase3, ResponseLimitSense3,
};
use fs_topopt::sdf3::response::{
    ProjectedResponseOptions3, ProjectedResponseStudy3, ResponseCase3, ResponseDesignOptions3,
    ResponseDesignStudy3, ResponseError3, ResponseOptions3, ResponseTarget3,
};
use fs_topopt::sdf3::{CutDensityStudy3, Sdf3Elasticity};
use fs_topopt::{EvaluationStop, SimpParams, SolveBudget, SolveControl, SolveProgress};
use std::{cell::Cell, ops::ControlFlow};

type Study = CutDensityStudy3<AdaptiveSolveSpace3>;
struct Full;
impl CutSdf3 for Full {
    fn value(&self, _: [f64; 3]) -> f64 { -1.0 }
    fn enclose(&self, _: [f64; 3], _: [f64; 3]) -> Interval {
        Interval::new(-1.0, -1.0)
    }
    fn derivative_enclose(&self, _: [f64; 3], _: [f64; 3], _: HeightAxis) -> Interval {
        Interval::new(0.0, 0.0)
    }
}
fn poll() -> ControlFlow<()> { ControlFlow::Continue(()) }
fn fixture(limit: f64, formulation: ElasticityFormulation3) -> (Study, Vec<f64>, Vec<f64>) {
    let mut callback = |_| poll();
    let mut quadrature = QuadratureControl3::new(Default::default(), &mut callback).unwrap();
    let op = AdaptiveElasticity3::build(
        HexCell::try_new([0.0; 3], [1.0; 3]).unwrap(),
        &Octree3::uniform(0, 2, 64).unwrap(), &Full,
        &IsotropicElastic::new(10.0, 0.0, limit).unwrap(),
        &|p| p[0] == 0.0,
        ElasticityOptions3 { formulation, ghost_gamma: 0.0, ..Default::default() },
        &mut quadrature,
    ).unwrap();
    // Manufacture a uniform axial solution with the actual elastic operator.
    // Under a uniform SIMP scale s, eps_xx = 0.005/s. The source limit 0.01
    // requires rho >= ((0.5-0.001)/0.999)^(1/3), approximately 0.7934.
    let u: Vec<_> = op.nodes().iter().flat_map(|p| [0.005*p[0], 0.0, 0.0]).collect();
    let q: Vec<_> = op.nodes().iter().flat_map(|p| [0.25*p[0], 0.0, 0.0]).collect();
    let mut force = vec![0.0; op.n()];
    op.apply(&u, &mut force);
    let params = SimpParams { e_min: 1e-3, penal: 3.0, beta: 0.0, eta: 0.5 };
    (CutDensityStudy3::new(AdaptiveSolveSpace3::jacobi(op, 100_000_000), 0.0, params), force, q)
}
fn optimizer() -> ProjectedAlOptions {
    ProjectedAlOptions { max_evaluations: 128, max_backtracks: 20, ..Default::default() }
}
fn response_options() -> ResponseOptions3 {
    ResponseOptions3 { volume_weight: 1.0, ..Default::default() }
}
fn projected_options() -> ProjectedResponseOptions3 {
    ProjectedResponseOptions3 { response: response_options(), volume_cap: 1.0,
        density_floor: 0.2, optimizer: optimizer(), ..Default::default() }
}
fn sqp_options() -> ResponseDesignOptions3 {
    ResponseDesignOptions3 { response: response_options(), volume_cap: 1.0,
        density_floor: 0.2, ..Default::default() }
}
fn minimum_options() -> MinimumMaterialOptions3 {
    MinimumMaterialOptions3 { density_floor: 0.2, optimizer: optimizer(), ..Default::default() }
}
fn joint_options() -> ResponseStressDesignOptions3 {
    ResponseStressDesignOptions3 { density_floor: 0.2, optimizer: optimizer(), ..Default::default() }
}
fn target(q: &[f64]) -> [ResponseTarget3<'_>; 1] {
    // Zero observation weight cannot suppress physical admission of the case.
    [ResponseTarget3 { q, target: 0.0, scale: 1.0, weight: 0.0 }]
}
fn limit(q: &[f64]) -> [ResponseLimit3<'_>; 1] {
    [ResponseLimit3 { response: BoundedResponse3::Displacement(q),
        sense: ResponseLimitSense3::AtMost, bound: 1.0, scale: 1.0 }]
}
fn check_endpoint(study: &Study, rho: &[f64], fields: &[Vec<f64>]) {
    assert_eq!(rho.len(), 1);
    assert!(rho[0] < 0.9 && rho[0] > 0.793);
    let scale = 0.001 + 0.999*rho[0].powi(3);
    assert!((study.operator().scales()[0] - scale).abs() < 1e-12);
    for u in fields {
        study.operator().admit_material_domain(u, poll).unwrap();
        for (p, u) in study.operator().elasticity().nodes().iter().zip(u.chunks_exact(3)) {
            assert!((u[0] - 0.005*p[0]/scale).abs() < 1e-10);
            assert!(u[1].abs() < 1e-10 && u[2].abs() < 1e-10);
        }
    }
}

#[test]
fn sqp_response_fitting_backtracks_real_material_exits() {
    for formulation in [ElasticityFormulation3::Displacement, ElasticityFormulation3::BBar] {
        let (mut study, force, q) = fixture(0.01, formulation);
        let targets = target(&q);
        let cases = [ResponseCase3 { force: &force, prescribed: None, targets: &targets }];
        let mut callback = |_| poll();
        let mut control = SolveControl::new(Default::default(), &mut callback);
        let mut session = ResponseDesignStudy3::new(
            &mut study, &cases, &[0.9], sqp_options(), &mut control).unwrap();
        let before = session.work();
        session.run(1).unwrap();
        assert_eq!(session.iterations(), 1);
        assert!(session.evaluations() > 2);
        assert!(session.work().linear_solves > before.linear_solves + 2);
        let accepted = session.accepted();
        check_endpoint(session.study(), &accepted.rho, &accepted.displacements);
        assert_eq!(session.point(), accepted.rho);
    }
}

#[test]
fn projected_response_fitting_charges_rejected_trials_and_keeps_accepted_fields() {
    for formulation in [ElasticityFormulation3::Displacement, ElasticityFormulation3::BBar] {
        let (mut study, force, q) = fixture(0.01, formulation);
        let targets = target(&q);
        let cases = [ResponseCase3 { force: &force, prescribed: None, targets: &targets }];
        let mut callback = |_| poll();
        let mut control = SolveControl::new(Default::default(), &mut callback);
        let mut session = ProjectedResponseStudy3::new(
            &mut study, &cases, &[0.9], projected_options(), &mut control).unwrap();
        let before = session.work();
        let report = session.run(1).unwrap();
        assert_eq!(report.work.iterations, 1);
        assert!(report.work.rejected_trials > 0);
        assert!(session.work().linear_solves > before.linear_solves + 2);
        let accepted = session.accepted();
        check_endpoint(session.study(), &accepted.rho, &accepted.displacements);
        assert_eq!(session.point(), accepted.rho);
        assert!(session.best_feasible().unwrap().volume_fraction < 0.9);
    }
}

#[test]
fn response_limits_and_joint_strength_backtrack_without_losing_their_caps() {
    for formulation in [ElasticityFormulation3::Displacement, ElasticityFormulation3::BBar] {
        let (mut study, force, q) = fixture(0.01, formulation);
        let limits = limit(&q);
        let cases = [ResponseLimitCase3 { force: &force, prescribed: None, limits: &limits }];
        let mut callback = |_| poll();
        let mut control = SolveControl::new(Default::default(), &mut callback);
        {
            let mut session = MinimumMaterialStudy3::new(
                &mut study, &cases, &[0.9], minimum_options(), &mut control).unwrap();
            let before = session.work();
            let report = session.run(1).unwrap();
            assert_eq!(report.work.iterations, 1);
            assert!(report.work.rejected_trials > 0);
            assert!(session.work().linear_solves > before.linear_solves + 2);
            let accepted = session.accepted();
            check_endpoint(session.study(), &accepted.rho, &accepted.displacements);
            assert!(accepted.feasible(1e-6));
            assert_eq!(session.history().last().unwrap().volume_fraction, accepted.volume_fraction);
        }
        let mut session = ResponseStressDesignStudy3::new(
            &mut study, &cases, &[0.9], joint_options(), &mut control).unwrap();
        let before = session.work();
        let report = session.run(1).unwrap();
        assert_eq!(report.work.iterations, 1);
        assert!(report.work.rejected_trials > 0);
        assert!(session.work().linear_solves > before.linear_solves + 2);
        let accepted = session.accepted();
        check_endpoint(session.study(), &accepted.response.rho, &accepted.response.displacements);
        assert_eq!(accepted.response.displacements, accepted.stress.displacements);
        assert!(accepted.feasible(1e-6));
        assert_eq!(session.point(), accepted.response.rho);
    }
}

#[test]
fn standalone_initialization_and_restored_endpoints_still_refuse_material_exits() {
    let (mut study, force, q) = fixture(0.01, ElasticityFormulation3::Displacement);
    let targets = target(&q);
    let cases = [ResponseCase3 { force: &force, prescribed: None, targets: &targets }];
    let limits = limit(&q);
    let bounded = [ResponseLimitCase3 { force: &force, prescribed: None, limits: &limits }];
    let mut callback = |_| poll();
    let mut control = SolveControl::new(Default::default(), &mut callback);
    let incoming = study.operator().scales().to_vec();
    assert!(matches!(study.evaluate_responses(&[0.2], &cases, response_options(), &mut control),
        Err(ResponseError3::MaterialDomain(_))));
    assert!(matches!(ResponseDesignStudy3::new(&mut study, &cases, &[0.2], sqp_options(), &mut control),
        Err(SqpError::Evaluation(ResponseError3::MaterialDomain(_)))));
    assert!(matches!(ProjectedResponseStudy3::new(&mut study, &cases, &[0.2], projected_options(), &mut control),
        Err(ProjectedAlError::Evaluation(ResponseError3::MaterialDomain(_)))));
    assert!(matches!(MinimumMaterialStudy3::new(&mut study, &bounded, &[0.2], minimum_options(), &mut control),
        Err(ProjectedAlError::Evaluation(ResponseError3::MaterialDomain(_)))));
    assert!(matches!(ResponseStressDesignStudy3::new(&mut study, &bounded, &[0.2], joint_options(), &mut control),
        Err(ProjectedAlError::Evaluation(ResponseStressError3::Response(ResponseError3::MaterialDomain(_))))));
    assert_eq!(study.operator().scales(), incoming);

    let projected = ProjectedResponseStudy3::new(
        &mut study, &cases, &[0.9], projected_options(), &mut control).unwrap().checkpoint();
    let mut minimum = MinimumMaterialStudy3::new(
        &mut study, &bounded, &[0.9], minimum_options(), &mut control).unwrap().checkpoint();
    let joint = ResponseStressDesignStudy3::new(
        &mut study, &bounded, &[0.9], joint_options(), &mut control).unwrap().checkpoint();
    // A saved state cannot authorize a displacement outside the supplied
    // original material's stricter domain. Restoration remains a fresh solve.
    let (mut changed, _, _) = fixture(0.002, ElasticityFormulation3::Displacement);
    let incoming = changed.operator().scales().to_vec();
    assert!(matches!(ProjectedResponseStudy3::restore(
        &mut changed, &cases, projected, projected_options(), &mut control),
        Err(ProjectedAlError::Evaluation(ResponseError3::MaterialDomain(_)))));
    assert!(matches!(MinimumMaterialStudy3::restore(
        &mut changed, &bounded, &mut minimum, minimum_options(), &mut control),
        Err(ProjectedAlError::Evaluation(ResponseError3::MaterialDomain(_)))));
    assert!(matches!(ResponseStressDesignStudy3::restore(
        &mut changed, &bounded, joint, joint_options(), &mut control),
        Err(ProjectedAlError::Evaluation(ResponseStressError3::Response(ResponseError3::MaterialDomain(_))))));
    assert_eq!(changed.operator().scales(), incoming);
}

#[test]
fn cancellation_and_work_exhaustion_are_not_unavailable_material_trials() {
    let (mut study, force, q) = fixture(0.01, ElasticityFormulation3::Displacement);
    let targets = target(&q);
    let cases = [ResponseCase3 { force: &force, prescribed: None, targets: &targets }];
    let stop = Cell::new(false);
    let mut callback = |progress: SolveProgress| {
        if stop.get() && progress.stage == "sdf3-material-domain" { ControlFlow::Break(()) }
        else { poll() }
    };
    let mut control = SolveControl::new(Default::default(), &mut callback);
    let initialization_iterations;
    {
        let mut session = ProjectedResponseStudy3::new(
            &mut study, &cases, &[0.9], projected_options(), &mut control).unwrap();
        initialization_iterations = session.work().linear_iterations;
        let before = session.accepted().clone();
        let before_scales = session.study().operator().scales().to_vec();
        stop.set(true);
        assert!(matches!(session.run(1),
            Err(ProjectedAlError::Evaluation(ResponseError3::Evaluation(EvaluationStop::Cancelled)))));
        assert_eq!(session.point(), before.rho);
        assert_eq!(session.study().operator().scales(), before_scales);
        assert_eq!(session.accepted().displacements, before.displacements);
        assert_eq!(session.optimizer_work().rejected_trials, 0);
    }

    let (mut study, _, _) = fixture(0.01, ElasticityFormulation3::Displacement);
    let mut callback = |_| poll();
    let mut control = SolveControl::new(
        SolveBudget { total_iterations: initialization_iterations, ..Default::default() },
        &mut callback);
    let mut session = ProjectedResponseStudy3::new(
        &mut study, &cases, &[0.9], projected_options(), &mut control).unwrap();
    let before = session.accepted().clone();
    let before_scales = session.study().operator().scales().to_vec();
    assert!(matches!(session.run(1),
        Err(ProjectedAlError::Evaluation(ResponseError3::Evaluation(EvaluationStop::TotalBudget { .. })))));
    assert_eq!(session.point(), before.rho);
    assert_eq!(session.study().operator().scales(), before_scales);
    assert_eq!(session.accepted().displacements, before.displacements);
    assert_eq!(session.optimizer_work().rejected_trials, 0);
}
