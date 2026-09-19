//! The actual optimizer, sensitivity chain, and goal-refinement consumer.
use std::ops::ControlFlow;
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_cutfem::elastic3::{ElasticityOptions3, adaptive::AdaptiveElasticity3};
use fs_cutfem::elastic3::adaptive::enrichment::precondition::{AdaptiveSolveOptions3, AdaptiveSolveSpace3};
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_topopt::{EvaluationStop, MultiLoadOcOptions, MultiLoadOcTermination, SimpParams, SolveBudget, SolveControl, SolveProgress};
use fs_topopt::pipeline::LoadCase;
use fs_topopt::sdf3::{CutDensityStudy3, Sdf3Elasticity, controlled_sdf3_optimality_criteria};
use fs_topopt::sdf3_goal::{GoalBodyLoad3, GoalPreconditioner3, GoalRefinementOptions3};

struct Slab;
impl CutSdf3 for Slab {
    fn value(&self, p: [f64; 3]) -> f64 { p[2] - 0.73 }
    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval { Interval::new(lo[2], hi[2]) - Interval::new(0.73, 0.73) }
    fn derivative_enclose(&self, _: [f64; 3], _: [f64; 3], a: HeightAxis) -> Interval {
        let d = if a == HeightAxis::Z { 1.0 } else { 0.0 }; Interval::new(d, d)
    }
}
fn y(_: [f64; 3]) -> [f64; 3] { [0.0, -1.0, 0.0] }
fn z(_: [f64; 3]) -> [f64; 3] { [0.0, 0.0, -1.0] }
fn build(level: u8) -> AdaptiveElasticity3 {
    let tree = Octree3::uniform(level, 4, 4096).unwrap();
    let mut p = |_| ControlFlow::Continue(());
    let mut q = QuadratureControl3::new(QuadratureOptions3 { depth: 1, ..Default::default() }, &mut p).unwrap();
    AdaptiveElasticity3::build(HexCell::try_new([0.0; 3], [1.0; 3]).unwrap(), &tree, &Slab,
        &IsotropicElastic::new(1.0, 0.3, 1.0).unwrap(), &|p| p[0] == 0.0, ElasticityOptions3::default(), &mut q).unwrap()
}
fn space() -> AdaptiveSolveSpace3 {
    AdaptiveSolveSpace3::two_level(build(1), &build(0), AdaptiveSolveOptions3::default(), || ControlFlow::Continue(())).unwrap()
}
fn forces(op: &AdaptiveElasticity3) -> (Vec<f64>, Vec<f64>) {
    (op.body_load(&y, || ControlFlow::Continue(())).unwrap(), op.body_load(&z, || ControlFlow::Continue(())).unwrap())
}
fn loads<'a>(f: &'a (Vec<f64>, Vec<f64>)) -> [LoadCase<'a>; 2] {
    [LoadCase { force: &f.0, weight: 0.3 }, LoadCase { force: &f.1, weight: 0.7 }]
}
fn params() -> SimpParams { SimpParams { beta: 2.0, penal: 3.0, ..Default::default() } }
fn close(a: &[f64], b: &[f64], tol: f64) {
    assert_eq!(a.len(), b.len()); let scale = b.iter().map(|v| v.abs()).fold(1e-30_f64, f64::max);
    assert!(a.iter().zip(b).all(|(a, b)| (a-b).abs() <= tol*scale));
}

#[test]
fn independent_loads_share_one_fresh_preparation_per_density() {
    let backend = space(); let f = forces(backend.elasticity()); let nc = backend.coarse_dofs();
    let mut study = CutDensityStudy3::new(backend, 0.15, params());
    let mut identity = CutDensityStudy3::new(build(1), 0.15, params());
    let mut preparations = 0;
    let mut callback = |p: SolveProgress| { if p.stage == "sdf3-preconditioner-start" { preparations += 1; } ControlFlow::Continue(()) };
    let mut c = SolveControl::new(SolveBudget::default(), &mut callback);
    let mut plain = |_| ControlFlow::Continue(()); let mut reference = SolveControl::new(SolveBudget::default(), &mut plain);
    for stage in 0..3 {
        let rho: Vec<_> = (0..study.cells()).map(|i| 0.35 + 0.02 * ((i+stage)%9) as f64).collect();
        let actual = study.evaluate(&rho, &loads(&f), &mut c).unwrap();
        let expected = identity.evaluate(&rho, &loads(&f), &mut reference).unwrap();
        close(&[actual.objective.compliance], &[expected.objective.compliance], 1e-8);
        close(&actual.objective.gradient, &expected.objective.gradient, 1e-7);
        for (a, b) in actual.objective.displacements.iter().zip(&expected.objective.displacements) { close(a, b, 1e-7); }
        assert_eq!(actual.projected_rho, expected.projected_rho);
        assert_eq!(c.work().preconditioner_operator_applications, (stage+1)*nc);
    }
    assert_eq!(reference.work().preconditioner_operator_applications, 0);
    assert_eq!(c.work().linear_solves, reference.work().linear_solves);
    drop(c); assert_eq!(preparations, 3);
}

#[test]
fn prepared_full_chain_coordinate_gradients_match_independent_resolves() {
    let backend = space(); let f = forces(backend.elasticity()); let mut study = CutDensityStudy3::new(backend, 0.15, params());
    let rho: Vec<_> = (0..study.cells()).map(|i| 0.35 + 0.03*i as f64).collect();
    let mut p = |_| ControlFlow::Continue(()); let mut c = SolveControl::new(SolveBudget::default(), &mut p);
    let base = study.evaluate(&rho, &loads(&f), &mut c).unwrap();
    for i in 0..rho.len() {
        let mut a = rho.clone(); let mut b = rho.clone(); a[i] += 1e-4; b[i] -= 1e-4;
        let plus = study.evaluate(&a, &loads(&f), &mut c).unwrap(); let minus = study.evaluate(&b, &loads(&f), &mut c).unwrap();
        let fd = (plus.objective.compliance-minus.objective.compliance)/2e-4;
        assert!((fd-base.objective.gradient[i]).abs() < 2e-4*fd.abs());
    }
}

#[test]
fn optimizer_descends_with_fixed_geometry_and_replays_accepted_fields() {
    let backend = space(); let f = forces(backend.elasticity()); let leaves = backend.elasticity().leaves().to_vec();
    let entries = backend.transfer_entries(); let mut study = CutDensityStudy3::new(backend, 0.15, params()); let rho = vec![0.5; study.cells()];
    let mut p = |_| ControlFlow::Continue(()); let mut c = SolveControl::new(SolveBudget::default(), &mut p);
    let report = controlled_sdf3_optimality_criteria(&mut study, &loads(&f), &rho,
        MultiLoadOcOptions { max_iterations: 3, ..Default::default() }, &mut c);
    assert!(report.history.len() > 1, "{report:?}");
    assert!(report.history.last().unwrap().compliance < report.history[0].compliance);
    for pair in report.history.windows(2) { assert!(pair[1].compliance <= pair[0].compliance); }
    assert!(report.history.iter().all(|row| row.volume_fraction <= 0.5+1e-8));
    let again = study.evaluate(&report.rho, &loads(&f), &mut c).unwrap();
    assert_eq!(again.objective.displacements, report.displacements);
    assert_eq!(again.projected_rho, report.projected_rho);
    assert_eq!(study.operator().elasticity().leaves(), leaves); assert_eq!(study.operator().transfer_entries(), entries);
}

#[test]
fn cancelled_trial_setup_is_charged_and_restores_the_accepted_operator() {
    let backend = space(); let nc = backend.coarse_dofs(); let f = forces(backend.elasticity());
    let mut study = CutDensityStudy3::new(backend, 0.15, params()); let rho = vec![0.5; study.cells()];
    let mut p = |s: SolveProgress| if s.work.preconditioner_operator_applications >= nc+3 { ControlFlow::Break(()) } else { ControlFlow::Continue(()) };
    let mut c = SolveControl::new(SolveBudget::default(), &mut p);
    let report = controlled_sdf3_optimality_criteria(&mut study, &loads(&f), &rho,
        MultiLoadOcOptions { max_iterations: 3, ..Default::default() }, &mut c);
    assert_eq!(report.termination, MultiLoadOcTermination::Cancelled);
    assert_eq!(report.work.preconditioner_operator_applications, nc+3); assert_eq!(report.history.len(), 1); assert_eq!(report.rho, rho);
    let scales = study.operator().scales().to_vec();
    let mut p = |_| ControlFlow::Continue(()); let mut c = SolveControl::new(SolveBudget::default(), &mut p);
    let retry = study.evaluate(&rho, &loads(&f), &mut c).unwrap();
    assert_eq!(retry.objective.displacements, report.displacements); assert_eq!(study.operator().scales(), scales);
}

#[test]
fn stops_after_setup_or_inside_solves_never_publish_a_new_density() {
    for target in ["sdf3-preconditioner-ready", "sdf3-elasticity", "sdf3-evaluation-publish"] {
        let backend = space(); let f = forces(backend.elasticity()); let mut study = CutDensityStudy3::new(backend, 0.15, params());
        let mut p = |_| ControlFlow::Continue(()); let mut c = SolveControl::new(SolveBudget::default(), &mut p);
        let baseline = study.evaluate(&vec![0.5; study.cells()], &loads(&f), &mut c).unwrap(); let scales = study.operator().scales().to_vec();
        let mut p = |s: SolveProgress| if s.stage == target && (target != "sdf3-elasticity" || s.solve_iterations > 0) { ControlFlow::Break(()) } else { ControlFlow::Continue(()) };
        let mut c = SolveControl::new(SolveBudget::default(), &mut p);
        assert!(matches!(study.evaluate(&vec![0.4; study.cells()], &loads(&f), &mut c), Err(EvaluationStop::Cancelled)));
        assert_eq!(study.operator().scales(), scales);
        assert_eq!(c.work().preconditioner_operator_applications, study.operator().coarse_dofs());
        let mut p = |_| ControlFlow::Continue(()); let mut c = SolveControl::new(SolveBudget::default(), &mut p);
        assert_eq!(study.evaluate(&vec![0.5; study.cells()], &loads(&f), &mut c).unwrap().objective.displacements, baseline.objective.displacements);
    }
}

#[test]
fn outer_budget_and_diagonal_refusal_preserve_previous_scales() {
    let backend = space(); let f = forces(backend.elasticity()); let mut study = CutDensityStudy3::new(backend, 0.15, params());
    let rho = vec![0.5; study.cells()]; let mut before_setup = 0;
    let mut p = |s: SolveProgress| { if s.stage == "sdf3-preconditioner-start" { before_setup = s.work.linear_iterations; } ControlFlow::Continue(()) };
    let mut c = SolveControl::new(SolveBudget::default(), &mut p); study.evaluate(&rho, &loads(&f), &mut c).unwrap(); drop(c);
    let scales = study.operator().scales().to_vec(); let mut p = |_| ControlFlow::Continue(());
    let mut c = SolveControl::new(SolveBudget { total_iterations: before_setup+1, ..Default::default() }, &mut p);
    assert!(matches!(study.evaluate(&rho, &loads(&f), &mut c), Err(EvaluationStop::TotalBudget { .. })));
    assert_eq!(c.work().linear_iterations, before_setup+1); assert_eq!(c.work().preconditioner_operator_applications, study.operator().coarse_dofs());
    assert_eq!(study.operator().scales(), scales);
    let mut bad = CutDensityStudy3::new(AdaptiveSolveSpace3::jacobi(build(1), 0), 0.15, params());
    let previous = bad.operator().scales().to_vec(); let mut p = |_| ControlFlow::Continue(()); let mut c = SolveControl::new(SolveBudget::default(), &mut p);
    assert!(matches!(bad.evaluate(&rho, &loads(&f), &mut c), Err(EvaluationStop::Breakdown { stage: "sdf3-preconditioner" })));
    assert_eq!(bad.operator().scales(), previous);
}

#[test]
fn preconditioned_accepted_fields_feed_the_existing_dwr_consumer() {
    let backend = space(); let f = forces(backend.elasticity()); let mut study = CutDensityStudy3::new(backend, 0.15, params());
    let mut p = |_| ControlFlow::Continue(()); let mut c = SolveControl::new(SolveBudget::default(), &mut p);
    let accepted = study.evaluate(&vec![0.5; study.cells()], &loads(&f), &mut c).unwrap();
    let scales = study.operator().scales().to_vec();
    let estimate = study.estimate_compliance_enrichment(build(2), &[GoalBodyLoad3 { density: &y, weight: 0.3 }, GoalBodyLoad3 { density: &z, weight: 0.7 }],
        &accepted.objective.displacements, GoalRefinementOptions3 { preconditioner: GoalPreconditioner3::Jacobi { max_contributions: 100_000_000 }, ..Default::default() }, &mut c).unwrap();
    close(&[estimate.coarse_value], &[accepted.objective.compliance], 1e-8);
    assert_eq!(estimate.cases.len(), 2); assert_eq!(study.operator().scales(), scales);
    assert!(!estimate.mark(0.5, 2, || ControlFlow::Continue(())).unwrap().marked.is_empty());
}
