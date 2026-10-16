use std::cell::Cell;
use fs_ascent::newton_continuation::{ContinuationError, ContinuationStop, ContinuationWork,
    NewtonContinuation, NewtonStage};
use fs_ascent::reverse_newton::{AugmentedNewtonBudget, AugmentedNewtonOptions, AugmentedNewtonStop};
use fs_ascent::trust::KrylovOptions;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey, TilePool};
use fs_opt::reverse::physics::{PhysicsBinding, PhysicsError, PhysicsHessian, PhysicsModel,
    PhysicsSample, PhysicsSignature};
use fs_opt::reverse::ReverseLimits;
use fs_opt::{ConstraintKind, Manifold, Problem, ProblemBuilder, ReverseProblem, Sense};
use fs_qty::Dims;

fn problem(target: f64, constrained: bool, name: &str) -> Problem {
    let mut b = ProblemBuilder::new();
    let var = b.var(name, Manifold::Rn { dim: 1 }, Dims::NONE).unwrap();
    let r = b.var_ref(var).unwrap(); let x = b.component(r, 0).unwrap();
    let t = b.konst(target, Dims::NONE).unwrap();
    let root = if constrained {
        let c = b.sub(t, x).unwrap();
        b.constraint(c, ConstraintKind::LeZero, "lower limit").unwrap();
        b.powi(x, 2).unwrap()
    } else { let r = b.sub(x, t).unwrap(); b.powi(r, 2).unwrap() };
    b.objective(root, Sense::Minimize, 1.0).unwrap(); b.finish()
}
fn limits() -> ReverseLimits { ReverseLimits { max_nodes: 128, max_scalar_slots: 1024 } }
fn budget() -> AugmentedNewtonBudget {
    AugmentedNewtonBudget { objective_evals: 2000, hessian_products: 4000, multiplier_updates: 100 }
}
fn stage<'a, 'p>(oracle: &'a ReverseProblem<'p>) -> NewtonStage<'a, 'p> {
    NewtonStage { oracle, tolerance: 1e-7, options: AugmentedNewtonOptions::default(), krylov: KrylovOptions::default() }
}

#[test]
fn split_schedule_matches_uninterrupted_physics_and_work_at_stage_boundary() {
    let problems: Vec<_> = [0.0, 1.0, 2.0].into_iter().map(|t| problem(t, false, "x")).collect();
    let oracles: Vec<_> = problems.iter().map(|p| ReverseProblem::new(p, limits()).unwrap()).collect();
    let stages: Vec<_> = oracles.iter().map(stage).collect();
    let mut whole = NewtonContinuation::new(&stages, &[0.0], None).unwrap();
    let zero = whole.run(0, budget(), None).unwrap();
    assert_eq!(zero.work, ContinuationWork::default()); assert!(zero.solution.is_none());
    let mut split = whole.clone();
    let final_report = whole.run(100, budget(), None).unwrap();
    assert_eq!(final_report.stop, ContinuationStop::Complete);
    assert!((final_report.solution.as_ref().unwrap().x[0] - 2.0).abs() < 1e-7);
    let boundary = split.run(1, budget(), None).unwrap();
    assert_eq!(boundary.stage, 1); assert_eq!(boundary.stop, ContinuationStop::IterationCap);
    let held = split.work(); split.run(0, budget(), None).unwrap(); assert_eq!(split.work(), held);
    let resumed = split.run(99, budget(), None).unwrap();
    assert_eq!(resumed, final_report);
    assert_eq!(split.completed_stages(), whole.completed_stages());
    assert_eq!(resumed.work.stage_initializations, 3);
    let reports = whole.completed_stages().iter().chain(resumed.solution.iter());
    assert_eq!(resumed.work.objective_evals, reports.clone().map(|r| r.objective_evals).sum::<usize>());
    assert_eq!(resumed.work.hessian_products, reports.map(|r| r.hessian_products).sum::<usize>());
}

#[test]
fn global_hessian_and_entry_limits_are_not_reissued_to_each_model() {
    let p = [problem(1.0, false, "x"), problem(2.0, false, "x")];
    let o: Vec<_> = p.iter().map(|p| ReverseProblem::new(p, limits()).unwrap()).collect();
    let stages: Vec<_> = o.iter().map(stage).collect();
    let mut study = NewtonContinuation::new(&stages, &[0.0], None).unwrap();
    let empty = study.run(10, AugmentedNewtonBudget { objective_evals: 0, ..budget() }, None).unwrap();
    assert_eq!(empty.stop, ContinuationStop::EntryBudget); assert_eq!(empty.work.objective_evals, 0);
    let stopped = study.run(100, AugmentedNewtonBudget { hessian_products: 2, ..budget() }, None).unwrap();
    assert_eq!(stopped.stage, 1);
    assert_eq!(stopped.stop, ContinuationStop::StageStopped(AugmentedNewtonStop::HessianBudget));
    assert_eq!(stopped.work.hessian_products, 2);
    let final_report = study.run(100, budget(), None).unwrap();
    assert_eq!(final_report.stop, ContinuationStop::Complete);
    assert_eq!(final_report.work.hessian_products, 4);
    assert_eq!(final_report.work.stage_initializations, 2);
}

#[derive(Debug)]
struct Model<'a> { calls: Cell<usize>, deny: Cell<bool>, cancel: Option<&'a CancelGate> }
impl PhysicsModel for Model<'_> {
    fn signature(&self) -> PhysicsSignature<'_> {
        PhysicsSignature { study: "q", parameter_count: 1, parameter_dims: Dims::NONE, residual_dims: Dims::NONE }
    }
    fn value_gradient(&self, x: &[f64], _: Option<&Cx<'_>>) -> Result<PhysicsSample, PhysicsError> {
        self.calls.set(self.calls.get() + 1);
        if self.deny.get() { return Err(PhysicsError::Budget { resource: "test solve", used: 0, limit: 0 }); }
        if let Some(gate) = self.cancel { gate.request(); }
        Ok(PhysicsSample { value: (x[0] - 1.0).powi(2), gradient: vec![2.0 * (x[0] - 1.0)] })
    }
    fn second_order(&self) -> Option<&dyn PhysicsHessian> { Some(self) }
}
impl PhysicsHessian for Model<'_> {
    fn hessian_vector_product(&self, _: &[f64], v: &[f64], _: Option<&Cx<'_>>) -> Result<Vec<f64>, PhysicsError> {
        Ok(vec![2.0 * v[0]])
    }
}
fn physical_problem() -> (Problem, fs_opt::NodeId) {
    let mut b = ProblemBuilder::new(); let v = b.var("x", Manifold::Rn { dim: 1 }, Dims::NONE).unwrap();
    let node = b.pde_residual("q", v, true, Dims::NONE).unwrap();
    b.objective(node, Sense::Minimize, 1.0).unwrap(); (b.finish(), node)
}

#[test]
fn failed_new_model_keeps_old_accepted_tape_and_charges_each_entry_attempt() {
    let p = problem(0.0, false, "x"); let (q, node) = physical_problem();
    let model = Model { calls: Cell::new(0), deny: Cell::new(true), cancel: None };
    let a = ReverseProblem::new(&p, limits()).unwrap();
    let b = ReverseProblem::new_with_physics(&q, limits(), &[PhysicsBinding { node, model: &model }], None).unwrap();
    let mut study = NewtonContinuation::new(&[stage(&a), stage(&b)], &[0.0], None).unwrap();
    for attempted in [2, 3] {
        assert!(matches!(study.run(10, budget(), None), Err(ContinuationError::Stage { stage: 1, .. })));
        assert_eq!(study.active_stage(), 0);
        assert_eq!(study.active().unwrap().inner().snapshot().x, vec![0.0]);
        assert_eq!(study.work().objective_evals, attempted);
        assert_eq!(study.work().stage_initializations, attempted);
        assert!(study.completed_stages().is_empty());
    }
    model.deny.set(false);
    let result = study.run(10, budget(), None).unwrap();
    assert_eq!(result.stop, ContinuationStop::Complete);
    assert_eq!(result.work.stage_initializations, 4);
    assert_eq!(result.work.objective_evals, 5); // first + two failures + next initial + trial
    assert_eq!(model.calls.get(), 4);
}

#[test]
fn cancellation_inside_stage_preparation_cannot_publish_the_new_model() {
    let pool = TilePool::for_host(1, 81); let gate = CancelGate::new_clock_free();
    let p = problem(0.0, false, "x"); let (q, node) = physical_problem();
    let model = Model { calls: Cell::new(0), deny: Cell::new(false), cancel: Some(&gate) };
    let a = ReverseProblem::new(&p, limits()).unwrap();
    let b = ReverseProblem::new_with_physics(&q, limits(), &[PhysicsBinding { node, model: &model }], None).unwrap();
    let mut study = NewtonContinuation::new(&[stage(&a), stage(&b)], &[0.0], None).unwrap();
    study.run(10, AugmentedNewtonBudget { objective_evals: 1, ..budget() }, None).unwrap();
    let before = study.active().unwrap().inner().snapshot();
    pool.arena_pool().scope(|arena| {
        let cx = Cx::new(&gate, arena, StreamKey { seed: 81, kernel_id: 7, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic);
        assert!(matches!(study.run(10, budget(), Some(&cx)), Err(ContinuationError::Stage { stage: 1, .. })));
    });
    assert_eq!(study.active_stage(), 0);
    assert_eq!(study.active().unwrap().inner().snapshot().x, before.x);
    assert_eq!(study.work().objective_evals, 2);
    assert_eq!(study.run(10, budget(), None).unwrap().stop, ContinuationStop::Complete);
}

#[test]
fn constrained_stages_reach_their_own_kkt_points_and_retain_combined_multiplier_work() {
    let p = [problem(1.0, true, "x"), problem(1.5, true, "x")];
    let o: Vec<_> = p.iter().map(|p| ReverseProblem::new(p, limits()).unwrap()).collect();
    let stages: Vec<_> = o.iter().map(stage).collect();
    let mut study = NewtonContinuation::new(&stages, &[0.0], None).unwrap();
    let capped = study.run(300, AugmentedNewtonBudget { multiplier_updates: 0, ..budget() }, None).unwrap();
    assert_eq!(capped.stop, ContinuationStop::StageStopped(AugmentedNewtonStop::MultiplierBudget));
    assert_eq!(capped.stage, 0); assert_eq!(capped.work.multiplier_attempts, 0);
    let result = study.run(600, budget(), None).unwrap();
    assert_eq!(result.stop, ContinuationStop::Complete);
    let last = result.solution.unwrap();
    assert!(last.kkt.within_tolerance(1e-7)); assert!((last.x[0] - 1.5).abs() < 1e-6);
    assert!((study.completed_stages()[0].x[0] - 1.0).abs() < 1e-6);
    assert_eq!(result.work.multiplier_attempts, study.completed_stages()[0].multiplier_attempts + last.multiplier_attempts);
}

#[test]
fn stage_schema_and_tolerances_are_admitted_before_any_initial_evaluation() {
    let a = problem(0.0, false, "x"); let b = problem(0.0, false, "different coordinate");
    let a = ReverseProblem::new(&a, limits()).unwrap(); let b = ReverseProblem::new(&b, limits()).unwrap();
    assert!(NewtonContinuation::new(&[], &[0.0], None).is_err());
    assert!(NewtonContinuation::new(&[stage(&a), stage(&b)], &[0.0], None).is_err());
    assert!(NewtonContinuation::new(&[stage(&a)], &[], None).is_err());
    assert!(NewtonContinuation::new(&[stage(&a)], &[f64::NAN], None).is_err());
    assert!(NewtonContinuation::new(&[NewtonStage { tolerance: f64::NAN, ..stage(&a) }], &[0.0], None).is_err());
}
